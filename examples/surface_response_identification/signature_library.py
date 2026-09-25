#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Signature-only process library from a version-2 preflight manifest.

One model per distinct feature hash of ``Identification.Features``, carrying the feature's
canonical ``Signature`` (the key-based matching contract of SURFACE-RESPONSE-IDENTIFICATION.md)
and placeholder matrix paths: enough for the patch dry run (``palace --surface-response-preflight``
builds the patches without reading any matrix), never for a field solve. A dry run with this
library must patch every feature of the manifest exactly once, i.e. cover the whole perimeter
minus the recorded exclusions.

    python3 -m surface_response_identification.signature_library MANIFEST OUTPUT.json
"""

import argparse
import json
import os
import sys

# Feature types the library cannot model (an UnclassifiedParallelPair is an exclusion
# described as a feature so that the partition stays exact).
UNMODELLED_TYPES = ("UnclassifiedParallelPair", "CurvedUnclassifiedParallelPair")

PAIR_TYPES = ("SameConductorGap", "DifferentConductorGap", "SameConductorStrip", "CurvedSameConductorGap", "CurvedDifferentConductorGap", "CurvedSameConductorStrip")
LONGITUDINAL_TYPES = ("IsolatedEdge", "CurvedEdge", "ParallelEdgeCluster") + PAIR_TYPES
VERTEX_TYPES = ("ConvexCorner", "ConcaveCorner", "Endpoint", "Junction")


def conductor_count(feature):
    """Canonical conductor labels of a feature (1 for single-conductor features)."""
    signature = feature["Signature"]
    if feature["Type"] == "SpatialEdgeCluster":
        return max(int(p["Conductor"]) for p in signature["Portions"])
    if feature["Type"] == "ParallelEdgeCluster" or feature["Type"] in PAIR_TYPES:
        return max(int(e["Conductor"]) for e in signature["Edges"])
    return 1


def build_signature_library(manifest, name="signature-only", matrix_directory="signature-only-matrices"):
    identification = manifest["Identification"]
    radius = float(identification["MatchingRadius"])
    models = []
    seen = set()
    for feature in sorted(identification["Features"], key=lambda f: (f["Type"], f["Hash"])):
        if feature["Type"] in UNMODELLED_TYPES or feature["Hash"] in seen:
            continue
        seen.add(feature["Hash"])
        model_name = f"{feature['Type']}-{feature['Hash'][:12]}"
        model = {
            "Name": model_name,
            "Topology": feature["Type"],
            "Signature": feature["Signature"],
            "FabricatedMatrix": f"{matrix_directory}/{model_name}-fabricated.csv",
            "ThinMatrix": f"{matrix_directory}/{model_name}-thin.csv",
            "BasisPoints": f"{matrix_directory}/{model_name}-basis-points.csv",
        }
        signature = feature["Signature"]
        if feature["Type"] in LONGITUDINAL_TYPES:
            # Longitudinal coupon depth: the matching radius (any positive depth; the dry run
            # records it with every patch weight).
            model["CouponDepth"] = radius
        # Version-1 geometry parameters derived from the signature (the library reader
        # validates them; the matching itself is by Signature).
        if feature["Type"] in PAIR_TYPES:
            model["Separation"] = float(signature["SeparationOverR"]) * radius
        elif feature["Type"] in ("ConvexCorner", "ConcaveCorner"):
            model["Angle"] = float(signature["AngleDegrees"])
            model["CornerRadius"] = float(signature["CornerRadiusOverR"]) * radius
        elif feature["Type"] == "Junction":
            angles, total = [], 0.0
            for difference in signature["ArmAnglesDegrees"]:
                angles.append(total)
                total += float(difference)
            model["ArmAngles"] = angles
        elif feature["Type"] == "ParallelEdgeCluster":
            model["Edges"] = [{"Offset": float(e["OffsetOverR"]) * radius, "GapDirection": int(e["GapSide"]), "Conductor": int(e["Conductor"])} for e in signature["Edges"]]
        conductors = conductor_count(feature)
        if conductors > 1:
            # One reference per canonical conductor label; positions are placeholders.
            model["ConductorReferences"] = [[0.0, 0.0, -radius * (k + 1)] for k in range(conductors)]
        if feature["Type"] != "SpatialEdgeCluster":
            law = json.loads(feature["Signature"]["Law"]) if "Law" in feature["Signature"] else {"Type": "PEC"}
            if law != {"Type": "PEC"}:
                model["BoundaryCondition"] = law
        models.append(model)
    return {"Version": 2, "Name": name, "MatchingRadius": radius, "TraceLiftVersion": 2, "Models": models}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("manifest")
    parser.add_argument("output")
    parser.add_argument("--name", default="signature-only")
    args = parser.parse_args(argv)
    with open(args.manifest) as source:
        manifest = json.load(source)
    library = build_signature_library(manifest, args.name)
    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w") as target:
        json.dump(library, target, indent=1)
    print(f"{len(library['Models'])} signature-only models -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
