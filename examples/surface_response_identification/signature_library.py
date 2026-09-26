#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Signature-only process library from a version-2 preflight manifest.

One model per COUPON of ``Identification.Features`` (decision 85(1), the library contract):
features of one type whose signature TOPOLOGY (the signature with every continuous parameter
removed) is equal and whose parameters agree within the recorded tolerance
(``SignatureParameterToleranceOverR`` = 1e-3 R for offsets / separations / radii,
``SignatureAngleToleranceDegrees`` = 1e-2 deg) are grouped by single linkage into one coupon
whose ``Signature`` is the group's representative (every parameter the midpoint of its range,
in the orientation nearest to the lexicographically smallest member — a set function; the same
rule as ``RepresentativeSignature`` in surfaceresponseidentification.cpp), with placeholder
matrix paths: enough for the patch dry run (``palace --surface-response-preflight`` builds the
patches without reading any matrix; the matcher accepts a model within the tolerance and takes
the nearest), never for a field solve. A dry run with this library must patch every feature of
the manifest exactly once, i.e. cover the whole perimeter minus the recorded exclusions.

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


# Continuous signature entries (surfaceresponseidentification.cpp IsLengthParameter /
# IsAngleParameter) and the recorded tolerances (Identification.Conventions).
LENGTH_KEYS = ("OffsetOverR", "SeparationOverR", "RadiusOverR", "CornerRadiusOverR")
ANGLE_KEYS = ("AngleDegrees", "ArmAnglesDegrees")
PARAMETER_TOLERANCE_OVER_R = 1.0e-3
ANGLE_TOLERANCE_DEGREES = 1.0e-2
LENGTH_QUANTUM_OVER_R = 1.0e-6
ANGLE_QUANTUM_DEGREES = 1.0e-6


def split_parameters(signature):
    """(topology key, lengths, angles) of a signature: the continuous entries replaced by null in
    the serialised topology key (SpatialEdgeCluster: the whole signature, no parameters)."""
    lengths, angles = [], []
    if isinstance(signature, dict) and signature.get("Type") == "SpatialEdgeCluster":
        return json.dumps(signature, sort_keys=True), lengths, angles

    def walk(node):
        if isinstance(node, dict):
            out = {}
            for key in sorted(node):
                value = node[key]
                if key in LENGTH_KEYS or key in ANGLE_KEYS:
                    target = lengths if key in LENGTH_KEYS else angles
                    target.extend(float(v) for v in (value if isinstance(value, list) else [value]))
                    out[key] = None
                else:
                    out[key] = walk(value)
            return out
        if isinstance(node, list):
            return [walk(v) for v in node]
        return node
    return json.dumps(walk(signature), sort_keys=True), lengths, angles


def substitute_parameters(signature, lengths, angles):
    """The signature with its continuous entries replaced in the same traversal order."""
    lengths, angles = list(lengths), list(angles)

    def walk(node):
        if isinstance(node, dict):
            out = {}
            for key in sorted(node):
                value = node[key]
                if key in LENGTH_KEYS or key in ANGLE_KEYS:
                    source, quantum = (lengths, LENGTH_QUANTUM_OVER_R) if key in LENGTH_KEYS else (angles, ANGLE_QUANTUM_DEGREES)
                    if isinstance(value, list):
                        out[key] = [round(source.pop(0) / quantum) * quantum for _ in value]
                    else:
                        out[key] = round(source.pop(0) / quantum) * quantum
                else:
                    out[key] = walk(value)
            return out
        if isinstance(node, list):
            return [walk(v) for v in node]
        return node
    out = walk(signature)
    assert not lengths and not angles
    return out


def mirror_translational(signature):
    """The mirror orientation of an ``Edges`` signature (reversed, gap sides negated, offsets from
    the top, conductors relabelled by first appearance); other signatures unchanged."""
    if not isinstance(signature, dict) or not signature.get("Edges"):
        return signature
    edges = signature["Edges"]
    span = float(edges[-1]["OffsetOverR"])
    labels, out = {}, []
    for edge in reversed(edges):
        e = dict(edge)
        e["OffsetOverR"] = round((span - float(edge["OffsetOverR"])) / LENGTH_QUANTUM_OVER_R) * LENGTH_QUANTUM_OVER_R
        e["GapSide"] = -int(edge["GapSide"])
        e["Conductor"] = labels.setdefault(int(edge["Conductor"]), len(labels) + 1)
        out.append(e)
    mirror = dict(signature)
    mirror["Edges"] = out
    return mirror


def signature_deviation(a, b):
    """max |difference| / tolerance over the parameters (both orientations of b; the smaller), or
    None when the topologies differ. Within tolerance iff <= 1."""
    ta, la, aa = split_parameters(a)
    best = None
    for candidate in (b, mirror_translational(b)):
        tb, lb, ab = split_parameters(candidate)
        if ta != tb or len(la) != len(lb) or len(aa) != len(ab):
            continue
        d = max([abs(x - y) / PARAMETER_TOLERANCE_OVER_R for x, y in zip(la, lb)] + [abs(x - y) / ANGLE_TOLERANCE_DEGREES for x, y in zip(aa, ab)] + [0.0])
        best = d if best is None else min(best, d)
    return best


def representative_signature(signatures):
    """The coupon's signature of a group of instances: the lexicographically smallest instance
    with every parameter replaced by the midpoint of its range over the group (each instance in
    the orientation nearest to that lead)."""
    lead = min(signatures, key=lambda s: json.dumps(s, sort_keys=True))
    t_lead, l_lead, a_lead = split_parameters(lead)
    lo_l, hi_l, lo_a, hi_a = list(l_lead), list(l_lead), list(a_lead), list(a_lead)
    for signature in signatures:
        best = None
        for candidate in (signature, mirror_translational(signature)):
            t, lens, angs = split_parameters(candidate)
            if t != t_lead or len(lens) != len(l_lead) or len(angs) != len(a_lead):
                continue
            d = max([abs(x - y) for x, y in zip(lens, l_lead)] + [abs(x - y) for x, y in zip(angs, a_lead)] + [0.0])
            if best is None or d < best[0]:
                best = (d, lens, angs)
        assert best is not None, "representative over different topologies"
        lo_l = [min(x, y) for x, y in zip(lo_l, best[1])]
        hi_l = [max(x, y) for x, y in zip(hi_l, best[1])]
        lo_a = [min(x, y) for x, y in zip(lo_a, best[2])]
        hi_a = [max(x, y) for x, y in zip(hi_a, best[2])]
    return substitute_parameters(lead, [0.5 * (x + y) for x, y in zip(lo_l, hi_l)], [0.5 * (x + y) for x, y in zip(lo_a, hi_a)])


def group_features(features):
    """Coupons of a feature list: per (type, topology key incl. both orientations) the distinct
    signatures grouped by single linkage at the tolerance (order independent); returns a list of
    (representative signature, member features, parameter spread)."""
    bases = {}
    for feature in features:
        signature = feature["Signature"]
        keys = [feature["Type"] + "|" + split_parameters(s)[0] for s in (signature, mirror_translational(signature))]
        bases.setdefault(min(keys), {}).setdefault(json.dumps(signature, sort_keys=True), (signature, []))[1].append(feature)
    groups = []
    for base_key in sorted(bases):
        distinct = [bases[base_key][k] for k in sorted(bases[base_key])]
        parent = list(range(len(distinct)))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i
        for i in range(len(distinct)):
            for j in range(i + 1, len(distinct)):
                d = signature_deviation(distinct[i][0], distinct[j][0])
                if d is not None and d <= 1.0:
                    parent[find(i)] = find(j)
        members = {}
        for i in range(len(distinct)):
            members.setdefault(find(i), []).append(i)
        for root in sorted(members, key=lambda r: json.dumps(distinct[r][0], sort_keys=True)):
            signatures = [distinct[i][0] for i in members[root]]
            representative = representative_signature(signatures)
            spread = max((signature_deviation(representative, s) or 0.0) for s in signatures)
            groups.append((representative, [f for i in members[root] for f in distinct[i][1]], spread))
    return groups


def signature_hash(signature):
    import hashlib
    return hashlib.sha256(json.dumps(signature, separators=(",", ":"), sort_keys=True).encode()).hexdigest()


def build_signature_library(manifest, name="signature-only", matrix_directory="signature-only-matrices"):
    identification = manifest["Identification"]
    radius = float(identification["MatchingRadius"])
    models = []
    modelled = [f for f in identification["Features"] if f["Type"] not in UNMODELLED_TYPES]
    for representative, members, spread in group_features(modelled):
        feature = {"Type": members[0]["Type"], "Signature": representative}
        model_name = f"{feature['Type']}-{signature_hash(representative)[:12]}"
        model = {
            "Name": model_name,
            "Topology": feature["Type"],
            "Signature": representative,
            "Instances": len(members),
            "ParameterSpread": spread,
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
