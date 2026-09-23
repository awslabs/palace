#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Derive a coupon's frozen semantic mesh contract (semantic-contract.json) from its
immutable source inputs; no number is authored by hand.

Rules (the ones the four-edge / ten-edge contracts were derived with, made
executable):

- materials and label families follow the frozen fabricated spatial-coupon producer
  (mesh_spatial_coupon.jl): substrate 1 / vacuum 2; matching surface 1; the
  substrate-vacuum interface of slot s is 3000 + s where it lies in the metal plane
  (un-etched) and 3100 + s where it is recessed (trench floor and walls); the
  metal-substrate (MS) surface of conductor c in slot s is 5000 + 100 s + c and the
  metal-vacuum (MA) surface 6000 + 100 s + c;
- slots and conductors come from Models[0].Edges of process-library.json
  (InterfaceSlot, Conductor) and must agree with the signature's Slot / Conductor
  columns;
- semantic corners are the plan-view boundary vertices classified Physical, in
  file order, at the vertex's Plane height;
- FeatureTopology is semantic_mesh_contract.derive_feature_topology (the finite
  oriented signature segments against the corners and the boundary classes);
- whether a slot's un-etched plane (3000 + s) exists is a producer outcome of the
  bound etch footprint (a device footprint or the producer-default collars) and
  the coupon box, not of the inputs alone: it is taken from the InterfaceAreas
  labels of a gmsh-build census of the same inputs (--build-census; a probe build
  with the provisional contract this tool writes without it).  Every census label
  must belong to a derived family and every family without an un-etched
  alternative must be present, otherwise the derivation fails closed.  The census
  digest and its labels are recorded under Derivation.

Two-pass workflow (by design): the un-etched plane set is a producer outcome, so a
contract is derived twice - (1) without --build-census (PROVISIONAL: the required
label set only), (2) a probe of the same inputs under the production options with that
provisional contract - the labels-only mesher pass of run_gmsh_only_case.py
--labels-only (decision 62(2): the physical surface groups are assigned on the CAD
entities before any mesh generation, so the pass stops there; register_case.py) or a
full census-only build - then this tool again with --build-census PROBE/build-census.json.
The probe census is NOT validated by the stage contract here (only its InterfaceAreas
labels are read, and only labels inside the derived families are accepted); the
production build that follows binds the final contract and validates its own census
against it (fail closed on a label mismatch).  Because the contract enters the build
through its SemanticCorners only, the label set of the labels-only pass (and the
gmsh-build.msh of a full probe) equals the production one.

The inputs are SOURCE_DIR/mesh-signature.csv, plan-view-boundary.csv and
process-library.json unless the case binds other file names (--signature,
--boundary, --process-library).  A case without a process library (the synthetic
repository fixtures) takes its slot / conductor pairs from the signature alone and
records that source under Derivation.SlotConductorSource.

usage: derive_semantic_contract.py SOURCE_DIR OUTPUT [--build-census CENSUS]
       [--signature PATH] [--boundary PATH] [--process-library PATH]
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from semantic_mesh_contract import derive_feature_topology, validate_semantic_contract

SIGNATURE = "mesh-signature.csv"
BOUNDARY = "plan-view-boundary.csv"
PROCESS_LIBRARY = "process-library.json"
MATCHING_SURFACE = 1
UNETCHED_BASE, ETCHED_BASE, THIN_BASE, MS_BASE, MA_BASE = 3000, 3100, 4000, 5000, 6000
COUPON_KINDS = ("fabricated", "thin")
THIN_RULE = ("thin kind (decision 66): the metal is a zero-thickness PEC sheet in the process plane, "
             "so there is no trench (no 3100 + slot label; the substrate-vacuum plane 3000 + slot is "
             "always present) and the sheet of conductor c in slot s is ONE physical surface "
             "4000 + 100 s + c carrying both sides - MS below and MA above (the producer's thin "
             "convention: make_config lists the 4000 family under both interface types)")
RULES = ("materials and physical label families follow the frozen fabricated spatial-coupon "
         "producer (substrate 1 / vacuum 2; matching surface 1; substrate-vacuum interface "
         "3000 + slot in the metal plane (un-etched) and 3100 + slot recessed; MS 5000 + 100 slot "
         "+ conductor; MA 6000 + 100 slot + conductor); slots/conductors come from "
         "Models[0].Edges and agree with the signature (from the signature alone when no "
         "process library is bound); semantic corners are plan-view vertices "
         "classified Physical; the un-etched plane labels present are those of the recorded "
         "gmsh-build census of the same inputs (a producer outcome of the etch footprint and "
         "the coupon box), every census label belonging to a derived family")


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


PROCESS_LIBRARY_PAIR_SOURCE = "process-library.json Models[0].Edges, equal to the signature's Slot / Conductor pairs"
SIGNATURE_PAIR_SOURCE = "mesh signature Slot / Conductor pairs (no process library is bound)"


SUPPORTED_TOPOLOGY = "SpatialEdgeCluster"


def process_library_pairs(process_library):
    """(slot, conductor) of every edge of the library's single SpatialEdgeCluster model.
    The label families assume one model whose Edges carry InterfaceSlot / Conductor;
    a multi-model library or another topology (e.g. an Arms model) is a contract
    error, not a lookup failure."""
    library = json.loads(Path(process_library).read_text())
    models = library.get("Models")
    if not isinstance(models, list) or len(models) != 1 or not isinstance(models[0], dict):
        raise ValueError(f"{process_library}: the process library must carry exactly one model "
                         f"(Models), found {len(models) if isinstance(models, list) else 'none'}")
    model = models[0]
    if model.get("Topology") != SUPPORTED_TOPOLOGY:
        raise ValueError(f"{process_library}: model topology {model.get('Topology')!r} is not "
                         f"{SUPPORTED_TOPOLOGY}; the label families are derived for edge clusters only")
    edges = model.get("Edges")
    if (not isinstance(edges, list) or not edges or
            any(not isinstance(edge, dict) or "InterfaceSlot" not in edge or "Conductor" not in edge
                for edge in edges)):
        raise ValueError(f"{process_library}: every Models[0].Edges entry must carry InterfaceSlot "
                         f"and Conductor")
    return {(int(edge["InterfaceSlot"]), int(edge["Conductor"])) for edge in edges}


def slot_conductor_pairs(signature, process_library):
    """Sorted (slot, conductor) pairs: Models[0].Edges of the process library when one is
    bound (it must agree with the signature), the signature's own pairs otherwise."""
    with Path(signature).open(newline="") as stream:
        signature_pairs = sorted({(int(row["Slot"]), int(row["Conductor"]))
                                  for row in csv.DictReader(stream)})
    if process_library is None:
        pairs = signature_pairs
    else:
        pairs = sorted(process_library_pairs(process_library))
        if pairs != signature_pairs:
            raise ValueError(f"Models[0].Edges slots/conductors {pairs} differ from the signature's "
                             f"{signature_pairs}")
    if not pairs:
        raise ValueError("the signature has no rows")
    if any(slot < 0 or slot > 99 or conductor < 1 or conductor > 99 for slot, conductor in pairs):
        raise ValueError("slot must lie in 0..99 and conductor in 1..99 for the label families")
    return pairs


def semantic_corners(boundary):
    with Path(boundary).open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    if not rows or not {"Class", "X", "Y", "Plane"} <= set(rows[0]):
        raise ValueError("plan-view boundary lacks Class / X / Y / Plane columns")
    if any(row["Class"] not in ("Physical", "Continuation") for row in rows):
        raise ValueError("plan-view boundary has an unknown vertex class")
    corners = [[float(row["X"]), float(row["Y"]), float(row["Plane"])]
               for row in rows if row["Class"] == "Physical"]
    if not corners:
        raise ValueError("plan-view boundary classifies no Physical vertex")
    return corners


def label_families(pairs, kind="fabricated"):
    """Attribute -> (role, adjacent materials, optional) of every derivable label of a
    coupon kind (fabricated: un-etched / etched planes, MS, MA; thin: the plane and the
    two-sided sheet, THIN_RULE)."""
    if kind not in COUPON_KINDS:
        raise ValueError(f"coupon kind must be one of {COUPON_KINDS}, not {kind!r}")
    slots = sorted({slot for slot, _ in pairs})
    families = {MATCHING_SURFACE: ("matching-surface", [1, 2], False)}
    if kind == "thin":
        for slot in slots:
            families[UNETCHED_BASE + slot] = (f"substrate-vacuum-slot-{slot}", [1, 2], False)
        for slot, conductor in pairs:
            families[THIN_BASE + 100 * slot + conductor] = (f"conductor-{conductor}-slot-{slot}-thin-sheet", [1, 2], False)
        return families
    for slot in slots:
        families[UNETCHED_BASE + slot] = (f"un-etched-substrate-vacuum-slot-{slot}", [1, 2], True)
        families[ETCHED_BASE + slot] = (f"etched-substrate-vacuum-slot-{slot}", [1, 2], False)
    for slot, conductor in pairs:
        families[MS_BASE + 100 * slot + conductor] = (f"conductor-{conductor}-slot-{slot}-ms", [1], False)
    for slot, conductor in pairs:
        families[MA_BASE + 100 * slot + conductor] = (f"conductor-{conductor}-slot-{slot}-ma", [2], False)
    return families


def census_labels(census_path):
    census = json.loads(Path(census_path).read_text())
    areas = census.get("InterfaceAreas")
    if not isinstance(areas, list) or not areas:
        raise ValueError("build census records no InterfaceAreas")
    return sorted({int(item["Attribute"]) for item in areas})


def derive(source, build_census=None, *, signature=None, boundary=None, process_library=None,
           kind="fabricated"):
    source = Path(source)
    signature = Path(signature) if signature is not None else source / SIGNATURE
    boundary = Path(boundary) if boundary is not None else source / BOUNDARY
    if process_library is None and (source / PROCESS_LIBRARY).exists():
        process_library = source / PROCESS_LIBRARY
    pairs = slot_conductor_pairs(signature, process_library)
    families = label_families(pairs, kind)
    required = {attribute for attribute, (_, _, optional) in families.items() if not optional}
    if build_census is not None:
        labels = census_labels(build_census)
        unknown = sorted(set(labels) - set(families))
        missing = sorted(required - set(labels))
        if unknown or missing:
            raise ValueError(f"build census labels outside the derived families {unknown} or "
                             f"required labels missing {missing}")
        present = set(labels)
    else:
        present = required
    boundary_labels = []
    for attribute in sorted(families):
        if attribute not in present:
            continue
        role, adjacent, _ = families[attribute]
        item = {"Attribute": attribute, "Role": role, "AdjacentMaterials": adjacent,
                "Protected": True}
        if attribute == MATCHING_SURFACE:
            item["AdjacentMaterialSets"] = [[1], [2]]
        boundary_labels.append(item)
    roles = [item["Role"] for item in boundary_labels]
    corners = semantic_corners(boundary)
    derivation = {"ProcessLibrarySHA256": (sha256(process_library) if process_library is not None
                                           else None),
                  "SlotConductorSource": (PROCESS_LIBRARY_PAIR_SOURCE if process_library is not None
                                          else SIGNATURE_PAIR_SOURCE),
                  "PlanViewBoundarySHA256": sha256(boundary),
                  "SignatureSHA256": sha256(signature),
                  "CouponKind": kind,
                  "Rules": RULES if kind == "fabricated" else RULES + "; " + THIN_RULE}
    if build_census is not None:
        derivation["BuildCensusSHA256"] = sha256(build_census)
        derivation["BuildCensusInterfaceLabels"] = labels
    elif any(optional for _, _, optional in families.values()):
        derivation["Provisional"] = ("un-etched plane labels unconfirmed: rebuild this contract "
                                     "with --build-census before freezing it")
    else:
        derivation["LabelSetRule"] = ("no optional label in the derived families (thin kind: the plane and the "
                                      "sheet are always present): the label set is final from the inputs alone")
    contract = {
        "Version": 1,
        "Derivation": derivation,
        "VolumeMaterials": [{"Attribute": 1, "Material": "substrate"},
                            {"Attribute": 2, "Material": "vacuum"}],
        "BoundaryLabels": boundary_labels,
        "SemanticCorners": corners,
        "ProtectedSupports": roles,
        "MetricSurfaceRoles": [role for role in roles if role.endswith(("-ms", "-ma", "-thin-sheet"))],
        "UnmatchedPolicy": "Error",
        "CutSurfaceRoles": ["matching-surface"],
        "FeatureTopology": derive_feature_topology(signature, boundary, corners),
    }
    return validate_semantic_contract(contract)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--build-census", type=Path)
    parser.add_argument("--signature", type=Path, help="default SOURCE_DIR/mesh-signature.csv")
    parser.add_argument("--boundary", type=Path, help="default SOURCE_DIR/plan-view-boundary.csv")
    parser.add_argument("--process-library", type=Path,
                        help="default SOURCE_DIR/process-library.json when it exists")
    parser.add_argument("--kind", choices=COUPON_KINDS, default="fabricated",
                        help="coupon kind: fabricated (default) or thin (the two-sided sheet family, decision 66)")
    args = parser.parse_args()
    contract = derive(args.source, args.build_census, signature=args.signature,
                      boundary=args.boundary, process_library=args.process_library, kind=args.kind)
    args.output.write_text(json.dumps(contract, indent=2) + "\n")
    print(f"{args.output}: {len(contract['BoundaryLabels'])} labels, "
          f"{len(contract['SemanticCorners'])} semantic corners"
          + (" (PROVISIONAL)" if "Provisional" in contract["Derivation"] else ""))


if __name__ == "__main__":
    main()
