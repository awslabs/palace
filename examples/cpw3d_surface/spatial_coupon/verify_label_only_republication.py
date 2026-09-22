#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Verify that a re-published production identity mesh differs from the previous
production identity in the radial MA shell labels only (supervisor decision 61a): the
$Nodes block byte-identical, every element identical in type / node list / tag count,
every tag identical except the (physical, elementary) pair of the MA elements whose new
physical label is a shell of the old label (shell_parent(new) == old, the shell label in
the bound census), the relabeled count equal to the census's, and - when the new root's
transform receipt is given - the receipt's ParentLabeledMeshSHA256 equal to the previous
identity's SHA-256 (the publication before the relabel IS the previous production mesh).

usage: verify_label_only_republication.py --previous IDENTITY.msh --current IDENTITY.msh
       --census IDENTITY.msh.radial-shells.json [--receipt identity-receipt.json] --out RECORD.json
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from mixed_mesh import SHELL_LABEL_STRIDE, shell_parent  # noqa: E402
from relabel_radial_ma_shells import read_msh22_binary  # noqa: E402


class RepublicationError(ValueError):
    """The two meshes differ beyond the MA shell labels."""


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def compare(previous_bytes, current_bytes, census):
    """The label-only record of the two MSH 2.2 binaries; raises RepublicationError."""
    before = read_msh22_binary(previous_bytes)
    after = read_msh22_binary(current_bytes)
    if previous_bytes[before["NodesSpan"][0]:before["NodesSpan"][1]] != current_bytes[after["NodesSpan"][0]:after["NodesSpan"][1]]:
        raise RepublicationError("the $Nodes block differs")
    if len(before["Elements"]) != len(after["Elements"]):
        raise RepublicationError(f"the element count differs ({len(before['Elements'])} vs {len(after['Elements'])})")
    shell_labels = {int(shell["Label"]): int(shell["Parent"]) for shell in census["Shells"]}
    relabeled = 0
    per_parent = {}
    for index, (old, new) in enumerate(zip(before["Elements"], after["Elements"])):
        if old[0] != new[0] or old[3] != new[3] or len(old[2]) != len(new[2]):
            raise RepublicationError(f"element {index} differs in type, node list or tag count")
        if old[2] == new[2]:
            continue
        new_label, old_label = new[2][0], old[2][0]
        if new_label < SHELL_LABEL_STRIDE or shell_labels.get(new_label) != old_label or shell_parent(new_label) != old_label:
            raise RepublicationError(f"element {index} changed its label {old_label} -> {new_label}, not a shell of the old label")
        if old[2][2:] != new[2][2:]:
            raise RepublicationError(f"MA element {index} changed a tag beyond the (physical, elementary) pair")
        relabeled += 1
        per_parent[old_label] = per_parent.get(old_label, 0) + 1
    expected = census.get("LabelOnly", {}).get("RelabeledElements")
    if expected is not None and expected != relabeled:
        raise RepublicationError(f"{relabeled} elements relabeled, the census records {expected}")
    old_physical = {(dimension, tag) for dimension, tag, _ in before["PhysicalNames"]}
    new_physical = {(dimension, tag) for dimension, tag, _ in after["PhysicalNames"]}
    return {"NodesBlockIdentical": True, "ElementsIdenticalApartFromShellLabels": True, "Elements": len(after["Elements"]),
            "RelabeledElements": relabeled, "RelabeledPerParent": {str(k): v for k, v in sorted(per_parent.items())},
            "CensusRelabeledElements": expected, "ShellLabels": sorted(shell_labels),
            "PhysicalNamesAdded": sorted(tag for dimension, tag in new_physical - old_physical),
            "PhysicalNamesRemoved": sorted(tag for dimension, tag in old_physical - new_physical)}


def verify(previous, current, census_path, receipt_path=None):
    previous, current, census_path = Path(previous), Path(current), Path(census_path)
    census = json.loads(census_path.read_text())
    current_bytes = current.read_bytes()
    if census.get("Mesh", {}).get("SHA256") != hashlib.sha256(current_bytes).hexdigest():
        raise RepublicationError("the census binds another mesh than the current identity")
    previous_bytes = previous.read_bytes()
    record = compare(previous_bytes, current_bytes, census)
    record.update({"Previous": {"Path": str(previous), "SHA256": hashlib.sha256(previous_bytes).hexdigest(), "Bytes": len(previous_bytes)},
                   "Current": {"Path": str(current), "SHA256": hashlib.sha256(current_bytes).hexdigest(), "Bytes": len(current_bytes)},
                   "Census": {"Path": str(census_path), "SHA256": sha256(census_path)}})
    if receipt_path is not None:
        receipt = json.loads(Path(receipt_path).read_text())
        parent = receipt.get("ParentLabeledMeshSHA256")
        record["Receipt"] = {"Path": str(receipt_path), "ParentLabeledMeshSHA256": parent,
                             "ParentLabeledEqualsPrevious": parent == record["Previous"]["SHA256"]}
        if parent != record["Previous"]["SHA256"]:
            raise RepublicationError("the receipt's ParentLabeledMeshSHA256 is not the previous identity's SHA-256")
    record["Verdict"] = "labels-only"
    record["Rule"] = ("the current identity equals the previous production identity in the $Nodes block and in every "
                      "element apart from the (physical, elementary) pair of the MA elements, whose new label is a shell "
                      "(10000 x ordinal + parent) of the old label bound in the census")
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--previous", type=Path, required=True)
    parser.add_argument("--current", type=Path, required=True)
    parser.add_argument("--census", type=Path, required=True)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args(argv)
    record = verify(args.previous, args.current, args.census, args.receipt)
    args.out.write_text(json.dumps(record, indent=2) + "\n")
    print(f"{record['Verdict']}: {record['RelabeledElements']} MA elements relabeled of {record['Elements']}; "
          f"{record['Previous']['SHA256'][:12]} -> {record['Current']['SHA256'][:12]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
