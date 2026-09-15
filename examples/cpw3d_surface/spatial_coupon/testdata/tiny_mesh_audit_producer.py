#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tiny deterministic producer used to exercise evidence provenance tests."""
import argparse
import json
from pathlib import Path

from semantic_mesh_contract import load_semantic_contract


def produce(contract_path, case, variant, output, mesh, *, dofs=120, features=2,
            subdivisions=2, invariant=4.0):
    contract = load_semantic_contract(contract_path)
    mesh.write_bytes((case + "\n" + variant + "\n").encode())
    report = {
        "CaseId": case, "Variant": variant,
        "ActualVolumeMaterials": contract["VolumeMaterials"],
        "ActualBoundaryAttributes": [item["Attribute"] for item in contract["BoundaryLabels"]],
        "ActualAdjacency": {str(item["Attribute"]): item["AdjacentMaterials"]
                            for item in contract["BoundaryLabels"]},
        "OwnershipClosure": {"UnmatchedPolicy": "Error", "Unmatched": 0,
                             "Overlaps": 0, "Exhaustive": True},
        "ActualSemanticCorners": contract["SemanticCorners"],
        "ProtectedSurfaces": {"Actual": contract["ProtectedSupports"], "Changed": 0},
        "AchievedAnisotropy": {"Samples": 8, "NormalTarget": 0.01,
                               "Transverse1P90": 0.012, "Transverse2P90": 0.013,
                               "TangentialP50": 0.04},
        "TraceDiagonal": {"GlobalDiagonalBands": 0},
        "Resources": {"ExitCode": 0, "Seconds": 0.01, "PeakRSSGiB": 0.01,
                      "Elements": 12},
        "MeshQuality": {"Samples": 12, "MinimumScaledJacobian": 0.2,
                        "MaximumJacobianCondition": 8.0},
        "Complexity": {"H1DOFs": dofs, "FeatureCount": features,
                       "CADSubdivisionCount": subdivisions},
        "ComparisonInvariants": {"Volume": invariant, "ProtectedArea": 2.0},
    }
    output.write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("contract", type=Path)
    parser.add_argument("case")
    parser.add_argument("variant")
    parser.add_argument("output", type=Path)
    parser.add_argument("mesh", type=Path)
    args = parser.parse_args()
    produce(args.contract, args.case, args.variant, args.output, args.mesh)


if __name__ == "__main__":
    main()
