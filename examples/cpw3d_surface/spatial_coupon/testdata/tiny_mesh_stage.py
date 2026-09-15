#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Deterministic stage executable for bounded mesh-DAG regression tests."""
import argparse
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import meshio
import numpy as np

from transform_coupon_source_contract import (read_transform, transform_semantic_contract,
                                                transformed_supports)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rewrite(source, output, binary):
    meshio.write(output, meshio.read(source), file_format="gmsh22", binary=binary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("metric", "adapt", "restore", "publish", "rigid"))
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--metric", type=Path)
    parser.add_argument("--mmg-seed", type=Path)
    parser.add_argument("--pins", type=Path)
    parser.add_argument("--recipe", type=Path)
    parser.add_argument("--fixed-triangles", type=Path)
    parser.add_argument("--ownership", type=Path)
    parser.add_argument("--source-local-output", type=Path)
    parser.add_argument("--semantic-contract", type=Path)
    parser.add_argument("--transformed-supports", type=Path)
    parser.add_argument("--ownership-quadrature", type=Path)
    parser.add_argument("--transform", type=Path)
    parser.add_argument("--source-directory", type=Path)
    parser.add_argument("--semantic-input", type=Path)
    parser.add_argument("--signature", type=Path)
    parser.add_argument("--boundary", type=Path)
    parser.add_argument("--mask", type=Path)
    parser.add_argument("--canonical-build-record", type=Path)
    parser.add_argument("--transformed-semantic-output", type=Path)
    parser.add_argument("--transformed-supports-output", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--ownership-runtime", type=Path)
    parser.add_argument("--ownership-auditor", type=Path)
    args = parser.parse_args()
    if args.stage == "metric":
        if (args.mmg_seed is None or args.pins is None or args.recipe is None or
                args.fixed_triangles is None):
            parser.error("metric requires MMG seed, pins, fixed triangles, and recipe")
        args.output.write_text(json.dumps({"SeedSHA256": digest(args.source),
                                          "Metric": [1.0, 0.0, 1.0]}) + "\n")
        rewrite(args.source, args.mmg_seed, False)
        args.pins.write_text("1\n")
        args.fixed_triangles.write_text("1\n")
        if args.semantic_contract is None or args.transformed_supports is None:
            parser.error("metric requires transformed semantic/support inputs")
        args.recipe.write_text(json.dumps({
            "SeedSHA256": digest(args.source),
            "SemanticContract": json.loads(args.semantic_contract.read_text()),
            "TransformedSupportsArtifact": str(args.transformed_supports.resolve()),
            "TransformedSupportsSHA256": digest(args.transformed_supports),
            "TransformedSupports": json.loads(args.transformed_supports.read_text()),
            "FarSize": 1.0,
            "FarFieldBudgetPolicy": {"Name": "seed-fraction-far-field-v1",
                "RequestedFarSize": 1.0, "Pressure": 1.0,
                "EffectiveFarSize": 1.0}}) + "\n")
    elif args.stage == "adapt":
        if (args.metric is None or not args.metric.is_file() or
                args.pins is None or not args.pins.is_file() or
                args.fixed_triangles is None or not args.fixed_triangles.is_file()):
            parser.error("adapt requires metric, pins, and fixed triangles")
        rewrite(args.source, args.output, True)
    elif args.stage == "restore":
        if (args.recipe is None or not args.recipe.is_file() or
                args.source_local_output is None):
            parser.error("restore requires --recipe and --source-local-output")
        rewrite(args.source, args.source_local_output, False)
        rewrite(args.source, args.output, False)
    elif args.stage == "publish":
        if args.ownership is None or args.ownership_quadrature is None:
            parser.error("publish requires ownership and quadrature ownership")
        rewrite(args.source, args.output, True)
        header = (
            "attribute,elements,area,ambiguous_area,ambiguous_fraction,"
            "unresolved_elements,unresolved_area,unresolved_fraction,quadrature_rule,"
            "quadrature_order,quadrature_points,quadrature_whole_measure,"
            "quadrature_owned_measure,quadrature_relative_closure,"
            "quadrature_closure_tolerance,quadrature_unmatched,quadrature_overlaps,"
            "quadrature_positive_weights\n")
        suffix = ",Gauss4,4,11,3,3,0,1e-12,0,0,1\n"
        args.ownership.write_text(
            header + "1,7,1,0,0,0,0,0" + suffix +
            "2,3,1,0,0,0,0,0" + suffix + "3,1,1,0,0,0,0,0" + suffix)
        args.ownership_quadrature.write_text("attribute,measure\n2,2\n3,1\n")
    else:
        required = (args.transform, args.source_directory, args.semantic_input,
                    args.signature, args.boundary, args.mask, args.canonical_build_record,
                    args.transformed_semantic_output, args.transformed_supports_output,
                    args.receipt, args.ownership, args.ownership_quadrature,
                    args.ownership_runtime, args.ownership_auditor)
        if any(value is None for value in required):
            parser.error("rigid publication inputs and outputs are required")
        matrix = read_transform(args.transform)
        mesh = meshio.read(args.source)
        mesh.points = mesh.points @ np.asarray(matrix)[:3, :3].T + np.asarray(matrix)[:3, 3]
        meshio.write(args.output, mesh, file_format="gmsh22", binary=True)
        semantic_hash = digest(args.semantic_input); transform_hash = digest(args.transform)
        semantic = transform_semantic_contract(
            json.loads(args.semantic_input.read_text()), matrix)
        semantic["SourceSemanticContractSHA256"] = semantic_hash
        semantic["CanonicalTransformSHA256"] = transform_hash
        args.transformed_semantic_output.write_text(json.dumps(semantic) + "\n")
        supports = transformed_supports(
            args.source_directory, matrix, signature=args.signature,
            boundary=args.boundary, mask=args.mask, transform_sha256=transform_hash,
            semantic_sha256=semantic_hash)
        args.transformed_supports_output.write_text(json.dumps(supports) + "\n")
        header = ("attribute,elements,ambiguous_fraction,unresolved_elements,"
                  "unresolved_fraction,quadrature_rule,quadrature_order,quadrature_points,"
                  "quadrature_whole_measure,quadrature_owned_measure,"
                  "quadrature_relative_closure,quadrature_closure_tolerance,"
                  "quadrature_unmatched,quadrature_overlaps,quadrature_positive_weights\n")
        suffix = ",1,0,0,0,Gauss4,4,11,3,3,0,1e-12,0,0,1\n"
        args.ownership.write_text(header + "1" + suffix + "2" + suffix + "3" + suffix)
        args.ownership_quadrature.write_text("attribute,measure\n2,2\n3,1\n")
        build = json.loads(args.canonical_build_record.read_text())
        expected = meshio.read(args.source).points @ np.asarray(matrix)[:3, :3].T + np.asarray(matrix)[:3, 3]
        receipt = {"CanonicalBuildId": build["CanonicalBuildId"],
                   "CanonicalBuildSHA256": build["CanonicalBuildSHA256"],
                   "OutputMeshSHA256": digest(args.output),
                   "TransformSHA256": transform_hash,
                   "MaximumCoordinateError": float(np.max(np.linalg.norm(mesh.points-expected, axis=1)))}
        args.receipt.write_text(json.dumps(receipt) + "\n")


if __name__ == "__main__":
    main()
