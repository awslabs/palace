#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Deterministic stage executable for bounded mesh-DAG regression tests."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import meshio
import numpy as np

from mesh_stage_contract import footprint_provenance, footprint_segments
from semantic_mesh_contract import cut_surface_attributes, material_interface_attributes
from transform_coupon_source_contract import (read_transform, transform_semantic_contract,
                                                transformed_supports)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def rewrite(source, output, binary):
    meshio.write(output, meshio.read(source), file_format="gmsh22", binary=binary)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("metric", "adapt", "restore", "publish", "rigid",
                                          "ownership"))
    parser.add_argument("source", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--metric", type=Path)
    parser.add_argument("--mmg-seed", type=Path)
    parser.add_argument("--pins", type=Path)
    parser.add_argument("--recipe", type=Path)
    parser.add_argument("--fixed-triangles", type=Path)
    parser.add_argument("--required-tetrahedra", type=Path)
    parser.add_argument("--ownership", type=Path)
    parser.add_argument("--source-local-output", type=Path)
    parser.add_argument("--semantic-contract", type=Path)
    # Canonical supports input for "metric"; transformed supports output for "rigid".
    parser.add_argument("--transformed-supports", type=Path)
    parser.add_argument("--ownership-quadrature", type=Path)
    parser.add_argument("--transform", type=Path)
    parser.add_argument("--semantic-input", type=Path)
    parser.add_argument("--signature", type=Path)
    parser.add_argument("--boundary", type=Path)
    parser.add_argument("--mask", type=Path)
    parser.add_argument("--process", type=Path)
    parser.add_argument("--canonical-build-record", type=Path)
    parser.add_argument("--transformed-semantic", type=Path)
    parser.add_argument("--receipt", type=Path)
    parser.add_argument("--ownership-runtime", type=Path)
    parser.add_argument("--ownership-auditor", type=Path)
    # Metric-stage corner isotropy prescription recorded in the recipe.
    parser.add_argument("--normal", type=float)
    parser.add_argument("--tangent", type=float)
    # Bound seed census whose simplified footprint edges the recipe records.
    parser.add_argument("--seed-census", type=Path)
    # The trace basis (and its dimensionless size ratio) only when the case binds one.
    parser.add_argument("--trace-basis-contract", type=Path)
    parser.add_argument("--trace-vertices", type=Path)
    parser.add_argument("--trace-triangles", type=Path)
    parser.add_argument("--process-library", type=Path)
    parser.add_argument("--trace-basis-size-ratio", type=float)
    # Label-restoration gates (the seed stage carries the same values).
    parser.add_argument("--maximum-corner-aspect", type=float)
    parser.add_argument("--minimum-scaled-jacobian", type=float)
    parser.add_argument("--maximum-jacobian-condition", type=float)
    parser.add_argument("--maximum-quality-displacement-over-normal", type=float)
    args = parser.parse_args()
    if args.stage == "metric":
        if (args.mmg_seed is None or args.pins is None or args.recipe is None or
                args.fixed_triangles is None or args.required_tetrahedra is None):
            parser.error("metric requires MMG seed, pins, fixed triangles, required "
                         "tetrahedra, and recipe")
        args.output.write_text(json.dumps({"SeedSHA256": digest(args.source),
                                          "Metric": [1.0, 0.0, 1.0]}) + "\n")
        rewrite(args.source, args.mmg_seed, False)
        args.pins.write_text("1\n")
        args.fixed_triangles.write_text("1\n")
        # The fixture seed's first tetrahedron stands for the corner balls.
        args.required_tetrahedra.write_text("1\n")
        if args.semantic_contract is None or args.transformed_supports is None:
            parser.error("metric requires transformed semantic/support inputs")
        if args.normal is None or args.tangent is None:
            parser.error("metric requires --normal and --tangent")
        if args.seed_census is None or not args.seed_census.is_file():
            parser.error("metric requires the bound --seed-census")
        semantic = json.loads(args.semantic_contract.read_text())
        census = json.loads(args.seed_census.read_text())
        basis_paths = (args.trace_basis_contract, args.trace_vertices, args.trace_triangles,
                       args.process_library)
        trace_record = None
        if any(path is not None for path in basis_paths):
            if any(path is None or not path.is_file() for path in basis_paths) or \
                    args.trace_basis_size_ratio is None:
                parser.error("metric needs all four bound trace basis files and the size ratio")
            # The recipe records the census's rule with the cut-surface size statistics.
            trace_record = {**census["TraceBasisSizing"],
                            "CutSurfaceSize": {"CutTriangles": 7, "Minimum": .001,
                                               "Median": .08, "Maximum": .08}}
            triangle = census["TraceBasisSizing"]["MeshFrameTriangles"][0]
            trace_edges = {"InputSHA256": census["TraceBasisSizing"]["InputSHA256"], "Count": 3,
                           "Segments": [[*triangle[a], *triangle[b]] for a, b in ((0, 1), (1, 2), (2, 0))]}
        elif args.trace_basis_size_ratio is not None:
            parser.error("a trace basis size ratio needs the bound trace basis")
        args.recipe.write_text(json.dumps({
            "SeedSHA256": digest(args.source),
            "SemanticContract": semantic,
            "FootprintSegments": {"Provenance": footprint_provenance(census),
                                  "EtchBoundary": census["EtchBoundary"],
                                  "Tolerance": census["FootprintCollinearTolerance"],
                                  "Polygons": len(census["FootprintPolygons"]),
                                  "Segments": footprint_segments(census)},
            # The fixture seed's junction lines are the census's CAD junction curves.
            "JunctionSegments": {
                "Segments": census["JunctionCurves"]["Segments"],
                "Count": len(census["JunctionCurves"]["Segments"]),
                "TotalLength": sum(math.dist(segment[:3], segment[3:])
                                   for segment in census["JunctionCurves"]["Segments"]),
                "CutSurfaceAttributes": sorted(cut_surface_attributes(semantic)),
                "MaterialInterfaceAttributes": sorted(material_interface_attributes(semantic))},
            "NormalSize": args.normal, "TangentialSize": args.tangent,
            "CornerIsotropyRadius": args.tangent,
            "TruePhysicalCorners": semantic["SemanticCorners"],
            "FixedSurfaceTriangles": 1,
            "ProtectedCornerBalls": {"Radius": args.tangent, "FrozenTriangles": 1,
                "PerCorner": [{"Point": corner, "FrozenTriangles": 1}
                              for corner in semantic["SemanticCorners"]]},
            "Tetrahedra": len(meshio.read(args.source).get_cells_type("tetra")),
            "RequiredTetrahedra": {"Count": 1, "CornerRadius": args.tangent,
                "PerCorner": [{"Point": corner, "Tetrahedra": 1}
                              for corner in semantic["SemanticCorners"]],
                "LayerRequiredReach": None, "PerSpan": [], "IndexBase": 1},
            "TransformedSupportsArtifact": str(args.transformed_supports.resolve()),
            "TransformedSupportsSHA256": digest(args.transformed_supports),
            "TransformedSupports": json.loads(args.transformed_supports.read_text()),
            "FarSize": 1.0,
            "FarFieldBudgetPolicy": {"Name": "seed-fraction-far-field-v1",
                "RequestedFarSize": 1.0, "Pressure": 1.0,
                "EffectiveFarSize": 1.0},
            **({"TraceBasisSizing": trace_record, "TraceBasisEdges": trace_edges}
               if trace_record is not None else {})}) + "\n")
    elif args.stage == "adapt":
        if (args.metric is None or not args.metric.is_file() or
                args.pins is None or not args.pins.is_file() or
                args.fixed_triangles is None or not args.fixed_triangles.is_file() or
                args.required_tetrahedra is None or not args.required_tetrahedra.is_file()):
            parser.error("adapt requires metric, pins, fixed triangles, and required "
                         "tetrahedra")
        rewrite(args.source, args.output, True)
    elif args.stage == "restore":
        if (args.recipe is None or not args.recipe.is_file() or
                args.source_local_output is None):
            parser.error("restore requires --recipe and --source-local-output")
        rewrite(args.source, args.source_local_output, False)
        rewrite(args.source, args.output, False)
        # The restorer's own report next to its output: it found the recipe's
        # required tetrahedra (the fixture adapter keeps the one listed cell).
        recipe = json.loads(args.recipe.read_text())
        args.output.with_suffix(".projection.json").write_text(json.dumps(
            {"RequiredTetrahedra": recipe["RequiredTetrahedra"]["Count"],
             "RequiredVertices": 4, "MinimumScaledJacobianAfter": 0.5,
             "RequiredMaximumJacobianCondition": 1.0}) + "\n")
    elif args.stage == "publish":
        if args.ownership is None or args.ownership_quadrature is None:
            parser.error("publish requires ownership and quadrature ownership")
        if any(value is None or not value.is_file()
               for value in (args.process, args.signature, args.boundary)):
            parser.error("publish requires explicit process, signature, and boundary")
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
    elif args.stage == "ownership":
        if any(value is None or not value.is_file()
               for value in (args.process, args.signature, args.boundary)):
            parser.error("ownership requires explicit process, signature, and boundary")
        quadrature = Path(str(args.output) + ".quadrature.csv")
        if args.output.exists() or quadrature.exists():
            parser.error("ownership outputs must be fresh")
        header = ("attribute,elements,ambiguous_fraction,unresolved_elements,"
                  "unresolved_fraction,quadrature_rule,quadrature_order,quadrature_points,"
                  "quadrature_whole_measure,quadrature_owned_measure,"
                  "quadrature_relative_closure,quadrature_closure_tolerance,"
                  "quadrature_unmatched,quadrature_overlaps,quadrature_positive_weights\n")
        suffix = ",1,0,0,0,Gauss4,4,11,3,3,0,1e-12,0,0,1\n"
        args.output.write_text(header + "1" + suffix + "2" + suffix + "3" + suffix)
        quadrature.write_text("attribute,measure\n2,2\n3,1\n")
    else:
        required = (args.transform, args.semantic_input, args.signature, args.boundary,
                    args.mask, args.process, args.canonical_build_record,
                    args.transformed_semantic, args.transformed_supports, args.receipt,
                    args.ownership, args.ownership_quadrature, args.ownership_runtime,
                    args.ownership_auditor)
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
        args.transformed_semantic.write_text(json.dumps(semantic) + "\n")
        supports = transformed_supports(
            args.signature.parent, matrix, signature=args.signature,
            boundary=args.boundary, mask=args.mask, transform_sha256=transform_hash,
            semantic_sha256=semantic_hash)
        args.transformed_supports.write_text(json.dumps(supports) + "\n")
        transform_csv = Path(str(args.receipt) + ".transform.csv")
        transform_csv.write_text(",".join(format(value, ".17g")
                                          for row in matrix for value in row) + "\n")
        command = [str(args.ownership_runtime.resolve()), str(args.ownership_auditor.resolve()),
                   "ownership", str(args.output.resolve()), str(args.ownership.resolve()),
                   "--transform", str(transform_csv.resolve()),
                   "--process", str(args.process.resolve()),
                   "--signature", str(args.signature.resolve()),
                   "--boundary", str(args.boundary.resolve())]
        subprocess.run(command, check=True)
        Path(str(args.ownership) + ".quadrature.csv").replace(args.ownership_quadrature)
        build = json.loads(args.canonical_build_record.read_text())
        expected = meshio.read(args.source).points @ np.asarray(matrix)[:3, :3].T + np.asarray(matrix)[:3, 3]
        receipt = {"CanonicalBuildId": build["CanonicalBuildId"],
                   "CanonicalBuildSHA256": build["CanonicalBuildSHA256"],
                   "OutputMeshSHA256": digest(args.output),
                   "TransformSHA256": transform_hash,
                   "SourceInputSHA256": {
                       "source-semantic-contract": semantic_hash,
                       "source-signature": digest(args.signature),
                       "source-boundary": digest(args.boundary),
                       "source-mask": digest(args.mask),
                       "source-process": digest(args.process)},
                   "MaximumCoordinateError": float(np.max(np.linalg.norm(mesh.points-expected, axis=1))),
                   "TransformedSemanticSHA256": digest(args.transformed_semantic),
                   "TransformedSupportsSHA256": digest(args.transformed_supports),
                   "OwnershipSHA256": digest(args.ownership),
                   "OwnershipQuadratureSHA256": digest(args.ownership_quadrature),
                   "OwnershipCommand": command,
                   "OwnershipRuntimeSHA256": digest(args.ownership_runtime),
                   "OwnershipAuditorSHA256": digest(args.ownership_auditor)}
        args.receipt.write_text(json.dumps(receipt) + "\n")

if __name__ == "__main__":
    main()
