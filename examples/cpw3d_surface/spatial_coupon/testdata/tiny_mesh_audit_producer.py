#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Write a real two-material Gmsh mesh for evidence-chain regression tests.

With `--prism-tubes true` (the Gmsh-only build of supervisor decision 38) the mesh
also carries one prism with a pyramid on a lateral face and two quadrangle boundary
faces, and the census records the prism tube design, the corner grading and the
per-type quality in the production mesher's schema."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import meshio
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from mixed_mesh import volume_quality  # noqa: E402


def tube_cells(scale, offset):
    """One right prism (bottom a b c, top a' b' c') with a pyramid on its lateral face
    (a b b' a'): six points, the prism, the pyramid, the boundary triangles and the two
    remaining lateral quadrangles (0-based indices from `offset`)."""
    base = 20.0
    points = scale * np.array([[base, 0., 0.], [base + 1., 0., 0.], [base, 1., 0.],
                               [base, 0., 1.], [base + 1., 0., 1.], [base, 1., 1.],
                               [base + .5, -.5, .5]])
    a, b, c, a1, b1, c1, apex = range(offset, offset + 7)
    wedge = [[a, b, c, a1, b1, c1]]
    pyramid = [[a, b, b1, a1, apex]]
    triangles = [[a, c, b], [a1, b1, c1], [a, b, apex], [b, b1, apex], [b1, a1, apex], [a1, a, apex]]
    quads = [[b, c, c1, b1], [c, a, a1, c1]]
    return points, wedge, pyramid, triangles, quads


def produce(output, transform, scale=1.0, prism_tubes=False):
    matrix = np.asarray(transform, dtype=float).reshape(4, 4)
    points = scale * np.array([[0., 0., 0.], [.08, 0., 0.], [0., .001, 0.],
                               [0., 0., .001], [0., 0., -.001], [10., 0., 0.],
                               [10.001, 0., 0.], [10., .001, 0.], [10., 0., .001]])
    tetrahedra = [[0, 1, 2, 3], [0, 2, 1, 4], [5, 6, 7, 8]]
    triangles = [[0, 1, 3], [0, 2, 3], [1, 2, 3],
                 [0, 2, 4], [0, 1, 4], [2, 1, 4], [0, 1, 2],
                 [5, 6, 7], [5, 8, 6], [5, 7, 8], [6, 8, 7]]
    extra_points = []
    for i in range(1, 6):
        for x in (i, i + 1):
            offset = len(points) + len(extra_points)
            extra_points.extend(([x, i, 0], [x + .08, i, 0],
                                 [x, i + .001, 0], [x, i, .001]))
            tetrahedra.append([offset, offset + 1, offset + 2, offset + 3])
            triangles.extend(([offset, offset + 1, offset + 2],
                              [offset, offset + 3, offset + 1],
                              [offset, offset + 2, offset + 3],
                              [offset + 1, offset + 3, offset + 2]))
    points = np.vstack((points, scale * np.asarray(extra_points)))
    tube = tube_cells(scale, len(points)) if prism_tubes else None
    if tube is not None:
        tube_points, wedge, pyramid, tube_triangles, quads = tube
        points = np.vstack((points, tube_points))
        triangles.extend(tube_triangles)
    homogeneous = np.column_stack((points, np.ones(len(points))))
    points = (homogeneous @ matrix.T)[:, :3]
    tetrahedra = np.asarray(tetrahedra)
    triangles = np.asarray(triangles)
    cells = [("triangle", triangles), ("tetra", tetrahedra)]
    physical = [np.array([1, 1, 1, 2, 2, 2, 3] + [1] * (len(triangles) - 7)),
                np.array([1, 7] + [1] * (len(tetrahedra) - 2))]
    if tube is not None:
        cells += [("quad", np.asarray(quads)), ("wedge", np.asarray(wedge)),
                  ("pyramid", np.asarray(pyramid))]
        physical += [np.ones(len(quads), dtype=int), np.ones(len(wedge), dtype=int),
                     np.ones(len(pyramid), dtype=int)]
    mesh = meshio.Mesh(points, cells,
                       cell_data={"gmsh:physical": physical,
                                  "gmsh:geometrical": [np.ones(len(block), dtype=int)
                                                       for _, block in cells]},
                       field_data={"substrate": np.array([1, 3]),
                                   "vacuum": np.array([7, 3]),
                                   "surface_1": np.array([1, 2]),
                                   "surface_2": np.array([2, 2]),
                                   "surface_3": np.array([3, 2])})
    meshio.write(output, mesh, file_format="gmsh22", binary=False)
    return mesh


def corner_grading(corner_size, ratio, normal, radius):
    """The production seeder's CornerGrading record (shell radii and sizes)."""
    radii, size = [], corner_size
    while size < normal:
        radii.append(radii[-1] + size if radii else size)
        size *= ratio
    return {"CornerSize": corner_size, "GrowthRatio": ratio, "NormalSize": normal, "Radius": radius,
            "Reach": (normal - corner_size) / (ratio - 1.0), "Rule": "fixture",
            "ShellRadii": radii + [radius],
            "ShellSizes": [corner_size * ratio**k for k in range(len(radii))] + [normal]}


def prism_tube_record(mesh, tubes):
    """The production seeder's PrismTubes census record for the fixture: one tube of
    one layer, three rings, the per-type quality measured on the written mesh."""
    quality = volume_quality(mesh)
    per_type = {name: {"Count": record["Samples"], "PositiveOrientation": record["PositiveOrientation"],
                       "NonpositiveCells": record["NonpositiveCells"],
                       "MinimumScaledJacobian": record["MinimumScaledJacobian"],
                       "MaximumJacobianCondition": record["MaximumJacobianCondition"]}
                for name, record in quality["ByType"].items()}
    per_type["Total"] = quality["Samples"]
    inner, ratio = tubes["edge_size"], tubes["ratio"]
    rings = 3
    return {"Rule": "fixture", "InnerSize": inner, "GrowthRatio": ratio,
            "TangentialSize": tubes["lc_tangent"], "NormalSize": tubes["lc_fine"],
            "FarSize": tubes["lc_far"], "FarGrowth": tubes["far_growth"],
            "Tubes": [{"Spacing": tubes["lc_tangent"], "Layers": 1, "Length": tubes["lc_tangent"],
                       "Edge": "top", "Conductor": 1}],
            "TubeCount": 1, "TotalTubeLength": tubes["lc_tangent"], "Layers": 1,
            "SpacingMinimum": tubes["lc_tangent"], "SpacingMaximum": tubes["lc_tangent"],
            "InnermostArc": inner * math.pi / 6, "MaximumPrismEdgeAspect": tubes["lc_tangent"] / inner,
            "Section": {"Rings": rings, "RingSizes": [inner * ratio**k for k in range(rings)],
                        "PyramidHeight": .5 * inner * ratio**(rings - 1),
                        "PyramidHeightOverOuterRing": .5, "SectorDegrees": 30.0},
            "Prisms": quality["ByType"].get("Prism", {}).get("Samples", 0),
            "Pyramids": quality["ByType"].get("Pyramid", {}).get("Samples", 0),
            "SizeLaws": {"NormalSize": tubes["lc_fine"], "FarSize": tubes["lc_far"],
                         "FarGrowth": tubes["far_growth"], "RadialGrowth": 1.0,
                         "ProtectedDistance": 2.0 * tubes["lc_fine"], "TubeRule": "fixture tube rule",
                         "BandRule": "fixture band rule", "Composition": "fixture composition",
                         "CornerExteriorGrowth": tubes["far_growth"],
                         "TraceBasisVolumeGrowth": tubes["far_growth"],
                         "CornerIsotropyRadius": tubes["radius"],
                         "JunctionVolumeRule": "fixture", "CornerExteriorRule": "fixture",
                         "TraceBasisVolumeRule": "fixture",
                         "Achieved": {"CornerExterior": {"Shells": [{"Cells": 1}]},
                                      "JunctionLines": {"Shells": [{"Cells": 1}]},
                                      "TraceApexes": ({"Shells": [{"Cells": 1}]}
                                                      if tubes["trace_basis"] else None)}},
            "Quality": per_type,
            "CapRegions": {"Caps": 1, "MinimumScaledJacobian": .5, "MaximumJacobianCondition": 2.0},
            "CutSurface": {"All": {"Elements": 7}},
            "Bands": {"Prescribed": tubes["lc_fine"],
                      "JunctionFirstLayer": {
                          "CutSurface": {"Elements": 2, "Prescribed": tubes["lc_fine"],
                                         "TransverseP50": tubes["lc_fine"],
                                         "TransverseP90": tubes["lc_fine"],
                                         "AchievedOverPrescribedP50": 1.0},
                          "Tetrahedra": {"Elements": 1, "Prescribed": tubes["lc_fine"],
                                         "TransverseP50": tubes["lc_fine"],
                                         "TransverseP90": tubes["lc_fine"],
                                         "AchievedOverPrescribedP50": 1.0}}},
            "BandCurves": {"Count": 0, "TotalLength": 0.0, "Spacing": tubes["lc_fine"],
                           "Segments": [], "Rule": "fixture"},
            "FarFieldBudgetPolicy": {"Name": "gmsh-only-fail-closed-cap", "Pressure": 1.0,
                                     "RequestedFarSize": tubes["lc_far"],
                                     "EffectiveFarSize": tubes["lc_far"],
                                     "Elements": quality["Samples"],
                                     "MaximumElements": tubes["max_elements"]}}


def junction_curves(scale):
    """The fixture's cut-surface/material-interface junction lines: the three edges of
    the interface triangle [0, 1, 2] shared with the cut triangles of the first tet."""
    a, b, c = [0.0, 0.0, 0.0], [.08 * scale, 0.0, 0.0], [0.0, .001 * scale, 0.0]
    segments = [[*a, *b], [*b, *c], [*c, *a]]
    return {"Count": len(segments),
            "TotalLength": sum(math.dist(segment[:3], segment[3:]) for segment in segments),
            "Segments": segments, "CurvedCurves": 0, "Rule": "fixture"}


def trace_basis_sizing(basis_paths, ratio, scale, slope=1.0):
    """Fixture record of the trace-basis cut-surface size rule (schema of the production
    seeder): the bound input digests, the dimensionless ratio and one mesh-frame basis
    triangle on the fixture's cut face."""
    if basis_paths is None:
        return None
    digests = {role: hashlib.sha256(path.read_bytes()).hexdigest()
               for role, path in zip(("BasisContract", "TraceVertices", "TraceTriangles",
                                      "ProcessLibrary"), basis_paths)}
    return {"Ratio": ratio, "RatioIsDimensionless": True, "Rule": "fixture",
            "InputSHA256": digests, "Lower": [0.0, 0.0, -.001 * scale],
            "Upper": [10.001 * scale, .001 * scale, .001 * scale],
            "Triangles": 1, "BasisEdgesBelowFarSize": 1, "MinimumRequestedSize": ratio * .001 * scale,
            "MeshSizeMinimum": min(ratio * .001 * scale, .1), "GradingSlope": slope,
            "MeshFrameTriangles": [[[0.0, 0.0, 0.0], [.08 * scale, 0.0, 0.0], [0.0, 0.0, .001 * scale]]]}


def corner_census(output, contract_path, radius, isotropic_size, etch_boundary=None,
                  scale=1.0, basis_paths=None, ratio=None, gates=None, tubes=None, mesh=None):
    """Recorded corner-ball census of the fixture seed (schema of the production seeder).

    `gates` is (MaximumCornerAspect, MinimumScaledJacobian, MaximumJacobianCondition,
    DisplacementBoundOverNormal): the seed-side required-region optimization record
    (decision 30; the condition gate, decision 34) is written when the production
    seeder's four gate options are passed."""
    contract = json.loads(contract_path.read_text())
    corners = contract["SemanticCorners"]
    quality = None
    if gates is not None:
        maximum_aspect, minimum_scaled, maximum_condition, ratio_bound = gates
        quality = {"MaximumCornerAspect": maximum_aspect, "CornerAspectTarget": .95 * maximum_aspect,
                   "MinimumScaledJacobian": minimum_scaled, "ScaledJacobianTarget": 2 * minimum_scaled,
                   "MaximumJacobianCondition": maximum_condition,
                   "RequiredMaximumJacobianConditionBefore": 1.0,
                   "RequiredMaximumJacobianConditionAfter": 1.0,
                   "RequiredCellsAboveConditionAfter": 0,
                   "DisplacementBoundOverNormal": ratio_bound, "LayerRequiredReach": None,
                   "RequiredTetrahedra": 1, "CornerAspectsBefore": [1.0 for _ in corners],
                   "CornerAspectsAfter": [1.0 for _ in corners], "CornerMoves": [0 for _ in corners],
                   "RequiredMinimumScaledJacobianBefore": 1.0,
                   "RequiredMinimumScaledJacobianAfter": 1.0,
                   "RequiredCellsBelowTargetBefore": 0, "RequiredCellsBelowTargetAfter": 0,
                   "RequiredCellsBelowGateAfter": 0, "RepairComponents": 0, "RepairMoves": 0,
                   "MovedVertices": 0, "MaximumDisplacement": 0.0,
                   "MaximumDisplacementOverBound": 0.0}
    tube_records = {}
    if tubes is not None:
        tube_records = {"PrismTubes": prism_tube_record(mesh, tubes),
                        "CornerGrading": corner_grading(tubes["corner_size"], tubes["ratio"],
                                                        isotropic_size, radius),
                        "EdgeLayer": None}
        if quality is not None:
            quality["CornerSize"] = tubes["corner_size"]
    output.write_text(json.dumps({
        "Version": 1, "Frame": "SourceLocal", "SemanticCorners": corners,
        "SemanticContract": str(contract_path),
        "SemanticContractSHA256": hashlib.sha256(contract_path.read_bytes()).hexdigest(),
        "CornerIsotropyRadius": radius, "IsotropicSize": isotropic_size, **tube_records,
        "EtchBoundary": "producer-default" if etch_boundary is None else str(etch_boundary),
        "EtchBoundarySHA256": None if etch_boundary is None else
                              hashlib.sha256(etch_boundary.read_bytes()).hexdigest(),
        # Simplified etch footprint polygons (schema of the production seeder): the
        # fixture footprint is one square whose duplicated vertex was merged.
        "FootprintCollinearTolerance": 1e-6,
        "FootprintSimplification": {"Rule": "fixture", "Polygons": 1, "RemovedVertices": 1,
                                    "MaximumRelativeDeviation": 0.0},
        "FootprintPolygons": [{"Conductor": 1, "Plane": 0.0, "Hole": False,
                               "Points": [[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]],
                               "Simplification": {"OriginalVertices": 5, "Vertices": 4,
                                                  "RemovedVertexCount": 1,
                                                  "RemovedVertexIndices": [3],
                                                  "MaximumDeviation": 0.0,
                                                  "MaximumDeviationLocalScale": 10.0,
                                                  "MaximumRelativeDeviation": 0.0,
                                                  "Tolerance": 1e-6}}],
        "JunctionCurves": junction_curves(scale),
        "TraceBasisSizing": trace_basis_sizing(basis_paths, ratio, scale,
                                               1.0 if tubes is None else tubes["far_growth"]),
        "SeedQualityOptimization": quality,
        "InterfaceAreaUnits": "um^2",
        # One row per contract boundary label, as written in the fixture seed.
        "InterfaceAreas": [{"Attribute": 1, "Name": "surface_1", "Triangles": 7,
                            "Area": 2.0},
                           {"Attribute": 2, "Name": "surface_2", "Triangles": 3, "Area": 2.0},
                           {"Attribute": 3, "Name": "surface_3", "Triangles": 1, "Area": 1.0}],
        "Corners": [{"Corner": index, "Point": corner, "BallEdges": 0}
                    for index, corner in enumerate(corners)],
        "LongitudinalFaces": [{"Surface": 1, "Triangles": 2, "InteriorNodes": 0,
                               "InteriorNodesAwayFromCorners": 0, "FullHeightTriangles": 0,
                               "FullHeightTrianglesAwayFromCorners": 0}]}, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path); parser.add_argument("transform", type=Path)
    parser.add_argument("--scale", type=float, default=1.0)
    # Bound source inputs: the fixture only requires them to exist, mirroring the
    # production seeder's positional signature and --mask/--boundary consumption.
    parser.add_argument("--signature", type=Path, required=True)
    parser.add_argument("--mask", type=Path, required=True)
    parser.add_argument("--boundary", type=Path, required=True)
    # Corner isotropy contract of the production seeder: contract corners, ball
    # radius, isotropic size and the recorded census artifact.
    parser.add_argument("--semantic-contract", type=Path, required=True)
    parser.add_argument("--corner-isotropy-radius", type=float, required=True)
    parser.add_argument("--lc-fine", type=float, required=True)
    parser.add_argument("--corner-census", type=Path, required=True)
    # The device etch footprint is consumed only when the case binds one.
    parser.add_argument("--etch-boundary", type=Path)
    # The trace basis (and its dimensionless size ratio) only when the case binds one.
    parser.add_argument("--trace-basis-contract", type=Path)
    parser.add_argument("--trace-vertices", type=Path)
    parser.add_argument("--trace-triangles", type=Path)
    parser.add_argument("--process-library", type=Path)
    parser.add_argument("--trace-basis-size-ratio", type=float)
    # The seed-side required-region gates (decisions 30 and 34), required together.
    parser.add_argument("--maximum-corner-aspect", type=float)
    parser.add_argument("--minimum-scaled-jacobian", type=float)
    parser.add_argument("--maximum-jacobian-condition", type=float)
    parser.add_argument("--maximum-quality-displacement-over-normal", type=float)
    # The Gmsh-only build (decision 38): the tube recipe options of the production mesher.
    parser.add_argument("--prism-tubes", choices=("true", "false"), default="false")
    parser.add_argument("--lc-tangent", type=float)
    parser.add_argument("--lc-far", type=float)
    parser.add_argument("--edge-size", type=float)
    parser.add_argument("--edge-growth-ratio", type=float, default=2.0)
    parser.add_argument("--corner-size", type=float)
    parser.add_argument("--far-growth", type=float)
    parser.add_argument("--max-elements", type=int)
    args = parser.parse_args()
    tubes = None
    if args.prism_tubes == "true":
        tubes = {"lc_tangent": args.lc_tangent, "lc_far": args.lc_far, "lc_fine": args.lc_fine,
                 "edge_size": args.edge_size, "ratio": args.edge_growth_ratio,
                 "corner_size": args.corner_size, "far_growth": args.far_growth,
                 "max_elements": args.max_elements}
        if any(value is None for value in tubes.values()) or args.edge_size != args.corner_size:
            parser.error("prism tubes need --lc-tangent, --lc-far, --edge-size == --corner-size, "
                         "--far-growth and --max-elements")
        tubes["radius"] = args.corner_isotropy_radius
        tubes["trace_basis"] = args.trace_basis_contract is not None
    gates = (args.maximum_corner_aspect, args.minimum_scaled_jacobian,
             args.maximum_jacobian_condition, args.maximum_quality_displacement_over_normal)
    if any(value is not None for value in gates):
        if any(value is None or value <= 0 for value in gates):
            parser.error("the four seed quality gate options are required together")
    else:
        gates = None
    if args.etch_boundary is not None and not args.etch_boundary.is_file():
        parser.error("the bound retained etch footprint must exist")
    basis_paths = (args.trace_basis_contract, args.trace_vertices, args.trace_triangles,
                   args.process_library)
    if any(path is not None for path in basis_paths):
        if any(path is None or not path.is_file() for path in basis_paths) or \
                args.trace_basis_size_ratio is None:
            parser.error("the trace basis needs all four bound files and its size ratio")
    elif args.trace_basis_size_ratio is not None:
        parser.error("a trace basis size ratio needs the bound trace basis")
    else:
        basis_paths = None
    if any(not path.is_file() for path in (args.signature, args.mask, args.boundary,
                                           args.semantic_contract)):
        parser.error("signature, mask, boundary, and semantic contract inputs must exist")
    if args.corner_census.exists():
        parser.error("corner census output must be fresh")
    mesh = produce(args.output, json.loads(args.transform.read_text()), args.scale, tubes is not None)
    corner_census(args.corner_census, args.semantic_contract, args.corner_isotropy_radius,
                  args.lc_fine, args.etch_boundary, args.scale, basis_paths,
                  args.trace_basis_size_ratio, gates, tubes, mesh)


if __name__ == "__main__":
    main()
