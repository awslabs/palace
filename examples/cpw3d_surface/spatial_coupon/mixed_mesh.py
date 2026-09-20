#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Mixed linear meshes (tetrahedra / prisms / pyramids; triangles / quadrangles) of
the Gmsh-only coupon build for the Python audits.

`simplicial_view` splits every prism into three and every pyramid into two
tetrahedra and every quadrangle into two triangles with the minimum-vertex
diagonal rule (the diagonal of a quadrangle face leaves its vertex of smallest
global index), which makes the split conforming across shared faces, so the
topology, adjacency, measure and protected-surface audits written for
tetrahedral meshes apply unchanged (planar-faced cells keep their volumes and
areas exactly).  Element quality is never measured on the view: `volume_quality`
evaluates the corner-frame Jacobians of the native cells (the same frames as the
mesher's census), and `element_counts` / `h1_dofs` count the native cells."""
import meshio
import numpy as np

VOLUME_KINDS = ("tetra", "wedge", "pyramid")
SURFACE_KINDS = ("triangle", "quad")
GMSH_NAMES = {"tetra": "Tetrahedron", "wedge": "Prism", "pyramid": "Pyramid",
              "triangle": "Triangle", "quad": "Quadrangle"}
VERTEX_COUNT = {"tetra": 4, "wedge": 6, "pyramid": 5, "triangle": 3, "quad": 4}
# Corner frames (vertex, three edge ends) whose edges span the local Jacobian;
# the tetrahedron keeps the production vertex-0 frame (0-based meshio order).
CORNER_FRAMES = {
    "tetra": [(0, 1, 2, 3)],
    "wedge": [(0, 1, 2, 3), (1, 2, 0, 4), (2, 0, 1, 5), (3, 5, 4, 0), (4, 3, 5, 1), (5, 4, 3, 2)],
    "pyramid": [(0, 1, 3, 4), (1, 2, 0, 4), (2, 3, 1, 4), (3, 0, 2, 4)],
}
QUALITY_QUANTILES = (0.0, 0.01, 0.05, 0.5, 0.95, 0.99, 1.0)


def _label_key(mesh):
    return "gmsh:physical" if "gmsh:physical" in mesh.cell_data else "medit:ref"


def _base_kind(cell_type):
    for kind in (*VOLUME_KINDS, *SURFACE_KINDS):
        if cell_type == kind or cell_type.startswith(kind):
            return kind
    return None


def cell_blocks(mesh, kinds):
    """[(kind, vertex connectivity, labels)] of every cell block whose kind is in
    `kinds`; only the vertex nodes of a higher-order block are kept."""
    key = _label_key(mesh)
    result = []
    for cell, labels in zip(mesh.cells, mesh.cell_data[key]):
        kind = _base_kind(cell.type)
        if kind in kinds:
            result.append((kind, np.asarray(cell.data)[:, :VERTEX_COUNT[kind]],
                           np.asarray(labels)))
    return result


def element_counts(mesh):
    """Native cell counts per Gmsh type name plus the volume and surface totals."""
    counts = {GMSH_NAMES[kind]: 0 for kind in (*VOLUME_KINDS, *SURFACE_KINDS)}
    for kind, connectivity, _ in cell_blocks(mesh, (*VOLUME_KINDS, *SURFACE_KINDS)):
        counts[GMSH_NAMES[kind]] += int(len(connectivity))
    counts["VolumeElements"] = sum(counts[GMSH_NAMES[kind]] for kind in VOLUME_KINDS)
    counts["SurfaceElements"] = sum(counts[GMSH_NAMES[kind]] for kind in SURFACE_KINDS)
    return counts


def is_mixed(mesh):
    counts = element_counts(mesh)
    return any(counts[GMSH_NAMES[kind]] for kind in ("wedge", "pyramid", "quad"))


def _quad_diagonal_from_first(quad):
    """True when the minimum-vertex diagonal of quad (a, b, c, d) is a-c."""
    return min(quad[0], quad[2]) < min(quad[1], quad[3])


def split_quad(quad):
    a, b, c, d = (int(v) for v in quad)
    if _quad_diagonal_from_first(quad):
        return [(a, b, c), (a, c, d)]
    return [(a, b, d), (b, c, d)]


def split_pyramid(cell):
    a, b, c, d, apex = (int(v) for v in cell)
    if _quad_diagonal_from_first(cell[:4]):
        return [(a, b, c, apex), (a, c, d, apex)]
    return [(a, b, d, apex), (b, c, d, apex)]


def split_wedge(cell):
    """Three tetrahedra of a prism (bottom 0 1 2, top 3 4 5) whose quadrangle faces
    are cut by the minimum-vertex diagonal: rotate the prism so its smallest
    vertex is at position 0, then the two faces at that vertex are cut from it and
    the third face (1 2 5 4) by its own smallest vertex (Dompierre et al.)."""
    cell = [int(v) for v in cell]
    smallest = int(np.argmin(cell))
    if smallest >= 3:
        # Swap the triangles (the volume is re-oriented by the caller).
        cell = cell[3:] + cell[:3]
        smallest -= 3
    rotation = [(smallest + i) % 3 for i in range(3)]
    cell = [cell[i] for i in rotation] + [cell[3 + i] for i in rotation]
    v0, v1, v2, v3, v4, v5 = cell
    if min(v1, v5) < min(v2, v4):
        return [(v0, v1, v2, v5), (v0, v1, v5, v4), (v0, v4, v5, v3)]
    return [(v0, v1, v2, v4), (v0, v4, v2, v5), (v0, v4, v5, v3)]


def _oriented(points, tetrahedra):
    xyz = points[tetrahedra]
    signed = np.einsum("ij,ij->i", np.cross(xyz[:, 1] - xyz[:, 0], xyz[:, 2] - xyz[:, 0]),
                       xyz[:, 3] - xyz[:, 0])
    flipped = tetrahedra.copy()
    negative = signed < 0
    flipped[negative, 2], flipped[negative, 3] = tetrahedra[negative, 3], tetrahedra[negative, 2]
    return flipped


def simplicial_view(mesh):
    """A meshio mesh of one triangle block and one tetrahedron block carrying the
    physical labels of the source cells (a tetrahedral mesh is returned with the
    same content); prisms, pyramids and quadrangles are split conformingly.  The
    split tetrahedra are oriented positively; native orientation is judged by
    `volume_quality`."""
    key = _label_key(mesh)
    points = np.asarray(mesh.points, dtype=float)
    triangles, triangle_labels, tetrahedra, tetra_labels = [], [], [], []
    for kind, connectivity, labels in cell_blocks(mesh, (*VOLUME_KINDS, *SURFACE_KINDS)):
        if kind == "triangle":
            triangles.append(connectivity); triangle_labels.append(labels)
        elif kind == "tetra":
            tetrahedra.append(connectivity); tetra_labels.append(labels)
        elif kind == "quad":
            split = [split_quad(q) for q in connectivity]
            triangles.append(np.asarray(split, dtype=np.int64).reshape(-1, 3))
            triangle_labels.append(np.repeat(labels, 2))
        else:
            splitter = split_wedge if kind == "wedge" else split_pyramid
            per_cell = 3 if kind == "wedge" else 2
            split = [splitter(c) for c in connectivity]
            cells = np.asarray(split, dtype=np.int64).reshape(-1, 4)
            tetrahedra.append(_oriented(points, cells)); tetra_labels.append(np.repeat(labels, per_cell))
    if not tetrahedra or not triangles:
        raise ValueError("Mesh lacks volume or surface elements")
    return meshio.Mesh(points,
                       [("triangle", np.concatenate(triangles).astype(np.int64)),
                        ("tetra", np.concatenate(tetrahedra).astype(np.int64))],
                       cell_data={key: [np.concatenate(triangle_labels),
                                        np.concatenate(tetra_labels)]},
                       point_data=dict(mesh.point_data), field_data=dict(mesh.field_data))


def _frame_quality(points, connectivity, frames):
    """Per cell: minimum corner determinant, maximum corner condition, minimum corner
    scaled Jacobian (|det| / product of the three edge lengths)."""
    determinant = np.full(len(connectivity), np.inf)
    condition = np.zeros(len(connectivity))
    scaled = np.full(len(connectivity), np.inf)
    xyz = points[connectivity]
    for vertex, a, b, c in frames:
        p0 = xyz[:, vertex]
        jacobian = np.stack((xyz[:, a] - p0, xyz[:, b] - p0, xyz[:, c] - p0), axis=2)
        det = np.linalg.det(jacobian)
        singular = np.linalg.svd(jacobian, compute_uv=False)
        lengths = np.linalg.norm(jacobian, axis=1)
        determinant = np.minimum(determinant, det)
        with np.errstate(divide="ignore", invalid="ignore"):
            condition = np.maximum(condition, singular[:, 0] / singular[:, -1])
            scaled = np.minimum(scaled, np.abs(det) / np.prod(lengths, axis=1))
    return determinant, condition, scaled


def volume_quality(mesh):
    """Native per-type Jacobian statistics with the production top-level semantics:
    PositiveOrientation over every type, MaximumJacobianCondition over every type,
    MinimumScaledJacobian and its quantiles over the tetrahedra (the prisms and
    pyramids are judged by orientation and condition; their scaled Jacobians are
    reported per type), Samples the number of volume elements of every type."""
    points = np.asarray(mesh.points, dtype=float)
    by_type = {}
    for kind, connectivity, _ in cell_blocks(mesh, VOLUME_KINDS):
        if not len(connectivity):
            continue
        determinant, condition, scaled = _frame_quality(points, connectivity, CORNER_FRAMES[kind])
        if not np.all(np.isfinite(condition)) or not np.all(np.isfinite(scaled)):
            raise ValueError(f"Invalid {GMSH_NAMES[kind]} Jacobian audit")
        name = GMSH_NAMES[kind]
        record = by_type.setdefault(name, {"Samples": 0, "determinant": [], "condition": [],
                                           "scaled": []})
        record["Samples"] += int(len(connectivity))
        record["determinant"].append(determinant); record["condition"].append(condition)
        record["scaled"].append(scaled)
    if "Tetrahedron" not in by_type:
        raise ValueError("Mesh lacks tetrahedra")
    result_by_type = {}
    for name, record in by_type.items():
        determinant = np.concatenate(record["determinant"])
        condition = np.concatenate(record["condition"])
        scaled = np.concatenate(record["scaled"])
        result_by_type[name] = {
            "Samples": record["Samples"],
            "PositiveOrientation": bool(np.all(determinant > 0.0)),
            "NonpositiveCells": int(np.sum(determinant <= 0.0)),
            "MinimumScaledJacobian": float(scaled.min()),
            "MaximumJacobianCondition": float(condition.max()),
            "ScaledJacobianQuantiles": np.quantile(scaled, QUALITY_QUANTILES).tolist(),
            "JacobianConditionQuantiles": np.quantile(condition, QUALITY_QUANTILES).tolist()}
    tetrahedra = result_by_type["Tetrahedron"]
    return {"Samples": sum(item["Samples"] for item in result_by_type.values()),
            "PositiveOrientation": all(item["PositiveOrientation"] for item in result_by_type.values()),
            "MinimumScaledJacobian": tetrahedra["MinimumScaledJacobian"],
            "MaximumJacobianCondition": max(item["MaximumJacobianCondition"]
                                            for item in result_by_type.values()),
            "ScaledJacobianQuantiles": tetrahedra["ScaledJacobianQuantiles"],
            "JacobianConditionQuantiles": tetrahedra["JacobianConditionQuantiles"],
            "ByType": result_by_type,
            "ScaledJacobianGate": "Tetrahedron",
            "Rule": ("corner-frame Jacobians of the native linear cells: orientation and "
                     "condition judge every type, the scaled Jacobian judges the tetrahedra")}


H1_ENTITY_NAMES = ("Vertices", "Edges", "TriangleFaces", "QuadFaces", "Tetrahedra", "Prisms", "Pyramids")


def h1_entity_counts(mesh):
    """The mesh entities an H1 space is counted on: unique vertices, edges,
    triangular and quadrangular faces of the volume cells and the cells by type."""
    edge_sets = []
    vertices = set()
    face_edges = {"tetra": ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)),
                  "wedge": ((0, 1), (1, 2), (2, 0), (3, 4), (4, 5), (5, 3), (0, 3), (1, 4), (2, 5)),
                  "pyramid": ((0, 1), (1, 2), (2, 3), (3, 0), (0, 4), (1, 4), (2, 4), (3, 4))}
    faces = {"tetra": (((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)), ()),
             "wedge": (((0, 1, 2), (3, 4, 5)), ((0, 1, 4, 3), (1, 2, 5, 4), (2, 0, 3, 5))),
             "pyramid": (((0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)), ((0, 1, 2, 3),))}
    cells = {"tetra": 0, "wedge": 0, "pyramid": 0}
    triangle_faces, quad_faces = [], []
    for kind, connectivity, _ in cell_blocks(mesh, VOLUME_KINDS):
        if not len(connectivity):
            continue
        cells[kind] += len(connectivity)
        vertices.update(np.unique(connectivity).tolist())
        edge_sets.append(np.sort(connectivity[:, face_edges[kind]].reshape(-1, 2), axis=1))
        tri, quad = faces[kind]
        if tri:
            triangle_faces.append(np.sort(connectivity[:, tri].reshape(-1, 3), axis=1))
        if quad:
            quad_faces.append(np.sort(connectivity[:, quad].reshape(-1, 4), axis=1))
    edges = len(np.unique(np.concatenate(edge_sets), axis=0)) if edge_sets else 0
    triangles = len(np.unique(np.concatenate(triangle_faces), axis=0)) if triangle_faces else 0
    quads = len(np.unique(np.concatenate(quad_faces), axis=0)) if quad_faces else 0
    return {"Vertices": len(vertices), "Edges": int(edges), "TriangleFaces": int(triangles),
            "QuadFaces": int(quads), "Tetrahedra": cells["tetra"], "Prisms": cells["wedge"],
            "Pyramids": cells["pyramid"]}


def h1_dofs_from_counts(counts, order):
    """Palace's H1 count on the entity counts: vertices, (p - 1) per edge, (p - 1)(p - 2)/2
    per triangle, (p - 1)^2 per quadrangle, (p - 1)(p - 2)(p - 3)/6 per tetrahedron,
    (p - 1)^2 (p - 2)/2 per prism and (p - 1)^3 per pyramid - the pyramid interior of
    the Fuentes H1 pyramid Palace's MFEM build selects; this closed form reproduces the
    Palace-printed H1 counts of the physics-09 / -11 / gallery-06b hybrid meshes at
    p3 / p4 / p5 exactly (the Bergot count (p - 1)(p - 2)(2p - 3)/6 does not)."""
    order = int(order)
    if order < 1:
        raise ValueError("H1 order must be positive")
    p = order
    return int(counts["Vertices"] + (p - 1) * counts["Edges"]
               + max((p - 1) * (p - 2) // 2, 0) * counts["TriangleFaces"]
               + (p - 1) ** 2 * counts["QuadFaces"]
               + max((p - 1) * (p - 2) * (p - 3) // 6, 0) * counts["Tetrahedra"]
               + max((p - 1) ** 2 * (p - 2) // 2, 0) * counts["Prisms"]
               + (p - 1) ** 3 * counts["Pyramids"])


def h1_dofs(mesh, order):
    """H1 (continuous Lagrange) degrees of freedom of order `order` on the native
    cells (h1_dofs_from_counts of h1_entity_counts)."""
    return h1_dofs_from_counts(h1_entity_counts(mesh), order)
