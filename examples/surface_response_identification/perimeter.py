# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Metal perimeter partition E of a Palace mesh, computed independently of the classifier.

The perimeter follows the definitions of palace/utils/metaledge.cpp ExtractMetalEdgeGeometry
and palace/utils/geodata.cpp GetBoundaryEdgeSegments:

* metal = the union of every boundary attribute with a metal boundary condition (PEC/Ground,
  AuxPEC, Terminal, PrescribedPotential (+TerminalAttributes), Conductivity, Impedance,
  RationalImpedance); lumped ports are not metal;
* a perimeter edge is a mesh edge of a metal boundary face with odd incidence among the
  metal faces (an edge shared by two metal faces of different attributes is an interior seam
  and is not perimeter), one segment per mesh edge;
* a perimeter edge which also lies on an exterior, non-metal, non-interface boundary
  attribute is a TRUNCATION edge (cut by the simulation domain), otherwise PHYSICAL;
* the perimeter is computed on the mesh as meshed (before Palace cracks the interior
  boundary elements), so metal sheets with the same material on both sides (airbridge
  spans and walls) keep their boundary edges here while the cracked-mesh odd-incidence
  rule of the classifier cancels them (both crack copies fall in the same side component);
* a degree-2 perimeter vertex is a CORNER when the turn exceeds the classifier's
  30 degree tolerance (metaledge.cpp corner_angle_tolerance_degrees), REGULAR otherwise;
  degree 1 is an ENDPOINT, degree >= 3 a JUNCTION; a physical chain is a maximal path
  through REGULAR vertices.

In addition (decision 73(3)) perimeter edges owned only by metal faces whose normal is not
parallel to the process normal are the excluded NONPLANAR class, planar edges off the primary
metal plane (the plane of largest metal area) are the excluded CROSS_LAYER class, and edges
owned by three or more metal faces (a wall standing on a sheet) are the NONMANIFOLD class.
"""

import math
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from .msh2 import ELEMENT_DIMENSION, QUAD, TETRAHEDRON_TYPES, TRIANGLE_TYPES

CORNER_ANGLE_TOLERANCE_DEGREES = 30.0
DIRECTION_QUANTUM = 1.0e-12
# Faces whose unit normal deviates from the process normal by more than this are non-planar
# metal (sidewalls, staples, TSV walls).
PLANAR_COSINE_TOLERANCE = 1.0e-6


def metal_attributes(config):
    """Attribute -> boundary-law type, following metaledge.cpp ExtractMetalEdgeGeometry."""
    boundaries = config.get("Boundaries", {})
    result = {}

    def add(attributes, law):
        for attribute in attributes or []:
            result.setdefault(int(attribute), law)

    for key in ("PEC", "Ground", "AuxPEC"):
        add(boundaries.get(key, {}).get("Attributes"), "PEC")
    for terminal in boundaries.get("Terminal", []):
        add(terminal.get("Attributes"), "PEC")
    for potential in boundaries.get("PrescribedPotential", []):
        add(potential.get("Attributes"), "PEC")
        add(potential.get("TerminalAttributes"), "PEC")
    for entry in boundaries.get("Conductivity", []):
        add(entry.get("Attributes"), "Conductivity")
    for entry in boundaries.get("Impedance", []):
        add(entry.get("Attributes"), "Impedance")
    for entry in boundaries.get("RationalImpedance", []):
        add(entry.get("Attributes"), "RationalImpedance")
    return result


def conductor_of_attribute(config, attribute):
    """Conductor index as GetConductor assigns it: 0 for every PEC attribute, the Terminal or
    PrescribedPotential index otherwise (surfaceresponseoperator.cpp GetConductor)."""
    boundaries = config.get("Boundaries", {})
    for key in ("PEC", "Ground", "AuxPEC"):
        if attribute in (boundaries.get(key, {}).get("Attributes") or []):
            return 0
    for terminal in boundaries.get("Terminal", []):
        if attribute in (terminal.get("Attributes") or []):
            return int(terminal["Index"])
    for potential in boundaries.get("PrescribedPotential", []):
        if attribute in (potential.get("Attributes") or []) or attribute in (
            potential.get("TerminalAttributes") or []
        ):
            return int(potential["Index"])
    return None


def interface_attributes(config):
    """(index, type) -> attribute list of every typed Boundaries.Postprocessing.Dielectric."""
    result = {}
    for entry in config.get("Boundaries", {}).get("Postprocessing", {}).get("Dielectric", []):
        kind = entry.get("Type", "Default")
        if kind in ("SA", "MS", "MA"):
            result[(int(entry["Index"]), kind)] = [int(a) for a in entry.get("Attributes", [])]
    return result


def target_interfaces(config):
    block = config.get("Solver", {}).get("SurfaceResponseCorrection")
    if block is None:
        block = config.get("Solver", {}).get("Electrostatic", {}).get("ResponseCorrection", {})
    return [int(i) for i in (block or {}).get("TargetInterfaces", [])]


@dataclass
class PerimeterVertex:
    point: np.ndarray
    edges: list = field(default_factory=list)
    kind: str = "REGULAR"  # REGULAR | CORNER | ENDPOINT | JUNCTION
    physical_kind: str = None  # same, counting PHYSICAL edges only
    turn_degrees: float = None  # deviation from straight, degree-2 vertices only
    convex: bool = None  # metal on the inside of the turn


@dataclass
class PerimeterEdge:
    vertices: tuple
    length: float
    attributes: tuple  # metal attributes of the faces owning this edge
    kind: str  # PHYSICAL | TRUNCATION | NONPLANAR | CROSS_LAYER | NONMANIFOLD
    conductors: tuple
    interfaces: tuple  # (index, type) of the typed dielectric interfaces coincident with it
    inward: np.ndarray  # in-plane unit vector from the edge into the metal
    plane: int = 0
    owners: int = 1  # number of metal faces owning the edge
    chain: int = -1


@dataclass
class Perimeter:
    vertices: list
    edges: list
    process_normal: np.ndarray
    planes: list  # plane offsets along the process normal, one per metal layer
    plane_areas: list = field(default_factory=list)
    primary_plane: int = 0
    chains: int = 0
    nonplanar_face_area: float = 0.0
    planar_face_area: float = 0.0
    scale: float = 1.0

    def edge_points(self, edge):
        return self.vertices[edge.vertices[0]].point, self.vertices[edge.vertices[1]].point

    def length(self, kind=None):
        return float(sum(e.length for e in self.edges if kind is None or e.kind == kind))


def _face_edges(corners):
    k = corners.shape[1]
    return [np.sort(np.stack([corners[:, i], corners[:, (i + 1) % k]], axis=1), axis=1) for i in range(k)]


def _face_normals(coordinates, corners):
    a, b, c = coordinates[corners[:, 0]], coordinates[corners[:, 1]], coordinates[corners[:, 2]]
    n = np.cross(b - a, c - a)
    area = 0.5 * np.linalg.norm(n, axis=1)
    with np.errstate(invalid="ignore", divide="ignore"):
        unit = n / (2.0 * area)[:, None]
    return unit, area


def _exterior_faces(mesh):
    """Set of sorted-node face keys which belong to exactly one volume element."""
    counts = defaultdict(int)
    for element_type in mesh.elements:
        if ELEMENT_DIMENSION[element_type] != 3:
            continue
        corners = mesh.corner_indices(element_type)
        if element_type in TETRAHEDRON_TYPES:
            faces = [corners[:, [0, 1, 2]], corners[:, [0, 1, 3]], corners[:, [0, 2, 3]], corners[:, [1, 2, 3]]]
        else:
            raise NotImplementedError("only tetrahedral volume elements are supported by the audit")
        for f in faces:
            for key in map(tuple, np.sort(f, axis=1)):
                counts[key] += 1
    return {key for key, count in counts.items() if count == 1}


def _canonical_node_map(mesh, tolerance):
    """Map every node index to a representative index, merging geometrically coincident
    nodes (duplicated crack nodes in a pre-cracked mesh)."""
    keys = np.round(mesh.coordinates / tolerance).astype(np.int64)
    representative = {}
    result = np.arange(len(mesh.coordinates))
    for index, key in enumerate(map(tuple, keys)):
        result[index] = representative.setdefault(key, index)
    return result


def extract_perimeter(mesh, config, process_normal=None, corner_tolerance_degrees=CORNER_ANGLE_TOLERANCE_DEGREES):
    metal = metal_attributes(config)
    if not metal:
        raise ValueError("the configuration names no metal boundary attributes")
    interfaces = interface_attributes(config)
    interface_attribute_set = {a for attributes in interfaces.values() for a in attributes}

    bbox = mesh.coordinates.max(axis=0) - mesh.coordinates.min(axis=0)
    extent = float(bbox.max())
    tolerance = 1.0e-10 * extent
    canonical = _canonical_node_map(mesh, tolerance)

    if mesh.has(QUAD):
        raise NotImplementedError("quadrilateral boundary faces are not supported by the audit")
    # First- and second-order triangles (Gmsh types 2 and 9); the perimeter is the straight
    # corner-to-corner edge graph, as in the classifier's mesh-vertex segments.
    triangle_types = [t for t in TRIANGLE_TYPES if mesh.has(t)]
    if not triangle_types:
        raise ValueError("the mesh has no triangular boundary faces")
    physical = np.concatenate([mesh.physical_tags(t) for t in triangle_types])
    corners_raw = np.concatenate([mesh.corner_indices(t) for t in triangle_types])
    corners = canonical[corners_raw]
    metal_mask = np.isin(physical, list(metal))
    metal_corners = corners[metal_mask]
    metal_tags = physical[metal_mask]
    if metal_corners.shape[0] == 0:
        raise ValueError(f"no boundary faces carry the metal attributes {sorted(metal)}")

    normals, areas = _face_normals(mesh.coordinates, metal_corners)
    if process_normal is None:
        # Dominant orientation by area: the process normal is the area-weighted principal
        # direction of the metal face normals (sign-invariant).
        m = (normals * areas[:, None]).T @ normals
        eigenvalues, eigenvectors = np.linalg.eigh(m)
        process_normal = eigenvectors[:, int(np.argmax(eigenvalues))]
    process_normal = np.asarray(process_normal, dtype=float)
    process_normal /= np.linalg.norm(process_normal)
    cosines = np.abs(normals @ process_normal)
    planar = cosines >= 1.0 - PLANAR_COSINE_TOLERANCE

    # Plane offsets of the planar metal along the process normal (one per metal layer).
    offsets = mesh.coordinates[metal_corners[planar]] @ process_normal
    offsets = offsets.mean(axis=1)
    plane_values = []
    plane_of_face = np.full(metal_corners.shape[0], -1)
    plane_tolerance = 1.0e-6 * extent
    planar_indices = np.flatnonzero(planar)
    for local, face_index in enumerate(planar_indices):
        value = offsets[local]
        for plane_index, plane_value in enumerate(plane_values):
            if abs(value - plane_value) <= plane_tolerance:
                plane_of_face[face_index] = plane_index
                break
        else:
            plane_values.append(float(value))
            plane_of_face[face_index] = len(plane_values) - 1

    plane_areas = [float(areas[plane_of_face == i].sum()) for i in range(len(plane_values))]
    primary_plane = int(np.argmax(plane_areas)) if plane_areas else -1

    # Perimeter edges: odd incidence among the metal faces.
    incidence = defaultdict(int)
    owners = defaultdict(list)
    for edge_key_array in _face_edges(metal_corners):
        for face_index, key in enumerate(map(tuple, edge_key_array)):
            incidence[key] += 1
            owners[key].append(face_index)
    perimeter_keys = sorted(key for key, count in incidence.items() if count % 2 == 1)

    # Truncation: exterior non-metal non-interface boundary faces.
    exterior = _exterior_faces(mesh) if any(ELEMENT_DIMENSION[t] == 3 for t in mesh.elements) else set()
    truncation_edges = set()
    interface_edges = defaultdict(set)
    other_mask = ~metal_mask
    for face_index in np.flatnonzero(other_mask):
        tag = int(physical[face_index])
        face = corners[face_index]
        face_key = tuple(np.sort(corners_raw[face_index]))
        is_truncation = tag not in interface_attribute_set and face_key in exterior
        for i in range(3):
            key = tuple(sorted((int(face[i]), int(face[(i + 1) % 3]))))
            if is_truncation:
                truncation_edges.add(key)
            for interface, attributes in interfaces.items():
                if tag in attributes:
                    interface_edges[key].add(interface)
    # Interfaces defined on metal attributes themselves (MS/MA on the PEC attribute).
    metal_interface = defaultdict(set)
    for interface, attributes in interfaces.items():
        for attribute in attributes:
            if attribute in metal:
                metal_interface[attribute].add(interface)

    vertex_index = {}
    vertices = []

    def vertex(node):
        if node not in vertex_index:
            vertex_index[node] = len(vertices)
            vertices.append(PerimeterVertex(point=mesh.coordinates[node].copy()))
        return vertex_index[node]

    edges = []
    for key in perimeter_keys:
        faces = owners[key]
        attributes = tuple(sorted({int(metal_tags[f]) for f in faces}))
        p0, p1 = mesh.coordinates[key[0]], mesh.coordinates[key[1]]
        length = float(np.linalg.norm(p1 - p0))
        if length <= 0.0:
            continue
        face_planar = [bool(planar[f]) for f in faces]
        plane = int(max((plane_of_face[f] for f in faces if plane_of_face[f] >= 0), default=-1))
        if len(faces) >= 3:
            kind = "NONMANIFOLD"
        elif not any(face_planar):
            kind = "NONPLANAR"
        elif plane != primary_plane:
            kind = "CROSS_LAYER"
        elif key in truncation_edges:
            kind = "TRUNCATION"
        else:
            kind = "PHYSICAL"
        # Inward direction: from the edge midpoint towards the owning face centroid,
        # projected onto the metal plane (planar faces only).
        inward = np.zeros(3)
        for f in faces:
            if planar[f]:
                centroid = mesh.coordinates[metal_corners[f]].mean(axis=0)
                d = centroid - 0.5 * (p0 + p1)
                d -= (d @ process_normal) * process_normal
                t = (p1 - p0) / length
                d -= (d @ t) * t
                norm = np.linalg.norm(d)
                if norm > 0:
                    inward += d / norm
        norm = np.linalg.norm(inward)
        inward = inward / norm if norm > 0 else inward
        edge_interfaces = set(interface_edges.get(key, set()))
        for attribute in attributes:
            edge_interfaces |= metal_interface.get(attribute, set())
        conductors = tuple(sorted({conductor_of_attribute(config, a) for a in attributes if conductor_of_attribute(config, a) is not None}))
        v0, v1 = vertex(key[0]), vertex(key[1])
        edge = PerimeterEdge(
            vertices=(v0, v1),
            length=length,
            attributes=attributes,
            kind=kind,
            conductors=conductors,
            interfaces=tuple(sorted(edge_interfaces)),
            inward=inward,
            plane=plane,
            owners=len(faces),
        )
        vertices[v0].edges.append(len(edges))
        vertices[v1].edges.append(len(edges))
        edges.append(edge)

    perimeter = Perimeter(
        vertices=vertices,
        edges=edges,
        process_normal=process_normal,
        planes=plane_values,
        plane_areas=plane_areas,
        primary_plane=primary_plane,
        nonplanar_face_area=float(areas[~planar].sum()),
        planar_face_area=float(areas[planar].sum()),
    )
    classify_vertices(perimeter, corner_tolerance_degrees)
    label_chains(perimeter)
    return perimeter


def _quantize(cosine):
    return round(cosine / DIRECTION_QUANTUM)


def classify_vertices(perimeter, corner_tolerance_degrees=CORNER_ANGLE_TOLERANCE_DEGREES):
    """metaledge.cpp ClassifyVertex: quantized direction cosines, 30 degree turn tolerance."""
    straight_dot = -math.cos(math.radians(corner_tolerance_degrees))
    quantized_straight = _quantize(straight_dot)
    for index, vertex in enumerate(perimeter.vertices):
        for physical in (False, True):
            edges = [
                e for e in vertex.edges if not physical or perimeter.edges[e].kind == "PHYSICAL"
            ]
            if not edges:
                kind = None
            elif len(edges) == 1:
                kind = "ENDPOINT"
            elif len(edges) > 2:
                kind = "JUNCTION"
            else:
                directions = []
                for e in edges:
                    edge = perimeter.edges[e]
                    other = edge.vertices[1] if edge.vertices[0] == index else edge.vertices[0]
                    d = perimeter.vertices[other].point - vertex.point
                    directions.append(d / np.linalg.norm(d))
                dot = float(directions[0] @ directions[1])
                kind = "REGULAR" if _quantize(dot) <= quantized_straight else "CORNER"
                if not physical:
                    vertex.turn_degrees = 180.0 - math.degrees(math.acos(max(-1.0, min(1.0, dot))))
                    bisector = directions[0] + directions[1]
                    inward = perimeter.edges[edges[0]].inward + perimeter.edges[edges[1]].inward
                    if np.linalg.norm(bisector) > 0 and np.linalg.norm(inward) > 0:
                        # Metal on the inside of the turn (bisector points into the metal)
                        # is a convex metal corner.
                        vertex.convex = bool(bisector @ inward > 0)
            if physical:
                vertex.physical_kind = kind
            else:
                vertex.kind = kind


def label_chains(perimeter):
    """Maximal PHYSICAL paths through REGULAR vertices (metaledge.cpp physical chains)."""
    chain = 0
    for seed, edge in enumerate(perimeter.edges):
        if edge.chain >= 0 or edge.kind != "PHYSICAL":
            continue
        stack = [seed]
        edge.chain = chain
        while stack:
            current = stack.pop()
            for v in perimeter.edges[current].vertices:
                vertex = perimeter.vertices[v]
                if vertex.physical_kind != "REGULAR":
                    continue
                for neighbor in vertex.edges:
                    other = perimeter.edges[neighbor]
                    if other.chain < 0 and other.kind == "PHYSICAL":
                        other.chain = chain
                        stack.append(neighbor)
        chain += 1
    perimeter.chains = chain


def segment_distance(p0, p1, q0, q1):
    """Closest distance between segments p0p1 and q0q1 (3D), by sampling the closed form
    for parallel and non-parallel cases."""
    u = p1 - p0
    v = q1 - q0
    w = p0 - q0
    a, b, c, d, e = u @ u, u @ v, v @ v, u @ w, v @ w
    denominator = a * c - b * b
    if denominator <= 1.0e-14 * a * c:
        s_candidates = [0.0, 1.0]
        best = math.inf
        for s in s_candidates:
            p = p0 + s * u
            t = min(1.0, max(0.0, ((p - q0) @ v) / c))
            best = min(best, float(np.linalg.norm(p - (q0 + t * v))))
        for t in (0.0, 1.0):
            q = q0 + t * v
            s = min(1.0, max(0.0, ((q - p0) @ u) / a))
            best = min(best, float(np.linalg.norm(q - (p0 + s * u))))
        return best
    s = (b * e - c * d) / denominator
    t = (a * e - b * d) / denominator
    s = min(1.0, max(0.0, s))
    t = ((p0 + s * u - q0) @ v) / c
    t = min(1.0, max(0.0, t))
    s = ((q0 + t * v - p0) @ u) / a
    s = min(1.0, max(0.0, s))
    return float(np.linalg.norm((p0 + s * u) - (q0 + t * v)))


def edge_interactions(perimeter, radius, kinds=("PHYSICAL",)):
    """Pairs of non-adjacent perimeter edges of different chains within 2R of each other:
    (edge a, edge b, distance, parallel cosine). A brute-force grid search; meshes of
    O(10^4) perimeter edges take seconds."""
    indices = [i for i, e in enumerate(perimeter.edges) if e.kind in kinds]
    if not indices:
        return []
    interaction = 2.0 * radius
    points = np.array([[*perimeter.edge_points(perimeter.edges[i])] for i in indices])  # (n, 2, 3)
    midpoints = points.mean(axis=1)
    half = 0.5 * np.linalg.norm(points[:, 1] - points[:, 0], axis=1)
    cell = interaction + 2.0 * float(half.max())
    grid = defaultdict(list)
    keys = np.floor(midpoints / cell).astype(np.int64)
    for local, key in enumerate(map(tuple, keys)):
        grid[key].append(local)
    results = []
    offsets = [(i, j, k) for i in (-1, 0, 1) for j in (-1, 0, 1) for k in (-1, 0, 1)]
    for local_a, key in enumerate(map(tuple, keys)):
        a = indices[local_a]
        edge_a = perimeter.edges[a]
        for offset in offsets:
            for local_b in grid.get((key[0] + offset[0], key[1] + offset[1], key[2] + offset[2]), []):
                if local_b <= local_a:
                    continue
                b = indices[local_b]
                edge_b = perimeter.edges[b]
                if set(edge_a.vertices) & set(edge_b.vertices):
                    continue
                if edge_a.chain == edge_b.chain and edge_a.chain >= 0:
                    continue
                if np.linalg.norm(midpoints[local_a] - midpoints[local_b]) > interaction + half[local_a] + half[local_b]:
                    continue
                p0, p1 = points[local_a]
                q0, q1 = points[local_b]
                distance = segment_distance(p0, p1, q0, q1)
                if distance <= interaction * (1.0 + 1.0e-9):
                    ta = (p1 - p0) / (2.0 * half[local_a])
                    tb = (q1 - q0) / (2.0 * half[local_b])
                    results.append((a, b, distance, abs(float(ta @ tb))))
    return results
