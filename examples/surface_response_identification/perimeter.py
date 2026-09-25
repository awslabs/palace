# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Metal perimeter partition E of a Palace mesh, computed independently of the classifier.

The perimeter follows the definitions of palace/utils/metaledge.cpp ExtractMetalEdgeGeometry
(phase 3 of the identification fix: crack-independent, derived from the distinct geometric
metal faces):

* metal = the union of every boundary attribute with a metal boundary condition (PEC/Ground,
  AuxPEC, Terminal, PrescribedPotential (+TerminalAttributes), Conductivity, Impedance,
  RationalImpedance); lumped ports are not metal;
* a mesh edge of the metal faces is classified by the in-plane inward directions of the
  distinct faces owning it: one direction class -> a one-sided perimeter edge; two opposite
  classes -> the metal continues (interior, not perimeter); two non-opposite classes -> a
  FOLD (the metal turns around the edge: staple, box edge; excluded class NonPlanar); three
  or more -> NONMANIFOLD (a wall standing on a sheet). Coincident nodes are merged, so a
  pre-cracked or an as-meshed mesh give the same perimeter;
* a one-sided edge on an exterior, non-metal, non-interface boundary attribute is a
  TRUNCATION edge (cut by the simulation domain); one whose face is not parallel to the
  process normal is NONPLANAR (excluded); one whose faces have the same material on both
  sides (an airbridge span, metal embedded in one dielectric) is EMBEDDED (excluded class
  UndeterminedProcessSide unless the target interfaces configure an EdgeFrameNormal);
  otherwise PHYSICAL;
* the parts of a PHYSICAL edge within 2R of metal off its own plane (a facing layer, a
  wall, a staple) are its CROSS_LAYER portions (excluded, decision 73(3)); a feature vertex
  within 2R of such metal is an excluded vertex;
* a degree-2 perimeter vertex is a CORNER when the turn exceeds the classifier's
  30 degree tolerance (metaledge.cpp corner_angle_tolerance_degrees), REGULAR otherwise;
  degree 1 is an ENDPOINT, degree >= 3 a JUNCTION; a physical chain is a maximal PHYSICAL
  path through REGULAR vertices;
* the conductor of an edge is its edge-connected metal component (the geometric identity the
  classifier uses); the Terminal / Ground labels are reported for comparison only.
"""

import math
from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from .msh2 import ELEMENT_DIMENSION, QUAD, TETRAHEDRON_TYPES, TRIANGLE_TYPES

CORNER_ANGLE_TOLERANCE_DEGREES = 30.0
DIRECTION_QUANTUM = 1.0e-12
# Faces whose unit normal deviates from the process normal by more than this are non-planar
# metal (sidewalls, staples, TSV walls); the classifier's kParallelCosineTolerance.
PLANAR_COSINE_TOLERANCE = 1.0e-8
# Metal off an edge's plane within this multiple of R excludes that part of the edge
# (Identification.Conventions CrossLayerReachOverR).
CROSS_LAYER_REACH_OVER_R = 2.0
# Edge kinds the classifier extracts as PHYSICAL segments (FOLD / NONMANIFOLD are its own
# exclusion types).
CLASSIFIER_PHYSICAL_KINDS = ("PHYSICAL", "TRUNCATION", "EMBEDDED", "NONPLANAR", "BOX")


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
    kind: str  # PHYSICAL | TRUNCATION | EMBEDDED | NONPLANAR | FOLD | NONMANIFOLD
    conductors: tuple  # boundary-condition labels (comparison only)
    interfaces: tuple  # (index, type) of the typed dielectric interfaces coincident with it
    inward: np.ndarray  # in-plane unit vector from the edge into the metal
    plane: int = 0
    owners: int = 1  # number of distinct metal faces owning the edge
    chain: int = -1
    component: int = -1  # edge-connected metal component (the classifier's conductor)
    cross_layer: list = field(default_factory=list)  # (s0, s1) portions within 2R of off-plane metal
    # Direction classes of the owning faces (ONE_SIDED | FOLD | NONMANIFOLD): a BOX edge where
    # two PEC box faces meet is a FOLD of the classifier (type FOLD, no PHYSICAL segment), so
    # it takes no part in the vertex census (classify_vertices).
    face_classes: str = "ONE_SIDED"

    @property
    def cross_layer_length(self):
        return float(sum(b - a for a, b in self.cross_layer))


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
    components: int = 0
    excluded_vertices: set = field(default_factory=set)  # feature vertices within 2R of off-plane metal

    def edge_points(self, edge):
        return self.vertices[edge.vertices[0]].point, self.vertices[edge.vertices[1]].point

    def length(self, kind=None):
        if kind == "CROSS_LAYER":
            return float(sum(e.cross_layer_length for e in self.edges if e.kind == "PHYSICAL"))
        return float(sum(e.length for e in self.edges if kind is None or e.kind == kind))


class BoxGrid:
    """Uniform grid over axis-aligned boxes (numpy): `query(lower, upper, margin)` returns the
    sorted indices of the items whose cells overlap the enlarged box — a superset of the items
    within `margin` of it; the caller applies its exact test. Acceleration only (the former
    all-items scans gave the same results); chip-scale meshes have 10^6 faces and 10^5 edges."""

    def __init__(self, lower, upper, cell):
        self.lower = np.asarray(lower, dtype=float)
        self.upper = np.asarray(upper, dtype=float)
        self.cell = float(cell)
        self.origin = self.lower.min(axis=0) if len(self.lower) else np.zeros(self.lower.shape[1] if self.lower.ndim == 2 else 3)
        self.cells = defaultdict(list)
        lo = np.floor((self.lower - self.origin) / self.cell).astype(np.int64)
        hi = np.floor((self.upper - self.origin) / self.cell).astype(np.int64)
        for index in range(len(self.lower)):
            ranges = [range(int(lo[index, d]), int(hi[index, d]) + 1) for d in range(self.lower.shape[1])]
            for key in _product(ranges):
                self.cells[key].append(index)

    def query(self, lower, upper, margin=0.0):
        lower = np.asarray(lower, dtype=float) - margin
        upper = np.asarray(upper, dtype=float) + margin
        lo = np.floor((lower - self.origin) / self.cell).astype(np.int64)
        hi = np.floor((upper - self.origin) / self.cell).astype(np.int64)
        found = []
        for key in _product([range(int(lo[d]), int(hi[d]) + 1) for d in range(len(lo))]):
            found.extend(self.cells.get(key, ()))
        return np.unique(np.asarray(found, dtype=np.int64))


def _product(ranges):
    if len(ranges) == 1:
        for a in ranges[0]:
            yield (a,)
        return
    for a in ranges[0]:
        for rest in _product(ranges[1:]):
            yield (a, *rest)


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


def _point_polygon_distance(p, polygon, normal):
    """Distance from a point to a planar (convex) polygon: the normal distance when the
    projection falls inside, else the distance to the nearest edge."""
    n = len(polygon)
    inside = True
    for i in range(n):
        a, b, c = polygon[i], polygon[(i + 1) % n], polygon[(i + 2) % n]
        interior = np.cross(normal, b - a)
        if float((c - b) @ interior) * float((p - a) @ interior) < 0.0:
            inside = False
            break
    if inside:
        return abs(float((p - polygon[0]) @ normal))
    best = math.inf
    for i in range(n):
        a, b = polygon[i], polygon[(i + 1) % n]
        d = b - a
        t = min(1.0, max(0.0, float((p - a) @ d) / float(d @ d)))
        best = min(best, float(np.linalg.norm(p - (a + t * d))))
    return best


def _sublevel_interval(f, lo, hi, level):
    """[s0, s1] on which the convex function f < level (golden section + bisection, the
    classifier's ConvexSublevelInterval), or None."""
    golden = 0.6180339887498949
    a, b = lo, hi
    c, d = b - golden * (b - a), a + golden * (b - a)
    fc, fd = f(c), f(d)
    for _ in range(90):
        if b - a <= 1.0e-14 * (hi - lo):
            break
        if fc < fd:
            b, d, fd = d, c, fc
            c = b - golden * (b - a)
            fc = f(c)
        else:
            a, c, fc = c, d, fd
            d = a + golden * (b - a)
            fd = f(d)
    s_min = 0.5 * (a + b)
    f_min = f(s_min)
    for candidate in (lo, hi):
        value = f(candidate)
        if value < f_min:
            f_min, s_min = value, candidate
    if not f_min < level:
        return None

    def bisect(inside, outside):
        if f(outside) < level:
            return outside
        for _ in range(100):
            if abs(outside - inside) <= 1.0e-15 * (hi - lo):
                break
            mid = 0.5 * (inside + outside)
            if f(mid) < level:
                inside = mid
            else:
                outside = mid
        return inside

    return (bisect(s_min, lo), bisect(s_min, hi))


def sample_quadratic_edge(p0, pm, p1):
    """Palace's sampling of a high-order edge (geodata.cpp SampleEdgeTransformation): the
    quadratic Lagrange edge through p0 (t = 0), pm (t = 1/2), p1 (t = 1) is bisected while the
    tangent turns by more than 7.5 degrees or the middle point deviates from the chord by more
    than 1e-3 of the endpoint distance (depth <= 12). Returns the sampled points p0 .. p1."""
    p0, pm, p1 = (np.asarray(x, dtype=float) for x in (p0, pm, p1))

    def evaluate(t):
        point = (1.0 - t) * (1.0 - 2.0 * t) * p0 + 4.0 * t * (1.0 - t) * pm + t * (2.0 * t - 1.0) * p1
        tangent = (4.0 * t - 3.0) * p0 + (4.0 - 8.0 * t) * pm + (4.0 * t - 1.0) * p1
        norm = np.linalg.norm(tangent)
        return t, point, (tangent / norm if norm > 1.0e-14 else None)

    def angle(a, b):
        if a[2] is None or b[2] is None:
            return 0.0
        return math.acos(max(-1.0, min(1.0, float(a[2] @ b[2]))))

    def chord_distance_squared(point, a, b):
        d = b - a
        length_squared = float(d @ d)
        if length_squared <= 1.0e-28:
            return float((point - a) @ (point - a))
        t = min(1.0, max(0.0, float((point - a) @ d) / length_squared))
        delta = point - (a + t * d)
        return float(delta @ delta)

    maximum_turn = math.radians(7.5)
    edge_scale = max(float(np.linalg.norm(p1 - p0)), 1.0e-14)
    chord_tolerance_squared = (1.0e-3 * edge_scale) ** 2
    first, last = evaluate(0.0), evaluate(1.0)
    points = [first[1]]

    def refine(left, right, depth):
        middle = evaluate(0.5 * (left[0] + right[0]))
        excessive_turn = angle(left, middle) > maximum_turn or angle(middle, right) > maximum_turn
        excessive_deviation = chord_distance_squared(middle[1], left[1], right[1]) > chord_tolerance_squared
        if depth < 12 and (excessive_turn or excessive_deviation):
            refine(left, middle, depth + 1)
            refine(middle, right, depth + 1)
        else:
            points.append(right[1])

    refine(first, last, 0)
    return points


def extract_perimeter(mesh, config, process_normal=None, corner_tolerance_degrees=CORNER_ANGLE_TOLERANCE_DEGREES, radius=None):
    metal = metal_attributes(config)
    if not metal:
        raise ValueError("the configuration names no metal boundary attributes")
    interfaces = interface_attributes(config)
    interface_attribute_set = {a for attributes in interfaces.values() for a in attributes}
    frame_normal_configured = any(
        entry.get("EdgeFrameNormal") for entry in config.get("Boundaries", {}).get("Postprocessing", {}).get("Dielectric", [])
    )

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
    # Mid-edge nodes of the second-order triangles (Gmsh order: corners, then the mid-edge
    # nodes of edges 0-1, 1-2, 2-0), keyed by the canonical corner pair.
    mid_node = {}
    if mesh.has(9):
        nodes9 = mesh.node_indices(9)
        canonical9 = canonical[nodes9[:, :3]]
        for local, (i, j) in enumerate(((0, 1), (1, 2), (2, 0))):
            for a, b, m in zip(canonical9[:, i], canonical9[:, j], nodes9[:, 3 + local]):
                mid_node[(min(int(a), int(b)), max(int(a), int(b)))] = int(m)
    metal_mask = np.isin(physical, list(metal))
    metal_corners = corners[metal_mask]
    metal_tags = physical[metal_mask]
    if metal_corners.shape[0] == 0:
        raise ValueError(f"no boundary faces carry the metal attributes {sorted(metal)}")

    # Distinct geometric faces (coincident copies count once).
    face_keys = [tuple(sorted(map(int, f))) for f in metal_corners]
    distinct = {}
    for index, key in enumerate(face_keys):
        distinct.setdefault(key, index)
    face_ids = np.array([distinct[key] for key in face_keys])
    unique_faces = sorted(set(face_ids.tolist()))
    metal_corners = metal_corners[unique_faces]
    metal_tags = metal_tags[unique_faces]

    normals, areas = _face_normals(mesh.coordinates, metal_corners)
    centroids = mesh.coordinates[metal_corners].mean(axis=1)
    # Metal faces on the bounding box of the mesh are a PEC simulation box, not process
    # metal: they neither vote for the process normal nor count as off-plane metal.
    lower_box, upper_box = mesh.coordinates.min(axis=0), mesh.coordinates.max(axis=0)
    face_points = mesh.coordinates[metal_corners]  # (n, 3, 3)
    on_box = np.zeros(metal_corners.shape[0], dtype=bool)
    for d in range(3):
        for bound in (lower_box[d], upper_box[d]):
            on_box |= np.all(np.abs(face_points[:, :, d] - bound) <= tolerance, axis=1)
    if process_normal is None:
        # Dominant orientation by area: the process normal is the area-weighted principal
        # direction of the metal face normals (sign-invariant); metaledge.cpp layer_normal.
        vote = ~on_box if not np.all(on_box) else np.ones_like(on_box)
        m = (normals[vote] * areas[vote, None]).T @ normals[vote]
        eigenvalues, eigenvectors = np.linalg.eigh(m)
        process_normal = eigenvectors[:, int(np.argmax(eigenvalues))]
    process_normal = np.asarray(process_normal, dtype=float)
    process_normal /= np.linalg.norm(process_normal)
    cosines = np.abs(normals @ process_normal)
    planar = cosines >= 1.0 - PLANAR_COSINE_TOLERANCE

    # Materials adjacent to every metal face (embedded sheets have one material).
    face_materials = defaultdict(set)
    for element_type in mesh.elements:
        if ELEMENT_DIMENSION[element_type] != 3:
            continue
        tet_corners = canonical[mesh.corner_indices(element_type)]
        tet_tags = mesh.physical_tags(element_type)
        if element_type not in TETRAHEDRON_TYPES:
            raise NotImplementedError("only tetrahedral volume elements are supported by the audit")
        for f in ([0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]):
            for key, tag in zip(map(tuple, np.sort(tet_corners[:, f], axis=1)), tet_tags):
                if key in distinct:
                    face_materials[key].add(int(tag))
    face_key_of = [tuple(sorted(map(int, f))) for f in metal_corners]

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

    # Edge incidence: the distinct faces owning every mesh edge, with their in-plane inward
    # directions; direction classes decide the kind (metaledge.cpp phase 3).
    owners = defaultdict(list)
    for edge_key_array in _face_edges(metal_corners):
        for face_index, key in enumerate(map(tuple, edge_key_array)):
            owners[key].append(face_index)

    def inward_of(face_index, key):
        p0, p1 = mesh.coordinates[key[0]], mesh.coordinates[key[1]]
        t = p1 - p0
        d = centroids[face_index] - 0.5 * (p0 + p1)
        d = d - (d @ t) / float(t @ t) * t
        return d / np.linalg.norm(d)

    same = _quantize(1.0 - 1.0e-8)
    opposite = _quantize(-1.0 + 1.0e-8)
    edge_kinds = {}
    edge_inwards = {}
    for key, faces in owners.items():
        classes = []
        inwards = [inward_of(f, key) for f in faces]
        for d in inwards:
            if not any(_quantize(float(d @ c)) >= same for c in classes):
                classes.append(d)
        if len(classes) == 1:
            kind = "ONE_SIDED"
        elif len(classes) == 2:
            if _quantize(float(classes[0] @ classes[1])) <= opposite:
                continue
            kind = "FOLD"
        else:
            kind = "NONMANIFOLD"
        edge_kinds[key] = kind
        edge_inwards[key] = inwards
    perimeter_keys = sorted(edge_kinds)

    # Edge-connected metal components (the classifier's conductor identity).
    parent = list(range(metal_corners.shape[0]))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for key, faces in owners.items():
        for f in faces[1:]:
            a, b = find(faces[0]), find(f)
            if a != b:
                parent[max(a, b)] = min(a, b)
    component_of_root = {}
    face_component = np.array([component_of_root.setdefault(find(f), len(component_of_root)) for f in range(metal_corners.shape[0])])

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
    extra_points = []  # sampled points of curved (second-order) edges

    def vertex(node):
        if node not in vertex_index:
            vertex_index[node] = len(vertices)
            point = mesh.coordinates[node] if node < len(mesh.coordinates) else extra_points[node - len(mesh.coordinates)]
            vertices.append(PerimeterVertex(point=np.array(point, dtype=float)))
        return vertex_index[node]

    def sampled_nodes(key):
        """Node indices along the (possibly curved) edge, as the classifier samples it."""
        m = mid_node.get(key)
        if m is None:
            return [key[0], key[1]]
        p0, pm, p1 = mesh.coordinates[key[0]], mesh.coordinates[m], mesh.coordinates[key[1]]
        sampled = sample_quadratic_edge(p0, pm, p1)
        if len(sampled) <= 2:
            return [key[0], key[1]]
        nodes = [key[0]]
        for point in sampled[1:-1]:
            extra_points.append(point)
            nodes.append(len(mesh.coordinates) + len(extra_points) - 1)
        nodes.append(key[1])
        return nodes

    edges = []
    for key in perimeter_keys:
        faces = owners[key]
        attributes = tuple(sorted({int(metal_tags[f]) for f in faces}))
        plane = int(max((plane_of_face[f] for f in faces if plane_of_face[f] >= 0), default=-1))
        materials = set()
        for f in faces:
            materials |= face_materials.get(face_key_of[f], set())
        if all(on_box[f] for f in faces) and key not in truncation_edges:
            kind = "BOX"
        elif edge_kinds[key] == "NONMANIFOLD":
            kind = "NONMANIFOLD"
        elif edge_kinds[key] == "FOLD":
            kind = "FOLD"
        elif not all(planar[f] for f in faces):
            kind = "NONPLANAR"
        elif key in truncation_edges:
            kind = "TRUNCATION"
        elif len(materials) == 1 and not frame_normal_configured:
            # One material on both sides (no adjacent volume elements at all leaves the
            # process side to the classifier's material scores: PHYSICAL).
            kind = "EMBEDDED"
        else:
            kind = "PHYSICAL"
        inward = np.zeros(3)
        for d in edge_inwards[key]:
            inward += d
        norm = np.linalg.norm(inward)
        inward = inward / norm if norm > 0 else inward
        edge_interfaces = set(interface_edges.get(key, set()))
        for attribute in attributes:
            edge_interfaces |= metal_interface.get(attribute, set())
        conductors = tuple(sorted({conductor_of_attribute(config, a) for a in attributes if conductor_of_attribute(config, a) is not None}))
        nodes = sampled_nodes(key)
        for a, b in zip(nodes[:-1], nodes[1:]):
            v0, v1 = vertex(a), vertex(b)
            length = float(np.linalg.norm(vertices[v1].point - vertices[v0].point))
            if length <= 0.0:
                continue
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
                component=int(face_component[faces[0]]),
                face_classes=edge_kinds[key],
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
        components=len(component_of_root),
    )
    classify_vertices(perimeter, corner_tolerance_degrees)
    label_chains(perimeter)
    if radius is not None:
        _cross_layer_zones(perimeter, mesh, metal_corners[~on_box], normals[~on_box], planar[~on_box], plane_of_face[~on_box], plane_values, radius)
    return perimeter


def _cross_layer_zones(perimeter, mesh, metal_corners, normals, planar, plane_of_face, plane_values, radius):
    """Portions of the PHYSICAL edges within 2R of metal off their plane, and the feature
    vertices within 2R of such metal (decision 73(3); Identifier::ClassifyPlanes)."""
    reach = CROSS_LAYER_REACH_OVER_R * radius
    polygons = mesh.coordinates[metal_corners]  # (n, 3, 3)
    lower = polygons.min(axis=1) - reach
    upper = polygons.max(axis=1) + reach
    face_offsets = np.where(plane_of_face >= 0, np.array([plane_values[i] if i >= 0 else 0.0 for i in plane_of_face]), 0.0)
    quantum = 1.0e-9 * radius
    # Faces by their (reach-enlarged) boxes: the candidates of an edge / vertex box are the
    # faces of the overlapping grid cells, then the former exact box and plane tests.
    grid = BoxGrid(lower, upper, 8.0 * radius) if len(lower) else None

    def candidates(box_lower, box_upper, offset):
        if grid is None:
            return np.zeros(0, dtype=np.int64)
        found = grid.query(box_lower, box_upper)
        near = np.all((lower[found] <= box_upper) & (upper[found] >= box_lower), axis=1)
        off_plane = ~planar[found] | ((np.abs(face_offsets[found] - offset) > quantum) & (np.abs(face_offsets[found] - offset) < reach - quantum))
        return found[near & off_plane]

    for edge in perimeter.edges:
        if edge.kind != "PHYSICAL":
            continue
        p0, p1 = perimeter.edge_points(edge)
        offset = float(0.5 * (p0 + p1) @ perimeter.process_normal)
        zones = []
        for f in candidates(np.minimum(p0, p1), np.maximum(p0, p1), offset):
            polygon = polygons[f]
            normal = normals[f]
            interval = _sublevel_interval(lambda s: _point_polygon_distance(p0 + s * (p1 - p0), polygon, normal), 0.0, 1.0, reach)
            if interval is not None:
                zones.append((interval[0] * edge.length, interval[1] * edge.length))
        merged = []
        for a, b in sorted(zones):
            if merged and a <= merged[-1][1] + quantum:
                merged[-1] = (merged[-1][0], max(merged[-1][1], b))
            else:
                merged.append((a, b))
        edge.cross_layer = merged
    for index, v in enumerate(perimeter.vertices):
        if v.physical_kind in (None, "REGULAR"):
            continue
        if not any(perimeter.edges[e].kind == "PHYSICAL" for e in v.edges):
            continue
        offset = float(v.point @ perimeter.process_normal)
        for f in candidates(v.point, v.point, offset):
            if _point_polygon_distance(v.point, polygons[f], normals[f]) < reach - quantum:
                perimeter.excluded_vertices.add(index)
                break


def _quantize(cosine):
    return round(cosine / DIRECTION_QUANTUM)


def classify_vertices(perimeter, corner_tolerance_degrees=CORNER_ANGLE_TOLERANCE_DEGREES):
    """metaledge.cpp ClassifyVertex: quantized direction cosines, 30 degree turn tolerance."""
    straight_dot = -math.cos(math.radians(corner_tolerance_degrees))
    quantized_straight = _quantize(straight_dot)
    # physical_kind counts the segments the classifier types PHYSICAL (one-sided, not on the
    # truncation boundary), i.e. the audit's PHYSICAL, EMBEDDED, NONPLANAR and one-sided BOX
    # kinds. Vertex-census rule (SURFACE-RESPONSE-IDENTIFICATION.md (a) 2): a vertex all of
    # whose incident edges are folds / non-manifold edges (the corner of a PEC box where three
    # box faces meet, the base corners of a bump) is a vertex of no one-sided perimeter and has
    # no record on either side; a BOX edge shared by two box faces is such a fold.
    physical_kinds = tuple(k for k in CLASSIFIER_PHYSICAL_KINDS if k != "TRUNCATION")
    for index, vertex in enumerate(perimeter.vertices):
        for physical in (False, True):
            edges = [
                e
                for e in vertex.edges
                if not physical
                or (perimeter.edges[e].kind in physical_kinds and perimeter.edges[e].face_classes == "ONE_SIDED")
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
    (edge a, edge b, distance, parallel cosine), in (a, b) order (the consumers aggregate)."""
    indices = [i for i, e in enumerate(perimeter.edges) if e.kind in kinds]
    if not indices:
        return []
    interaction = 2.0 * radius
    points = np.array([[*perimeter.edge_points(perimeter.edges[i])] for i in indices])  # (n, 2, 3)
    midpoints = points.mean(axis=1)
    half = 0.5 * np.linalg.norm(points[:, 1] - points[:, 0], axis=1)
    # Edges by their boxes (a cell of 2R; a long straight edge occupies many cells): the
    # candidates of an edge are the edges of the cells within the interaction distance of
    # its box, then the former exact tests; the pairs are listed by (a, b).
    grid = BoxGrid(points.min(axis=1), points.max(axis=1), interaction)
    results = []
    for local_a in range(len(indices)):
        a = indices[local_a]
        edge_a = perimeter.edges[a]
        for local_b in grid.query(points[local_a].min(axis=0), points[local_a].max(axis=0), interaction * (1.0 + 1.0e-6)):
            local_b = int(local_b)
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


def chain_vertex_sequences(perimeter):
    """Ordered vertex sequences of every PHYSICAL chain (a maximal path through REGULAR
    vertices); a closed chain without any non-regular vertex starts at its lowest vertex."""
    sequences = []
    by_chain = defaultdict(list)
    for index, edge in enumerate(perimeter.edges):
        if edge.kind == "PHYSICAL" and edge.chain >= 0:
            by_chain[edge.chain].append(index)
    for chain, edge_indices in sorted(by_chain.items()):
        adjacency = defaultdict(list)
        for e in edge_indices:
            a, b = perimeter.edges[e].vertices
            adjacency[a].append(e)
            adjacency[b].append(e)
        ends = sorted(v for v, es in adjacency.items() if perimeter.vertices[v].physical_kind != "REGULAR" or len(es) == 1)
        start = ends[0] if ends else min(adjacency)
        sequence = [start]
        used = set()
        current = start
        while True:
            following = [e for e in adjacency[current] if e not in used]
            if not following:
                break
            e = following[0]
            used.add(e)
            a, b = perimeter.edges[e].vertices
            current = b if a == current else a
            sequence.append(current)
            if current == start or (perimeter.vertices[current].physical_kind != "REGULAR" and current != start):
                break
        sequences.append((chain, sequence, sequence[-1] == start and len(sequence) > 1))
    return sequences


def rounded_runs(perimeter, radius):
    """Fillet arcs with the classifier's rounded-corner reading (DetectRoundedCorners): a
    chain's runs are its maximal collinear chord sequences (a collinear mesh vertex inserted
    by refinement or a spline discretisation does not break an arc); a run shorter than R with
    sub-threshold turns at both ends is an arc chord, a longer run is an arm; a maximal arc
    sequence bounded by two arms is a rounded corner when the tangent distances from the
    virtual corner are both below R and equal within 5 % and the fillet radius
    (a + b) / 2 tan(turn / 2) lies in (0, R). DS-SCT-001: the previous reading split the
    arcs at collinear vertices and accepted one-chord arms (29 rounded runs against the
    classifier's 18 fillets)."""
    runs = []
    turn_threshold = 1.0e-6 * 180.0 / math.pi
    quantum = 1.0e-8 * radius
    for chain, sequence, cycle in chain_vertex_sequences(perimeter):
        points = [perimeter.vertices[v].point for v in sequence]
        n = len(sequence) - (1 if cycle else 0)
        if n < 3:
            continue
        # Turn at every chain vertex (None at the ends of an open chain and at non-regular
        # vertices); chord i joins vertex i to vertex i + 1 (cyclic for a closed chain).
        turns = []
        for i in range(n):
            v = sequence[i]
            if perimeter.vertices[v].physical_kind != "REGULAR" or (not cycle and (i == 0 or i == n - 1)):
                turns.append(None)
                continue
            previous = points[(i - 1) % n]
            following = points[(i + 1) % n]
            d0 = points[i] - previous
            d1 = following - points[i]
            cosine = float(d0 @ d1 / (np.linalg.norm(d0) * np.linalg.norm(d1)))
            turns.append(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))
        chords = n if cycle else n - 1
        # Runs: maximal chord sequences whose interior vertices are collinear. A run is
        # (first chord, chord count); boundary k of the run list is the vertex before run k.
        boundaries = [i for i in range(n) if turns[i] is None or turns[i] > turn_threshold]
        if cycle and not boundaries:
            continue  # a closed polyline without any turn: a curve without arms
        if not cycle:
            boundaries = [0] + [b for b in boundaries if 0 < b < n - 1] + [n - 1]
        run_list = []
        for k in range(len(boundaries) - (0 if cycle else 1)):
            first = boundaries[k]
            last = boundaries[(k + 1) % len(boundaries)]
            count = (last - first) % n if cycle else last - first
            if count == 0:
                count = n
            run_list.append((first, count))
        m = len(run_list)
        if m < 3:
            continue

        def run_points(k):
            first, count = run_list[k]
            return points[first % n], points[(first + count) % n] if cycle else points[first + count]

        def run_length(k):
            a, b = run_points(k)
            return float(np.linalg.norm(b - a))

        def run_tangent(k):
            a, b = run_points(k)
            return (b - a) / np.linalg.norm(b - a)

        def turn_at_boundary(k):
            # Boundary k: the vertex starting run k (a sub-threshold turn between two runs).
            if not cycle and (k == 0 or k >= m):
                return False
            vertex = run_list[k % m][0]
            return turns[vertex] is not None and turns[vertex] > turn_threshold

        is_arc = [run_length(k) < radius - quantum and turn_at_boundary(k) and turn_at_boundary(k + 1) for k in range(m)]
        if all(is_arc):
            continue
        start = is_arc.index(False)
        step = 0
        while step < m:
            k = (start + step) % m
            if not is_arc[k]:
                step += 1
                continue
            count = 0
            while step + count < m and is_arc[(start + step + count) % m]:
                count += 1
            before = (k + m - 1) % m
            after = (k + count) % m
            has_arms = (cycle or (k > 0 and k + count < m)) and not is_arc[before] and not is_arc[after] and turn_at_boundary(k) and turn_at_boundary(k + count)
            arc_vertices = [run_list[(k + c) % m][0] for c in range(1, count)] + [run_list[k][0], run_list[after][0]]
            total = sum(turns[v] for v in arc_vertices if turns[v] is not None)
            record = {"Chain": chain, "Vertices": count + 1, "TotalTurnDegrees": total, "Rounded": False, "Radius": None, "AngleDegrees": None}
            if has_arms:
                ta = run_tangent(before)
                tb = run_tangent(after)
                arm_a_end = run_points(before)[1]
                arm_b_start = run_points(after)[0]
                turn = math.acos(max(-1.0, min(1.0, float(ta @ tb))))
                if math.sin(turn) > 1.0e-9:
                    # Virtual corner X = arm_a_end + a ta = arm_b_start - b tb, solved in the
                    # (ta, y) plane with y = n x ta.
                    w = arm_b_start - arm_a_end
                    y = np.cross(perimeter.process_normal, ta)
                    wx, wy = float(w @ ta), float(w @ y)
                    tbx, tby = float(tb @ ta), float(tb @ y)
                    if abs(tby) > 1.0e-12:
                        b = wy / tby
                        a = wx - b * tbx
                        if a > 0.0 and b > 0.0:
                            fillet = 0.5 * (a + b) / math.tan(0.5 * turn)
                            if a < radius - quantum and b < radius - quantum and abs(a - b) <= 0.05 * max(a, b) and 0.0 < fillet < radius - quantum:
                                record.update({"Rounded": True, "Radius": fillet, "AngleDegrees": 180.0 - math.degrees(turn)})
            runs.append(record)
            step += count
    return runs
