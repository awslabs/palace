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
* a degree-2 perimeter vertex is a CORNER when the turn exceeds the classifier's joint
  noise threshold (metaledge.hpp kCornerTurnToleranceDegrees = 1 deg; USER decision 117(4):
  every sharper joint is a corner unless a fitted arc absorbs it), REGULAR otherwise;
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

# The geometric joint noise threshold (metaledge.hpp kJointNoiseSagittaOverRadius; manifest
# Identification.Conventions.JointNoiseSagittaOverR; USER decision 121 (B), 2026-09-28,
# replacing the 1 deg angular threshold of 117(4)): a two-segment vertex turning by t between
# two straight pieces (collinear mesh edges merged) is a straight continuation of its chain
# when the sagitta (c / 2) tan(t / 4) it implies on the SHORTER adjacent piece c is below this
# multiple of R; otherwise it is a CORNER, which a fitted arc absorbs or which is a corner
# feature. The same constant is the arc rule's mesh-coarseness diagnostic (SAGITTA_OVER_R).
JOINT_NOISE_SAGITTA_OVER_R = 0.05
# The arc rule's joint-turn cap (Conventions.ArcMaxJointTurnDegrees; USER decision 122): a
# joint turning this much or more is never a joint of an arc (regular polygons such as squares
# and hexagons stay corners).
ARC_MAX_JOINT_TURN_DEGREES = 50.0
DIRECTION_QUANTUM = 1.0e-12


def implied_joint_sagitta(turn_radians, shorter_piece):
    """metaledge.hpp ImpliedJointSagitta: (c / 2) tan(t / 4)."""
    return 0.5 * shorter_piece * math.tan(0.25 * turn_radians)


def joint_is_noise(turn_radians, shorter_piece, radius):
    """metaledge.hpp JointIsNoise: the implied sagitta below JOINT_NOISE_SAGITTA_OVER_R x R on
    a 1e-9 relative grid."""
    threshold = JOINT_NOISE_SAGITTA_OVER_R * radius
    grid = 1.0e-9 * threshold
    return round(implied_joint_sagitta(turn_radians, shorter_piece) / grid) < round(threshold / grid)
# Faces whose unit normal deviates from the process normal by more than this are non-planar
# metal (sidewalls, staples, TSV walls); the classifier's kParallelCosineTolerance.
PLANAR_COSINE_TOLERANCE = 1.0e-8
# Metal off an edge's plane within this multiple of R excludes that part of the edge
# (Identification.Conventions CrossLayerReachOverR).
CROSS_LAYER_REACH_OVER_R = 2.0
# Edge kinds the classifier extracts as PHYSICAL segments (FOLD / NONMANIFOLD are its own
# exclusion types).
CLASSIFIER_PHYSICAL_KINDS = ("PHYSICAL", "TRUNCATION", "PORT", "EMBEDDED", "NONPLANAR", "BOX")
# Cuts of the classifier's physical-edge graph: a chain stops there and its end is no feature.
CUT_KINDS = ("TRUNCATION", "PORT")


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


def port_attributes(config):
    """Port boundary attributes of the configuration (LumpedPort elements, WavePort), metal
    attributes excluded: the metal perimeter bordering them is the Port exclusion (decision
    82(5)); metaledge.cpp ExtractMetalEdgeGeometry."""
    boundaries = config.get("Boundaries", {})
    metal = metal_attributes(config)
    result = set()
    for port in boundaries.get("LumpedPort", []):
        for element in port.get("Elements", [port]):
            result |= {int(a) for a in element.get("Attributes", [])}
    for port in boundaries.get("WavePort", []):
        result |= {int(a) for a in port.get("Attributes", [])}
    return {a for a in result if a not in metal}


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


def extract_perimeter(mesh, config, process_normal=None, radius=None):
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
    port_edges = set()
    ports = port_attributes(config)
    interface_edges = defaultdict(set)
    other_mask = ~metal_mask
    for face_index in np.flatnonzero(other_mask):
        tag = int(physical[face_index])
        face = corners[face_index]
        face_key = tuple(np.sort(corners_raw[face_index]))
        is_truncation = tag not in interface_attribute_set and tag not in ports and face_key in exterior
        for i in range(3):
            key = tuple(sorted((int(face[i]), int(face[(i + 1) % 3]))))
            if tag in ports:
                port_edges.add(key)
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
        if edge_kinds[key] == "ONE_SIDED" and key in port_edges:
            kind = "PORT"  # the port is not metal: a cut, the Port exclusion
        elif all(on_box[f] for f in faces) and key not in truncation_edges:
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
    if radius is None:
        raise ValueError("extract_perimeter needs the matching radius: the joint noise rule is geometric")
    classify_vertices(perimeter, radius)
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


def classify_vertices(perimeter, radius):
    """metaledge.cpp ClassifyVertex: quantized direction cosines; a two-segment vertex is a
    collinear continuation (no joint) within the direction quantum, else REGULAR iff the
    geometric joint noise rule holds (joint_is_noise on the shorter adjacent straight piece,
    collinear edges merged), else CORNER."""
    quantized_collinear = _quantize(1.0 - DIRECTION_QUANTUM)
    # physical_kind counts the segments the classifier types PHYSICAL (one-sided, not on the
    # truncation boundary), i.e. the audit's PHYSICAL, EMBEDDED, NONPLANAR and one-sided BOX
    # kinds. Vertex-census rule (SURFACE-RESPONSE-IDENTIFICATION.md (a) 2): a vertex all of
    # whose incident edges are folds / non-manifold edges (the corner of a PEC box where three
    # box faces meet, the base corners of a bump) is a vertex of no one-sided perimeter and has
    # no record on either side; a BOX edge shared by two box faces is such a fold.
    physical_kinds = tuple(k for k in CLASSIFIER_PHYSICAL_KINDS if k not in CUT_KINDS)

    def edges_at(index, physical):
        return [
            e
            for e in perimeter.vertices[index].edges
            if not physical
            or (perimeter.edges[e].kind in physical_kinds and perimeter.edges[e].face_classes == "ONE_SIDED")
        ]

    def unit(from_index, to_index):
        d = perimeter.vertices[to_index].point - perimeter.vertices[from_index].point
        norm = float(np.linalg.norm(d))
        return d / norm, norm

    def piece_length(index, edge, physical):
        """The straight piece leaving the vertex along the edge through collinear two-edge
        vertices (metaledge.cpp PieceLength)."""
        length = 0.0
        current, e = index, edge
        for _ in range(len(perimeter.edges)):
            other = _other_vertex(perimeter, e, current)
            d_in, piece = unit(current, other)
            length += piece
            following = edges_at(other, physical)
            if len(following) != 2 or other == index:
                break
            nxt = following[1] if following[0] == e else following[0]
            d_out, _ = unit(other, _other_vertex(perimeter, nxt, other))
            if _quantize(float(d_in @ d_out)) < quantized_collinear:
                break
            current, e = other, nxt
        return length

    for index, vertex in enumerate(perimeter.vertices):
        for physical in (False, True):
            edges = edges_at(index, physical)
            if not edges:
                kind = None
            elif len(edges) == 1:
                kind = "ENDPOINT"
            elif len(edges) > 2:
                kind = "JUNCTION"
            else:
                directions = [unit(index, _other_vertex(perimeter, e, index))[0] for e in edges]
                dot = float(directions[0] @ directions[1])
                if _quantize(-dot) >= quantized_collinear:
                    kind = "REGULAR"
                else:
                    turn = math.acos(max(-1.0, min(1.0, -dot)))
                    shorter = min(piece_length(index, edges[0], physical), piece_length(index, edges[1], physical))
                    kind = "REGULAR" if joint_is_noise(turn, shorter, radius) else "CORNER"
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


# The arc rule's constants (Identification.Conventions ArcFitToleranceOverR / SagittaOverR,
# USER decision 117(4)): every joint of an arc lies on its circle within the signature
# parameter tolerance and every chord's sagitta is below SAGITTA_OVER_R x R.
SIGNATURE_PARAMETER_TOLERANCE_OVER_R = 1.0e-3
ARC_FIT_TOLERANCE_OVER_R = SIGNATURE_PARAMETER_TOLERANCE_OVER_R
# The mesh-coarseness diagnostic of the arc rule (Conventions.SagittaOverR): an arc whose
# largest chord sagitta reaches this multiple of R is listed in the manifest's
# MeshCoarsenessWarning; it is NOT a membership test (USER decision 122).
SAGITTA_OVER_R = JOINT_NOISE_SAGITTA_OVER_R


def arc_groups(perimeter, radius):
    """The classifier's arc rule (surfaceresponseidentification.cpp DetectArcs; the
    concyclicity form of USER decisions 121 / 122): along every path of PHYSICAL edges through
    vertices with exactly two of them (regular or corner), the joints are the non-collinear
    vertices; from every unconsumed joint the largest range of at least three following joints
    turning the same way, each by less than ARC_MAX_JOINT_TURN_DEGREES, and totalling at most
    180 deg is an arc iff its joint vertices lie on one circle within ARC_FIT_TOLERANCE_OVER_R
    x R — whatever the chord sagitta, which is recorded per arc (MaxChordSagittaOverR; at or
    above SAGITTA_OVER_R it is the mesh-coarseness diagnostic). The circle is the one tangent
    to both arms at the end joints when every joint lies on it (tangent-length radius);
    otherwise, for a radius >= R over at least four joints, the least-squares circle of the
    joints, each arm meeting the circle's tangent at its end joint within the geometric joint
    noise rule (on the shorter of the arm piece and the first chord) or lying on the circle as
    a chord. No piece-length rule enters. A closed path of joints turning one way through
    360 deg on one circle (fit tolerance, every joint below the cap) is one arc. Both
    traversal directions are scanned and the set absorbing more joints wins (then fewer arcs,
    then the smaller invariant serialisation). Radius below R: one rounded corner (any turn);
    otherwise a bend of that exact radius. Corner vertices inside an arc are absorbed (no
    corner of their own). Returns the arcs with their absorbed corner vertices."""
    quantum = 1.0e-8 * radius
    joint_turn_cap = math.radians(ARC_MAX_JOINT_TURN_DEGREES)
    fit_tolerance = ARC_FIT_TOLERANCE_OVER_R * radius
    incident = defaultdict(list)
    for index, edge in enumerate(perimeter.edges):
        if edge.kind == "PHYSICAL" and edge.chain >= 0:
            incident[edge.vertices[0]].append(index)
            incident[edge.vertices[1]].append(index)

    def continues(v):
        return len(incident[v]) == 2 and perimeter.vertices[v].physical_kind in ("REGULAR", "CORNER")

    def other(e, v):
        a, b = perimeter.edges[e].vertices
        return b if a == v else a

    visited = set()
    arcs = []
    normal = np.asarray(perimeter.process_normal, dtype=float)

    def least_squares_circle(points, tangent):
        """Algebraic (Kasa) circle through the points in the plane of the process normal, in
        the frame (u = tangent, v = normal x u) at the first point; None when degenerate."""
        o = points[0]
        u = tangent - normal * float(tangent @ normal)
        if np.linalg.norm(u) <= 1.0e-12:
            return None
        u = u / np.linalg.norm(u)
        v = np.cross(normal, u)
        xy = np.array([[float((p - o) @ u), float((p - o) @ v)] for p in points])
        a = np.column_stack([xy[:, 0], xy[:, 1], np.ones(len(xy))])
        b = -(xy[:, 0] ** 2 + xy[:, 1] ** 2)
        try:
            (d, e, f), *_ = np.linalg.lstsq(a, b, rcond=None)
        except np.linalg.LinAlgError:
            return None
        cx, cy = -0.5 * d, -0.5 * e
        r2 = cx * cx + cy * cy - f
        if not r2 > 0.0:
            return None
        return o + cx * u + cy * v, math.sqrt(r2)

    for seed in range(len(perimeter.edges)):
        if seed in visited or perimeter.edges[seed].kind != "PHYSICAL" or perimeter.edges[seed].chain < 0:
            continue
        # back to the start (or around a loop)
        s, v = seed, perimeter.edges[seed].vertices[0]
        seen = {seed}
        while continues(v):
            previous = incident[v][0] if incident[v][1] == s else incident[v][1]
            if previous in seen:
                break
            seen.add(previous)
            s, v = previous, other(previous, v)
        path_edges, path_vertices = [], []
        while True:
            visited.add(s)
            path_edges.append(s)
            path_vertices.append(v)
            v = other(s, v)
            if not continues(v):
                path_vertices.append(v)
                break
            nxt = incident[v][0] if incident[v][1] == s else incident[v][1]
            if nxt in visited:
                path_vertices.append(v)
                break
            s = nxt
        closed = path_vertices[0] == path_vertices[-1] and len(path_edges) > 2 and continues(path_vertices[0])
        chain = perimeter.edges[path_edges[0]].chain
        found = None
        for direction in (1, -1):
            if direction > 0:
                edges_d, vertices_d = path_edges, path_vertices
            else:
                edges_d, vertices_d = path_edges[::-1], path_vertices[::-1]
            candidate = _scan_arc_path(perimeter, edges_d, vertices_d, closed, normal, radius, quantum, joint_turn_cap,
                                       fit_tolerance, least_squares_circle)
            if found is None or _arc_set_score(candidate, perimeter, path_vertices, path_edges, closed, normal, radius) > \
                    _arc_set_score(found, perimeter, path_vertices, path_edges, closed, normal, radius):
                found = candidate
        for arc in found:
            arc["Chain"] = chain
            rho, total = arc["Radius"], arc["Turn"]
            corner = arc["Tangent"] and rho < radius - quantum
            arcs.append({
                "Joints": arc["Joints"],
                "AbsorbedCorners": [v for v in arc["Joints"] if perimeter.vertices[v].physical_kind == "CORNER"],
                "Radius": rho,
                "RadiusOverR": rho / radius,
                "TurnDegrees": math.degrees(total),
                "AngleDegrees": 180.0 - math.degrees(total),
                "MaxChordSagittaOverR": arc["Sagitta"] / radius,
                "Coarse": arc["Sagitta"] / radius >= SAGITTA_OVER_R,
                "Rounded": corner,
                "Chain": chain,
                "Vertices": len(arc["Joints"]),
            })
    return arcs


def _arc_set_score(arc_list, perimeter, path_vertices, path_edges, closed, normal, radius):
    """The classifier's tie-break between the two traversal directions: more joints absorbed,
    fewer arcs, then the smaller serialisation of (radius, turn, joints, centre distance from
    the path centroid, first joint's distance from the nearer path end) on the signature grid
    (invariant under translation, rotation, mirroring and reversal)."""
    distinct = path_vertices[:-1] if closed else path_vertices
    centroid = np.mean([perimeter.vertices[v].point for v in distinct], axis=0)
    position = [0.0]
    for e in path_edges:
        position.append(position[-1] + perimeter.edges[e].length)
    total = position[-1]
    index_of = {v: k for k, v in reversed(list(enumerate(path_vertices)))}
    keys = []
    for arc in arc_list:
        x_first = position[index_of[arc["Joints"][0]]]
        from_end = 0.0 if closed else min(x_first, total - x_first)
        keys.append("%.6f,%.6f,%d,%.6f,%.6f" % (round(arc["Radius"] / radius, 6), round(math.degrees(arc["Turn"]), 6), len(arc["Joints"]),
                                                round(float(np.linalg.norm(arc["Centre"] - centroid)) / radius, 6), round(from_end / radius, 6)))
    serial = ";".join(sorted(keys))
    absorbed = sum(len(arc["Joints"]) for arc in arc_list)
    # Larger tuple wins: more absorbed, fewer arcs (negated), then the smaller serial (the
    # classifier compares the serial of the other side, so a smaller serial wins).
    return (absorbed, -len(arc_list), _Reversed(serial))


class _Reversed:
    """Orders strings in reverse (the smaller serialisation wins a max-comparison)."""

    def __init__(self, value):
        self.value = value

    def __lt__(self, other):
        return self.value > other.value

    def __gt__(self, other):
        return self.value < other.value

    def __eq__(self, other):
        return self.value == other.value


def _scan_arc_path(perimeter, path_edges, path_vertices, closed, normal, radius, quantum, joint_turn_cap,
                   fit_tolerance, least_squares_circle):
    """One traversal direction of the classifier's greedy arc scan (DetectArcs ScanPath)."""
    n = len(path_edges)
    points = [perimeter.vertices[x].point for x in path_vertices]
    position = [0.0]
    for k in range(n):
        position.append(position[-1] + float(np.linalg.norm(points[k + 1] - points[k])))
    length = position[n]

    def direction(k):
        d = points[k + 1] - points[k]
        return d / np.linalg.norm(d)

    joints = []  # (vertex, path index, in, out, |turn|, sign)
    for k in range(0 if closed else 1, n):
        d_in = direction((k - 1) % n)
        d_out = direction(k)
        dot = float(max(-1.0, min(1.0, d_in @ d_out)))
        if _quantize(dot) >= _quantize(1.0 - DIRECTION_QUANTUM):
            continue
        sign = 1 if float(np.cross(d_in, d_out) @ normal) >= 0.0 else -1
        joints.append((path_vertices[k], k, d_in, d_out, math.acos(dot), sign))
    found = []
    if len(joints) < 2:
        return found
    m = len(joints)
    if closed:
        gaps = []
        for j in range(m):
            gap = position[joints[j][1]] - position[joints[j - 1][1]]
            gaps.append(gap + length if gap <= 0.0 else gap)
        best = int(np.argmax(gaps))
        joints = joints[best:] + joints[:best]

    def point(j):
        return perimeter.vertices[joints[j][0]].point

    def on_circle(indices, centre, rho):
        return all(abs(float(np.linalg.norm(point(j) - centre)) - rho) < fit_tolerance - quantum for j in indices)

    def max_sagitta(indices, rho, cyclic):
        worst = 0.0
        count = len(indices)
        for q in range(count if cyclic else count - 1):
            chord = float(np.linalg.norm(point(indices[q]) - point(indices[(q + 1) % count])))
            if chord >= 2.0 * rho:
                return math.inf
            worst = max(worst, rho - math.sqrt(rho * rho - 0.25 * chord * chord))
        return worst

    def below_cap(j):
        return _quantize(math.cos(joint_turn_cap)) < _quantize(math.cos(joints[j][4]))

    def piece_before(j):
        idx = joints[j][1]
        if j == 0:
            return (position[idx] - position[joints[m - 1][1]] + length) % length if closed else position[idx]
        return position[idx] - position[joints[j - 1][1]]

    def piece_after(j):
        idx = joints[j][1]
        if j + 1 == m:
            return (position[joints[0][1]] - position[idx] + length) % length if closed else length - position[idx]
        return position[joints[j + 1][1]] - position[idx]

    def end_consistent(j, arm_direction, neighbour, arm_edge, centre, rho, first, arm_piece, first_chord):
        """The kink between the arm and the circle's tangent at an end joint of a least-squares
        bend is noise under the geometric rule on the shorter of the arm piece and the first
        chord, or the arm's far vertex lies on the circle (a chord arm: an arc starting at a
        corner on its circle)."""
        at = point(j)
        tangent = np.cross(normal, at - centre)
        tangent = tangent / np.linalg.norm(tangent)
        along = float(tangent @ (point(neighbour) - at))
        if first != (along > 0.0):
            tangent = -tangent
        kink = math.acos(max(-1.0, min(1.0, float(arm_direction @ tangent))))
        if joint_is_noise(kink, min(arm_piece, first_chord), radius):
            return True
        far = perimeter.vertices[_other_vertex(perimeter, arm_edge, joints[j][0])].point
        return abs(float(np.linalg.norm(far - centre)) - rho) < fit_tolerance - quantum

    def try_fit(i, count):
        indices = [(i + j) % m for j in range(count)]
        first, last = joints[i], joints[indices[-1]]
        ta, tb = first[2], last[3]
        Ta, Tb = point(i), point(indices[-1])
        turn = sum(joints[j][4] for j in indices)
        angle = math.acos(max(-1.0, min(1.0, float(ta @ tb))))
        if abs(angle - turn) > 1.0e-6 and abs(2.0 * math.pi - angle - turn) > 1.0e-6:
            return None
        na = first[5] * np.cross(normal, ta)
        na = na / np.linalg.norm(na)
        # the circle tangent to both arms at the end joints
        rho, centre, tangent_circle = 0.0, None, False
        if math.sin(turn) > 1.0e-9 and turn < math.pi - 1.0e-9:
            w = Tb - Ta
            wx, wy = float(w @ ta), float(w @ na)
            tbx, tby = float(tb @ ta), float(tb @ na)
            if abs(tby) > 1.0e-12:
                b = wy / tby
                a = wx - b * tbx
                if a > 0.0 and b > 0.0:
                    rho = 0.5 * (a + b) / math.tan(0.5 * turn)
                    centre = Ta + rho * na
                    tangent_circle = rho > 0.0
        else:
            w = Tb - Ta
            across = float(w @ na)
            if across > 0.0:
                rho = 0.5 * across
                centre = Ta + rho * na
                tangent_circle = True
        if tangent_circle and on_circle(indices, centre, rho):
            sagitta = max_sagitta(indices, rho, False)
            if math.isfinite(sagitta):
                return {"Radius": rho, "Turn": turn, "Centre": centre, "Tangent": True, "Sagitta": sagitta}
            return None
        # least-squares bend (radius >= R) over at least four joints, arms consistent at the ends
        if count < 4:
            return None
        ls = least_squares_circle([point(j) for j in indices], ta)
        if ls is None:
            return None
        centre, rho = ls
        if rho < radius - quantum or not on_circle(indices, centre, rho):
            return None
        if not end_consistent(i, ta, indices[1], path_edges[(first[1] - 1) % n], centre, rho, True, piece_before(i), piece_after(i)):
            return None
        if not end_consistent(indices[-1], tb, indices[-2], path_edges[last[1] % n], centre, rho, False, piece_after(indices[-1]), piece_before(indices[-1])):
            return None
        sagitta = max_sagitta(indices, rho, False)
        if not math.isfinite(sagitta):
            return None
        return {"Radius": rho, "Turn": turn, "Centre": centre, "Tangent": False, "Sagitta": sagitta}

    consumed = [False] * m
    # A closed path of joints turning one way through 360 deg on one circle is one arc (a round
    # pad or hole); the joint-turn cap tells a circle from a polygon (a square or hexagonal
    # hole is corners, an octagon is a circle).
    if closed and m >= 3:
        total = sum(j[4] for j in joints)
        same_sign = all(j[5] == joints[0][5] for j in joints)
        if same_sign and all(below_cap(j) for j in range(m)) and abs(total - 2.0 * math.pi) < 1.0e-6:
            fit = least_squares_circle([point(j) for j in range(m)], joints[0][2])
            if fit is not None:
                centre, rho = fit
                if on_circle(list(range(m)), centre, rho):
                    sagitta = max_sagitta(list(range(m)), rho, True)
                    if math.isfinite(sagitta):
                        found.append({"Joints": [j[0] for j in joints], "Radius": rho, "Turn": 2.0 * math.pi, "Centre": centre, "Tangent": False, "Sagitta": sagitta})
                        consumed = [True] * m
    for i in range(m):
        if consumed[i] or not below_cap(i):
            continue
        best_count, best = 0, None
        turn = joints[i][4]
        for count in range(2, m + 1):
            k = (i + count - 1) % m
            if (not closed and i + count - 1 >= m) or k == i or consumed[k] or joints[k][5] != joints[i][5] or not below_cap(k):
                break
            turn += joints[k][4]
            if turn > math.pi + 1.0e-9:
                break
            if count < 3:
                continue
            fit = try_fit(i, count)
            if fit is not None:
                best, best_count = fit, count
        if best_count == 0:
            continue
        members = [joints[(i + j) % m][0] for j in range(best_count)]
        for j in range(best_count):
            consumed[(i + j) % m] = True
        found.append(dict(best, Joints=members))
    return found


def _other_vertex(perimeter, edge, vertex):
    a, b = perimeter.edges[edge].vertices
    return b if a == vertex else a


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
