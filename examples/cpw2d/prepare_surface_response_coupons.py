#!/usr/bin/env python3

"""Plan, build, qualify, cache, and merge surface-response coupons."""

import argparse
import concurrent.futures
import copy
import csv
import hashlib
import json
import math
import shlex
import shutil
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IDENTIFICATION_TOOLS = ROOT.parent / "surface_response_identification"
if str(IDENTIFICATION_TOOLS) not in sys.path:
    sys.path.insert(0, str(IDENTIFICATION_TOOLS))
import local_pieces as local_pieces_module  # noqa: E402
from signature_library import refuse_duplicate_signatures  # noqa: E402

STRAIGHT_BUILDER = ROOT / "build_surface_process_library.py"
STRAIGHT_QUALIFIER = ROOT / "qualify_surface_response.py"
CLUSTER_MESH = ROOT / "mesh" / "mesh_edge_cluster_coupon.jl"
CLUSTER_GENERATOR = ROOT / "generate_edge_cluster_response.py"
CORNER_ROOT = ROOT.parent / "cpw3d_surface" / "corner_coupon"
CORNER_MESH = CORNER_ROOT / "mesh_corner_coupon.jl"
CORNER_GENERATOR = CORNER_ROOT / "generate_corner_response.py"
CORNER_FINALIZER = CORNER_ROOT / "finalize_corner_response.py"
CORNER_CONVERGENCE = CORNER_ROOT / "run_probe_convergence.py"
# The trace RESOLVABILITY gate of every corner coupon mesh (block (b) family 4 round 2,
# supervisor decision 328 (1), 2026-10-05): every free hat of the trace basis must hold at
# least MINIMUM_ACTIVE_NODES_PER_ORDER_SQUARED x p^2 active boundary nodes of the order-p
# mesh, fail-closed before any solve; the mesher sizes the mesh at the knots, every free
# knot at min(its along-ring gap, its ring's across-ring spacing) (decisions 787 / 791 R-A;
# mesh_corner_coupon.jl --trace-mesh) so that it holds, whatever the basis layout.
CORNER_RESOLVABILITY = CORNER_ROOT / "trace_resolvability.py"
# Knots per ring of the corner trace basis rules (generate_corner_response.RULES: 2 crossings
# + MetalInteriorKnots + FreeKnots; legacy = the recorded MetalRingsOnly rule).
CORNER_TRACE_BASIS_RING_SIZES = {
    "legacy": 8,
    "all-rings-follow-metal": 16,
    # The metal-graded refinement (the concave-MS remedy, decisions 590 / 629 O-9; its DECIDE
    # step): MetalKnotGrading (1/9, 1/3), MetalInteriorKnots 9, RingSize 20, 220 knots. Named
    # here so that a corner build can request it explicitly; NOT the default (the adoption is
    # a separate decision on the decide step's measurements).
    "all-rings-follow-metal-metal-graded": 20,
}
# The CONVEXITY-CONDITIONAL selector (the 633 (1) follow-up; S7-DESIGN 3.2 R2, decision 739):
# a corner build requesting it resolves PER NODE to the metal-graded rule on a CONCAVE node
# (the metal-side chord of decision 590 is a property of the metal-interior knots under the
# MS plane at a concave node) and to the refined default on a CONVEX node (whose metal interior
# is the other side of the box: the 633 (b) convex MS drop cannot recur). Convexity is a
# geometric class of every corner node, so the rule stays generic; the resolved rule is the
# one the coupon spec, the generator and the records carry, so a convex node built under the
# selector is the SAME coupon spec as the default build (its basis files byte-identical by
# construction: test_corner_trace_basis_concave_graded_selector). Not the default.
CORNER_TRACE_BASIS_CONCAVE_GRADED_SELECTOR = "all-rings-follow-metal-concave-graded"
CORNER_TRACE_BASIS_SELECTORS = {
    CORNER_TRACE_BASIS_CONCAVE_GRADED_SELECTOR: {
        "concave": "all-rings-follow-metal-metal-graded",
        "convex": "all-rings-follow-metal",
    },
}
# The default rule of every corner coupon (corner-basis refinement 2026-09-30, USER decision
# 161 (1)). The box trace basis is independent of CornerRadius: the fillet lies inside the
# matching box (its tangency points at r / tan(angle / 2) < R from the apex) and the arms
# are straight where they cross the box, so a rounded corner's basis files are byte-identical
# to the sharp corner's at the same angle (fillet-basis design 2026-10-07 section 1.2, the
# generator test test_corner_basis_files_independent_of_corner_radius). Until decision 511
# a rounded corner was forced onto the legacy MetalRingsOnly rule, which fails the held-out
# self-check on every 90-degree corner (S7's rounded corners: -16.9 / +205 %, decision 310);
# the refined rule now applies to rounded corners too (qualified on the S7 classes and the
# radius sweep of that design) and the runtime's load-time rule check covers them.
CORNER_TRACE_BASIS_DEFAULT = "all-rings-follow-metal"
LIBRARY_COMBINER = CORNER_ROOT / "combine_process_libraries.py"
SPATIAL_ROOT = ROOT.parent / "cpw3d_surface" / "spatial_coupon"
SPATIAL_FIELD_MESH = SPATIAL_ROOT / "mesh_spatial_coupon.jl"
SPATIAL_SWEPT_MESH = SPATIAL_ROOT / "mesh_spatial_coupon_swept.jl"
# Retain the historical name for plan manifests and callers that inspect it. Actual
# execution resolves the explicitly selected implementation with spatial_mesh_path().
SPATIAL_MESH = SPATIAL_FIELD_MESH
SPATIAL_GENERATOR = SPATIAL_ROOT / "generate_spatial_response.py"

PAIRED_TOPOLOGIES = {
    "SameConductorGap": "same-conductor-gaps",
    "DifferentConductorGap": "different-conductor-gaps",
    "SameConductorStrip": "same-conductor-strips",
}
CORNER_TOPOLOGIES = {"ConvexCorner", "ConcaveCorner"}
SPATIAL_TOPOLOGIES = {
    "SpatialEdgeCluster",
    "Endpoint",
    "Junction",
}


def load_json(path):
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid JSON in {path}: {error}") from error


def write_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n")


def run(command, check=True, cwd=None):
    command = [str(value) for value in command]
    location = f" (cwd={cwd})" if cwd else ""
    print("+ " + shlex.join(command) + location, flush=True)
    return subprocess.run(command, check=check, cwd=cwd).returncode


def slug(value):
    return (
        f"{value:.12g}".replace("-", "m").replace(".", "p")
        if isinstance(value, float)
        else str(value).lower()
    )


def coupon_signature(requirement):
    geometry = copy.deepcopy(requirement.get("Geometry") or {})
    if "PlanViewBoundary" in geometry:
        geometry.pop("PlanViewFacets", None)
    return {
        "Topology": requirement["Topology"],
        "Geometry": geometry,
        "BoundaryCondition": requirement.get("BoundaryCondition", {}),
        "Interfaces": requirement.get("Interfaces", []),
    }


def coupon_id(requirement):
    topology = requirement["Topology"]
    geometry = requirement.get("Geometry") or {}
    parts = [topology.lower()]
    for key in (
        "Separation",
        "AngleDegrees",
        "CornerRadius",
        "ArmCount",
        "EdgeCount",
    ):
        if key in geometry:
            parts.append(f"{key.lower()}-{slug(geometry[key])}")
    digest = hashlib.sha256(
        json.dumps(
            coupon_signature(requirement),
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()[:12]
    return "_".join(parts + [digest])


def verified_boundary_condition(condition):
    if not isinstance(condition, dict):
        raise ValueError(
            "boundary-law metadata must be an object with dimensional named parameters"
        )
    boundary_type = condition.get("Type")
    allowed = {
        "PEC": {"Type"},
        "Conductivity": {
            "Type",
            "Conductivity",
            "Permeability",
            "Thickness",
            "External",
        },
        "Impedance": {"Type", "Rs", "Ls", "Cs"},
        "RationalImpedance": {"Type", "Numerator", "Denominator"},
    }
    if boundary_type not in allowed:
        raise ValueError(f"unknown boundary-law type {boundary_type!r}")
    unknown = set(condition) - allowed[boundary_type]
    if unknown:
        raise ValueError(
            f"{boundary_type} boundary law has unsupported metadata "
            + ", ".join(sorted(unknown))
        )

    def finite(name, default=None):
        if name not in condition:
            if default is None:
                raise ValueError(f"{boundary_type} boundary law is missing {name}")
            return default
        value = condition[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{boundary_type} boundary-law {name} is not numeric")
        value = float(value)
        if not math.isfinite(value):
            raise ValueError(f"{boundary_type} boundary-law {name} is not finite")
        return value

    if boundary_type == "Conductivity":
        if finite("Conductivity") <= 0.0:
            raise ValueError("Conductivity boundary law must have positive Conductivity")
        if finite("Permeability", 1.0) <= 0.0:
            raise ValueError("Conductivity boundary law must have positive Permeability")
        if finite("Thickness", 0.0) < 0.0:
            raise ValueError("Conductivity boundary law cannot have negative Thickness")
        if "External" in condition and not isinstance(condition["External"], bool):
            raise ValueError("Conductivity boundary-law External must be boolean")
    elif boundary_type == "Impedance":
        values = [finite(name, 0.0) for name in ("Rs", "Ls", "Cs")]
        if not any(value != 0.0 for value in values):
            raise ValueError("Impedance boundary law must have a nonzero Rs, Ls, or Cs")
    elif boundary_type == "RationalImpedance":
        for name in ("Numerator", "Denominator"):
            coefficients = condition.get(name)
            if not isinstance(coefficients, list) or not coefficients:
                raise ValueError(
                    f"RationalImpedance boundary-law {name} must be a nonempty array"
                )
            if any(
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(float(value))
                for value in coefficients
            ):
                raise ValueError(
                    f"RationalImpedance boundary-law {name} must be finite and numeric"
                )
            if not any(float(value) != 0.0 for value in coefficients):
                raise ValueError(
                    f"RationalImpedance boundary-law {name} must be nonzero"
                )
    return copy.deepcopy(condition)


def coupon_boundary_conditions(coupon):
    conditions = [coupon.get("BoundaryCondition")]
    geometry = coupon.get("Geometry", {})
    for entries in (geometry.get("Edges", []), geometry.get("Arms", [])):
        conditions.extend(
            entry.get("BoundaryCondition", coupon.get("BoundaryCondition"))
            for entry in entries
            if isinstance(entry, dict)
        )
    return conditions


def require_verified_boundary_conditions(coupon):
    for condition in coupon_boundary_conditions(coupon):
        verified_boundary_condition(condition)


def uses_finite_impedance(coupon):
    return any(
        isinstance(condition, dict) and condition.get("Type") != "PEC"
        for condition in coupon_boundary_conditions(coupon)
    )


def preparation(requirement, spatial_mesh=SPATIAL_MESH):
    topology = requirement["Topology"]
    geometry = requirement.get("Geometry") or {}
    try:
        require_verified_boundary_conditions(requirement)
    except ValueError as error:
        return {
            "Method": "Unsupported",
            "Reason": f"Unverified boundary-law metadata: {error}",
        }
    if topology == "IsolatedEdge" or topology in PAIRED_TOPOLOGIES:
        return {
            "Method": "StraightEdgeBuilder",
            "Driver": str(STRAIGHT_BUILDER.relative_to(ROOT.parent.parent)),
        }
    if topology == "ParallelEdgeCluster":
        return {
            "Method": "ParallelClusterCoupon",
            "MeshGenerator": str(CLUSTER_MESH.relative_to(ROOT.parent.parent)),
            "ResponseGenerator": str(
                CLUSTER_GENERATOR.relative_to(ROOT.parent.parent)
            ),
        }
    if topology in CORNER_TOPOLOGIES:
        angle = float(geometry.get("AngleDegrees", math.nan))
        # (0, 180]: 180 with a zero corner radius is the corner family's straight anchor.
        if math.isfinite(angle) and 0.0 < angle <= 180.0 and (angle < 180.0 or float(geometry.get("CornerRadius", 0.0)) == 0.0):
            return {
                "Method": "CornerCoupon",
                "MeshGenerator": str(CORNER_MESH.relative_to(ROOT.parent.parent)),
                "ResponseGenerator": str(
                    CORNER_GENERATOR.relative_to(ROOT.parent.parent)
                ),
            }
        return {
            "Method": "Unsupported",
            "Reason": "Corner angles must lie in (0, 180] degrees (180 only for the sharp straight anchor)",
        }
    if topology in SPATIAL_TOPOLOGIES:
        return {
            "Method": "SpatialCoupon",
            "MeshGenerator": str(
                spatial_mesh.relative_to(ROOT.parent.parent)
            ),
            "ResponseGenerator": str(
                SPATIAL_GENERATOR.relative_to(ROOT.parent.parent)
            ),
        }
    return {"Method": "Unsupported", "Reason": "Unknown coupon topology"}


def plan_from_manifest(
    manifest_path,
    manifest,
    library_path,
    library,
    include_matched,
    spatial_mesh=SPATIAL_MESH,
    one_reference=False,
):
    """`one_reference`: every coupon planned in the one-reference form (the opt-out of the piece-form default,
    plan_local_pieces_form / LOCAL_PIECES_PLANNER_RULE)."""
    coupons = {}
    signatures = {}
    matching_radius = float(manifest["Library"]["MatchingRadius"])
    try:
        parameters = process_parameters(library)
    except ValueError:
        parameters = None   # no version-3 Fabrication metadata: a version-2 spatial coupon's pieces are not derivable
                            # (plan_local_pieces_form: Unsupported by name, decision 1145)
    for requirement in manifest["Requirements"]:
        if not include_matched and requirement.get("Status") != "Missing":
            continue
        # A version-2 record without geometry parameters (an isolated edge) carries null.
        requirement["Geometry"] = requirement.get("Geometry") or {}
        geometry = requirement["Geometry"]
        facets = geometry.get("PlanViewFacets")
        boundary = geometry.get("PlanViewBoundary")
        if facets and boundary:
            process_axis = (
                1 if requirement["Topology"] == "SpatialEdgeCluster" else 2
            )
            if unclassified_plan_view_boundary(boundary) != (
                canonical_plan_view_boundary(
                    facets, matching_radius, process_axis
                )
            ):
                raise ValueError(
                    "Palace PlanViewBoundary does not match its exported facets"
                )
        identifier = coupon_id(requirement)
        coupon = {
            "Id": identifier,
            "Topology": requirement["Topology"],
            "Geometry": requirement.get("Geometry") or {},
            "Interfaces": requirement.get("Interfaces", []),
            "BoundaryCondition": requirement.get("BoundaryCondition", {}),
            "DeviceOccurrences": requirement.get("Count", 0),
            "DeviceEdgeLength": requirement.get("TotalEdgeLength", 0.0),
            "CoverageStatus": requirement.get("Status"),
            "Preparation": preparation(requirement, spatial_mesh),
        }
        # Version-2 records (SURFACE-RESPONSE-IDENTIFICATION.md (d)): the canonical
        # Signature is the library model's key; the group statistics are recorded. A
        # mirror-formed cluster (round-3 DESIGN 4.4, decision 557: impl-B5 CONTRACT.md) carries
        # MirrorFormed + MirrorFormedContract, consumed by device_coupons. KeyText (decision
        # 605 (2) F-4 (c); absent on record-era manifests) is the exact text palace hashed
        # into Hash, carried as an opaque string for the key checks downstream.
        for key in ("Signature", "Hash", "KeyText", "Instances", "DistinctSignatures", "ParameterSpread", "ExactParameters",
                    "NearKeys", "SpanCapAllowance", "MirrorFormed", "MirrorFormedContract"):
            if key in requirement:
                coupon[key] = copy.deepcopy(requirement[key])
        if uses_finite_impedance(coupon):
            coupon["Preparation"]["BoundaryLawQualification"] = "Missing"
        if "SelectedModels" in requirement:
            coupon["SelectedModels"] = requirement["SelectedModels"]
        if "Reason" in requirement:
            coupon["CoverageReason"] = requirement["Reason"]
        signature = coupon_signature(requirement)
        if identifier in coupons:
            if signatures[identifier] != signature:
                raise ValueError(f"Coupon identifier collision for {identifier}")
            existing = coupons[identifier]
            for field in (
                "Topology",
                "Interfaces",
                "BoundaryCondition",
                "CoverageStatus",
                "Preparation",
                "SelectedModels",
                "CoverageReason",
                "Signature",
            ):
                if existing.get(field) != coupon.get(field):
                    raise ValueError(
                        f"Inconsistent repeated coupon {identifier}: {field}"
                    )
            existing["DeviceOccurrences"] += coupon["DeviceOccurrences"]
            existing["DeviceEdgeLength"] += coupon["DeviceEdgeLength"]
        else:
            coupons[identifier] = coupon
            signatures[identifier] = signature
    coupons = sorted(coupons.values(), key=lambda coupon: coupon["Id"])
    for coupon in coupons:
        plan_local_pieces_form(coupon, matching_radius, parameters, one_reference)
    return {
        "Version": 2,
        "SourceManifest": str(manifest_path),
        "ProcessLibrary": str(library_path),
        "Fabrication": library.get("Fabrication"),
        "MatchingRadius": manifest["Library"]["MatchingRadius"],
        "Coupons": coupons,
        "LocalPiecesForm": {"Rule": LOCAL_PIECES_PLANNER_RULE, "OneReference": bool(one_reference),
                            "Coupons": [coupon["Id"] for coupon in coupons if coupon_local_pieces(coupon)]},
        "Summary": {
            "CouponCount": len(coupons),
            "LocalPiecesForm": sum(coupon_local_pieces(coupon) for coupon in coupons),
            "AutomaticallyPreparatable": sum(
                coupon["Preparation"]["Method"] != "Unsupported"
                for coupon in coupons
            ),
            "Unsupported": sum(
                coupon["Preparation"]["Method"] == "Unsupported"
                for coupon in coupons
            ),
        },
    }


def metadata_value(fabrication, names, required=True, default=None):
    for name in names:
        if name in fabrication:
            return fabrication[name]
    if required:
        raise ValueError(
            "Fabrication metadata is missing " + " or ".join(names)
        )
    return default


def process_parameters(library):
    fabrication = library.get("Fabrication")
    if not isinstance(fabrication, dict):
        raise ValueError(
            "Automatic coupon execution requires version-3 Fabrication metadata"
        )
    if fabrication.get("LengthUnit", "um") != "um":
        raise ValueError("Automatic coupon generation currently requires microns")
    layers = fabrication.get("InterfaceLayers")
    if not isinstance(layers, dict) or any(
        interface not in layers for interface in ("SA", "MS", "MA")
    ):
        raise ValueError("Fabrication metadata must define SA, MS, and MA layers")
    result = {
        "metal_thickness": float(
            metadata_value(fabrication, ("MetalThickness",))
        ),
        "overetch": float(
            metadata_value(fabrication, ("OveretchDepth", "Overetch"))
        ),
        "sidewall_angle": float(
            metadata_value(
                fabrication, ("SidewallAngleDegrees", "SidewallAngle")
            )
        ),
        "top_radius": float(
            metadata_value(
                fabrication,
                ("TopRoundingRadius", "TopRadius", "TopRounding"),
            )
        ),
        "bottom_radius": float(
            metadata_value(
                fabrication,
                ("BottomRoundingRadius", "BottomRadius", "TrenchRounding"),
            )
        ),
        "substrate_permittivity": float(
            metadata_value(fabrication, ("SubstratePermittivity",))
        ),
    }
    for interface in ("SA", "MS", "MA"):
        result[f"{interface.lower()}_thickness"] = float(
            metadata_value(layers[interface], ("Thickness",))
        )
        result[f"{interface.lower()}_permittivity"] = float(
            metadata_value(layers[interface], ("Permittivity",))
        )
    if any(not math.isfinite(value) for value in result.values()):
        raise ValueError("Fabrication metadata contains a nonfinite value")
    return result


def fingerprint(data):
    return hashlib.sha256(
        json.dumps(data, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def tool_fingerprint(paths):
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def julia_command(args):
    command = [args.julia]
    if args.julia_project:
        command.append(f"--project={args.julia_project.expanduser().resolve()}")
    return command


def palace_command(args, config):
    command = [args.palace]
    command.extend(["--serial"] if args.ranks == 1 else ["-np", args.ranks])
    command.append(config)
    return command


def material_options(parameters):
    options = ["--substrate-permittivity", parameters["substrate_permittivity"]]
    for interface in ("sa", "ms", "ma"):
        options.extend(
            [
                f"--{interface}-thickness",
                parameters[f"{interface}_thickness"],
                f"--{interface}-permittivity",
                parameters[f"{interface}_permittivity"],
            ]
        )
    return options


def process_resolution(parameters, fine_size, minimum_elements):
    if fine_size <= 0.0 or minimum_elements <= 0.0:
        raise ValueError("Process resolution parameters must be positive")
    features = {
        "MetalThickness": parameters["metal_thickness"],
        "OveretchDepth": parameters["overetch"],
    }
    features = {
        name: float(value)
        for name, value in features.items()
        if float(value) > 0.0
    }
    if not features:
        raise ValueError("Fabrication metadata has no positive geometric feature")
    limiting_name, limiting_size = min(features.items(), key=lambda item: item[1])
    elements = limiting_size / fine_size
    report = {
        "FineSize": fine_size,
        "MinimumElementsPerFeature": minimum_elements,
        "LimitingFeature": limiting_name,
        "LimitingFeatureSize": limiting_size,
        "ElementsAcrossLimitingFeature": elements,
        "RoundingMinimumCirclePoints": 24,
        "Passed": elements >= minimum_elements,
    }
    if not report["Passed"]:
        raise ValueError(
            f"Coupon fine size {fine_size:g} resolves {limiting_name} "
            f"({limiting_size:g}) with only {elements:.3g} elements; "
            f"at least {minimum_elements:g} are required"
        )
    return report


def canonical_plan_view_boundary(facets, matching_radius, process_axis=1):
    if not facets:
        return []
    if process_axis not in (0, 1, 2):
        raise ValueError("Plan-view process axis is invalid")
    tolerance = 1.0e-9 * matching_radius
    # Block (b) DESIGN A3 (1): two consecutive boundary sides turning by at most this angle
    # (radians) are one side (generate_spatial_response.JUNCTION_TANGENT_ANGLE, the mesher's
    # constant); exactly collinear sides merge as before.
    junction_tangent_angle = 1.0e-4

    def subtract(first, second):
        return tuple(a - b for a, b in zip(first, second))

    def dot(first, second):
        return sum(a * b for a, b in zip(first, second))

    def cross(first, second):
        return (
            first[1] * second[2] - first[2] * second[1],
            first[2] * second[0] - first[0] * second[2],
            first[0] * second[1] - first[1] * second[0],
        )

    def near_collinear(first, second):
        # Both directions leave the vertex: a smooth joint has them anti-parallel.
        turn = cross(first, second)
        if turn == (0, 0, 0):
            return True
        return dot(first, second) < 0 and math.sqrt(
            dot(turn, turn)
        ) <= junction_tangent_angle * math.sqrt(dot(first, first)) * math.sqrt(
            dot(second, second)
        )

    def quantize(point):
        def llround(value):
            scaled = value / tolerance
            return (
                math.floor(scaled + 0.5)
                if scaled >= 0.0
                else math.ceil(scaled - 0.5)
            )

        return tuple(llround(value) for value in point)

    groups = {}
    for facet in facets:
        conductor = int(facet["Conductor"])
        points = facet["Points"]
        if conductor <= 0 or len(points) < 3:
            raise ValueError("Plan-view facet is invalid")
        ring = []
        for point in points:
            if len(point) != 3 or not all(
                math.isfinite(float(value)) for value in point
            ):
                raise ValueError("Plan-view facet is invalid")
            key = quantize(tuple(float(value) for value in point))
            if not ring or ring[-1] != key:
                ring.append(key)
        if len(ring) > 1 and ring[0] == ring[-1]:
            ring.pop()
        if len(ring) < 3:
            continue
        plane = ring[0][process_axis]
        if any(point[process_axis] != plane for point in ring):
            raise ValueError("Plan-view facet is not on one process plane")
        groups.setdefault((conductor, plane), []).append(ring)

    result = []
    for (conductor, _), group_polygons in sorted(groups.items()):
        polygons = []
        seen = set()
        for ring in group_polygons:
            sequences = []
            for reverse in (False, True):
                ordered = list(reversed(ring)) if reverse else ring
                sequences.extend(
                    tuple(ordered[index:] + ordered[:index])
                    for index in range(len(ordered))
                )
            key = min(sequences)
            if key not in seen:
                seen.add(key)
                polygons.append(ring)

        vertices = [point for polygon in polygons for point in polygon]
        counts = {}
        for polygon in polygons:
            for begin, end in zip(polygon, polygon[1:] + polygon[:1]):
                direction = subtract(end, begin)
                length_squared = dot(direction, direction)
                if length_squared == 0:
                    continue
                split = [begin, end]
                for point in vertices:
                    if point == begin or point == end:
                        continue
                    offset = subtract(point, begin)
                    if cross(direction, offset) != (0, 0, 0):
                        continue
                    coordinate = dot(offset, direction)
                    if 0 < coordinate < length_squared:
                        split.append(point)
                split = sorted(
                    set(split),
                    key=lambda point: dot(subtract(point, begin), direction),
                )
                for first, second in zip(split, split[1:]):
                    edge = tuple(sorted((first, second)))
                    counts[edge] = counts.get(edge, 0) + 1

        if any(count > 2 for count in counts.values()):
            raise ValueError("Plan-view facets form a nonmanifold surface")
        boundary = {edge for edge, count in counts.items() if count % 2 == 1}
        boundary_degree = {}
        for first, second in boundary:
            boundary_degree[first] = boundary_degree.get(first, 0) + 1
            boundary_degree[second] = boundary_degree.get(second, 0) + 1
        if any(degree % 2 for degree in boundary_degree.values()):
            raise ValueError("Plan-view facet union has an open boundary")
        while True:
            adjacency = {}
            for first, second in boundary:
                adjacency.setdefault(first, set()).add(second)
                adjacency.setdefault(second, set()).add(first)
            merge = None
            for vertex, neighbors in adjacency.items():
                if len(neighbors) != 2:
                    continue
                first, second = tuple(neighbors)
                if near_collinear(
                    subtract(first, vertex),
                    subtract(second, vertex),
                ):
                    merge = (vertex, first, second)
                    break
            if merge is None:
                break
            vertex, first, second = merge
            boundary.remove(tuple(sorted((first, vertex))))
            boundary.remove(tuple(sorted((vertex, second))))
            boundary.add(tuple(sorted((first, second))))
        boundary = sorted(boundary)
        if not boundary:
            raise ValueError("Plan-view facets have no union boundary")
        result.append(
            {
                "Conductor": conductor,
                "Segments": [[list(first), list(second)] for first, second in boundary],
            }
        )
    result.sort(
        key=lambda entry: json.dumps(
            entry, sort_keys=True, separators=(",", ":")
        )
    )
    return result


def unclassified_plan_view_boundary(boundary):
    return [
        {
            "Conductor": component["Conductor"],
            "Segments": component["Segments"],
        }
        for component in boundary
    ]


def model_matches_coupon(model, coupon):
    if model.get("Topology") != coupon["Topology"]:
        return False
    geometry = coupon.get("Geometry") or {}

    def close(a, b):
        return math.isclose(float(a), float(b), rel_tol=1.0e-9, abs_tol=1.0e-9)

    if coupon["Topology"] in PAIRED_TOPOLOGIES:
        return close(model.get("Separation", math.nan), geometry["Separation"])
    if coupon["Topology"] == "ParallelEdgeCluster":
        model_edges = model.get("Edges") or []
        coupon_edges = geometry.get("Edges") or []
        if len(model_edges) != len(coupon_edges):
            return False
        for first, second in zip(model_edges, coupon_edges):
            offset = second["Offset"][0] if isinstance(second["Offset"], list) else second["Offset"]
            direction = second["GapDirection"][0] if isinstance(second["GapDirection"], list) else second["GapDirection"]
            if not close(first["Offset"], offset) or int(round(float(first["GapDirection"]))) != int(round(float(direction))):
                return False
        return True
    if coupon["Topology"] in CORNER_TOPOLOGIES:
        model_connectivity = (model.get("TraceBasis") or {}).get("ConnectivityAngleDegrees")
        coupon_connectivity = geometry.get("ConnectivityAngleDegrees")
        if (model_connectivity is None) != (coupon_connectivity is None) or (
            model_connectivity is not None
            and not close(model_connectivity, coupon_connectivity)
        ):
            return False
        return close(model.get("Angle", math.nan), geometry["AngleDegrees"]) and close(
            model.get("CornerRadius", 0.0), geometry.get("CornerRadius", 0.0)
        )
    if coupon["Topology"] == "Junction":
        return list(model.get("ArmAngles", [])) == list(geometry.get("ArmAnglesDegrees", []))
    if coupon["Topology"] == "SpatialEdgeCluster":
        # OI-PLANNER-SKIPPED-SIGNATURE-STAMP (decisions 1117 / 1137; RG-7: the merged library re-stamped v3-j's
        # spatialedgecluster_edgecount-10_13c7cc02fd9b with the SKIPPED stub cluster's 6-edge Signature and the binary
        # refused the library): a spatial model never matches a coupon by topology alone. A model carrying a Signature
        # (a version-2 model, the matcher's exact key) matches the coupon of the SAME Signature only; a model without
        # one matches by its Edges (count, Point, GapDirection, Conductor) against a coupon that carries Edges.
        if "Signature" in model:
            return "Signature" in coupon and json.dumps(model["Signature"], sort_keys=True) == json.dumps(
                coupon["Signature"], sort_keys=True)
        model_edges = model.get("Edges") or []
        coupon_edges = geometry.get("Edges") or []
        if not coupon_edges or len(model_edges) != len(coupon_edges):
            return False
        for first, second in zip(model_edges, coupon_edges):
            if int(first.get("Conductor", 1)) != int(second.get("Conductor", 1)):
                return False
            for key in ("Point", "GapDirection"):
                if len(first.get(key) or []) != 3 or len(second.get(key) or []) != 3 or not all(
                        close(a, b) for a, b in zip(first[key], second[key])):
                    return False
        return True
    return True


def built_coupons_for_stamping(built, skipped, failures):
    """The coupons that may stamp the merged library (stamp_library_signatures): the ones built by this execution, with
    every Skipped / failed Id refused by name should one appear (a contract check: a skipped requirement never stamps a
    library model; OI-PLANNER-SKIPPED-SIGNATURE-STAMP, decisions 1117 / 1137)."""
    excluded = {identifier for entry in list(skipped) + list(failures) for identifier in entry.get("Ids", [])}
    offending = sorted(coupon["Id"] for coupon in built if coupon["Id"] in excluded)
    if offending:
        raise ValueError(f"skipped / failed coupon(s) {offending} were listed as built: a skipped requirement never stamps "
                         "a library model (OI-PLANNER-SKIPPED-SIGNATURE-STAMP)")
    return list(built)


def stamp_library_signatures(path, coupons):
    """Every model of a written library that corresponds to exactly one version-2 coupon
    carries that coupon's canonical Signature (the matcher's key: exact for a
    SpatialEdgeCluster, within the signature tolerance otherwise) and the coupon's group
    statistics; a model matching several signature coupons is a contract error (fail
    closed), a model matching none (a version-1 coupon, a base library) is left as written.
    `coupons` = the coupons BUILT by this execution (execute: never a Skipped, Unsupported or
    failed coupon; OI-PLANNER-SKIPPED-SIGNATURE-STAMP, decisions 1117 / 1137: a skipped
    requirement never stamps a library model).
    Two models carrying one Signature after the stamping (a base model and this run's model of
    one key; OI-LIBRARY-DUPLICATE-SIGNATURE, decisions 1103 / 1126 (5)) are refused by name
    BEFORE the library is written: a library holds at most one model per signature."""
    library = load_json(path)
    stamped = 0
    for model in library.get("Models", []):
        matches = [
            coupon for coupon in coupons if "Signature" in coupon and model_matches_coupon(model, coupon)
        ]
        signatures = {json.dumps(coupon["Signature"], sort_keys=True) for coupon in matches}
        if len(signatures) > 1:
            # model_matches_coupon compares the geometry only (the boundary law is stamped
            # per coupon by stamp_library_boundary_conditions): name each coupon's law so
            # that two coupons of one geometry under different laws are read as the cause.
            raise ValueError(
                f"model {model.get('Name')} matches {len(signatures)} version-2 coupon signatures: "
                + ", ".join(
                    f"{coupon['Id']} (BoundaryCondition {json.dumps(coupon.get('BoundaryCondition') or {'Type': 'PEC'}, sort_keys=True)})"
                    for coupon in sorted(matches, key=lambda c: c["Id"])
                )
            )
        if not matches:
            continue
        coupon = matches[0]
        model["Signature"] = copy.deepcopy(coupon["Signature"])
        for key in ("Instances", "DistinctSignatures", "ParameterSpread"):
            if key in coupon:
                model[key] = coupon[key]
        stamped += 1
    refuse_duplicate_signatures(library.get("Models", []), library=str(path))
    write_json(path, library)
    return stamped


def stamp_library_boundary_conditions(path, coupons):
    library = load_json(path)
    models = []
    for model in library.get("Models", []):
        matches = [
            coupon for coupon in coupons if model_matches_coupon(model, coupon)
        ]
        if not matches:
            models.append(model)
            continue
        unique = {}
        for coupon in matches:
            require_verified_boundary_conditions(coupon)
            signature = json.dumps(
                coupon_boundary_conditions(coupon),
                sort_keys=True,
                separators=(",", ":"),
            )
            unique.setdefault(signature, coupon)
        for coupon in unique.values():
            stamped = copy.deepcopy(model)
            if coupon["Topology"] == "SpatialEdgeCluster":
                source_edges = coupon["Geometry"].get("Edges", [])
                target_edges = stamped.get("Edges", [])
                if len(source_edges) != len(target_edges):
                    raise ValueError(
                        f"{coupon['Id']} does not match its generated spatial model"
                    )
                for source, target in zip(source_edges, target_edges):
                    target["BoundaryCondition"] = verified_boundary_condition(
                        source.get(
                            "BoundaryCondition", coupon["BoundaryCondition"]
                        )
                    )
                stamped.pop("BoundaryCondition", None)
            else:
                stamped["BoundaryCondition"] = verified_boundary_condition(
                    coupon["BoundaryCondition"]
                )
            if uses_finite_impedance(coupon):
                stamped["BoundaryLawQualification"] = {
                    "Version": 1,
                    "Status": "Unqualified",
                    "Calibration": "QuasiElectrostatic",
                    "FrequencyUniversal": False,
                }
            else:
                stamped.pop("BoundaryLawQualification", None)
            if len(unique) > 1:
                digest = hashlib.sha256(
                    json.dumps(
                        coupon_boundary_conditions(coupon),
                        sort_keys=True,
                        separators=(",", ":"),
                    ).encode()
                ).hexdigest()[:8]
                stamped["Name"] = f"{stamped['Name']}-bc-{digest}"
            models.append(stamped)
    library["Models"] = models
    write_json(path, library)


def stamp_geometry_coverage_only(path):
    library = load_json(path)
    for model in library.get("Models", []):
        model["BoundaryLawQualification"] = {
            "Version": 1,
            "Status": "Unqualified",
            "Calibration": "GeometryCoverageOnly",
            "FrequencyUniversal": False,
        }
    write_json(path, library)


def pad_missing_surface_response_matrices(root, interface_count, matching_radius):
    """Add exact-zero blocks for interfaces absent from a generated coupon mesh."""
    for kind in ("thin", "fabricated"):
        postpro = root / "postpro" / f"spatial_{kind}"
        domain_path = postpro / "domain-response-matrix.csv"
        surface_path = postpro / "surface-response-matrix.csv"
        with domain_path.open(newline="") as stream:
            domain_rows = list(csv.DictReader(stream, skipinitialspace=True))
        basis_size = max(
            max(int(float(row["basis_i"])), int(float(row["basis_j"])))
            for row in domain_rows
        )
        rows = []
        fieldnames = [
            "interface",
            "edge",
            "R (m)",
            "basis_i",
            "basis_j",
            "Q_ij (J)",
            "Q_ij normal (J)",
            "Q_ij tangential (J)",
            "Q_total_ij (J)",
            "Q_total_ij normal (J)",
            "Q_total_ij tangential (J)",
        ]
        if surface_path.is_file():
            with surface_path.open(newline="") as stream:
                reader = csv.DictReader(stream, skipinitialspace=True)
                fieldnames = [name.strip() for name in reader.fieldnames]
                rows = [
                    {name.strip(): value for name, value in row.items()}
                    for row in reader
                ]
        present = {int(float(row["interface"])) for row in rows}
        for interface in range(1, interface_count + 1):
            if interface in present:
                continue
            for i in range(1, basis_size + 1):
                for j in range(i, basis_size + 1):
                    rows.append(
                        {
                            "interface": interface,
                            "edge": 1,
                            "R (m)": matching_radius * 1.0e-6,
                            "basis_i": i,
                            "basis_j": j,
                            "Q_ij (J)": 0.0,
                            "Q_ij normal (J)": 0.0,
                            "Q_ij tangential (J)": 0.0,
                            "Q_total_ij (J)": 0.0,
                            "Q_total_ij normal (J)": 0.0,
                            "Q_total_ij tangential (J)": 0.0,
                        }
                    )
        with surface_path.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)


def normalize_parallel_edges(coupon, matching_radius):
    raw_edges = coupon.get("Geometry", {}).get("Edges")
    if not isinstance(raw_edges, list) or len(raw_edges) < 3:
        raise ValueError(
            f"{coupon['Id']} has no complete ParallelEdgeCluster edge signature"
        )
    tolerance = max(1.0e-10 * matching_radius, 1.0e-12)
    edges = []
    for index, raw in enumerate(raw_edges):
        offset = raw.get("Offset")
        direction = raw.get("GapDirection")
        if isinstance(offset, list):
            if len(offset) != 2 or abs(float(offset[1])) > tolerance:
                raise ValueError(
                    f"{coupon['Id']} edge {index + 1} is not coplanar"
                )
            offset = offset[0]
        if isinstance(direction, list):
            if (
                len(direction) != 2
                or abs(float(direction[1])) > 1.0e-10
                or not math.isclose(abs(float(direction[0])), 1.0, abs_tol=1.0e-10)
            ):
                raise ValueError(
                    f"{coupon['Id']} edge {index + 1} is not a canonical "
                    "parallel cross-section edge"
                )
            direction = direction[0]
        edge = {
            "Offset": float(offset),
            "GapDirection": int(round(float(direction))),
            "Conductor": int(raw["Conductor"]),
        }
        edges.append(edge)
    if edges[0]["Conductor"] == 0:
        for edge in edges:
            edge["Conductor"] += 1
    labels = []
    for edge in edges:
        conductor = edge["Conductor"]
        if conductor not in labels:
            labels.append(conductor)
        if conductor != labels.index(conductor) + 1:
            raise ValueError(
                f"{coupon['Id']} conductor labels are not canonical"
            )
    if (
        abs(edges[0]["Offset"]) > tolerance
        or any(
            second["Offset"] - first["Offset"] <= tolerance
            for first, second in zip(edges, edges[1:])
        )
        or any(abs(edge["GapDirection"]) != 1 for edge in edges)
    ):
        raise ValueError(f"{coupon['Id']} has invalid parallel-edge geometry")
    # A stack is bounded per consecutive link (two consecutive edges interact iff their
    # separation is below 2R; decision 82(2)): its total span may exceed 2R (the 4-edge
    # 2 / 2 / 2 um flux line at R 1.9 spans 3.16 R).
    if any(
        second["Offset"] - first["Offset"] > 2.0 * matching_radius + tolerance
        for first, second in zip(edges, edges[1:])
    ):
        raise ValueError(
            f"{coupon['Id']} has consecutive edges farther apart than the 2R interaction "
            "distance"
        )
    return edges


def write_cluster_signature(root, edges, local_pieces=False):
    """The JSON signature (the key's labels, read by the response generator) and the mesher's CSV; with
    local_pieces the CSV's Conductor column is the LOCAL piece (R0-b: the mesher tags the metal per piece through
    its existing per-label tagging, no mesher change), the JSON unchanged."""
    data = {"Edges": edges}
    json_path = root / "cluster-signature.json"
    write_json(json_path, data)
    csv_path = root / "cluster-signature.csv"
    labels = [edge["Conductor"] for edge in edges]
    if local_pieces:
        labels = local_pieces_module.local_pieces_of_parallel_edges(
            [edge["Offset"] for edge in edges], [edge["GapDirection"] for edge in edges])
    lines = ["Offset,GapDirection,Conductor"]
    lines.extend(
        f"{edge['Offset']:.17g},{edge['GapDirection']},{label}"
        for edge, label in zip(edges, labels)
    )
    csv_path.write_text("\n".join(lines) + "\n")
    return json_path, csv_path


def straight_coupon_directory(root, coupon):
    if coupon["Topology"] == "IsolatedEdge":
        return root / "coupons" / "isolated-edge"
    topology = {
        "SameConductorGap": "same-conductor-gap",
        "DifferentConductorGap": "different-conductor-gap",
        "SameConductorStrip": "same-conductor-strip",
    }[coupon["Topology"]]
    width = float(coupon["Geometry"]["Separation"])
    return root / "coupons" / topology / f"w{slug(width)}um"


def straight_spec(coupons, args, parameters, local_pieces=None):
    """The straight-edge planner cache spec (its fingerprint names the cache root); `local_pieces` None = the
    family's planned form (straight_family_local_pieces)."""
    if local_pieces is None:
        local_pieces = straight_family_local_pieces(coupons)
    resolution = process_resolution(
        parameters, args.straight_lc_fine, args.min_process_feature_elements
    )
    spec = {
        "Version": 1,
        "Family": "straight",
        "Fabrication": parameters,
        "MatchingRadius": args.matching_radius,
        "Coupons": [
            {
                "Topology": coupon["Topology"],
                "Geometry": coupon["Geometry"],
                "BoundaryCondition": coupon["BoundaryCondition"],
            }
            for coupon in coupons
        ],
        "Orders": args.orders,
        "Mesh": {
            "FineSize": args.straight_lc_fine,
            "FarSize": args.straight_lc_far,
            "Order": args.mesh_order,
            "EdgeDistances": sorted(getattr(args, "straight_edge_distances", None) or [0.2]),
        },
        "ProcessResolution": resolution,
        "Response": {
            "BasisSize": args.basis_size,
            "Samples": args.samples,
            "CouponDepth": args.coupon_depth,
        },
        "ToolFingerprint": tool_fingerprint(
            (STRAIGHT_BUILDER, STRAIGHT_QUALIFIER, Path(__file__))
        ),
    }
    add_local_pieces_form(spec, local_pieces)
    return spec


def parallel_cluster_model_name(key, edges, local_pieces):
    """parallel-edge-cluster-<key[:12]>; the piece form carries the -pieces<P> suffix (local_pieces_model_name)."""
    name = f"parallel-edge-cluster-{key[:12]}"
    if local_pieces:
        pieces = local_pieces_module.local_pieces_of_parallel_edges(
            [edge["Offset"] for edge in edges], [edge["GapDirection"] for edge in edges])
        name = local_pieces_module.local_pieces_model_name(name, pieces)
    return name


def add_local_pieces_form(spec, local_pieces):
    """R0-b (DESIGN A.3 (b); review d1099 MAJOR-1): the piece-labelled form of a cross-section is not its one-reference
    form, so the planner cache key must carry the form; the key is added ONLY in the piece form so that every
    one-reference spec (and its cache root) stays byte-identical."""
    if local_pieces:
        spec["LocalPieces"] = True
    return spec


# THE PLANNER DEFAULT (decision 1136 (4); the 1134 result: a one-reference same-conductor CPW coupon is wrong by
# -57 / -25 / -38 % under a strip-vs-ground field, the three-piece coupon within 1 pt): every coupon whose
# cross-section / plan view holds a LABEL with >= 2 local metal pieces (a one-label multi-piece or a multi-piece-label
# class: ParallelEdgeCluster by local_pieces_of_parallel_edges, SameConductorGap (its two sides), SpatialEdgeCluster
# by cluster_signature_geometry.local_pieces_of_arrangement) is PLANNED in the piece form (Preparation.LocalPieces
# true, read by the builders here and by device_coupons); a class with one piece per label keeps the one-reference
# form (identical outputs). --one-reference opts every coupon out (the ES record form of v3-j, decisions 994 / 1136 (1)).
# A coupon of a class WITH a piece form whose pieces are NOT derivable (the refusals by name of the piece rules, or no
# version-3 Fabrication for a version-2 spatial record) is Unsupported, never the one-reference form silently (1145).
LOCAL_PIECES_PLANNER_RULE = ("decision 1136 (4): the piece form (--local-pieces) is the DEFAULT for every coupon whose "
                             "cross-section or plan view holds a conductor label with >= 2 local metal pieces "
                             "(ParallelEdgeCluster / SameConductorGap / SpatialEdgeCluster); one piece per label = the "
                             "one-reference form; a piece form that is not derivable = Unsupported (decision 1145); "
                             "--one-reference opts out (the ES record form)")


def local_pieces_per_label(coupon, matching_radius, parameters=None):
    """{label: piece count} of a coupon's cross-section / plan view, or None for a class without a piece form
    (IsolatedEdge, strips, DifferentConductorGap, corners, curved, Endpoint / Junction, a version-1 spatial record).
    ValueError when the class has a piece form that is not derivable from this coupon (fail closed, by name)."""
    topology = coupon["Topology"]
    if topology == "SameConductorGap":
        return {1: 2}
    if topology == "ParallelEdgeCluster":
        edges = normalize_parallel_edges(coupon, matching_radius)
        pieces = local_pieces_module.local_pieces_of_parallel_edges(
            [edge["Offset"] for edge in edges], [edge["GapDirection"] for edge in edges])
        per_label = {}
        for edge, piece in zip(edges, pieces):
            per_label.setdefault(int(edge["Conductor"]), set()).add(piece)
        return {label: len(found) for label, found in sorted(per_label.items())}
    geometry = coupon.get("Geometry") or {}
    if topology == "SpatialEdgeCluster" and "Signature" in geometry and not geometry.get("Edges"):
        if parameters is None:
            raise ValueError("the spatial pieces need the process parameters (version-3 Fabrication metadata)")
        if str(SPATIAL_ROOT) not in sys.path:
            sys.path.insert(0, str(SPATIAL_ROOT))
        import cluster_signature_geometry  # noqa: E402  (imports this module; resolved at call time)
        try:
            built, _ = cluster_signature_geometry.cluster_coupon(
                copy.deepcopy(coupon), matching_radius, parameters["metal_thickness"], parameters["overetch"],
                local_pieces=True)
        except cluster_signature_geometry.SignatureGeometryError as error:
            raise ValueError(str(error)) from error
        per_label = {}
        for row in built["Geometry"]["Edges"]:
            per_label.setdefault(int(row["Conductor"]), set()).add(int(row["LocalPiece"]))
        return {label: len(found) for label, found in sorted(per_label.items())}
    return None


def plan_local_pieces_form(coupon, matching_radius, parameters=None, one_reference=False):
    """Stamps coupon["Preparation"] with LocalPieces (bool), LocalPiecesReason and PiecesPerLabel (the rule above).
    A class WITH a piece form whose pieces are not derivable from this coupon (local_pieces_per_label's ValueError:
    the spatial arrangement refusals, lane C's parallel-edge refusals, a library without version-3 Fabrication for a
    version-2 spatial record) FAILS CLOSED: Preparation.Method = Unsupported with the refusal as the Reason (decision
    1145 Q1 (a): the one-reference form of such a coupon is the model class decision 1134 found wrong by tens of
    percent, so it is never planned silently; --one-reference builds it explicitly)."""
    preparation = coupon["Preparation"]
    if one_reference:
        preparation["LocalPieces"] = False
        preparation["LocalPiecesReason"] = "--one-reference: the one-reference form on every coupon (the ES record form)"
        return preparation
    try:
        per_label = local_pieces_per_label(coupon, matching_radius, parameters)
    except ValueError as error:
        preparation["LocalPieces"] = False
        preparation["LocalPiecesReason"] = f"the pieces are not derivable from this coupon: {error}"
        if preparation["Method"] != "Unsupported":
            # the builder keys (Driver / MeshGenerator / ResponseGenerator) go: the coupon is not built (the
            # Unsupported shape of preparation(): Method + Reason)
            coupon["Preparation"] = preparation = {
                "Method": "Unsupported",
                "Reason": (f"the piece form is the default (decision 1136 (4)) and its pieces are not derivable from "
                           f"this coupon: {error}; --one-reference builds the one-reference form"),
                **{key: value for key, value in preparation.items()
                   if key in ("BoundaryLawQualification", "LocalPieces", "LocalPiecesReason")}}
        return preparation
    if per_label is None:
        preparation["LocalPieces"] = False
        preparation["LocalPiecesReason"] = f"{coupon['Topology']}: a class without a piece form (one piece per label)"
        return preparation
    preparation["PiecesPerLabel"] = {str(label): count for label, count in per_label.items()}
    multi = sorted(label for label, count in per_label.items() if count >= 2)
    preparation["LocalPieces"] = bool(multi)
    preparation["LocalPiecesReason"] = (
        f"label(s) {multi} hold >= 2 local metal pieces: the piece form ({LOCAL_PIECES_PLANNER_RULE})" if multi
        else "one local metal piece per label: the one-reference form")
    return preparation


def coupon_local_pieces(coupon):
    """The planned form of a coupon (plan_local_pieces_form); a plan written before the rule carries none = the
    one-reference form."""
    return bool((coupon.get("Preparation") or {}).get("LocalPieces", False))


def straight_family_local_pieces(coupons):
    """The straight family is built by one builder call: its --local-pieces (the SameConductorGap's two sides) is on
    when any SameConductorGap coupon of the family is planned in the piece form."""
    return any(coupon["Topology"] == "SameConductorGap" and coupon_local_pieces(coupon) for coupon in coupons)


# The held-out escalation rule (decisions 1116 / 1126 (6) / 1128 Q3; the evidence: RG-7's one-reference 3 / 3 / 3 coupon
# FAILED its held-out check at the planner defaults (basis 96 / samples 1200: fabricated MS +11.09 % vs the 10 % bar, PBS 83667)
# and QUALIFIED at 1.5x (basis 144 / samples 1800: MS +6.40 %, PBS 83740); T2 ran on that library). On a held-out-ONLY failure
# (the matrix + convergence checks pass) ONE retry at 1.5x --basis-size / --samples; a second failure = Unqualified, fail closed.
# Only the two 2D families carry the lever: ParallelClusterCoupon (per coupon) and StraightEdgeBuilder (family-wide); corner /
# spatial coupons are not escalated. The escalated build has its own cache root and model Name (the fingerprint carries
# BasisSize / Samples: d54065837782 at the defaults -> 9c177156962b at 1.5x on RG-7). Default ON (--no-heldout-escalation).
HELDOUT_ESCALATION_FACTOR = 1.5
HELDOUT_ESCALATION_FAMILIES = ("StraightEdgeBuilder", "ParallelClusterCoupon")
HELDOUT_ESCALATION_NOT_APPLICABLE = ("CornerCoupon", "SpatialCoupon")
HELDOUT_ESCALATION_RULE = ("decision 1116: on a held-out-only qualification failure ONE retry at 1.5x --basis-size / --samples; a "
                           "second failure = Unqualified (fail closed); the two attempts' records point at each other")


class HeldOutQualificationFailure(RuntimeError):
    """A coupon family failed its qualification ONLY on the held-out checks (the matrix / convergence checks passed); the
    escalation rule applies. qualification_path = the attempt's planner record (already written)."""

    def __init__(self, message, qualification_path):
        super().__init__(message)
        self.qualification_path = Path(qualification_path)


def heldout_only_failure(report):
    """True when a qualify_surface_response.py report fails ONLY on its HeldoutChecks."""
    if report.get("Passed"):
        return False
    other = list(report.get("MatrixChecks") or []) + list(report.get("ConvergenceChecks") or [])
    heldout = list(report.get("HeldoutChecks") or [])
    return all(check.get("Passed") for check in other) and any(not check.get("Passed") for check in heldout)


def escalated_args(args, factor=HELDOUT_ESCALATION_FACTOR):
    """A copy of the planner arguments at factor x --basis-size (rounded up to an even number: the generators require it) and
    factor x --samples (at least the basis size); everything else identical."""
    escalated = copy.copy(args)
    basis = int(math.ceil(args.basis_size * factor))
    escalated.basis_size = basis + (basis % 2)
    escalated.samples = max(int(math.ceil(args.samples * factor)), escalated.basis_size)
    return escalated


def record_heldout_escalation(qualification_path, block):
    """Merge `block` into the planner record's HeldOutEscalation (the coupon's record of the rule's course)."""
    qualification = load_json(qualification_path)
    qualification["HeldOutEscalation"] = {**(qualification.get("HeldOutEscalation") or {}), **block}
    write_json(qualification_path, qualification)


def with_heldout_escalation(builder, args, label):
    """Run builder(args); on a HeldOutQualificationFailure retry ONCE with escalated_args (when the rule is on); a second
    held-out failure = Unqualified, fail closed, both records cross-referenced. -> (library, qualification path)."""
    if not getattr(args, "heldout_escalation", True):
        return builder(args)
    try:
        return builder(args)
    except HeldOutQualificationFailure as first:
        escalated = escalated_args(args)
        settings = {"BasisSize": args.basis_size, "Samples": args.samples}
        escalated_settings = {"BasisSize": escalated.basis_size, "Samples": escalated.samples}
        print(f"Held-out qualification failed for {label} at basis {args.basis_size} / samples {args.samples}: ONE escalation "
              f"at {HELDOUT_ESCALATION_FACTOR}x (basis {escalated.basis_size} / samples {escalated.samples}; {HELDOUT_ESCALATION_RULE})",
              flush=True)
        record_heldout_escalation(first.qualification_path, {"Rule": HELDOUT_ESCALATION_RULE, "Attempt": 1, "Factor": HELDOUT_ESCALATION_FACTOR,
                                                              "Settings": settings, "Outcome": "HeldOutFailure", "Escalation": escalated_settings})
        try:
            library, qualification_path = builder(escalated)
        except HeldOutQualificationFailure as second:
            record_heldout_escalation(second.qualification_path, {"Rule": HELDOUT_ESCALATION_RULE, "Attempt": 2, "Factor": HELDOUT_ESCALATION_FACTOR,
                                                                   "Settings": escalated_settings, "Outcome": "HeldOutFailure", "Unqualified": True,
                                                                   "FirstAttempt": {"Qualification": str(first.qualification_path), **settings}})
            record_heldout_escalation(first.qualification_path, {"SecondAttempt": str(second.qualification_path), "Unqualified": True})
            raise RuntimeError(
                f"{label}: Unqualified (fail closed, {HELDOUT_ESCALATION_RULE}): the held-out check failed at basis {args.basis_size} / "
                f"samples {args.samples} ({first.qualification_path}) and at the {HELDOUT_ESCALATION_FACTOR}x escalation basis "
                f"{escalated.basis_size} / samples {escalated.samples} ({second.qualification_path})") from second
        record_heldout_escalation(qualification_path, {"Rule": HELDOUT_ESCALATION_RULE, "Attempt": 2, "Factor": HELDOUT_ESCALATION_FACTOR,
                                                       "Settings": escalated_settings, "Outcome": "Qualified",
                                                       "FirstAttempt": {"Qualification": str(first.qualification_path), **settings}})
        record_heldout_escalation(first.qualification_path, {"EscalatedTo": str(qualification_path)})
        return library, qualification_path


def raise_recorded_heldout_failure(qualification, key, qualification_path, args, family):
    """A cache root whose record is a held-out failure at THESE settings (the same fingerprint) is not re-solved: the recorded
    failure is raised again so that the escalation (or the fail-closed verdict) follows from the record; --force rebuilds."""
    if (
        qualification
        and qualification.get("Fingerprint") == key
        and qualification.get("HeldOutFailure")
        and not args.force
    ):
        raise HeldOutQualificationFailure(
            f"{family} held-out qualification failed at these settings (recorded, not re-solved; --force rebuilds): {qualification_path}",
            qualification_path,
        )


def build_straight(coupons, args, parameters, cache):
    for coupon in coupons:
        require_verified_boundary_conditions(coupon)
    separations = {
        topology: sorted(
            {
                float(coupon["Geometry"]["Separation"])
                for coupon in coupons
                if coupon["Topology"] == topology
            }
        )
        for topology in PAIRED_TOPOLOGIES
    }
    isolated = any(coupon["Topology"] == "IsolatedEdge" for coupon in coupons)
    local_pieces = straight_family_local_pieces(coupons)
    spec = straight_spec(coupons, args, parameters, local_pieces)
    key = fingerprint(spec)
    root = cache / f"straight-{key}"
    qualification_path = root / "qualification.json"
    final_library = root / f"p{max(args.orders)}" / "library" / "process-library.json"
    qualification = (
        load_json(qualification_path) if qualification_path.is_file() else None
    )
    if (
        qualification
        and qualification.get("Fingerprint") == key
        and (
            qualification.get("Passed")
            or (
                getattr(args, "coverage_only", False)
                and qualification.get("GeometryComplete")
            )
        )
        and final_library.is_file()
        and not args.force
    ):
        print(f"Reusing completed straight-edge cache {root}")
        return final_library, qualification_path
    raise_recorded_heldout_failure(qualification, key, qualification_path, args, "Straight-edge")

    root.mkdir(parents=True, exist_ok=True)
    write_json(root / "coupon-spec.json", spec)
    order_roots = {}
    for order in args.orders:
        order_root = root / f"p{order}"
        order_roots[order] = order_root
        command = [
            sys.executable,
            STRAIGHT_BUILDER,
            "--output",
            order_root,
            "--name",
            f"{args.name}-straight-p{order}",
            "--palace",
            args.palace,
            "--julia",
            args.julia,
            "--ranks",
            args.ranks,
            "--matching-radius",
            args.matching_radius,
            "--metal-thickness",
            parameters["metal_thickness"],
            "--overetch",
            parameters["overetch"],
            "--sidewall-angle",
            parameters["sidewall_angle"],
            "--top-radius",
            parameters["top_radius"],
            "--bottom-radius",
            parameters["bottom_radius"],
            "--lc-fine",
            args.straight_lc_fine,
            "--lc-far",
            args.straight_lc_far,
            "--mesh-order",
            args.mesh_order,
            "--basis-size",
            args.basis_size,
            "--samples",
            args.samples,
            "--order",
            order,
            "--coupon-depth",
            args.coupon_depth,
            *material_options(parameters),
        ]
        if getattr(args, "straight_edge_distances", None):
            command.extend(["--edge-distances", *[str(d) for d in args.straight_edge_distances]])
        if args.julia_project:
            command.extend(["--julia-project", args.julia_project])
        if not isolated:
            command.append("--skip-isolated")
        if local_pieces:
            command.append("--local-pieces")
        for topology, option in PAIRED_TOPOLOGIES.items():
            command.append(f"--{option}")
            command.extend(separations[topology])
        run(command)

    if getattr(args, "coverage_only", False):
        stamp_library_boundary_conditions(final_library, coupons)
        stamp_geometry_coverage_only(final_library)
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "straight",
            "Library": str(final_library),
            "GeometryComplete": True,
            "Passed": False,
            "Reason": "Geometry-coverage generation does not perform convergence qualification",
        }
        write_json(qualification_path, qualification)
        return final_library, qualification_path

    previous = order_roots[args.orders[-2]]
    current = order_roots[args.orders[-1]]
    checks = []
    reports = []
    passed = True
    for coupon in coupons:
        old_coupon = straight_coupon_directory(previous, coupon)
        new_coupon = straight_coupon_directory(current, coupon)
        report = new_coupon / "qualification.json"
        code = run(
            [
                sys.executable,
                STRAIGHT_QUALIFIER,
                new_coupon,
                "--previous",
                old_coupon,
                "--output",
                report,
                "--max-heldout-error",
                args.max_heldout_error,
                "--max-fabricated-matrix-change",
                args.max_fabricated_matrix_change,
                "--max-fabricated-energy-change",
                args.max_fabricated_energy_change,
                "--max-domain-defect-change",
                args.max_domain_defect_change,
            ],
            check=False,
        )
        result = load_json(report)
        checks.append(
            {
                "Coupon": coupon["Id"],
                "Report": str(report),
                "Passed": result.get("Passed", False),
            }
        )
        passed = passed and code == 0 and result.get("Passed", False)
        reports.append(result)
    # held-out-ONLY when every failing coupon of the family fails on its held-out checks alone
    heldout_failure = not passed and all(r.get("Passed") or heldout_only_failure(r) for r in reports)
    qualification = {
        "Version": 1,
        "Fingerprint": key,
        "Family": "straight",
        "Library": str(final_library),
        "Checks": checks,
        "Passed": passed,
        "HeldOutFailure": heldout_failure,
    }
    write_json(qualification_path, qualification)
    if heldout_failure:
        raise HeldOutQualificationFailure(f"Straight-edge held-out qualification failed: {qualification_path}", qualification_path)
    if not passed:
        raise RuntimeError(f"Straight-edge qualification failed: {qualification_path}")
    stamp_library_boundary_conditions(final_library, coupons)
    return final_library, qualification_path


def cluster_complete(root):
    matrix_names = ("domain-response-matrix.csv", "surface-response-matrix.csv")
    heldout_names = ("domain-E.csv", "surface-Q.csv", "surface-Q-edge.csv")
    return all(
        (root / "postpro" / f"cluster_{kind}" / name).is_file()
        for kind in ("thin", "fabricated")
        for name in matrix_names
    ) and all(
        (root / "postpro" / f"heldout_cluster_{kind}" / name).is_file()
        for kind in ("thin", "fabricated")
        for name in heldout_names
    )


def run_probe_convergence(calibration, output, args):
    command = [
        sys.executable,
        CORNER_CONVERGENCE,
        calibration,
        "--output",
        output,
        "--palace",
        args.palace,
        "--orders",
        *args.orders,
        "--ranks",
        args.ranks,
        "--max-fabricated-matrix-change",
        args.max_fabricated_matrix_change,
        "--max-fabricated-energy-change",
        args.max_fabricated_energy_change,
        "--max-domain-defect-change",
        args.max_domain_defect_change,
    ]
    if args.force:
        command.append("--force")
    code = run(command, check=False)
    report = output / "probe-convergence.json"
    result = (
        load_json(report)
        if report.is_file()
        else {
            "Version": 1,
            "Passed": False,
            "Failure": "Probe runner did not write a convergence report",
        }
    )
    return code, report, result


def run_mesh_convergence(calibrations, output, args):
    command = [
        sys.executable,
        CORNER_CONVERGENCE,
        "--output",
        output,
        "--palace",
        args.palace,
        "--fixed-order",
        max(args.orders),
        "--ranks",
        args.ranks,
        "--max-fabricated-matrix-change",
        args.max_fabricated_matrix_change,
        "--max-fabricated-energy-change",
        args.max_fabricated_energy_change,
        "--max-domain-defect-change",
        args.max_domain_defect_change,
    ]
    for name, calibration in calibrations:
        command.extend(("--case", f"{name}={calibration}"))
    if args.force:
        command.append("--force")
    code = run(command, check=False)
    report = output / "probe-convergence.json"
    result = (
        load_json(report)
        if report.is_file()
        else {
            "Version": 1,
            "Study": "MeshResolution",
            "Passed": False,
            "Failure": "Probe runner did not write a mesh-convergence report",
        }
    )
    return code, report, result


def generate_parallel_cluster_meshes(
    root, csv_signature, args, parameters, spec, factor
):
    mesh_root = (
        root
        if math.isclose(factor, 1.0)
        else root / "mesh-calibrations" / f"h-{slug(factor)}"
    )
    mesh_root.mkdir(parents=True, exist_ok=True)
    meshes = {}
    for kind in ("thin", "fabricated"):
        mesh = mesh_root / f"cluster_{kind}.msh"
        meshes[kind] = mesh
        if mesh.is_file() and not args.force:
            continue
        run(
            [
                *julia_command(args),
                CLUSTER_MESH,
                csv_signature,
                kind,
                mesh,
                "--radius",
                args.matching_radius,
                "--metal-thickness",
                parameters["metal_thickness"],
                "--overetch",
                parameters["overetch"],
                "--sidewall-angle",
                parameters["sidewall_angle"],
                "--top-radius",
                parameters["top_radius"],
                "--bottom-radius",
                parameters["bottom_radius"],
                "--lc-fine",
                factor * args.cluster_lc_fine,
                "--lc-far",
                args.cluster_lc_far,
                "--mesh-order",
                spec["Mesh"]["Order"],
            ]
        )
    return mesh_root, meshes


def parallel_cluster_spec(coupon, edges, args, parameters, local_pieces=None):
    """The parallel-cluster planner cache spec; `local_pieces` None = the coupon's planned form (coupon_local_pieces)."""
    if local_pieces is None:
        local_pieces = coupon_local_pieces(coupon)
    """The parallel-edge-cluster planner cache spec (its fingerprint names the cache root and the model Name)."""
    resolution = process_resolution(
        parameters, args.cluster_lc_fine, args.min_process_feature_elements
    )
    spec = {
        "Version": 1,
        "Family": "parallel-edge-cluster",
        "Edges": edges,
        "BoundaryCondition": coupon["BoundaryCondition"],
        "Fabrication": parameters,
        "MatchingRadius": args.matching_radius,
        "Orders": args.orders,
        "Mesh": {
            "FineSize": args.cluster_lc_fine,
            "FarSize": args.cluster_lc_far,
            "Order": max(2, args.mesh_order),
            "HRefinementFactors": sorted(
                set(getattr(args, "cluster_h_factors", (2.0, 1.0))),
                reverse=True,
            ),
        },
        "ProcessResolution": resolution,
        "Response": {
            "BasisSize": args.basis_size,
            "Samples": args.samples,
            "CouponDepth": args.coupon_depth,
            "EdgeOffsetTolerance": args.edge_offset_tolerance,
        },
        "ToolFingerprint": tool_fingerprint(
            (
                CLUSTER_MESH,
                CLUSTER_GENERATOR,
                STRAIGHT_QUALIFIER,
                Path(__file__),
            )
        ),
    }
    add_local_pieces_form(spec, local_pieces)
    return spec


def build_parallel_cluster(coupon, args, parameters, cache):
    require_verified_boundary_conditions(coupon)
    edges = normalize_parallel_edges(coupon, args.matching_radius)
    local_pieces = coupon_local_pieces(coupon)
    spec = parallel_cluster_spec(coupon, edges, args, parameters, local_pieces)
    key = fingerprint(spec)
    root = cache / f"parallel-cluster-{key}"
    qualification_path = root / "qualification.json"
    final_library = root / f"p{max(args.orders)}" / "process-library.json"
    qualification = (
        load_json(qualification_path) if qualification_path.is_file() else None
    )
    if (
        qualification
        and qualification.get("Fingerprint") == key
        and (
            qualification.get("Passed")
            or (
                getattr(args, "coverage_only", False)
                and qualification.get("GeometryComplete")
            )
            or (
                getattr(args, "probe_study_only", False)
                and qualification.get("ProbeStudyComplete")
            )
        )
        and (final_library.is_file() or getattr(args, "probe_study_only", False))
        and not args.force
    ):
        print(f"Reusing completed parallel-cluster cache {root}")
        return (
            None if getattr(args, "probe_study_only", False) else final_library,
            qualification_path,
        )
    raise_recorded_heldout_failure(qualification, key, qualification_path, args, "Parallel-cluster")

    root.mkdir(parents=True, exist_ok=True)
    write_json(root / "coupon-spec.json", spec)
    json_signature, csv_signature = write_cluster_signature(root, edges, local_pieces)
    _, meshes = generate_parallel_cluster_meshes(
        root, csv_signature, args, parameters, spec, 1.0
    )

    final_order = max(args.orders)
    final_root = root / f"p{final_order}"
    final_root.mkdir(parents=True, exist_ok=True)
    run(
        [
            sys.executable,
            CLUSTER_GENERATOR,
            json_signature,
            "--output",
            final_root,
            "--thin-mesh",
            meshes["thin"],
            "--fabricated-mesh",
            meshes["fabricated"],
            "--radius",
            args.matching_radius,
            "--metal-thickness",
            parameters["metal_thickness"],
            "--overetch-depth",
            parameters["overetch"],
            "--sidewall-angle",
            parameters["sidewall_angle"],
            "--top-rounding",
            parameters["top_radius"],
            "--trench-rounding",
            parameters["bottom_radius"],
            "--basis-size",
            args.basis_size,
            "--samples",
            args.samples,
            "--order",
            final_order,
            "--coupon-depth",
            args.coupon_depth,
            "--edge-offset-tolerance",
            args.edge_offset_tolerance,
            *(["--local-pieces"] if local_pieces else []),
            "--library-name",
            f"{args.name}-parallel-cluster-p{final_order}",
            "--model-name",
            parallel_cluster_model_name(key, edges, local_pieces),
            *material_options(parameters),
        ]
    )
    if getattr(args, "coverage_only", False):
        for name in ("cluster_thin", "cluster_fabricated"):
            run(
                palace_command(args, final_root / f"{name}.json"),
                cwd=final_root,
            )
        stamp_library_boundary_conditions(final_library, [coupon])
        stamp_geometry_coverage_only(final_library)
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "parallel-edge-cluster",
            "Library": str(final_library),
            "GeometryComplete": True,
            "Passed": False,
            "Reason": "Geometry-coverage generation does not perform convergence qualification",
        }
        write_json(qualification_path, qualification)
        return final_library, qualification_path

    convergence_root = root / "probe-convergence"
    convergence_code, convergence_report, convergence = run_probe_convergence(
        final_root, convergence_root, args
    )
    if convergence_code != 0 or not convergence.get("Passed", False):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "parallel-edge-cluster",
            "Library": str(final_library),
            "ResponseReport": None,
            "ConvergenceReport": str(convergence_report),
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        raise RuntimeError(
            f"Parallel-cluster probe convergence failed: {qualification_path}"
        )

    mesh_calibrations = []
    for factor in spec["Mesh"]["HRefinementFactors"]:
        if math.isclose(factor, 1.0):
            calibration = convergence_root / f"p{final_order}"
        else:
            mesh_root, factor_meshes = generate_parallel_cluster_meshes(
                root, csv_signature, args, parameters, spec, factor
            )
            calibration = prepare_probe_mesh_calibration(
                final_root, mesh_root / "calibration", factor_meshes
            )
        mesh_calibrations.append((f"h-{slug(factor)}", calibration))
    mesh_convergence_root = root / "mesh-convergence"
    (
        mesh_convergence_code,
        mesh_convergence_report,
        mesh_convergence,
    ) = run_mesh_convergence(mesh_calibrations, mesh_convergence_root, args)
    if (
        mesh_convergence_code != 0
        or not mesh_convergence.get("Passed", False)
    ):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "parallel-edge-cluster",
            "Library": str(final_library),
            "ResponseReport": None,
            "ConvergenceReport": str(convergence_report),
            "MeshConvergenceReport": str(mesh_convergence_report),
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        raise RuntimeError(
            f"Parallel-cluster mesh convergence failed: {qualification_path}"
        )
    if getattr(args, "probe_study_only", False):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "parallel-edge-cluster",
            "Library": None,
            "ResponseReport": None,
            "ConvergenceReport": str(convergence_report),
            "MeshConvergenceReport": str(mesh_convergence_report),
            "ProbeStudyComplete": True,
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        return None, qualification_path

    if cluster_complete(final_root) and not args.force:
        print(f"Reusing completed p{final_order} parallel-cluster responses")
    else:
        for name in (
            "cluster_thin",
            "cluster_fabricated",
            "heldout_cluster_thin",
            "heldout_cluster_fabricated",
        ):
            run(
                palace_command(args, final_root / f"{name}.json"),
                cwd=final_root,
            )
        if not cluster_complete(final_root):
            raise RuntimeError(
                f"Palace did not complete parallel-cluster responses in {final_root}"
            )

    response_report = final_root / "qualification.json"
    response_code = run(
        [
            sys.executable,
            STRAIGHT_QUALIFIER,
            final_root,
            "--output",
            response_report,
            "--max-heldout-error",
            args.max_heldout_error,
        ],
        check=False,
    )
    response = load_json(response_report)
    passed = response_code == 0 and response.get("Passed", False)
    heldout_failure = not passed and heldout_only_failure(response)
    qualification = {
        "Version": 1,
        "Fingerprint": key,
        "Family": "parallel-edge-cluster",
        "Library": str(final_library),
        "ResponseReport": str(response_report),
        "ConvergenceReport": str(convergence_report),
        "MeshConvergenceReport": str(mesh_convergence_report),
        "Passed": passed,
        "HeldOutFailure": heldout_failure,
    }
    write_json(qualification_path, qualification)
    if heldout_failure:
        raise HeldOutQualificationFailure(
            f"Parallel-cluster held-out qualification failed: {qualification_path}", qualification_path
        )
    if not passed:
        raise RuntimeError(
            f"Parallel-cluster qualification failed: {qualification_path}"
        )
    stamp_library_boundary_conditions(final_library, [coupon])
    return final_library, qualification_path


def check_corner_trace_resolvability(root, mesh_root, meshes, args):
    """The trace resolvability gate on a corner coupon's meshes (CORNER_RESOLVABILITY at the
    final solve order; report mesh_root/trace-resolvability.json); fail-closed."""
    report = mesh_root / "trace-resolvability.json"
    code = run(
        [
            sys.executable,
            CORNER_RESOLVABILITY,
            root,
            *(option for kind in ("thin", "fabricated") for option in ("--mesh", meshes[kind])),
            "--order",
            max(args.orders),
            "--radius",
            args.matching_radius,
            "--report",
            report,
        ],
        check=False,
    )
    # The gate exits 0 only when every mesh passed (its report records the counts per hat).
    return code == 0, report


def generate_corner_meshes(
    root, topology, angle, corner_radius, args, parameters, spec, factor
):
    """The corner coupon's meshes at lc_fine x factor, sized at the trace basis knots (the
    generator's trace mesh in `root`, written before the meshes: `--trace-mesh`; the mesher
    fails closed without it) and gated for trace resolvability
    (check_corner_trace_resolvability; a failing mesh root raises)."""
    mesh_root = (
        root
        if math.isclose(factor, 1.0)
        else root / "mesh-calibrations" / f"h-{slug(factor)}"
    )
    mesh_root.mkdir(parents=True, exist_ok=True)
    meshes = {}
    for kind in ("thin", "fabricated"):
        mesh = mesh_root / f"corner_{kind}.msh"
        meshes[kind] = mesh
        if mesh.is_file() and not args.force:
            continue
        run(
            [
                *julia_command(args),
                CORNER_MESH,
                f"{topology}-{kind}",
                mesh,
                "--trace-mesh",
                root,
                "--radius",
                args.matching_radius,
                "--angle",
                angle,
                "--corner-radius",
                corner_radius,
                "--metal-thickness",
                parameters["metal_thickness"],
                "--overetch",
                parameters["overetch"],
                "--sidewall-angle",
                parameters["sidewall_angle"],
                "--top-radius",
                parameters["top_radius"],
                "--bottom-radius",
                parameters["bottom_radius"],
                "--lc-fine",
                factor * args.corner_lc_fine,
                "--lc-far",
                args.corner_lc_far,
                "--mesh-order",
                spec["Mesh"]["Order"],
            ]
        )
    passed, report = check_corner_trace_resolvability(root, mesh_root, meshes, args)
    if not passed:
        raise RuntimeError(
            f"Corner coupon meshes in {mesh_root} do not resolve the trace basis: {report}"
        )
    return mesh_root, meshes


def corner_resolvability_constant():
    """MINIMUM_ACTIVE_NODES_PER_ORDER_SQUARED of the trace resolvability gate (recorded in the
    coupon spec)."""
    if str(CORNER_ROOT) not in sys.path:
        sys.path.insert(0, str(CORNER_ROOT))
    import trace_resolvability

    return trace_resolvability.MINIMUM_ACTIVE_NODES_PER_ORDER_SQUARED


def corner_trace_basis(coupon_id, args, topology="convex"):
    """The trace basis rule of a corner coupon: --corner-trace-basis when given, else the
    refined default. The rule does not depend on the corner radius (the box basis never
    sees the fillet; CORNER_TRACE_BASIS_DEFAULT): a rounded corner resolves exactly as the
    sharp corner of its angle, and an explicit legacy request is honoured on either. A
    convexity-conditional SELECTOR (CORNER_TRACE_BASIS_SELECTORS) resolves to the concrete
    rule of the node's topology ("convex" / "concave")."""
    requested = getattr(args, "corner_trace_basis", None)
    if requested is None:
        return CORNER_TRACE_BASIS_DEFAULT
    if requested in CORNER_TRACE_BASIS_SELECTORS:
        if topology not in CORNER_TRACE_BASIS_SELECTORS[requested]:
            raise ValueError(
                f"{coupon_id}: the selector {requested!r} has no rule for topology {topology!r}"
            )
        return CORNER_TRACE_BASIS_SELECTORS[requested][topology]
    if requested not in CORNER_TRACE_BASIS_RING_SIZES:
        raise ValueError(f"{coupon_id}: unknown corner trace basis rule {requested!r}")
    return requested


def build_corner(coupon, args, parameters, cache):
    require_verified_boundary_conditions(coupon)
    resolution = process_resolution(
        parameters, args.corner_lc_fine, args.min_process_feature_elements
    )
    geometry = coupon["Geometry"]
    angle = float(geometry["AngleDegrees"])
    if not math.isfinite(angle) or not 0.0 < angle <= 180.0:
        raise ValueError(f"{coupon['Id']} requests invalid corner angle {angle:g}")
    topology = "convex" if coupon["Topology"] == "ConvexCorner" else "concave"
    corner_radius = float(geometry.get("CornerRadius", 0.0))
    if angle == 180.0 and corner_radius > 0.0:
        raise ValueError(f"{coupon['Id']}: the straight anchor (180 deg) has no corner radius")
    # Corner family segment connectivity (corner-qualification block 2026-09-29): the node's
    # band triangulation follows the rule's layout at this angle (generate_corner_response.py
    # --connectivity-angle, record TraceBasis ConnectivityAngleDegrees); part of the coupon
    # spec and Id (two coupons at one angle, one per side of a knot-corner passage).
    connectivity_angle = geometry.get("ConnectivityAngleDegrees")
    if connectivity_angle is not None:
        connectivity_angle = float(connectivity_angle)
        if not 0.0 < connectivity_angle < 180.0 or corner_radius > 0.0:
            raise ValueError(
                f"{coupon['Id']} ConnectivityAngleDegrees must lie in (0, 180) on a sharp corner"
            )
    # The trace basis rule (corner-basis refinement 2026-09-30, USER decision 161 (1)): the
    # refined AllRingsFollowMetal layout by default on every corner, sharp or rounded (no
    # events: no connectivity angle; the rule is independent of CornerRadius, decision 511),
    # the recorded MetalRingsOnly layout as "legacy" on request; part of the coupon spec and
    # Id. The ring size is the rule's (2 + MetalInteriorKnots + FreeKnots).
    trace_basis = corner_trace_basis(coupon["Id"], args, topology)
    corner_ring_size = CORNER_TRACE_BASIS_RING_SIZES[trace_basis]
    if trace_basis != "legacy" and connectivity_angle is not None:
        raise ValueError(
            f"{coupon['Id']} ConnectivityAngleDegrees is a legacy (MetalRingsOnly) trace "
            "basis option: the refined rule has no events"
        )
    if (
        not math.isfinite(corner_radius)
        or not 0.0 <= corner_radius < args.matching_radius
    ):
        raise ValueError(
            f"{coupon['Id']} corner radius must be finite and lie in "
            f"[0, {args.matching_radius:g})"
        )
    tangent_distance = 0.0 if angle == 180.0 else corner_radius / math.tan(0.5 * math.radians(angle))
    if corner_radius > 0.0 and tangent_distance >= args.matching_radius:
        raise ValueError(
            f"{coupon['Id']} rounded-corner tangency distance "
            f"{tangent_distance:g} must be smaller than matching radius "
            f"{args.matching_radius:g}"
        )
    geometry_order = max(2, args.mesh_order)
    spec = {
        "Version": 3,
        "Family": "corner",
        "Topology": topology,
        "AngleDegrees": angle,
        "CornerRadius": corner_radius,
        **({"ConnectivityAngleDegrees": connectivity_angle} if connectivity_angle is not None else {}),
        "BoundaryCondition": coupon["BoundaryCondition"],
        "Fabrication": parameters,
        "MatchingRadius": args.matching_radius,
        "Orders": args.orders,
        "Mesh": {
            "FineSize": args.corner_lc_fine,
            "FarSize": args.corner_lc_far,
            "Order": geometry_order,
            "HRefinementFactors": sorted(
                set(getattr(args, "corner_h_factors", (2.0, 1.0))),
                reverse=True,
            ),
        },
        "ProcessResolution": resolution,
        "Response": {
            "RingSize": corner_ring_size,
            "TraceBasis": trace_basis,
            # Decision 328: the mesh is sized at the trace knots and gated (the recipe change
            # of family-4 round 2; the fingerprint below covers the mesher and the gate).
            # Decision 787 (R-A) / 791: the knot size is min(along-ring gap, across-ring
            # spacing), named by the label; the fingerprint stays the byte discriminator.
            "TraceResolvability": {
                "MeshSizing": "KnotGapAndRingSpacing",
                "MinimumActiveNodesPerOrderSquared": corner_resolvability_constant(),
            },
        },
        "ToolFingerprint": tool_fingerprint(
            (
                CORNER_MESH,
                CORNER_GENERATOR,
                CORNER_FINALIZER,
                CORNER_CONVERGENCE,
                CORNER_RESOLVABILITY,
                Path(__file__),
            )
        ),
    }
    key = fingerprint(spec)
    root = cache / f"corner-{key}"
    library_path = root / "process-library.json"
    qualification_path = root / "qualification.json"
    qualification = (
        load_json(qualification_path) if qualification_path.is_file() else None
    )
    if (
        qualification
        and qualification.get("Fingerprint") == key
        and (
            qualification.get("Passed")
            or (
                getattr(args, "coverage_only", False)
                and qualification.get("GeometryComplete")
            )
            or (
                getattr(args, "probe_study_only", False)
                and qualification.get("ProbeStudyComplete")
            )
        )
        and (library_path.is_file() or getattr(args, "probe_study_only", False))
        and not args.force
    ):
        print(f"Reusing completed corner cache {root}")
        return (
            None if getattr(args, "probe_study_only", False) else library_path,
            qualification_path,
        )

    root.mkdir(parents=True, exist_ok=True)
    write_json(root / "coupon-spec.json", spec)
    # The generator first (it reads no mesh: the configs name the mesh files), then the
    # meshes sized at its trace knots and gated for resolvability (decision 328).
    meshes = {kind: root / f"corner_{kind}.msh" for kind in ("thin", "fabricated")}
    run(
        [
            sys.executable,
            CORNER_GENERATOR,
            "--output",
            root,
            "--thin-mesh",
            meshes["thin"],
            "--fabricated-mesh",
            meshes["fabricated"],
            "--radius",
            args.matching_radius,
            "--angle",
            angle,
            *(
                ["--connectivity-angle", connectivity_angle]
                if connectivity_angle is not None
                else []
            ),
            "--corner-radius",
            corner_radius,
            "--trace-basis",
            trace_basis,
            "--ring-size",
            corner_ring_size,
            "--order",
            max(args.orders),
            "--metal-thickness",
            parameters["metal_thickness"],
            "--overetch-depth",
            parameters["overetch"],
            "--sidewall-angle",
            parameters["sidewall_angle"],
            "--top-rounding",
            parameters["top_radius"],
            "--trench-rounding",
            parameters["bottom_radius"],
            "--topology",
            topology,
            *material_options(parameters),
        ]
    )
    try:
        _, generated = generate_corner_meshes(
            root,
            topology,
            angle,
            corner_radius,
            args,
            parameters,
            spec,
            1.0,
        )
    except RuntimeError as failure:
        write_json(
            qualification_path,
            {
                "Version": 1,
                "Fingerprint": key,
                "Family": "corner",
                "Library": str(library_path),
                "TraceResolvabilityReport": str(root / "trace-resolvability.json"),
                "Passed": False,
                "Reason": str(failure),
            },
        )
        raise
    if generated != meshes:
        raise RuntimeError(f"{coupon['Id']}: the corner meshes are not the generator's")
    if getattr(args, "coverage_only", False):
        for name in ("thin", "fabricated"):
            run(palace_command(args, root / f"{name}.json"), cwd=root)
        run([sys.executable, CORNER_FINALIZER, root])
        stamp_library_boundary_conditions(library_path, [coupon])
        stamp_geometry_coverage_only(library_path)
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "corner",
            "Library": str(library_path),
            "GeometryComplete": True,
            "Passed": False,
            "Reason": "Geometry-coverage generation does not perform convergence qualification",
        }
        write_json(qualification_path, qualification)
        return library_path, qualification_path

    convergence_root = root / "convergence"
    convergence_code, convergence_report, convergence = run_probe_convergence(
        root, convergence_root, args
    )
    if convergence_code != 0 or not convergence.get("Passed", False):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "corner",
            "Library": str(library_path),
            "HeldoutReport": None,
            "ConvergenceReport": str(convergence_report),
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        raise RuntimeError(f"Corner probe convergence failed: {qualification_path}")

    mesh_calibrations = []
    for factor in spec["Mesh"]["HRefinementFactors"]:
        if math.isclose(factor, 1.0):
            calibration = convergence_root / f"p{max(args.orders)}"
        else:
            try:
                mesh_root, factor_meshes = generate_corner_meshes(
                    root,
                    topology,
                    angle,
                    corner_radius,
                    args,
                    parameters,
                    spec,
                    factor,
                )
            except RuntimeError as failure:
                write_json(
                    qualification_path,
                    {
                        "Version": 1,
                        "Fingerprint": key,
                        "Family": "corner",
                        "Library": str(library_path),
                        "ConvergenceReport": str(convergence_report),
                        "TraceResolvabilityReport": str(
                            root
                            / "mesh-calibrations"
                            / f"h-{slug(factor)}"
                            / "trace-resolvability.json"
                        ),
                        "Passed": False,
                        "Reason": str(failure),
                    },
                )
                raise
            calibration = prepare_probe_mesh_calibration(
                root, mesh_root / "calibration", factor_meshes
            )
        mesh_calibrations.append((f"h-{slug(factor)}", calibration))
    mesh_convergence_root = root / "mesh-convergence"
    (
        mesh_convergence_code,
        mesh_convergence_report,
        mesh_convergence,
    ) = run_mesh_convergence(mesh_calibrations, mesh_convergence_root, args)
    if (
        mesh_convergence_code != 0
        or not mesh_convergence.get("Passed", False)
    ):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "corner",
            "Library": str(library_path),
            "HeldoutReport": None,
            "ConvergenceReport": str(convergence_report),
            "MeshConvergenceReport": str(mesh_convergence_report),
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        raise RuntimeError(
            f"Corner mesh convergence failed: {qualification_path}"
        )
    if getattr(args, "probe_study_only", False):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "corner",
            "Library": None,
            "HeldoutReport": None,
            "ConvergenceReport": str(convergence_report),
            "MeshConvergenceReport": str(mesh_convergence_report),
            "ProbeStudyComplete": True,
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        return None, qualification_path

    for name in ("thin", "fabricated"):
        run(palace_command(args, root / f"{name}.json"), cwd=root)
    run([sys.executable, CORNER_FINALIZER, root])
    for name in ("heldout-thin", "heldout-fabricated"):
        run(palace_command(args, root / f"{name}.json"), cwd=root)
    heldout_report = root / "heldout-qualification.json"
    heldout_code = run(
        [
            sys.executable,
            CORNER_FINALIZER,
            root,
            "--report",
            heldout_report,
            "--require-heldout",
            "--max-heldout-error",
            args.max_heldout_error,
        ],
        check=False,
    )
    heldout = load_json(heldout_report)
    passed = heldout_code == 0 and heldout.get("Passed", False)
    qualification = {
        "Version": 1,
        "Fingerprint": key,
        "Family": "corner",
        "Library": str(library_path),
        "HeldoutReport": str(heldout_report),
        "ConvergenceReport": str(convergence_report),
        "MeshConvergenceReport": str(mesh_convergence_report),
        "Passed": passed,
    }
    write_json(qualification_path, qualification)
    if not passed:
        raise RuntimeError(f"Corner qualification failed: {qualification_path}")
    stamp_library_boundary_conditions(library_path, [coupon])
    return library_path, qualification_path


def spatial_complete(root):
    matrix_names = ("domain-response-matrix.csv", "surface-response-matrix.csv")
    heldout_names = ("domain-E.csv", "surface-Q.csv", "surface-Q-edge.csv")
    return all(
        (root / "postpro" / f"spatial_{kind}" / name).is_file()
        for kind in ("thin", "fabricated")
        for name in matrix_names
    ) and all(
        (root / "postpro" / f"heldout_spatial_{kind}" / name).is_file()
        for kind in ("thin", "fabricated")
        for name in heldout_names
    )


def spatial_mesh_path(args):
    mode = getattr(args, "spatial_mesher", "field")
    if mode == "field":
        return SPATIAL_FIELD_MESH
    if mode == "swept":
        return SPATIAL_SWEPT_MESH
    raise ValueError(f"Unknown spatial mesher {mode!r}")


def spatial_spec(coupon, args, parameters):
    mesher = spatial_mesh_path(args)
    resolution = process_resolution(
        parameters, args.spatial_lc_fine, args.min_process_feature_elements
    )
    geometry = copy.deepcopy(coupon["Geometry"])
    facets = geometry.pop("PlanViewFacets", [])
    if facets:
        process_axis = 1 if coupon["Topology"] == "SpatialEdgeCluster" else 2
        canonical_boundary = canonical_plan_view_boundary(
            facets, args.matching_radius, process_axis
        )
        if "PlanViewBoundary" in geometry:
            exported_boundary = geometry["PlanViewBoundary"]
            if (
                unclassified_plan_view_boundary(exported_boundary)
                != canonical_boundary
            ):
                raise ValueError(
                    "Palace PlanViewBoundary does not match its exported facets"
                )
        else:
            geometry["PlanViewBoundary"] = canonical_boundary
        geometry["MaskRegularization"] = {
            "Version": 1,
            "PhysicalBoundary": "TaperAndRound",
            "ContinuationBoundary": "Vertical",
        }
    return {
        "Version": 2,
        "Family": "spatial",
        "Coupon": {
            "Topology": coupon["Topology"],
            "Geometry": geometry,
            "Interfaces": coupon["Interfaces"],
            "BoundaryCondition": coupon["BoundaryCondition"],
        },
        "Fabrication": parameters,
        "MatchingRadius": args.matching_radius,
        "Orders": args.orders,
        "Mesh": {
            "Mesher": getattr(args, "spatial_mesher", "field"),
            "FineSize": args.spatial_lc_fine,
            "TangentialSize": getattr(args, "spatial_lc_tangent", 0.0),
            "FarSize": args.spatial_lc_far,
            "ProcessCoreWidth": getattr(args, "spatial_process_core_width", 0.0),
            "ProcessFineWidth": getattr(args, "spatial_process_fine_width", 0.0),
            "ProcessGradingPower": getattr(args, "spatial_process_grading_power", 1.7),
            "NormalGrowthRatio": getattr(args, "spatial_normal_growth_ratio", 1.4),
            "MaxNodes": getattr(args, "spatial_max_nodes", 500_000),
            "MaxElements": getattr(args, "spatial_max_elements", 2_000_000),
            "Order": max(2, args.mesh_order),
            "HRefinementFactors": sorted(
                set(getattr(args, "spatial_h_factors", (2.0, 1.0))),
                reverse=True,
            ),
        },
        "ProcessResolution": resolution,
        "Response": {"RingSize": args.spatial_ring_size},
        "RegularizationPolicy": geometry.get("MaskRegularization"),
        "ToolFingerprint": tool_fingerprint(
            (
                mesher,
                SPATIAL_GENERATOR,
                STRAIGHT_QUALIFIER,
                CORNER_CONVERGENCE,
                Path(__file__),
            )
        ),
    }


def generate_spatial_meshes(
    root,
    signature,
    mask,
    boundary,
    args,
    parameters,
    spec,
    factor,
):
    mesher = spatial_mesh_path(args)
    mesh_root = (
        root
        if math.isclose(factor, 1.0)
        else root / "mesh-calibrations" / f"h-{slug(factor)}"
    )
    mesh_root.mkdir(parents=True, exist_ok=True)
    meshes = {}
    for kind in ("thin", "fabricated"):
        mesh = mesh_root / f"spatial_{kind}.msh"
        meshes[kind] = mesh
        if mesh.is_file() and not args.force:
            continue
        run(
            [
                *julia_command(args),
                mesher,
                signature,
                kind,
                mesh,
                *(["--mask", mask] if mask.is_file() else []),
                *(["--boundary", boundary] if boundary.is_file() else []),
                "--radius",
                args.matching_radius,
                "--metal-thickness",
                parameters["metal_thickness"],
                "--overetch",
                parameters["overetch"],
                "--sidewall-angle",
                parameters["sidewall_angle"],
                "--top-radius",
                parameters["top_radius"],
                "--bottom-radius",
                parameters["bottom_radius"],
                "--lc-fine",
                factor * args.spatial_lc_fine,
                "--lc-tangent",
                getattr(args, "spatial_lc_tangent", 0.0),
                "--lc-far",
                args.spatial_lc_far,
                "--process-core-width",
                getattr(args, "spatial_process_core_width", 0.0),
                "--process-fine-width",
                getattr(args, "spatial_process_fine_width", 0.0),
                "--process-grading-power",
                getattr(args, "spatial_process_grading_power", 1.7),
                "--normal-growth-ratio",
                getattr(args, "spatial_normal_growth_ratio", 1.4),
                "--max-nodes",
                getattr(args, "spatial_max_nodes", 500_000),
                "--max-elements",
                getattr(args, "spatial_max_elements", 2_000_000),
                "--mesh-order",
                spec["Mesh"]["Order"],
            ]
        )
    return mesh_root, meshes


def prepare_probe_mesh_calibration(source, destination, meshes):
    destination.mkdir(parents=True, exist_ok=True)
    for name in ("probe-manifest.json", "process-library.json"):
        shutil.copy2(source / name, destination / name)
    for kind in ("thin", "fabricated"):
        config = load_json(source / f"probe-{kind}.json")
        config["Model"]["Mesh"] = str(meshes[kind])
        config["Problem"]["Output"] = str(
            destination / "postpro" / f"probe-{kind}"
        )
        write_json(destination / f"probe-{kind}.json", config)
    return destination


def build_spatial(coupon, args, parameters, cache):
    require_verified_boundary_conditions(coupon)
    spec = spatial_spec(coupon, args, parameters)
    key = fingerprint(spec)
    root = cache / f"spatial-{key}"
    library_path = root / "process-library.json"
    qualification_path = root / "qualification.json"
    qualification = (
        load_json(qualification_path) if qualification_path.is_file() else None
    )
    if (
        qualification
        and qualification.get("Fingerprint") == key
        and (
            qualification.get("Passed")
            or (
                getattr(args, "coverage_only", False)
                and qualification.get("GeometryComplete")
            )
            or (
                getattr(args, "probe_study_only", False)
                and qualification.get("ProbeStudyComplete")
            )
        )
        and (library_path.is_file() or getattr(args, "probe_study_only", False))
        and not args.force
    ):
        print(f"Reusing completed spatial cache {root}")
        return (
            None if getattr(args, "probe_study_only", False) else library_path,
            qualification_path,
        )

    root.mkdir(parents=True, exist_ok=True)
    coupon_path = root / "coupon.json"
    generated_coupon = copy.deepcopy(coupon)
    canonical_geometry = spec["Coupon"]["Geometry"]
    if "PlanViewBoundary" in canonical_geometry:
        generated_coupon["Geometry"]["PlanViewBoundary"] = canonical_geometry[
            "PlanViewBoundary"
        ]
        generated_coupon["Geometry"]["MaskRegularization"] = canonical_geometry[
            "MaskRegularization"
        ]
    write_json(coupon_path, generated_coupon)
    write_json(root / "coupon-spec.json", spec)
    common_generator = [
        sys.executable,
        SPATIAL_GENERATOR,
        coupon_path,
        "--output",
        root,
        "--radius",
        args.matching_radius,
        "--metal-thickness",
        parameters["metal_thickness"],
        "--overetch-depth",
        parameters["overetch"],
        "--sidewall-angle",
        parameters["sidewall_angle"],
        "--top-rounding",
        parameters["top_radius"],
        "--trench-rounding",
        parameters["bottom_radius"],
        "--ring-size",
        args.spatial_ring_size,
        "--order",
        max(args.orders),
        "--model-name",
        f"{coupon['Topology'].lower()}-{key[:12]}",
        *material_options(parameters),
    ]
    try:
        run([*common_generator, "--signature-only"])
    except subprocess.CalledProcessError:
        failure_path = root / "generation-failure.json"
        if failure_path.is_file():
            failure = load_json(failure_path)
            raise RuntimeError(failure.get("Reason", "Spatial geometry failed")) from None
        raise
    signature = root / "mesh-signature.csv"
    mask = root / "plan-view-mask.csv"
    boundary = root / "plan-view-boundary.csv"
    _, meshes = generate_spatial_meshes(
        root,
        signature,
        mask,
        boundary,
        args,
        parameters,
        spec,
        1.0,
    )
    run(
        [
            *common_generator,
            "--thin-mesh",
            meshes["thin"],
            "--fabricated-mesh",
            meshes["fabricated"],
        ]
    )
    if getattr(args, "coverage_only", False):
        for name in ("spatial_thin", "spatial_fabricated"):
            run(palace_command(args, root / f"{name}.json"), cwd=root)
        interface_count = len(
            {
                (int(interface.get("Slot", 0)), interface["Type"])
                for interface in coupon.get("Interfaces", [])
            }
        )
        pad_missing_surface_response_matrices(
            root, interface_count, args.matching_radius
        )
        stamp_library_boundary_conditions(library_path, [coupon])
        stamp_geometry_coverage_only(library_path)
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "spatial",
            "Topology": coupon["Topology"],
            "Library": str(library_path),
            "GeometryComplete": True,
            "Passed": False,
            "Reason": "Geometry-coverage generation does not perform convergence qualification",
        }
        write_json(qualification_path, qualification)
        return library_path, qualification_path

    convergence_root = root / "probe-convergence"
    convergence_code, convergence_report, convergence = run_probe_convergence(
        root, convergence_root, args
    )
    if convergence_code != 0 or not convergence.get("Passed", False):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "spatial",
            "Topology": coupon["Topology"],
            "Library": str(library_path),
            "HeldoutReport": None,
            "ConvergenceReport": str(convergence_report),
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        raise RuntimeError(f"Spatial probe convergence failed: {qualification_path}")

    mesh_calibrations = []
    for factor in spec["Mesh"]["HRefinementFactors"]:
        if math.isclose(factor, 1.0):
            calibration = convergence_root / f"p{max(args.orders)}"
        else:
            mesh_root, factor_meshes = generate_spatial_meshes(
                root,
                signature,
                mask,
                boundary,
                args,
                parameters,
                spec,
                factor,
            )
            calibration = prepare_probe_mesh_calibration(
                root, mesh_root / "calibration", factor_meshes
            )
        mesh_calibrations.append((f"h-{slug(factor)}", calibration))
    mesh_convergence_root = root / "mesh-convergence"
    (
        mesh_convergence_code,
        mesh_convergence_report,
        mesh_convergence,
    ) = run_mesh_convergence(mesh_calibrations, mesh_convergence_root, args)
    if (
        mesh_convergence_code != 0
        or not mesh_convergence.get("Passed", False)
    ):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "spatial",
            "Topology": coupon["Topology"],
            "Library": str(library_path),
            "HeldoutReport": None,
            "ConvergenceReport": str(convergence_report),
            "MeshConvergenceReport": str(mesh_convergence_report),
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        raise RuntimeError(
            f"Spatial mesh convergence failed: {qualification_path}"
        )
    if getattr(args, "probe_study_only", False):
        qualification = {
            "Version": 1,
            "Fingerprint": key,
            "Family": "spatial",
            "Topology": coupon["Topology"],
            "Library": None,
            "HeldoutReport": None,
            "ConvergenceReport": str(convergence_report),
            "MeshConvergenceReport": str(mesh_convergence_report),
            "ProbeStudyComplete": True,
            "Passed": False,
        }
        write_json(qualification_path, qualification)
        return None, qualification_path

    if spatial_complete(root) and not args.force:
        print(f"Reusing completed spatial responses in {root}")
    else:
        for name in (
            "spatial_thin",
            "spatial_fabricated",
            "heldout_spatial_thin",
            "heldout_spatial_fabricated",
        ):
            run(palace_command(args, root / f"{name}.json"), cwd=root)
        if not spatial_complete(root):
            raise RuntimeError(f"Palace did not complete spatial responses in {root}")

    heldout_report = root / "heldout-qualification.json"
    heldout_code = run(
        [
            sys.executable,
            STRAIGHT_QUALIFIER,
            root,
            "--output",
            heldout_report,
            "--max-heldout-error",
            args.max_heldout_error,
        ],
        check=False,
    )
    heldout = load_json(heldout_report)
    passed = heldout_code == 0 and heldout.get("Passed", False)
    qualification = {
        "Version": 1,
        "Fingerprint": key,
        "Family": "spatial",
        "Topology": coupon["Topology"],
        "Library": str(library_path),
        "HeldoutReport": str(heldout_report),
        "ConvergenceReport": str(convergence_report),
        "MeshConvergenceReport": str(mesh_convergence_report),
        "Passed": passed,
    }
    write_json(qualification_path, qualification)
    if not passed:
        raise RuntimeError(f"Spatial qualification failed: {qualification_path}")
    stamp_library_boundary_conditions(library_path, [coupon])
    return library_path, qualification_path


def execute(plan, library_path, args):
    parameters = process_parameters(load_json(library_path))
    cache = args.cache.expanduser().resolve()
    cache.mkdir(parents=True, exist_ok=True)
    missing = [
        coupon for coupon in plan["Coupons"] if coupon["CoverageStatus"] == "Missing"
    ]
    unsupported = [
        coupon
        for coupon in missing
        if coupon["Preparation"]["Method"] == "Unsupported"
    ]
    # Families left to another builder (--skip-methods, e.g. SpatialCoupon: the prism-tube
    # library of coupon_library.py build / qualify) are recorded, never failures.
    skip_methods = set(getattr(args, "skip_methods", None) or [])
    skipped = [
        {"Ids": [coupon["Id"]], "Topologies": [coupon["Topology"]], "Method": coupon["Preparation"]["Method"],
         "Reason": "left to another builder (--skip-methods)"}
        for coupon in missing
        if coupon["Preparation"]["Method"] in skip_methods
    ]
    missing = [coupon for coupon in missing if coupon["Preparation"]["Method"] not in skip_methods]
    libraries = [library_path]
    qualifications = []
    failures = [
        {
            "Ids": [coupon["Id"]],
            "Topologies": [coupon["Topology"]],
            "Method": "Unsupported",
            "Reason": coupon["Preparation"]["Reason"],
        }
        for coupon in unsupported
    ]

    # The coupons whose builder produced a library in THIS execution: the only ones that stamp
    # a Signature on the merged library (OI-PLANNER-SKIPPED-SIGNATURE-STAMP, decisions 1117 /
    # 1137: never a Skipped, Unsupported or failed coupon).
    built = []

    def evaluate(task):
        coupons, method, builder = task
        try:
            generated, qualification = builder()
            return generated, qualification, None, coupons
        except Exception as error:
            failure = {
                "Ids": [coupon["Id"] for coupon in coupons],
                "Topologies": sorted(
                    {coupon["Topology"] for coupon in coupons}
                ),
                "Method": method,
                "Reason": f"{type(error).__name__}: {error}",
            }
            return None, None, failure, coupons

    def collect(result):
        generated, qualification, failure, coupons = result
        if failure:
            failures.append(failure)
            print(
                "Coupon qualification failed for "
                + ", ".join(failure["Ids"])
                + ": "
                + failure["Reason"],
                file=sys.stderr,
                flush=True,
            )
            return
        if generated is not None:
            libraries.append(generated)
            built.extend(coupons)
        qualifications.append(qualification)

    # Build the straight-edge family first. Besides normally being inexpensive, this
    # avoids concurrent first-use Julia package precompilation before independent corner
    # and spatial coupon workers start.
    straight = [
        coupon
        for coupon in missing
        if coupon["Preparation"]["Method"] == "StraightEdgeBuilder"
    ]
    if straight:
        collect(
            evaluate(
                (
                    straight,
                    "StraightEdgeBuilder",
                    lambda: with_heldout_escalation(
                        lambda a: build_straight(straight, a, parameters, cache), args, "straight-edge family"
                    ),
                )
            )
        )

    tasks = []
    for coupon in missing:
        method = coupon["Preparation"]["Method"]
        if method == "ParallelClusterCoupon":
            builder = lambda coupon=coupon: with_heldout_escalation(
                lambda a, coupon=coupon: build_parallel_cluster(coupon, a, parameters, cache), args, coupon["Id"]
            )
        elif method == "CornerCoupon":
            builder = lambda coupon=coupon: build_corner(
                coupon, args, parameters, cache
            )
        elif method == "SpatialCoupon":
            builder = lambda coupon=coupon: build_spatial(
                coupon, args, parameters, cache
            )
        else:
            continue
        tasks.append(([coupon], method, builder))

    coupon_jobs = getattr(args, "coupon_jobs", 1)
    if coupon_jobs == 1:
        results = map(evaluate, tasks)
        for result in results:
            collect(result)
    else:
        print(
            f"Executing {len(tasks)} independent coupon families with "
            f"{coupon_jobs} concurrent workers",
            flush=True,
        )
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=coupon_jobs
        ) as executor:
            # executor.map preserves deterministic task order for library combination even
            # though subprocesses run concurrently.
            for result in executor.map(evaluate, tasks):
                collect(result)

    destination = args.output / "library"
    run(
        [
            sys.executable,
            LIBRARY_COMBINER,
            "--output",
            destination,
            "--name",
            args.name,
            *(["--allow-empty"] if len(libraries) == 1 else []),
            *libraries,
        ]
    )
    coverage_only = getattr(args, "coverage_only", False)
    probe_study_only = getattr(args, "probe_study_only", False)
    stamped = 0
    if not probe_study_only and (destination / "process-library.json").is_file():
        stamped = stamp_library_signatures(destination / "process-library.json", built_coupons_for_stamping(built, skipped, failures))
    geometry_complete = not failures and not probe_study_only
    probe_study_complete = not failures and probe_study_only
    qualification_manifest = {
        "Version": 1,
        "SourceManifest": plan["SourceManifest"],
        "SourceLibrary": str(library_path),
        "GeneratedLibraries": [str(path) for path in libraries[1:]],
        "QualificationReports": [str(path) for path in qualifications],
        "Failures": failures,
        "Skipped": skipped,
        "SignatureModels": stamped,
        "HeldOutEscalation": {"Rule": HELDOUT_ESCALATION_RULE, "Enabled": getattr(args, "heldout_escalation", True),
                              "Factor": HELDOUT_ESCALATION_FACTOR, "Families": list(HELDOUT_ESCALATION_FAMILIES),
                              "NotApplicable": list(HELDOUT_ESCALATION_NOT_APPLICABLE)},
        "GeometryComplete": geometry_complete,
        "ProbeStudyComplete": probe_study_complete,
        "QualificationRequired": not coverage_only and not probe_study_only,
        "BoundaryLawPhysics": {
            "Complete": not any(uses_finite_impedance(coupon) for coupon in missing),
            "UnqualifiedCoupons": [
                coupon["Id"] for coupon in missing if uses_finite_impedance(coupon)
            ],
        },
        "Passed": (geometry_complete and not coverage_only) or probe_study_complete,
    }
    qualification_path = destination / "qualification-manifest.json"
    write_json(qualification_path, qualification_manifest)
    plan["Execution"] = {
        "Complete": geometry_complete and not coverage_only,
        "GeometryComplete": geometry_complete,
        "ProbeStudyComplete": probe_study_complete,
        "QualificationRequired": not coverage_only and not probe_study_only,
        "Library": str(destination / "process-library.json"),
        "QualificationManifest": str(qualification_path),
        "Failures": failures,
        "Skipped": skipped,
    }
    if probe_study_only:
        return probe_study_complete
    return geometry_complete if coverage_only else not failures


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument(
        "--include-matched",
        action="store_true",
        help="Include exact and interpolated requirements as validation coupons",
    )
    parser.add_argument("--execute", action="store_true")
    parser.add_argument(
        "--coverage-only",
        action="store_true",
        help=(
            "Generate every required response model without convergence or held-out "
            "qualification; entries are explicitly marked GeometryCoverageOnly/Unqualified"
        ),
    )
    parser.add_argument(
        "--probe-study-only",
        action="store_true",
        help=(
            "Run FEM-order and mesh-resolution probes for every coupon, but stop before "
            "assembling full corner/spatial response matrices"
        ),
    )
    parser.add_argument(
        "--skip-methods",
        nargs="+",
        default=[],
        choices=("StraightEdgeBuilder", "ParallelClusterCoupon", "CornerCoupon", "SpatialCoupon"),
        help=(
            "coupon families left to another builder (recorded under Skipped, never a failure): "
            "SpatialCoupon when the prism-tube library (coupon_library.py build / qualify) builds them"
        ),
    )
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--palace", type=Path)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--julia-project", type=Path)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument(
        "--coupon-jobs",
        type=int,
        default=1,
        help=(
            "Number of independent corner/spatial coupon families to execute concurrently; "
            "total MPI ranks can reach --coupon-jobs times --ranks"
        ),
    )
    parser.add_argument("--name", default="qualified-fabrication-process")
    parser.add_argument("--orders", type=int, nargs="+", default=[2, 3])
    parser.add_argument("--mesh-order", type=int, default=1)
    parser.add_argument("--straight-lc-fine", type=float, default=0.002)
    parser.add_argument("--straight-lc-far", type=float, default=0.05)
    parser.add_argument("--cluster-lc-fine", type=float, default=0.002)
    parser.add_argument("--cluster-lc-far", type=float, default=0.05)
    parser.add_argument(
        "--cluster-h-factors",
        type=float,
        nargs="+",
        default=[2.0, 1.0],
        help=(
            "Coarse-to-fine multipliers on --cluster-lc-fine used for the "
            "mesh-resolution convergence gate"
        ),
    )
    parser.add_argument(
        "--straight-edge-distances",
        type=float,
        nargs="+",
        default=None,
        help=(
            "localized-energy radii (um) of the StraightEdgeBuilder interfaces, ending at 0.2 "
            "(default 0.2 alone); the 3D thin ring radii give the 2D MA per shell (decision 66 part C)"
        ),
    )
    parser.add_argument("--corner-lc-fine", type=float, default=0.02)
    parser.add_argument("--corner-lc-far", type=float, default=0.3)
    parser.add_argument(
        "--corner-h-factors",
        type=float,
        nargs="+",
        default=[2.0, 1.0],
        help=(
            "Coarse-to-fine multipliers on --corner-lc-fine used for the "
            "mesh-resolution convergence gate"
        ),
    )
    parser.add_argument(
        "--spatial-mesher",
        choices=("field", "swept"),
        default="field",
        help=(
            "Spatial coupon mesher. 'field' is the existing experimental "
            "unstructured mode; 'swept' selects the explicit swept-prism path"
        ),
    )
    parser.add_argument("--spatial-lc-fine", type=float, default=0.02)
    parser.add_argument(
        "--spatial-lc-tangent",
        type=float,
        default=0.0,
        help=(
            "Independent longitudinal size. The field mesher accepts 0 for its "
            "isotropic fallback; the swept mesher requires a positive value"
        ),
    )
    parser.add_argument("--spatial-lc-far", type=float, default=0.3)
    parser.add_argument(
        "--spatial-process-core-width",
        type=float,
        default=0.0,
        help=(
            "Width over which the mesh grades from fine to far size; 0 selects "
            "max(2*metal thickness, 4*overetch, 8*lc-fine)"
        ),
    )
    parser.add_argument(
        "--spatial-process-fine-width",
        type=float,
        default=0.0,
        help="Optional flat inner band retained at lc-fine; 0 grades immediately",
    )
    parser.add_argument(
        "--spatial-process-grading-power",
        type=float,
        default=1.7,
        help="Power-law exponent for edge-normal grading from fine to far size",
    )
    parser.add_argument(
        "--spatial-normal-growth-ratio",
        type=float,
        default=1.4,
        help="Geometric growth ratio for explicit swept normal coordinates",
    )
    parser.add_argument(
        "--spatial-max-nodes",
        type=int,
        help=(
            "Spatial mesh node budget; defaults to 500k for field mode and "
            "3M for explicit swept meshes"
        ),
    )
    parser.add_argument("--spatial-max-elements", type=int, default=2_000_000)
    parser.add_argument(
        "--spatial-h-factors",
        type=float,
        nargs="+",
        default=[2.0, 1.0],
        help=(
            "Coarse-to-fine multipliers on --spatial-lc-fine used for the "
            "mesh-resolution convergence gate"
        ),
    )
    parser.add_argument("--min-process-feature-elements", type=float, default=2.0)
    parser.add_argument("--basis-size", type=int, default=96)
    parser.add_argument("--samples", type=int, default=1200)
    parser.add_argument("--ring-size", type=int, default=8)
    parser.add_argument(
        "--corner-trace-basis",
        choices=sorted((*CORNER_TRACE_BASIS_RING_SIZES, *CORNER_TRACE_BASIS_SELECTORS)),
        default=None,
        help="trace basis rule of the corner coupons (generate_corner_response.py "
        "--trace-basis; the ring size is the rule's, --ring-size does not apply to corners); "
        f"default: {CORNER_TRACE_BASIS_DEFAULT} on every corner, sharp or rounded (the box "
        "basis does not depend on CornerRadius); "
        f"{CORNER_TRACE_BASIS_CONCAVE_GRADED_SELECTOR} = the convexity-conditional selector "
        "(the metal-graded rule on a concave node, the default on a convex node)",
    )
    parser.add_argument("--spatial-ring-size", type=int, default=16)
    parser.add_argument("--coupon-depth", type=float, default=1055.0)
    parser.add_argument("--edge-offset-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--local-pieces", action="store_true",
                        help="R0-b (identification fix round DESIGN A.3 (b); decisions 1064 (4) / 1096 (8) / 1136 (4)): the piece form - a "
                             "coupon's conductor labels are its LOCAL metal pieces (one reference per piece; the parallel-cluster "
                             "signature CSV and the same-conductor-gap mesher carry Conductor := piece; Edges[].LocalPiece / LocalPieces "
                             "on the library entries; the straight / parallel-cluster cache keys carry the form and the model Name "
                             "the -pieces<P> suffix). THE DEFAULT since decision 1136 (4) for every coupon whose cross-section / plan "
                             "view holds a label with >= 2 pieces (Preparation.LocalPieces in the plan); this flag names the default "
                             "explicitly and is refused with --one-reference")
    parser.add_argument("--one-reference", action="store_true",
                        help="opt out of the piece-form default (decision 1136 (4)): every coupon in the one-reference form (the ES "
                             "record form of v3-j; every output byte-identical to the pre-1136 planner)")
    parser.add_argument("--max-heldout-error", type=float, default=10.0)
    parser.add_argument("--no-heldout-escalation", action="store_true",
                        help=f"switch the held-out escalation rule OFF ({HELDOUT_ESCALATION_RULE}; default ON for "
                             f"{' / '.join(HELDOUT_ESCALATION_FAMILIES)})")
    parser.add_argument("--max-fabricated-matrix-change", type=float, default=5.0)
    parser.add_argument("--max-fabricated-energy-change", type=float, default=10.0)
    parser.add_argument("--max-domain-defect-change", type=float, default=5.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    args.heldout_escalation = not args.no_heldout_escalation
    args.output = args.output.expanduser().resolve()
    args.cache = (
        args.cache.expanduser().resolve()
        if args.cache
        else args.output / "cache"
    )
    args.orders = sorted(set(args.orders))
    if args.coverage_only and args.probe_study_only:
        parser.error("--coverage-only and --probe-study-only are mutually exclusive")
    if args.local_pieces and args.one_reference:
        parser.error("--local-pieces (the default form) and --one-reference (the opt-out) are mutually exclusive")
    if args.execute and args.palace is None:
        parser.error("--palace is required with --execute")
    if args.ranks <= 0:
        parser.error("--ranks must be positive")
    if args.coupon_jobs <= 0:
        parser.error("--coupon-jobs must be positive")
    if any(order <= 0 for order in args.orders) or (
        not args.coverage_only and len(args.orders) < 2
    ):
        parser.error(
            "--orders requires positive FEM orders and at least two unless "
            "--coverage-only is selected"
        )
    if args.basis_size <= 0 or args.basis_size % 2:
        parser.error("--basis-size must be a positive even number")
    if args.ring_size < 8 or args.ring_size % 8:
        parser.error("--ring-size must be a multiple of eight and at least eight")
    if args.spatial_ring_size < 8 or args.spatial_ring_size % 4:
        parser.error(
            "--spatial-ring-size must be a multiple of four and at least eight"
        )
    if args.samples < args.basis_size:
        parser.error("--samples must be at least --basis-size")
    if args.edge_offset_tolerance < 0.0:
        parser.error("--edge-offset-tolerance must be nonnegative")
    if args.min_process_feature_elements <= 0.0:
        parser.error("--min-process-feature-elements must be positive")
    if args.spatial_lc_fine <= 0.0:
        parser.error("--spatial-lc-fine must be positive")
    if args.spatial_lc_far < args.spatial_lc_fine:
        parser.error("--spatial-lc-far must not be smaller than fine size")
    if args.spatial_lc_tangent < 0.0:
        parser.error("--spatial-lc-tangent must be nonnegative")
    if args.spatial_mesher == "swept":
        if args.spatial_lc_tangent <= 0.0:
            parser.error("the swept spatial mesher requires --spatial-lc-tangent > 0")
    elif args.spatial_lc_tangent and not (
        args.spatial_lc_fine <= args.spatial_lc_tangent <= args.spatial_lc_far
    ):
        parser.error("--spatial-lc-tangent must lie between fine and far sizes")
    if args.spatial_process_core_width < 0.0:
        parser.error("--spatial-process-core-width must be nonnegative")
    if args.spatial_process_fine_width < 0.0:
        parser.error("--spatial-process-fine-width must be nonnegative")
    if args.spatial_process_grading_power <= 0.0:
        parser.error("--spatial-process-grading-power must be positive")
    if args.spatial_normal_growth_ratio <= 1.0:
        parser.error("--spatial-normal-growth-ratio must exceed one")
    if args.spatial_max_nodes is None:
        args.spatial_max_nodes = (
            3_000_000 if args.spatial_mesher == "swept" else 500_000
        )
    if args.spatial_max_nodes <= 0 or args.spatial_max_elements <= 0:
        parser.error("spatial mesh budgets must be positive")
    args.cluster_h_factors = sorted(
        set(args.cluster_h_factors), reverse=True
    )
    if (
        len(args.cluster_h_factors) < 2
        or any(factor < 1.0 for factor in args.cluster_h_factors)
        or not math.isclose(args.cluster_h_factors[-1], 1.0)
    ):
        parser.error(
            "--cluster-h-factors requires at least two values >= 1 ending at 1"
        )
    args.corner_h_factors = sorted(
        set(args.corner_h_factors), reverse=True
    )
    if (
        len(args.corner_h_factors) < 2
        or any(factor < 1.0 for factor in args.corner_h_factors)
        or not math.isclose(args.corner_h_factors[-1], 1.0)
    ):
        parser.error(
            "--corner-h-factors requires at least two values >= 1 ending at 1"
        )
    args.spatial_h_factors = sorted(
        set(args.spatial_h_factors), reverse=True
    )
    if (
        len(args.spatial_h_factors) < 2
        or any(factor < 1.0 for factor in args.spatial_h_factors)
        or not math.isclose(args.spatial_h_factors[-1], 1.0)
    ):
        parser.error(
            "--spatial-h-factors requires at least two values >= 1 ending at 1"
        )
    if (
        args.spatial_mesher == "field"
        and args.spatial_lc_tangent
        and args.spatial_h_factors[0] * args.spatial_lc_fine
        > args.spatial_lc_tangent
    ):
        parser.error(
            "the field spatial mesher requires --spatial-lc-tangent to be at "
            "least the coarsest h-refined normal size"
        )
    return args


def main():
    args = parse_args()
    manifest_path = args.manifest.expanduser().resolve()
    manifest = load_json(manifest_path)
    # Version 2 (geometry identification, SURFACE-RESPONSE-IDENTIFICATION.md (d)) keeps every
    # version-1 field; its Requirements are derived from the identification features.
    if manifest.get("Version") not in (1, 2) or not isinstance(
        manifest.get("Requirements"), list
    ):
        raise ValueError(f"{manifest_path} is not a version-1 or version-2 requirements manifest")
    args.matching_radius = float(manifest["Library"]["MatchingRadius"])

    library_path = Path(manifest["Library"]["Path"]).expanduser()
    if not library_path.is_absolute():
        library_path = (manifest_path.parent / library_path).resolve()
    library = load_json(library_path)
    plan = plan_from_manifest(
        manifest_path,
        manifest,
        library_path,
        library,
        args.include_matched,
        spatial_mesh_path(args),
        one_reference=args.one_reference,
    )
    args.output.mkdir(parents=True, exist_ok=True)
    destination = args.output / "coupon-plan.json"
    write_json(destination, plan)
    print(destination)
    if args.execute:
        complete = False
        try:
            complete = execute(plan, library_path, args)
        finally:
            write_json(destination, plan)
        print(plan["Execution"]["Library"])
        if not complete:
            description = (
                "geometry generation"
                if args.coverage_only
                else "probe study"
                if args.probe_study_only
                else "qualification"
            )
            print(
                f"One or more coupon requirements failed {description}; "
                "see coupon-plan.json",
                file=sys.stderr,
            )
            raise SystemExit(1)
    elif plan["Summary"]["Unsupported"]:
        print(
            f"Warning: {plan['Summary']['Unsupported']} coupon requirement(s) need "
            "a generic spatial coupon generator"
        )


if __name__ == "__main__":
    main()
