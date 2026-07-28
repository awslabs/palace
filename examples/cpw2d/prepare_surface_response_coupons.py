#!/usr/bin/env python3

"""Plan, build, qualify, cache, and merge surface-response coupons."""

import argparse
import copy
import hashlib
import json
import math
import shlex
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
STRAIGHT_BUILDER = ROOT / "build_surface_process_library.py"
STRAIGHT_QUALIFIER = ROOT / "qualify_surface_response.py"
CLUSTER_MESH = ROOT / "mesh" / "mesh_edge_cluster_coupon.jl"
CLUSTER_GENERATOR = ROOT / "generate_edge_cluster_response.py"
CORNER_ROOT = ROOT.parent / "cpw3d_surface" / "corner_coupon"
CORNER_MESH = CORNER_ROOT / "mesh_corner_coupon.jl"
CORNER_GENERATOR = CORNER_ROOT / "generate_corner_response.py"
CORNER_FINALIZER = CORNER_ROOT / "finalize_corner_response.py"
CORNER_CONVERGENCE = CORNER_ROOT / "run_probe_convergence.py"
LIBRARY_COMBINER = CORNER_ROOT / "combine_process_libraries.py"
SPATIAL_ROOT = ROOT.parent / "cpw3d_surface" / "spatial_coupon"
SPATIAL_MESH = SPATIAL_ROOT / "mesh_spatial_coupon.jl"
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


def run(command, check=True):
    command = [str(value) for value in command]
    print("+ " + shlex.join(command), flush=True)
    return subprocess.run(command, check=check).returncode


def slug(value):
    return (
        f"{value:.12g}".replace("-", "m").replace(".", "p")
        if isinstance(value, float)
        else str(value).lower()
    )


def coupon_id(requirement):
    topology = requirement["Topology"]
    geometry = requirement.get("Geometry", {})
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
            {
                "Topology": topology,
                "Geometry": geometry,
                "BoundaryCondition": requirement.get("BoundaryCondition", {}),
                "Interfaces": requirement.get("Interfaces", []),
            },
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    ).hexdigest()[:12]
    return "_".join(parts + [digest])


def preparation(requirement):
    topology = requirement["Topology"]
    geometry = requirement.get("Geometry", {})
    boundary_conditions = [requirement.get("BoundaryCondition", {})]
    for entries in (geometry.get("Edges", []), geometry.get("Arms", [])):
        boundary_conditions.extend(
            entry.get("BoundaryCondition", {})
            for entry in entries
            if isinstance(entry, dict)
        )
    finite_impedance = next(
        (
            condition.get("Type", "unknown")
            for condition in boundary_conditions
            if condition and condition.get("Type") != "PEC"
        ),
        None,
    )
    if finite_impedance:
        return {
            "Method": "Unsupported",
            "Reason": (
                "Automatic coupon generation is not qualified for "
                f"{finite_impedance} metal"
            ),
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
        if math.isclose(angle, 90.0, abs_tol=1.0e-8):
            return {
                "Method": "CornerCoupon",
                "MeshGenerator": str(CORNER_MESH.relative_to(ROOT.parent.parent)),
                "ResponseGenerator": str(
                    CORNER_GENERATOR.relative_to(ROOT.parent.parent)
                ),
            }
        return {
            "Method": "Unsupported",
            "Reason": "The current corner coupon generator supports 90-degree corners",
        }
    if topology in SPATIAL_TOPOLOGIES:
        return {
            "Method": "SpatialCoupon",
            "MeshGenerator": str(
                SPATIAL_MESH.relative_to(ROOT.parent.parent)
            ),
            "ResponseGenerator": str(
                SPATIAL_GENERATOR.relative_to(ROOT.parent.parent)
            ),
        }
    return {"Method": "Unsupported", "Reason": "Unknown coupon topology"}


def plan_from_manifest(manifest_path, manifest, library_path, library, include_matched):
    coupons = []
    for requirement in manifest["Requirements"]:
        if not include_matched and requirement.get("Status") != "Missing":
            continue
        coupon = {
            "Id": coupon_id(requirement),
            "Topology": requirement["Topology"],
            "Geometry": requirement.get("Geometry", {}),
            "Interfaces": requirement.get("Interfaces", []),
            "BoundaryCondition": requirement.get("BoundaryCondition", {}),
            "DeviceOccurrences": requirement.get("Count", 0),
            "DeviceEdgeLength": requirement.get("TotalEdgeLength", 0.0),
            "CoverageStatus": requirement.get("Status"),
            "Preparation": preparation(requirement),
        }
        if "SelectedModels" in requirement:
            coupon["SelectedModels"] = requirement["SelectedModels"]
        if "Reason" in requirement:
            coupon["CoverageReason"] = requirement["Reason"]
        coupons.append(coupon)
    coupons.sort(key=lambda coupon: coupon["Id"])
    return {
        "Version": 2,
        "SourceManifest": str(manifest_path),
        "ProcessLibrary": str(library_path),
        "Fabrication": library.get("Fabrication"),
        "MatchingRadius": manifest["Library"]["MatchingRadius"],
        "Coupons": coupons,
        "Summary": {
            "CouponCount": len(coupons),
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
                if (
                    cross(
                        subtract(first, vertex),
                        subtract(second, vertex),
                    )
                    == (0, 0, 0)
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


def require_pec(coupon):
    condition = coupon.get("BoundaryCondition", {})
    if condition.get("Type") != "PEC":
        raise ValueError(
            f"{coupon['Id']} uses {condition.get('Type', 'unknown')} metal; "
            "preflight currently lacks the dimensional boundary-law metadata needed "
            "to stamp a generated coupon safely"
        )


def require_spatial_pec(coupon):
    require_pec(coupon)
    for edge in coupon.get("Geometry", {}).get("Edges", []):
        condition = edge.get("BoundaryCondition", coupon["BoundaryCondition"])
        if condition.get("Type") != "PEC":
            raise ValueError(
                f"{coupon['Id']} contains a non-PEC spatial edge; automatic "
                "finite-impedance coupon calibration is not yet qualified"
            )


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
    if edges[-1]["Offset"] > 2.0 * matching_radius + tolerance:
        raise ValueError(
            f"{coupon['Id']} spans more than the 2R interaction diameter"
        )
    return edges


def write_cluster_signature(root, edges):
    data = {"Edges": edges}
    json_path = root / "cluster-signature.json"
    write_json(json_path, data)
    csv_path = root / "cluster-signature.csv"
    lines = ["Offset,GapDirection,Conductor"]
    lines.extend(
        f"{edge['Offset']:.17g},{edge['GapDirection']},{edge['Conductor']}"
        for edge in edges
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


def build_straight(coupons, args, parameters, cache):
    for coupon in coupons:
        require_pec(coupon)
    resolution = process_resolution(
        parameters, args.straight_lc_fine, args.min_process_feature_elements
    )
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
        and qualification.get("Passed")
        and final_library.is_file()
        and not args.force
    ):
        print(f"Reusing qualified straight-edge cache {root}")
        return final_library, qualification_path

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
        if args.julia_project:
            command.extend(["--julia-project", args.julia_project])
        if not isolated:
            command.append("--skip-isolated")
        for topology, option in PAIRED_TOPOLOGIES.items():
            command.append(f"--{option}")
            command.extend(separations[topology])
        run(command)

    previous = order_roots[args.orders[-2]]
    current = order_roots[args.orders[-1]]
    checks = []
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
    qualification = {
        "Version": 1,
        "Fingerprint": key,
        "Family": "straight",
        "Library": str(final_library),
        "Checks": checks,
        "Passed": passed,
    }
    write_json(qualification_path, qualification)
    if not passed:
        raise RuntimeError(f"Straight-edge qualification failed: {qualification_path}")
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


def build_parallel_cluster(coupon, args, parameters, cache):
    require_pec(coupon)
    resolution = process_resolution(
        parameters, args.cluster_lc_fine, args.min_process_feature_elements
    )
    edges = normalize_parallel_edges(coupon, args.matching_radius)
    spec = {
        "Version": 1,
        "Family": "parallel-edge-cluster",
        "Edges": edges,
        "Fabrication": parameters,
        "MatchingRadius": args.matching_radius,
        "Orders": args.orders,
        "Mesh": {
            "FineSize": args.cluster_lc_fine,
            "FarSize": args.cluster_lc_far,
            "Order": args.mesh_order,
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
        and qualification.get("Passed")
        and final_library.is_file()
        and not args.force
    ):
        print(f"Reusing qualified parallel-cluster cache {root}")
        return final_library, qualification_path

    root.mkdir(parents=True, exist_ok=True)
    write_json(root / "coupon-spec.json", spec)
    json_signature, csv_signature = write_cluster_signature(root, edges)
    meshes = {}
    for kind in ("thin", "fabricated"):
        mesh = root / f"cluster_{kind}.msh"
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
                args.cluster_lc_fine,
                "--lc-far",
                args.cluster_lc_far,
                "--mesh-order",
                args.mesh_order,
            ]
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
            "--library-name",
            f"{args.name}-parallel-cluster-p{final_order}",
            "--model-name",
            f"parallel-edge-cluster-{key[:12]}",
            *material_options(parameters),
        ]
    )
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

    if cluster_complete(final_root) and not args.force:
        print(f"Reusing completed p{final_order} parallel-cluster responses")
    else:
        for name in (
            "cluster_thin",
            "cluster_fabricated",
            "heldout_cluster_thin",
            "heldout_cluster_fabricated",
        ):
            run(palace_command(args, final_root / f"{name}.json"))
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
    qualification = {
        "Version": 1,
        "Fingerprint": key,
        "Family": "parallel-edge-cluster",
        "Library": str(final_library),
        "ResponseReport": str(response_report),
        "ConvergenceReport": str(convergence_report),
        "Passed": passed,
    }
    write_json(qualification_path, qualification)
    if not passed:
        raise RuntimeError(
            f"Parallel-cluster qualification failed: {qualification_path}"
        )
    return final_library, qualification_path


def build_corner(coupon, args, parameters, cache):
    require_pec(coupon)
    resolution = process_resolution(
        parameters, args.corner_lc_fine, args.min_process_feature_elements
    )
    geometry = coupon["Geometry"]
    angle = float(geometry["AngleDegrees"])
    if not math.isclose(angle, 90.0, abs_tol=1.0e-8):
        raise ValueError(f"{coupon['Id']} requests unsupported corner angle {angle:g}")
    topology = "convex" if coupon["Topology"] == "ConvexCorner" else "concave"
    corner_radius = float(geometry.get("CornerRadius", 0.0))
    spec = {
        "Version": 1,
        "Family": "corner",
        "Topology": topology,
        "CornerRadius": corner_radius,
        "Fabrication": parameters,
        "MatchingRadius": args.matching_radius,
        "Orders": args.orders,
        "Mesh": {
            "FineSize": args.corner_lc_fine,
            "FarSize": args.corner_lc_far,
            "Order": args.mesh_order,
        },
        "ProcessResolution": resolution,
        "Response": {"RingSize": args.ring_size},
        "ToolFingerprint": tool_fingerprint(
            (
                CORNER_MESH,
                CORNER_GENERATOR,
                CORNER_FINALIZER,
                CORNER_CONVERGENCE,
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
        and qualification.get("Passed")
        and library_path.is_file()
        and not args.force
    ):
        print(f"Reusing qualified corner cache {root}")
        return library_path, qualification_path

    root.mkdir(parents=True, exist_ok=True)
    write_json(root / "coupon-spec.json", spec)
    meshes = {}
    for kind in ("thin", "fabricated"):
        mesh = root / f"corner_{kind}.msh"
        meshes[kind] = mesh
        run(
            [
                *julia_command(args),
                CORNER_MESH,
                f"{topology}-{kind}",
                mesh,
                "--radius",
                args.matching_radius,
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
                args.corner_lc_fine,
                "--lc-far",
                args.corner_lc_far,
                "--mesh-order",
                args.mesh_order,
            ]
        )
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
            "--corner-radius",
            corner_radius,
            "--ring-size",
            args.ring_size,
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

    for name in ("thin", "fabricated"):
        run(palace_command(args, root / f"{name}.json"))
    run([sys.executable, CORNER_FINALIZER, root])
    for name in ("heldout-thin", "heldout-fabricated"):
        run(palace_command(args, root / f"{name}.json"))
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
        "Passed": passed,
    }
    write_json(qualification_path, qualification)
    if not passed:
        raise RuntimeError(f"Corner qualification failed: {qualification_path}")
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


def spatial_spec(coupon, args, parameters):
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
            "FineSize": args.spatial_lc_fine,
            "FarSize": args.spatial_lc_far,
            "Order": args.mesh_order,
        },
        "ProcessResolution": resolution,
        "Response": {"RingSize": args.spatial_ring_size},
        "RegularizationPolicy": geometry.get("MaskRegularization"),
        "ToolFingerprint": tool_fingerprint(
            (
                SPATIAL_MESH,
                SPATIAL_GENERATOR,
                STRAIGHT_QUALIFIER,
                CORNER_CONVERGENCE,
                Path(__file__),
            )
        ),
    }


def build_spatial(coupon, args, parameters, cache):
    require_spatial_pec(coupon)
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
        and qualification.get("Passed")
        and library_path.is_file()
        and not args.force
    ):
        print(f"Reusing qualified spatial cache {root}")
        return library_path, qualification_path

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
    meshes = {}
    for kind in ("thin", "fabricated"):
        mesh = root / f"spatial_{kind}.msh"
        meshes[kind] = mesh
        if mesh.is_file() and not args.force:
            continue
        run(
            [
                *julia_command(args),
                SPATIAL_MESH,
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
                args.spatial_lc_fine,
                "--lc-far",
                args.spatial_lc_far,
                "--mesh-order",
                args.mesh_order,
            ]
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

    if spatial_complete(root) and not args.force:
        print(f"Reusing completed spatial responses in {root}")
    else:
        for name in (
            "spatial_thin",
            "spatial_fabricated",
            "heldout_spatial_thin",
            "heldout_spatial_fabricated",
        ):
            run(palace_command(args, root / f"{name}.json"))
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
        "Passed": passed,
    }
    write_json(qualification_path, qualification)
    if not passed:
        raise RuntimeError(f"Spatial qualification failed: {qualification_path}")
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

    def attempt(coupons, method, builder):
        try:
            generated, qualification = builder()
        except Exception as error:
            failure = {
                "Ids": [coupon["Id"] for coupon in coupons],
                "Topologies": sorted(
                    {coupon["Topology"] for coupon in coupons}
                ),
                "Method": method,
                "Reason": f"{type(error).__name__}: {error}",
            }
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
        libraries.append(generated)
        qualifications.append(qualification)

    straight = [
        coupon
        for coupon in missing
        if coupon["Preparation"]["Method"] == "StraightEdgeBuilder"
    ]
    if straight:
        attempt(
            straight,
            "StraightEdgeBuilder",
            lambda: build_straight(straight, args, parameters, cache),
        )
    for coupon in missing:
        if coupon["Preparation"]["Method"] != "ParallelClusterCoupon":
            continue
        attempt(
            [coupon],
            "ParallelClusterCoupon",
            lambda coupon=coupon: build_parallel_cluster(
                coupon, args, parameters, cache
            ),
        )
    for coupon in missing:
        if coupon["Preparation"]["Method"] != "CornerCoupon":
            continue
        attempt(
            [coupon],
            "CornerCoupon",
            lambda coupon=coupon: build_corner(
                coupon, args, parameters, cache
            ),
        )
    for coupon in missing:
        if coupon["Preparation"]["Method"] != "SpatialCoupon":
            continue
        attempt(
            [coupon],
            "SpatialCoupon",
            lambda coupon=coupon: build_spatial(
                coupon, args, parameters, cache
            ),
        )

    destination = args.output / "library"
    run(
        [
            sys.executable,
            LIBRARY_COMBINER,
            "--output",
            destination,
            "--name",
            args.name,
            *libraries,
        ]
    )
    qualification_manifest = {
        "Version": 1,
        "SourceManifest": plan["SourceManifest"],
        "SourceLibrary": str(library_path),
        "GeneratedLibraries": [str(path) for path in libraries[1:]],
        "QualificationReports": [str(path) for path in qualifications],
        "Failures": failures,
        "Passed": not failures,
    }
    qualification_path = destination / "qualification-manifest.json"
    write_json(qualification_path, qualification_manifest)
    plan["Execution"] = {
        "Complete": not failures,
        "Library": str(destination / "process-library.json"),
        "QualificationManifest": str(qualification_path),
        "Failures": failures,
    }
    return not failures


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
    parser.add_argument("--cache", type=Path)
    parser.add_argument("--palace", type=Path)
    parser.add_argument("--julia", default="julia")
    parser.add_argument("--julia-project", type=Path)
    parser.add_argument("--ranks", type=int, default=1)
    parser.add_argument("--name", default="qualified-fabrication-process")
    parser.add_argument("--orders", type=int, nargs="+", default=[2, 3])
    parser.add_argument("--mesh-order", type=int, default=1)
    parser.add_argument("--straight-lc-fine", type=float, default=0.002)
    parser.add_argument("--straight-lc-far", type=float, default=0.05)
    parser.add_argument("--cluster-lc-fine", type=float, default=0.002)
    parser.add_argument("--cluster-lc-far", type=float, default=0.05)
    parser.add_argument("--corner-lc-fine", type=float, default=0.02)
    parser.add_argument("--corner-lc-far", type=float, default=0.3)
    parser.add_argument("--spatial-lc-fine", type=float, default=0.02)
    parser.add_argument("--spatial-lc-far", type=float, default=0.3)
    parser.add_argument("--min-process-feature-elements", type=float, default=2.0)
    parser.add_argument("--basis-size", type=int, default=96)
    parser.add_argument("--samples", type=int, default=1200)
    parser.add_argument("--ring-size", type=int, default=8)
    parser.add_argument("--spatial-ring-size", type=int, default=16)
    parser.add_argument("--coupon-depth", type=float, default=1055.0)
    parser.add_argument("--edge-offset-tolerance", type=float, default=1.0e-3)
    parser.add_argument("--max-heldout-error", type=float, default=10.0)
    parser.add_argument("--max-fabricated-matrix-change", type=float, default=5.0)
    parser.add_argument("--max-fabricated-energy-change", type=float, default=10.0)
    parser.add_argument("--max-domain-defect-change", type=float, default=5.0)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    args.output = args.output.expanduser().resolve()
    args.cache = (
        args.cache.expanduser().resolve()
        if args.cache
        else args.output / "cache"
    )
    args.orders = sorted(set(args.orders))
    if args.execute and args.palace is None:
        parser.error("--palace is required with --execute")
    if args.ranks <= 0:
        parser.error("--ranks must be positive")
    if len(args.orders) < 2 or any(order <= 0 for order in args.orders):
        parser.error("--orders requires at least two positive FEM orders")
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
    return args


def main():
    args = parse_args()
    manifest_path = args.manifest.expanduser().resolve()
    manifest = load_json(manifest_path)
    if manifest.get("Version") != 1 or not isinstance(
        manifest.get("Requirements"), list
    ):
        raise ValueError(f"{manifest_path} is not a version-1 requirements manifest")
    args.matching_radius = float(manifest["Library"]["MatchingRadius"])

    library_path = Path(manifest["Library"]["Path"]).expanduser()
    if not library_path.is_absolute():
        library_path = (manifest_path.parent / library_path).resolve()
    library = load_json(library_path)
    plan = plan_from_manifest(
        manifest_path, manifest, library_path, library, args.include_matched
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
            print(
                "One or more coupon requirements failed qualification; "
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
