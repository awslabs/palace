#!/usr/bin/env python3

"""Generate response inputs for endpoint, junction, and spatial-edge coupons."""

import argparse
import json
import math
from pathlib import Path

import numpy as np


INTERFACE_DEFAULTS = {
    "SA": (0.002, 4.0),
    "MS": (0.0003, 11.47),
    "MA": (0.03, 10.0),
}
INTERFACE_LOSS = {"SA": 0.002, "MS": 0.0003, "MA": 0.03}


def load_json(path):
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError as error:
        raise ValueError(f"Invalid JSON in {path}: {error}") from error


def unit(values, name):
    vector = np.asarray(values, dtype=float)
    if vector.shape != (3,) or not np.all(np.isfinite(vector)):
        raise ValueError(f"{name} must be a finite three-vector")
    norm = np.linalg.norm(vector)
    if not math.isclose(norm, 1.0, rel_tol=0.0, abs_tol=1.0e-9):
        raise ValueError(f"{name} must be a unit vector")
    return vector / norm


def frame_from_geometry(topology, geometry):
    entries = (
        geometry.get("Edges", [])
        if topology == "SpatialEdgeCluster"
        else geometry.get("Arms", [])
    )
    if not entries:
        raise ValueError(f"{topology} has no complete edge geometry")
    normal = unit(entries[0]["ProcessNormal"], "ProcessNormal")
    gap = unit(entries[0]["GapDirection"], "GapDirection")
    if abs(np.dot(normal, gap)) > 1.0e-9:
        raise ValueError("ProcessNormal and GapDirection must be orthogonal")
    # Local z is the fabrication normal. Local x/y form a right-handed process plane.
    axis_x = gap
    axis_y = np.cross(normal, gap)
    axis_y /= np.linalg.norm(axis_y)
    return np.vstack((axis_x, axis_y, normal))


def normalize_geometry(coupon, radius):
    topology = coupon["Topology"]
    geometry = coupon.get("Geometry", {})
    frame = frame_from_geometry(topology, geometry)
    entries = (
        geometry["Edges"]
        if topology == "SpatialEdgeCluster"
        else geometry["Arms"]
    )
    edges = []
    for index, entry in enumerate(entries):
        point = np.asarray(entry.get("Point", [0.0, 0.0, 0.0]), dtype=float)
        gap = unit(entry["GapDirection"], f"edge {index + 1} GapDirection")
        normal = unit(entry["ProcessNormal"], f"edge {index + 1} ProcessNormal")
        if abs(np.dot(gap, normal)) > 1.0e-9:
            raise ValueError(f"edge {index + 1} has a nonorthogonal frame")
        if topology == "SpatialEdgeCluster":
            tangent = np.cross(gap, normal)
        else:
            tangent = unit(entry["Direction"], f"arm {index + 1} Direction")
            expected = np.cross(gap, normal)
            if abs(abs(np.dot(tangent, expected)) - 1.0) > 1.0e-9:
                raise ValueError(f"arm {index + 1} direction is inconsistent")
        interval = [float(value) for value in entry["Interval"]]
        if (
            len(interval) != 2
            or not all(math.isfinite(value) for value in interval)
            or interval[0] > 0.0
            or interval[1] < 0.0
            or interval[0] >= interval[1]
        ):
            raise ValueError(f"edge {index + 1} has an invalid Interval")
        local_point = frame @ point
        local_gap = frame @ gap
        local_normal = frame @ normal
        local_tangent = frame @ tangent
        if (
            abs(local_gap[2]) > 1.0e-8
            or abs(local_tangent[2]) > 1.0e-8
            or abs(abs(local_normal[2]) - 1.0) > 1.0e-8
        ):
            raise ValueError(
                "Automatic spatial coupons currently require parallel or antiparallel "
                "fabrication planes"
            )
        conductor = int(entry.get("Conductor", 1))
        slot = int(entry.get("InterfaceSlot", 0))
        if conductor <= 0 or slot < 0:
            raise ValueError("Conductor and InterfaceSlot labels must be nonnegative")
        boundary = entry.get(
            "BoundaryCondition", coupon.get("BoundaryCondition", {"Type": "PEC"})
        )
        if boundary.get("Type") != "PEC":
            raise ValueError(
                "Automatic spatial coupon calibration currently supports PEC only"
            )
        edges.append(
            {
                "Point": local_point.tolist(),
                "GapDirection": local_gap.tolist(),
                "ProcessNormal": local_normal.tolist(),
                "Tangent": local_tangent.tolist(),
                "Interval": interval,
                "Conductor": conductor,
                "InterfaceSlot": slot,
                "BoundaryCondition": boundary,
                "VertexArm": topology != "SpatialEdgeCluster",
            }
        )

    labels = []
    for edge in edges:
        conductor = edge["Conductor"]
        if conductor not in labels:
            labels.append(conductor)
        if conductor != labels.index(conductor) + 1:
            raise ValueError("Conductor labels are not canonical")
    if topology != "SpatialEdgeCluster" and labels != [1]:
        raise ValueError("Endpoint and junction coupons require one conductor")
    if any(np.linalg.norm(edge["Point"]) > 8.0 * radius for edge in edges):
        raise ValueError("Spatial coupon geometry is too large for its matching radius")

    facets = []
    for index, entry in enumerate(geometry.get("PlanViewFacets", [])):
        conductor = int(entry.get("Conductor", 0))
        points = np.asarray(entry.get("Points", []), dtype=float)
        if (
            conductor <= 0
            or points.ndim != 2
            or points.shape[0] < 3
            or points.shape[1] != 3
            or not np.all(np.isfinite(points))
        ):
            raise ValueError(f"plan-view facet {index + 1} is invalid")
        local = (frame @ points.T).T
        if np.ptp(local[:, 2]) > 1.0e-8 * radius:
            raise ValueError(
                f"plan-view facet {index + 1} is not on a process plane"
            )
        facets.append(
            {
                "Conductor": conductor,
                "Plane": float(np.mean(local[:, 2])),
                "Points": local[:, :2].tolist(),
            }
        )
    expected = {edge["Conductor"] for edge in edges}
    found = {facet["Conductor"] for facet in facets}
    if facets and found != expected:
        raise ValueError(
            "Plan-view facets must cover every spatial-coupon conductor"
        )
    return frame, edges, facets


def extended_interval(edge, radius):
    begin, end = edge["Interval"]
    extension = 2.0 * radius
    if edge["VertexArm"]:
        if abs(begin) <= 1.0e-10 * radius:
            end += extension
        elif abs(end) <= 1.0e-10 * radius:
            begin -= extension
    else:
        tolerance = 1.0e-10 * radius
        if begin <= -radius + tolerance:
            begin -= extension
        if end >= radius - tolerance:
            end += extension
    return begin, end


def metal_polygon(edge, radius):
    point = np.asarray(edge["Point"])[:2]
    tangent = np.asarray(edge["Tangent"])[:2]
    gap = np.asarray(edge["GapDirection"])[:2]
    begin, end = extended_interval(edge, radius)
    return np.asarray(
        (
            point + begin * tangent,
            point + end * tangent,
            point + end * tangent - 3.0 * radius * gap,
            point + begin * tangent - 3.0 * radius * gap,
        )
    )


def convex_polygons_overlap(first, second, tolerance):
    for polygon in (first, second):
        for start, end in zip(polygon, np.roll(polygon, -1, axis=0)):
            delta = end - start
            axis = np.asarray((-delta[1], delta[0]))
            norm = np.linalg.norm(axis)
            if norm <= tolerance:
                continue
            axis /= norm
            first_projection = first @ axis
            second_projection = second @ axis
            if (
                np.max(first_projection) <= np.min(second_projection) + tolerance
                or np.max(second_projection)
                <= np.min(first_projection) + tolerance
            ):
                return False
    return True


def validate_plan_view_geometry(edges, radius, facets):
    if facets:
        return
    tolerance = 1.0e-10 * radius
    polygons = [metal_polygon(edge, radius) for edge in edges]
    for first in range(len(edges)):
        for second in range(first + 1, len(edges)):
            a = edges[first]
            b = edges[second]
            same_layer = (
                abs(a["Point"][2] - b["Point"][2]) <= tolerance
                and np.dot(a["ProcessNormal"], b["ProcessNormal"])
                > 1.0 - 1.0e-10
            )
            if (
                same_layer
                and a["Conductor"] != b["Conductor"]
                and convex_polygons_overlap(
                    polygons[first], polygons[second], tolerance
                )
            ):
                raise ValueError(
                    "The edge-only spatial signature reconstructs overlapping "
                    f"plan-view metal for conductors {a['Conductor']} and "
                    f"{b['Conductor']}. Exact coupon generation requires additional "
                    "plan-view conductor boundaries."
                )


def coupon_bounds(edges, radius, metal_thickness, overetch):
    points = []
    for edge in edges:
        point = np.asarray(edge["Point"])
        tangent = np.asarray(edge["Tangent"])
        gap = np.asarray(edge["GapDirection"])
        begin, end = extended_interval(edge, radius)
        for coordinate in (begin, end):
            boundary = point + coordinate * tangent
            for side in (-1.0, 1.0):
                # coupon_bounds in mesh_spatial_coupon.jl adds one more radius
                # of padding, placing the transverse matching boundary 2R away.
                points.append(boundary + side * radius * gap)
    points = np.asarray(points)
    padding = radius
    lower = np.min(points, axis=0) - padding
    upper = np.max(points, axis=0) + padding
    lower[2] = min(
        lower[2],
        min(edge["Point"][2] for edge in edges) - radius - overetch,
    )
    upper[2] = max(
        upper[2],
        max(edge["Point"][2] for edge in edges) + radius + metal_thickness,
    )
    return lower, upper


def rectangle_perimeter_point(bounds, z, coordinate):
    xmin, ymin, _ = bounds[0]
    xmax, ymax, _ = bounds[1]
    width = xmax - xmin
    height = ymax - ymin
    perimeter = 2.0 * (width + height)
    coordinate %= perimeter
    if coordinate < width:
        return (xmin + coordinate, ymin, z)
    if coordinate < width + height:
        return (xmax, ymin + coordinate - width, z)
    if coordinate < 2.0 * width + height:
        return (xmax - coordinate + width + height, ymax, z)
    return (xmin, ymax - coordinate + 2.0 * width + height, z)


def rectangle_perimeter_coordinate(bounds, point, tolerance):
    xmin, ymin, _ = bounds[0]
    xmax, ymax, _ = bounds[1]
    width = xmax - xmin
    height = ymax - ymin
    x, y = point
    if abs(y - ymin) <= tolerance:
        return x - xmin
    if abs(x - xmax) <= tolerance:
        return width + y - ymin
    if abs(y - ymax) <= tolerance:
        return 2.0 * width + height - (x - xmin)
    if abs(x - xmin) <= tolerance:
        return 2.0 * (width + height) - (y - ymin)
    raise ValueError("Matching feature does not lie on the coupon perimeter")


def matching_perimeter_coordinates(
    bounds,
    levels,
    ring_size,
    edges,
    radius,
    metal_thickness,
    sidewall_angle,
    facets,
):
    if ring_size < 8 or ring_size % 4:
        raise ValueError("ring size must be a multiple of four and at least eight")
    xmin, ymin, _ = bounds[0]
    xmax, ymax, _ = bounds[1]
    width = xmax - xmin
    height = ymax - ymin
    perimeter = 2.0 * (width + height)
    coordinates = list(
        np.arange(ring_size, dtype=float) * perimeter / ring_size
    )
    tolerance = 1.0e-10 * max(width, height)
    pullback = metal_thickness / math.tan(math.radians(sidewall_angle))
    for facet in facets:
        for point in facet["Points"]:
            if (
                abs(point[0] - xmin) <= tolerance
                or abs(point[0] - xmax) <= tolerance
                or abs(point[1] - ymin) <= tolerance
                or abs(point[1] - ymax) <= tolerance
            ):
                coordinates.append(
                    rectangle_perimeter_coordinate(bounds, point, tolerance)
                )
    for level in levels:
        for edge in edges:
            sign = 1.0 if edge["ProcessNormal"][2] > 0.0 else -1.0
            height_from_plane = sign * (level - edge["Point"][2])
            if (
                height_from_plane < -tolerance
                or height_from_plane > metal_thickness + tolerance
            ):
                continue
            shift = pullback * np.clip(
                height_from_plane / metal_thickness, 0.0, 1.0
            )
            point = np.asarray(edge["Point"])[:2]
            tangent = np.asarray(edge["Tangent"])[:2]
            gap = np.asarray(edge["GapDirection"])[:2]
            begin, end = extended_interval(edge, radius)
            polygon = [
                point + begin * tangent - 3.0 * radius * gap,
                point + end * tangent - 3.0 * radius * gap,
                point + end * tangent - shift * gap,
                point + begin * tangent - shift * gap,
            ]
            for first, second in zip(polygon, polygon[1:] + polygon[:1]):
                delta = second - first
                intersections = []
                if abs(delta[0]) > tolerance:
                    for x in (xmin, xmax):
                        fraction = (x - first[0]) / delta[0]
                        if -tolerance <= fraction <= 1.0 + tolerance:
                            y = first[1] + fraction * delta[1]
                            if ymin - tolerance <= y <= ymax + tolerance:
                                intersections.append((x, np.clip(y, ymin, ymax)))
                if abs(delta[1]) > tolerance:
                    for y in (ymin, ymax):
                        fraction = (y - first[1]) / delta[1]
                        if -tolerance <= fraction <= 1.0 + tolerance:
                            x = first[0] + fraction * delta[0]
                            if xmin - tolerance <= x <= xmax + tolerance:
                                intersections.append((np.clip(x, xmin, xmax), y))
                coordinates.extend(
                    rectangle_perimeter_coordinate(bounds, value, tolerance)
                    for value in intersections
                )
    coordinates = sorted(value % perimeter for value in coordinates)
    unique = []
    for coordinate in coordinates:
        if not unique or coordinate - unique[-1] > tolerance:
            unique.append(coordinate)
    if len(unique) > 1 and perimeter - unique[-1] + unique[0] <= tolerance:
        unique.pop()
    return np.asarray(unique)


def rectangle_ring(bounds, z, coordinates):
    return np.asarray(
        [
            rectangle_perimeter_point(bounds, z, coordinate)
            for coordinate in coordinates
        ]
    )


def connect_rings(triangles, first, second, size):
    for index in range(size):
        following = (index + 1) % size
        triangles.append((first + index, first + following, second + following))
        triangles.append((first + index, second + following, second + index))


def cap_ring(triangles, points, offset, size, reverse):
    remaining = list(range(offset, offset + size))
    scale = np.ptp(points[remaining, :2], axis=0)
    area_tolerance = 1.0e-14 * max(scale[0] * scale[1], 1.0)
    while len(remaining) > 2:
        for position, current in enumerate(remaining):
            previous = remaining[position - 1]
            following = remaining[(position + 1) % len(remaining)]
            first = points[current, :2] - points[previous, :2]
            second = points[following, :2] - points[current, :2]
            signed_area = first[0] * second[1] - first[1] * second[0]
            if signed_area <= area_tolerance:
                continue
            triangle = (previous, current, following)
            triangles.append(tuple(reversed(triangle)) if reverse else triangle)
            remaining.pop(position)
            break
        else:
            polygon = points[remaining, :2]
            twice_area = np.sum(
                polygon[:, 0] * np.roll(polygon[:, 1], -1)
                - polygon[:, 1] * np.roll(polygon[:, 0], -1)
            )
            if abs(twice_area) <= 2.0 * area_tolerance:
                break
            raise ValueError("Unable to triangulate matching-surface cap")


def build_matching_surface(
    bounds,
    levels,
    ring_size,
    edges,
    radius,
    metal_thickness,
    sidewall_angle,
    facets,
):
    levels = sorted(set(float(value) for value in levels))
    coordinates = matching_perimeter_coordinates(
        bounds,
        levels,
        ring_size,
        edges,
        radius,
        metal_thickness,
        sidewall_angle,
        facets,
    )
    ring_size = len(coordinates)
    rings = [rectangle_ring(bounds, level, coordinates) for level in levels]
    points = np.vstack(rings)
    triangles = []
    for ring in range(len(rings) - 1):
        connect_rings(
            triangles, ring * ring_size, (ring + 1) * ring_size, ring_size
        )
    cap_ring(triangles, points, 0, ring_size, True)
    cap_ring(
        triangles,
        points,
        (len(rings) - 1) * ring_size,
        ring_size,
        False,
    )
    return points, np.asarray(triangles, dtype=int), [ring_size] * len(rings)


def points_in_polygon(points, polygon, tolerance):
    polygon = np.asarray(polygon)
    inside = np.zeros(len(points), dtype=bool)
    boundary = np.zeros(len(points), dtype=bool)
    previous = polygon[-1]
    for current in polygon:
        edge = current - previous
        relative = points - previous
        cross = edge[0] * relative[:, 1] - edge[1] * relative[:, 0]
        projection = relative @ edge
        edge_norm_squared = edge @ edge
        boundary |= (
            (np.abs(cross) <= tolerance * max(np.linalg.norm(edge), 1.0))
            & (projection >= -tolerance)
            & (projection <= edge_norm_squared + tolerance)
        )
        crossing = (previous[1] > points[:, 1]) != (
            current[1] > points[:, 1]
        )
        denominator = current[1] - previous[1]
        if abs(denominator) > tolerance:
            intersection = previous[0] + (
                (points[:, 1] - previous[1]) * edge[0] / denominator
            )
            inside ^= crossing & (intersection > points[:, 0])
        previous = current
    return inside | boundary


def points_in_plan_view_mask(points, facets, conductor, plane, tolerance):
    inside = np.zeros(len(points), dtype=bool)
    for facet in facets:
        if (
            facet["Conductor"] == conductor
            and abs(facet["Plane"] - plane) <= tolerance
        ):
            inside |= points_in_polygon(points, facet["Points"], tolerance)
    return inside


def conductor_at_points(
    points,
    edges,
    radius,
    metal_thickness,
    sidewall_angle,
    facets,
):
    labels = np.zeros(len(points), dtype=int)
    pullback = metal_thickness / math.tan(math.radians(sidewall_angle))
    width = 3.0 * radius
    tolerance = 1.0e-10 * radius
    for edge in edges:
        point = np.asarray(edge["Point"])
        tangent = np.asarray(edge["Tangent"])
        gap = np.asarray(edge["GapDirection"])
        sign = 1.0 if edge["ProcessNormal"][2] > 0.0 else -1.0
        height = sign * (points[:, 2] - point[2])
        active_height = (height >= -tolerance) & (
            height <= metal_thickness + tolerance
        )
        shift = pullback * np.clip(height / metal_thickness, 0.0, 1.0)
        delta = points - point
        longitudinal = delta @ tangent
        transverse = delta @ gap
        begin, end = extended_interval(edge, radius)
        active = (
            active_height
            & (longitudinal >= begin - tolerance)
            & (longitudinal <= end + tolerance)
            & (transverse <= -shift + tolerance)
            & (transverse >= -width - tolerance)
        )
        if facets:
            active &= points_in_plan_view_mask(
                points[:, :2],
                facets,
                edge["Conductor"],
                edge["Point"][2],
                tolerance,
            )
        conflict = active & (labels != 0) & (labels != edge["Conductor"])
        if np.any(conflict):
            raise ValueError("Different coupon conductors overlap on the matching surface")
        labels[active] = edge["Conductor"]
    return labels


def write_surface_trace(path, points, triangles, values):
    rows = []
    for triangle_index, triangle in enumerate(triangles, start=1):
        for vertex in triangle:
            rows.append((*points[vertex], values[vertex], triangle_index))
    np.savetxt(
        path,
        np.asarray(rows),
        delimiter=",",
        header="x,y,z,V,triangle",
        comments="",
        fmt=("%.16e", "%.16e", "%.16e", "%.16e", "%d"),
    )


def write_basis(output, points, triangles, active):
    trace_directory = output / "traces"
    trace_directory.mkdir(parents=True, exist_ok=True)
    for stale in trace_directory.glob("basis-*.csv"):
        stale.unlink()
    paths = []
    for basis, vertex in enumerate(active, start=1):
        values = np.zeros(len(points))
        values[vertex] = 1.0
        path = trace_directory / f"basis-{basis:04d}.csv"
        write_surface_trace(path, points, triangles, values)
        paths.append(path)
    zero = output / "zero-trace.csv"
    write_surface_trace(zero, points, triangles, np.zeros(len(points)))
    return paths, zero


def open_paths(contour_groups, labels, active, conductor_count):
    active_lookup = {vertex: index + 1 for index, vertex in enumerate(active)}
    paths = []
    offset = 0
    for size in contour_groups:
        ring = np.arange(offset, offset + size)
        free = [vertex for vertex in ring if labels[vertex] == 0]
        if free:
            # Split cyclic free arcs at every constrained run. Fully free rings remain
            # one path and are anchored below when the conductor graph is canonicalized.
            if len(free) == size:
                paths.append([active_lookup[vertex] for vertex in free])
            else:
                start = next(
                    index
                    for index in range(size)
                    if labels[ring[index]] == 0
                    and labels[ring[index - 1]] != 0
                )
                run = []
                for step in range(size):
                    vertex = ring[(start + step) % size]
                    if labels[vertex] == 0:
                        run.append(active_lookup[vertex])
                    elif run:
                        paths.append(run)
                        run = []
                if run:
                    paths.append(run)
        offset += size
    if len(paths) < conductor_count - 1:
        raise ValueError("Trace mesh has too few open arcs to connect all conductors")
    result = []
    for index, indices in enumerate(paths):
        end = index + 2 if index < conductor_count - 1 else 2
        result.append(
            {
                "Indices": indices,
                "StartConductor": 1,
                "EndConductor": end,
            }
        )
    assigned = sorted(index for path in result for index in path["Indices"])
    if assigned != list(range(1, len(active) + 1)):
        raise ValueError("Open contour paths do not partition the active trace")
    return result


def canonical_points(local_points, frame):
    return local_points @ frame


def dielectric(
    index,
    attributes,
    interface_type,
    radius,
    thickness,
    permittivity,
    edge_frame_normal,
    edge_attributes,
):
    return {
        "Index": index,
        "Attributes": attributes,
        "Type": interface_type,
        "Thickness": thickness,
        "Permittivity": permittivity,
        "LossTan": INTERFACE_LOSS[interface_type],
        "LocalizeEdgeEnergy": True,
        "EdgeAttributes": edge_attributes,
        "EdgeExcludeAttributes": [1],
        "EdgeDistances": [radius],
        "EdgeFrameNormal": edge_frame_normal,
    }


def conductor_attributes(edges, fabricated, conductor):
    if fabricated:
        return [5000 + conductor, 6000 + conductor]
    return [4000 + conductor]


def interface_attributes(edges, fabricated, slot, interface_type):
    conductors = sorted(
        {edge["Conductor"] for edge in edges}
    )
    if interface_type == "SA":
        return [3000 + slot, 3100 + slot] if fabricated else [3000 + slot]
    if fabricated:
        base = 5000 if interface_type == "MS" else 6000
        return [base + conductor for conductor in conductors]
    return [4000 + conductor for conductor in conductors]


def edge_attributes(edges, fabricated, slot):
    return [3000 + slot, 3100 + slot] if fabricated else [3000 + slot]


def make_config(
    output,
    name,
    mesh,
    traces,
    zero_trace,
    terminal_conductors,
    fabricated,
    order,
    radius,
    substrate_permittivity,
    interface_layers,
    interfaces,
    edges,
    response_matrix=True,
):
    potentials = [
        {"Index": index, "Attributes": [1], "DataFile": str(trace)}
        for index, trace in enumerate(traces, start=1)
    ]
    for conductor in terminal_conductors:
        potentials.append(
            {
                "Index": len(potentials) + 1,
                "Attributes": [1],
                "TerminalAttributes": conductor_attributes(
                    edges, fabricated, conductor
                ),
                "DataFile": str(zero_trace),
            }
        )
    dielectric_data = []
    for index, interface in enumerate(interfaces, start=1):
        interface_type = interface["Type"]
        slot_normals = {
            tuple(edge["ProcessNormal"])
            for edge in edges
            if edge["InterfaceSlot"] == interface["Slot"]
        }
        if len(slot_normals) != 1:
            raise ValueError(
                f"Interface slot {interface['Slot']} has inconsistent process normals"
            )
        dielectric_data.append(
            dielectric(
                index,
                interface_attributes(
                    edges, fabricated, interface["Slot"], interface_type
                ),
                interface_type,
                radius,
                *interface_layers[interface_type],
                list(next(iter(slot_normals))),
                edge_attributes(edges, fabricated, interface["Slot"]),
            )
        )
    conductor_count = max(edge["Conductor"] for edge in edges)
    ground = [
        attribute
        for conductor in range(1, conductor_count + 1)
        for attribute in conductor_attributes(edges, fabricated, conductor)
    ]
    return {
        "Problem": {
            "Type": "Electrostatic",
            "Verbose": 1,
            "Output": str(output / "postpro" / name),
            "OutputFormats": {"Paraview": False, "GridFunction": False},
        },
        "Model": {
            "Mesh": str(mesh),
            "L0": 1.0e-6,
            "Refinement": {"MaxIts": 0},
        },
        "Domains": {
            "Materials": [
                {"Attributes": [1], "Permittivity": substrate_permittivity},
                {"Attributes": [2], "Permittivity": 1.0},
            ],
            "Postprocessing": {
                "Energy": [
                    {"Index": 1, "Attributes": [1]},
                    {"Index": 2, "Attributes": [2]},
                ]
            },
        },
        "Boundaries": {
            "Ground": {"Attributes": ground},
            "PrescribedPotential": potentials,
            "Postprocessing": {"Dielectric": dielectric_data},
        },
        "Solver": {
            "Order": order,
            "Electrostatic": {
                "Save": 0,
                "ResponseMatrix": response_matrix,
                "AggregateResponseMatrix": response_matrix,
            },
            "Linear": {
                "Type": "BoomerAMG",
                "KSPType": "CG",
                "Tol": 1.0e-10,
                "MaxIts": 1000,
                "EstimatorTol": 5.0e-1,
                "EstimatorMaxIts": 5,
                "EstimatorMG": True,
            },
        },
    }


def model_interfaces(coupon):
    interfaces = []
    seen = set()
    for entry in coupon.get("Interfaces", []):
        key = (int(entry.get("Slot", 0)), entry["Type"])
        if key in seen:
            continue
        seen.add(key)
        interfaces.append({"Slot": key[0], "Type": key[1]})
    interfaces.sort(key=lambda item: (item["Slot"], item["Type"]))
    if not interfaces:
        raise ValueError("Spatial coupon has no target interfaces")
    if any(item["Type"] not in INTERFACE_DEFAULTS for item in interfaces):
        raise ValueError("Spatial coupons require explicit SA, MS, or MA interfaces")
    return interfaces


def write_library(
    output,
    coupon,
    radius,
    basis_path,
    contour_groups,
    zero_indices,
    paths,
    references,
    interface_entries,
    fabrication,
    model_name,
):
    topology = coupon["Topology"]
    geometry = coupon["Geometry"]
    model = {
        "Name": model_name,
        "Topology": topology,
        "FabricatedMatrix": "postpro/spatial_fabricated/domain-response-matrix.csv",
        "ThinMatrix": "postpro/spatial_thin/domain-response-matrix.csv",
        "FabricatedSurfaceMatrix":
            "postpro/spatial_fabricated/surface-response-matrix.csv",
        "ThinSurfaceMatrix": "postpro/spatial_thin/surface-response-matrix.csv",
        "BasisPoints": str(basis_path.relative_to(output)),
        "Interfaces": [
            {
                "Slot": interface["Slot"],
                "Type": interface["Type"],
                "Coupon": index,
            }
            for index, interface in enumerate(interface_entries, start=1)
        ],
    }
    if topology == "Endpoint":
        model["Reference"] = references[0]
    elif topology == "Junction":
        model["Reference"] = references[0]
        model["ArmAngles"] = geometry["ArmAnglesDegrees"]
        model["ArmAngleTolerance"] = 1.0e-6
    elif topology == "SpatialEdgeCluster":
        model["Edges"] = []
        for source in geometry["Edges"]:
            edge = dict(source)
            edge["BoundaryCondition"] = "PEC"
            model["Edges"].append(edge)
        if "PlanViewBoundary" in geometry:
            model["PlanViewBoundary"] = geometry["PlanViewBoundary"]
        model["EdgePositionTolerance"] = 1.0e-6 * radius
        model["EdgeAngleTolerance"] = 1.0e-6
        model["ConductorReferences"] = references
    else:
        raise ValueError(f"Unsupported spatial topology {topology}")
    if paths:
        model["OpenContourPaths"] = paths
    else:
        model["ContourGroups"] = contour_groups
        if zero_indices:
            model["ZeroTraceIndices"] = zero_indices
        model["BoundaryCondition"] = "PEC"
    library = {
        "Version": 3,
        "Name": model_name,
        "MatchingRadius": radius,
        "Fabrication": fabrication,
        "Models": [model],
    }
    path = output / "process-library.json"
    path.write_text(json.dumps(library, indent=2) + "\n")
    return path


def reference_points(coupon, edges, frame, radius):
    conductor_count = max(edge["Conductor"] for edge in edges)
    references = []
    for conductor in range(1, conductor_count + 1):
        edge = next(item for item in edges if item["Conductor"] == conductor)
        local = np.asarray(edge["Point"]) - 0.5 * radius * np.asarray(
            edge["GapDirection"]
        )
        references.append((local @ frame).tolist())
    return references


def write_mesh_signature(path, edges):
    header = (
        "Index,Slot,Conductor,Px,Py,Pz,Gx,Gy,Gz,Tx,Ty,Tz,Nz,S0,S1,VertexArm"
    )
    lines = [header]
    for index, edge in enumerate(edges, start=1):
        values = [
            index,
            edge["InterfaceSlot"],
            edge["Conductor"],
            *edge["Point"],
            *edge["GapDirection"],
            *edge["Tangent"],
            edge["ProcessNormal"][2],
            *edge["Interval"],
            int(edge["VertexArm"]),
        ]
        lines.append(",".join(f"{value:.17g}" for value in values))
    path.write_text("\n".join(lines) + "\n")


def write_plan_view_mask(path, facets):
    lines = ["Facet,Conductor,Plane,X,Y"]
    for facet_index, facet in enumerate(facets, start=1):
        for point in facet["Points"]:
            lines.append(
                ",".join(
                    f"{value:.17g}"
                    for value in (
                        facet_index,
                        facet["Conductor"],
                        facet["Plane"],
                        *point,
                    )
                )
            )
    path.write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("coupon", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--thin-mesh", type=Path)
    parser.add_argument("--fabricated-mesh", type=Path)
    parser.add_argument("--radius", type=float, required=True)
    parser.add_argument("--metal-thickness", type=float, required=True)
    parser.add_argument("--overetch-depth", type=float, required=True)
    parser.add_argument("--sidewall-angle", type=float, required=True)
    parser.add_argument("--top-rounding", type=float, required=True)
    parser.add_argument("--trench-rounding", type=float, required=True)
    parser.add_argument("--substrate-permittivity", type=float, default=11.47)
    parser.add_argument("--ring-size", type=int, default=32)
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--model-name", required=True)
    parser.add_argument(
        "--signature-only",
        action="store_true",
        help="Validate the coupon and write only mesh-signature.csv",
    )
    for interface, (thickness, permittivity) in INTERFACE_DEFAULTS.items():
        parser.add_argument(
            f"--{interface.lower()}-thickness", type=float, default=thickness
        )
        parser.add_argument(
            f"--{interface.lower()}-permittivity",
            type=float,
            default=permittivity,
        )
    args = parser.parse_args()
    if args.radius <= 0.0 or args.metal_thickness <= 0.0:
        parser.error("radius and metal thickness must be positive")
    if not 0.0 <= args.overetch_depth < args.radius:
        parser.error("overetch depth must lie in [0, radius)")
    if not 0.0 < args.sidewall_angle <= 90.0:
        parser.error("sidewall angle must lie in (0, 90]")
    if args.order <= 0:
        parser.error("order must be positive")

    coupon_path = args.coupon.expanduser().resolve()
    coupon = load_json(coupon_path)
    output = args.output.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    failure_path = output / "generation-failure.json"
    try:
        frame, edges, facets = normalize_geometry(coupon, args.radius)
        validate_plan_view_geometry(edges, args.radius, facets)
        interfaces = model_interfaces(coupon)
    except ValueError as error:
        failure_path.write_text(
            json.dumps(
                {
                    "Version": 1,
                    "Stage": "SpatialGeometry",
                    "Reason": str(error),
                },
                indent=2,
            )
            + "\n"
        )
        parser.error(str(error))
    failure_path.unlink(missing_ok=True)
    write_mesh_signature(output / "mesh-signature.csv", edges)
    mask_path = output / "plan-view-mask.csv"
    if facets:
        write_plan_view_mask(mask_path, facets)
    else:
        mask_path.unlink(missing_ok=True)
    if args.signature_only:
        print(output / "mesh-signature.csv")
        return
    if args.thin_mesh is None or args.fabricated_mesh is None:
        parser.error("--thin-mesh and --fabricated-mesh are required")

    lower, upper = coupon_bounds(
        edges, args.radius, args.metal_thickness, args.overetch_depth
    )
    bounds = np.vstack((lower, upper))
    levels = [lower[2], upper[2]]
    for edge in edges:
        plane = edge["Point"][2]
        sign = 1.0 if edge["ProcessNormal"][2] > 0.0 else -1.0
        levels.extend(
            (
                plane,
                plane + sign * args.metal_thickness,
                plane - sign * args.overetch_depth,
            )
        )
    points, triangles, contour_groups = build_matching_surface(
        bounds,
        levels,
        args.ring_size,
        edges,
        args.radius,
        args.metal_thickness,
        args.sidewall_angle,
        facets,
    )
    labels = conductor_at_points(
        points,
        edges,
        args.radius,
        args.metal_thickness,
        args.sidewall_angle,
        facets,
    )
    conductor_count = max(edge["Conductor"] for edge in edges)
    if conductor_count == 1:
        active = np.arange(len(points))
        zero_indices = (np.flatnonzero(labels) + 1).tolist()
        paths = []
    else:
        active = np.flatnonzero(labels == 0)
        zero_indices = []
        paths = open_paths(
            contour_groups, labels, active.tolist(), conductor_count
        )
    traces, zero_trace = write_basis(
        output, points, triangles, active.tolist()
    )
    canonical = canonical_points(points[active], frame)
    basis_path = output / "basis-points.csv"
    np.savetxt(
        basis_path,
        canonical,
        delimiter=",",
        header="x,y,z",
        comments="",
        fmt="%.16e",
    )

    scale = np.maximum(upper - lower, np.finfo(float).tiny)
    normalized = (points - 0.5 * (lower + upper)) / scale
    heldout_values = (
        0.35
        + 0.20 * normalized[:, 0]
        - 0.15 * normalized[:, 1]
        + 0.10 * normalized[:, 2]
        + 0.08 * normalized[:, 0] * normalized[:, 1]
        + 0.06 * normalized[:, 2] ** 2
    )
    heldout_values[labels != 0] = 0.0
    heldout_trace = output / "heldout-trace.csv"
    write_surface_trace(heldout_trace, points, triangles, heldout_values)
    heldout_coefficients = np.concatenate(
        (heldout_values[active], np.ones(conductor_count - 1))
    )
    np.savetxt(
        output / "heldout_coefficients.csv",
        heldout_coefficients,
        delimiter=",",
        header="coefficient_V",
        comments="",
        fmt="%.16e",
    )

    probes = (
        ("common", "cutoff", np.ones(len(points))),
        ("x-linear", "cutoff*x/Lx", normalized[:, 0]),
        ("y-linear", "cutoff*y/Ly", normalized[:, 1]),
        ("z-linear", "cutoff*z/Lz", normalized[:, 2]),
        ("xy-mixed", "cutoff*x*y/(Lx*Ly)", normalized[:, 0] * normalized[:, 1]),
        ("z-quadratic", "cutoff*z^2/Lz^2", normalized[:, 2] ** 2),
    )
    probe_paths = []
    probe_metadata = []
    for index, (name, expression, values) in enumerate(probes, start=1):
        values = values.copy()
        values[labels != 0] = 0.0
        path = output / f"probe-{index:02d}.csv"
        write_surface_trace(path, points, triangles, values)
        probe_paths.append(path)
        probe_metadata.append({"Name": name, "Expression": expression})
    for conductor in range(2, conductor_count + 1):
        probe_metadata.append(
            {
                "Name": f"conductor-{conductor}",
                "Expression": f"V_{conductor}-V_1",
            }
        )
    (output / "probe-manifest.json").write_text(
        json.dumps({"Version": 1, "Probes": probe_metadata}, indent=2) + "\n"
    )

    interface_layers = {
        interface: (
            getattr(args, f"{interface.lower()}_thickness"),
            getattr(args, f"{interface.lower()}_permittivity"),
        )
        for interface in INTERFACE_DEFAULTS
    }
    terminal_conductors = list(range(2, conductor_count + 1))
    meshes = {
        "thin": args.thin_mesh.expanduser().resolve(),
        "fabricated": args.fabricated_mesh.expanduser().resolve(),
    }
    for mesh in meshes.values():
        if not mesh.is_file():
            raise FileNotFoundError(mesh)
    for kind, mesh in meshes.items():
        fabricated = kind == "fabricated"
        name = f"spatial_{kind}"
        config = make_config(
            output,
            name,
            mesh,
            traces,
            zero_trace,
            terminal_conductors,
            fabricated,
            args.order,
            args.radius,
            args.substrate_permittivity,
            interface_layers,
            interfaces,
            edges,
        )
        (output / f"{name}.json").write_text(json.dumps(config, indent=2) + "\n")
        heldout_name = f"heldout_{name}"
        heldout = make_config(
            output,
            heldout_name,
            mesh,
            [heldout_trace],
            zero_trace,
            [],
            fabricated,
            args.order,
            args.radius,
            args.substrate_permittivity,
            interface_layers,
            interfaces,
            edges,
            False,
        )
        heldout_terminals = [
            attribute
            for conductor in terminal_conductors
            for attribute in conductor_attributes(edges, fabricated, conductor)
        ]
        if heldout_terminals:
            heldout["Boundaries"]["PrescribedPotential"][0][
                "TerminalAttributes"
            ] = heldout_terminals
        (output / f"{heldout_name}.json").write_text(
            json.dumps(heldout, indent=2) + "\n"
        )
        probe_name = f"probe-{kind}"
        probe = make_config(
            output,
            probe_name,
            mesh,
            probe_paths,
            zero_trace,
            terminal_conductors,
            fabricated,
            args.order,
            args.radius,
            args.substrate_permittivity,
            interface_layers,
            interfaces,
            edges,
        )
        (output / f"{probe_name}.json").write_text(
            json.dumps(probe, indent=2) + "\n"
        )

    fabrication = {
        "LengthUnit": "um",
        "MetalThickness": args.metal_thickness,
        "OveretchDepth": args.overetch_depth,
        "SidewallAngleDegrees": args.sidewall_angle,
        "TopRoundingRadius": args.top_rounding,
        "BottomRoundingRadius": args.trench_rounding,
        "SubstratePermittivity": args.substrate_permittivity,
        "InterfaceLayers": {
            interface: {
                "Thickness": values[0],
                "Permittivity": values[1],
            }
            for interface, values in interface_layers.items()
        },
    }
    library = write_library(
        output,
        coupon,
        args.radius,
        basis_path,
        contour_groups,
        zero_indices,
        paths,
        reference_points(coupon, edges, frame, args.radius),
        interfaces,
        fabrication,
        args.model_name,
    )
    print(output / "mesh-signature.csv")
    print(output / "spatial_thin.json")
    print(output / "spatial_fabricated.json")
    print(output / "heldout_spatial_thin.json")
    print(output / "heldout_spatial_fabricated.json")
    print(output / "probe-thin.json")
    print(output / "probe-fabricated.json")
    print(library)


if __name__ == "__main__":
    main()
