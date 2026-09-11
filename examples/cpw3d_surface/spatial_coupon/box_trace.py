# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Complete and validate box-boundary traces without dropping retained nodes."""
import collections
import numpy as np
from generate_spatial_response import cap_ring, connect_rings, rectangle_perimeter_coordinate


def complete_box_trace(retained_points):
    original = np.asarray(retained_points, dtype=float)
    if original.ndim != 2 or original.shape[1] != 3 or not np.isfinite(original).all():
        raise ValueError("Invalid retained trace points")
    lower, upper = original.min(axis=0), original.max(axis=0)
    bounds = np.asarray([lower, upper])
    tolerance = 1e-10 * max(upper - lower)
    if np.any(upper <= lower):
        raise ValueError("Degenerate matching box")
    xy = sorted({tuple(point[:2]) for point in original})
    snapped = {}
    for corner in ((lower[0], lower[1]), (upper[0], lower[1]),
                   (upper[0], upper[1]), (lower[0], upper[1])):
        nearby = [point for point in xy if max(abs(a-b) for a, b in zip(point, corner)) <= tolerance]
        if len(nearby) > 1:
            raise ValueError("Distinct retained basis nodes collapse at a box corner")
        if nearby:
            old = nearby[0]
            xy[xy.index(old)] = corner
            snapped[old] = corner
        else:
            xy.append(corner)
    xy = sorted(xy, key=lambda point: rectangle_perimeter_coordinate(bounds, point, tolerance))
    levels = sorted(set(original[:, 2]))
    points = np.asarray([(x, y, z) for z in levels for x, y in xy])
    lookup = {tuple(point): i for i, point in enumerate(points)}
    old_to_vertex = [lookup[(*snapped.get(tuple(point[:2]), tuple(point[:2])), point[2])]
                     for point in original]
    if len(set(old_to_vertex)) != len(original):
        raise ValueError("Retained trace nodes were collapsed")
    triangles = []
    size = len(xy)
    for ring in range(len(levels)-1):
        connect_rings(triangles, ring*size, (ring+1)*size, size)
    cap_ring(triangles, points, 0, size, True)
    cap_ring(triangles, points, (len(levels)-1)*size, size, False)
    triangles = np.asarray(triangles, dtype=int)
    report = validate_box_trace(points, triangles, bounds)
    report.update(OriginalVertices=len(original), AddedVertices=len(points)-len(original),
                  MaximumRetainedNodeDisplacement=float(np.max(np.abs(points[old_to_vertex]-original))),
                  RetainedVertexMap=old_to_vertex)
    return points, triangles, [size]*len(levels), report


def validate_box_trace(points, triangles, bounds=None):
    points, triangles = np.asarray(points, dtype=float), np.asarray(triangles, dtype=int)
    bounds = np.asarray(bounds if bounds is not None else [points.min(axis=0), points.max(axis=0)])
    lower, upper = bounds
    tolerance = 1e-9 * max(upper-lower)
    if len({tuple(point) for point in points}) != len(points):
        raise ValueError("Coincident distinct trace nodes")
    areas = np.zeros((3,2))
    edges = collections.defaultdict(list)
    used = set()
    for triangle in triangles:
        if len(set(triangle)) != 3 or min(triangle) < 0 or max(triangle) >= len(points):
            raise ValueError("Invalid trace triangle connectivity")
        xyz = points[triangle]
        area = np.linalg.norm(np.cross(xyz[1]-xyz[0], xyz[2]-xyz[0]))/2
        if not area > 0:
            raise ValueError("Zero-area trace triangle")
        faces = [(d,side) for d in range(3) for side in (0,1)
                 if np.all(np.abs(xyz[:,d]-bounds[side,d]) <= tolerance)]
        if len(faces) != 1:
            raise ValueError(f"Trace triangle does not lie on one matching-box face: {xyz.tolist()}")
        areas[faces[0]] += area
        used.update(triangle)
        for i in range(3):
            a,b = int(triangle[i]),int(triangle[(i+1)%3])
            edges[min(a,b),max(a,b)].append(1 if a < b else -1)
    if used != set(range(len(points))):
        raise ValueError("Unused trace node (possibly a dropped cap boundary node)")
    if any(len(signs) != 2 or sum(signs) != 0 for signs in edges.values()):
        raise ValueError("Trace surface is open, nonconforming, or inconsistently oriented")
    sizes = upper-lower
    expected = np.asarray([[np.prod(np.delete(sizes,d))]*2 for d in range(3)])
    worst = float(np.max(np.abs(areas/expected-1)))
    if worst > 1e-9:
        raise ValueError(f"Matching-box face coverage differs by {worst}")
    return {"Vertices":len(points), "Triangles":len(triangles), "OffBoxTriangles":0,
            "ClosedOrientedSurface":True, "MaximumRelativeFaceAreaDifference":worst,
            "Lower":lower.tolist(), "Upper":upper.tolist()}
