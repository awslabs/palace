#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The corner coupon's trace RESOLVABILITY gate (block (b) family 4 round 2, supervisor
decision 328 (1)): the order-p boundary node lattice is MFEM's (closed Gauss-Lobatto), the
MSH 2.2 reader finds the matching-surface triangles, the counter reproduces the failure of the
concave 48.75-degree coupon (the z = +-R cap rings' inner hats on a 300-nm far mesh hold too
few — or no — nodes: FAIL before) and passes on a mesh graded at the knots (PASS after)."""

import importlib.util
import json
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


GENERATOR = load("generate_corner_response")
GATE = load("trace_resolvability")
MESHER = ROOT / "mesh_corner_coupon.jl"
RADIUS, THICKNESS, OVERETCH = 1.9, 0.1, 0.05
FAR_SIZE = 0.3  # the recipe's lc_far: the mesh size at the z = +-R faces


def julia_with_gmsh():
    """The julia executable when Gmsh.jl is importable (the mesher test), else None."""
    julia = shutil.which("julia")
    if julia is None:
        return None
    probe = subprocess.run(
        [julia, "-e", "import Gmsh"], capture_output=True, text=True, timeout=240
    )
    return julia if probe.returncode == 0 else None


def write_trace(directory, angle, topology="concave"):
    surface = GENERATOR.build_surface(
        RADIUS, 16, THICKNESS, OVERETCH, angle_degrees=angle, topology=topology,
        rule=GENERATOR.REFINED_RULE,
    )
    GENERATOR.write_trace_mesh(directory, surface)
    (directory / "process-library.json").write_text(
        json.dumps(
            {"Models": [{"ZeroTraceIndices": [i + 1 for i in surface.zero_trace_indices()]}]}
        )
    )
    return surface


def face_grid_msh(path, breakpoints_of_face):
    """An ASCII MSH 2.2 file whose physical surface 1 is the matching box |x|, |y|, |z| <= R
    triangulated face by face on tensor grids (`breakpoints_of_face(axis, sign)` gives the
    two in-face coordinate arrays), first-order triangles."""
    nodes = []
    index_of = {}

    def node(point):
        key = tuple(np.round(point, 12))
        if key not in index_of:
            index_of[key] = len(nodes) + 1
            nodes.append(key)
        return index_of[key]

    triangles = []
    for axis in range(3):
        for sign in (-1.0, 1.0):
            u_axis, v_axis = [a for a in range(3) if a != axis]
            u, v = breakpoints_of_face(axis, sign)
            for i in range(len(u) - 1):
                for j in range(len(v) - 1):
                    corners = []
                    for du, dv in ((0, 0), (1, 0), (1, 1), (0, 1)):
                        point = [0.0, 0.0, 0.0]
                        point[axis] = sign * RADIUS
                        point[u_axis] = u[i + du]
                        point[v_axis] = v[j + dv]
                        corners.append(node(point))
                    triangles.append((corners[0], corners[1], corners[2]))
                    triangles.append((corners[0], corners[2], corners[3]))
    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat", "$Nodes", str(len(nodes))]
    lines += [f"{i + 1} {x:.17g} {y:.17g} {z:.17g}" for i, (x, y, z) in enumerate(nodes)]
    lines += ["$EndNodes", "$Elements", str(len(triangles) + 1)]
    # A tetrahedron of another physical tag: the reader must skip it.
    lines.append("1 4 2 1 1 1 2 3 4")
    lines += [
        f"{i + 2} 2 2 1 7 {a} {b} {c}" for i, (a, b, c) in enumerate(triangles)
    ]
    lines += ["$EndElements", ""]
    path.write_text("\n".join(lines))
    return len(triangles)


def uniform_breakpoints(h):
    count = int(np.ceil(2.0 * RADIUS / h))
    return np.linspace(-RADIUS, RADIUS, count + 1)


class LatticeTest(unittest.TestCase):
    def test_closed_gauss_lobatto_points_are_mfems(self):
        np.testing.assert_allclose(GATE.gauss_lobatto_closed_points(1), [0.0, 1.0])
        np.testing.assert_allclose(GATE.gauss_lobatto_closed_points(2), [0.0, 0.5, 1.0])
        np.testing.assert_allclose(
            GATE.gauss_lobatto_closed_points(3),
            [0.0, 0.5 * (1.0 - 1.0 / np.sqrt(5.0)), 0.5 * (1.0 + 1.0 / np.sqrt(5.0)), 1.0],
        )
        np.testing.assert_allclose(
            GATE.gauss_lobatto_closed_points(4),
            [0.0, 0.5 * (1.0 - np.sqrt(3.0 / 7.0)), 0.5, 0.5 * (1.0 + np.sqrt(3.0 / 7.0)), 1.0],
        )

    def test_triangle_lattice_count_and_symmetry(self):
        for p in (1, 2, 3, 4, 5):
            lattice = GATE.lagrange_triangle_lattice(p)
            self.assertEqual(len(lattice), (p + 1) * (p + 2) // 2)
            # Inside the reference triangle, vertices present, symmetric under the vertex
            # permutation (x, y) -> (y, 1 - x - y).
            self.assertTrue(np.all(lattice >= -1e-15))
            self.assertTrue(np.all(lattice.sum(axis=1) <= 1.0 + 1e-15))
            for vertex in ((0.0, 0.0), (1.0, 0.0), (0.0, 1.0)):
                self.assertTrue(np.any(np.all(np.isclose(lattice, vertex), axis=1)))
            permuted = np.column_stack([lattice[:, 1], 1.0 - lattice.sum(axis=1)])
            for point in permuted:
                self.assertTrue(np.any(np.all(np.isclose(lattice, point, atol=1e-12), axis=1)))

    def test_required_count_scales_with_order_squared(self):
        self.assertEqual(GATE.required_active_nodes(3), 54)
        self.assertEqual(GATE.required_active_nodes(4), 96)


class GateTest(unittest.TestCase):
    def test_fails_on_the_far_mesh_and_passes_graded_at_the_knots(self):
        # The concave 48.75-degree node (decision 325): on a uniform 300-nm surface mesh at
        # p4 the cap rings' inner free hats (7.4 nm apart) hold a few nodes or none (FAIL,
        # the wide hats hold hundreds); on a surface mesh graded like the mesher's every free
        # hat holds >= 6 p^2 nodes (PASS); the zero set is never counted against.
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            surface = write_trace(directory, 48.75)
            zero = GATE.read_zero_trace_indices(directory)
            self.assertEqual(zero, [i + 1 for i in surface.zero_trace_indices()])
            far = directory / "far.msh"
            face_grid_msh(far, lambda axis, sign: (uniform_breakpoints(FAR_SIZE),) * 2)
            audit = GATE.audit_mesh(far, directory, 4, RADIUS, zero)
            self.assertFalse(audit["Passed"])
            counts = np.asarray(audit["ActiveNodes"])
            # The five inner free hats of the two cap rings (1-based 148-152, 164-168).
            cap = np.concatenate([counts[147:152], counts[163:168]])
            self.assertLess(cap.max(), 20)
            self.assertTrue(set(range(148, 153)) <= set(audit["FailingFreeHats"]))
            self.assertTrue(set(range(164, 169)) <= set(audit["FailingFreeHats"]))
            # A wide hat (the graded R/3 knot of the z = -R ring, 1-based 2) is resolved.
            self.assertGreater(counts[1], 96)
            # The zero set is not gated.
            for index in zero:
                self.assertNotIn(index, audit["FailingFreeHats"])

            # PASS after: the box triangulated at 50 nm (the cap faces and the side faces'
            # far parts) and 10 nm across the process band z in [-3 d, t + 6 d] of the side
            # faces — a surface mesh that samples every hat of the 48.75-degree layout as the
            # mesher's knot-gap sizing does (the real coupon reads 270-290 on the cap hats).
            def graded(axis, sign):
                if axis == 2:
                    return (uniform_breakpoints(0.05),) * 2
                band = np.arange(-3.0 * OVERETCH, THICKNESS + 6.0 * OVERETCH + 1e-12, 0.01)
                return uniform_breakpoints(0.05), np.union1d(uniform_breakpoints(0.05), band)

            fine = directory / "graded.msh"
            face_grid_msh(fine, graded)
            audit = GATE.audit_mesh(fine, directory, 4, RADIUS, zero)
            self.assertTrue(audit["Passed"], audit["FailingFreeHats"])
            self.assertGreaterEqual(audit["MinimumActiveNodesFreeHat"], 96)
            counts = np.asarray(audit["ActiveNodes"])
            self.assertGreater(min(counts[147:152].min(), counts[163:168].min()), 150)

    def test_reader_skips_other_tags_and_requires_the_surface(self):
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            path = directory / "far.msh"
            count = face_grid_msh(path, lambda axis, sign: (uniform_breakpoints(1.0),) * 2)
            coordinates, triangles = GATE.read_msh2_boundary_triangles(path, 1)
            self.assertEqual(len(triangles), count)
            self.assertEqual(set(len(t) for t in triangles), {3})
            with self.assertRaises(ValueError):
                GATE.read_msh2_boundary_triangles(path, 5)
            nodes = GATE.boundary_lattice_nodes(coordinates, triangles, 2, RADIUS)
            # p2 on a 4 x 4 grid per face: the box's (2 n + 1)^2 lattice per face, shared
            # edges and corners counted once: 6 x 81 - 12 x 9 + 8.
            self.assertEqual(len(nodes), 6 * 81 - 12 * 9 + 8)


class MesherTest(unittest.TestCase):
    def test_mesher_knot_gap_sizing_passes_the_gate(self):
        # End to end at a cheap resolution (lc_fine 0.1 / lc_far 0.6, ~5 s, < 1 GB): the
        # concave 48.75-degree THIN coupon meshed WITHOUT the trace mesh fails the gate at p4
        # (cap hats with a node or none), meshed WITH `--trace-mesh` (the knot-gap size
        # fields) every free hat holds >= 6 p^2 nodes.
        julia = julia_with_gmsh()
        if julia is None:
            self.skipTest("julia with Gmsh.jl is not available")
        with tempfile.TemporaryDirectory() as directory:
            directory = Path(directory)
            write_trace(directory, 48.75)
            zero = GATE.read_zero_trace_indices(directory)
            results = {}
            for name, extra in (("base", []), ("knots", ["--trace-mesh", str(directory)])):
                mesh = directory / f"{name}.msh"
                subprocess.run(
                    [
                        julia, str(MESHER), "concave-thin", str(mesh), "--radius", str(RADIUS),
                        "--angle", "48.75", "--corner-radius", "0", "--metal-thickness",
                        str(THICKNESS), "--overetch", str(OVERETCH), "--sidewall-angle", "90",
                        "--top-radius", "0", "--bottom-radius", "0", "--lc-fine", "0.1",
                        "--lc-far", "0.6", "--mesh-order", "2", *extra,
                    ],
                    check=True, capture_output=True, text=True, timeout=240,
                )
                results[name] = GATE.audit_mesh(mesh, directory, 4, RADIUS, zero)
            self.assertFalse(results["base"]["Passed"])
            self.assertLess(results["base"]["MinimumActiveNodesFreeHat"], 10)
            self.assertTrue(results["knots"]["Passed"], results["knots"]["FailingFreeHats"])
            self.assertGreaterEqual(results["knots"]["MinimumActiveNodesFreeHat"], 96)


if __name__ == "__main__":
    unittest.main()
