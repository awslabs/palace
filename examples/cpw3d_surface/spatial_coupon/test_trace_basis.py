# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import csv
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import numpy as np

from trace_basis import (basis_statistics, cut_surface_size_report, load_trace_basis,
                         point_triangle_distances, process_frame, requested_sizes,
                         trace_basis_sizes, transform_trace_basis)

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
FOUR = HERE / "testdata" / "four-edge-9d2cb9bbb3fe"
TEN = HERE / "testdata" / "ten-edge-6791f1c84123"


def load(source):
    return load_trace_basis(source / "basis-contract.json", source / "trace-vertices.csv",
                            source / "trace-triangles.csv", source / "process-library.json")


def rotation(axis, angle):
    axis = np.asarray(axis, dtype=float); axis /= np.linalg.norm(axis)
    cross = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
    return np.eye(3) * np.cos(angle) + (1 - np.cos(angle)) * np.outer(axis, axis) + np.sin(angle) * cross


class TraceBasisTest(unittest.TestCase):
    def test_frame_is_the_campaign_producer_frame_and_places_the_basis_on_the_box(self):
        from generate_spatial_response import frame_from_geometry
        library = json.loads((FOUR / "process-library.json").read_text())
        model = library["Models"][0]
        expected = frame_from_geometry(model["Topology"], {"Edges": model["Edges"]})
        np.testing.assert_array_equal(process_frame(library, model["Name"]), expected)
        basis = load(FOUR)
        contract = json.loads((FOUR / "basis-contract.json").read_text())
        self.assertEqual(len(basis["Points"]), contract["Geometry"]["Vertices"])
        self.assertEqual(len(basis["Triangles"]), contract["Geometry"]["Triangles"])
        np.testing.assert_array_equal(basis["Lower"], contract["Geometry"]["Lower"])
        np.testing.assert_array_equal(basis["Upper"], contract["Geometry"]["Upper"])
        # Every vertex on the box, every triangle in one face plane, every source index once.
        on_face = (np.abs(basis["Points"] - basis["Lower"]) < 1e-9) | (np.abs(basis["Points"] - basis["Upper"]) < 1e-9)
        self.assertTrue(np.all(on_face.any(axis=1)))
        self.assertTrue(np.all(on_face[basis["Triangles"]].all(axis=1).any(axis=1)))
        self.assertEqual(sorted(basis["Basis"][basis["Basis"] > 0]), list(range(1, contract["Sources"] + 1)))
        ten = load(TEN)
        self.assertEqual((len(ten["Points"]), len(ten["Triangles"])), (260, 516))
        for role in ("BasisContract", "TraceVertices", "TraceTriangles", "ProcessLibrary"):
            self.assertEqual(len(basis["InputSHA256"][role]), 64)

    def test_loading_fails_closed_on_frame_count_and_connectivity_tampering(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in ("basis-contract.json", "trace-vertices.csv", "trace-triangles.csv",
                         "process-library.json"):
                shutil.copy(FOUR / name, root / name)
            load(root)
            # A vertex moved off the box (the frame would no longer place it on a face).
            with (FOUR / "trace-vertices.csv").open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            rows[5]["x"] = repr(float(rows[5]["x"]) + 0.3)
            with (root / "trace-vertices.csv").open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
            with self.assertRaisesRegex(ValueError, "box surface|box face plane"):
                load(root)
            shutil.copy(FOUR / "trace-vertices.csv", root / "trace-vertices.csv")
            # A dropped triangle contradicts the contract's geometry counts.
            lines = (FOUR / "trace-triangles.csv").read_text().splitlines()
            (root / "trace-triangles.csv").write_text("\n".join(lines[:-1]) + "\n")
            with self.assertRaisesRegex(ValueError, "geometry counts"):
                load(root)
            shutil.copy(FOUR / "trace-triangles.csv", root / "trace-triangles.csv")
            # A process library without the contract's model has no frame.
            library = json.loads((FOUR / "process-library.json").read_text())
            library["Models"][0]["Name"] = "other"
            (root / "process-library.json").write_text(json.dumps(library))
            with self.assertRaisesRegex(ValueError, "exactly once"):
                load(root)
            # The wrong frame (ten-edge library model name is different): rejected the same way.
            shutil.copy(TEN / "process-library.json", root / "process-library.json")
            with self.assertRaisesRegex(ValueError, "exactly once"):
                load(root)
            shutil.copy(FOUR / "process-library.json", root / "process-library.json")
            contract = json.loads((FOUR / "basis-contract.json").read_text())
            contract["Geometry"]["OffBoxTriangles"] = 1
            (root / "basis-contract.json").write_text(json.dumps(contract))
            with self.assertRaisesRegex(ValueError, "Unsupported"):
                load(root)

    def test_point_triangle_distance_is_exact_in_every_region(self):
        triangle = np.array([[0., 0., 0.], [2., 0., 0.], [0., 1., 0.]])
        query = np.array([[.5, .25, 3.],      # above the interior: plane distance
                          [-1., 0., 0.],      # beyond vertex a along the edge line
                          [1., -2., 0.],      # below edge ab
                          [3., 3., 0.],       # beyond the hypotenuse
                          [-1., -1., 1.]])    # vertex region, off plane
        distances = point_triangle_distances(query, triangle)
        np.testing.assert_allclose(distances[:3], [3., 1., 2.])
        # Hypotenuse from (2,0) to (0,1): the foot of (3,3) is (1.6, 0.2) -> sqrt(9.8);
        # (-1,-1,1) is nearest to vertex a at sqrt(3).
        np.testing.assert_allclose(distances[3], np.sqrt(9.8))
        np.testing.assert_allclose(distances[4], np.sqrt(3.))
        rng = np.random.default_rng(7)
        tri = rng.normal(size=(3, 3)); points = 3 * rng.normal(size=(50, 3))
        s, t = np.meshgrid(np.linspace(0, 1, 401), np.linspace(0, 1, 401), indexing="ij")
        keep = s + t <= 1
        samples = tri[0] + s[keep, None] * (tri[1] - tri[0]) + t[keep, None] * (tri[2] - tri[0])
        brute = np.array([np.linalg.norm(samples - p, axis=1).min() for p in points])
        exact = point_triangle_distances(points, tri)
        # The exact distance is never above the sampled one and within the sampling step.
        self.assertTrue(np.all(exact <= brute + 1e-12))
        step = np.max(np.linalg.norm(tri[[1, 2, 2]] - tri[[0, 0, 1]], axis=1)) / 400
        self.assertLess(np.max(brute - exact), 2 * step)

    def test_size_rule_is_per_triangle_graded_and_rigid_covariant(self):
        basis = load(FOUR)
        statistics = basis_statistics(basis, 1.0, 0.16)
        self.assertEqual(statistics["UniqueEdges"], 234)
        self.assertEqual(statistics["BasisEdgesBelowFarSize"], 46)
        self.assertAlmostEqual(statistics["MinimumBasisEdge"], 0.0217069561, places=9)
        self.assertEqual(statistics["TrianglesBelowFarSize"], 76)
        requested = requested_sizes(basis["Points"], basis["Triangles"], 1.0)
        narrow = int(np.argmin(requested))
        centroid = basis["Points"][basis["Triangles"][narrow]].mean(axis=0)
        normal = np.zeros(3); face = np.flatnonzero(np.abs(basis["Points"][basis["Triangles"][narrow]] - centroid).max(axis=0) < 1e-12)
        normal[face[0]] = 1.0
        inward = -np.sign(centroid[face[0]] - 0.5 * (basis["Lower"][face[0]] + basis["Upper"][face[0]]))
        query = np.array([centroid, centroid + inward * normal * 0.05, centroid + inward * normal * 5.0,
                          0.5 * (basis["Lower"] + basis["Upper"])])
        sizes = trace_basis_sizes(query, basis, 1.0, 0.16, 1.0)
        self.assertAlmostEqual(sizes[0], requested[narrow])
        self.assertAlmostEqual(sizes[1], requested[narrow] + 0.05)  # growth 1 along the normal
        self.assertEqual(sizes[2], 0.16); self.assertEqual(sizes[3], 0.16)
        # Doubling the ratio doubles the surface size where it stays below the far size.
        self.assertAlmostEqual(trace_basis_sizes(query[:1], basis, 2.0, 0.16, 1.0)[0], 2 * requested[narrow])
        # A wide-hat point (a bottom-face triangle with 2 um edges) keeps the far size.
        wide = int(np.argmax(requested))
        far_point = basis["Points"][basis["Triangles"][wide]].mean(axis=0)
        self.assertEqual(trace_basis_sizes(far_point, basis, 1.0, 0.16, 1.0)[0], 0.16)
        # Rigid covariance of sizes and of the placed basis.
        q = rotation([1., 2., -1.], 0.9); shift = np.array([3., -2., 1.])
        matrix = np.eye(4); matrix[:3, :3] = q; matrix[:3, 3] = shift
        placed = transform_trace_basis(basis, matrix.reshape(-1))
        np.testing.assert_allclose(trace_basis_sizes(query @ q.T + shift, placed, 1.0, 0.16, 1.0),
                                   sizes, rtol=1e-12, atol=1e-12)
        with self.assertRaises(ValueError):
            requested_sizes(basis["Points"], basis["Triangles"], 0.0)
        with self.assertRaises(ValueError):
            trace_basis_sizes(query, basis, 1.0, 0.16, 0.0)

    def test_cut_surface_report_measures_compliance_across_the_shortest_edge(self):
        # One box face z=0 over [0,4]x[0,1] as two basis triangles split by the diagonal;
        # shortest edges are the x=0 / x=4 sides (length 1), so extents are measured along y.
        points = np.array([[0., 0., 0.], [4., 0., 0.], [4., 1., 0.], [0., 1., 0.]])
        basis = {"Points": points, "Triangles": np.array([[0, 1, 2], [0, 2, 3]]),
                 "Lower": np.array([0., 0., 0.]), "Upper": np.array([4., 1., 3.])}
        # A cut mesh of the face with elements 0.5 x 0.25 (aligned rectangles split).
        xs, ys = np.linspace(0, 4, 9), np.linspace(0, 1, 5)
        grid = np.array([[x, y, 0.] for y in ys for x in xs])
        triangles = []
        for j in range(4):
            for i in range(8):
                a = j * 9 + i; b = a + 1; c = a + 10; d = a + 9
                triangles += [[a, b, c], [a, c, d]]
        report = cut_surface_size_report(grid, np.array(triangles), basis, ratio=0.5, far_size=10.)
        self.assertEqual(report["CutTriangles"], 64)
        self.assertAlmostEqual(report["Minimum"], np.hypot(.5, .25))
        rows = report["BasisTrianglesBelowFarSize"]
        self.assertEqual([row["RequestedSize"] for row in rows], [0.5, 0.5])
        # Every cut triangle is centred in a basis triangle (those centred on the shared
        # diagonal count for both).
        self.assertGreaterEqual(sum(row["CutTriangles"] for row in rows), 64)
        self.assertTrue(all(row["CutTriangles"] >= 24 for row in rows))
        for row in rows:
            self.assertAlmostEqual(row["MaximumExtentAlongShortestEdge"], 0.25)
            self.assertAlmostEqual(row["ExtentOverRequested"], 0.5)
        self.assertAlmostEqual(report["MaximumExtentOverRequested"], 0.5)


class TraceBasisJuliaTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.julia = shutil.which(os.environ.get("JULIA", "julia"))
        if cls.julia is None:
            raise unittest.SkipTest("Julia is not available")
        cls.project = REPO / "test" / "examples"

    def test_julia_seed_rule_matches_the_python_metric_rule(self):
        """The seed's Julia implementation of the trace basis (frame, triangles, per-triangle
        rule with the process-band slope) equals the Python metric implementation."""
        basis = load(FOUR)
        rng = np.random.default_rng(3)
        query = np.vstack((basis["Points"][basis["Triangles"]].mean(axis=1),
                           basis["Lower"] + rng.random((40, 3)) * (basis["Upper"] - basis["Lower"])))
        far, slope, ratio, lc = 0.6, 0.25, 4.0, 0.55
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            np.savetxt(root / "query.csv", query, delimiter=",")
            j = json.dumps  # Julia string literals with double quotes
            code = f"""
include(joinpath({j(str(HERE))}, "mesh_spatial_coupon.jl"))
using DelimitedFiles
basis = read_trace_basis({j(str(FOUR / 'basis-contract.json'))}, {j(str(FOUR / 'trace-vertices.csv'))},
                         {j(str(FOUR / 'trace-triangles.csv'))}, {j(str(FOUR / 'process-library.json'))})
record = prepare_trace_basis_sizing!(basis, {ratio}, {far}, {slope})
query = readdlm({j(str(root / 'query.csv'))}, ',')
sizes = [trace_basis_size(query[i, 1], query[i, 2], query[i, 3], {lc}) for i in axes(query, 1)]
writedlm({j(str(root / 'sizes.csv'))}, sizes)
open({j(str(root / 'record.json'))}, "w") do io; write_json(io, record); end
"""
            result = subprocess.run([self.julia, "--startup-file=no", f"--project={self.project}",
                                     "-e", code], cwd=REPO, capture_output=True, text=True,
                                    check=False, timeout=600)
            self.assertEqual(result.returncode, 0, f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}")
            julia_sizes = np.loadtxt(root / "sizes.csv")
            record = json.loads((root / "record.json").read_text())
        expected = np.minimum(lc, trace_basis_sizes(query, basis, ratio, far, slope))
        np.testing.assert_allclose(julia_sizes, expected, rtol=1e-12, atol=1e-12)
        np.testing.assert_array_equal(np.asarray(record["MeshFrameTriangles"]),
                                      basis["Points"][basis["Triangles"]])
        np.testing.assert_array_equal(np.asarray(record["Frame"]), basis["Frame"])
        self.assertEqual(record["InputSHA256"], basis["InputSHA256"])
        statistics = basis_statistics(basis, ratio, far)
        for key in ("UniqueEdges", "BasisEdgesBelowFarSize", "TrianglesBelowFarSize", "Triangles", "Vertices"):
            self.assertEqual(record[key], statistics[key], key)
        self.assertAlmostEqual(record["MinimumRequestedSize"], statistics["MinimumRequestedSize"])
        self.assertTrue(record["RatioIsDimensionless"]); self.assertEqual(record["Ratio"], ratio)


if __name__ == "__main__":
    unittest.main()
