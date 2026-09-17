# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Fail-closed publication tests for the optional native MMG adapter."""
import os
from pathlib import Path
import subprocess
import tempfile
import unittest


EXECUTABLE = os.environ.get("EDGE_METRIC_ADAPTER_EXE")


@unittest.skipUnless(EXECUTABLE and Path(EXECUTABLE).is_file(),
                     "EDGE_METRIC_ADAPTER_EXE is not available")
class EdgeMetricAdapterPublicationTest(unittest.TestCase):
    def invoke(self, root, output):
        return subprocess.run(
            [EXECUTABLE, str(root / "missing.meshb"), str(root / "missing.f64"),
             str(root / "missing.txt"), str(output), ".01", ".1", "1.2"],
            capture_output=True, text=True, timeout=10)

    def test_existing_solution_file_is_not_mutated(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "accepted.meshb"
            solution = Path(str(output) + ".sol")
            solution.write_bytes(b"immutable solution")
            result = self.invoke(root, output)
            self.assertEqual(result.returncode, 2)
            self.assertIn("Refusing to overwrite", result.stderr)
            self.assertEqual(solution.read_bytes(), b"immutable solution")
            self.assertFalse(output.exists())

    def test_unsavable_existing_solution_directory_is_not_mutated(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            output = root / "accepted.meshb"
            solution = Path(str(output) + ".sol")
            solution.mkdir()
            result = self.invoke(root, output)
            self.assertEqual(result.returncode, 2)
            self.assertIn("Refusing to overwrite", result.stderr)
            self.assertTrue(solution.is_dir())
            self.assertFalse(output.exists())

    def test_existing_mesh_and_rejected_attempt_are_not_mutated(self):
        for suffix in ("", ".rejected.meshb"):
            with self.subTest(suffix=suffix), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                output = root / "accepted.meshb"
                protected = Path(str(output) + suffix)
                protected.write_bytes(b"immutable attempt")
                result = self.invoke(root, output)
                self.assertEqual(result.returncode, 2)
                self.assertEqual(protected.read_bytes(), b"immutable attempt")

    @staticmethod
    def cube(root):
        """A Kuhn-subdivided unit cube (6 tets, 12 boundary triangles), an isotropic
        metric, one pin: the smallest complete adapter input."""
        import numpy as np
        import meshio
        points = np.array([[x, y, z] for z in (0., 1.) for y in (0., 1.) for x in (0., 1.)])
        tets = np.array([[0, 1, 3, 7], [0, 1, 5, 7], [0, 2, 3, 7], [0, 2, 6, 7],
                         [0, 4, 5, 7], [0, 4, 6, 7]])
        # Positively oriented like a Gmsh seed, so MMG keeps the vertex order.
        for cell in tets:
            edges = points[cell[1:]] - points[cell[0]]
            if np.linalg.det(edges) < 0:
                cell[2], cell[3] = cell[3], cell[2]
        faces = {}
        for cell in tets:
            for face in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
                key = tuple(sorted(cell[list(face)]))
                faces[key] = faces.get(key, 0) + 1
        triangles = np.array([face for face, count in faces.items() if count == 1])
        mesh = meshio.Mesh(points, [("triangle", triangles), ("tetra", tets)],
                           cell_data={"medit:ref": [np.full(len(triangles), 10001),
                                                    np.full(len(tets), 1)]})
        meshio.write(root / "seed.mesh", mesh, file_format="medit")
        np.tile(np.array([4., 0., 0., 4., 0., 4.]), len(points)).astype("<f8").tofile(
            root / "metric.f64")
        (root / "pins.txt").write_text("1\n")
        return len(tets)

    def test_required_tetrahedra_flag_fails_closed_and_keeps_listed_cells(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary); tets = self.cube(root)
            base = [EXECUTABLE, str(root / "seed.mesh"), str(root / "metric.f64"),
                    str(root / "pins.txt"), str(root / "out.meshb"), ".5", "2", "1.3"]
            def attempt(required, *extra):
                (root / "required.txt").write_text(required)
                result = subprocess.run(
                    [*base, *extra, "--required-tetrahedra", str(root / "required.txt")],
                    capture_output=True, text=True, timeout=60)
                return result
            for required, message in ((f"{tets + 1}\n", "Invalid required tetrahedron index"),
                                      ("0\n", "Invalid required tetrahedron index"),
                                      ("", "malformed required tetrahedron list"),
                                      ("1 x\n", "malformed required tetrahedron list")):
                with self.subTest(required=required):
                    result = attempt(required)
                    self.assertEqual(result.returncode, 2)
                    self.assertIn(message, result.stderr)
                    self.assertFalse((root / "out.meshb").exists())
            result = subprocess.run([*base, "--required-tetrahedra"], capture_output=True,
                                    text=True, timeout=60)
            self.assertEqual(result.returncode, 2)
            self.assertIn("exactly one file", result.stderr)
            result = attempt("1\n", "--required-tetrahedra", str(root / "required.txt"))
            self.assertEqual(result.returncode, 2)
            self.assertIn("exactly one file", result.stderr)
            # A valid list runs MMG and the listed seed cells survive verbatim.
            result = attempt("1\n2\n")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Required tetrahedra: 2", result.stdout)
            self.assertIn("required tetrahedra=2", result.stdout)
            import sys
            sys.path.insert(0, str(Path(__file__).resolve().parent))
            from mesh_array_io import read_mesh
            import numpy as np
            adapted = read_mesh(root / "out.meshb")
            flags = np.concatenate([flag for block, flag in
                                    zip(adapted.cells, adapted.cell_data["medit:required"])
                                    if block.type == "tetra"])
            cells = np.concatenate([block.data for block in adapted.cells if block.type == "tetra"])
            self.assertEqual(int(flags.sum()), 2)
            kept = {tuple(map(tuple, np.round(adapted.points[cell], 12))) for cell in cells[flags == 1]}
            points = np.array([[x, y, z] for z in (0., 1.) for y in (0., 1.) for x in (0., 1.)])
            seed_cells = np.array([[0, 1, 3, 7], [0, 1, 5, 7]])
            for cell in seed_cells:
                if np.linalg.det(points[cell[1:]] - points[cell[0]]) < 0:
                    cell[2], cell[3] = cell[3], cell[2]
            self.assertEqual(kept, {tuple(map(tuple, points[cell])) for cell in seed_cells})


if __name__ == "__main__":
    unittest.main()
