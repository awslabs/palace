#!/usr/bin/env python3

"""ma_shells_2d.py: per-shell MA and the spatial coupons' sharp-edge estimator on a 2D
StraightEdgeBuilder response run (decision 66 part C), and the --edge-distances /
MA edge-point rule of the 2D response generators."""

import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

CPW2D = Path(__file__).resolve().parent
sys.path.insert(0, str(CPW2D))
import ma_shells_2d  # noqa: E402


def load(name):
    spec = importlib.util.spec_from_file_location(name, CPW2D / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RADII_UM = [0.002, 0.006, 0.014, 0.030, 0.062, 0.2]
LENGTH = 1.0e-3  # the edge length carrying the rings (arbitrary, cancels in the deficit)


def ring_energy(alpha, inner, outer, c=1.0):
    """c L (r_k^(1+a) - r_(k-1)^(1+a)) / (1 + a): the exact power-law ring energy."""
    return c * LENGTH * (outer ** (1.0 + alpha) - inner ** (1.0 + alpha)) / (1.0 + alpha)


def write_run(directory, *, sources, edges, ring1_factor=1.0):
    """A synthetic postpro directory: domain energies and, per edge (alpha, scale), the
    cumulative localized MA energies of an exact power law whose ring 1 carries only
    `ring1_factor` of the law (the unresolved singularity); Q_total = the sum of every
    ring plus a far remainder."""
    directory.mkdir(parents=True, exist_ok=True)
    with (directory / "domain-response-matrix.csv").open("w") as stream:
        stream.write("basis_i, basis_j, Q_ij (J)\n")
        for i, energy in sources.items():
            stream.write(f"{i}, {i}, {energy}\n")
    columns = ["Q_ij (J)", "Q_ij normal (J)", "Q_ij tangential (J)", "Q_total_ij (J)", "Q_total_ij normal (J)",
               "Q_total_ij tangential (J)"]
    with (directory / "surface-response-matrix.csv").open("w") as stream:
        stream.write("interface, edge, R (m), basis_i, basis_j, " + ", ".join(columns) + "\n")
        for i in sources:
            for interface in (1, 2, 3):
                for edge, (alpha, scale) in enumerate(edges, start=1):
                    rings = [ring_energy(alpha, inner, outer, scale) for inner, outer in zip([0.0] + RADII_UM, RADII_UM)]
                    rings[0] *= ring1_factor
                    total = math.fsum(rings) + 0.1 * rings[-1]
                    cumulative = 0.0
                    for radius, ring in zip(RADII_UM, rings):
                        cumulative += ring
                        value = cumulative if interface == 3 else 0.5 * cumulative
                        row = [value, value, 0.0, total if interface == 3 else 0.5 * total, total, 0.0]
                        stream.write(f"{interface}, {edge}, {radius * 1e-6:.12e}, {i}, {i}, " + ", ".join(f"{v:.16e}" for v in row) + "\n")


class MAShells2DTest(unittest.TestCase):
    def test_fabricated_rings_kinds_and_estimator(self):
        # Edge 1 (bottom corner) alpha -0.33, edge 2 (top corner) alpha -2/3; ring 1 holds
        # 70% of the law: the Consistent estimator (top Theory@2, bottom Fit2-4) recovers
        # the remainder of each kind exactly on exact power laws.
        with tempfile.TemporaryDirectory() as temp:
            postpro = Path(temp) / "postpro"
            write_run(postpro, sources={1: 2.0e-20, 2: 3.0e-20}, edges=[(-0.33, 2.0), (-2.0 / 3.0, 1.0)], ring1_factor=0.7)
            record = ma_shells_2d.analyze(postpro, kind="fabricated")
        self.assertEqual(record["RingRadii"], RADII_UM[:-1])
        self.assertEqual(record["Cutoff"], 0.002)
        source = record["PerSource"]["1"]
        self.assertEqual(set(source["Kinds"]), {"bottom", "top"})
        bottom, top = source["Kinds"]["bottom"], source["Kinds"]["top"]
        self.assertEqual(bottom["Estimator"], "Fit2-4")
        self.assertEqual(top["Estimator"], "Theory@2")
        self.assertAlmostEqual(top["Rings"]["1"]["Q"] / (0.7 * ring_energy(-2.0 / 3.0, 0.0, 0.002)), 1.0, places=9)
        self.assertAlmostEqual(top["Rings"]["2"]["Q"] / ring_energy(-2.0 / 3.0, 0.002, 0.006), 1.0, places=9)
        self.assertEqual(len(top["Rings"]), 6)  # 5 tube rings + the far ring to the coupon radius
        expected_top = 0.3 * ring_energy(-2.0 / 3.0, 0.0, 0.002)
        expected_bottom = 0.3 * ring_energy(-0.33, 0.0, 0.002, 2.0)
        self.assertAlmostEqual(top["Remainder"] / expected_top, 1.0, places=6)
        self.assertAlmostEqual(bottom["Remainder"] / expected_bottom, 1.0, places=3)
        self.assertAlmostEqual(source["Q_MA_tail"] / (expected_top + expected_bottom), 1.0, places=3)
        self.assertAlmostEqual(source["Q_MA_sharp"], source["Q_MA_raw"] + source["Q_MA_tail"])
        self.assertAlmostEqual(source["Alpha"], -2.0 / 3.0, places=3)
        self.assertAlmostEqual(source["Ring1Factor"], 0.7, places=6)
        self.assertGreater(source["Deficit"]["Consistent"], 0.0)
        self.assertEqual(record["Summary"]["Sources"], 2)
        self.assertEqual(record["Summary"]["Estimator"], "Consistent")
        self.assertAlmostEqual(record["Summary"]["p_MA"]["SharpMedian"] / record["Summary"]["p_MA"]["RawMedian"],
                               1.0 + record["Summary"]["Deficit"]["Consistent"]["Median"], places=6)
        self.assertIn("deficit (Consistent) median", ma_shells_2d.markdown(record))

    def test_thin_reports_cutoff_without_tail(self):
        with tempfile.TemporaryDirectory() as temp:
            postpro = Path(temp) / "postpro"
            write_run(postpro, sources={1: 1.0e-20}, edges=[(-1.0 + 1e-3, 1.0)])
            record = ma_shells_2d.analyze(postpro, kind="thin")
        source = record["PerSource"]["1"]
        self.assertEqual(list(source["Kinds"]), ["sheet"])
        self.assertNotIn("Q_MA_tail", source)
        self.assertEqual(source["Cutoff"], 0.002)
        self.assertAlmostEqual(record["Summary"]["LocalSlopeMedian"], -1.0, places=1)
        self.assertIn("no extrapolation", record["Rule"])

    def test_fabricated_needs_bottom_top_pairs_and_ring_radii(self):
        with tempfile.TemporaryDirectory() as temp:
            postpro = Path(temp) / "postpro"
            write_run(postpro, sources={1: 1.0e-20}, edges=[(-0.5, 1.0)])
            with self.assertRaisesRegex(ValueError, "bottom, top"):
                ma_shells_2d.analyze(postpro, kind="fabricated")
            # A historical run (0.2 alone) carries no rings.
            text = (postpro / "surface-response-matrix.csv").read_text().splitlines()
            keep = [text[0]] + [line for line in text[1:] if "2.000000000000e-07" in line]
            (postpro / "surface-response-matrix.csv").write_text("\n".join(keep) + "\n")
            with self.assertRaisesRegex(ValueError, "--edge-distances"):
                ma_shells_2d.analyze(postpro, kind="fabricated")

    def test_cli_writes_records(self):
        with tempfile.TemporaryDirectory() as temp:
            postpro = Path(temp) / "postpro"
            write_run(postpro, sources={1: 1.0e-20}, edges=[(-0.33, 1.0), (-2.0 / 3.0, 1.0)], ring1_factor=0.8)
            out_json, out_md = Path(temp) / "shells.json", Path(temp) / "shells.md"
            ma_shells_2d.main(["--postpro", str(postpro), "--kind", "fabricated", "--out-json", str(out_json),
                               "--out-md", str(out_md)])
            record = json.loads(out_json.read_text())
            self.assertEqual(record["Kind"], "fabricated")
            self.assertIn("MA shells - fabricated", out_md.read_text())


class EdgeDistancesOptionTest(unittest.TestCase):
    def test_ma_edge_points_follow_the_shells(self):
        for name, foot, sidewall in (("generate_edge_response", [2], [4]), ("generate_edge_pair_response", [2, 7], [5, 9])):
            module = load(name)
            self.assertEqual(module.ma_edge_attributes(foot, sidewall), foot)
            module.EDGE_DISTANCES[:] = RADII_UM
            try:
                self.assertEqual(module.ma_edge_attributes(foot, sidewall), sidewall)
            finally:
                module.EDGE_DISTANCES[:] = module.DEFAULT_EDGE_DISTANCES

    def test_isolated_config_carries_the_rings_on_every_interface(self):
        module = load("generate_edge_response")
        module.EDGE_DISTANCES[:] = RADII_UM
        try:
            layers = {"SA": (0.002, 4.0), "MS": (0.002, 11.47), "MA": (0.002, 10.0)}
            fabricated = module.make_config(Path("/out"), "edge_fabricated", Path("/m.msh"), [Path("/t.csv")], True, 2,
                                            1055.0, 11.45, layers)
            thin = module.make_config(Path("/out"), "edge_thin", Path("/m.msh"), [Path("/t.csv")], False, 2, 1055.0, 11.45,
                                      layers)
        finally:
            module.EDGE_DISTANCES[:] = module.DEFAULT_EDGE_DISTANCES
        dielectrics = {entry["Type"]: entry for entry in fabricated["Boundaries"]["Postprocessing"]["Dielectric"]}
        self.assertEqual(dielectrics["MA"]["EdgeAttributes"], [4])
        self.assertEqual(dielectrics["MS"]["EdgeAttributes"], [2])
        self.assertTrue(all(entry["EdgeDistances"] == RADII_UM for entry in dielectrics.values()))
        thin_ma = next(entry for entry in thin["Boundaries"]["Postprocessing"]["Dielectric"] if entry["Type"] == "MA")
        self.assertEqual(thin_ma["EdgeAttributes"], [2])

    def test_edge_distances_must_end_at_the_coupon_radius(self):
        for script, extra in (("generate_edge_response.py", []), ("generate_edge_pair_response.py", ["--separation", "1.0"])):
            with tempfile.TemporaryDirectory() as temp:
                command = [sys.executable, str(CPW2D / script), "--library-name", "x", "--output", temp, *extra,
                           "--edge-distances", "0.002", "0.1"]
                result = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(result.returncode, 0)
            self.assertIn("coupon radius 0.2", result.stderr)


if __name__ == "__main__":
    unittest.main()
