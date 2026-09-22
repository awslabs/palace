# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Sharp-edge MA extrapolation in qualify (decision 61a; qualify/ma_tail.py): the tail
of synthetic power-law shells, the modelled reference side, the sharp offsets in the
comparison, the p-sequence observable and the gate quantities; the stored radial run
(/tmp/library-radial-ma-01, RADIAL-MA.md) reproduces its recorded deficits."""
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "qualify"))
import compare_matrices  # noqa: E402
import gates  # noqa: E402
import ma_tail  # noqa: E402
import p_sequence  # noqa: E402
import radial_ma_profile  # noqa: E402

RADIAL_RUN = Path("/tmp/library-radial-ma-01")
RADIAL_CASE = "two-edge-calib-radial-ma-shells"
RADII = [0.00025, 0.00075, 0.00175, 0.00375, 0.00775, 0.01575, 0.03175]
SURFACE_HEADER = ("interface,edge,R (m),basis_i,basis_j,Q_ij (J),Q_ij normal (J),Q_ij tangential (J),Q_total_ij (J),"
                  "Q_total_ij normal (J),Q_total_ij tangential (J)")


def shell_map(radii=RADII):
    """Interface index -> shell record as case_inputs.expand_radial_shells writes it:
    index 3 far, 4..10 top rings 1..7, 11..17 bottom rings 1..7; MS at 2."""
    out = {3: {"Kind": "far", "Ring": 0, "InnerRadius": radii[-1], "OuterRadius": None}}
    index = 4
    for kind in ("top", "bottom"):
        for ring, outer in enumerate(radii, start=1):
            out[index] = {"Kind": kind, "Ring": ring, "InnerRadius": 0.0 if ring == 1 else radii[ring - 2], "OuterRadius": outer}
            index += 1
    return out


def ring_energy(alpha, inner, outer, c=1.0):
    return c * radial_ma_profile.ring_energy_factor(alpha, inner, outer)


def write_reducer(directory, sources, *, ring1_factor, alpha_top=-2.0 / 3.0, alpha_bottom=-1.0 / 3.0, far=0.3, e=2.0,
                  bottom_scale=0.2):
    """A reducer directory whose top-edge shells follow r^alpha_top exactly with ring 1
    resolved at `ring1_factor` of the law, the bottom edge r^alpha_bottom (ring 1 exact)."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    domain = ["basis_i,basis_j,Q_ij (J)"]
    surface = [SURFACE_HEADER]
    shells = shell_map()
    for i in sources:
        domain.append(f"{i},{i},{e}")
        for index, shell in shells.items():
            if shell["Kind"] == "far":
                q = far
            elif shell["Kind"] == "top":
                q = ring_energy(alpha_top, shell["InnerRadius"], shell["OuterRadius"]) * (ring1_factor if shell["Ring"] == 1 else 1.0)
            else:
                q = bottom_scale * ring_energy(alpha_bottom, shell["InnerRadius"], shell["OuterRadius"])
            q *= i  # source-dependent scale
            surface.append(f"{index},1,2e-06,{i},{i},{q},{q},0,{q},{q},0")
        ms = 0.1 * i
        surface.append(f"2,1,2e-06,{i},{i},{ms},{ms},0,{ms},{ms},0")
    (directory / "domain-response-matrix.csv").write_text("\n".join(domain) + "\n")
    (directory / "surface-response-matrix.csv").write_text("\n".join(surface) + "\n")


def write_reference(directory, sources, p_ma, e=2.0):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    domain = ["basis_i,basis_j,Q_ij (J)"] + [f"{i},{i},{e}" for i in sources]
    surface = [SURFACE_HEADER]
    for i in sources:
        q = p_ma[i] * e
        surface.append(f"1,1,2e-06,{i},{i},{q},{q},0,{q},{q},0")
        surface.append(f"2,1,2e-06,{i},{i},{0.1 * i},{0.1 * i},0,{0.1 * i},{0.1 * i},0")
    (directory / "domain-response-matrix.csv").write_text("\n".join(domain) + "\n")
    (directory / "surface-response-matrix.csv").write_text("\n".join(surface) + "\n")


INTERFACE_TYPES = {index: "MA" for index in shell_map()} | {2: "MS"}
REFERENCE_TYPES = {1: "MA", 2: "MS"}


class MATailRuleTest(unittest.TestCase):
    def setUp(self):
        self.tmp = Path(tempfile.mkdtemp())
        self.addCleanup(lambda: __import__("shutil").rmtree(self.tmp, True))
        self.sources = [1, 2, 3]

    def test_tail_recovers_the_exact_law_and_the_ring1_factor(self):
        write_reducer(self.tmp / "run", self.sources, ring1_factor=0.7)
        tails = ma_tail.tails(self.tmp / "run", shell_map(), INTERFACE_TYPES)
        self.assertEqual(sorted(tails), self.sources)
        for i in self.sources:
            entry = tails[i]
            # Theory@2 on an exact -2/3 law: the ring-1 model is the exact ring-1 energy, the
            # tail is the 30% the resolved ring 1 misses; the bottom edge (exact ring 1) adds nothing.
            exact_ring1 = ring_energy(-2.0 / 3.0, 0.0, RADII[0]) * i
            self.assertAlmostEqual(entry["Ring1Factor"], 0.7, places=12)
            self.assertAlmostEqual(entry["Kinds"]["top"]["Remainder"], 0.3 * exact_ring1, places=14)
            self.assertAlmostEqual(entry["Kinds"]["bottom"]["Remainder"], 0.0, places=12)
            self.assertAlmostEqual(entry["Q_MA_tail"], 0.3 * exact_ring1, places=14)
            self.assertAlmostEqual(entry["Q_MA_sharp"], entry["Q_MA_raw"] + entry["Q_MA_tail"])
            self.assertAlmostEqual(entry["p_MA_sharp"], entry["Q_MA_sharp"] / 2.0)
            self.assertAlmostEqual(entry["Alpha"], -2.0 / 3.0, places=3)
            self.assertAlmostEqual(entry["Deficit"]["Consistent"], entry["Q_MA_tail"] / entry["Q_MA_raw"])
            # Fit2-4 on the exact law agrees with Theory@2 (the estimator spread vanishes).
            self.assertAlmostEqual(entry["Deficit"]["Fit2-4"], entry["Deficit"]["Consistent"], places=6)
        summary = ma_tail.summary(tails, self.sources, strongest=[3])
        self.assertAlmostEqual(summary["Ring1Factor"]["Median"], 0.7, places=12)
        self.assertAlmostEqual(summary["Alpha"]["Median"], -2.0 / 3.0, places=3)
        self.assertEqual(summary["Deficit"]["Consistent"]["At"], {})
        self.assertAlmostEqual(summary["Deficit"]["Consistent"]["StrongestMedian"], tails[3]["Deficit"]["Consistent"])
        # The shell map must be the MA interfaces of the run config.
        with self.assertRaisesRegex(ValueError, "shell interfaces"):
            ma_tail.tails(self.tmp / "run", shell_map(), INTERFACE_TYPES | {99: "MA"})

    def test_reference_side_modelled_extrapolated_and_raw(self):
        write_reducer(self.tmp / "run", self.sources, ring1_factor=0.7)
        tails = ma_tail.tails(self.tmp / "run", shell_map(), INTERFACE_TYPES)
        ref_pma = {i: tails[i]["p_MA_raw"] * 0.99 for i in self.sources}
        # Modelled: deficit_ref = deficit_run x (0.5 / 0.25)^(1/3) from a 0.5 nm edge size.
        modelled = ma_tail.reference_side(ref_pma, tails, reference_edge_size=0.0005, run_inner_radius=RADII[0])
        self.assertTrue(modelled["Modelled"] and not modelled["Extrapolated"])
        for i in self.sources:
            expected = tails[i]["Deficit"]["Consistent"] * 2.0 ** (1.0 / 3.0)
            self.assertAlmostEqual(modelled["PerSource"][i]["ModelledDeficit"], expected)
            self.assertAlmostEqual(modelled["PerSource"][i]["p_MA_sharp"], ref_pma[i] * (1 + expected))
        self.assertIn("reference unextrapolated; sharp-edge deficit modelled at", modelled["Annotation"])
        self.assertIn("0.5 nm edge size via the eps^(1/3) law", modelled["Annotation"])
        # Extrapolated: a shelled reference carries its own tails.
        extrapolated = ma_tail.reference_side(ref_pma, tails, reference_tails=tails)
        self.assertTrue(extrapolated["Extrapolated"])
        self.assertAlmostEqual(extrapolated["PerSource"][2]["p_MA_sharp"], tails[2]["p_MA_sharp"])
        # Raw: no sizing and no edge size - the raw reference stands in, annotated.
        raw = ma_tail.reference_side(ref_pma, tails)
        self.assertFalse(raw["Modelled"] or raw["Extrapolated"])
        self.assertEqual(raw["PerSource"][1]["p_MA_sharp"], ref_pma[1])
        self.assertIn("not modelled", raw["Annotation"])

    def test_sharp_offsets_flow_through_comparison_p_sequence_and_gates(self):
        # The run's raw MA is 2.68% below its sharp value (ring-1 factor tuned so); the
        # reference's raw MA equals the run's SHARP value: raw-vs-raw shows the deficit,
        # sharp-vs-modelled shows the modelled reference deficit alone.
        write_reducer(self.tmp / "run", self.sources, ring1_factor=0.7)
        tails = ma_tail.tails(self.tmp / "run", shell_map(), INTERFACE_TYPES)
        ref_pma = {i: tails[i]["p_MA_sharp"] for i in self.sources}
        write_reference(self.tmp / "ref", self.sources, ref_pma)
        side = ma_tail.reference_side(ref_pma, tails, reference_edge_size=0.0005, run_inner_radius=RADII[0])
        summary = compare_matrices.write_comparison(self.tmp / "ref", self.tmp / "run", self.tmp / "cmp", zero_trace_indices=[],
                                                    interface_types=INTERFACE_TYPES, reference_interface_types=REFERENCE_TYPES,
                                                    run_ma_tails=tails, reference_ma_side=side)
        for i in self.sources:
            record = summary["PerSource"][str(i)]
            deficit = tails[i]["Deficit"]["Consistent"]
            self.assertAlmostEqual(record["p_MA_rel"], 1.0 / (1.0 + deficit) - 1.0, places=12)
            self.assertAlmostEqual(record["p_MA_sharp_rel"], 1.0 / (1.0 + deficit * 2.0 ** (1.0 / 3.0)) - 1.0, places=12)
        # Without the tails the comparison carries the raw offsets only.
        plain = compare_matrices.write_comparison(self.tmp / "ref", self.tmp / "run", self.tmp / "plain", zero_trace_indices=[],
                                                  interface_types=INTERFACE_TYPES, reference_interface_types=REFERENCE_TYPES)
        self.assertNotIn("p_MA_sharp_rel", plain["PerSource"]["1"])
        # p-sequence: the sharp observable per order from that order's own shells; a
        # higher order resolving ring 1 better changes p_MA but not p_MA_sharp (exact law).
        write_reducer(self.tmp / "high", self.sources, ring1_factor=0.75)
        write_reducer(self.tmp / "low", self.sources, ring1_factor=0.65)
        runs = {"low": self.tmp / "low", "main": self.tmp / "run", "high": self.tmp / "high", "ref": self.tmp / "ref"}
        ma_tails = {key: ma_tail.tails(path, shell_map(), INTERFACE_TYPES) for key, path in runs.items() if key != "ref"}
        sequence = p_sequence.p_sequence(runs, [2], INTERFACE_TYPES, REFERENCE_TYPES, ma_tails=ma_tails, reference_ma_side=side)
        seq_raw, seq_sharp = sequence[2]["p_MA"]["seq"], sequence[2]["p_MA_sharp"]["seq"]
        self.assertGreater(abs(seq_raw["d_high"]), 1e-3)
        self.assertLess(abs(seq_sharp["d_high"]), 1e-12)
        self.assertAlmostEqual(sequence[2]["p_MA_sharp"]["values"]["ref"], side["PerSource"][2]["p_MA_sharp"])
        self.assertAlmostEqual(sequence[2]["Q_MA_sharp"]["values"]["main"], tails[2]["Q_MA_sharp"])
        # Without tails the sharp observables are absent (a run without shells).
        self.assertNotIn("p_MA_sharp", p_sequence.p_sequence(runs, [2], INTERFACE_TYPES, REFERENCE_TYPES)[2])
        # Gates: p_MA on p_MA_sharp_rel with the recorded rule and the p-sequence on p_MA_sharp.
        table, digest = gates.load_gates()
        self.assertEqual(table["Gates"]["p_MA"]["Quantity"], "p_MA_sharp")
        self.assertEqual(table["Gates"]["PSequenceControls"]["MAObservable"], "p_MA_sharp")
        record = gates.evaluate(table, comparison=summary, classes={i: "wide hats, bottom/top faces" for i in self.sources},
                                ref_pma=ref_pma, p_sequence_summary=sequence, reference_order=4, gates_sha256=digest,
                                interfaces=["MA", "MS"])
        ma_gate = record["Gates"]["p_MA"]
        self.assertEqual(ma_gate["Quantity"], "p_MA_sharp")
        self.assertEqual(ma_gate["QuantityRule"], table["Gates"]["p_MA"]["QuantityRule"])
        self.assertAlmostEqual(ma_gate["Free"]["signed_median"], summary["PerSource"]["2"]["p_MA_sharp_rel"])
        self.assertAlmostEqual(ma_gate["RawFree"]["signed_median"], summary["PerSource"]["2"]["p_MA_rel"])
        controls = record["Gates"]["PSequenceControls"]
        self.assertEqual(controls["MAObservable"], "p_MA_sharp")
        self.assertEqual(set(controls["Controls"]["2"]), {"E", "p_MA_sharp", "p_MS"})
        self.assertEqual(controls["NotApplicableObservables"], ["p_SA"])
        self.assertTrue(controls["Passed"])
        # A run without shells is gated on the raw p_MA and says so.
        raw_record = gates.evaluate(table, comparison=plain, classes={i: "wide hats, bottom/top faces" for i in self.sources},
                                    ref_pma=ref_pma, p_sequence_summary=p_sequence.p_sequence(runs, [2], INTERFACE_TYPES, REFERENCE_TYPES),
                                    reference_order=4, gates_sha256=digest, interfaces=["MA", "MS"])
        self.assertEqual(raw_record["Gates"]["p_MA"]["Quantity"], "p_MA")
        self.assertIn("no radial MA shells", raw_record["Gates"]["p_MA"]["QuantityRule"])
        self.assertIsNone(raw_record["Gates"]["p_MA"]["RawFree"])
        self.assertEqual(raw_record["Gates"]["PSequenceControls"]["MAObservable"], "p_MA")
        self.assertEqual(set(raw_record["Gates"]["PSequenceControls"]["Controls"]["2"]), {"E", "p_MA", "p_MS"})
        # The written record and CSV.
        tail_record = {"Rule": ma_tail.RULE, "Orders": {"p4": {"Prefix": "x-p4", "Strongest": [3],
                                                              "Summary": ma_tail.summary(tails, self.sources, strongest=[3])}},
                       "Reference": side}
        ma_tail.write_record(self.tmp / "ma-tail.json", self.tmp / "ma-tail.md", tail_record, "test")
        ma_tail.write_csv(self.tmp / "ma-sharp.csv", tails)
        text = (self.tmp / "ma-tail.md").read_text()
        self.assertIn("ring-1 factor", text)
        self.assertIn("reference unextrapolated", text)
        self.assertEqual(len((self.tmp / "ma-sharp.csv").read_text().splitlines()), 1 + len(self.sources))


@unittest.skipUnless((RADIAL_RUN / "library-qualification.json").is_file(), "the stored radial run is not available")
class StoredRadialRunTest(unittest.TestCase):
    """The decision-56 run reproduces RADIAL-MA.md: Consistent deficit median 2.68% (p5), 53 / 58
    4.44%; Fit2-4 1.89% / 4.75%; alpha -0.648; ring-1 factor 0.706 (0.687-0.710); p4 3.07% / 5.11%."""

    def test_reproduces_radial_ma_md(self):
        case = next(c for c in json.loads((RADIAL_RUN / "library-qualification.json").read_text())["Cases"] if c["Case"] == RADIAL_CASE)
        shells = ma_tail.shell_map_of(case["Inputs"]["RadialShells"])
        types = {int(k): v for k, v in case["Inputs"]["Interfaces"].items()}
        zero = set(case["Sources"]["ZeroTrace"])
        results = Path(case["Root"]) / "results" / "main"
        expected = {5: (0.0268, 0.0444, 0.0189, 0.0475, 0.706), 4: (0.0307, 0.0511, 0.0221, 0.0540, 0.665)}
        for order, (consistent, at53, fit, fit53, ring1) in expected.items():
            tails = ma_tail.tails(results / f"{RADIAL_CASE}-p{order}" / "reducer", shells, types)
            free = [i for i in tails if i not in zero]
            self.assertEqual(len(free), 78)
            summary = ma_tail.summary(tails, free)
            self.assertAlmostEqual(summary["Deficit"]["Consistent"]["Median"], consistent, places=4)
            self.assertAlmostEqual(summary["Deficit"]["Consistent"]["At"]["53"], at53, places=4)
            self.assertAlmostEqual(summary["Deficit"]["Consistent"]["At"]["58"], at53, places=4)
            self.assertAlmostEqual(summary["Deficit"]["Fit2-4"]["Median"], fit, places=4)
            self.assertAlmostEqual(summary["Deficit"]["Fit2-4"]["At"]["53"], fit53, places=4)
            self.assertAlmostEqual(summary["Alpha"]["Median"], -0.648, places=3)
            self.assertAlmostEqual(summary["Ring1Factor"]["Median"], ring1, places=3)
            if order == 5:
                self.assertAlmostEqual(summary["Ring1Factor"]["Range"][0], 0.687, places=3)
                self.assertAlmostEqual(summary["Ring1Factor"]["Range"][1], 0.710, places=3)
                self.assertAlmostEqual(summary["Alpha"]["At"]["53"][0], -0.672, places=3)
        # The modelled 0.5 nm reference deficit at 53 / 58: 4.44% x 2^(1/3) = 5.6% (RADIAL-MA.md).
        tails = ma_tail.tails(results / f"{RADIAL_CASE}-p5" / "reducer", shells, types)
        side = ma_tail.reference_side({53: 1.0, 58: 1.0}, tails, reference_edge_size=0.0005,
                                      run_inner_radius=case["Inputs"]["RadialShells"]["RingRadii"][0])
        self.assertAlmostEqual(side["PerSource"][53]["ModelledDeficit"], 0.0559, places=3)


CALIBRATION_BUILD = Path("/tmp/coupon-calibration-radial-ma-20260921/library-build.json")
REFERENCE_CAMPAIGN = Path("/tmp/library-acceptance-01-reference")
BINARY_SHA256 = "b28f089ae12c25863493566b2b8ca11af2c8ffb0e273e7aa67a2b42046eacf27"


@unittest.skipUnless(CALIBRATION_BUILD.is_file() and (RADIAL_RUN / RADIAL_CASE / "results" / "main").is_dir()
                     and REFERENCE_CAMPAIGN.is_dir(), "the stored radial run, its build record and the reference are not available")
class StoredRadialRunAnalysisTest(unittest.TestCase):
    """qualify's analysis of the stored radial run (results read only) gates p_MA on p_MA_sharp
    against the 0.5 nm reference modelled by the eps^(1/3) law and the p-sequence on p_MA_sharp;
    the raw offsets are RADIAL-MA.md's (median +1.25%, 53 +3.50%), the sharp-vs-modelled ones
    smaller; process-library carries MA_raw and MA_sharp per source."""

    def test_analysis_gates_on_the_sharp_ma(self):
        import argparse
        import estimate_stages
        import qualify_library
        manifest_path = HERE / "geometry-independence-calibration-sizing.json"
        manifest = json.loads(manifest_path.read_text())
        manifest["Path"] = str(manifest_path)
        build = json.loads(CALIBRATION_BUILD.read_text())
        case = next(item for item in build["Cases"] if item["Case"] == RADIAL_CASE)
        args = argparse.Namespace(reference=REFERENCE_CAMPAIGN, control_source=[1, 2, 3, 5, 6, 20, 25, 49, 53, 58], control_count=8,
                                  stage_prefix=RADIAL_CASE, orders=[4, 5], controls=[3, 5], frozen_binary_sha256=BINARY_SHA256,
                                  max_jobs=1, job_policy=None, fixed_jobs=None, reference_edge_size_nm=0.5)
        profile = json.loads((HERE / "qualify" / "cluster-profile.json").read_text())
        table, digest = gates.load_gates()
        results = RADIAL_RUN / RADIAL_CASE / "results"
        before = {path: path.stat().st_mtime_ns for path in results.rglob("*") if path.is_file()}
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            record, context = qualify_library.prepare_case(case, manifest_path=manifest_path, manifest=manifest, args=args, root=root,
                                                           remote={"Host": "h", "Root": "/r"}, profile=profile,
                                                           cost_model=estimate_stages.load_cost_model(), gates=table, gates_digest=digest)
            gate_record = qualify_library.analyze_case(record, context, results, gates=table, gates_digest=digest, profile=profile)
            self.assertEqual({path: path.stat().st_mtime_ns for path in results.rglob("*") if path.is_file()}, before)
            ma = gate_record["Gates"]["p_MA"]
            self.assertEqual(ma["Quantity"], "p_MA_sharp")
            self.assertEqual(gate_record["GatedOrder"], 5)
            self.assertAlmostEqual(ma["RawFree"]["signed_median"], 0.0125, places=3)
            self.assertAlmostEqual(ma["RawFree"]["worst_offset"], 0.0350, places=3)
            self.assertEqual(ma["RawFree"]["worst_source"], 53)
            self.assertLess(ma["Free"]["signed_median"], ma["RawFree"]["signed_median"])
            self.assertLess(ma["Free"]["worst_offset"], ma["RawFree"]["worst_offset"])
            self.assertGreater(ma["Free"]["within_2pct"], ma["RawFree"]["within_2pct"])
            controls = gate_record["Gates"]["PSequenceControls"]
            self.assertEqual(controls["MAObservable"], "p_MA_sharp")
            self.assertTrue(controls["Passed"])
            self.assertIn("p_MA_sharp", controls["Controls"]["53"])
            tail = record["MATail"]
            self.assertAlmostEqual(tail["Orders"]["p5"]["Summary"]["Deficit"]["Consistent"]["Median"], 0.0268, places=4)
            self.assertAlmostEqual(tail["Orders"]["p5"]["Summary"]["Deficit"]["Consistent"]["At"]["53"], 0.0444, places=4)
            self.assertTrue(tail["Reference"]["Modelled"])
            self.assertIn("reference unextrapolated; sharp-edge deficit modelled at", tail["Reference"]["Annotation"])
            self.assertEqual(record["Qualification"]["MAQuantity"], "p_MA_sharp")
            comparison = root / RADIAL_CASE / "comparison"
            for name in ("ma-tail.json", "ma-tail.md", f"ma-sharp-{RADIAL_CASE}-p4.csv", f"ma-sharp-{RADIAL_CASE}-p5.csv"):
                self.assertTrue((comparison / name).is_file(), name)
            offsets = json.loads((comparison / "ma-ms-offsets.json").read_text())
            gated = offsets["Distributions"][f"{RADIAL_CASE}-p5-vs-reference"]
            self.assertEqual(gated["p_MA_rel:strongest20"]["n"], 20)
            self.assertEqual(gated["p_MA_sharp_rel:strongest20"]["n"], 20)
            library = qualify_library.process_library_entries([record], {RADIAL_CASE: context}, manifest_path=manifest_path,
                                                              manifest=manifest, root=root)
            model = library["Models"][0]["MA"]
            self.assertEqual(len(model["MA_sharp"]), 78)
            self.assertGreater(model["MA_sharp"]["53"], model["MA_raw"]["53"])
            self.assertEqual(model["Order"], 4)


if __name__ == "__main__":
    unittest.main()
