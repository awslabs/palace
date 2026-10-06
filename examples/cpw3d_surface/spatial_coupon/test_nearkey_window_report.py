# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""nearkey_window_report: the reused share per Type and the reuse error budget (MA of record on MA_sharp) from
a run's model-energy CSV + ModelCatalog + config and the library's reused models - a synthetic run (always),
and the three measured windows reproducing the DESIGN v2 4.4 budgets (S3p 0.057 / 0.034 / 0.113, S4 0.257 /
0.170 / 0.346, C3 0.794 / 0.539 / 0.717 points) when the local evidence mirror is mounted."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nearkey_window_report as report  # noqa: E402

S2 = Path("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/stage2-20261004")
DEMO_LIBRARY = S2 / "nearkey-reuse-impl" / "work" / "demo" / "s2-r1p9-v3-b1-nearkey-demo.json"
DERIVED = S2 / "nearkey-sensitivity" / "derived-sens.json"
RUNS = S2 / "nearkey-sensitivity" / "runs"
DESIGN_BUDGETS = {"sct002-S3p": (0.057, 0.034, 0.113), "sct002-S4": (0.257, 0.170, 0.346), "ctx003-C3": (0.794, 0.539, 0.717)}


def write_run(postpro, energies, classes):
    """A synthetic postpro: model-energy rows (evaluations 0 / 1 / 2; the fixed trace carries `energies`),
    palace.json ModelCatalog and config_resolved.json Dielectric entries."""
    postpro = Path(postpro)
    postpro.mkdir(parents=True)
    k_max = max(classes)
    header = "source,evaluation,model,patch count,patch weight,domain correction (J)," + ",".join(
        f"fabricated surface energy[{k}] (J)" for k in range(1, k_max + 1))
    lines = [header]
    for index, per_interface in energies.items():
        for evaluation in (0, 1, 2):
            scale = 1.0 if evaluation == 0 else 1.1
            values = ",".join(f"{scale * per_interface.get(k, 0.0):+.12e}" for k in range(1, k_max + 1))
            lines.append(f"1.00e+00,{evaluation:.2e},{index:.2e},1.00e+00,+1.0e+00,0.0,{values}")
    (postpro / "surface-response-model-energy.csv").write_text("\n".join(lines) + "\n")
    (postpro / "palace.json").write_text(json.dumps({"SurfaceResponse": {"ModelCatalog": [
        {"Index": index, "Name": f"model-{index}"} for index in energies]}}))
    (postpro / "config_resolved.json").write_text(json.dumps({"Boundaries": {"Postprocessing": {"Dielectric": [
        {"Index": k, "Type": cls} for k, cls in classes.items()]}}}))


def reused_entry(name, bounds, sharp):
    return {"Name": name, "QualificationStatus": "ReusedResponse", "ReuseMode": "Fallback", "ReusedFrom": {"Donor": "donor"},
            "NearKey": {"W": -0.01}, "PredictedReuseError": {**{T: {"Bound": b} for T, b in zip(report.TYPES, bounds)}, "MA_sharp": {"Bound": sharp}}}


class Synthetic(unittest.TestCase):
    def test_share_and_budget(self):
        with tempfile.TemporaryDirectory() as tmp:
            classes = {1: "SA", 2: "MS", 3: "MA", 4: "MS", 5: "MA"}
            energies = {1: {1: 2.0, 2: 1.0, 3: 0.1, 4: 1.0, 5: 0.1},     # reused: SA 2, MS 2, MA 0.2
                        2: {1: 6.0, 2: 3.0, 3: 0.3, 4: 3.0, 5: 0.3},     # exact: SA 6, MS 6, MA 0.6
                        3: {1: 2.0, 2: 2.0, 3: 0.2, 4: 2.0, 5: 0.2}}     # reused: SA 2, MS 4, MA 0.4
            write_run(Path(tmp) / "postpro", energies, classes)
            library = {"Name": "lib", "NearKeyReuse": {"RuleVersion": "nearkey-reuse-rule-v1"},
                       "Models": [reused_entry("model-1", (0.005, 0.0025, 0.01), 0.0105), {"Name": "model-2", "QualificationStatus": "Qualified"},
                                  reused_entry("model-3", (0.01, 0.01, 0.02), 0.021)]}
            out = report.window_report(library, Path(tmp) / "postpro", {"SA": 20.0, "MS": 20.0, "MA": 2.0})
            self.assertEqual(out["Counts"], {"Reused": 2, "Models": 3})
            self.assertAlmostEqual(out["ReusedShare"]["SA"], 0.2)
            self.assertAlmostEqual(out["ReusedShare"]["MS"], 0.3)
            self.assertAlmostEqual(out["ReusedShare"]["MA"], 0.3)
            # budget SA = 100 x (0.005 x 0.1 + 0.01 x 0.1) = 0.15 points; MS = 100 x (0.0025 x 0.1 + 0.01 x 0.2) = 0.225; MA raw = 100 x (0.01 x 0.1 + 0.02 x 0.2) = 0.5
            self.assertAlmostEqual(out["ReuseBudgetPoints"]["SA"], 0.15)
            self.assertAlmostEqual(out["ReuseBudgetPoints"]["MS"], 0.225)
            self.assertAlmostEqual(out["ReuseBudgetPoints"]["MA"], 0.5)
            self.assertAlmostEqual(out["ReuseBudgetPoints"]["MA_sharp"], 100 * (0.0105 * 0.1 + 0.021 * 0.2))
            self.assertEqual(out["BudgetOfRecord"]["MA_sharp"], out["ReuseBudgetPoints"]["MA_sharp"])
            self.assertEqual(out["ReuseBudgetLarge"], [])
            self.assertTrue(out["ShareConvention"].startswith("A:"))
            large = report.window_report(library, Path(tmp) / "postpro", {"SA": 1.0, "MS": 20.0, "MA": 2.0})
            self.assertEqual(large["ReuseBudgetLarge"], ["SA"])
            summed = report.window_report(library, Path(tmp) / "postpro")
            self.assertAlmostEqual(summed["ReusedShare"]["SA"], 0.4)
            self.assertTrue(summed["ShareConvention"].startswith("summed-model"))
            text = report.markdown(out)
            self.assertIn("model-1", text)
            self.assertIn("**window**", text)
            with self.assertRaises(report.WindowReportError):
                report.window_report(library, Path(tmp) / "postpro", {"SA": 0.0, "MS": 1.0, "MA": 1.0})
            no_prediction = {"Name": "lib", "Models": [{"Name": "model-1", "QualificationStatus": "ReusedResponse"}]}
            with self.assertRaises(report.WindowReportError):
                report.window_report(no_prediction, Path(tmp) / "postpro")


@unittest.skipUnless(DEMO_LIBRARY.is_file() and DERIVED.is_file() and RUNS.is_dir(), "the local evidence mirror is not mounted")
class MeasuredWindows(unittest.TestCase):
    def test_design_budgets_on_the_exact_runs(self):
        library = json.loads(DEMO_LIBRARY.read_text())
        derived = json.loads(DERIVED.read_text())
        for window, expected in DESIGN_BUDGETS.items():
            postpro = RUNS / f"{window}-thin-P4-exact" / "postpro"
            if not postpro.is_dir():
                self.skipTest(f"{postpro} not mirrored")
            out = report.window_report(library, postpro, derived["windows"][window]["reference"]["E_surf"])
            self.assertEqual(out["Counts"]["Reused"], 2, window)
            for T, value in zip(report.TYPES, expected):
                self.assertAlmostEqual(out["ReuseBudgetPoints"][T], value, places=3, msg=f"{window} {T}")
            self.assertGreater(out["ReuseBudgetPoints"]["MA_sharp"], out["ReuseBudgetPoints"]["MA"])
            self.assertLess(out["ReuseBudgetPoints"]["MA_sharp"], 1.03 * out["ReuseBudgetPoints"]["MA"] + 0.01)
            self.assertEqual(out["ReuseBudgetLarge"], [])


if __name__ == "__main__":
    unittest.main()
