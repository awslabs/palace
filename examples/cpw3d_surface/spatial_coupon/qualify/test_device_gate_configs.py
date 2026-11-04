#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests of device_gate_configs.py (fix round R-4; decisions 1174 / 1238 / 1239) and of the gate reader on the Maxwell-form
record: (1) the ES default reproduces the (F) chains' write_c0_configs v2 recipe BYTE-IDENTICALLY (the recipe body is kept
here verbatim as the reference); (2) the driven form's leaves; (3) the refusals (form / type mismatch, no excited port, no
frequency, AMR on an ES base, overwrite); (4) spatial_qualification.gate_solve_record + resolution_aware_gate on a palace.json
of the Maxwell path (the shared FinishConductorConsistencyGate record form: Maxwell true / Frequency per record) read against
the same bars as the electrostatic record."""

import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import device_gate_configs as dgc  # noqa: E402
import spatial_qualification as sq  # noqa: E402


def es_base():
    """A device ES config in the registration form of record (the spatial-rebuild <W>-<KEY>-k2d.json shape)."""
    return {
        "Problem": {"Type": "Electrostatic", "Verbose": 1, "Output": "/cluster/registration/key/preflight/postpro",
                    "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": "/cluster/mesh/window.msh", "L0": 1e-6,
                  "Refinement": {"MaxIts": 0, "UniformLevels": 0, "SerialUniformLevels": 0}},
        "Domains": {"Materials": [{"Attributes": [1], "Permittivity": 1.0}]},
        "Boundaries": {"Ground": {"Attributes": [2]}, "Terminal": [{"Index": 1, "Attributes": [3]}],
                       "Postprocessing": {"Dielectric": [{"Index": 1, "Type": "SA", "Attributes": [4]}]}},
        "Solver": {"Order": 2,
                   "Electrostatic": {"Save": 0,
                                     "ResponseCorrection": {"Library": "/cluster/library/base/process-library-abs.json",
                                                            "TargetInterfaces": [1], "UnmatchedPolicy": "Warn",
                                                            "CorrectionMode": "Both", "TraceCoupling": "Collocated"}},
                   "Linear": {"Type": "BoomerAMG", "KSPType": "CG", "Tol": 1e-10, "MaxIts": 1000}},
    }


def driven_base():
    """A device Driven config in the F-39 port form (the RG-7 / mxthin shape: Samples, LumpedPort, the AMR ladder)."""
    return {
        "Problem": {"Type": "Driven", "Verbose": 2, "Output": "/cluster/runs/W-corr/postpro",
                    "OutputFormats": {"Paraview": False, "GridFunction": False}},
        "Model": {"Mesh": "/cluster/mesh/window.msh", "L0": 1e-6,
                  "Refinement": {"MaxIts": 15, "UniformLevels": 0, "SerialUniformLevels": 0, "Tol": 0.01,
                                 "Nonconformal": True, "SaveAdaptIterations": True}},
        "Domains": {"Materials": [{"Attributes": [1], "Permittivity": 1.0}]},
        "Boundaries": {"PEC": {"Attributes": [2]},
                       "LumpedPort": [{"Index": 1, "R": 50.0, "Excitation": True,
                                       "Elements": [{"Attributes": [5], "Direction": [1.0, 0.0, 0.0], "Length": 3.0}]}],
                       "Postprocessing": {"Dielectric": [{"Index": 1, "Type": "SA", "Attributes": [4]}]}},
        "Solver": {"Order": 4,
                   "Driven": {"Samples": [{"Type": "Point", "Freq": [5.2, 5.7, 6.2], "SaveStep": 0}]},
                   "Linear": {"Type": "STRUMPACK", "KSPType": "GMRES", "Tol": 1e-8, "MaxIts": 500},
                   "SurfaceResponseCorrection": {"Library": "/cluster/library/E/process-library-abs.json",
                                                 "TargetInterfaces": [1], "UnmatchedPolicy": "Warn",
                                                 "CorrectionMode": "PostprocessOnly", "TraceCoupling": "Collocated",
                                                 "Diagnostics": {"ExclusionEnergy": True}, "PatchEnergy": True}},
    }


def chain_recipe_v2(base, lib, outdir):
    """The body of the (F) chains' write_c0_configs v2 (rerun-3 R3 -> spatial-rebuild write_c0_configs_sr.py), verbatim in
    its assignments and serialisation, returning {file name -> text} instead of writing."""
    c = json.loads(json.dumps(base))
    c["Solver"]["Electrostatic"]["ResponseCorrection"]["Library"] = lib
    assert c["Model"]["Refinement"]["MaxIts"] == 0
    outdir = Path(outdir)
    c["Problem"]["Output"] = f"{outdir}/preflight/postpro"
    pre = json.loads(json.dumps(c))
    pre["Solver"]["Electrostatic"]["ResponseCorrection"]["GeometryFailurePolicy"] = "RecordRaw"
    files = {"preflight.json": json.dumps(pre, indent=1) + "\n"}
    for order in (3, 4, 5):
        d = json.loads(json.dumps(c))
        d["Solver"]["Order"] = order
        d["Solver"]["Electrostatic"]["Save"] = 0
        d["Solver"]["Electrostatic"]["ResponseCorrection"]["TraceCoupling"] = "SurfaceMortar"
        d["Problem"]["Output"] = f"{outdir}/c0-p{order}/postpro"
        files[f"c0-p{order}.json"] = json.dumps(d, indent=1) + "\n"
    return files


LIB = "/cluster/library/interim/process-library-abs.json"


class EsDefaultByteIdentity(unittest.TestCase):
    def test_default_form_is_es_and_byte_identical_to_the_chain_recipe(self):
        with tempfile.TemporaryDirectory() as tmp:
            base_path = Path(tmp) / "window-key-k2d.json"
            base_path.write_text(json.dumps(es_base(), indent=1) + "\n")
            out = Path(tmp) / "registration" / "key"
            record = dgc.write_configs(base_path, form="es", library=LIB, out=str(out))
            expected = chain_recipe_v2(es_base(), LIB, str(out))
            self.assertEqual(sorted(record["Files"]), sorted(expected))
            for name, text in expected.items():
                self.assertEqual((out / name).read_text(), text, name)
            self.assertEqual(record["Form"], "es")
            self.assertEqual(record["ProblemType"], "Electrostatic")
            # A base with a missing Save leaf: the recipe appends it after the same keys (insertion order preserved).
            base = es_base()
            del base["Solver"]["Electrostatic"]["Save"]
            got = dgc.device_gate_configs(base, form="es", library=LIB, out=str(out))
            for name, text in chain_recipe_v2(base, LIB, str(out)).items():
                self.assertEqual(dgc.dump_config(got[name]), text, name)

    def test_cli_default_form_and_record(self):
        with tempfile.TemporaryDirectory() as tmp:
            base_path = Path(tmp) / "base.json"
            base_path.write_text(json.dumps(es_base()))
            out = os.path.join(tmp, "out")
            self.assertEqual(dgc.main(["--base", str(base_path), "--library", LIB, "--out", out]), 0)
            record = json.loads(Path(out, dgc.RECORD_NAME).read_text())
            self.assertEqual(record["Form"], "es")
            self.assertEqual(record["Orders"], [3, 4, 5])
            self.assertIsNone(record["FrequencyGHz"])
            self.assertTrue(record["Verification"].startswith("the form of record"))
            self.assertEqual(sorted(p.name for p in Path(out).glob("*.json")),
                             sorted(["preflight.json", "c0-p3.json", "c0-p4.json", "c0-p5.json", dgc.RECORD_NAME]))


class DrivenForm(unittest.TestCase):
    def test_driven_leaves(self):
        base = driven_base()
        got = dgc.device_gate_configs(base, form="driven", library=LIB, out="/cluster/gate/key", frequency=5.2)
        self.assertEqual(sorted(got), ["c0-p3.json", "c0-p4.json", "c0-p5.json", "preflight.json"])
        pre = got["preflight.json"]
        self.assertEqual(pre["Problem"]["Type"], "Driven")
        self.assertEqual(pre["Problem"]["Output"], "/cluster/gate/key/preflight/postpro")
        self.assertEqual(pre["Solver"]["SurfaceResponseCorrection"]["Library"], LIB)
        self.assertEqual(pre["Solver"]["SurfaceResponseCorrection"]["GeometryFailurePolicy"], "RecordRaw")
        self.assertEqual(pre["Solver"]["Driven"]["Samples"], [{"Type": "Point", "Freq": [5.2], "SaveStep": 0}])
        self.assertEqual(pre["Model"]["Refinement"]["MaxIts"], 0)
        self.assertNotIn("Electrostatic", pre["Solver"])
        for order in (3, 4, 5):
            c0 = got[f"c0-p{order}.json"]
            self.assertEqual(c0["Solver"]["Order"], order)
            self.assertEqual(c0["Problem"]["Output"], f"/cluster/gate/key/c0-p{order}/postpro")
            self.assertEqual(c0["Solver"]["SurfaceResponseCorrection"]["TraceCoupling"], "Collocated")
            self.assertNotIn("GeometryFailurePolicy", c0["Solver"]["SurfaceResponseCorrection"])
            self.assertEqual(c0["Solver"]["Driven"]["Samples"][0]["Freq"], [5.2])
            self.assertEqual(c0["Model"]["Refinement"]["MaxIts"], 0)
            # The port, the PEC, the dielectric postprocessing and the linear solver of the device config are kept.
            self.assertEqual(c0["Boundaries"]["LumpedPort"], base["Boundaries"]["LumpedPort"])
            self.assertEqual(c0["Solver"]["Linear"], base["Solver"]["Linear"])
            self.assertEqual(c0["Solver"]["SurfaceResponseCorrection"]["PatchEnergy"], True)
        # The base is not mutated.
        self.assertEqual(base, driven_base())

    def test_driven_drops_the_range_and_sweep_leaves_by_name(self):
        """Confirmation d1243 MINOR-2: a base carrying the legacy range or the sweep-level Save / AdaptiveTol / Restart leaves
        yields the one-sample Driven block (Palace would add a second range, verify Save against the one sample, or run an
        adaptive / restarted sweep over it)."""
        base = driven_base()
        base["Solver"]["Driven"].update({"MinFreq": 5.0, "MaxFreq": 6.5, "FreqStep": 0.5, "SaveStep": 1,
                                         "Save": [5.7], "AdaptiveTol": 1e-3, "Restart": 2})
        got = dgc.device_gate_configs(base, form="driven", library=LIB, out="/cluster/gate/key", frequency=5.2)
        self.assertEqual(dgc.DRIVEN_LEAVES_DROPPED,
                         ("MinFreq", "MaxFreq", "FreqStep", "SaveStep", "Save", "AdaptiveTol", "Restart"))
        for name, config in got.items():
            with self.subTest(file=name):
                self.assertEqual(config["Solver"]["Driven"],
                                 {"Samples": [{"Type": "Point", "Freq": [5.2], "SaveStep": 0}]})
        # A base without them is unchanged by the drop (the test_driven_leaves base).
        plain = dgc.device_gate_configs(driven_base(), form="driven", library=LIB, out="/cluster/gate/key", frequency=5.2)
        self.assertEqual(plain["c0-p3.json"]["Solver"]["Driven"], got["c0-p3.json"]["Solver"]["Driven"])
        self.assertEqual(sorted(base["Solver"]["Driven"]),
                         ["AdaptiveTol", "FreqStep", "MaxFreq", "MinFreq", "Restart", "Samples", "Save", "SaveStep"])

    def test_driven_record_is_unverified(self):
        with tempfile.TemporaryDirectory() as tmp:
            base_path = Path(tmp) / "driven.json"
            base_path.write_text(json.dumps(driven_base()))
            record = dgc.write_configs(base_path, form="driven", library=LIB, out=os.path.join(tmp, "g"), frequency=5.7,
                                       orders=(4,))
            self.assertEqual(record["FrequencyGHz"], 5.7)
            self.assertEqual(record["Orders"], [4])
            self.assertEqual(record["Files"], ["c0-p4.json", "preflight.json"])
            self.assertTrue(record["Verification"].startswith("UNVERIFIED"))
            self.assertEqual(record["Leaves"]["Correction"], "Solver.SurfaceResponseCorrection")


class Refusals(unittest.TestCase):
    def test_form_type_mismatch(self):
        with self.assertRaisesRegex(dgc.DeviceGateConfigError, "Problem.Type Driven"):
            dgc.device_gate_configs(es_base(), form="driven", library=LIB, out="/o", frequency=5.2)
        with self.assertRaisesRegex(dgc.DeviceGateConfigError, "Problem.Type Electrostatic"):
            dgc.device_gate_configs(driven_base(), form="es", library=LIB, out="/o")

    def test_driven_needs_a_frequency_and_an_excited_port(self):
        with self.assertRaisesRegex(dgc.DeviceGateConfigError, "--frequency"):
            dgc.device_gate_configs(driven_base(), form="driven", library=LIB, out="/o")
        base = driven_base()
        base["Boundaries"]["LumpedPort"][0]["Excitation"] = False
        with self.assertRaisesRegex(dgc.DeviceGateConfigError, "Excitation true"):
            dgc.device_gate_configs(base, form="driven", library=LIB, out="/o", frequency=5.2)

    def test_es_refuses_amr_and_an_unknown_form(self):
        base = es_base()
        base["Model"]["Refinement"]["MaxIts"] = 3
        with self.assertRaisesRegex(dgc.DeviceGateConfigError, "MaxIts 0"):
            dgc.device_gate_configs(base, form="es", library=LIB, out="/o")
        with self.assertRaisesRegex(dgc.DeviceGateConfigError, "unknown form"):
            dgc.device_gate_configs(es_base(), form="eigen", library=LIB, out="/o")

    def test_never_overwrites(self):
        with tempfile.TemporaryDirectory() as tmp:
            base_path = Path(tmp) / "base.json"
            base_path.write_text(json.dumps(es_base()))
            out = os.path.join(tmp, "out")
            dgc.write_configs(base_path, form="es", library=LIB, out=out)
            with self.assertRaisesRegex(dgc.DeviceGateConfigError, "refusing to overwrite"):
                dgc.write_configs(base_path, form="es", library=LIB, out=out)
            self.assertEqual(dgc.main(["--base", str(base_path), "--library", LIB, "--out", out]), 2)


def maxwell_record(model, max_ratio, frequency, excluded=False):
    """One ConductorConsistency record as FinishConductorConsistencyGate writes it on the Maxwell path (the leaves the
    reader uses plus the path markers)."""
    return {"Patch": 0, "Model": model, "ModelIndex": 1, "Source": 1, "MaxRatio": max_ratio, "Excluded": excluded,
            "Maxwell": True, "Frequency": frequency, "FloorApplied": False, "Amplitude": 1.0, "State": 1.0,
            "ExcitationPotential": 1.0, "Normalization": 1.0, "MaxDeviation": max_ratio,
            "ProbeKind": "TranslationalEdge", "SelfConsistentPreGate": False, "RecordOnly": False}


class MaxwellFormGateRecord(unittest.TestCase):
    """The Driven c0's palace.json carries the SAME SurfaceResponse.Diagnostics.ConductorConsistency form as the electrostatic
    one (FinishConductorConsistencyGate is shared by both ApplyConductorConsistencyGate overloads), with Maxwell true and the
    Frequency per record: gate_solve_record / resolution_aware_gate read it against the same 1e-4 bar."""

    MODEL = "parallel-edge-cluster-95b13e7fdc7c-pieces3"

    def palace_json(self, tmp, name, records, count, order, problem_type="Driven"):
        run = Path(tmp) / name
        run.mkdir(parents=True)
        (run / "palace.json").write_text(json.dumps({"Problem": {"Type": problem_type}, "SurfaceResponse": {"Diagnostics": {
            "ConductorConsistency": {"Count": count, "Tolerance": 0.02, "AmplitudeFloor": 0.001, "TestedPatches": len(records),
                                     "Records": records}}}}))
        (run / "config_resolved.json").write_text(json.dumps({"Problem": {"Type": problem_type}, "Solver": {"Order": order}}))
        return run / "palace.json"

    def test_passing_maxwell_records_at_three_frequencies(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [maxwell_record(self.MODEL, r, f) for f, r in ((5.2, 2.2e-7), (5.7, 6.3e-7), (6.2, 4.1e-7))]
            records += [maxwell_record("same-conductor-gap-3um", 1.1e-4, 5.2)]  # another model's ratio does not count
            solves = [sq.gate_solve_record(self.palace_json(tmp, f"c0-p{order}", records, 0, order)) for order in (4, 5)]
            self.assertEqual([s["Order"] for s in solves], [4, 5])
            gate = sq.resolution_aware_gate(solves, self.MODEL)
            self.assertTrue(gate["Passed"])
            self.assertEqual(gate["Orders"], [4, 5])
            self.assertEqual(gate["ProductionOrdersMissing"], [])
            self.assertAlmostEqual(gate["MaxRatio"], 6.3e-7)
            self.assertEqual((gate["Probed"], gate["Count"], gate["Excluded"]), (6, 0, 0))
            self.assertEqual(gate["Tolerance"], sq.GATE_MAX_RATIO)

    def test_a_maxwell_ratio_over_the_bar_at_one_frequency_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            records = [maxwell_record(self.MODEL, 2.2e-7, 5.2), maxwell_record(self.MODEL, 3.0e-4, 5.7)]
            gate = sq.resolution_aware_gate([sq.gate_solve_record(self.palace_json(tmp, "c0-p4", records, 0, 4))], self.MODEL)
            self.assertFalse(gate["Passed"])
            self.assertAlmostEqual(gate["MaxRatio"], 3.0e-4)
            # The solve-time exclusion (MaxRatio > 0.02, the C++ Tolerance) fails too, as on the electrostatic path.
            records = [maxwell_record(self.MODEL, 0.03, 5.2, excluded=True)]
            gate = sq.resolution_aware_gate([sq.gate_solve_record(self.palace_json(tmp, "c0-p5", records, 1, 5))], self.MODEL)
            self.assertFalse(gate["Passed"])
            self.assertEqual((gate["Count"], gate["Excluded"]), (1, 1))

    def test_same_bars_on_both_paths(self):
        """The electrostatic record (no Maxwell / Frequency leaves) and the Maxwell record of the same ratios grade alike."""
        with tempfile.TemporaryDirectory() as tmp:
            es_records = [{k: v for k, v in maxwell_record(self.MODEL, r, 0.0).items() if k not in ("Maxwell", "Frequency")}
                          for r in (2.2e-7, 6.3e-7)]
            mx_records = [maxwell_record(self.MODEL, r, f) for r, f in ((2.2e-7, 5.2), (6.3e-7, 5.2))]
            es = sq.resolution_aware_gate([sq.gate_solve_record(self.palace_json(tmp, "es", es_records, 0, 4, "Electrostatic"))],
                                          self.MODEL)
            mx = sq.resolution_aware_gate([sq.gate_solve_record(self.palace_json(tmp, "mx", mx_records, 0, 4))], self.MODEL)
            for key in ("Passed", "MaxRatio", "Probed", "Count", "Excluded", "Orders", "Tolerance"):
                self.assertEqual(es[key], mx[key], key)


if __name__ == "__main__":
    unittest.main()
