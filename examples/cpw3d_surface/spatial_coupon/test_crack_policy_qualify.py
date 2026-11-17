#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""P3-39 + ruling 1597 (upstream integration P3): the crack policy in the spatial qualification.
(2) --reuse-main: a stored main stage is reused only under the same crack policy; a stored record without
Inputs.CrackPolicy reads as legacy, so a legacy run still reuses a pre-integration stored main and a default run
stops. (3a) --reference: a reference config without the key reads as legacy; under the default policy the
cross-form pair is RECORDED (ConfigEqualsDerived.CrackPolicy) and (a) LABELLED in the qualification outcome: the
model's Qualification in the written library, the library totals and the LIBRARY line. The device adapter passes
the policy to the generator and records it."""

import json
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "qualify"))
sys.path.insert(0, str(HERE))
import case_inputs  # noqa: E402
import device_coupons  # noqa: E402
import qualify_library  # noqa: E402
import test_write_process_library as writer_fixtures  # noqa: E402  (module import: its TestCase is not collected here)
from crack_policy import CRACK_KEY  # noqa: E402

LABEL = "reference crack policy: legacy (assumed: no record) vs run: default"


def sha256(path):
    return qualify_library.sha256(path)


class ReusedMainAcrossCrackPolicies(unittest.TestCase):
    """A stored qualify root as test_qualify_dry_run.ReusedMainTest builds it, with and without Inputs.CrackPolicy."""

    CASE = "spatial-8-edge-05c322f6cda8"
    PREFIX = "spatial-8-edge-05c322f6cda8-p4"

    def store(self, stored_config, stored_policy=None):
        root = Path(tempfile.mkdtemp(prefix="p3f-reuse-")) / f"fab-{self.CASE}"
        self.addCleanup(shutil.rmtree, root.parent, True)
        case_root = root / self.CASE
        (case_root / "inputs").mkdir(parents=True)
        (case_root / "inputs" / "run-config.json").write_text(json.dumps(stored_config, indent=2) + "\n")
        (case_root / "main" / "jobs" / "worker-1").mkdir(parents=True)
        (case_root / "main" / "jobs" / "worker-1" / "plan.json").write_text(json.dumps(
            {"PinnedSHA256": {**{f"/old/root/inputs/traces/{name}": digest for name, digest in self.traces.items()},
                              "/old/root/mesh/identity-cdb48405d063.msh": "c" * 64}}))
        reducer = case_root / "results" / "main" / self.PREFIX / "reducer"
        reducer.mkdir(parents=True)
        (reducer / "domain-response-matrix.csv").write_text("basis_i,basis_j,Q_ij (J)\n1,1,1.5\n")
        (reducer / "surface-response-matrix.csv").write_text("interface,edge,R (m),basis_i,basis_j,Q_ij (J)\n2,0,2e-6,1,1,0.5\n")
        digests = {name: sha256(reducer / name) for name in qualify_library.REUSED_MAIN_CSVS}
        inputs = {"ConfigSHA256": sha256(case_root / "inputs" / "run-config.json")}
        if stored_policy is not None:
            inputs["CrackPolicy"] = stored_policy
        stored = {"Case": self.CASE, "Status": "failed", "StoppedBy": None,
                  "Mesh": {"Local": "/old/identity.msh", "SHA256": "c" * 64, "Verified": True}, "Inputs": inputs,
                  "Stages": [{"Prefix": self.PREFIX, "Role": "main", "Order": 4}],
                  "Controls": {"Indices": [1, 2, 3]},
                  "Jobs": [{"Name": "reducer", "Kind": "reducer", "Nodes": 1, "Ranks": 192, "Submission": {"Job": "57843.h"}}],
                  "Nodes": {"MultiNode": False, "Main": 1, "Fixed": 1},
                  "ResultDigests": {f"main/{self.PREFIX}/reducer/{name}": {"Local": digest, "Remote": digest, "OK": True}
                                    for name, digest in digests.items()},
                  "Qualification": {"Verdict": "Failed", "UnjudgedTypes": []},
                  "Cost": {"MainStage": {"NodeHours": 1.2}, "JobNodeHours": 2.3}}
        (root / "library-qualification.json").write_text(json.dumps({"Version": 1, "ToolCommit": "7b95205da0",
                                                                     "Gates": {"SHA256": "g" * 64}, "Cases": [stored]}))
        return root

    def setUp(self):
        self.traces = {f"basis-{k:04d}.csv": f"{k:064x}" for k in (1, 2, 3)}
        self.pre_integration = {"Problem": {"Type": "Electrostatic", "Output": "/old/root/main"},
                                "Model": {"Mesh": "/old/root/mesh/identity-cdb48405d063.msh"},
                                "Solver": {"Order": 4, "Linear": {"Tol": 1e-10}},
                                "Boundaries": {"PrescribedPotential": [
                                    {"Index": k, "DataFile": f"/old/root/inputs/traces/basis-{k:04d}.csv"} for k in (1, 2, 3)]}}

    def derived(self, policy):
        """The run config derived now: the producer's Model under the policy (legacy writes the key true)."""
        config = json.loads(json.dumps(self.pre_integration))
        config["Model"]["Mesh"] = "/new/root/mesh/identity-cdb48405d063.msh"
        config["Problem"]["Output"] = "/new/root/main"
        if policy == "legacy":
            config["Model"][CRACK_KEY] = True
        for entry in config["Boundaries"]["PrescribedPotential"]:
            entry["DataFile"] = entry["DataFile"].replace("/old/root", "/new/root")
        return config

    def load(self, root, policy, **kwargs):
        return qualify_library.load_reused_main(root, self.CASE, identity_sha256="c" * 64, run_config=self.derived(policy),
                                                trace_digests=dict(self.traces), main_prefix=self.PREFIX, **kwargs)

    def stop(self, root, policy, *messages):
        with self.assertRaises(qualify_library.CaseStop) as context:
            self.load(root, policy, crack_policy=policy)
        self.assertEqual(context.exception.record["Kind"], "ReuseMain")
        for message in messages:
            self.assertIn(message, context.exception.record["Message"])

    def test_a_legacy_run_still_reuses_a_pre_integration_stored_main(self):
        root = self.store(self.pre_integration)   # no Inputs.CrackPolicy, no key in the stored config
        self.assertEqual(qualify_library.stored_crack_policy(json.loads(
            (root / "library-qualification.json").read_text())["Cases"][0]), "legacy")
        reused = self.load(root, "legacy", crack_policy="legacy")
        self.assertTrue(reused["Mesh"]["Identical"] and reused["Reducer"]["VerifiedAgainstStoredRecord"])
        self.assertEqual(self.load(root, "legacy")["Case"], self.CASE)   # the keyword's default = the shared default

    def test_a_default_run_never_reuses_a_cracked_stored_main(self):
        root = self.store(self.pre_integration)
        # config_differences alone would accept it (both configs lack the key): the policy check stops it
        self.assertEqual(case_inputs.config_differences(self.derived("default"), self.pre_integration), [])
        self.stop(root, "default", "crack policy legacy (assumed: no record)", "this run's is default",
                  "never reused across crack policies")
        recorded = self.store(self.derived("legacy") | {"Problem": self.pre_integration["Problem"]}, "legacy")
        self.stop(recorded, "default", "crack policy legacy,", "this run's is default")

    def test_stored_default_main_is_reused_by_default_runs_only(self):
        root = self.store(self.pre_integration, "default")
        self.assertEqual(self.load(root, "default", crack_policy="default")["Case"], self.CASE)
        self.stop(root, "legacy", "crack policy default", "this run's is legacy")

    def test_a_stored_record_with_an_unknown_policy_fails_closed(self):
        with self.assertRaises(ValueError):
            qualify_library.stored_crack_policy({"Inputs": {"CrackPolicy": "cracked"}})


class ReferenceCrossFormLabel(unittest.TestCase):
    def test_the_pair_is_recorded_and_labelled_only_when_cross_form(self):
        absent, explicit = {"Model": {"Mesh": "m"}}, {"Model": {"Mesh": "m", CRACK_KEY: True}}
        same = qualify_library.reference_crack_policy("legacy", absent)
        self.assertEqual((same["Run"], same["Reference"], same["CrossForm"], same["Label"]), ("legacy", "legacy", False, None))
        self.assertEqual(same["ReferenceBasis"], "assumed: no record")
        cross = qualify_library.reference_crack_policy("default", absent)
        self.assertEqual((cross["CrossForm"], cross["Label"]), (True, LABEL))
        self.assertEqual(cross["Rule"], qualify_library.CROSS_FORM_REFERENCE_RULE)
        self.assertEqual(qualify_library.reference_crack_policy("default", explicit)["Label"],
                         f"reference crack policy: legacy (Model.{CRACK_KEY} true) vs run: default")
        self.assertIsNone(qualify_library.reference_crack_policy("legacy", explicit)["Label"])
        with self.assertRaises(qualify_library.CaseStop):
            qualify_library.reference_crack_policy("default", {"Model": {CRACK_KEY: False}})
        record = {"Reference": {"ConfigEqualsDerived": {"CrackPolicy": cross}}}
        self.assertEqual(qualify_library.cross_form_reference_label(record), LABEL)
        self.assertIsNone(qualify_library.cross_form_reference_label({"Reference": {"Rule": "--reference none"}}))

    def test_the_label_reaches_the_library_and_its_summary_line(self):
        tmp = Path(tempfile.mkdtemp(prefix="p3f-label-"))
        self.addCleanup(shutil.rmtree, tmp, True)
        root = tmp / "root"
        records, contexts, entries = [], {}, []
        for case_id, name, label in (("case-a", "model_a", LABEL), ("case-b", "model_b", None)):
            source = writer_fixtures.write_source(tmp / "sources" / case_id, name)
            reducer = root / case_id / "results" / "main" / f"{case_id}-p4" / "reducer"
            writer_fixtures.write_reducer(reducer)
            qualification = {"Verdict": "PendingQualification", "Path": str(root / case_id / "qualification.json"),
                             "ReferenceAnchor": None}
            if label:
                qualification["CrossFormReference"] = label
            records.append({"Case": case_id, "Root": str(root / case_id), "MainReducer": str(reducer), "Status": "pending",
                            "Mesh": {"Local": "/nowhere/identity.msh", "SHA256": "0" * 64}, "Qualification": qualification})
            contexts[case_id] = {"layout": [{"Role": "main", "Prefix": f"{case_id}-p4", "Order": 4}], "ma_tail": None}
            entries.append({"Id": case_id, "Source": {"Directory": str(source.parent),
                                                      "Files": {"ProcessLibrary": {"Name": "process-library.json"}}}})
        library = qualify_library.process_library_entries(records, contexts, manifest_path=tmp / "manifest.json",
                                                          manifest={"RepositoryRoot": ".", "Cases": entries}, root=root)
        models = {model["Name"]: model for model in library["Models"]}
        self.assertEqual(models["model_a"]["Qualification"]["CrossFormReference"], LABEL)
        self.assertNotIn("CrossFormReference", models["model_b"]["Qualification"])
        counts = qualify_library.cross_form_reference_counts(records)
        self.assertEqual(counts, {LABEL: 1})
        totals = {"CouponsAttempted": 2, "CouponsQualified": 0, "CouponsPending": 2, "CouponsFailed": 0, "CouponsPlanned": 0,
                  "CouponsSkipped": 0, "CouponsUnsupported": 0, "NodeHours": 1.0, "JobsSubmitted": 2, "MaxJobs": 40,
                  "UserJobCap": 100}
        line = qualify_library.library_summary_line({"Root": "/r", "Library": dict(totals, CouponsCrossFormReference=counts)})
        self.assertIn(f"pending 2 failed 0 planned 0 skipped 0 unsupported 0 cross-form-reference 1 ({LABEL}) node-h 1.000", line)
        # without a cross-form coupon the line is the pre-P3 line, byte for byte
        plain = qualify_library.library_summary_line({"Root": "/r", "Library": dict(totals, CouponsCrossFormReference={})})
        self.assertEqual(plain, "LIBRARY attempted 2 qualified 0 pending 2 failed 0 planned 0 skipped 0 unsupported 0 "
                                "node-h 1.000 jobs 2/40 (cap 100); record /r/library-qualification.json")


class DeviceAdapterPassesThePolicy(unittest.TestCase):
    def test_generator_command_carries_the_policy(self):
        coupon = {"Id": "spatial-x", "Geometry": {}}
        parameters = {"metal_thickness": 0.1, "overetch": 0.05, "sidewall_angle": 80.0, "top_radius": 0.01,
                      "bottom_radius": 0.01}
        with tempfile.TemporaryDirectory() as directory, \
                mock.patch.object(device_coupons, "coupon_geometry", return_value={"Geometry": {}}), \
                mock.patch.object(device_coupons.planner, "material_options", return_value=[]), \
                mock.patch.object(device_coupons, "stamp_signature_model"), \
                mock.patch.object(device_coupons.subprocess, "run", return_value=subprocess.CompletedProcess([], 0)) as run:
            for policy in ("legacy", "default"):
                command = device_coupons.generate_sources(coupon, Path(directory) / policy, radius=2.0,
                                                          parameters=parameters, ring_size=16, crack_policy=policy)
                self.assertEqual(command[-2:], ["--crack-policy", policy])
                self.assertEqual(run.call_args[0][0], command)
                self.assertIn("--basis-only", command)   # no Palace config is written by this run
            default = device_coupons.generate_sources(coupon, Path(directory) / "unset", radius=2.0, parameters=parameters,
                                                      ring_size=16)
            self.assertEqual(default[-1], device_coupons.CRACK_POLICY_DEFAULT)
            with self.assertRaises(ValueError):
                device_coupons.generate_sources(coupon, Path(directory) / "x", radius=2.0, parameters=parameters,
                                                ring_size=16, crack_policy="none")


if __name__ == "__main__":
    unittest.main()
