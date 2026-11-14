#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""P3-39 (P3-F review m-1): prepare_case passes the run's crack policy on, data-free.
The policy of --crack-policy reaches case_inputs.derive, and the policy derive records
(Inputs.CrackPolicy) reaches load_reused_main (--controls-only --reuse-main) and
reference_crack_policy (--reference). The case's sources are mocked: no build, mesh, or run
data is read; each leg stops at the mocked call it checks."""

import hashlib
import json
import sys
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / "qualify"))
sys.path.insert(0, str(HERE))
import qualify_library  # noqa: E402


class Reached(Exception):
    """Raised by the mocked call under test: prepare_case stops there."""


class PrepareCaseCrackPolicy(unittest.TestCase):
    CASE = "coupon-a"

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        mesh = self.root / "identity.msh"
        mesh.write_text("identity mesh\n")
        self.case_record = {"Case": self.CASE, "Status": "OK", "Passed": True, "Root": None,
                            "Variants": {"identity": {
                                "Path": str(mesh),
                                "SHA256": hashlib.sha256(mesh.read_bytes()).hexdigest()}}}
        self.reference_config = {"Problem": {"Type": "Electrostatic"}}
        reference_path = self.root / "reference-config.json"
        reference_path.write_text(json.dumps(self.reference_config))
        self.reference = {"Config": str(reference_path), "ReferenceOrder": 3, "LinearTol": 1e-10}
        case = {"Source": {"Files": {"BasisContract": {"SHA256": "b" * 64},
                                     "Signature": {"Name": "signature.csv"}}}}
        self.derive = mock.Mock(side_effect=self.derived)
        patches = [
            mock.patch.object(qualify_library, "manifest_case", return_value=case),
            mock.patch.object(qualify_library, "case_kind", return_value="fabricated"),
            mock.patch.object(qualify_library, "source_directory", return_value=self.root),
            mock.patch.object(qualify_library.case_inputs, "physics_run_parameters",
                              return_value={"Order": 2, "LinearTol": 1e-8}),
            mock.patch.object(qualify_library.case_inputs, "derive", self.derive),
            mock.patch.object(qualify_library.case_inputs, "interface_types", return_value={}),
            mock.patch.object(qualify_library.case_inputs, "config_differences", return_value=[]),
            mock.patch.object(qualify_library.locate_sources, "read_signature_rows",
                              return_value=[]),
            mock.patch.object(qualify_library.locate_sources, "signature_layers",
                              return_value=[]),
            mock.patch.object(qualify_library.locate_sources, "check_layers", return_value=None),
            mock.patch.object(qualify_library.reference_campaign, "resolve",
                              return_value=self.reference),
        ]
        for patch in patches:
            patch.start()
            self.addCleanup(patch.stop)
        self.addCleanup(self.tmp.cleanup)

    @staticmethod
    def derived(case, directory, *, crack_policy, **_):
        """case_inputs.derive's pass-through: Inputs.CrackPolicy = the policy it was given."""
        return ({"Problem": {"Type": "Electrostatic"}},
                {"CrackPolicy": crack_policy, "Interfaces": {}, "InterfaceTypes": {},
                 "Sources": []})

    def prepare(self, policy, **options):
        args = SimpleNamespace(crack_policy=policy, controls_only=False, reuse_main=None,
                               stage_prefix=None, control_amplitudes=None, reference=None,
                               build_gate_override=None)
        for name, value in options.items():
            setattr(args, name, value)
        with self.assertRaises(Reached):
            qualify_library.prepare_case(self.case_record, manifest_path=self.root / "m.json",
                                         manifest={}, args=args, root=self.root / "run",
                                         remote=None, profile=None, cost_model=None,
                                         gates=None, gates_digest=None)
        self.assertEqual(self.derive.call_args.kwargs["crack_policy"], policy)

    def test_reused_main_gets_the_run_policy(self):
        for policy in ("legacy", "default"):
            with self.subTest(policy=policy), mock.patch.object(
                    qualify_library, "load_reused_main", side_effect=Reached) as reused:
                self.prepare(policy, controls_only=True, reuse_main=str(self.root / "stored"))
                self.assertEqual(reused.call_args.kwargs["crack_policy"], policy)

    def test_reference_label_gets_the_run_policy(self):
        for policy in ("legacy", "default"):
            with self.subTest(policy=policy), mock.patch.object(
                    qualify_library, "reference_crack_policy", side_effect=Reached) as label:
                self.prepare(policy, reference=str(self.root / "campaign"))
                self.assertEqual(label.call_args.args, (policy, self.reference_config))


if __name__ == "__main__":
    unittest.main()
