# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
import copy
import csv
import hashlib
import json
from pathlib import Path
from types import SimpleNamespace
import tempfile
import unittest

from run_general_mesh_suite import audit_manifest_evidence, run_manifest


GATES = {
    "CornerTolerance": 1e-8,
    "MaximumNormalFactor": 2,
    "MinimumAchievedAspect": 1.5,
    "RotationTolerance": 1e-8,
    "MaximumSeconds": 10,
    "MaximumRSSGiB": 1,
    "MaximumElements": 1000,
}


def passing_evidence():
    return {
        "ExactLabelsMaterials": {
            "ExpectedVolumeAttributes": [1, 2], "ActualVolumeAttributes": [2, 1],
            "ExpectedBoundaryAttributes": [1, 3100, 5001, 6001],
            "ActualBoundaryAttributes": [6001, 5001, 3100, 1],
        },
        "OwnershipClosure": {
            "UnmatchedPolicy": "Error", "Unmatched": 0, "Overlaps": 0,
            "Exhaustive": True,
        },
        "SemanticCorners": {"Expected": [[0, 0, 0]], "Actual": [[0, 0, 0]]},
        "ProtectedSurfaces": {
            "Expected": ["matching", "physical"],
            "Actual": ["physical", "matching"], "Changed": 0,
        },
        "AchievedAnisotropy": {
            "Samples": 4, "NormalTarget": .01, "Transverse1P90": .012,
            "Transverse2P90": .013, "TangentialP50": .04,
        },
        "RotationCovariance": {
            "ComparedVariant": "rotate-z-0.63", "MaximumRelativeInvariantError": 1e-10,
        },
        "TraceDiagonal": {"GlobalDiagonalBands": 0},
        "Resources": {"ExitCode": 0, "Seconds": 2, "PeakRSSGiB": .2, "Elements": 100},
    }


class GeneralMeshManifestTest(unittest.TestCase):
    def make_case(self, root, name, rows):
        directory = root / name
        directory.mkdir()
        signature = directory / "signature.csv"
        with signature.open("w", newline="") as stream:
            writer = csv.writer(stream)
            writer.writerow(["Index", "Slot", "Conductor"])
            writer.writerows((i + 1, i % 2, i % 3 + 1) for i in range(rows))
        digest = hashlib.sha256(signature.read_bytes()).hexdigest()
        return {
            "Id": name, "Variants": ["identity", "rotated"],
            "Source": {
                "Directory": name, "SignatureRole": "Signature",
                "SignatureColumns": ["Index", "Slot", "Conductor"],
                "Files": {"Signature": {"Name": "signature.csv", "SHA256": digest}},
            },
        }

    def args(self, manifest, output, *, audit=None):
        return SimpleNamespace(manifest=manifest, input=[], root=output,
                               preflight_only=audit is None, audit_root=audit)

    def test_discovers_counts_without_fixed_edge_or_source_assumptions(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            cases = [self.make_case(root, "short", 1), self.make_case(root, "long", 10)]
            manifest = root / "suite.json"
            manifest.write_text(json.dumps({"Version": 1, "RepositoryRoot": ".",
                                            "Gates": GATES, "Cases": cases}))
            output = root / "out"
            self.assertTrue(run_manifest(self.args(manifest, output)))
            summary = json.loads((output / "summary.json").read_text())
            self.assertEqual([case["DiscoveredEdgeCount"] for case in summary["Cases"]], [1, 10])
            audits = root / "audits"
            audits.mkdir()
            for case in cases:
                (audits / f"{case['Id']}.json").write_text(json.dumps(passing_evidence()))
            audited_output = root / "audited"
            self.assertTrue(run_manifest(self.args(manifest, audited_output, audit=audits)))
            audited = json.loads((audited_output / "summary.json").read_text())
            self.assertTrue(audited["Passed"])
            self.assertTrue(all(not case["GateFailures"] for case in audited["Cases"]))

    def test_required_six_and_ten_inputs_fail_closed_before_audits(self):
        for missing in ("six", "ten"):
            with self.subTest(missing=missing), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)
                cases = [self.make_case(root, name, count)
                         for name, count in (("six", 6), ("ten", 10))]
                for case in cases:
                    if case["Id"] == missing:
                        case["Source"]["Directory"] = None
                manifest = root / "suite.json"
                manifest.write_text(json.dumps({"Version": 1, "RepositoryRoot": ".",
                                                "Gates": GATES, "Cases": cases}))
                output = root / "out"
                args = self.args(manifest, output, audit=root / "absent-audits")
                self.assertFalse(run_manifest(args))
                summary = json.loads((output / "summary.json").read_text())
                self.assertFalse(summary["PreflightPassed"])
                failed = next(case for case in summary["Cases"] if case["Id"] == missing)
                self.assertIn("unavailable", failed["Error"])
                self.assertNotIn("AuditEvidence", failed)

    def test_every_mesh_gate_is_fail_closed(self):
        self.assertEqual(audit_manifest_evidence(passing_evidence(), GATES), [])
        self.assertEqual(set(audit_manifest_evidence({}, GATES)), {
            "exact-labels-materials", "ownership-exhaustive-closure", "semantic-corners",
            "protected-surfaces", "achieved-anisotropy", "rotation-covariance",
            "trace-diagonal-overrefinement", "bounded-resources",
        })
        mutations = {
            "exact-labels-materials": ("ExactLabelsMaterials", "ActualVolumeAttributes", [1]),
            "ownership-exhaustive-closure": ("OwnershipClosure", "Unmatched", 1),
            "semantic-corners": ("SemanticCorners", "Actual", [[1, 0, 0]]),
            "protected-surfaces": ("ProtectedSurfaces", "Changed", 1),
            "achieved-anisotropy": ("AchievedAnisotropy", "Samples", 0),
            "rotation-covariance": ("RotationCovariance", "MaximumRelativeInvariantError", 1),
            "trace-diagonal-overrefinement": ("TraceDiagonal", "GlobalDiagonalBands", 1),
            "bounded-resources": ("Resources", "Seconds", 11),
        }
        for gate, (section, key, value) in mutations.items():
            with self.subTest(gate=gate):
                evidence = copy.deepcopy(passing_evidence())
                evidence[section][key] = value
                self.assertIn(gate, audit_manifest_evidence(evidence, GATES))


if __name__ == "__main__":
    unittest.main()
