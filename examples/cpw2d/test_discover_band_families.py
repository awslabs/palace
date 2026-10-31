#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""The discovery's near-match band families (identification fix round DESIGN sections 4 / 5 / A.5; decisions 1075
Q1, 1079 Q4, 1097 Q2 / Q3): within one pass a Missing pair / ParallelEdgeCluster key within the binary's near-match
band of an already-planned placeholder gets NO placeholder (ServedByPlaceholder; order-independent representatives),
the binary then reads it NearMatched / GeometricNearMatched, and the final manifest keeps that status with
ServedByRequirement (counted apart from Exact / Interpolated / Missing; the planner skips it). A fake Palace
implements the C++ bands on its pass library; a band disagreement fails closed by name."""

import json
import stat
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path

CPW2D = Path(__file__).parent
SPATIAL = CPW2D.parent / "cpw3d_surface" / "spatial_coupon"
for p in (str(CPW2D), str(SPATIAL)):
    if p not in sys.path:
        sys.path.insert(0, p)
import discover_surface_response_requirements as DISCOVERY  # noqa: E402
import identification_constants as ic  # noqa: E402
import prepare_surface_response_coupons as planner  # noqa: E402

R = 1.9
LAW = "{\"Type\":\"PEC\"}"


def strip_requirement(hash_char, separation_over_R, length, count=2):
    return {"Topology": "SameConductorStrip", "Hash": hash_char * 64, "Count": count, "Instances": 1,
            "TotalEdgeLength": length, "BoundaryCondition": {"Type": "PEC"},
            "Interfaces": [{"Slot": 0, "Target": 1, "Type": "SA"}],
            "Geometry": {"EdgeCount": 2, "Separation": separation_over_R * R},
            "Signature": {"Type": "SameConductorStrip", "SeparationOverR": separation_over_R, "Edges": [
                {"Conductor": 1, "GapSide": -1, "Interfaces": ["SA"], "Law": LAW, "OffsetOverR": 0.0},
                {"Conductor": 1, "GapSide": 1, "Interfaces": ["SA"], "Law": LAW, "OffsetOverR": separation_over_R}]}}


def cluster_requirement(hash_char, offsets, gap_sides, length, count=3):
    edges = [{"Conductor": 1, "GapSide": g, "Interfaces": ["SA"], "Law": LAW, "OffsetOverR": o} for o, g in zip(offsets, gap_sides)]
    return {"Topology": "ParallelEdgeCluster", "Hash": hash_char * 64, "Count": count, "Instances": 1,
            "TotalEdgeLength": length, "BoundaryCondition": {"Type": "PEC"},
            "Interfaces": [{"Slot": 0, "Target": 1, "Type": "SA"}], "Geometry": {"EdgeCount": len(offsets)},
            "Signature": {"Type": "ParallelEdgeCluster", "Edges": edges}}


# C (the longest strip) and A are two families 0.048 R apart; B lies 1.8e-3 R from A (inside the 5e-3 R band);
# P1 is a 4-edge cluster and P2 its MIRROR orientation with one offset moved by 4e-3 R (inside the 1e-2 R band).
REQUIREMENTS = [
    strip_requirement("a", 1.5520, 54.0),
    strip_requirement("b", 1.5538, 40.0),
    strip_requirement("c", 1.6000, 100.0),
    cluster_requirement("d", [0.0, 1.0, 2.0, 3.5], [-1, 1, -1, 1], 30.0),
    cluster_requirement("e", [0.0, 1.5, 2.504, 3.5], [-1, 1, -1, 1], 20.0),
]


def fake_device(tmp, requirements, near_match=True):
    """A device + an executable fake Palace reading the pass library: a requirement is Exact on its own placeholder,
    else (near_match) NearMatched / GeometricNearMatched to a placeholder of its family within the C++ band,
    else Missing. The manifest carries the binary's Summary.Counts form (every status by name)."""
    req_path = tmp / "fake-requirements.json"
    req_path.write_text(json.dumps(requirements))
    library = tmp / "process-library.json"
    library.write_text(json.dumps({"Version": 3, "Name": "fake-seed", "MatchingRadius": R, "Models": []}))
    config = tmp / "device.json"
    config.write_text(json.dumps({"Problem": {"Type": "Electrostatic", "Output": str(tmp / "postpro")},
                                  "Solver": {"Electrostatic": {"ResponseCorrection": {"Library": str(library), "UnmatchedPolicy": "Warn"}}}}))
    palace = tmp / "fake-palace.py"
    palace.write_text(textwrap.dedent(f"""\
        #!{sys.executable}
        import json, pathlib, sys
        sys.path.insert(0, {str(CPW2D.parent / "surface_response_identification")!r})
        import signature_library
        assert sys.argv[1] == "--surface-response-preflight"
        config = json.load(open(sys.argv[2]))
        response = config["Solver"]["Electrostatic"]["ResponseCorrection"]
        library = json.load(open(response["Library"]))
        models = {{m["Name"]: m for m in library["Models"]}}
        requirements = json.load(open({str(req_path)!r}))
        counts, lengths = {{}}, {{}}
        def near(req):
            if not {near_match!r}:
                return None
            sig = req["Signature"]
            for name, m in models.items():
                ms = m.get("Signature") or {{}}
                if ms.get("Type") != sig["Type"]:
                    continue
                if sig["Type"] == "SameConductorStrip":
                    if abs(ms["SeparationOverR"] - sig["SeparationOverR"]) <= {ic.GDS_SCATTER_OVER_RADIUS!r}:
                        return name, "NearMatched"
                else:
                    for cand in (ms, signature_library.mirror_translational(ms)):
                        if [e["GapSide"] for e in cand["Edges"]] == [e["GapSide"] for e in sig["Edges"]] and \\
                           max(abs(a["OffsetOverR"] - b["OffsetOverR"]) for a, b in zip(cand["Edges"], sig["Edges"])) <= {ic.GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS!r}:
                            return name, "GeometricNearMatched"
            return None
        for req in requirements:
            own = "__preflight_placeholder_" + req["Hash"][:16]
            if own in models:
                req["Status"] = "Exact"; req["SelectedModels"] = [{{"Name": own, "Weight": 1.0}}]
            else:
                hit = near(req)
                if hit:
                    req["Status"] = hit[1]; req["SelectedModels"] = [{{"Name": hit[0], "Weight": 1.0}}]
                else:
                    req["Status"] = "Missing"
            counts[req["Status"]] = counts.get(req["Status"], 0) + req["Count"]
            lengths[req["Status"]] = lengths.get(req["Status"], 0.0) + req["TotalEdgeLength"]
        out = pathlib.Path(config["Problem"]["Output"]); out.mkdir(parents=True, exist_ok=True)
        manifest = {{"Version": 2, "Complete": counts.get("Missing", 0) == 0, "Requirements": requirements,
                    "Library": {{"Name": library["Name"], "Path": response["Library"], "MatchingRadius": library["MatchingRadius"]}},
                    "Summary": {{"Counts": counts, "TotalEdgeLengths": lengths}}}}
        (out / "surface-response-requirements.json").write_text(json.dumps(manifest))
        """))
    palace.chmod(palace.stat().st_mode | stat.S_IXUSR)
    return config, palace


class BandFamilyDiscoveryTest(unittest.TestCase):
    def test_family_reading(self):
        a, b, c, d, e = REQUIREMENTS
        (ka, na), = DISCOVERY.band_family_orientations(a)
        (kb, nb), = DISCOVERY.band_family_orientations(b)
        self.assertEqual(ka, kb)
        self.assertAlmostEqual(DISCOVERY.band_deviation_over_R(ka, na, nb), 0.0018, places=12)
        self.assertEqual(DISCOVERY.band_over_R(ka), ic.GDS_SCATTER_OVER_RADIUS)
        reps = {}
        DISCOVERY.register_representative(d, reps)
        hit = DISCOVERY.serving_representative(e, reps)
        self.assertEqual(hit[0], d["Hash"])
        self.assertAlmostEqual(hit[1], 0.004, places=9)   # through the mirror orientation
        self.assertEqual(hit[2], ic.GEOMETRIC_NEAR_MATCH_CLAIMS_OVER_RADIUS)
        self.assertIsNone(DISCOVERY.serving_representative(c, {ka: [(a["Hash"], na)]}))
        self.assertIsNone(DISCOVERY.band_family({"Topology": "IsolatedEdge", "Hash": "f" * 64, "Signature": {"Type": "IsolatedEdge"}}))

    def test_one_coupon_per_band_family_and_the_statuses_kept(self):
        with tempfile.TemporaryDirectory(prefix="discover-band-") as tmp:
            tmp = Path(tmp)
            config, palace = fake_device(tmp, REQUIREMENTS)
            manifest = DISCOVERY.discover(config, tmp / "out", palace)
            history = json.loads((tmp / "out" / "closure-history.json").read_text())
            self.assertEqual(len(history["Passes"]), 2)
            added = sorted(p["Id"] for p in history["Passes"][0]["AddedPlaceholders"])
            self.assertEqual(added, sorted(h * 16 for h in "acd"))   # the representatives: C, A, P1 (B and P2 served)
            self.assertEqual(history["PlaceholderCount"], 3)
            served = {r["Hash"]: r for r in manifest["ServedByPlaceholder"]}
            self.assertEqual(sorted(served), ["b" * 64, "e" * 64])
            self.assertEqual(served["b" * 64]["Representative"], "a" * 64)
            self.assertEqual(served["e" * 64]["Representative"], "d" * 64)
            self.assertAlmostEqual(served["b" * 64]["DeviationOverR"], 0.0018, places=12)
            self.assertEqual(manifest["Summary"]["ServedByPlaceholder"], {"Count": 2, "TotalEdgeLength": 60.0})
            by_hash = {r["Hash"]: r for r in manifest["Requirements"]}
            self.assertEqual(by_hash["b" * 64]["Status"], "NearMatched")
            self.assertEqual(by_hash["b" * 64]["ServedByRequirement"], "a" * 64)
            self.assertNotIn("SelectedModels", by_hash["b" * 64])
            self.assertIn("served by the coupon planned", by_hash["b" * 64]["Reason"])
            self.assertEqual(by_hash["e" * 64]["Status"], "GeometricNearMatched")
            self.assertEqual(by_hash["e" * 64]["ServedByRequirement"], "d" * 64)
            for h in "acd":
                self.assertEqual(by_hash[h * 64]["Status"], "Missing")   # the restored representatives
            self.assertEqual(manifest["Summary"]["Counts"],
                             {"Exact": 0, "Interpolated": 0, "Missing": 7, "NearMatched": 2, "GeometricNearMatched": 3})
            self.assertFalse(manifest["Complete"])
            # The planner plans the Missing set only: the served keys get no coupon of their own.
            library = json.loads((tmp / "process-library.json").read_text())
            plan = planner.plan_from_manifest(tmp / "out" / "surface-response-requirements.json", manifest,
                                              tmp / "process-library.json", library, include_matched=False)
            self.assertEqual({c["Hash"][:1] for c in plan["Coupons"]}, {"a", "c", "d"})

    def test_representatives_are_order_independent(self):
        with tempfile.TemporaryDirectory(prefix="discover-perm-") as tmp:
            tmp = Path(tmp)
            results = []
            for order in (REQUIREMENTS, list(reversed(REQUIREMENTS)), [REQUIREMENTS[i] for i in (1, 4, 0, 3, 2)]):
                root = tmp / f"o{len(results)}"
                root.mkdir()
                config, palace = fake_device(root, order)
                manifest = DISCOVERY.discover(config, root / "out", palace)
                history = json.loads((root / "out" / "closure-history.json").read_text())
                results.append((sorted(p["Id"] for p in history["Passes"][0]["AddedPlaceholders"]),
                                {r["Hash"]: r["Representative"] for r in manifest["ServedByPlaceholder"]}))
            self.assertEqual(results[0], results[1])
            self.assertEqual(results[0], results[2])

    def test_band_disagreement_fails_closed(self):
        """A served key the binary still reports Missing with its representative present = this reader's band
        disagrees with the C++: RuntimeError by name, never a silent second pass."""
        with tempfile.TemporaryDirectory(prefix="discover-disagree-") as tmp:
            tmp = Path(tmp)
            config, palace = fake_device(tmp, REQUIREMENTS, near_match=False)
            with self.assertRaisesRegex(RuntimeError, "band disagreement.*bbbbbbbbbbbb"):
                DISCOVERY.discover(config, tmp / "out", palace)

    def test_restore_counts_every_status_by_name(self):
        manifest = {"Library": {"Path": "x", "Name": "n"}, "Summary": {}, "Requirements": [
            {"Hash": "1" * 64, "Status": "DomainBoundary", "Count": 2, "TotalEdgeLength": 1.0},
            {"Hash": "2" * 64, "Status": "NearMatched", "Count": 1, "TotalEdgeLength": 2.0,
             "SelectedModels": [{"Name": "same-conductor-strip-2um"}]},
            {"Hash": "3" * 64, "Status": "Exact", "Count": 3, "TotalEdgeLength": 3.0}]}
        restored = DISCOVERY.restore_source_status(manifest, {"__SourcePath": "lib.json", "Name": "lib"}, {})
        self.assertEqual(restored["Summary"]["Counts"], {"Exact": 3, "Interpolated": 0, "Missing": 0, "DomainBoundary": 2, "NearMatched": 1})
        self.assertTrue(restored["Complete"])
        self.assertEqual(restored["Requirements"][1]["SelectedModels"], [{"Name": "same-conductor-strip-2um"}])   # a real model: kept


if __name__ == "__main__":
    unittest.main()
