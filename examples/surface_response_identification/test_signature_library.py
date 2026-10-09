# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Canonical chording of SpatialEdgeCluster arc portions (option A)."""

import math
import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from surface_response_identification import signature_library as L  # noqa: E402


def arc_portion(center, radius, theta0, theta1, gap_radial, conductor=1):
    p0 = (center[0] + radius * math.cos(theta0), center[1] + radius * math.sin(theta0))
    p1 = (center[0] + radius * math.cos(theta1), center[1] + radius * math.sin(theta1))
    tm = 0.5 * (theta0 + theta1)
    m = (center[0] + radius * math.cos(tm), center[1] + radius * math.sin(tm))
    a, b = sorted([p0, p1])
    return {"P": [a[0], a[1], b[0], b[1]], "Arc": [center[0], center[1], m[0], m[1]], "GapRadial": gap_radial, "Conductor": conductor, "Interfaces": ["MetalAir"], "Law": "{}"}


class ClusterChordingTest(unittest.TestCase):
    R = 2.0

    def test_straight_portion_is_one_edge(self):
        signature = {"Portions": [{"P": [0.0, 0.0, 1.5, 0.0], "Gap": [0.0, 1.0], "Conductor": 1}]}
        edges = L.cluster_plan_view_edges(signature, self.R)
        self.assertEqual(len(edges), 1)
        self.assertEqual(edges[0]["P0"], (0.0, 0.0))
        self.assertEqual(edges[0]["P1"], (3.0, 0.0))
        self.assertEqual(edges[0]["Gap"], (0.0, 1.0))

    def test_quarter_arc_chord_count_and_geometry(self):
        # A 90 deg arc of radius 0.5 R: 90 / 5 = 18 chords by the angular step; the length
        # rule (pi / 4 R over 0.25 R = 4 chords) is weaker.
        signature = {"Portions": [arc_portion((0.0, 0.0), 0.5, 0.0, 0.5 * math.pi, +1)]}
        edges = L.cluster_plan_view_edges(signature, self.R)
        self.assertEqual(len(edges), 18)
        for k, edge in enumerate(edges):
            for point in (edge["P0"], edge["P1"]):
                self.assertAlmostEqual(math.hypot(*point), 0.5 * self.R, places=12)
            self.assertEqual(edge["Chord"], k)
            # Outward gap (metal inside the circle): the gap direction is the outward radial.
            mid = (0.5 * (edge["P0"][0] + edge["P1"][0]), 0.5 * (edge["P0"][1] + edge["P1"][1]))
            self.assertGreater(edge["Gap"][0] * mid[0] + edge["Gap"][1] * mid[1], 0.0)
        # The chording runs from the lexicographically smaller end to the larger one.
        ends = sorted([edges[0]["P0"], edges[-1]["P1"]])
        self.assertAlmostEqual(ends[0][1], 0.5 * self.R, places=12)
        self.assertAlmostEqual(ends[1][0], 0.5 * self.R, places=12)

    def test_length_rule_refines_a_large_radius(self):
        # A 20 deg arc of radius 10 R: 4 chords by the angle, 10 R x 0.349 / 0.25 R = 14 by length.
        signature = {"Portions": [arc_portion((3.0, -2.0), 10.0, 0.3, 0.3 + math.radians(20.0), -1)]}
        edges = L.cluster_plan_view_edges(signature, self.R)
        self.assertEqual(len(edges), 14)
        # Inward gap (metal outside the circle).
        edge = edges[5]
        mid = (0.5 * (edge["P0"][0] + edge["P1"][0]) - 3.0 * self.R, 0.5 * (edge["P0"][1] + edge["P1"][1]) + 2.0 * self.R)
        self.assertLess(edge["Gap"][0] * mid[0] + edge["Gap"][1] * mid[1], 0.0)

    def test_side_of_the_arc_from_the_midpoint(self):
        # The same two end points, the midpoint on the other side: the long way round.
        short = arc_portion((0.0, 0.0), 1.0, -0.25 * math.pi, 0.25 * math.pi, +1)
        long_way = dict(short, Arc=[0.0, 0.0, -1.0 * self.R / self.R, 0.0])
        n_short = len(L.cluster_plan_view_edges({"Portions": [short]}, self.R))
        n_long = len(L.cluster_plan_view_edges({"Portions": [long_way]}, self.R))
        self.assertEqual(n_short, 18)
        self.assertEqual(n_long, 54)

    def test_closed_circle(self):
        portion = arc_portion((0.0, 0.0), 0.3, 0.0, 2.0 * math.pi, +1)
        edges = L.cluster_plan_view_edges({"Portions": [portion]}, self.R)
        self.assertEqual(len(edges), 72)
        self.assertAlmostEqual(edges[0]["P0"][0], edges[-1]["P1"][0], places=12)
        self.assertAlmostEqual(edges[0]["P0"][1], edges[-1]["P1"][1], places=12)

    def test_chording_is_mesh_independent_by_construction(self):
        # Two serialisations of one arc (the builder's input) give identical edges.
        a = arc_portion((1.0, 1.0), 0.75, 0.1, 1.3, +1)
        b = dict(a)
        self.assertEqual(L.cluster_plan_view_edges({"Portions": [a]}, self.R), L.cluster_plan_view_edges({"Portions": [b]}, self.R))




class SpatialSupportContractTest(unittest.TestCase):
    """Contract v3 (decision 282): the claims-derived box and the Context entries."""

    # The 4-edge lead-end cluster of the C++ unit case SurfaceResponseIdentificationSpatialSupportContract
    # (pad edge claimed over 6.873 R, a 2-um lead ending 1 um above it; R = 2): the C++ record's ClaimsBox.
    SIGNATURE = {"Type": "SpatialEdgeCluster", "EdgeCount": 4, "Portions": [
        {"Conductor": 1, "Gap": [-1.0, 0.0], "Interfaces": ["MS"], "Law": "{}", "P": [-0.218559, -0.5, -0.218559, 0.5]},
        {"Conductor": 1, "Gap": [0.0, -1.0], "Interfaces": ["MS"], "Law": "{}", "P": [-0.218559, -0.5, 2.281441, -0.5]},
        {"Conductor": 1, "Gap": [0.0, 1.0], "Interfaces": ["MS"], "Law": "{}", "P": [-0.218559, 0.5, 2.281441, 0.5]},
        {"Conductor": 2, "Gap": [1.0, 0.0], "Interfaces": ["MS"], "Law": "{}", "P": [-0.718559, -3.436492, -0.718559, 3.436492]}],
        "Vertices": [{"P": [-0.218559, -0.5], "TurnDegrees": 90.0, "Type": "ConvexCorner"}, {"P": [-0.218559, 0.5], "TurnDegrees": 90.0, "Type": "ConvexCorner"}]}
    CPP_BOX = [-3.218559, -6.436492, 5.281441, 6.436492]

    def test_box_matches_the_cpp_record(self):
        box = L.cluster_support_box(self.SIGNATURE)
        for value, expected in zip(box, self.CPP_BOX):
            self.assertAlmostEqual(value, expected, places=9)

    def test_box_rule_on_one_free_portion(self):
        # A lone 2 R portion: both ends free (|half| = 1 >= R) -> continued by 2R, widened by R, padded by R.
        signature = {"Portions": [{"P": [-1.0, 0.0, 1.0, 0.0], "Gap": [0.0, 1.0], "Conductor": 1}]}
        self.assertEqual(L.cluster_support_box(signature), [-4.0, -2.0, 4.0, 2.0])
        # A connected short end (a vertex there) keeps its length: no continuation on that side.
        signature["Vertices"] = [{"P": [-0.5, 0.0], "Type": "ConvexCorner", "TurnDegrees": 90.0}]
        signature["Portions"][0]["P"] = [-0.5, 0.0, 0.5, 0.0]
        self.assertEqual(L.cluster_support_box(signature), [-1.5, -2.0, 4.0, 2.0])

    def test_legacy_contract_alias(self):
        signature = dict(self.SIGNATURE)
        signature["Box"] = list(self.CPP_BOX)
        signature["Context"] = [{"Chain": True, "Conductor": 1, "Gap": [0.0, -1.0], "Interfaces": ["MS"], "Law": "{}", "P": [2.281441, -0.5, 5.281441, -0.5]}]
        digest = L.context_digest(signature)
        self.assertEqual(len(digest), 64)
        self.assertEqual(L.context_digest(self.SIGNATURE), "")
        feature = {"Id": 7, "Hash": L.signature_hash(signature), "Signature": signature, "SpatialSupport": {"ContextDigest": digest, "ClaimsKey": L.signature_hash(self.SIGNATURE)}}
        alias = L.legacy_contract_alias(feature, "unit test")
        self.assertEqual(alias["Key"], feature["Hash"])
        self.assertEqual(alias["ContextDigest"], digest)
        self.assertEqual(alias["ClaimsKey"], L.signature_hash(self.SIGNATURE))
        self.assertEqual(alias["Context"]["Box"], self.CPP_BOX)
        with self.assertRaises(ValueError):
            L.legacy_contract_alias({"Id": 8, "Hash": "x", "Signature": self.SIGNATURE}, "claims only")
        with self.assertRaises(ValueError):
            L.legacy_contract_alias(dict(feature, SpatialSupport={"ContextDigest": "0" * 64}), "digest mismatch")

    def test_context_entries_follow_the_claims(self):
        signature = dict(self.SIGNATURE)
        signature["Box"] = list(self.CPP_BOX)
        signature["Context"] = [{"Chain": True, "Conductor": 1, "Gap": [0.0, -1.0], "Interfaces": ["MS"], "Law": "{}", "P": [2.281441, -0.5, 5.281441, -0.5]},
                                {"Chain": False, "Conductor": 3, "Gap": [0.0, 1.0], "Interfaces": ["MS"], "Law": "{}", "P": [-3.0, 5.0, -2.0, 5.0]}]
        claims = L.cluster_plan_view_edges(signature, 2.0)
        self.assertEqual(len(claims), 4)
        self.assertTrue(all("Context" not in e for e in claims))
        edges = L.cluster_plan_view_edges(signature, 2.0, include_context=True)
        self.assertEqual(len(edges), 6)
        self.assertEqual([e.get("Context", False) for e in edges], [False] * 4 + [True] * 2)
        self.assertEqual([e["Chain"] for e in edges[4:]], [True, False])
        self.assertEqual([e["Portion"] for e in edges[4:]], [0, 1])
        self.assertEqual(edges[4]["P0"], (2.0 * 2.281441, -1.0))
        self.assertEqual(L.conductor_count({"Type": "SpatialEdgeCluster", "Signature": signature}), 3)
        self.assertEqual(L.conductor_count({"Type": "SpatialEdgeCluster", "Signature": self.SIGNATURE}), 2)


if __name__ == "__main__":
    unittest.main()


class ClusterQuantumNearMatchTest(unittest.TestCase):
    """Block (b) DESIGN section 4 (decision 303): the Python mirror of the C++ cluster quantum
    near-match on the three stage-2 loop-end keys (S1p 284d6c2b5b66, S2p 9e103a0f291c, S4
    20ac3e14a128: one topology key, 17 / 8 of 290 numbers one quantum apart)."""

    FIXTURE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "cpw3d_surface",
                           "spatial_coupon", "testdata", "census-b-signatures", "signatures.json")

    @classmethod
    def setUpClass(cls):
        import json
        cls.census = json.load(open(cls.FIXTURE))
        cls.s1p = cls.census["284d6c2b5b66"]["Signature"]
        cls.s2p = cls.census["9e103a0f291c"]["Signature"]
        cls.s4 = cls.census["20ac3e14a128"]["Signature"]

    def test_hashes_match_the_cpp_keys(self):
        for prefix in ("284d6c2b5b66", "9e103a0f291c", "20ac3e14a128"):
            self.assertTrue(L.signature_hash(self.census[prefix]["Signature"]).startswith(prefix))

    def test_one_topology_key_and_the_quantum_differences(self):
        t1, l1, a1 = L.split_parameters(self.s1p)
        t2, _, _ = L.split_parameters(self.s2p)
        t4, _, _ = L.split_parameters(self.s4)
        self.assertEqual(t1, t2)
        self.assertEqual(t1, t4)
        self.assertEqual(len(l1) + len(a1), 290)
        worst, paths = L.cluster_quantum_difference(self.s1p, self.s2p)
        self.assertAlmostEqual(worst, 1.0, places=6)
        self.assertEqual(len(paths), 17)
        self.assertIn("Vertices[5].TurnDegrees", paths)
        self.assertIn("Portions[0].Arc[1]", paths)
        worst, paths = L.cluster_quantum_difference(self.s1p, self.s4)
        self.assertAlmostEqual(worst, 1.0, places=6)
        self.assertEqual(len(paths), 8)
        # Normalised by the half-quantum inclusive threshold 4.5 quanta (decision 317 MINOR-1).
        self.assertAlmostEqual(L.signature_deviation(self.s1p, self.s2p), 1.0 / 4.5, places=6)
        self.assertAlmostEqual(L.signature_deviation(self.s2p, self.s4), 1.0 / 4.5, places=6)

    def test_grouping_and_the_representative(self):
        import json
        features = [{"Id": k, "Type": "SpatialEdgeCluster", "Signature": s} for k, s in enumerate((self.s1p, self.s2p, self.s4))]
        groups = L.group_features(features)
        self.assertEqual(len(groups), 1)
        representative, members, spread = groups[0]
        self.assertEqual(json.dumps(representative, sort_keys=True), min(json.dumps(s, sort_keys=True) for s in (self.s1p, self.s2p, self.s4)))
        self.assertEqual(len(members), 3)
        self.assertLessEqual(spread, 1.0)
        # The census keys are Unboxable placeholders (the span cap; step 2 boxes them under an
        # allowance): the library build is exercised on boxable copies.
        boxable = [{**f, "Signature": {k: v for k, v in f["Signature"].items() if k != "Unboxable"}} for f in features]
        for f in boxable:
            self.assertNotIn("Unboxable", f["Signature"])
        manifest = {"Identification": {"MatchingRadius": 1.9, "Features": boxable}}
        library = L.build_signature_library(manifest)
        representative = L.group_features(boxable)[0][0]
        self.assertEqual(len(library["Models"]), 1)
        model = library["Models"][0]
        representative_hash = L.signature_hash(representative)
        self.assertEqual(model["Name"], f"SpatialEdgeCluster-{representative_hash[:12]}")
        self.assertEqual(model["NearKeys"], sorted({L.signature_hash(f["Signature"]) for f in boxable} - {representative_hash}))
        self.assertEqual(len(model["NearKeys"]), 2)

    def test_five_quanta_and_a_permutation_stay_apart(self):
        import copy

        def shifted(quanta):
            out = copy.deepcopy(self.s1p)
            out["Portions"][20]["P"][0] += quanta * 1.0e-6
            return out
        five, four, eight, nine = shifted(5), shifted(4), shifted(8), shifted(9)
        # Exactly 4 / 8 ON-GRID quanta compute to 4 / 8 +- 1e-9: the half-quantum inclusive
        # rule (decision 317 MINOR-1) matches at 4 and refuses the duplicate at 8; 5 / 9 are apart.
        self.assertAlmostEqual(L.cluster_quantum_difference(self.s1p, four)[0], 4.0, places=6)
        self.assertNotEqual(L.cluster_quantum_difference(self.s1p, four)[0], 4.0)
        self.assertTrue(L.within_cluster_quantum_near_match(L.cluster_quantum_difference(self.s1p, four)[0]))
        self.assertFalse(L.within_cluster_quantum_near_match(L.cluster_quantum_difference(self.s1p, five)[0]))
        self.assertTrue(L.cluster_quantum_duplicate(L.cluster_quantum_difference(self.s1p, eight)[0]))
        self.assertFalse(L.cluster_quantum_duplicate(L.cluster_quantum_difference(self.s1p, nine)[0]))
        self.assertGreater(L.signature_deviation(self.s1p, five), 1.0)
        self.assertLessEqual(L.signature_deviation(self.s1p, four), 1.0)
        self.assertEqual(L.cluster_quantum_difference(self.s1p, five)[1], ["Portions[20].P[0]"])
        # The matcher (grouping): 4 quanta group into one case, 5 stay two.
        self.assertEqual(len(L.group_features([{"Id": 0, "Type": "SpatialEdgeCluster", "Signature": self.s1p},
                                               {"Id": 1, "Type": "SpatialEdgeCluster", "Signature": four}])), 1)
        permuted = copy.deepcopy(self.s1p)
        permuted["Portions"][0], permuted["Portions"][1] = permuted["Portions"][1], permuted["Portions"][0]
        self.assertIsNone(L.signature_deviation(self.s1p, permuted))
        features = [{"Id": 0, "Type": "SpatialEdgeCluster", "Signature": self.s1p}, {"Id": 1, "Type": "SpatialEdgeCluster", "Signature": five}]
        self.assertEqual(len(L.group_features(features)), 2)
        # Non-cluster signatures keep the parameter-tolerance comparator.
        a = {"Type": "SameConductorGap", "SeparationOverR": 1.5}
        b = {"Type": "SameConductorGap", "SeparationOverR": 1.5001}
        self.assertAlmostEqual(L.signature_deviation(a, b), 0.1, places=9)


if __name__ == "__main__":
    unittest.main()


# The three-lexeme KeyText fixture (decision 605 (2) F-4 (c); the R2-failure diagnosis F-4): the
# nlohmann (Grisu2) key text of a signature whose doubles 7.378134999999999 / -2.953708 /
# -2.732016 are written as 17-digit lexemes, pinned together with the C++ unit test
# (test-signaturekeytext.cpp "SignatureKeyAndHash key text carries the hashed lexemes").
KEYTEXT_FIXTURE = '{"EdgeCount":2,"Portions":[{"P":[7.3781349999999994,-2.9537079999999998,-2.7320159999999998,0.5]}],"Type":"SpatialEdgeCluster"}'
KEYTEXT_FIXTURE_HASH = "95423956626c93e2982f74f01c7a119a57a8cd44495a054243071caffa1b3779"
KEYTEXT_FIXTURE_FLOAT_HASH = "5e30744232102a28c1997259012c6465418308055fcdb63cb055f9f1fd99da83"


class SignatureKeyTextTest(unittest.TestCase):
    """verify_signature_key: KeyText path (sha256(KeyText) == key and json.loads(KeyText) == Signature), the
    legacy float-hash path for a record without KeyText, and fail-closed by name on every inconsistency."""

    def setUp(self):
        import hashlib
        import json
        self.signature = json.loads(KEYTEXT_FIXTURE)
        self.assertEqual(hashlib.sha256(KEYTEXT_FIXTURE.encode()).hexdigest(), KEYTEXT_FIXTURE_HASH)

    def test_keytext_path_reproduces_the_cpp_digest_and_the_float_hash_does_not(self):
        # T13: the verbatim text hashes to the C++ key; Python's shortest-repr re-dump does not.
        self.assertEqual(L.signature_hash(self.signature), KEYTEXT_FIXTURE_FLOAT_HASH)
        self.assertNotEqual(L.signature_hash(self.signature), KEYTEXT_FIXTURE_HASH)
        record = L.verify_signature_key(self.signature, KEYTEXT_FIXTURE_HASH, KEYTEXT_FIXTURE)
        self.assertEqual(record, {"Path": L.KEY_PATH_KEYTEXT, "KeyText": KEYTEXT_FIXTURE, "SignatureHash": KEYTEXT_FIXTURE_HASH,
                                  "RequirementKey": KEYTEXT_FIXTURE_HASH, "Equal": True})
        # the >= 12-hex prefix form of the key (the 438 (3) rule) and no key at all
        self.assertEqual(L.verify_signature_key(self.signature, KEYTEXT_FIXTURE_HASH[:12], KEYTEXT_FIXTURE)["SignatureHash"],
                         KEYTEXT_FIXTURE_HASH)
        self.assertEqual(L.verify_signature_key(self.signature, None, KEYTEXT_FIXTURE)["SignatureHash"], KEYTEXT_FIXTURE_HASH)
        # the same key WITHOUT KeyText takes the legacy path and fails by name: the record-era defect, unchanged
        with self.assertRaises(L.SignatureKeyError) as context:
            L.verify_signature_key(self.signature, KEYTEXT_FIXTURE_HASH, None)
        self.assertIn("legacy float-hash path", str(context.exception))
        self.assertIn("!= the requirement key", str(context.exception))

    def test_legacy_path_without_keytext_is_the_float_hash(self):
        signature = {"Type": "SpatialEdgeCluster", "EdgeCount": 2, "Portions": [{"P": [0.5, -1.25, 2.0, 0.0]}]}
        key = L.signature_hash(signature)
        record = L.verify_signature_key(signature, key, None)
        self.assertEqual(record, {"Path": L.KEY_PATH_LEGACY, "KeyText": None, "SignatureHash": key, "RequirementKey": key,
                                  "Equal": True})
        self.assertEqual(L.verify_signature_key(signature, key[:12], None)["Path"], L.KEY_PATH_LEGACY)
        self.assertEqual(L.verify_signature_key(signature, None, None)["SignatureHash"], key)
        with self.assertRaises(L.SignatureKeyError):
            L.verify_signature_key(signature, "0" * 64, None)
        with self.assertRaises(L.SignatureKeyError):
            L.verify_signature_key(signature, key[:11], None)  # a prefix shorter than 12 hex is not a key

    def test_inconsistent_keytext_fails_closed_by_name(self):
        import json
        # sha256(KeyText) != the key
        with self.assertRaises(L.SignatureKeyError) as context:
            L.verify_signature_key(self.signature, "0" * 64, KEYTEXT_FIXTURE)
        self.assertIn("sha256(KeyText)", str(context.exception))
        # the text hashes to the key but names another signature
        other = dict(self.signature, EdgeCount=3)
        with self.assertRaises(L.SignatureKeyError) as context:
            L.verify_signature_key(other, KEYTEXT_FIXTURE_HASH, KEYTEXT_FIXTURE)
        self.assertIn("json.loads(KeyText) != Signature", str(context.exception))
        # the float re-dump offered as KeyText: it parses equal but is not the hashed text
        shortest = json.dumps(self.signature, separators=(",", ":"), sort_keys=True)
        with self.assertRaises(L.SignatureKeyError) as context:
            L.verify_signature_key(self.signature, KEYTEXT_FIXTURE_HASH, shortest)
        self.assertIn("sha256(KeyText)", str(context.exception))
        # not JSON / not a string / empty
        import hashlib
        broken = KEYTEXT_FIXTURE[:-1]
        with self.assertRaises(L.SignatureKeyError) as context:
            L.verify_signature_key(self.signature, hashlib.sha256(broken.encode()).hexdigest(), broken)
        self.assertIn("not JSON", str(context.exception))
        for bad in (self.signature, 12, ""):
            with self.assertRaises(L.SignatureKeyError):
                L.verify_signature_key(self.signature, KEYTEXT_FIXTURE_HASH, bad)

    def test_feature_signature_key_and_the_library_model_name(self):
        # A manifest feature carrying KeyText names its signature-only model by the recorded (C++) key;
        # without KeyText (a record-era manifest) the legacy re-derivation is unchanged. The fixture is the
        # S1p census signature (boxable copy) with ONE lexeme of its key text written in a parse-equal
        # 17-digit form (as nlohmann writes some doubles): the float re-dump then hashes elsewhere.
        import hashlib
        import json
        census = json.load(open(ClusterQuantumNearMatchTest.FIXTURE))
        signature = {k: v for k, v in census["284d6c2b5b66"]["Signature"].items() if k != "Unboxable"}
        float_text = json.dumps(signature, separators=(",", ":"), sort_keys=True)
        self.assertIn("90.0,", float_text)
        key_text = float_text.replace("90.0,", "90.000000000000000,", 1)
        self.assertNotEqual(key_text, float_text)
        self.assertEqual(json.loads(key_text), signature)
        key = hashlib.sha256(key_text.encode()).hexdigest()
        float_hash = L.signature_hash(signature)
        self.assertNotEqual(key, float_hash)
        feature = {"Id": 1, "Type": "SpatialEdgeCluster", "Hash": key, "KeyText": key_text, "Signature": signature}
        self.assertEqual(L.feature_signature_key(feature), key)
        legacy = {k: v for k, v in feature.items() if k != "KeyText"}
        self.assertEqual(L.feature_signature_key(legacy), float_hash)
        with self.assertRaises(L.SignatureKeyError):
            L.feature_signature_key({**feature, "Hash": float_hash})
        manifest = {"Identification": {"MatchingRadius": 1.9, "Features": [feature]}}
        self.assertEqual(L.build_signature_library(manifest)["Models"][0]["Name"], f"SpatialEdgeCluster-{key[:12]}")
        legacy_manifest = {"Identification": {"MatchingRadius": 1.9, "Features": [legacy]}}
        self.assertEqual(L.build_signature_library(legacy_manifest)["Models"][0]["Name"],
                         f"SpatialEdgeCluster-{float_hash[:12]}")
        # two features of one key: still the recorded key, no NearKeys; an inconsistent KeyText on a member
        # fails the library build by name
        manifest["Identification"]["Features"].append({**feature, "Id": 2})
        model = L.build_signature_library(manifest)["Models"][0]
        self.assertEqual((model["Name"], model["Instances"], model.get("NearKeys")), (f"SpatialEdgeCluster-{key[:12]}", 2, None))
        manifest["Identification"]["Features"][1]["Hash"] = "0" * 64
        with self.assertRaises(L.SignatureKeyError):
            L.build_signature_library(manifest)
        # a group of several distinct signatures (the S1p / S2p / S4 near-match group) is re-derived as before:
        # the representative is not any member's signature; the members' NearKeys use their recorded keys
        others = [{k: v for k, v in census[prefix]["Signature"].items() if k != "Unboxable"} for prefix in ("9e103a0f291c", "20ac3e14a128")]
        group = [feature] + [{"Id": 2 + i, "Type": "SpatialEdgeCluster", "Hash": L.signature_hash(o), "Signature": o}
                             for i, o in enumerate(others)]
        model = L.build_signature_library({"Identification": {"MatchingRadius": 1.9, "Features": group}})["Models"][0]
        representative = L.group_features(group)[0][0]
        self.assertEqual(model["Name"], f"SpatialEdgeCluster-{L.signature_hash(representative)[:12]}")
        self.assertEqual(model["NearKeys"], sorted({key} | {L.signature_hash(o) for o in others}))

if __name__ == "__main__":
    unittest.main()
