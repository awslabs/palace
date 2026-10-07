# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""nearkey_predictor: the RuleVersion v1 predictor reproduces the DESIGN v2 per-pair bounds (the six
measured pairs, the four consistency combinations, the tabulated bounds at |W| = 1.0-2.5 %, the
domain corner and validation pair 5), the MA_sharp form and the Option-A policy decisions
(B2 / B4 / B5 Default with the SHIPPED rule since validation pair 5 (decision 459), B3 / C1 / C2 Fallback-only,
Option B B4 only, |W| 15 % refused, DefaultActive false -> fallback-only, fail closed on every missing record);
the rule file's and the shipped activation record's BYTE sha256 are pinned (decision 440 MINOR-A / 460)."""
import atexit
import copy
import hashlib
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import nearkey_predictor as predictor  # noqa: E402

EVIDENCE = Path("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/stage2-20261004/nearkey-reuse-design/evidence")
# decision 459 / 460: the validation-pair-5 DefaultActivation record of the evidence tree and its byte-identical copy shipped beside
# the rule file (the rule's RecordPath); both sha256s and the rule file's own BYTE sha256 are pinned here (decision 440 MINOR-A)
PAIR5_EVIDENCE_RECORD = Path("/Users/simlap/bedrock-tests/coupon-accuracy-assessment-20260913/stage2-20261004/nearkey-validation-pair5/records/"
                             "default-activation-pair5.json")
SHIPPED_RECORD = HERE / "nearkey-default-activation-pair5.json"
SHIPPED_RECORD_SHA256 = "cdc8c059901357564fc5ba724f26013741aac269f45aac5d2572569eb0654606"
RULE_FILE_SHA256 = "7d40555bd2dee115840072f6a11376e6538a27980eee598ac0554bb3d83bfdac"
# DESIGN v2 section 2.2 / evidence calibration.json Data + transplant-fidelity-energy.json (2 x max_X T2e_X is the
# transplant term) + ma-tail.json (the donor's shell tail): the recorded inputs of the six measured pairs.
PAIRS = {
    "B4": {"W": 0.003123125975976694, "S": 0.00020145052418115068, "T2eMax": 8.219448546808002e-08, "T2": True,
           "tail": 0.0213810413806117, "r": {"SA": -0.0008414719204352661, "MS": -0.00045955586489099254, "MA": -0.0006728247146722266}},
    "B2": {"W": -0.009527170077628535, "S": 0.0007573071691547667, "T2eMax": 1.755583837093352e-06, "T2": True,
           "tail": 0.022262609816823264, "r": {"SA": 0.0022163908763748186, "MS": 0.0007679103460078718, "MA": -0.0009982347106820555}},
    "B3": {"W": 0.014625330182121218, "S": 0.001268858969695948, "T2eMax": 0.0007953881034261512, "T2": False,
           "tail": 0.022262609816823264, "r": {"SA": -0.003309229228237176, "MS": -0.0018355484905907549, "MA": -0.005005444604594067}},
    "B5": {"W": 0.015825391576225978, "S": 0.00101035359560976, "T2eMax": 1.7179236635988526e-06, "T2": True,
           "tail": 0.0213810413806117, "r": {"SA": -0.0033612962206749364, "MS": -0.0016866130019945746, "MA": -0.004498547595610747}},
    "C1": {"W": -0.06064641022620769, "S": 0.0073833166876292955, "T2eMax": 6.603141345834059e-06, "T2": True,
           "tail": 0.02185151742987057, "r": {"SA": 0.01148419074978424, "MS": 0.0076924154900468444, "MA": 0.009543470326705439}},
    "C2": {"W": -0.11280045501005607, "S": 0.006391354806027371, "T2eMax": 8.391083511532648e-06, "T2": True,
           "tail": 0.02071196223047833, "r": {"SA": 0.023750688778903628, "MS": 0.015082181268420536, "MA": 0.025068248990486763}},
}
# DESIGN 2.4 "Per pair (governing form; |r| / bound)": bounds in % (3 decimals) and the governing form.
DESIGN_PAIR_BOUNDS = {"B4": (0.088, 0.054, 0.401, "Default"), "B2": (0.267, 0.150, 0.635, "Default"),
                      "B3": (0.568, 0.385, 0.981, "Default"), "B5": (0.440, 0.244, 0.860, "Default"),
                      "C1": (1.468, 0.923, 1.834, "Fallback"), "C2": (2.643, 1.677, 3.071, "Fallback")}
# DESIGN 2.5: the four consistency combinations (W, S, T2e, r_pre) and their bounds (%).
COMBOS = {
    "B3<-B2": ({"W": 0.024384818573612677, "S": 0.002027140225541266, "T2eMax": 0.0007936918077775882}, (0.839, 0.531, 1.335),
               {"SA": -0.005586194690442525, "MS": -0.002667295014184057, "MA": -0.00404518325423342}),
    "B2<-B3": ({"W": -0.023804353726724403, "S": 0.0020245361649786076, "T2eMax": 8.867212532372998e-07}, (0.665, 0.364, 1.156),
               {"SA": 0.005717607328518115, "MS": 0.003252646711990659, "MA": 0.004147410257180795}),
    "B5<-B4": ({"W": 0.012662718335688617, "S": 0.0008088582240332809, "T2eMax": 2.8474566639532766e-06}, (0.352, 0.197, 0.746),
               {"SA": -0.0025594052250899058, "MS": -0.0012728772982094627, "MA": -0.0038923543949905826}),
    "B4<-B5": ({"W": -0.012504378907618716, "S": 0.0008086781941683067, "T2eMax": 4.0555300831476113e-07}, (0.348, 0.194, 0.740),
               {"SA": 0.0025478807499681455, "MS": 0.0012673224046175768, "MA": 0.0038997837220227094}),
}
STOP = {"Path": "registration/stop-record.json", "SHA256": "0" * 64}


_ACTIVATION_FIXTURE = {}


def activation_fixture():
    """The synthetic DefaultActivation record of these tests (also imported by test_nearkey_reuse): written
    once per process under a temporary directory that is removed at interpreter exit - a module-level
    mkdtemp leaked one nearkey-activation-* directory per run.  {Directory, Record, SHA256}."""
    if not _ACTIVATION_FIXTURE:
        directory = Path(tempfile.mkdtemp(prefix="nearkey-activation-"))
        atexit.register(shutil.rmtree, directory, True)
        record = directory / "pair-5-record.json"
        record.write_text(json.dumps({"ValidationPair": "pair-5", "r_pre": {"SA": -0.0046, "MS": -0.0022, "MA": -0.0055},
                                      "InsideBound": True}))
        _ACTIVATION_FIXTURE.update(Directory=directory, Record=record, SHA256=hashlib.sha256(record.read_bytes()).hexdigest())
    return _ACTIVATION_FIXTURE


def activated(rule, record_path=None, record_sha=None):
    """The rule with a synthetic DefaultActivation record (a readable validation-pair record that re-hashes to RecordSHA256)."""
    active = copy.deepcopy(rule)
    active["Policy"]["DefaultActive"] = True
    fixture = activation_fixture()
    active["Policy"]["DefaultActivation"] = {"ValidationPair": "pair-5", "RecordPath": str(record_path or fixture["Record"]),
                                             "RecordSHA256": record_sha or fixture["SHA256"], "Decision": "synthetic (test)"}
    return active


def deactivated(rule):
    """The rule as it shipped before decision 459 (DefaultActive false, no DefaultActivation record): the fallback-only state."""
    inactive = copy.deepcopy(rule)
    inactive["Policy"]["DefaultActive"] = False
    inactive["Policy"]["DefaultActivation"] = None
    return inactive


class RuleFile(unittest.TestCase):
    def test_rule_pins_design_v2(self):
        rule = predictor.load_rule()
        self.assertEqual(rule["Design"]["SHA256"], "e9c10d84af99b9ecb5cb29c840ec0f6a3a16ea8f249993e458b5b15ec4351f50")
        self.assertEqual(rule["Policy"]["Option"], "A")
        self.assertEqual(rule["Policy"]["DefaultBound"], {"SA": 0.005, "MS": 0.005, "MA_sharp": 0.010})
        self.assertEqual(rule["Policy"]["FallbackMaxBound"], 0.035)
        self.assertEqual(rule["Policy"]["Domain"]["WDefaultMax"], 0.025)
        # decision 459: DEFAULT reuse is ACTIVE by the validation-pair-5 record (the policy limits unchanged: Option A)
        self.assertTrue(rule["Policy"]["DefaultActive"])
        activation = rule["Policy"]["DefaultActivation"]
        self.assertEqual(activation["ValidationPair"], "pair-5")
        self.assertEqual(activation["RecordPath"], "nearkey-default-activation-pair5.json")
        self.assertEqual(activation["RecordSHA256"], SHIPPED_RECORD_SHA256)
        self.assertEqual((activation["Decision"], activation["UserDecision"], activation["Verdict"]), (459, 428, "PASS"))
        self.assertEqual(activation["Bound_pct"], {"SA": 0.551, "MS": 0.301, "MA_sharp": 1.007})
        for T, r in activation["Measured_r_pct"].items():
            self.assertLess(abs(r), activation["Bound_pct"][T], T)
        self.assertLess(activation["Measured_r_pct"]["SA"], 0)
        self.assertLess(activation["Measured_r_pct"]["MS"], 0)
        keys = sorted(k[:16] for k in predictor.admissible_structure_keys(rule))
        self.assertEqual(keys, ["2d27be59ff1b0774", "4682966a36a26319", "52c670f3cfe5f26b", "d44841a19b9546cf"])
        c = rule["Predictor"]["Coefficients"]
        for T, a_def, rho_def, a_full, rho_full, b, phi_def, phi_full in (
                ("SA", 0.2352, 0.15, 0.2063, 0.10, 0.093, 0.002, 0.021), ("MS", 0.1150, 0.30, 0.1314, 0.10, 0.004, 0.007, 0.042),
                ("MA", 0.2806, 0.25, 0.2088, 0.15, 0.160, 0.288, 0.258)):
            self.assertAlmostEqual(c[T]["Default"]["a"], a_def, places=4)
            self.assertAlmostEqual(c[T]["Default"]["rho"], rho_def, places=9)
            self.assertAlmostEqual(c[T]["Fallback"]["a"], a_full, places=4)
            self.assertAlmostEqual(c[T]["Fallback"]["rho"], rho_full, places=9)
            self.assertAlmostEqual(c[T]["Default"]["b"], b, places=3)
            self.assertAlmostEqual(c[T]["Fallback"]["b"], b, places=3)
            self.assertAlmostEqual(100 * c[T]["Default"]["phi"], phi_def, places=3)
            self.assertAlmostEqual(100 * c[T]["Fallback"]["phi"], phi_full, places=3)

    def test_rule_file_byte_sha256_is_pinned(self):
        """Decision 440 MINOR-A / 459: the shipped rule file's BYTE sha256 (= RuleFileSHA256 on every prediction and library header)."""
        self.assertEqual(hashlib.sha256(predictor.RULE_FILE.read_bytes()).hexdigest(), RULE_FILE_SHA256)
        self.assertEqual(predictor.load_rule()["_sha256"], RULE_FILE_SHA256)

    def test_shipped_activation_record_is_pinned(self):
        """Decision 460: the activation record shipped beside the rule file re-hashes to the rule's RecordSHA256 and names pair 5 PASS."""
        self.assertEqual(hashlib.sha256(SHIPPED_RECORD.read_bytes()).hexdigest(), SHIPPED_RECORD_SHA256)
        record = json.loads(SHIPPED_RECORD.read_text())
        self.assertEqual((record["Record"], record["ValidationPair"], record["Verdict"]), ("DefaultActivation", "pair-5", "PASS"))
        self.assertEqual(record["RuleVersion"], predictor.RULE_VERSION)
        for T, bound in (("SA", 0.551), ("MS", 0.301), ("MA_sharp", 1.007)):
            scored = record["Scoring"]["PerType"][T]
            self.assertEqual(scored["Bound_pct"], bound)
            self.assertTrue(scored["InsideBound"])
            self.assertLess(abs(scored["r_pct"]), bound)
        self.assertTrue(record["Scoring"]["PerType"]["SA"]["SignRight"])
        self.assertTrue(record["Scoring"]["PerType"]["MS"]["SignRight"])
        self.assertEqual(len(record["Records"]["QExactModelFiles"]["Files"]), 7)
        active, activation = predictor.default_active(predictor.load_rule())
        self.assertTrue(active)
        self.assertEqual(activation["RecordSHA256"], SHIPPED_RECORD_SHA256)
        self.assertEqual(predictor.activation_record_path(predictor.load_rule(), activation), SHIPPED_RECORD)

    @unittest.skipUnless(PAIR5_EVIDENCE_RECORD.is_file(), "the validation-pair-5 evidence tree is not mounted")
    def test_shipped_activation_record_equals_the_evidence_record(self):
        self.assertEqual(SHIPPED_RECORD.read_bytes(), PAIR5_EVIDENCE_RECORD.read_bytes())
        self.assertEqual(hashlib.sha256(PAIR5_EVIDENCE_RECORD.read_bytes()).hexdigest(), SHIPPED_RECORD_SHA256)

    @unittest.skipUnless((EVIDENCE / "calibration.json").is_file(), "the design evidence tree is not mounted")
    def test_rule_matches_calibration_record(self):
        rule = predictor.load_rule()
        payload = (EVIDENCE / "calibration.json").read_bytes()
        import hashlib
        self.assertEqual(hashlib.sha256(payload).hexdigest(), rule["CalibrationRecord"]["SHA256"])
        calibration = json.loads(payload)
        for T in predictor.TYPES:
            d, f = calibration["DefaultRegime"][T], calibration["Final"][T]
            self.assertEqual(rule["Predictor"]["Coefficients"][T]["Default"]["a"], d["a_def"])
            self.assertEqual(rule["Predictor"]["Coefficients"][T]["Fallback"]["a"], f["a"])
            self.assertAlmostEqual(rule["Predictor"]["Coefficients"][T]["Default"]["phi"], d["phi_def_pct"] / 100, places=15)

    def test_malformed_rule_fails_closed(self):
        rule = json.loads(predictor.RULE_FILE.read_text())
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "rule.json"
            bad = copy.deepcopy(rule)
            bad["RuleVersion"] = "nearkey-reuse-rule-v2"
            path.write_text(json.dumps(bad))
            with self.assertRaises(predictor.NearKeyRuleError):
                predictor.load_rule(path)
            bad = copy.deepcopy(rule)
            bad["AdmissibleStructureKeys"].append({"StructureKey": "a" * 64})
            path.write_text(json.dumps(bad))
            with self.assertRaises(predictor.NearKeyRuleError):
                predictor.load_rule(path)
            bad = copy.deepcopy(rule)
            bad["Predictor"]["Coefficients"]["SA"]["Default"]["a"] = -1.0
            path.write_text(json.dumps(bad))
            with self.assertRaises(predictor.NearKeyRuleError):
                predictor.load_rule(path)


class Bounds(unittest.TestCase):
    def setUp(self):
        self.rule = predictor.load_rule()

    def prediction(self, data):
        return predictor.predict(self.rule, w=data["W"], s=data["S"], t2e_max=data["T2eMax"], tail_donor=data.get("tail", 0.0))

    def test_six_pairs_reproduce_design_bounds_and_cover_the_measurements(self):
        for pair, data in PAIRS.items():
            p = self.prediction(data)
            sa, ms, ma, regime = DESIGN_PAIR_BOUNDS[pair]
            self.assertEqual(p["Regime"], regime, pair)
            for T, expected in zip(predictor.TYPES, (sa, ms, ma)):
                self.assertAlmostEqual(100 * p[T]["Bound"], expected, places=3, msg=f"{pair} {T}")
                ratio = abs(data["r"][T]) / p[T]["Bound"]
                self.assertTrue(0.155 <= ratio <= 0.955, f"{pair} {T} |r| / bound {ratio}")   # DESIGN: 0.16-0.95 (rounded)
                if not (pair == "B2" and T == "MA"):   # the one MA sign flip of the data (DESIGN 2.3)
                    self.assertGreater(p[T]["Central"] * data["r"][T], 0, f"{pair} {T} sign")
            self.assertEqual(p["SA"]["Terms"]["Transplant"], 2 * data["T2eMax"])
        self.assertEqual(self.prediction(PAIRS["C1"])["SA"]["GoverningForm"], "Fallback")

    def test_consistency_combinations(self):
        for name, (inputs, expected, r_pre) in COMBOS.items():
            p = predictor.predict(self.rule, w=inputs["W"], s=inputs["S"], t2e_max=inputs["T2eMax"], tail_donor=0.0)
            for T, value in zip(predictor.TYPES, expected):
                self.assertAlmostEqual(100 * p[T]["Bound"], value, places=3, msg=f"{name} {T}")
                self.assertLessEqual(abs(r_pre[T]) / p[T]["Bound"], 0.895, f"{name} {T}")   # DESIGN: <= 0.89 (rounded)

    @unittest.skipUnless((EVIDENCE / "calibration.json").is_file(), "the design evidence tree is not mounted")
    def test_bounds_equal_the_calibration_record_to_1e_9(self):
        calibration = json.loads((EVIDENCE / "calibration.json").read_text())
        crossval = json.loads((EVIDENCE / "crossval.json").read_text())["Combos"]
        for pair, data in PAIRS.items():
            p = self.prediction(data)
            for T in predictor.TYPES:
                self.assertAlmostEqual(100 * p[T]["Bound"], calibration["DefaultRegime"][T]["pairs"][pair]["bound_rule_pct"], delta=1e-9)
        for name, combo in crossval.items():
            p = predictor.predict(self.rule, w=combo["W"], s=combo["S"], t2e_max=combo["Tests"]["T2e"], tail_donor=0.0)
            for T in predictor.TYPES:
                self.assertAlmostEqual(100 * p[T]["Bound"], calibration["DefaultRegime"][T]["combos"][name]["bound_rule_pct"], delta=1e-9)

    def test_tabulated_bounds_at_s_equal_0p1_w(self):
        # DESIGN 2.4: |W| 1.0 -> 0.282 / 0.157 / 0.655; 1.5 -> 0.422 / 0.232 / 0.838; 1.7 -> 0.478 / 0.262 / 0.912;
        # 2.0 -> 0.561 / 0.307 / 1.022; 2.5 -> 0.701 / 0.382 / 1.205; the domain corner (11.28 %, 0.75 %) -> 2.652 / 1.676 / 3.088.
        table = {0.010: (0.282, 0.157, 0.655), 0.015: (0.422, 0.232, 0.838), 0.017: (0.478, 0.262, 0.912),
                 0.020: (0.561, 0.307, 1.022), 0.025: (0.701, 0.382, 1.205)}
        for w, expected in table.items():
            p = predictor.predict(self.rule, w=w, s=0.1 * w, t2e_max=0.0, tail_donor=0.0)
            for T, value in zip(predictor.TYPES, expected):
                self.assertAlmostEqual(100 * p[T]["Bound"], value, places=3, msg=f"W {w} {T}")
        corner = predictor.predict(self.rule, w=0.1128, s=0.0075, t2e_max=0.0, tail_donor=0.0)
        for T, value in zip(predictor.TYPES, (2.652, 1.676, 3.088)):
            self.assertAlmostEqual(100 * corner[T]["Bound"], value, places=3)

    def test_regime_boundary_is_monotone(self):
        below = predictor.predict(self.rule, w=0.025, s=0.0025, t2e_max=0.0, tail_donor=0.0)
        above = predictor.predict(self.rule, w=0.0250001, s=0.0025, t2e_max=0.0, tail_donor=0.0)
        self.assertEqual(below["Regime"], "Default")
        self.assertEqual(above["Regime"], "Fallback")
        for T in predictor.TYPES:
            self.assertGreaterEqual(above[T]["Bound"], below[T]["Bound"] - 1e-15)
        # just above the boundary the default form at WDefaultMax governs (the full-range form is lower there)
        self.assertEqual(above["SA"]["GoverningForm"], "DefaultAtWDefaultMax")
        wide = predictor.predict(self.rule, w=-0.03, s=0.003, t2e_max=0.0, tail_donor=0.0)
        self.assertEqual(wide["SA"]["GoverningForm"], "Fallback")
        self.assertGreater(wide["SA"]["Bound"], above["SA"]["Bound"])

    def test_validation_pair_5_prediction(self):
        # DESIGN 5.1 item 5 / P-NK2 (recomputed by the confirmation review MINOR-4): W +1.961 %, S 0.2 % ->
        # central SA -0.461 / MS -0.225 / MA -0.550 %, bounds 0.551 / 0.301 / 1.008 %.
        w = (0.208 - 0.204) / 0.204
        p = predictor.predict(self.rule, w=w, s=0.002, t2e_max=0.0, tail_donor=0.0)
        for T, central, bound in (("SA", -0.461, 0.551), ("MS", -0.225, 0.301), ("MA", -0.550, 1.008)):
            self.assertAlmostEqual(100 * p[T]["Central"], central, places=3)
            self.assertAlmostEqual(100 * p[T]["Bound"], bound, places=3)
        # pair 4 (W +5.25 %, S 0.53 %): the full form governs: -1.083 / -0.690 / -1.096, bounds 1.262 / 0.803 / 1.604
        p4 = predictor.predict(self.rule, w=(0.321 - 0.305) / 0.305, s=0.0053, t2e_max=0.0, tail_donor=0.0)
        for T, central, bound in (("SA", -1.083, 1.262), ("MS", -0.690, 0.803), ("MA", -1.096, 1.604)):
            self.assertAlmostEqual(100 * p4[T]["Central"], central, places=2)
            self.assertAlmostEqual(100 * p4[T]["Bound"], bound, places=2)

    def test_ma_sharp_form(self):
        p = self.prediction(PAIRS["C2"])
        self.assertAlmostEqual(p["MA_sharp"]["Bound"], p["MA"]["Bound"] * (1 + PAIRS["C2"]["tail"]) + 0.00013512593462183808, places=15)
        self.assertTrue(p["MA_sharp"]["TailFromDonor"])
        with self.assertRaises(predictor.NearKeyRuleError):
            predictor.predict(self.rule, w=0.01, s=0.001, t2e_max=0.0, tail_donor=-0.01)

    def test_implied_default_tolerance_is_sa_binding_near_1p7_percent(self):
        # DESIGN 3 / the confirmation review: SA 1.78 %, MS 3.29 %, MA_sharp 1.84 % (t_donor ~ 2.2 %) -> SA-binding -> |W| <= 1.7 %;
        # Option B (0.5 % on MA_sharp too) -> 0.51 % (B4 at 0.312 % still admitted)
        tolerance = predictor.implied_default_tolerance(self.rule, tail_donor=0.0222)
        self.assertEqual(tolerance["Binding"], "SA")
        self.assertAlmostEqual(100 * tolerance["SA"], 1.78, places=2)
        self.assertAlmostEqual(100 * tolerance["MS"], 3.29, places=2)
        self.assertAlmostEqual(100 * tolerance["MA_sharp"], 1.84, places=1)
        self.assertLess(tolerance["SA"], 0.0179)
        self.assertGreater(tolerance["SA"], 0.017)   # => the implied Option-A tolerance |W| <= 1.7 %
        strict = predictor.implied_default_tolerance(self.rule, tail_donor=0.0222, limits={"SA": 0.005, "MS": 0.005, "MA_sharp": 0.005})
        self.assertEqual(strict["Binding"], "MA_sharp")
        self.assertAlmostEqual(100 * strict["MA_sharp"], 0.51, places=2)

    def test_bad_inputs_fail_closed(self):
        for kwargs in ({"w": float("nan"), "s": 0.0, "t2e_max": 0.0}, {"w": 0.01, "s": -1e-3, "t2e_max": 0.0},
                       {"w": 0.01, "s": 0.0, "t2e_max": -1e-9}, {"w": "0.01", "s": 0.0, "t2e_max": 0.0}):
            with self.assertRaises(predictor.NearKeyRuleError):
                predictor.predict(self.rule, tail_donor=0.0, **kwargs)


class Policy(unittest.TestCase):
    def setUp(self):
        self.rule = predictor.load_rule()
        self.active = activated(self.rule)

    def decide(self, rule, pair, mode, **overrides):
        data = PAIRS[pair]
        p = predictor.predict(rule, w=data["W"], s=data["S"], t2e_max=data["T2eMax"], tail_donor=data["tail"])
        kwargs = {"requested_mode": mode, "t2_passed": data["T2"], "gates_passed": True, "in_domain": True}
        kwargs.update(overrides)
        return predictor.policy_decision(rule, p, **kwargs)

    def test_option_a_scoring_of_the_six_pairs_once_activated(self):
        for pair in ("B4", "B2", "B5"):
            decision = self.decide(self.active, pair, "default")
            self.assertTrue(decision["Allowed"], (pair, decision["Reasons"]))
            self.assertEqual(decision["Mode"], "Default")
            self.assertEqual(decision["DefaultActivation"]["ValidationPair"], "pair-5")
        for pair, reason in (("B3", "T2TraceGateFailed"), ("B3", "BoundAbovePolicyLimit: SA"), ("C1", "OutsideDefaultRegime"),
                             ("C2", "OutsideDefaultRegime")):
            decision = self.decide(self.active, pair, "default")
            self.assertFalse(decision["Allowed"], pair)
            self.assertTrue(any(r.startswith(reason) for r in decision["Reasons"]), (pair, reason, decision["Reasons"]))
        for pair in ("B3", "C1", "C2", "B2", "B4", "B5"):
            decision = self.decide(self.active, pair, "fallback", stop_record=STOP, approval="supervisor decision NNN")
            self.assertTrue(decision["Allowed"], (pair, decision["Reasons"]))
            self.assertEqual(decision["Mode"], "Fallback")

    def test_shipped_rule_allows_default_reuse_of_b4_b2_b5_and_refuses_b3_c1_c2(self):
        """Decision 459 / DESIGN v2 section 3 with the SHIPPED rule (no synthetic record): B4 / B2 / B5 are Default reuses now, B3 (T2
        failed, SA bound 0.568 > 0.5 %), C1 and C2 (|W| > 2.5 %) stay refused in Default and remain Fallback-only."""
        shipped = self.rule
        self.assertTrue(predictor.is_shipped_rule(shipped))
        for pair in ("B4", "B2", "B5"):
            decision = self.decide(shipped, pair, "default")
            self.assertTrue(decision["Allowed"], (pair, decision["Reasons"]))
            self.assertEqual(decision["Mode"], "Default")
            self.assertEqual(decision["DefaultActivation"], shipped["Policy"]["DefaultActivation"])
            self.assertEqual(decision["DefaultActivation"]["Decision"], 459)
        for pair, reason in (("B3", "BoundAbovePolicyLimit: SA"), ("B3", "T2TraceGateFailed"), ("C1", "OutsideDefaultRegime"),
                             ("C2", "OutsideDefaultRegime")):
            decision = self.decide(shipped, pair, "default")
            self.assertFalse(decision["Allowed"], pair)
            self.assertEqual(decision["Mode"], "None")
            self.assertTrue(any(r.startswith(reason) for r in decision["Reasons"]), (pair, reason, decision["Reasons"]))
            self.assertFalse(any(r.startswith("DefaultNotActive") for r in decision["Reasons"]), decision["Reasons"])
            self.assertTrue(self.decide(shipped, pair, "fallback", stop_record=STOP, approval="x")["Allowed"], pair)

    def test_relative_record_path_resolves_against_the_rule_file_directory(self):
        """Decision 460: a relative RecordPath is read beside the rule file that carries it (the shipped copy for the shipped rule; a
        rule loaded from another directory looks there and fails closed when the record is absent); absolute paths are read as given."""
        shipped = self.rule
        self.assertEqual(predictor.activation_record_path(shipped, shipped["Policy"]["DefaultActivation"]), SHIPPED_RECORD)
        self.assertEqual(predictor.activation_record_path(shipped, {"RecordPath": str(activation_fixture()["Record"])}),
                         activation_fixture()["Record"])
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "rule.json"
            path.write_text(json.dumps({k: v for k, v in shipped.items() if not k.startswith("_")}))
            elsewhere = predictor.load_rule(path)
            self.assertEqual(predictor.activation_record_path(elsewhere, elsewhere["Policy"]["DefaultActivation"]),
                             (Path(tmp) / "nearkey-default-activation-pair5.json").resolve())
            with self.assertRaises(predictor.NearKeyRuleError):
                predictor.default_active(elsewhere)
            (Path(tmp) / "nearkey-default-activation-pair5.json").write_bytes(SHIPPED_RECORD.read_bytes())
            self.assertTrue(predictor.default_active(elsewhere)[0])

    def test_shipped_record_mismatch_still_fails_closed(self):
        """Decision 438 (5) on the SHIPPED rule: a RecordSHA256 that does not re-hash, a RecordPath to another file or a missing file
        keeps default reuse inactive (NearKeyRuleError), and so does the flag without a record."""
        tampered = copy.deepcopy(self.rule)
        tampered["Policy"]["DefaultActivation"]["RecordSHA256"] = "0" * 64
        with self.assertRaises(predictor.NearKeyRuleError) as context:
            predictor.default_active(tampered)
        self.assertIn("default reuse stays inactive", str(context.exception))
        other_file = copy.deepcopy(self.rule)
        other_file["Policy"]["DefaultActivation"]["RecordPath"] = str(activation_fixture()["Record"])   # readable, re-hashes to another sha
        with self.assertRaises(predictor.NearKeyRuleError):
            predictor.default_active(other_file)
        missing = copy.deepcopy(self.rule)
        missing["Policy"]["DefaultActivation"]["RecordPath"] = "nearkey-default-activation-pair6.json"
        with self.assertRaises(predictor.NearKeyRuleError):
            predictor.default_active(missing)
        for broken in (tampered, other_file, missing):
            with self.assertRaises(predictor.NearKeyRuleError):
                self.decide(broken, "B4", "default")
        # the information-only shape used by the non-Default evaluations (decision 460 (B))
        refused = predictor.refused_decision(self.rule, "default", "DefaultNotActive: record unreadable (test)")
        self.assertEqual((refused["Mode"], refused["Allowed"], refused["Reasons"]), ("None", False, ["DefaultNotActive: record unreadable (test)"]))
        self.assertEqual(refused["PolicyLimit"]["Default"], self.rule["Policy"]["DefaultBound"])

    def test_option_b_admits_b4_only(self):
        strict = copy.deepcopy(self.active)
        strict["Policy"]["DefaultBound"]["MA_sharp"] = 0.005
        admitted = [pair for pair in PAIRS if self.decide(strict, pair, "default")["Allowed"]]
        self.assertEqual(admitted, ["B4"])

    def test_default_inactive_until_a_default_activation_record_exists(self):
        inactive = deactivated(self.rule)   # the rule as shipped before decision 459
        for pair in ("B4", "B2", "B5"):
            decision = self.decide(inactive, pair, "default")
            self.assertFalse(decision["Allowed"])
            self.assertEqual(decision["Mode"], "None")
            self.assertTrue(any(r.startswith("DefaultNotActive") for r in decision["Reasons"]), decision["Reasons"])
            self.assertTrue(self.decide(inactive, pair, "fallback", stop_record=STOP, approval="x")["Allowed"])
        flag_only = copy.deepcopy(inactive)
        flag_only["Policy"]["DefaultActive"] = True   # no record: fail closed
        with self.assertRaises(predictor.NearKeyRuleError):
            self.decide(flag_only, "B4", "default")
        # decision 438 (5): the activation record is verified against the validation-pair record: a missing path, an unreadable
        # path or a sha mismatch fails closed
        partial = copy.deepcopy(inactive)
        partial["Policy"]["DefaultActive"] = True
        partial["Policy"]["DefaultActivation"] = {"ValidationPair": "pair-5", "RecordSHA256": activation_fixture()["SHA256"]}
        with self.assertRaises(predictor.NearKeyRuleError):
            predictor.default_active(partial)
        with self.assertRaises(predictor.NearKeyRuleError):
            predictor.default_active(activated(inactive, record_path=activation_fixture()["Directory"] / "missing.json"))
        with self.assertRaises(predictor.NearKeyRuleError):
            predictor.default_active(activated(inactive, record_sha="0" * 64))
        active, activation = predictor.default_active(self.active)
        self.assertTrue(active)
        self.assertEqual(activation["RecordSHA256"], activation_fixture()["SHA256"])

    def test_only_the_shipped_rule_file_is_the_shipped_rule(self):
        self.assertTrue(predictor.is_shipped_rule(self.rule))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "rule.json"
            path.write_text(predictor.RULE_FILE.read_text())
            self.assertFalse(predictor.is_shipped_rule(predictor.load_rule(path)))

    def test_fallback_needs_stop_record_approval_and_the_cap(self):
        self.assertFalse(self.decide(self.rule, "C2", "fallback", approval="x")["Allowed"])
        self.assertFalse(self.decide(self.rule, "C2", "fallback", stop_record=STOP)["Allowed"])
        self.assertFalse(self.decide(self.rule, "C2", "fallback", stop_record={"Path": "x"}, approval="x")["Allowed"])
        wide = predictor.predict(self.rule, w=0.15, s=0.0075, t2e_max=0.0, tail_donor=0.02)
        decision = predictor.policy_decision(self.rule, wide, requested_mode="fallback", t2_passed=True, gates_passed=True,
                                             in_domain=False, in_domain_reasons=["|W| 0.15 > 0.1128"], stop_record=STOP, approval="x")
        self.assertFalse(decision["Allowed"])
        self.assertTrue(any(r.startswith("OutsideCalibratedDomain") for r in decision["Reasons"]))
        self.assertTrue(any(r.startswith("BoundAboveFallbackCap") for r in decision["Reasons"]))
        failed = self.decide(self.rule, "B4", "fallback", stop_record=STOP, approval="x", gates_passed=False)
        self.assertFalse(failed["Allowed"])
        self.assertIn("TransplantGateFailed", failed["Reasons"])

    def test_off_mode_is_not_a_reuse_mode(self):
        p = predictor.predict(self.rule, w=0.001, s=0.0, t2e_max=0.0, tail_donor=0.0)
        with self.assertRaises(predictor.NearKeyRuleError):
            predictor.policy_decision(self.rule, p, requested_mode="off", t2_passed=True, gates_passed=True, in_domain=True)


if __name__ == "__main__":
    unittest.main()
