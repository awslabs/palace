#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""The near-key reuse PREDICTOR and POLICY of RuleVersion nearkey-reuse-rule-v1 (USER decisions
420 / 428 = Option A; DESIGN v2 sections 2, 2.6 and 3). Pure arithmetic on the recorded inputs of
a reused coupon; the coefficients, regime boundary, policy limits and admissible structure keys
are read from the pinned rule file `nearkey-reuse-rule-v1.json` (its sha256 is recorded on every
prediction so that a later RuleVersion can re-score a stored model from its recorded W / S / T2e).

    central   r_c,T    = -a_T W                    (W = (w_donor - w_exact) / w_exact of the narrowest claim strip;
                                                     a narrower donor reads MORE energy)
    bound     |r|_T    = (1 + rho_T) a_T |W| + b_T S + phi_T + 2 max_X T2e_X
    regimes   DEFAULT  |W| <= WDefaultMax (2.5 %): the group-B per-pair-ratio calibration;
              FALLBACK |W| > WDefaultMax: max(the full-range form, the default form at WDefaultMax) - monotone
    MA_sharp  Bound_MA x (1 + t_donor) + TailDifferenceBound                                  (decision 422)
    policy    DEFAULT  iff bound <= 0.5 % SA, <= 0.5 % MS, <= 1.0 % MA_sharp AND |W| <= WDefaultMax AND the rule's
                       DefaultActive is true (validation pair 5 recorded) AND the T2 trace gate passed (<= 1e-3);
              FALLBACK iff every Type's bound <= 3.5 % inside the calibrated domain, the exact coupon is unbuildable
                       (a STOP record) and an approval is recorded; else NO reuse (the feature stays Missing / F2).
"""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
RULE_FILE = HERE / "nearkey-reuse-rule-v1.json"
RULE_VERSION = "nearkey-reuse-rule-v1"
TYPES = ("SA", "MS", "MA")
REGIME_DEFAULT = "Default"
REGIME_FALLBACK = "Fallback"
MODE_DEFAULT = "Default"
MODE_FALLBACK = "Fallback"
MODE_NONE = "None"


class NearKeyRuleError(ValueError):
    """The rule file or a prediction input is not usable (fail closed)."""


def load_rule(path=RULE_FILE):
    """The pinned rule with its sha256 (recorded as PredictedReuseError.RuleFileSHA256)."""
    path = Path(path)
    payload = path.read_bytes()
    rule = json.loads(payload)
    if rule.get("RuleVersion") != RULE_VERSION:
        raise NearKeyRuleError(f"{path}: RuleVersion {rule.get('RuleVersion')!r} is not {RULE_VERSION}")
    for T in TYPES:
        for regime in (REGIME_DEFAULT, REGIME_FALLBACK):
            coefficients = rule["Predictor"]["Coefficients"][T][regime]
            for name in ("a", "rho", "b", "phi"):
                value = coefficients[name]
                if not isinstance(value, (int, float)) or value < 0:
                    raise NearKeyRuleError(f"{path}: coefficient {T}.{regime}.{name} = {value!r} is not a non-negative number")
    keys = [entry["StructureKey"] for entry in rule["AdmissibleStructureKeys"]]
    if len(keys) != 4 or len(set(keys)) != 4 or any(len(k) != 64 for k in keys):
        raise NearKeyRuleError(f"{path}: RuleVersion v1 pins exactly four 64-hex structure keys, found {keys}")
    rule["_sha256"] = hashlib.sha256(payload).hexdigest()
    rule["_path"] = str(path)
    return rule


def admissible_structure_keys(rule):
    return {entry["StructureKey"]: entry for entry in rule["AdmissibleStructureKeys"]}


def _bound_form(coefficients, w_abs, s, transplant):
    terms = {"Width": (1.0 + coefficients["rho"]) * coefficients["a"] * w_abs, "Vertex": coefficients["b"] * s,
             "Floor": coefficients["phi"], "Transplant": transplant}
    return sum(terms.values()), terms


def predict_type(rule, T, w, s, t2e_max):
    """One Type's prediction: {Regime, Central, Bound, Terms, Coefficients}. The regime is set by |W|
    against WDefaultMax; beyond it the bound is max(fallback form, default form at WDefaultMax)."""
    if T not in TYPES:
        raise NearKeyRuleError(f"unknown Type {T!r}")
    for name, value in (("W", w), ("S", s), ("T2e", t2e_max)):
        if not isinstance(value, (int, float)) or value != value:
            raise NearKeyRuleError(f"{name} = {value!r} is not a finite number")
    if s < 0 or t2e_max < 0:
        raise NearKeyRuleError(f"S and T2e must be non-negative (S {s}, T2e {t2e_max})")
    w_default_max = float(rule["Predictor"]["WDefaultMax"])
    transplant = 2.0 * float(t2e_max)
    coefficients = rule["Predictor"]["Coefficients"][T]
    if abs(w) <= w_default_max:
        regime = REGIME_DEFAULT
        used = coefficients[REGIME_DEFAULT]
        bound, terms = _bound_form(used, abs(w), s, transplant)
        governing = REGIME_DEFAULT
    else:
        regime = REGIME_FALLBACK
        used = coefficients[REGIME_FALLBACK]
        bound_fallback, terms_fallback = _bound_form(used, abs(w), s, transplant)
        bound_boundary, terms_boundary = _bound_form(coefficients[REGIME_DEFAULT], w_default_max, s, transplant)
        if bound_fallback >= bound_boundary:
            bound, terms, governing = bound_fallback, terms_fallback, REGIME_FALLBACK
        else:
            bound, terms, governing = bound_boundary, terms_boundary, "DefaultAtWDefaultMax"
    return {"Regime": regime, "GoverningForm": governing, "Central": -used["a"] * w, "Bound": bound, "Terms": terms,
            "Coefficients": {name: used[name] for name in ("a", "rho", "b", "phi")}}


def ma_sharp_bound(rule, bound_ma, tail_donor):
    """DESIGN 2.6: Bound_MA x (1 + t_donor) + TailDifferenceBound (absolute, on the MA energy)."""
    if not isinstance(tail_donor, (int, float)) or tail_donor != tail_donor or tail_donor < 0:
        raise NearKeyRuleError(f"the donor's MA tail t_donor = {tail_donor!r} must be a non-negative number")
    tail_bound = float(rule["Predictor"]["MASharp"]["TailDifferenceBound"])
    return bound_ma * (1.0 + tail_donor) + tail_bound, tail_bound


def predict(rule, *, w, s, t2e_max, tail_donor):
    """The PredictedReuseError block without the Mode / PolicyLimit (policy_decision adds them)."""
    types = {T: predict_type(rule, T, w, s, t2e_max) for T in TYPES}
    regimes = {types[T]["Regime"] for T in TYPES}
    assert len(regimes) == 1, regimes
    sharp, tail_bound = ma_sharp_bound(rule, types["MA"]["Bound"], tail_donor)
    return {"RuleVersion": rule["RuleVersion"], "RuleFileSHA256": rule["_sha256"],
            "CalibrationRecordSHA256": rule["CalibrationRecord"]["SHA256"], "Regime": regimes.pop(),
            "Inputs": {"W": w, "S": s, "T2eMax": t2e_max, "TailDonor": tail_donor},
            "Coefficients": {T: types[T]["Coefficients"] for T in TYPES},
            **{T: {key: types[T][key] for key in ("Central", "Bound", "Terms", "GoverningForm")} for T in TYPES},
            "MA_sharp": {"Bound": sharp, "TailDonor": tail_donor, "TailDifferenceBound": tail_bound, "TailFromDonor": True}}


def default_active(rule):
    """The rule's DefaultActive flag is honoured only with a DefaultActivation {ValidationPair, RecordPath,
    RecordSHA256, Decision} record whose RecordPath is readable and re-hashes to RecordSHA256 (decision 438 (5):
    a missing / unreadable / mismatching validation-pair record fails closed)."""
    policy = rule["Policy"]
    activation = policy.get("DefaultActivation")
    if not policy.get("DefaultActive"):
        return False, None
    if not (isinstance(activation, dict) and activation.get("ValidationPair") and activation.get("RecordSHA256")
            and activation.get("RecordPath")):
        raise NearKeyRuleError("the rule says DefaultActive true without a DefaultActivation {ValidationPair, RecordPath, RecordSHA256} "
                               "record")
    record = Path(activation["RecordPath"])
    if not record.is_file():
        raise NearKeyRuleError(f"DefaultActivation.RecordPath {record} is not readable: the validation-pair record cannot be verified")
    actual = hashlib.sha256(record.read_bytes()).hexdigest()
    if actual != activation["RecordSHA256"]:
        raise NearKeyRuleError(f"DefaultActivation.RecordSHA256 {activation['RecordSHA256'][:16]}… != the validation-pair record's "
                               f"{actual[:16]}… ({record}): default reuse stays inactive")
    return True, activation


def is_shipped_rule(rule):
    """True iff the rule was loaded from the repository's pinned rule file (decision 438 (5): only the shipped,
    test-pinned rule activates default reuse; a --nearkey-rule override is refused in Default mode)."""
    return Path(rule.get("_path", "")).resolve() == RULE_FILE.resolve()


def policy_decision(rule, prediction, *, requested_mode, t2_passed, gates_passed, in_domain, in_domain_reasons=(),
                    stop_record=None, approval=None):
    """The per-coupon policy (DESIGN section 3, Option A) -> {Mode, Allowed, Reasons, PolicyLimit,
    DefaultActivation}. `requested_mode` is the build's --nearkey-reuse setting ("default" or
    "fallback"); `in_domain` is the detection's verdict (every item of DESIGN 1.2 inside its limit);
    `stop_record` / `approval` are the fallback's STOP record of the exact key and the approval text."""
    policy = rule["Policy"]
    limits = policy["DefaultBound"]
    reasons = []
    bounds = {T: prediction[T]["Bound"] for T in TYPES}
    sharp = prediction["MA_sharp"]["Bound"]
    result = {"PolicyLimit": {"Default": dict(limits), "FallbackMaxPct": 100.0 * policy["FallbackMaxBound"]},
              "RequestedMode": requested_mode, "DefaultActivation": None, "Mode": MODE_NONE, "Allowed": False}
    if requested_mode not in (MODE_DEFAULT.lower(), MODE_FALLBACK.lower()):
        raise NearKeyRuleError(f"--nearkey-reuse must be 'default' or 'fallback' to reuse anything, not {requested_mode!r}")
    if not gates_passed:
        reasons.append("TransplantGateFailed")
    if not in_domain:
        reasons.append("OutsideCalibratedDomain" + (f": {'; '.join(in_domain_reasons)}" if in_domain_reasons else ""))
    if requested_mode == MODE_DEFAULT.lower():
        active, activation = default_active(rule)
        if not active:
            reasons.append("DefaultNotActive: no DefaultActivation record (validation pair 5 not yet measured inside its bound)")
        result["DefaultActivation"] = activation
        if prediction["Regime"] != REGIME_DEFAULT:
            reasons.append(f"OutsideDefaultRegime: |W| {abs(prediction['Inputs']['W']):.5f} > {rule['Predictor']['WDefaultMax']}")
        if not t2_passed:
            reasons.append("T2TraceGateFailed: default reuse needs T2 <= the trace tolerance")
        over = [f"{T} {100 * bounds[T]:.4f} % > {100 * limits[T]:.2f} %" for T in ("SA", "MS") if bounds[T] > limits[T]]
        if sharp > limits["MA_sharp"]:
            over.append(f"MA_sharp {100 * sharp:.4f} % > {100 * limits['MA_sharp']:.2f} %")
        if over:
            reasons.append("BoundAbovePolicyLimit: " + "; ".join(over))
        if not reasons:
            result.update({"Mode": MODE_DEFAULT, "Allowed": True})
    else:
        if not (isinstance(stop_record, dict) and stop_record.get("Path") and stop_record.get("SHA256")):
            reasons.append("NoStopRecord: fallback reuse needs the registration / build STOP record of the exact key")
        if not (isinstance(approval, str) and approval.strip()):
            reasons.append("NoApproval: fallback reuse needs a recorded approval")
        cap = policy["FallbackMaxBound"]
        over = [f"{T} {100 * bounds[T]:.4f} %" for T in TYPES if bounds[T] > cap]
        if sharp > cap:
            over.append(f"MA_sharp {100 * sharp:.4f} %")
        if over:
            reasons.append(f"BoundAboveFallbackCap {100 * cap:.1f} %: " + "; ".join(over))
        if not reasons:
            result.update({"Mode": MODE_FALLBACK, "Allowed": True})
    result["Reasons"] = reasons
    return result


def implied_default_tolerance(rule, s_over_w=0.10, tail_donor=0.0, limits=None):
    """Information: the |W| at which each Type's bound (default form, S = s_over_w |W|, no transplant
    term, MA in the MA_sharp form with the given donor tail) meets the policy limit; the binding Type
    is the minimum (DESIGN 3: SA-binding at 1.78 % -> |W| <= 1.7 % under Option A)."""
    limits = limits or rule["Policy"]["DefaultBound"]
    out = {}
    for T in TYPES:
        c = rule["Predictor"]["Coefficients"][T][REGIME_DEFAULT]
        slope = (1.0 + c["rho"]) * c["a"] + c["b"] * s_over_w
        if T == "MA":
            tail_bound = rule["Predictor"]["MASharp"]["TailDifferenceBound"]
            out["MA_sharp"] = ((limits["MA_sharp"] - tail_bound) / (1.0 + tail_donor) - c["phi"]) / slope
        else:
            out[T] = (limits[T] - c["phi"]) / slope
    out["Binding"] = min((T for T in out), key=out.get)
    return out
