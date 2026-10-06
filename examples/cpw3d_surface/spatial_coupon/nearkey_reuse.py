#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Near-key REUSE of a library coupon for a Missing SpatialEdgeCluster requirement (USER decisions 420 /
428 = Option A; supervisor decision 431; DESIGN v2 sections 1-4): detection -> transplant -> predictor ->
policy -> the reused model directory and its records. Python phase 1: the runtime C++ is unchanged (a
reused model carries the exact feature's Signature and basis, so the runtime matches and places it like
an exact coupon).

    nearkey_reuse.py evaluate   --library LIB.json --exact-basis DIR [--exact-signature-from NAME] --output DIR
                                [--mode default|fallback|measurement-only] [--donor NAME ...]
                                [--donor-matrices-root DIR] [--donor-shelled NAME=PATH ...] [--donor-record NAME=PATH ...]
                                [--stop-record PATH --approval TEXT] [--rule PATH]
    nearkey_reuse.py assemble   --library LIB.json --reused-model MODEL.json ... --output LIB.json --name NAME
                                [--note TEXT] [--measurement-only] [--keep-exact-names]

The model directory `<output>/<exact case>-reused/` = the generated basis files + the five transplanted
matrices (fabricated / thin domain + surface, the fabricated SHELLED surface) + interpolation-map-P.csv +
transplant-record.json; the model entry (reused-model.json) = the exact feature's generated entry with
the donor's response and the records of DESIGN 4.3: QualificationStatus ReusedResponse, ReuseMode,
ReusedFrom, NearKey, PredictedReuseError, TransplantTests, FallbackReason / FallbackStopRecord / Approval
or DefaultActivation, TailFromDonor, Note. `evaluate` exits 0 when a reuse was written, 2 when every
candidate was refused (the refusal recorded in <output>/nearkey-refusal.json), 1 on an error.
"""
import argparse
import copy
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import nearkey_detection as detection  # noqa: E402
import nearkey_predictor as predictor  # noqa: E402
import nearkey_transplant as transplant  # noqa: E402

STATUS_REUSED = detection.STATUS_REUSED
MODES = ("default", "fallback", "measurement-only")
MODE_MEASUREMENT = "measurement-only"
BASIS_FILES = ("trace-vertices.csv", "trace-triangles.csv", "basis-points.csv", "basis-contract.json", "zero-trace.csv")
REUSED_MODEL_RECORD = "reused-model.json"
REFUSAL_RECORD = "nearkey-refusal.json"


class NearKeyReuseError(ValueError):
    """A fail-closed stop of the reuse (the missing piece is named)."""


def sha256(path):
    return transplant.sha256(path)


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, text=True, stderr=subprocess.DEVNULL).strip()
    except (subprocess.CalledProcessError, OSError):
        return None


def donor_tail(model):
    """The donor's per-coupon shell tail t = SUM_s MA_sharp_s / SUM_s MA_raw_s - 1 (DESIGN 2.6) from its
    library `MA` record; fail closed without one."""
    ma = model.get("MA")
    if not isinstance(ma, dict) or not isinstance(ma.get("MA_raw"), dict) or not isinstance(ma.get("MA_sharp"), dict):
        raise NearKeyReuseError(f"donor {model.get('Name')}: no per-source MA shell record (MA_raw / MA_sharp): the MA_sharp bound "
                                "needs t_donor")
    raw, sharp = ma["MA_raw"], ma["MA_sharp"]
    total_raw = sum(float(v) for v in raw.values())
    if total_raw <= 0 or set(raw) != set(sharp):
        raise NearKeyReuseError(f"donor {model.get('Name')}: unusable MA shell record")
    return sum(float(v) for v in sharp.values()) / total_raw - 1.0


def donor_stored_state2(model, record_path=None):
    """The donor's (F) record state-2 MatrixIdentity Predicted energies per class, when the record is
    readable (the library's SpatialQualification.Record or an explicit local path); None otherwise."""
    path = record_path or (model.get("SpatialQualification") or {}).get("Record")
    if not path or not Path(path).is_file():
        return None, None
    record = json.loads(Path(path).read_text())
    trace = next((t for t in record.get("Traces", []) if t.get("Name") == "state-2"), None)
    if trace is None or "MatrixIdentity" not in trace:
        return None, sha256(path)
    predicted = {cls: value["Predicted"] for cls, value in trace["MatrixIdentity"].items()
                 if isinstance(value, dict) and "Predicted" in value}
    return predicted, sha256(path)


def load_exact_basis(basis_dir, name, model_entry):
    return transplant.TraceBasis(basis_dir, name, support_box=model_entry.get("SupportBox"), interfaces=model_entry.get("Interfaces"))


def load_donor(model, library_root, *, matrices_root=None, shelled_path=None, require_shelled=True):
    """The donor's basis + matrices. `matrices_root` = a local mirror holding <model directory name>/
    (the CSVs of the model's relative paths); `shelled_path` overrides the FabricatedSurfaceMatrixShelled
    path (a cluster path in the libraries of record). A missing shelled matrix is a fail-closed stop when
    `require_shelled` (Default / Fallback); the measurement-only demo records it as NotTransplanted."""
    paths = transplant.model_matrix_paths(model, library_root)
    if matrices_root is not None:
        directory = Path(matrices_root) / Path(model["FabricatedMatrix"]).parent.name
        for role, _, name in transplant.MATRIX_ROLES[:4]:
            paths[role] = str(directory / Path(model[{"fab_domain": "FabricatedMatrix", "fab_surface": "FabricatedSurfaceMatrix",
                                                      "thin_domain": "ThinMatrix", "thin_surface": "ThinSurfaceMatrix"}[role]]).name)
        basis_dir = directory
    else:
        basis_dir = Path(paths["fab_domain"]).parent
    if shelled_path is not None:
        paths["fab_surface_shelled"] = str(shelled_path)
    shelled_missing = paths["fab_surface_shelled"] is None or not Path(paths["fab_surface_shelled"]).is_file()
    if shelled_missing and require_shelled:
        raise NearKeyReuseError(f"donor {model['Name']}: the fabricated SHELLED surface matrix is not readable "
                                f"({paths['fab_surface_shelled']}): the MA_sharp transplant (decisions 422 / 431) cannot be made")
    donor = transplant.TraceBasis(basis_dir, model["Name"], support_box=model.get("SupportBox"), interfaces=model.get("Interfaces"))
    donor.load_matrices(paths)
    return donor, paths, shelled_missing


def choose_donor(evaluations):
    """DESIGN 1.3: among the admissible candidates the smallest max_T bound; ties by |W|, S, Name."""
    admissible = [e for e in evaluations if e["Admissible"]]
    if not admissible:
        return None
    return min(admissible, key=lambda e: (max(e["Prediction"][T]["Bound"] for T in predictor.TYPES), abs(e["NearKey"]["W"]),
                                           e["NearKey"]["S"], e["Donor"]))


def evaluate_candidate(exact, exact_signature, model, near_key, *, rule, library_root, mode, matrices_root=None, shelled_path=None,
                       record_path=None, stop_record=None, approval=None, radius):
    """One qualifying candidate: transplant + prediction + policy; returns the evaluation record (the
    transplant result kept under "Result" for the writer)."""
    require_shelled = mode != MODE_MEASUREMENT
    donor, paths, shelled_missing = load_donor(model, library_root, matrices_root=matrices_root, shelled_path=shelled_path,
                                               require_shelled=require_shelled)
    stored, record_sha = donor_stored_state2(model, record_path)
    result = transplant.transplant(exact, donor, radius, gates=rule["Policy"]["Gates"], domain_limits=rule["Policy"]["Domain"],
                                   donor_stored_state2=stored)
    tail = donor_tail(model)
    prediction = predictor.predict(rule, w=near_key["W"], s=near_key["S"], t2e_max=result["T2eMax"], tail_donor=tail)
    override = detection.donor_override_record(model)
    default_admits_donor = near_key["DefaultAdmissible"]
    in_domain = near_key["Qualifies"]
    decisions = {}
    for requested in ("default", "fallback"):
        reasons = []
        if requested == "default" and not default_admits_donor:
            reasons.append("DonorBuildGateOverride: " + (override["Kind"] if override else "donor status")
                           + " refuses default reuse (decision 431)")
        decision = predictor.policy_decision(rule, prediction, requested_mode=requested, t2_passed=result["T2Passed"],
                                             gates_passed=result["GatesPassed"], in_domain=in_domain,
                                             in_domain_reasons=near_key["Refusals"],
                                             stop_record=stop_record if requested == "fallback" else None,
                                             approval=approval if requested == "fallback" else None)
        if reasons:
            decision["Reasons"] = reasons + decision["Reasons"]
            decision["Allowed"], decision["Mode"] = False, predictor.MODE_NONE
        decisions[requested] = decision
    if mode == MODE_MEASUREMENT:
        # information: what the policy would say with a STOP record + approval present
        as_if = predictor.policy_decision(rule, prediction, requested_mode="fallback", t2_passed=result["T2Passed"],
                                          gates_passed=result["GatesPassed"], in_domain=in_domain, in_domain_reasons=near_key["Refusals"],
                                          stop_record={"Path": "(measurement-only: as if a STOP record existed)", "SHA256": "-"},
                                          approval="(measurement-only: as if approved)")
        decisions["FallbackIfStopRecorded"] = as_if
        admissible = result["GatesPassed"]
    else:
        admissible = decisions[mode]["Allowed"]
    return {"Donor": model["Name"], "DonorEntry": model, "DonorBasis": donor, "DonorPaths": paths, "ShelledMissing": shelled_missing,
            "DonorRecordSHA256": record_sha or (model.get("SpatialQualification") or {}).get("RecordSHA256"),
            "NearKey": near_key, "Result": result, "Prediction": prediction, "TailDonor": tail, "Override": override,
            "Decisions": decisions, "Admissible": admissible, "DefaultAdmitsDonor": default_admits_donor}


def candidate_summary(evaluation, chosen):
    summary = {"Model": evaluation["Donor"], "W": evaluation["NearKey"]["W"], "S": evaluation["NearKey"]["S"],
               "Bound": {T: evaluation["Prediction"][T]["Bound"] for T in predictor.TYPES},
               "MA_sharpBound": evaluation["Prediction"]["MA_sharp"]["Bound"], "T2eMax": evaluation["Result"]["T2eMax"],
               "GatesPassed": evaluation["Result"]["GatesPassed"], "T2Passed": evaluation["Result"]["T2Passed"],
               "Chosen": chosen}
    if not chosen:
        summary["RejectedReason"] = "not the smallest bound" if evaluation["Admissible"] else \
            "; ".join(r for d in evaluation["Decisions"].values() for r in d.get("Reasons", [])) or "transplant gates failed"
    return summary


def reuse_requirement(*, exact_signature, exact_basis_dir, exact_model_entry, library, library_path, rule, mode, output,
                      requirement_key=None, donors=None, matrices_root=None, shelled_paths=None, record_paths=None, stop_record=None,
                      approval=None, exact_interfaces=None, exact_boundary_condition=None, log=print, generator=None):
    """The reuse decision for one requirement; writes the model directory and returns the record
    {"Reused": bool, "Model": entry | None, "Refused": {...} | None, "Candidates": [...]}."""
    if mode not in MODES:
        raise NearKeyReuseError(f"--nearkey-reuse mode {mode!r} is not one of {MODES}")
    if mode == "default":
        active, _ = predictor.default_active(rule)
        if not active:
            raise NearKeyReuseError("--nearkey-reuse default: the rule file carries no DefaultActivation record (validation pair 5): "
                                    "default reuse is not active; use fallback (with a STOP record + approval) or off")
    if mode == "fallback":
        if not (isinstance(stop_record, dict) and stop_record.get("Path") and stop_record.get("SHA256")):
            raise NearKeyReuseError("--nearkey-reuse fallback needs the registration / build STOP record of the exact key")
        if not (isinstance(approval, str) and approval.strip()):
            raise NearKeyReuseError("--nearkey-reuse fallback needs --nearkey-fallback-approval TEXT")
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    library_root = Path(library_path).resolve().parent
    radius = float(library["MatchingRadius"])
    exact_name = exact_model_entry["Name"]
    key_hash = requirement_key or detection.signature_library.signature_hash(exact_signature)
    stamped = exact_model_entry.get("Signature")
    signature_hash = detection.signature_library.signature_hash
    if stamped is not None and signature_hash(stamped) != signature_hash(exact_signature):
        raise NearKeyReuseError(f"{exact_name}: the generated basis's Signature differs from the requirement's (hash mismatch)")
    records = detection.candidate_donors(exact_signature, library, rule=rule, exact_interfaces=exact_interfaces,
                                         exact_boundary_condition=exact_boundary_condition)
    record = {"Requirement": exact_name, "FeatureKey": key_hash, "Mode": mode, "RuleVersion": rule["RuleVersion"],
              "RuleFileSHA256": rule["_sha256"],
              "Candidates": [], "Reused": False, "Model": None, "Refused": None}
    if "*" in records:
        # the requirement's structure key is not calibrated: refused before any basis is read (DESIGN 1.2 item 0)
        record["Refused"] = {"Reason": records["*"]["Refusals"][0], "StructureKey": records["*"]["StructureKey"], "Candidates": [],
                             "Rule": "DESIGN 1.4: the feature stays Missing -> uncovered (F2 raw energy, decision 398)"}
        (output / REFUSAL_RECORD).write_text(json.dumps(record, indent=1) + "\n")
        log(f"{exact_name}: near-key reuse REFUSED: {record['Refused']['Reason']}")
        return record
    exact = load_exact_basis(exact_basis_dir, exact_name, exact_model_entry)
    by_name = {m["Name"]: m for m in library["Models"]}
    candidates, evaluations = [], []
    for name, near_key in records.items():
        if name == exact_name:
            continue
        if donors is not None and name not in donors:
            if near_key["Qualifies"]:
                candidates.append({"Model": name, "W": near_key.get("W"), "S": near_key.get("S"), "Chosen": False,
                                   "RejectedReason": "not in the --donor restriction (recorded)"})
            continue
        if not near_key["Qualifies"]:
            candidates.append({"Model": name, "W": near_key.get("W"), "S": near_key.get("S"), "Chosen": False,
                               "RejectedReason": "; ".join(near_key["Refusals"])})
            continue
        log(f"{exact_name}: candidate donor {name} (W {100 * near_key['W']:+.3f} %, S {100 * near_key['S']:.4f} %): transplant")
        try:
            evaluation = evaluate_candidate(exact, exact_signature, by_name[name], near_key, rule=rule,
                                            library_root=library_root, mode=mode,
                                            matrices_root=matrices_root, shelled_path=(shelled_paths or {}).get(name),
                                            record_path=(record_paths or {}).get(name), stop_record=stop_record, approval=approval,
                                            radius=radius)
        except (NearKeyReuseError, transplant.TransplantError) as error:
            # a candidate whose matrices / records cannot be read is not a donor (fail closed per candidate, recorded)
            candidates.append({"Model": name, "W": near_key.get("W"), "S": near_key.get("S"), "Chosen": False,
                               "RejectedReason": f"NotLoadable: {error}"})
            log(f"{exact_name}: candidate donor {name} rejected: {error}")
            continue
        evaluations.append(evaluation)
    chosen = choose_donor(evaluations)
    for evaluation in evaluations:
        candidates.append(candidate_summary(evaluation, evaluation is chosen))
    record.update({"Candidates": candidates, "Reused": chosen is not None})
    if chosen is None:
        reason = "no admissible candidate donor" if evaluations else "no candidate donor qualifies"
        record["Refused"] = {"Reason": reason, "Candidates": candidates,
                             "Rule": "DESIGN 1.4: the feature stays Missing -> uncovered (F2 raw energy, decision 398)"}
        (output / REFUSAL_RECORD).write_text(json.dumps(record, indent=1) + "\n")
        log(f"{exact_name}: near-key reuse REFUSED: {reason}")
        return record
    entry, model_dir = write_reused_model(exact=exact, exact_basis_dir=exact_basis_dir, exact_model_entry=exact_model_entry,
                                          evaluation=chosen,
                                         library=library, library_path=library_path, rule=rule, mode=mode, output=output, key_hash=key_hash,
                                         candidates=candidates, stop_record=stop_record, approval=approval, generator=generator)
    record["Model"] = entry
    record["ModelDirectory"] = str(model_dir)
    log(f"{exact_name}: near-key reuse {entry['ReuseMode']} from {chosen['Donor']} (W {100 * chosen['NearKey']['W']:+.3f} %; bounds "
        + " / ".join(f"{100 * chosen['Prediction'][T]['Bound']:.3f}" for T in predictor.TYPES)
        + f" %, MA_sharp {100 * chosen['Prediction']['MA_sharp']['Bound']:.3f} %)")
    return record


def write_reused_model(*, exact, exact_basis_dir, exact_model_entry, evaluation, library, library_path, rule, mode, output, key_hash,
                       candidates, stop_record, approval, generator):
    """The model directory and the reused model entry (DESIGN 4.1 item 3 / 4.3)."""
    exact_name = exact_model_entry["Name"]
    case = Path(exact_basis_dir).name
    model_dir = Path(output) / f"{case}-reused"
    if model_dir.exists():
        shutil.rmtree(model_dir)
    model_dir.mkdir(parents=True)
    copied = {}
    for name in BASIS_FILES:
        source = Path(exact_basis_dir) / name
        if source.is_file():
            shutil.copyfile(source, model_dir / name)
            copied[name] = sha256(model_dir / name)
    donor = evaluation["DonorBasis"]
    result = evaluation["Result"]
    written, map_record = transplant.write_matrices(model_dir, donor, result)
    tests = transplant.test_summary(result)
    near_key = dict(evaluation["NearKey"])
    near_key.update({"FeatureKey": key_hash, "Knots": {k: tests["Knots"][k] for k in ("Exact", "Donor", "Displaced",
                                                                                      "DonorOrphans", "ExactOrphans",
                                                                                      "MaxMatchedDisplacementOverR", "KnotsInterpolated")},
                     "KnotMaxOverR": tests["Knots"]["MaxMatchedDisplacementOverR"], "Candidates": candidates})
    prediction = dict(evaluation["Prediction"])
    reuse_mode = {"default": "Default", "fallback": "Fallback", MODE_MEASUREMENT: "MeasurementOnly"}[mode]
    prediction["Mode"] = reuse_mode
    prediction["PolicyLimit"] = {"Default": dict(rule["Policy"]["DefaultBound"]), "FallbackMax": rule["Policy"]["FallbackMaxBound"]}
    prediction["PolicyDecisions"] = evaluation["Decisions"]
    donor_entry = evaluation["DonorEntry"]
    transplant_record = {"Tool": "examples/cpw3d_surface/spatial_coupon/nearkey_transplant.py", "ToolSHA256": sha256(transplant.__file__),
                         "Commit": git_commit(), "Requirement": exact_name, "FeatureKey": key_hash, "Donor": donor_entry["Name"],
                         "Sizes": result["Sizes"], "Affine": result["Affine"], "Map": result["Map"], "Tests": result["Tests"],
                         "DonorMatrixFiles": {role: {"Path": path, "SHA256": sha256(path) if path and Path(path).is_file() else None}
                                              for role, path in evaluation["DonorPaths"].items()},
                         "Written": written, "MapCSV": map_record, "CopiedBasisFiles": copied}
    record_path = model_dir / "transplant-record.json"
    record_path.write_text(json.dumps(transplant_record, indent=1) + "\n")
    entry = copy.deepcopy(exact_model_entry)
    entry["Name"] = exact_name if mode == MODE_MEASUREMENT else f"{exact_name}-reused"
    for role, key, name in transplant.MATRIX_ROLES:
        if role in written:
            entry[key] = str(model_dir / name)
    if "fab_surface_shelled" not in written:
        entry["FabricatedSurfaceMatrixShelled"] = donor_entry.get("FabricatedSurfaceMatrixShelled")
        entry["ShelledMatrix"] = {"Status": "NotTransplanted: donor file unavailable (measurement-only marking, decision 431)",
                                  "DonorPath": donor_entry.get("FabricatedSurfaceMatrixShelled")}
    else:
        entry["ShelledMatrix"] = {"Status": "Transplanted", "Rule": "P^T Q_D,shelled P per radial shell (decisions 422 / 431)"}
    entry["MA"] = copy.deepcopy(donor_entry.get("MA"))
    entry["TailFromDonor"] = True
    entry["QualificationStatus"] = STATUS_REUSED
    entry["StatusProvisional"] = False
    entry["LibraryQualified"] = False
    entry["ReuseMode"] = reuse_mode
    entry["ReusedFrom"] = {"Donor": donor_entry["Name"],
                           "DonorKey": (detection.signature_library.signature_hash(donor_entry["Signature"])
                                        if donor_entry.get("Signature") else None),
                           "DonorLibrary": {"Name": library.get("Name"), "Version": library.get("Version"), "Path": str(library_path),
                                            "SHA256": sha256(library_path)},
                           "DonorRecordSHA256": evaluation["DonorRecordSHA256"],
                           "DonorMatrixSHA256": {role: value["SHA256"] for role, value in transplant_record["DonorMatrixFiles"].items()},
                           "MapCSVSHA256": map_record["SHA256"], "TransplantRecordSHA256": sha256(record_path),
                           "Tool": transplant_record["Tool"], "ToolSHA256": transplant_record["ToolSHA256"],
                           "BasisGenerator": generator or {"Command": "generate_spatial_response.py --basis-only (the registration's "
                                                                      "bound basis)",
                                                           "BasisFilesSHA256": copied},
                           "GeneratedSignatureHashEqualsKey": True, "DonorBuildGateOverride": evaluation["Override"]}
    entry["NearKey"] = near_key
    entry["PredictedReuseError"] = prediction
    entry["TransplantTests"] = tests
    if mode == "fallback":
        entry["FallbackReason"] = "the exact coupon cannot be built (STOP record recorded)"
        entry["FallbackStopRecord"] = stop_record
        entry["Approval"] = approval
    elif mode == "default":
        entry["DefaultActivation"] = rule["Policy"]["DefaultActivation"]
    else:
        entry["MeasurementOnly"] = True
        entry["ExactModelName"] = exact_name
    entry["Note"] = (f"REUSED RESPONSE ({reuse_mode}; {rule['RuleVersion']}; decisions 420 / 428 / 431): the response matrices of the "
                     f"donor "
                     f"{donor_entry['Name']} transplanted onto this feature's own generated basis (W {100 * near_key['W']:+.3f} %, S "
                     f"{100 * near_key['S']:.4f} %); predicted |r| bound SA / MS / MA "
                     + " / ".join(f"{100 * prediction[T]['Bound']:.3f}" for T in predictor.TYPES)
                     + f" %, MA_sharp {100 * prediction['MA_sharp']['Bound']:.3f} % (TailFromDonor). Never Qualified."
                     + (" MEASUREMENT-ONLY: never a library of record, never for any window verdict." if mode == MODE_MEASUREMENT else ""))
    (model_dir / REUSED_MODEL_RECORD).write_text(json.dumps(entry, indent=1) + "\n")
    return entry, model_dir


def library_header(rule, reused_models):
    """The library header block NearKeyReuse (DESIGN 4.3)."""
    counts = {"Default": 0, "Fallback": 0, "MeasurementOnly": 0}
    for model in reused_models:
        counts[model["ReuseMode"]] = counts.get(model["ReuseMode"], 0) + 1
    policy = rule["Policy"]
    return {"RuleVersion": rule["RuleVersion"], "RuleFileSHA256": rule["_sha256"],
            "AdmissibleStructureKeys": [entry["StructureKey"] for entry in rule["AdmissibleStructureKeys"]],
            "Policy": {"Option": policy["Option"], "DefaultBoundPct": {k: 100 * v for k, v in policy["DefaultBound"].items()},
                       "FallbackMaxPct": 100 * policy["FallbackMaxBound"],
                       "Domain": {"WMaxPct": 100 * policy["Domain"]["WMax"], "WDefaultMaxPct": 100 * policy["Domain"]["WDefaultMax"],
                                  "SMaxPct": 100 * policy["Domain"]["SMax"], "ShiftMaxOverR": policy["Domain"]["ShiftMaxOverR"]},
                       "DefaultActive": bool(policy.get("DefaultActive")), "DefaultActivation": policy.get("DefaultActivation")},
            "Models": [model["Name"] for model in reused_models], "Counts": counts}


def assemble_library(library, reused_models, *, rule, name, note=None, measurement_only=False, library_path=None):
    """The base library with the reused models added (or, under a measurement-only assembly keeping the
    exact names, replacing the exact models of the same Name) and the NearKeyReuse header."""
    out = copy.deepcopy(library)
    out["Name"] = name
    replaced = 0
    by_name = {m["Name"]: n for n, m in enumerate(out["Models"])}
    for model in reused_models:
        if model.get("QualificationStatus") != STATUS_REUSED:
            raise NearKeyReuseError(f"{model.get('Name')}: not a ReusedResponse model")
        if model["Name"] in by_name:
            if not measurement_only:
                raise NearKeyReuseError(f"{model['Name']}: a model of that Name is already in the library (only a measurement-only "
                                        f"assembly replaces the exact model)")
            out["Models"][by_name[model["Name"]]] = model
            replaced += 1
        else:
            out["Models"].append(model)
    out["NearKeyReuse"] = library_header(rule, reused_models)
    out["NearKeyReuse"]["Replaced"] = replaced
    if measurement_only:
        out["MeasurementOnly"] = True
        out["Note"] = ("MEASUREMENT-ONLY near-key reuse demo (decisions 428 / 431; nearkey-reuse-impl lane): NEVER a library of record, "
                       "never "
                       "for production or any window verdict. " + (note or "") + " " + str(library.get("Note", "")))
    elif note:
        out["Note"] = note + " " + str(library.get("Note", ""))
    if library_path is not None:
        out["DerivedFrom"] = {"Path": str(library_path), "SHA256": sha256(library_path), "Name": library.get("Name"),
                              "Rule": "near-key reused models added by nearkey_reuse.py assemble (DESIGN v2 section 4)"}
    return out


def parse_assignments(values, option):
    out = {}
    for value in values or ():
        name, separator, path = str(value).partition("=")
        if not separator or not name or not path:
            raise NearKeyReuseError(f"{option} expects NAME=PATH, not {value!r}")
        out[name] = path
    return out


def stop_record_from_path(path, requirement_name):
    """The fallback STOP record of the exact key: an existing JSON file naming the requirement (its Name or
    its case id) in its text; recorded with its sha256 (fail closed otherwise)."""
    path = Path(path)
    if not path.is_file():
        raise NearKeyReuseError(f"STOP record {path} does not exist")
    text = path.read_text()
    try:
        json.loads(text)
    except json.JSONDecodeError as error:
        raise NearKeyReuseError(f"STOP record {path} is not JSON: {error}") from error
    if requirement_name.split("_")[-1] not in text:
        raise NearKeyReuseError(f"STOP record {path} does not name the requirement {requirement_name}")
    return {"Path": str(path), "SHA256": sha256(path)}


def command_evaluate(args):
    rule = predictor.load_rule(args.rule)
    library = json.loads(Path(args.library).read_text())
    basis_dir = Path(args.exact_basis)
    entry = None
    if (basis_dir / "process-library.json").is_file():
        generated = json.loads((basis_dir / "process-library.json").read_text())
        if len(generated.get("Models", [])) != 1:
            raise NearKeyReuseError(f"{basis_dir}: the basis directory's process-library.json carries "
                                    f"{len(generated.get('Models', []))} models, not one")
        entry = generated["Models"][0]
    elif not args.exact_signature_from:
        raise NearKeyReuseError(f"{basis_dir}: no process-library.json (the generated model entry) and no --exact-signature-from")
    if args.exact_signature_from:
        source = next((m for m in library["Models"] if m["Name"] == args.exact_signature_from), None)
        if source is None:
            raise NearKeyReuseError(f"--exact-signature-from {args.exact_signature_from}: no such model in the library")
        signature = source["Signature"]
        entry = copy.deepcopy(source)
        for key in ("BasisPoints", "TraceMesh"):
            entry.pop(key, None)
        entry.update({"BasisPoints": "basis-points.csv", "TraceMesh": {"Vertices": "trace-vertices.csv",
                                                                       "Triangles": "trace-triangles.csv"}})
        for key in ("FabricatedSurfaceMatrixShelled", "CouponMesh", "ThinCouponMesh", "Qualification", "ThinQualification",
                    "SpatialQualification", "CombinedFrom", "SourceProcessLibrary", "BuildGateOverride", "StatusProvisionalTool"):
            entry.pop(key, None)
    else:
        signature = entry.get("Signature")
        if signature is None:
            raise NearKeyReuseError(f"{basis_dir}: the generated model carries no Signature")
    stop = stop_record_from_path(args.stop_record, entry["Name"]) if args.stop_record else None
    record = reuse_requirement(exact_signature=signature, exact_basis_dir=basis_dir, exact_model_entry=entry, library=library,
                               library_path=args.library, rule=rule, mode=args.mode, output=args.output, donors=args.donor or None,
                               matrices_root=args.donor_matrices_root, shelled_paths=parse_assignments(args.donor_shelled,
                                                                                                       "--donor-shelled"),
                               record_paths=parse_assignments(args.donor_record,
                                                              "--donor-record"), stop_record=stop, approval=args.approval)
    (Path(args.output) / f"{entry['Name']}-nearkey-reuse.json").write_text(json.dumps(record, indent=1) + "\n")
    return 0 if record["Reused"] else 2


def command_assemble(args):
    rule = predictor.load_rule(args.rule)
    library = json.loads(Path(args.library).read_text())
    reused = [json.loads(Path(path).read_text()) for path in args.reused_model]
    out = assemble_library(library, reused, rule=rule, name=args.name, note=args.note, measurement_only=args.measurement_only,
                           library_path=args.library)
    Path(args.output).write_text(json.dumps(out, indent=1) + "\n")
    print(f"library {args.output} {sha256(args.output)} models {len(out['Models'])} reused {len(reused)} "
          f"({'MEASUREMENT-ONLY' if args.measurement_only else 'of record candidate'})")
    return 0


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    commands = parser.add_subparsers(dest="command", required=True)
    evaluate = commands.add_parser("evaluate", help="detect, transplant, predict and decide one requirement")
    evaluate.add_argument("--library", required=True,
                          help="the process library (the donor pool; its Models' paths resolve against its directory)")
    evaluate.add_argument("--exact-basis", required=True,
                          help="the requirement's --basis-only directory (trace mesh, basis points, process-library.json)")
    evaluate.add_argument("--exact-signature-from",
                          help="take the requirement's Signature / geometric entry from this library model (the demo: "
                                                         "the exact model exists in the library; its own entry is refused as a quantum "
                                                         "match)")
    evaluate.add_argument("--output", required=True)
    evaluate.add_argument("--mode", choices=MODES, default="fallback")
    evaluate.add_argument("--donor", action="append", help="restrict the candidate donors to these model Names (recorded)")
    evaluate.add_argument("--donor-matrices-root", help="local mirror of the donors' model directories (<model dir name>/<csv>)")
    evaluate.add_argument("--donor-shelled", action="append", metavar="NAME=PATH",
                          help="the donor's fabricated SHELLED surface matrix file")
    evaluate.add_argument("--donor-record", action="append", metavar="NAME=PATH",
                          help="the donor's (F) spatial-qualification.json (T4 vs the stored energies)")
    evaluate.add_argument("--stop-record", help="fallback: the registration / build STOP record of the exact key (JSON)")
    evaluate.add_argument("--approval", help="fallback: the recorded approval text")
    evaluate.add_argument("--rule", default=str(predictor.RULE_FILE))
    evaluate.set_defaults(func=command_evaluate)
    assemble = commands.add_parser("assemble", help="add reused models to a library with the NearKeyReuse header")
    assemble.add_argument("--library", required=True)
    assemble.add_argument("--reused-model", action="append", required=True, help="a reused-model.json written by evaluate (repeatable)")
    assemble.add_argument("--output", required=True)
    assemble.add_argument("--name", required=True)
    assemble.add_argument("--note")
    assemble.add_argument("--measurement-only", action="store_true",
                          help="a measurement-only library (the exact models of the same Name replaced)")
    assemble.add_argument("--rule", default=str(predictor.RULE_FILE))
    assemble.set_defaults(func=command_assemble)
    return parser


def main(argv=None):
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except (NearKeyReuseError, predictor.NearKeyRuleError, detection.NearKeyDetectionError, transplant.TransplantError) as error:
        print(f"NEARKEY_REUSE_FAILED: {error}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
