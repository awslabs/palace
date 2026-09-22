#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Register a coupon source directory as a Gmsh-only manifest case (supervisor
decision 48, `coupon-library build` step 1): the only hand-edited JSON of the mesh
path becomes a derived record.

A source directory carries the producer's frozen inputs - mesh-signature.csv,
plan-view-boundary.csv, plan-view-mask.csv, process.toml, optionally
process-library.json, the trace basis (basis-contract.json, trace-vertices.csv,
trace-triangles.csv; all three, and only with a process library), provenance.json and
retained-etch.csv - and either its own mesh-recipe.json or a --mesh-recipe repository
path.  The etch footprint is an explicit declaration (--footprint bound: the
directory's retained-etch.csv is frozen; --footprint producer-default: no such file
exists); registration fails closed when it is missing or contradicts the directory.

Registration:

1. SHA-256 of every source file; the recipe scope classes exhibited by the inputs
   (mesh_stage_contract.scope_classes_of_case_inputs) - a class the recipe guards
   stops the registration as "unsupported-class" with the guard id, no probe built;
2. the two-pass contract derivation of derive_semantic_contract.py, orchestrated:
   a PROVISIONAL contract from the inputs, a labels-only probe of a staging copy
   (run_gmsh_only_case.py --labels-only under the production recipe, decision 62(2):
   headroom gate, source validation, the mesher with the production gmsh-build options
   stopping right after the CAD-entity labelling - the label set is the only thing the
   derivation learns from a build and it exists before any mesh is generated; a mesher
   ScopeGuard stop is recorded as unsupported-class), then the final contract from the
   probe's label census (--build-census), written to SOURCE_DIR/semantic-contract.json;
   the production build that follows fails closed when its own census labels differ
   from the contract (mesh_stage_contract.validate_gmsh_build_census);
3. the manifest case: the source file digests, the production Variants (identity and
   rotate-z) / TransformComparison / SignatureColumns shared by every existing case
   (fail closed when they differ), InventoryStatus, Features (default: the exhibited
   scope classes), FixtureVersion and a Provenance statement naming the commit, the
   footprint declaration and the probe root;
4. refreeze_manifest_tools.refreeze over the written manifest (the calibration
   mirrors too when it is the production manifest).

Idempotent by content: a case whose recorded source digests equal the directory's is
reused (Status "reused", nothing written); a case with any changed source becomes the
next FixtureVersion and the previous entry's changed bindings are kept under
RetiredFixtures.Entries (never a silent edit).  The outcome is recorded in
WORK/register-case.json.

usage: register_case.py CASE_ID SOURCE_DIR --footprint {bound,producer-default}
       --inventory-status STATUS [--manifest PATH] [--mesh-recipe REPOSITORY_PATH]
       [--feature NAME ...] [--provenance TEXT] [--work DIR] [--julia PATH]
"""
import argparse
import copy
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from derive_semantic_contract import derive as derive_contract  # noqa: E402
from general_mesh_manifest import (CASE_KIND_KEY, CASE_KINDS, FABRICATED_CASE_KEY,  # noqa: E402
                                   PRODUCER_DEFAULT_ETCH_FOOTPRINT, TRACE_BASIS_ROLES, case_kind,
                                   preflight_recipe_scope, validate_manifest, validate_thin_recipe)
from refreeze_manifest_tools import CALIBRATION_MANIFESTS, PRODUCTION_MANIFEST, refreeze  # noqa: E402
from run_gmsh_only_case import BUILD_SUMMARY  # noqa: E402

REQUIRED_SOURCE_FILES = {"Signature": "mesh-signature.csv", "Boundary": "plan-view-boundary.csv",
                         "Mask": "plan-view-mask.csv", "Process": "process.toml"}
OPTIONAL_SOURCE_FILES = {"ProcessLibrary": "process-library.json", "BasisContract": "basis-contract.json",
                         "TraceVertices": "trace-vertices.csv", "TraceTriangles": "trace-triangles.csv",
                         "Provenance": "provenance.json"}
RETAINED_ETCH_FILE = "retained-etch.csv"
SEMANTIC_CONTRACT_FILE = "semantic-contract.json"
# The thin case of a source directory (decision 66) shares every frozen source of its
# fabricated case and derives its own contract (the thin label families) into this file.
THIN_SEMANTIC_CONTRACT_FILE = "semantic-contract-thin.json"
THIN_CASE_SUFFIX = "-thin"
CONTRACT_FILES = {"fabricated": SEMANTIC_CONTRACT_FILE, "thin": THIN_SEMANTIC_CONTRACT_FILE}


def thin_case_id(fabricated_case_id):
    return fabricated_case_id + THIN_CASE_SUFFIX
MESH_RECIPE_FILE = "mesh-recipe.json"
FOOTPRINT_BOUND = "bound"
FOOTPRINT_DECLARATIONS = (FOOTPRINT_BOUND, PRODUCER_DEFAULT_ETCH_FOOTPRINT)
# The production manifest's inventory vocabulary; "Calibration" belongs to the labeled
# calibration manifests only and is never registered here.  "DeviceDerived" is a case
# whose sources the device adapter (device_coupons.py) produced from a device's
# discovery closure (decision 52).
INVENTORY_STATUSES = ("RepositoryFixture", "RepositoryAssessmentFixture", "RemoteVerified",
                      "LocalImmutableCalibration", "DeviceDerived")
# Case fields every production case shares (the production placement variants, the
# physical covariance comparison, the signature preflight columns).
SHARED_CASE_FIELDS = ("Variants", "TransformComparison")
SHARED_SOURCE_FIELDS = ("SignatureColumns",)
RETIRED_FIXTURES_RULE = ("superseded fixture versions are never edited silently: every re-bound or "
                         "re-derived frozen input keeps its old SHA256 here with the evidence that "
                         "retired it; a retired version is not buildable and has no evidentiary value "
                         "(supervisor decisions 46 and 48)")
REGISTRATION_RECORD = "register-case.json"
STATUS_REGISTERED, STATUS_REUSED, STATUS_UNSUPPORTED, STATUS_FAILED = (
    "registered", "reused", "unsupported-class", "failed")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


class RegistrationError(ValueError):
    """A fail-closed registration outcome (recorded, never a manifest edit)."""


def source_files(directory, footprint):
    """Role -> path of the source files the directory binds, validated against the
    footprint declaration and the all-or-none trace basis rule."""
    directory = Path(directory)
    if footprint not in FOOTPRINT_DECLARATIONS:
        raise RegistrationError(f"the etch footprint must be declared as one of "
                                f"{list(FOOTPRINT_DECLARATIONS)}, not {footprint!r}")
    paths = {}
    for role, name in REQUIRED_SOURCE_FILES.items():
        path = directory / name
        if not path.is_file():
            raise RegistrationError(f"source directory lacks {name} ({role})")
        paths[role] = path
    for role, name in OPTIONAL_SOURCE_FILES.items():
        if (directory / name).is_file():
            paths[role] = directory / name
    basis = {role for role in TRACE_BASIS_ROLES if role in paths}
    if basis and basis != set(TRACE_BASIS_ROLES):
        raise RegistrationError(f"a trace basis binds all of {list(TRACE_BASIS_ROLES)} (the build passes "
                                f"the process library with it), found {sorted(basis)}")
    retained = directory / RETAINED_ETCH_FILE
    if footprint == FOOTPRINT_BOUND:
        if not retained.is_file():
            raise RegistrationError(f"--footprint {FOOTPRINT_BOUND} declared but {RETAINED_ETCH_FILE} "
                                    f"is absent")
        paths["RetainedEtch"] = retained
    elif retained.is_file():
        raise RegistrationError(f"--footprint {PRODUCER_DEFAULT_ETCH_FOOTPRINT} declared but "
                                f"{RETAINED_ETCH_FILE} exists: declare --footprint {FOOTPRINT_BOUND} "
                                f"or remove the file")
    return paths


def mesh_recipe_binding(directory, repository, mesh_recipe):
    """The frozen mesh recipe: the directory's mesh-recipe.json (Name binding) or an
    explicit repository path (RepositoryPath binding); exactly one."""
    local = Path(directory) / MESH_RECIPE_FILE
    if mesh_recipe is None:
        if not local.is_file():
            raise RegistrationError(f"no {MESH_RECIPE_FILE} in the source directory and no --mesh-recipe")
        return {"Name": MESH_RECIPE_FILE, "SHA256": sha256(local)}, local
    if local.is_file():
        raise RegistrationError(f"--mesh-recipe given but the source directory carries {MESH_RECIPE_FILE}")
    path = Path(mesh_recipe)
    path = path if path.is_absolute() else repository / path
    if not path.is_file():
        raise RegistrationError(f"mesh recipe {path} does not exist")
    return {"RepositoryPath": os.path.relpath(path.resolve(), repository), "SHA256": sha256(path)}, path


def shared_case_fields(manifest):
    """The fields every existing production case shares; fail closed when they differ,
    so a registered case never introduces a second placement or comparison rule."""
    cases = manifest["Cases"]
    if not cases:
        raise RegistrationError("the manifest declares no case to take the shared fields from")
    shared = {}
    for field in SHARED_CASE_FIELDS:
        values = {json.dumps(case[field], sort_keys=True) for case in cases}
        if len(values) != 1:
            raise RegistrationError(f"existing cases disagree on {field}: cannot register")
        shared[field] = copy.deepcopy(cases[0][field])
    for field in SHARED_SOURCE_FIELDS:
        values = {json.dumps(case["Source"][field], sort_keys=True) for case in cases}
        if len(values) != 1:
            raise RegistrationError(f"existing cases disagree on Source.{field}: cannot register")
        shared[field] = copy.deepcopy(cases[0]["Source"][field])
    return shared


def recorded_directory(directory, repository):
    """Source.Directory as the manifest records it: repository-relative inside the
    repository, absolute elsewhere."""
    directory = Path(directory).resolve()
    try:
        return str(directory.relative_to(repository))
    except ValueError:
        return str(directory)


def file_entries(paths, recipe_entry):
    entries = {role: {"Name": path.name, "SHA256": sha256(path)} for role, path in paths.items()}
    entries["MeshRecipe"] = recipe_entry
    return entries


def source_digests(entries):
    """The content identity of a case's sources: every role but the derived contract."""
    return {role: entry["SHA256"] for role, entry in entries.items() if role != "SemanticContract"}


def case_directory(case, repository):
    directory = Path(case["Source"]["Directory"])
    return directory if directory.is_absolute() else repository / directory


PROBE_MODE = "--labels-only"


def run_probe_build(python, julia, manifest_path, case_id, root, log):
    """The labels-only probe of the staging case (decision 62(2)); returns its build summary."""
    command = [python, str(HERE / "run_gmsh_only_case.py"), case_id, "--manifest", str(manifest_path),
               "--root", str(root), PROBE_MODE]
    if julia is not None:
        command += ["--julia", str(julia)]
    with open(log, "w") as stream:
        result = subprocess.run(command, stdout=stream, stderr=subprocess.STDOUT)
    summary = root / BUILD_SUMMARY
    if summary.is_file():
        return json.loads(summary.read_text())
    return {"Status": STATUS_FAILED, "Stage": "probe-launch", "ReturnCode": result.returncode,
            "ScopeGuard": None, "Message": f"probe build wrote no {BUILD_SUMMARY}; see {log}"}


def derive_two_pass(manifest, repository, case, paths, recipe_path, work, probe):
    """Provisional contract -> labels-only probe of a staging copy -> final contract
    from the probe's label census.  Returns (contract, probe root, probe summary);
    raises RegistrationError with the probe summary attached when the probe stopped."""
    staging = work / "probe-source"
    staging.mkdir(parents=True)
    staged = {}
    for role, path in paths.items():
        shutil.copyfile(path, staging / path.name)
        staged[role] = staging / path.name
    process_library = staged.get("ProcessLibrary")
    kind = case_kind(case)
    provisional = derive_contract(staging, None, signature=staged["Signature"],
                                  boundary=staged["Boundary"], process_library=process_library, kind=kind)
    (staging / CONTRACT_FILES[kind]).write_text(json.dumps(provisional, indent=2) + "\n")
    staged["SemanticContract"] = staging / CONTRACT_FILES[kind]
    recipe_entry = dict(case["Source"]["Files"]["MeshRecipe"])
    if "Name" in recipe_entry:
        shutil.copyfile(recipe_path, staging / recipe_entry["Name"])
    probe_case = copy.deepcopy(case)
    probe_case["Source"]["Directory"] = str(staging)
    probe_case["Source"]["Files"] = file_entries(staged, recipe_entry)
    probe_manifest = copy.deepcopy(manifest)
    probe_manifest["RepositoryRoot"] = str(repository)
    probe_manifest["Cases"] = [item for item in manifest["Cases"] if item["Id"] != case["Id"]] + [probe_case]
    if kind == "thin":
        # The probe's fabricated case shares the staging sources (validate_case_kind).
        fabricated = copy.deepcopy(next(item for item in probe_manifest["Cases"] if item["Id"] == case[FABRICATED_CASE_KEY]))
        fabricated["Source"]["Directory"] = str(staging)
        fabricated["Source"]["Files"] = {role: (entry if role == "SemanticContract" else probe_case["Source"]["Files"][role])
                                         for role, entry in fabricated["Source"]["Files"].items()}
        probe_manifest["Cases"] = [fabricated if item["Id"] == fabricated["Id"] else item for item in probe_manifest["Cases"]]
    probe_manifest_path = work / "probe-manifest.json"
    probe_manifest_path.write_text(json.dumps(probe_manifest, indent=2) + "\n")
    probe_root = work / "probe"
    summary = probe(probe_manifest_path, case["Id"], probe_root, work / "probe.log")
    if summary["Status"] != "built":
        error = RegistrationError(summary.get("Message") or f"probe build {summary['Status']}")
        error.summary, error.probe_root = summary, probe_root
        raise error
    contract = derive_contract(staging, probe_root / "build-census.json", signature=staged["Signature"],
                               boundary=staged["Boundary"], process_library=process_library, kind=kind)
    return contract, probe_root, summary


def retired_entry(old_case, new_entries, version, reason, evidence):
    """The RetiredFixtures entry of a superseded case: every binding whose digest
    changed or vanished, the derived contract always."""
    old_files = old_case["Source"]["Files"]
    changed = {role: entry for role, entry in old_files.items()
               if role == "SemanticContract" or new_entries.get(role, {}).get("SHA256") != entry["SHA256"]}
    return {"Id": old_case["Id"], "RetiredVersion": version,
            "Retired": f"{time.strftime('%Y-%m-%d')} (register_case.py)",
            "RetiredBinding": changed, "Reason": reason, "Evidence": evidence}


class PreparedRegistration:
    """The outcome of the manifest-free part of a registration (prepare_registration):
    the record so far and, for a case to register, the manifest case entry, its derived
    contract, the probe root / summary and the source paths the commit needs."""

    def __init__(self, record, work, *, case=None, contract=None, probe_root=None, summary=None, directory=None,
                 footprint=None, scope=None, provenance=None, features=None, error=None):
        self.record, self.work = record, work
        self.case, self.contract, self.probe_root, self.summary = case, contract, probe_root, summary
        self.directory, self.footprint, self.scope, self.provenance, self.features = directory, footprint, scope, provenance, features
        self.error = error

    @property
    def done(self):
        """True when the outcome is final (reused / unsupported / failed): nothing to commit."""
        return self.case is None


def finish_record(record, work, status, **fields):
    record.update(Status=status, **fields)
    (work / REGISTRATION_RECORD).write_text(json.dumps(record, indent=2) + "\n")
    return record


def load_production_manifest(manifest_path):
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("Pipeline") != "gmsh-only" or "Calibration" in manifest:
        raise RegistrationError("cases are registered into the Gmsh-only production manifest only")
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    validate_manifest(manifest, manifest_path)
    return manifest_path, manifest, repository


def prepare_registration(case_id, directory, *, footprint, inventory_status, manifest_path, mesh_recipe=None,
                         features=None, provenance=None, work=None, python=sys.executable, julia=None, probe=None,
                         kind="fabricated", fabricated_case=None):
    """Steps 1-2 of a registration (the manifest is read, never written): the source
    digests, the scope classes, the reuse-by-content check against the recorded case,
    then the two-pass contract derivation with the labels-only probe.  Independent per
    case, so device_coupons runs it for several coupons at once; commit_registration
    appends the outcome to the manifest serially.  Returns a PreparedRegistration; a
    final outcome (reused / unsupported / failed) carries `error` (the RegistrationError
    of a stop) instead of a case to commit."""
    manifest_path, manifest, repository = load_production_manifest(manifest_path)
    if inventory_status not in INVENTORY_STATUSES:
        raise RegistrationError(f"--inventory-status must be one of {list(INVENTORY_STATUSES)}")
    if not case_id or any(character in case_id for character in "/ \t\n"):
        raise RegistrationError("the case id must be a nonempty token")
    if kind not in CASE_KINDS or (kind == "thin") != (fabricated_case is not None):
        raise RegistrationError(f"the kind must be one of {list(CASE_KINDS)}; a thin case names its fabricated case")
    if kind == "thin":
        # Decision 66: the thin case pairs with a registered fabricated case of the same
        # source directory under a manifest recording the thin convention (ThinRecipe).
        if validate_thin_recipe(manifest.get("ProductionRecipe") or {}) is None:
            raise RegistrationError("the manifest's ProductionRecipe carries no ThinRecipe: no thin convention to build")
        paired = next((item for item in manifest["Cases"] if item["Id"] == fabricated_case), None)
        if paired is None or case_kind(paired) != "fabricated":
            raise RegistrationError(f"the fabricated case {fabricated_case!r} is not registered")
        if case_directory(paired, repository).resolve() != Path(directory).resolve():
            raise RegistrationError(f"the fabricated case {fabricated_case} was registered from another directory "
                                    f"({paired['Source']['Directory']})")
    directory = Path(directory).resolve()
    if not directory.is_dir():
        raise RegistrationError(f"source directory {directory} does not exist")
    commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=repository, text=True).strip()
    work = Path(work) if work is not None else Path(
        f"/tmp/coupon-register-{case_id}-{commit}-{time.strftime('%Y%m%d-%H%M%S')}")
    work.mkdir(parents=True, exist_ok=True)
    record = {"Case": case_id, "Kind": kind, "FabricatedCase": fabricated_case,
              "Manifest": str(manifest_path), "SourceDirectory": str(directory),
              "Footprint": footprint, "InventoryStatus": inventory_status, "Commit": commit,
              "Work": str(work), "Status": None, "FixtureVersion": None, "SourceSHA256": None,
              "Scope": None, "StoppedBy": None, "ProbeRoot": None, "ContractSHA256": None,
              "RetiredVersion": None, "Message": None}

    def finish(status, **fields):
        return finish_record(record, work, status, **fields)

    try:
        paths = source_files(directory, footprint)
        recipe_entry, recipe_path = mesh_recipe_binding(directory, repository, mesh_recipe)
        shared = shared_case_fields(manifest)
    except RegistrationError as error:
        finish(STATUS_FAILED, Message=str(error))
        raise
    entries = file_entries(paths, recipe_entry)
    record["SourceSHA256"] = source_digests(entries)
    scope = preflight_recipe_scope(manifest, paths, {"Source": {"Files": entries}, CASE_KIND_KEY: kind})
    record["Scope"] = scope
    existing = next((item for item in manifest["Cases"] if item["Id"] == case_id), None)
    if existing is not None and case_kind(existing) != kind:
        finish(STATUS_FAILED, Message=f"{case_id} is registered as a {case_kind(existing)} case: a case never changes kind")
        return PreparedRegistration(record, work, error=RegistrationError(record["Message"]))
    if existing is not None:
        old_entries = existing["Source"]["Files"]
        if source_digests(old_entries) == source_digests(entries):
            # Same content: the recorded case stands (whatever directory it was registered
            # from) provided its frozen contract is intact.
            old_directory = case_directory(existing, repository)
            contract = old_directory / old_entries["SemanticContract"]["Name"]
            if not contract.is_file() or sha256(contract) != old_entries["SemanticContract"]["SHA256"]:
                finish(STATUS_FAILED, Message=f"the frozen contract of {case_id} ({contract}) differs from the "
                                              f"recorded digest while every source is unchanged: restore it")
                return PreparedRegistration(record, work, error=RegistrationError(record["Message"]))
            finish(STATUS_REUSED, FixtureVersion=existing.get("FixtureVersion"),
                   ContractSHA256=old_entries["SemanticContract"]["SHA256"],
                   Message=f"every recorded source digest equals the directory's: case reused "
                           f"(recorded directory {existing['Source']['Directory']})")
            return PreparedRegistration(record, work)
    if scope["UnsupportedClasses"]:
        guard = scope["UnsupportedClasses"][0]
        finish(STATUS_UNSUPPORTED, StoppedBy={"Kind": "ScopeGuard", "Id": guard, "Stage": "inputs"},
               Message=f"unsupported class {guard} (from the inputs; no probe built)")
        return PreparedRegistration(record, work, error=RegistrationError(record["Message"]))
    case = {"Id": case_id, "InventoryStatus": inventory_status, "Frozen": True,
            "Features": list(features) if features else list(scope["ExhibitedClasses"]),
            "Variants": shared["Variants"], "TransformComparison": shared["TransformComparison"],
            "Source": {"Directory": recorded_directory(directory, repository), "SignatureRole": "Signature",
                       "SignatureColumns": shared["SignatureColumns"], "Files": dict(entries)},
            "FixtureVersion": None, "Provenance": None}
    if kind == "thin":
        case[CASE_KIND_KEY] = kind
        case[FABRICATED_CASE_KEY] = fabricated_case
    if footprint != FOOTPRINT_BOUND:
        case["Source"]["EtchFootprint"] = PRODUCER_DEFAULT_ETCH_FOOTPRINT
    probe = probe if probe is not None else (
        lambda probe_manifest, probe_case_id, root, log: run_probe_build(python, julia, probe_manifest,
                                                                          probe_case_id, root, log))
    try:
        contract, probe_root, summary = derive_two_pass(manifest, repository, case, paths, recipe_path, work, probe)
    except RegistrationError as error:
        summary = getattr(error, "summary", None)
        if summary is not None and summary["Status"] == STATUS_UNSUPPORTED:
            finish(STATUS_UNSUPPORTED, ProbeRoot=str(error.probe_root),
                   StoppedBy={"Kind": "ScopeGuard", "Id": summary["ScopeGuard"], "Stage": summary["Stage"]},
                   Message=str(error))
        else:
            stopped = None
            if summary is not None:
                stopped = {"Kind": "HeadroomGate" if summary.get("Stage") == "headroom-gate" else "Stage",
                           "Id": summary.get("Stage"), "ReturnCode": summary.get("ReturnCode")}
            finish(STATUS_FAILED, ProbeRoot=str(getattr(error, "probe_root", "")) or None, StoppedBy=stopped,
                   Message=str(error))
        return PreparedRegistration(record, work, error=error)
    return PreparedRegistration(record, work, case=case, contract=contract, probe_root=probe_root, summary=summary,
                                directory=directory, footprint=footprint, scope=scope, provenance=provenance,
                                features=features)


def commit_registration(prepared, *, manifest_path, refreeze_calibration=None):
    """Step 3-4 of a registration: the contract written to the source directory, the
    case appended to (or superseding its version in) the manifest read afresh, the
    manifest validated and written, the tools refrozen.  Serial: one manifest
    read-modify-write per call."""
    if prepared.done:
        if prepared.error is not None:
            raise prepared.error
        return prepared.record
    record, work, case, directory = prepared.record, prepared.work, prepared.case, prepared.directory
    case_id, commit = record["Case"], record["Commit"]
    manifest_path, manifest, repository = load_production_manifest(manifest_path)
    entries = case["Source"]["Files"]
    existing = next((item for item in manifest["Cases"] if item["Id"] == case_id), None)
    version = 1 if existing is None else int(existing.get("FixtureVersion") or 1) + 1
    case["FixtureVersion"] = version
    contract_file = CONTRACT_FILES[case_kind(case)]
    contract_path = directory / contract_file
    contract_path.write_text(json.dumps(prepared.contract, indent=2) + "\n")
    case["Source"]["Files"]["SemanticContract"] = {"Name": contract_file, "SHA256": sha256(contract_path)}
    case["Provenance"] = (
        f"Registered by register_case.py at {commit} ({time.strftime('%Y-%m-%d')}), fixture version {version}: "
        + (f"THIN case of {case[FABRICATED_CASE_KEY]} (decision 66: the same frozen sources, the thin label "
           f"families, the ThinRecipe options), " if case_kind(case) == "thin" else "")
        + f"source directory {case['Source']['Directory']}, etch footprint {prepared.footprint}, recipe scope classes "
        f"{prepared.scope['ExhibitedClasses']}; the contract is derived by derive_semantic_contract.py from the inputs "
        f"and the label census of a labels-only production-option probe pass (run_gmsh_only_case.py --labels-only, "
        f"decision 62(2); {prepared.probe_root}, commit {prepared.summary.get('Commit')})"
        + (f"; {prepared.provenance}" if prepared.provenance else "")
        + (f"; version {version - 1} retired in RetiredFixtures" if existing is not None else ""))
    if existing is not None:
        old_entries = existing["Source"]["Files"]
        changed = sorted(role for role in set(old_entries) | set(entries)
                         if role != "SemanticContract" and
                         old_entries.get(role, {}).get("SHA256") != entries.get(role, {}).get("SHA256"))
        retired = manifest.setdefault("RetiredFixtures", {"Rule": RETIRED_FIXTURES_RULE, "Entries": []})
        retired["Entries"].append(retired_entry(
            existing, case["Source"]["Files"], version - 1,
            f"source roles changed: {changed}; re-derived and re-registered by register_case.py as "
            f"FixtureVersion {version} at {commit}",
            f"{work / REGISTRATION_RECORD}; probe root {prepared.probe_root}"))
        record["RetiredVersion"] = version - 1
        manifest["Cases"] = [case if item["Id"] == case_id else item for item in manifest["Cases"]]
    else:
        manifest["Cases"].append(case)
    validate_manifest(manifest, manifest_path)
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n")
    # refreeze_manifest_tools freezes the stage tools relative to the manifest's own
    # directory (the tool directory): a manifest elsewhere carries the mirrored digests
    # and is not refrozen here (recorded).
    if manifest_path.parent == HERE:
        calibration = (list(CALIBRATION_MANIFESTS) if refreeze_calibration is None and manifest_path == PRODUCTION_MANIFEST
                       else list(refreeze_calibration or []))
        changes, mirror_stale = refreeze(manifest_path, *calibration, check_only=False)
        refrozen = {"Changes": [list(change) for change in changes], "MirrorUpdated": mirror_stale}
    else:
        refrozen = f"skipped: the manifest is not in {HERE} (stage tools are frozen relative to it)"
    return finish_record(record, work, STATUS_REGISTERED, FixtureVersion=version, ProbeRoot=str(prepared.probe_root),
                         ContractSHA256=case["Source"]["Files"]["SemanticContract"]["SHA256"],
                         Refreeze=refrozen,
                         Message=f"registered {case_id} as FixtureVersion {version}")


def register(case_id, directory, *, footprint, inventory_status, manifest_path, mesh_recipe=None,
             features=None, provenance=None, work=None, python=sys.executable, julia=None,
             probe=None, refreeze_calibration=None, kind="fabricated", fabricated_case=None):
    """Register (or reuse / supersede) the case; returns the registration record
    (prepare_registration then commit_registration).  `kind` thin with
    `fabricated_case` registers the thin counterpart of a registered fabricated case
    (decision 66)."""
    prepared = prepare_registration(case_id, directory, footprint=footprint, inventory_status=inventory_status,
                                    manifest_path=manifest_path, mesh_recipe=mesh_recipe, features=features,
                                    provenance=provenance, work=work, python=python, julia=julia, probe=probe,
                                    kind=kind, fabricated_case=fabricated_case)
    return commit_registration(prepared, manifest_path=manifest_path, refreeze_calibration=refreeze_calibration)


def add_arguments(parser):
    parser.add_argument("case_id")
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("--footprint", required=True, choices=FOOTPRINT_DECLARATIONS,
                        help="the etch footprint declaration (mandatory; fail closed without it)")
    parser.add_argument("--inventory-status", required=True, choices=INVENTORY_STATUSES)
    parser.add_argument("--manifest", type=Path, default=PRODUCTION_MANIFEST)
    parser.add_argument("--mesh-recipe", help="repository path of the frozen mesh recipe when the source "
                                              "directory carries no mesh-recipe.json")
    parser.add_argument("--feature", action="append", default=None,
                        help="case feature (repeatable; default: the exhibited recipe scope classes)")
    parser.add_argument("--provenance", help="text appended to the Provenance statement (e.g. the origin)")
    parser.add_argument("--work", type=Path, help="work directory (default /tmp/coupon-register-<case>-<commit>-<ts>)")
    parser.add_argument("--julia", default=shutil.which("julia"))
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--thin-of", metavar="FABRICATED_CASE_ID",
                        help="register CASE_ID as the thin counterpart of this registered fabricated case of the same "
                             "source directory (decision 66; the manifest must carry ProductionRecipe.ThinRecipe)")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    args = parser.parse_args(argv)
    try:
        record = register(args.case_id, args.source_dir, footprint=args.footprint,
                          inventory_status=args.inventory_status, manifest_path=args.manifest,
                          mesh_recipe=args.mesh_recipe, features=args.feature, provenance=args.provenance,
                          work=args.work, python=args.python, julia=args.julia,
                          kind="thin" if args.thin_of else "fabricated", fabricated_case=args.thin_of)
    except RegistrationError as error:
        print(f"REGISTRATION_FAILED {args.case_id}: {error}", file=sys.stderr)
        return 1
    print(f"{record['Status'].upper()} {args.case_id}: {record['Message']} ({record['Work']}/{REGISTRATION_RECORD})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
