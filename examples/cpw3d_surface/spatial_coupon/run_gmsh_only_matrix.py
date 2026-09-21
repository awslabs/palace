#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Build a set of Gmsh-only manifest cases as a job pool and consolidate one
library-build.json (supervisor decision 48, `coupon-library build` step 2).

Every selected case runs the unchanged per-case driver run_gmsh_only_case.py - its
headroom gate (pre-build element estimate against MaximumElements), the canonical DAG,
the rigid placements, the consolidated audits and verify_canonical_case_entries.py -
under the manifest's stage bounds (MaximumSeconds / MaximumRSSGiB for every stage,
audits and verification included), --jobs cases at a time (default 2).  A case fails
closed on its own (headroom gate, ScopeGuard, a stage stop, verification) and the
others continue; the exit status is nonzero unless every case passed.

library-build.json records per case: the recipe scope classes of its inputs, the
Status (built / unsupported-class / failed) and the exact guard / gate / stage that
stopped it, the calibration label and options when the manifest is a labeled
calibration manifest (Library.Manifest.Kind calibration, never mixed with production
records), CanonicalBuildId, the identity and rotate-z mesh SHA256 and paths, elements
by type, H1 DOFs at --h1-order with the entity counts (vertices, edges, faces, cells by
type: Palace's H1 closed form at any order), the pre-build estimate against the actual count and the
cap (EstimateOverCap), wall seconds and peak GiB of every bounded stage, the
verification verdict, and HeadroomFlags: any measure at or above HEADROOM_FRACTION of
its bound (elements or estimate vs the cap, a stage's seconds or peak RSS vs its
limit).  Library totals: cases attempted / built / passed / unsupported / failed, the
flagged cases, wall clock.  Nothing here is an evidence tool: every number is read from
the per-case root's records.

usage: run_gmsh_only_matrix.py [--manifest PATH] [--case ID ...] [--jobs N] [--root DIR]
       [--output PATH] [--h1-order P] [--julia PATH] [--python PATH]
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from estimate_build_cost import actual_counts, case_paths  # noqa: E402
from general_mesh_manifest import preflight_recipe_scope  # noqa: E402
from run_gmsh_only_case import BUILD_SUMMARY  # noqa: E402

LIBRARY_BUILD_RECORD = "library-build.json"
LIBRARY_BUILD_VERSION = 1
HEADROOM_FRACTION = 0.9
HEADROOM_RULE = (f"a case is flagged when its element count or pre-build estimate reaches "
                 f"{HEADROOM_FRACTION} x MaximumElements, or any bounded stage's wall seconds or peak "
                 f"process-tree RSS reaches {HEADROOM_FRACTION} x its limit (the ten-edge production root "
                 f"sits at 0.87 of the cap and 0.96 of the memory bound: the margin rule the throughput "
                 f"plan lacked)")
STATUS_BUILT, STATUS_UNSUPPORTED, STATUS_FAILED, STATUS_UNBUILT = "built", "unsupported-class", "failed", "registered-unbuilt"
HEADROOM_GATE_STAGE = "headroom-gate"
VERIFICATION_STAGE = "per-entry-verification"
STAGE_REPORT_SUFFIXES = (".audit.log.json", ".log.json")
ELEMENT_TYPES = ("Tetrahedron", "Prism", "Pyramid")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    path = Path(path)
    return json.loads(path.read_text()) if path.is_file() else None


def stage_reports(root):
    """Stage name -> wall seconds, peak GiB, return code, stop reason and limits of every
    bounded stage report (run_bounded_mesher.py `<stage>.log.json`) under the root."""
    stages = {}
    for path in sorted(Path(root).glob("*.log.json")):
        name = path.name
        for suffix in STAGE_REPORT_SUFFIXES:
            if name.endswith(suffix):
                name = name[:-len(suffix)]
                break
        report = json.loads(path.read_text())
        stages[name] = {"Seconds": report["Seconds"],
                        "PeakGiB": report["PeakProcessTreeRSSBytes"] / 2**30,
                        "ReturnCode": report["ReturnCode"], "StopReason": report["StopReason"],
                        "Limits": report["Limits"], "Report": str(path)}
    return stages


def stopped_by(summary, stages, verification):
    """The exact guard / gate / stage / verification that stopped a case, or None."""
    if summary is None:
        return {"Kind": "Driver", "Id": "run_gmsh_only_case.py",
                "Message": f"no {BUILD_SUMMARY} was written (the driver stopped before its first stage)"}
    status, stage = summary["Status"], summary.get("Stage")
    if status == STATUS_BUILT:
        return None
    if status == STATUS_UNSUPPORTED:
        return {"Kind": "ScopeGuard", "Id": summary["ScopeGuard"], "Stage": stage, "Message": summary["Message"]}
    if stage == HEADROOM_GATE_STAGE:
        return {"Kind": "HeadroomGate", "Id": "MaximumElements", "Stage": stage, "Message": summary["Message"]}
    if stage == VERIFICATION_STAGE:
        return {"Kind": "Verification", "Id": stage, "Stage": stage, "Message": summary["Message"],
                "Failures": verification["Failures"] if verification else None}
    report = stages.get(stage) or {}
    return {"Kind": "Stage", "Id": stage, "Stage": stage, "Message": summary["Message"],
            "StopReason": report.get("StopReason"), "ReturnCode": summary.get("ReturnCode")}


def headroom_flags(elements, estimate, stages, gates):
    flags = []
    cap = gates["MaximumElements"]
    if elements is not None and elements >= HEADROOM_FRACTION * cap:
        flags.append(f"elements {elements} >= {HEADROOM_FRACTION} x MaximumElements {cap}")
    if estimate is not None and estimate["EstimateOverCap"] >= HEADROOM_FRACTION:
        flags.append(f"estimate {estimate['EstimatedElements']:.0f} >= {HEADROOM_FRACTION} x MaximumElements {cap}")
    for name, stage in stages.items():
        if stage["Seconds"] >= HEADROOM_FRACTION * stage["Limits"]["Seconds"]:
            flags.append(f"stage {name} {stage['Seconds']:.0f} s >= {HEADROOM_FRACTION} x {stage['Limits']['Seconds']} s")
        if stage["PeakGiB"] >= HEADROOM_FRACTION * stage["Limits"]["MemoryGiB"]:
            flags.append(f"stage {name} {stage['PeakGiB']:.2f} GiB >= {HEADROOM_FRACTION} x "
                         f"{stage['Limits']['MemoryGiB']} GiB")
    return flags


def h1_record_of_mesh(mesh_path, order):
    """H1 DOFs at `order` with the entity counts they are computed from (the qualify
    step estimates the other orders from the same counts without re-reading the mesh)."""
    from mesh_array_io import read_mesh
    from mixed_mesh import h1_dofs_from_counts, h1_entity_counts
    counts = h1_entity_counts(read_mesh(mesh_path))
    return {"Order": order, "DOFs": h1_dofs_from_counts(counts, order), "EntityCounts": counts}


def case_record(manifest, manifest_path, case, root, *, h1_order, driver_return_code, wall_seconds):
    """The library-build record of one case from its root's records alone."""
    root = Path(root)
    case_id = case["Id"]
    paths = case_paths(manifest, manifest_path, case)
    scope = preflight_recipe_scope(manifest, paths, case)
    summary = read_json(root / BUILD_SUMMARY)
    estimate = read_json(root / "build-cost-estimate.json")
    verification = read_json(root / f"{VERIFICATION_STAGE}.json")
    stages = stage_reports(root)
    elements = None
    census = root / "build-census.json"
    if census.is_file():
        counts = actual_counts(census)
        elements = {**counts, "Total": sum(counts.values())}
    estimate_record = None
    if estimate is not None:
        estimate_record = {key: estimate[key] for key in ("EstimatedElements", "EstimateOverCap",
                                                          "MaximumElements", "Passed")}
        estimate_record["ActualElements"] = elements["Total"] if elements else None
        estimate_record["EstimateOverActual"] = (estimate["EstimatedElements"] / elements["Total"]
                                                 if elements else None)
    variants = {}
    for variant in case["Variants"]:
        evidence = read_json(root / "audits" / f"{case_id}--{variant['Id']}.json")
        mesh = root / f"{variant['Id']}.msh"
        if evidence is not None:
            variants[variant["Id"]] = {"Path": evidence["Mesh"]["Path"], "SHA256": evidence["Mesh"]["SHA256"]}
        elif mesh.is_file():
            variants[variant["Id"]] = {"Path": str(mesh), "SHA256": sha256(mesh)}
        else:
            variants[variant["Id"]] = None
    built = summary is not None and summary["Status"] == STATUS_BUILT
    passed = bool(built and verification is not None and verification["Passed"])
    h1 = None
    identity = root / "identity.msh"
    if built and identity.is_file():
        h1 = h1_record_of_mesh(identity, h1_order)
    status = summary["Status"] if summary is not None else STATUS_FAILED
    canonical = None
    if verification is not None:
        canonical = verification["SharedCanonicalBuildId"]
    elif (root / "canonical-build.json").is_file():
        canonical = read_json(root / "canonical-build.json").get("CanonicalBuildId")
    calibration = case.get("Calibration") if "Calibration" in manifest else None
    return {"Case": case_id, "InventoryStatus": case["InventoryStatus"], "FixtureVersion": case.get("FixtureVersion"),
            "Calibration": ({"Label": calibration["Label"], "BaseCase": calibration.get("BaseCase"),
                             "BuildCommandOptions": calibration.get("BuildCommandOptions"),
                             "ProductionValues": calibration["ProductionValues"]} if calibration else None),
            "Scope": scope, "Status": status, "Passed": passed,
            "StoppedBy": stopped_by(summary, stages, verification) if not passed else None,
            "CanonicalBuildId": canonical, "Variants": variants,
            "Elements": elements, "H1": h1, "Estimate": estimate_record, "Stages": stages,
            "Verification": (None if verification is None else
                             {"Passed": verification["Passed"], "Failures": verification["Failures"],
                              "Path": str(root / f"{VERIFICATION_STAGE}.json")}),
            "HeadroomFlags": headroom_flags(elements["Total"] if elements else None, estimate_record, stages,
                                            manifest["Gates"]),
            "Root": str(root), "DriverReturnCode": driver_return_code, "WallSeconds": wall_seconds}


def run_case(manifest_path, manifest, case, root, *, python, julia, h1_order):
    """One job of the pool: the per-case driver under the manifest bounds, then the record."""
    root = Path(root)
    command = [python, str(HERE / "run_gmsh_only_case.py"), case["Id"], "--manifest", str(manifest_path),
               "--root", str(root), "--audit-memory-gib", repr(float(manifest["Gates"]["MaximumRSSGiB"]))]
    if julia is not None:
        command += ["--julia", str(julia)]
    start = time.time()
    root.parent.mkdir(parents=True, exist_ok=True)
    with open(root.parent / f"{case['Id']}.driver.log", "w") as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
    return case_record(manifest, manifest_path, case, root, h1_order=h1_order,
                       driver_return_code=result.returncode, wall_seconds=time.time() - start)


def unbuilt_record(manifest, manifest_path, case, estimate, rank, limit):
    """The record of a selected case the --build-limit left unbuilt: its scope, its
    pre-build estimate and the rule that ranked it (nothing is built or judged)."""
    paths = case_paths(manifest, manifest_path, case)
    return {"Case": case["Id"], "InventoryStatus": case["InventoryStatus"], "FixtureVersion": case.get("FixtureVersion"),
            "Calibration": None, "Scope": preflight_recipe_scope(manifest, paths, case), "Status": STATUS_UNBUILT,
            "Passed": False,
            "StoppedBy": {"Kind": "BuildLimit", "Id": "--build-limit", "Stage": None,
                          "Message": f"ranked {rank} by the pre-build element estimate; --build-limit {limit} builds the "
                                     f"{limit} smallest selected cases"},
            "CanonicalBuildId": None, "Variants": {variant["Id"]: None for variant in case["Variants"]},
            "Elements": None, "H1": None,
            "Estimate": {key: estimate[key] for key in ("EstimatedElements", "EstimateOverCap", "MaximumElements", "Passed")}
            | {"ActualElements": None, "EstimateOverActual": None},
            "Stages": {}, "Verification": None, "HeadroomFlags": [], "Root": None, "DriverReturnCode": None, "WallSeconds": 0.0}


def library_totals(records, *, jobs, wall_seconds, manifest, manifest_path, commit):
    return {"CasesAttempted": len(records),
            "CasesBuilt": sum(record["Status"] == STATUS_BUILT for record in records),
            "CasesPassed": sum(record["Passed"] for record in records),
            "CasesUnsupported": sum(record["Status"] == STATUS_UNSUPPORTED for record in records),
            "CasesFailed": sum(record["Status"] == STATUS_FAILED for record in records),
            "CasesUnbuilt": sum(record["Status"] == STATUS_UNBUILT for record in records),
            "FlaggedCases": [record["Case"] for record in records if record["HeadroomFlags"]],
            "WallClockSeconds": wall_seconds, "Jobs": jobs, "Commit": commit,
            "Bounds": {key: manifest["Gates"][key] for key in ("MaximumElements", "MaximumSeconds", "MaximumRSSGiB")},
            "HeadroomFraction": HEADROOM_FRACTION, "HeadroomRule": HEADROOM_RULE,
            "Manifest": {"Path": str(manifest_path), "SHA256": sha256(manifest_path),
                         "Kind": "calibration" if "Calibration" in manifest else "production",
                         "ProductionManifest": (str((manifest_path.parent / manifest["Calibration"]["ProductionManifest"]).resolve())
                                                if "Calibration" in manifest else None)}}


def run_matrix(manifest_path, case_ids=None, *, root, jobs=2, output=None, python=sys.executable,
               julia=None, h1_order=4, build_limit=None, extra=None):
    """Build the selected cases (default: every case) as a pool of `jobs` drivers and
    write the library-build record; returns it.  `build_limit` builds only the N
    selected cases with the smallest pre-build element estimate and records the others
    as registered-unbuilt with their estimates; `extra` (a dict) is merged into the
    record (the device adapter's record under Device)."""
    manifest_path = Path(manifest_path).resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("Pipeline") != "gmsh-only":
        raise ValueError("the matrix runs a Gmsh-only manifest (the production manifest or a labeled calibration manifest)")
    if jobs < 1:
        raise ValueError("--jobs must be positive")
    by_id = {case["Id"]: case for case in manifest["Cases"]}
    selected = list(by_id) if not case_ids else list(case_ids)
    unknown = [case_id for case_id in selected if case_id not in by_id]
    if unknown or len(set(selected)) != len(selected):
        raise ValueError(f"unknown or repeated cases {unknown or selected}")
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=repository, text=True).strip()
    root = Path(root)
    root.mkdir(parents=True, exist_ok=True)
    output = Path(output) if output is not None else root / LIBRARY_BUILD_RECORD
    start = time.time()
    unbuilt = []
    if build_limit is not None:
        if build_limit < 0:
            raise ValueError("--build-limit must be nonnegative")
        from estimate_build_cost import gate as estimate_gate
        estimates = {case_id: estimate_gate(manifest, manifest_path, by_id[case_id]) for case_id in selected}
        ranked = sorted(selected, key=lambda case_id: (estimates[case_id]["EstimatedElements"], case_id))
        unbuilt = [(case_id, rank) for rank, case_id in enumerate(ranked, start=1) if rank > build_limit]
        selected = ranked[:build_limit]
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(run_case, manifest_path, manifest, by_id[case_id], root / case_id,
                               python=python, julia=julia, h1_order=h1_order) for case_id in selected]
        records = [future.result() for future in futures]
    records += [unbuilt_record(manifest, manifest_path, by_id[case_id], estimates[case_id], rank, build_limit)
                for case_id, rank in unbuilt]
    record = {"Version": LIBRARY_BUILD_VERSION, "Command": "coupon-library build",
              "Root": str(root), "Cases": records,
              "Library": library_totals(records, jobs=jobs, wall_seconds=time.time() - start,
                                        manifest=manifest, manifest_path=manifest_path, commit=commit)}
    if build_limit is not None:
        record["Library"]["BuildLimit"] = {"Limit": build_limit, "Rule": "the selected cases ranked by the pre-build "
                                           "element estimate (estimate_build_cost.gate); the smallest Limit built, the "
                                           "others recorded registered-unbuilt with their estimates",
                                           "Ranked": ranked}
    if extra:
        record.update(extra)
    output.write_text(json.dumps(record, indent=2) + "\n")
    return record


def add_arguments(parser):
    parser.add_argument("--manifest", type=Path, default=HERE / "geometry-independence-suite.json")
    parser.add_argument("--case", action="append", default=None, help="case id (repeatable; default: every case)")
    parser.add_argument("--jobs", type=int, default=2, help="concurrent case drivers (default 2)")
    parser.add_argument("--root", type=Path, help="library root (default /tmp/coupon-library-build-<commit>-<ts>)")
    parser.add_argument("--output", type=Path, help=f"record path (default ROOT/{LIBRARY_BUILD_RECORD})")
    parser.add_argument("--h1-order", type=int, default=4, help="H1 DOF order recorded per built case")
    parser.add_argument("--build-limit", type=int, help="build only this many of the selected cases (the smallest by the "
                                                        "pre-build element estimate); the others are recorded "
                                                        "registered-unbuilt with their estimates")
    parser.add_argument("--julia", default=shutil.which("julia"))
    parser.add_argument("--python", default=sys.executable)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    add_arguments(parser)
    args = parser.parse_args(argv)
    return run_build(args)


def run_build(args, extra=None):
    """Run the matrix from parsed arguments; prints one line per case and the totals."""
    manifest_path = args.manifest.resolve()
    root = args.root
    if root is None:
        repository = (manifest_path.parent / json.loads(manifest_path.read_text())["RepositoryRoot"]).resolve()
        commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=repository, text=True).strip()
        root = Path(f"/tmp/coupon-library-build-{commit}-{time.strftime('%Y%m%d-%H%M%S')}")
    record = run_matrix(manifest_path, args.case, root=root, jobs=args.jobs, output=args.output,
                        python=args.python, julia=args.julia, h1_order=args.h1_order,
                        build_limit=getattr(args, "build_limit", None), extra=extra)
    for case in record["Cases"]:
        stopped = case["StoppedBy"]
        print(f"{case['Case']}: {case['Status']} passed={case['Passed']}"
              + (f" elements={case['Elements']['Total']}" if case["Elements"] else "")
              + (f" stopped-by {stopped['Kind']}[{stopped['Id']}]" if stopped else "")
              + (f" FLAGS {case['HeadroomFlags']}" if case["HeadroomFlags"] else ""), flush=True)
    totals = record["Library"]
    print(f"LIBRARY attempted {totals['CasesAttempted']} built {totals['CasesBuilt']} passed {totals['CasesPassed']} "
          f"unsupported {totals['CasesUnsupported']} failed {totals['CasesFailed']} unbuilt {totals['CasesUnbuilt']} "
          f"wall {totals['WallClockSeconds']:.0f} s "
          f"flagged {totals['FlaggedCases']}; record {Path(args.output) if args.output else Path(root) / LIBRARY_BUILD_RECORD}")
    return 0 if totals["CasesPassed"] + totals["CasesUnbuilt"] == totals["CasesAttempted"] else 1


if __name__ == "__main__":
    sys.exit(main())
