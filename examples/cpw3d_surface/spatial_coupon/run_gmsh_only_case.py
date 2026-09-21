#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Build one Gmsh-only manifest case - a production case or a labeled calibration
case (supervisor decision 41) - through the Gmsh-only canonical DAG (supervisor
decision 38) and its evidence chain:

  headroom gate (estimate_build_cost.py: the pre-build element estimate of a production
  case must not exceed MaximumElements; recorded as build-cost-estimate.json)
  -> canonical-source-validation -> gmsh-build -> canonical-gmsh-publication
  -> canonical build record -> proper-rigid-publication per variant
  -> consolidated audits + normalization per variant -> per-entry verification.

Every stage command is generated from the manifest case (immutable inputs, process,
mesh recipe, the build options: ProductionRecipe.BuildCommandOptions of a production
manifest, or a calibration case's Calibration.ProductionValues overridden by its
Calibration.BuildCommandOptions; --trace-basis-size-ratio is one of the options and
has no default here - it is passed, at the manifest's value, exactly when the case
binds a trace basis) and run under run_bounded_mesher.py with the
manifest's stage bounds (audits: --audit-memory-gib).
Nothing here is an evidence tool: the evidence is the bounded stage reports, the
audit records and the verification report written under --root.  The run's outcome is
recorded in --root/build-summary.json (Status built / unsupported-class / failed): a
mesher stop at a recipe scope guard ("ScopeGuard[<id>]" in the gmsh-build log,
supervisor decision 48) is recorded as an unsupported class with the guard id,
distinctly from any other failure.

usage: run_gmsh_only_case.py CASE_ID [--manifest PATH] [--root DIR] [--julia PATH]
       [--python PATH] [--stages-only] [--audits-only]
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import tomllib

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from estimate_build_cost import gate as estimate_gate  # noqa: E402
from general_mesh_manifest import case_gates  # noqa: E402
from mesh_stage_contract import scope_guard_in_text  # noqa: E402
BUILD_SUMMARY = "build-summary.json"
TRACE_BASIS = {"BasisContract": ("source-basis-contract", "--trace-basis-contract"),
               "TraceVertices": ("source-trace-vertices", "--trace-vertices"),
               "TraceTriangles": ("source-trace-triangles", "--trace-triangles"),
               "ProcessLibrary": ("source-process-library", "--process-library")}
CANONICAL_STAGES = ("canonical-source-validation", "gmsh-build", "canonical-gmsh-publication")
TRACE_BASIS_RATIO_OPTION = "--trace-basis-size-ratio"
STAGE_STEMS = {"canonical-source-validation": "canonical-source", "gmsh-build": "gmsh-build",
               "canonical-gmsh-publication": "canonical-publish"}
AUDIT_KINDS = ("bounded-run", "mesh-topology-quality", "mesh-complexity", "mesh-invariants",
               "variant-transform")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value):
    return repr(float(value)) if isinstance(value, float) else str(value)


def case_build_options(manifest, case):
    """The gmsh-build recipe options of a case: a production case executes the
    manifest's ProductionRecipe.BuildCommandOptions; a case of a labeled calibration
    manifest executes its Calibration.ProductionValues overridden by its
    Calibration.BuildCommandOptions (decision 41).  Returns (options without the trace
    basis ratio, trace basis ratio or None, label); the ratio is taken from the
    options alone (no default; production 0.5 since decision 42) and is passed with
    the bound trace basis only.  A case binding a trace basis under a manifest that
    declares no ratio fails closed."""
    calibration = case.get("Calibration") if "Calibration" in manifest else None
    if calibration is not None:
        options = dict(calibration["ProductionValues"], **calibration["BuildCommandOptions"])
        label = (f"CALIBRATION build {case['Id']} ({calibration['Label']}): options "
                 f"{calibration['BuildCommandOptions']} against production {calibration['ProductionValues']}")
    else:
        options = dict(manifest["ProductionRecipe"]["BuildCommandOptions"])
        label = f"Gmsh-only production build {case['Id']}"
    ratio = options.pop(TRACE_BASIS_RATIO_OPTION, None)
    if "BasisContract" in case["Source"]["Files"]:
        if ratio is None:
            raise ValueError(f"{case['Id']} binds a trace basis but the manifest declares no "
                             f"{TRACE_BASIS_RATIO_OPTION}")
        ratio = float(ratio)
    else:
        ratio = None
    return options, ratio, label


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case_id")
    parser.add_argument("--manifest", type=Path, default=HERE / "geometry-independence-suite.json")
    parser.add_argument("--root", type=Path)
    parser.add_argument("--julia", default=shutil.which("julia"))
    parser.add_argument("--python", default=sys.executable)
    parser.add_argument("--julia-project", default="test/examples",
                        help="Julia project (relative to the repository root)")
    parser.add_argument("--audit-memory-gib", default="16")
    parser.add_argument("--stages-only", action="store_true")
    parser.add_argument("--audits-only", action="store_true")
    args = parser.parse_args()
    manifest_path = args.manifest.resolve()
    manifest = json.loads(manifest_path.read_text())
    if manifest.get("Pipeline") != "gmsh-only":
        parser.error("the manifest must be a Gmsh-only (production or labeled calibration) manifest")
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    os.chdir(repository)
    case = next(item for item in manifest["Cases"] if item["Id"] == args.case_id)
    recipe, trace_basis_ratio, label = case_build_options(manifest, case)
    # The manifest Gates are the canonical cache key (gates.json); the mesher's quality
    # bounds are the gates that judge THIS case - a labeled calibration-only per-case
    # deviation (element cap, Jacobian condition bound) applies to that case alone.
    gates = manifest["Gates"]
    judged = case_gates(manifest, case)
    commit = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], text=True).strip()
    root = args.root or Path(f"/tmp/coupon-gmsh-only-{args.case_id}-{commit}-{time.strftime('%Y%m%d-%H%M%S')}")
    root.mkdir(parents=True, exist_ok=args.audits_only)
    source = case["Source"]
    directory = Path(source["Directory"])
    directory = directory if directory.is_absolute() else repository / directory
    paths = {}
    for role, entry in source["Files"].items():
        path = repository / entry["RepositoryPath"] if entry.get("RepositoryPath") else directory / entry["Name"]
        if sha256(path) != entry["SHA256"]:
            raise SystemExit(f"immutable {role} hash mismatch: {path}")
        paths[role] = path
    process = tomllib.loads(paths["Process"].read_text())
    mesh_recipe = json.loads(paths["MeshRecipe"].read_text())
    normal = mesh_recipe["NormalSizeOverThickness"] * process["MetalThickness"]
    tangent = normal * mesh_recipe["TangentialSizeOverNormalSize"]
    far = mesh_recipe["FarSizeOverRadius"] * process["Radius"]
    python, julia = args.python, args.julia
    if julia is None:
        parser.error("--julia is required")
    tools = {name: HERE / name for name in (
        "transform_coupon_source_contract.py", "mesh_spatial_coupon.jl",
        "relabel_frozen_interface_mesh.jl", "publish_rigid_coupon_mesh.py",
        "audit_rigid_coupon_ownership.jl", "general_mesh_audit_producer.py",
        "normalize_general_mesh_evidence.py", "verify_canonical_case_entries.py",
        "run_bounded_mesher.py", "canonical_mesh_build.py")}
    seconds, memory = str(int(gates["MaximumSeconds"])), number(gates["MaximumRSSGiB"])
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", JULIA_NUM_THREADS="1")
    S = {role: str(path) for role, path in paths.items()}
    variants = [variant["Id"] for variant in case["Variants"]]

    def write_summary(status, *, stage=None, return_code=0, scope_guard=None, message=None):
        """The machine-readable outcome of this run (decision 48): Status "built",
        "unsupported-class" (the mesher stopped at a recipe scope guard; ScopeGuard is its
        id) or "failed" (any other stage failure)."""
        (root / BUILD_SUMMARY).write_text(json.dumps({
            "Case": args.case_id, "Commit": commit, "Root": str(root), "Status": status,
            "Stage": stage, "ReturnCode": return_code, "ScopeGuard": scope_guard,
            "Message": message}, indent=2) + "\n")

    def launch(name, *command, memory_gib=memory, check=True):
        with open(root / f"{name}.launch.stdout", "w") as out, open(root / f"{name}.launch.stderr", "w") as err:
            result = subprocess.run([python, str(tools["run_bounded_mesher.py"]), "--seconds", seconds,
                                     "--memory-gib", memory_gib, *command], env=env, stdout=out, stderr=err)
        if check and result.returncode != 0:
            log = root / f"{name}.log"
            guard = scope_guard_in_text(log.read_text(errors="replace")) if log.is_file() else None
            if guard is not None:
                write_summary("unsupported-class", stage=name, return_code=result.returncode,
                              scope_guard=guard, message=f"unsupported class {guard}")
                raise SystemExit(f"UNSUPPORTED_CLASS {guard}: stage {name} stopped at the recipe scope "
                                 f"guard ScopeGuard[{guard}]; see {log} and {root}/{BUILD_SUMMARY}")
            write_summary("failed", stage=name, return_code=result.returncode,
                          message=f"stage {name} failed rc={result.returncode}")
            raise SystemExit(f"stage {name} failed rc={result.returncode}; see {root}/{name}.launch.stderr")
        return result.returncode

    if not args.audits_only:
        # Headroom gate (fail closed before any build): the pre-build element estimate of the
        # case - a calibration case with its own labeled options and the production model -
        # against the unchanged MaximumElements; recorded in the root.
        cost = estimate_gate(manifest, manifest_path, case)
        (root / "build-cost-estimate.json").write_text(json.dumps(cost, indent=2) + "\n")
        print(f"ESTIMATE {args.case_id}: {cost['EstimatedElements']:.0f} elements "
              f"({cost['EstimateOverCap']:.3f} of the cap {cost['MaximumElements']})", flush=True)
        if not cost["Passed"]:
            write_summary("failed", stage="headroom-gate", return_code=1,
                          message=f"pre-build element estimate {cost['EstimatedElements']:.0f} exceeds the cap")
            raise SystemExit(f"pre-build element estimate {cost['EstimatedElements']:.0f} exceeds the cap "
                             f"{cost['MaximumElements']}: not building (see {root}/build-cost-estimate.json)")
        (root / ("CALIBRATION.txt" if "Calibration" in manifest else "PRODUCTION.txt")).write_text(
            f"{label} at {commit} (decision 38): build options {recipe}, {TRACE_BASIS_RATIO_OPTION} "
            f"{trace_basis_ratio if trace_basis_ratio is not None else 'not passed (no trace basis)'}; "
            f"NormalSize {normal} CornerIsotropyRadius {tangent} FarSize {far}; "
            f"process {process}; manifest {manifest_path}\n")
        (root / "canonical-transform.json").write_text("[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]\n")
        (root / "input-hashes.json").write_text(json.dumps({k: v["SHA256"] for k, v in source["Files"].items()}, indent=2))
        (root / "gates.json").write_text(json.dumps(gates, indent=2))
        if judged != gates:
            (root / "case-gates.json").write_text(json.dumps(
                {"Rule": "the gates that judge this case (general_mesh_manifest.case_gates: the manifest Gates with the "
                         "case's labeled calibration-only deviations); gates.json stays the canonical cache key",
                 "Gates": judged, "Deviations": {key: value for key, value in judged.items() if gates.get(key) != value}},
                indent=2) + "\n")
        (root / "canonical-tool-hashes.json").write_text(json.dumps(
            {f"{stage}/{role}": digest for stage in CANONICAL_STAGES
             for role, digest in manifest["StageToolSHA256"][stage].items()}, indent=2))
        for variant in case["Variants"]:
            (root / f"{variant['Id']}-transform.json").write_text(json.dumps(variant["Transform"]) + "\n")
        launch("canonical-source", "--log", f"{root}/canonical-source.log", "--stage", "canonical-source-validation",
               "--input", f"source-semantic-contract={S['SemanticContract']}", "--input", f"source-signature={S['Signature']}",
               "--input", f"source-boundary={S['Boundary']}", "--input", f"source-mask={S['Mask']}",
               "--input", f"canonical-transform={root}/canonical-transform.json",
               "--artifact", f"canonical-semantic-contract={root}/canonical-semantic.json",
               "--artifact", f"canonical-supports={root}/canonical-supports.json",
               "--tool", f"runtime={python}", "--tool", f"source-validator={tools['transform_coupon_source_contract.py']}", "--",
               python, str(tools["transform_coupon_source_contract.py"]), str(directory), f"{root}/canonical-transform.json",
               f"{root}/canonical-semantic.json", f"{root}/canonical-supports.json", "--semantic-input", S["SemanticContract"],
               "--signature", S["Signature"], "--boundary", S["Boundary"], "--mask", S["Mask"])
        basis_inputs, basis_options = [], []
        if "BasisContract" in paths:
            for role, (name, option) in TRACE_BASIS.items():
                basis_inputs += ["--input", f"{name}={S[role]}"]; basis_options += [option, S[role]]
            basis_options += [TRACE_BASIS_RATIO_OPTION, number(trace_basis_ratio)]
        etch_inputs, etch_options = [], []
        if "RetainedEtch" in paths:
            etch_inputs = ["--input", f"source-retained-etch={S['RetainedEtch']}"]
            etch_options = ["--etch-boundary", S["RetainedEtch"]]
        recipe_options = [token for option, value in recipe.items() for token in (option, number(value))]
        with open(paths["Signature"], newline="") as stream:
            slots = {int(row["Slot"]) for row in csv.DictReader(stream)}
        # Multi-slot coupons need the mesher's surface-partition postprocessor.
        ownership_report = (["--interface-ownership-report", f"{root}/build-interface-ownership.csv"]
                            if len(slots) > 1 else [])
        launch("gmsh-build", "--log", f"{root}/gmsh-build.log", "--stage", "gmsh-build",
               "--input", f"source-signature={S['Signature']}", "--input", f"source-boundary={S['Boundary']}",
               "--input", f"source-mask={S['Mask']}", "--input", f"canonical-semantic-contract={root}/canonical-semantic.json",
               *etch_inputs, *basis_inputs,
               "--artifact", f"gmsh-mesh={root}/gmsh-build.msh", "--artifact", f"build-census={root}/build-census.json",
               "--tool", f"runtime={julia}", "--tool", f"mesher={tools['mesh_spatial_coupon.jl']}", "--",
               julia, "--startup-file=no", f"--project={args.julia_project}", str(tools["mesh_spatial_coupon.jl"]),
               S["Signature"], "fabricated", f"{root}/gmsh-build.msh", "--mask", S["Mask"], "--boundary", S["Boundary"],
               "--radius", number(process["Radius"]), "--metal-thickness", number(process["MetalThickness"]),
               "--overetch", number(process["Overetch"]), "--sidewall-angle", number(process["SidewallAngle"]),
               "--top-radius", number(process["TopRounding"]), "--bottom-radius", number(process["TrenchRounding"]),
               "--lc-fine", number(normal), "--lc-far", number(far), "--mesh-order", str(mesh_recipe["GeometryOrder"]),
               "--max-nodes", str(int(judged["MaximumElements"])), "--max-elements", str(int(judged["MaximumElements"])),
               "--semantic-contract", f"{root}/canonical-semantic.json", "--corner-isotropy-radius", number(tangent),
               "--corner-census", f"{root}/build-census.json", *ownership_report, *etch_options, *basis_options,
               "--prism-tubes", "true", *recipe_options,
               "--maximum-corner-aspect", number(judged["MaximumCornerAspect"]),
               "--minimum-scaled-jacobian", number(judged["MinimumScaledJacobian"]),
               "--maximum-jacobian-condition", number(judged["MaximumJacobianCondition"]),
               "--maximum-quality-displacement-over-normal", "0.75")
        launch("canonical-publish", "--log", f"{root}/canonical-publish.log", "--stage", "canonical-gmsh-publication",
               "--input", f"gmsh-mesh={root}/gmsh-build.msh", "--input", f"source-process={S['Process']}",
               "--input", f"source-signature={S['Signature']}", "--input", f"source-boundary={S['Boundary']}",
               "--artifact", f"canonical-candidate-mesh={root}/canonical.msh",
               "--artifact", f"canonical-ownership-partition={root}/canonical.msh.interface-partition.csv",
               "--artifact", f"canonical-ownership-quadrature-partition={root}/canonical.msh.interface-partition.csv.quadrature.csv",
               "--tool", f"runtime={julia}", "--tool", f"publisher={tools['relabel_frozen_interface_mesh.jl']}", "--",
               julia, "--startup-file=no", f"--project={args.julia_project}", str(tools["relabel_frozen_interface_mesh.jl"]),
               str(directory), "fabricated", f"{root}/gmsh-build.msh", f"{root}/canonical.msh", "--process", S["Process"],
               "--signature", S["Signature"], "--boundary", S["Boundary"])
        subprocess.run([python, str(tools["canonical_mesh_build.py"]), f"{root}/input-hashes.json", f"{root}/gates.json",
                        f"{root}/canonical-tool-hashes.json", f"{root}/canonical-build.json",
                        *[token for stage in CANONICAL_STAGES
                          for token in ("--stage-report", f"{stage}={root}/{STAGE_STEMS[stage]}.log.json")]],
                       cwd=HERE, check=True, env=env)
        for variant in variants:
            transform = f"{root}/{variant}-transform.json"; output = f"{root}/{variant}.msh"
            ownership = f"{output}.interface-partition.csv"
            launch(f"{variant}-publication", "--log", f"{root}/{variant}-publication.log", "--stage", "proper-rigid-publication",
                   "--input", f"canonical-candidate-mesh={root}/canonical.msh", "--input", f"canonical-build-record={root}/canonical-build.json",
                   "--input", f"placement-transform={transform}", "--input", f"source-semantic-contract={S['SemanticContract']}",
                   "--input", f"source-signature={S['Signature']}", "--input", f"source-boundary={S['Boundary']}",
                   "--input", f"source-mask={S['Mask']}", "--input", f"source-process={S['Process']}",
                   "--artifact", f"candidate-mesh={output}", "--artifact", f"transformed-semantic-contract={root}/{variant}-semantic.json",
                   "--artifact", f"transformed-supports={root}/{variant}-supports.json", "--artifact", f"ownership-partition={ownership}",
                   "--artifact", f"ownership-quadrature-partition={ownership}.quadrature.csv",
                   "--artifact", f"transform-receipt={root}/{variant}-receipt.json",
                   "--tool", f"runtime={python}", "--tool", f"rigid-publisher={tools['publish_rigid_coupon_mesh.py']}",
                   "--tool", f"ownership-runtime={julia}", "--tool", f"ownership-auditor={tools['audit_rigid_coupon_ownership.jl']}", "--",
                   python, str(tools["publish_rigid_coupon_mesh.py"]), f"{root}/canonical.msh", transform, output,
                   f"{root}/{variant}-receipt.json", "--semantic-input", S["SemanticContract"], "--signature", S["Signature"],
                   "--boundary", S["Boundary"], "--mask", S["Mask"], "--process", S["Process"],
                   "--canonical-build-record", f"{root}/canonical-build.json", "--transformed-semantic", f"{root}/{variant}-semantic.json",
                   "--transformed-supports", f"{root}/{variant}-supports.json", "--ownership", ownership,
                   "--ownership-quadrature", f"{ownership}.quadrature.csv", "--ownership-runtime", julia,
                   "--ownership-auditor", str(tools["audit_rigid_coupon_ownership.jl"]), "--kind", "fabricated")
        (root / "stages.done").touch()
        print(f"STAGES_DONE {root}", flush=True)
    if args.stages_only:
        write_summary("built", stage="stages-only")
        print(root); return
    (root / "audits").mkdir(exist_ok=True)
    stage_reports = [token for stage in CANONICAL_STAGES
                     for token in ("--stage-report", f"{stage}={root}/{STAGE_STEMS[stage]}.log.json")]
    for variant in variants:
        start = time.time()
        code = launch(f"{variant}-variant-audits", "--log", f"{root}/{variant}-variant-audits.audit.log",
                      *[token for kind in AUDIT_KINDS for token in ("--artifact", f"{root}/{variant}-{kind}.json")], "--",
                      python, str(tools["general_mesh_audit_producer.py"]), "variant-audits", args.case_id, variant,
                      f"{root}/{variant}.msh", f"{root}/input-hashes.json", f"{root}/{variant}-transform.json", f"{root}/{variant}-",
                      "--contract", S["SemanticContract"], "--recipe", S["MeshRecipe"], "--process", S["Process"],
                      "--signature", S["Signature"], "--identity-mesh", f"{root}/canonical.msh",
                      "--identity-seed-mesh", f"{root}/gmsh-build.msh", *stage_reports,
                      "--stage-report", f"proper-rigid-publication={root}/{variant}-publication.log.json",
                      memory_gib=args.audit_memory_gib, check=False)
        print(f"{variant} variant-audits rc={code} wall {time.time() - start:.0f} s", flush=True)
        if code != 0:
            continue
        with open(root / f"{variant}-normalize.stdout", "w") as out, open(root / f"{variant}-normalize.stderr", "w") as err:
            subprocess.run([python, str(tools["normalize_general_mesh_evidence.py"]), str(manifest_path), args.case_id, variant,
                            f"{root}/{variant}.msh", f"{root}/audits/{args.case_id}--{variant}.json",
                            *[token for kind in AUDIT_KINDS for token in ("--audit-record", f"{kind}={root}/{variant}-{kind}.json")]],
                           env=env, stdout=out, stderr=err)
    (root / "audits.done").touch()
    code = launch("per-entry-verification", "--log", f"{root}/per-entry-verification.audit.log",
                  "--artifact", f"{root}/per-entry-verification.json", "--",
                  python, str(tools["verify_canonical_case_entries.py"]), str(manifest_path), f"{root}/audits",
                  f"{root}/per-entry-verification.json", args.case_id, memory_gib=args.audit_memory_gib, check=False)
    print(f"verification rc={code}", flush=True)
    write_summary("built" if code == 0 else "failed", stage="per-entry-verification", return_code=code,
                  message=None if code == 0 else "per-entry verification failed")
    print(f"DONE {root}", flush=True)


if __name__ == "__main__":
    main()
