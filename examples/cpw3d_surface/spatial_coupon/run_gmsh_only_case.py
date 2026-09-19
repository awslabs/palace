#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Build one production manifest case through the Gmsh-only canonical DAG
(supervisor decision 38) and its evidence chain:

  canonical-source-validation -> gmsh-build -> canonical-gmsh-publication
  -> canonical build record -> proper-rigid-publication per variant
  -> consolidated audits + normalization per variant -> per-entry verification.

Every stage command is generated from the manifest case (immutable inputs, process,
mesh recipe, ProductionRecipe.BuildCommandOptions) and run under
run_bounded_mesher.py with the manifest's stage bounds (audits: --audit-memory-gib).
Nothing here is an evidence tool: the evidence is the bounded stage reports, the
audit records and the verification report written under --root.

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
TRACE_BASIS = {"BasisContract": ("source-basis-contract", "--trace-basis-contract"),
               "TraceVertices": ("source-trace-vertices", "--trace-vertices"),
               "TraceTriangles": ("source-trace-triangles", "--trace-triangles"),
               "ProcessLibrary": ("source-process-library", "--process-library")}
CANONICAL_STAGES = ("canonical-source-validation", "gmsh-build", "canonical-gmsh-publication")
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
        parser.error("the manifest must be the Gmsh-only production manifest")
    repository = (manifest_path.parent / manifest["RepositoryRoot"]).resolve()
    os.chdir(repository)
    case = next(item for item in manifest["Cases"] if item["Id"] == args.case_id)
    recipe = manifest["ProductionRecipe"]["BuildCommandOptions"]
    gates = manifest["Gates"]
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

    def launch(name, *command, memory_gib=memory, check=True):
        with open(root / f"{name}.launch.stdout", "w") as out, open(root / f"{name}.launch.stderr", "w") as err:
            result = subprocess.run([python, str(tools["run_bounded_mesher.py"]), "--seconds", seconds,
                                     "--memory-gib", memory_gib, *command], env=env, stdout=out, stderr=err)
        if check and result.returncode != 0:
            raise SystemExit(f"stage {name} failed rc={result.returncode}; see {root}/{name}.launch.stderr")
        return result.returncode

    if not args.audits_only:
        (root / "PRODUCTION.txt").write_text(
            f"Gmsh-only production build {args.case_id} at {commit} (decision 38): build options {recipe}; "
            f"NormalSize {normal} TangentialSize {tangent} FarSize {far}; process {process}; manifest {manifest_path}\n")
        (root / "canonical-transform.json").write_text("[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1]\n")
        (root / "input-hashes.json").write_text(json.dumps({k: v["SHA256"] for k, v in source["Files"].items()}, indent=2))
        (root / "gates.json").write_text(json.dumps(gates, indent=2))
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
            basis_options += ["--trace-basis-size-ratio", "1.0"]
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
               "--max-nodes", str(int(gates["MaximumElements"])), "--max-elements", str(int(gates["MaximumElements"])),
               "--semantic-contract", f"{root}/canonical-semantic.json", "--corner-isotropy-radius", number(tangent),
               "--corner-census", f"{root}/build-census.json", *ownership_report, *etch_options, *basis_options,
               "--prism-tubes", "true", *recipe_options,
               "--maximum-corner-aspect", number(gates["MaximumCornerAspect"]),
               "--minimum-scaled-jacobian", number(gates["MinimumScaledJacobian"]),
               "--maximum-jacobian-condition", number(gates["MaximumJacobianCondition"]),
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
    print(f"DONE {root}", flush=True)


if __name__ == "__main__":
    main()
