#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Preflight matrix for one geometry: ranks x libraries x uniform refinement x cracking.

Every cell runs the geometry-only preflight (palace --surface-response-preflight), then the
identification audit against the mesh; the summary compares the canonical manifest digests
across cells (A3 library independence, A4 rank determinism, A5 refinement invariance, and
the CrackInternalBoundaryElements dependence).

    python3 -m surface_response_identification.preflight_matrix --mesh M.msh2 \\
        --ground 5 6 7 --terminal 9 --sa 8 --library seed=seed.json --library full=full.json \\
        --ranks 1 2 4 6 --uniform-levels 0 1 --crack true false --output DIR

Local limits (decision 69 incident): ranks <= 6, no oversubscription; the runner refuses
more ranks than cores.
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import time

if __package__ in (None, ""):
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    __package__ = "surface_response_identification"

from . import audit, manifest as M  # noqa: E402
from .preflight_config import preflight_config  # noqa: E402

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DEFAULT_PALACE = os.path.join(REPO, "build", "bin", "palace")


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def run_cell(palace, config_path, ranks, log_path, timeout):
    command = [palace, "-np", str(ranks), "--surface-response-preflight", config_path]
    environment = dict(os.environ, OMP_NUM_THREADS="1")
    for name in list(environment):
        if name.startswith("PALACE_RESPONSE_"):
            del environment[name]
    started = time.time()
    with open(log_path, "w") as log:
        completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=environment, timeout=timeout)
    return completed.returncode, time.time() - started


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mesh", required=True)
    parser.add_argument("--ground", type=int, nargs="+", required=True)
    parser.add_argument("--terminal", type=int, nargs="+", action="append", default=[])
    parser.add_argument("--sa", type=int, nargs="*", default=[])
    parser.add_argument("--ms", type=int, nargs="*")
    parser.add_argument("--ma", type=int, nargs="*")
    parser.add_argument("--substrate", type=int, nargs="+", default=[1])
    parser.add_argument("--vacuum", type=int, nargs="+", default=[2])
    parser.add_argument("--library", action="append", required=True, help="label=path (repeatable)")
    parser.add_argument("--ranks", type=int, nargs="+", default=[1])
    parser.add_argument("--uniform-levels", type=int, nargs="+", default=[0])
    parser.add_argument("--crack", nargs="+", default=["true"], choices=["true", "false"])
    parser.add_argument("--l0", type=float, default=1.0e-6)
    parser.add_argument("--palace", default=DEFAULT_PALACE)
    parser.add_argument("--timeout", type=float, default=1100.0, help="seconds per cell")
    parser.add_argument("--max-local-ranks", type=int, default=6)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    cores = os.cpu_count() or 1
    for ranks in args.ranks:
        if ranks > args.max_local_ranks or ranks > cores:
            raise SystemExit(f"{ranks} ranks exceeds the local limit ({args.max_local_ranks}, cores {cores})")
    libraries = dict(entry.split("=", 1) for entry in args.library)
    os.makedirs(args.output, exist_ok=True)
    cells = []
    for label, library in libraries.items():
        for levels in args.uniform_levels:
            for crack in args.crack:
                for ranks in args.ranks:
                    name = f"{label}-u{levels}-crack{crack}-np{ranks}"
                    directory = os.path.join(args.output, name)
                    os.makedirs(directory, exist_ok=True)
                    config = preflight_config(
                        args.mesh, args.ground, args.terminal, args.sa, library, os.path.join(directory, "postpro"), ms=args.ms, ma=args.ma, l0=args.l0,
                        uniform_levels=levels, crack=crack == "true", substrate_attributes=args.substrate, vacuum_attributes=args.vacuum,
                    )
                    config_path = os.path.join(directory, "config.json")
                    with open(config_path, "w") as target:
                        json.dump(config, target, indent=2)
                    log_path = os.path.join(directory, "palace.log")
                    manifest_path = os.path.join(directory, "postpro", "surface-response-requirements.json")
                    cell = {"Name": name, "Library": label, "UniformLevels": levels, "Crack": crack == "true", "Ranks": ranks, "Directory": directory}
                    if os.path.exists(manifest_path) and not os.environ.get("PREFLIGHT_MATRIX_RERUN"):
                        cell["ExitCode"] = 0
                        cell["Seconds"] = None
                        cell["Reused"] = True
                    else:
                        try:
                            cell["ExitCode"], cell["Seconds"] = run_cell(args.palace, config_path, ranks, log_path, args.timeout)
                        except subprocess.TimeoutExpired:
                            cell["ExitCode"], cell["Seconds"] = "timeout", args.timeout
                    if os.path.exists(manifest_path):
                        manifest = M.load_manifest(manifest_path)
                        cell["ManifestSha256"] = sha256(manifest_path)
                        cell["DigestFull"] = M.canonical_digest(manifest)
                        cell["DigestGeometryOnly"] = M.canonical_digest(manifest, geometry_only_counts=True)
                        cell["GeometryDigest"] = manifest.get("Identification", {}).get("GeometryDigest")
                        cell["Summary"] = manifest.get("Summary")
                        cell["Statistics"] = manifest.get("Statistics", {}).get("Geometry")
                        cell["Log"] = M.parse_palace_log(log_path) if os.path.exists(log_path) else None
                        if levels == 0:
                            audit_args = argparse.Namespace(
                                mesh=args.mesh, config=config_path, manifest=manifest_path, log=log_path if os.path.exists(log_path) else None, compare=None, radius=None,
                                corner_tolerance=audit.P.CORNER_ANGLE_TOLERANCE_DEGREES, output_prefix=os.path.join(directory, "audit"),
                            )
                            result = audit.run_audit(audit_args)
                            with open(os.path.join(directory, "audit.json"), "w") as target:
                                json.dump(result, target, indent=1, default=str)
                            with open(os.path.join(directory, "audit.md"), "w") as target:
                                target.write(audit.render_markdown(result))
                            cell["Gates"] = {g["Gate"]: g["Status"] for g in result["Gates"]}
                            cell["GapBound"] = {k: v for k, v in result["GapBound"].items() if k != "OmittedByClassFromLog"}
                        else:
                            cell["Gates"] = "audit skipped on the refined mesh (the audit reads the unrefined MSH); compare digests"
                    else:
                        cell["Error"] = "no manifest"
                    cells.append(cell)
                    print(json.dumps({k: cell.get(k) for k in ("Name", "ExitCode", "Seconds", "DigestFull", "Gates")}, default=str), flush=True)

    # Cross-cell comparisons.
    by_name = {c["Name"]: c for c in cells}
    comparisons = []
    reference = next((c for c in cells if c.get("DigestFull")), None)
    for cell in cells:
        if not cell.get("DigestFull") or cell is reference:
            continue
        a = M.load_manifest(os.path.join(reference["Directory"], "postpro", "surface-response-requirements.json"))
        b = M.load_manifest(os.path.join(cell["Directory"], "postpro", "surface-response-requirements.json"))
        diff = M.diff_manifests(a, b)
        diff_geometry = M.diff_manifests(a, b, geometry_only_counts=True)
        # Version 2: the identification's GeometryDigest is the A3 / A4 / A5 identity.
        digest_identical = reference.get("GeometryDigest") == cell.get("GeometryDigest") if reference.get("GeometryDigest") else None
        comparisons.append(
            {
                "Reference": reference["Name"],
                "Cell": cell["Name"],
                "GeometryDigestIdentical": digest_identical,
                "FullIdentical": diff["Identical"],
                "GeometryOnlyIdentical": diff_geometry["Identical"],
                "Added": len(diff["Added"]),
                "Removed": len(diff["Removed"]),
                "Changed": len(diff["Changed"]),
                "AddedGeometry": [e["Topology"] for e in diff["Added"]],
                "RemovedGeometry": [e["Topology"] for e in diff["Removed"]],
                "ChangedGeometry": [(e["A"]["Topology"], e["A"]["Count"], e["B"]["Count"], round(e["A"]["TotalEdgeLength"], 6), round(e["B"]["TotalEdgeLength"], 6)) for e in diff["Changed"]],
                "Diff": diff if not diff["Identical"] else None,
            }
        )
    summary = {
        "Mesh": os.path.abspath(args.mesh),
        "MeshSha256": sha256(args.mesh),
        "Palace": args.palace,
        "Libraries": {k: {"Path": os.path.abspath(v), "Sha256": sha256(v)} for k, v in libraries.items()},
        "Cells": cells,
        "Comparisons": comparisons,
    }
    with open(os.path.join(args.output, "matrix.json"), "w") as target:
        json.dump(summary, target, indent=1, default=str)
    lines = [f"# Preflight matrix: {os.path.basename(args.mesh)}", "", f"mesh sha256 {summary['MeshSha256'][:16]}", ""]
    lines.append("| Cell | exit | s | manifest sha | digest full | digest geom | Exact/Missing | gates |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for c in cells:
        s = c.get("Summary") or {}
        counts = s.get("Counts", {})
        gates = c.get("Gates")
        gate_text = ", ".join(f"{k.split('-', 1)[-1]}={v[0]}" for k, v in gates.items()) if isinstance(gates, dict) else str(gates)[:40]
        lines.append(f"| {c['Name']} | {c.get('ExitCode')} | {c.get('Seconds') if c.get('Seconds') is None else round(c['Seconds'], 1)} | {c.get('ManifestSha256', '')[:12]} | {(c.get('GeometryDigest') or c.get('DigestFull') or '')[:12]} | {c.get('DigestGeometryOnly', '')[:12]} | {counts.get('Exact')}/{counts.get('Missing')} | {gate_text} |")
    lines.append("")
    lines.append("| Reference | Cell | geometry digest identical | full identical | geometry-only identical | added | removed | changed |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for c in comparisons:
        lines.append(f"| {c['Reference']} | {c['Cell']} | {c['GeometryDigestIdentical']} | {c['FullIdentical']} | {c['GeometryOnlyIdentical']} | {c['Added']} {c['AddedGeometry']} | {c['Removed']} {c['RemovedGeometry']} | {c['Changed']} {c['ChangedGeometry'][:6]} |")
    with open(os.path.join(args.output, "matrix.md"), "w") as target:
        target.write("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
