#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Hash-bound Linux source bundle of the decision-62(4) executable (the recorded
coupon-accuracy-assessment-20260913/experiments/freeze_linux_source.py procedure,
parametrized by the checkout and the output root; no dependency installs / builds).

usage: freeze_linux_source_perf62.py REPO OUT_ROOT [--name source-freeze-perf62]
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

TOOLS = ("estimate_archived_fields.py", "test_estimate_archived_fields.py", "run_bounded_mesher.py")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("repo", type=Path)
    parser.add_argument("out_root", type=Path)
    parser.add_argument("--name", default="source-freeze-perf62")
    parser.add_argument("--schema", type=Path, default=None,
                        help="embedded_schema.hpp of the checkout's build (default REPO/build/generated/embedded_schema.hpp)")
    args = parser.parse_args()
    repo = args.repo.resolve()
    out = args.out_root.resolve() / args.name
    assert not out.exists(), out
    out.mkdir(parents=True)
    for src in sorted(p for p in (repo / "palace").rglob("*") if p.is_file()):
        dst = out / src.relative_to(repo)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
    (out / "generated").mkdir()
    shutil.copyfile(args.schema or repo / "build/generated/embedded_schema.hpp", out / "generated/embedded_schema.hpp")
    (out / "tools").mkdir()
    for name in TOOLS:
        shutil.copyfile(repo / "examples/cpw3d_surface/spatial_coupon" / name, out / "tools" / name)
    (out / "tracked.diff").write_bytes(subprocess.check_output(["git", "diff"], cwd=repo))
    (out / "status.txt").write_bytes(subprocess.check_output(["git", "status", "--short"], cwd=repo))
    assert not subprocess.check_output(["git", "diff", "--cached"], cwd=repo)
    hashes = {str(p.relative_to(out)): hashlib.sha256(p.read_bytes()).hexdigest()
              for p in sorted(out.rglob("*")) if p.is_file()}
    (out / "manifest.json").write_text(json.dumps({
        "HEAD": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repo, text=True).strip(),
        "Branch": subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=repo, text=True).strip(),
        "Files": hashes, "NoStagedFiles": True}, indent=2) + "\n")
    tar_path = out.parent / f"{args.name}.tar.gz"
    with tarfile.open(tar_path, "x:gz") as tar:
        tar.add(out, arcname=args.name)
    print(hashlib.sha256(tar_path.read_bytes()).hexdigest(), tar_path)


if __name__ == "__main__":
    main()
