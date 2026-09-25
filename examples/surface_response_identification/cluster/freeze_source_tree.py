#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Hash-bound freeze of the CMake source tree needed to build the palace executable.

The freeze carries every file `palace/CMakeLists.txt` reaches at configure time: `palace/`
(sources, `palace/cmake/`), `cmake/` (EmbedSchema), `scripts/schema/` (embedded JSON schemas)
`test/unit/` (`add_subdirectory(../test/unit)` at configure time; the unit tests are not
built) and this `cluster/` recipe itself (so the freeze carries the build script that consumes it). Working-tree contents are frozen (tracked files plus untracked, non-ignored files
under those directories), the tracked diff and status are recorded, and `manifest.json`
carries HEAD, `git describe`, and the SHA-256 of every frozen file. The tarball is what the
cluster build (`build_from_cmake_tree.py`) verifies before configuring.

    python3 freeze_source_tree.py REPO OUT_ROOT --name source-freeze-<label>
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tarfile

FROZEN_DIRECTORIES = ("palace", "cmake", "scripts/schema", "test/unit", "examples/surface_response_identification/cluster")


def git(repo, *args):
    return subprocess.check_output(["git", *args], cwd=repo, text=True)


def frozen_files(repo):
    tracked = git(repo, "ls-files", "-z", "--", *FROZEN_DIRECTORIES).split("\0")
    untracked = git(repo, "ls-files", "-z", "--others", "--exclude-standard", "--", *FROZEN_DIRECTORIES).split("\0")
    files = sorted({p for p in tracked + untracked if p and (repo / p).is_file()})
    return files, sorted(p for p in untracked if p)


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("repo", type=Path)
    parser.add_argument("out_root", type=Path)
    parser.add_argument("--name", required=True, help="freeze directory / tarball name, e.g. source-freeze-sct002")
    args = parser.parse_args()
    repo = args.repo.resolve()
    out = args.out_root.resolve() / args.name
    if out.exists() or out.with_suffix(".tar.gz").exists():
        raise SystemExit(f"{out} or its tarball exists: refusing to overwrite a freeze")
    if git(repo, "diff", "--cached", "--name-only").strip():
        raise SystemExit("staged changes present: commit or unstage before freezing")
    files, untracked = frozen_files(repo)
    out.mkdir(parents=True)
    for relative in files:
        target = out / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(repo / relative, target)
    (out / "tracked.diff").write_text(git(repo, "diff", "--", *FROZEN_DIRECTORIES))
    (out / "status.txt").write_text(git(repo, "status", "--short", "--", *FROZEN_DIRECTORIES))
    describe = git(repo, "describe", "--tags", "--always", "--dirty").strip()
    manifest = {
        "HEAD": git(repo, "rev-parse", "HEAD").strip(),
        "Branch": git(repo, "rev-parse", "--abbrev-ref", "HEAD").strip(),
        "Describe": describe,
        "FrozenDirectories": list(FROZEN_DIRECTORIES),
        "UntrackedFiles": untracked,
        "TrackedDiffEmpty": not (out / "tracked.diff").read_text().strip(),
        "Files": {p: sha256(out / p) for p in files},
        "NoStagedFiles": True,
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    tar_path = out.parent / f"{args.name}.tar.gz"
    with tarfile.open(tar_path, "x:gz") as tar:
        tar.add(out, arcname=args.name)
    print(json.dumps({"Tarball": str(tar_path), "SHA256": sha256(tar_path), "HEAD": manifest["HEAD"],
                      "Describe": describe, "Files": len(files), "TrackedDiffEmpty": manifest["TrackedDiffEmpty"]}))


if __name__ == "__main__":
    main()
