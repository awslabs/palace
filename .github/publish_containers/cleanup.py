#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Delete dev-channel containers whose branch was merged to main or deleted.

`Publish Containers` writes a moving `dev-<slug>` channel for a branch on
`workflow_dispatch` or from a `push-containers`-labeled PR. This workflow
removes everything that is no longer needed (e.g., branches merged to main
or deleted).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field

from authorize import DEV_BRANCH_RE, GitHubApi, slugify

# `dev-<slug>-<arch_label>`. Arch labels are [a-z0-9_] with no hyphens (see
# publish.IMAGE_NAME_RE), so the greedy selector group backtracks to the LAST
# hyphen and the split is unambiguous even when the slug contains hyphens
# (`dev-gbozzola-foo-x86_64_v3` -> `dev-gbozzola-foo` + `x86_64_v3`).
DEV_ECR_TAG_RE = re.compile(r"^(?P<selector>dev-[A-Za-z0-9_.-]+)-(?P<arch>[a-z0-9_]+)\Z")

# `dev/<slug>/<name>` — exactly three components, since neither the slug nor the
# object name may contain `/`. That component count is what rejects
# traversal-shaped keys like `dev/../main/x86_64_v3.sif`.
DEV_S3_KEY_RE = re.compile(r"^dev/(?P<slug>[A-Za-z0-9_.-]+)/(?P<name>[A-Za-z0-9_.-]+)\Z")

# S3 keys are opaque strings, so `dev/../x` is a real object inside the dev
# namespace rather than an escape. Reject these anyway: they cannot have come
# from publish.py, and refusing them keeps the plan legible.
_RESERVED_SLUGS = frozenset({".", ".."})


@dataclass(frozen=True)
class Plan:
    """What to delete. `s3` maps bucket name -> keys."""

    ecr_tags: list[str] = field(default_factory=list)
    s3: dict[str, list[str]] = field(default_factory=dict)

    def is_empty(self) -> bool:
        return not self.ecr_tags and not any(self.s3.values())


def ecr_tag_selector(tag: str) -> str | None:
    """Selector an ECR tag belongs to, or None if it is not a dev tag."""
    m = DEV_ECR_TAG_RE.match(tag)
    return m.group("selector") if m else None


def s3_key_selector(key: str) -> str | None:
    """Selector an S3 key belongs to, or None if it is not a dev key."""
    m = DEV_S3_KEY_RE.match(key)
    if not m or m.group("slug") in _RESERVED_SLUGS:
        return None
    return f"dev-{m.group('slug')}"


def orphan_ecr_tags(tags, live) -> list[str]:
    """Dev ECR tags whose selector no live branch claims. Pure."""
    return [t for t in tags if (s := ecr_tag_selector(t)) is not None and s not in live]


def orphan_s3_keys(keys, live) -> list[str]:
    """Dev S3 keys whose selector no live branch claims. Pure."""
    return [k for k in keys if (s := s3_key_selector(k)) is not None and s not in live]


def live_selectors(api) -> set[str]:
    """Selectors a branch still claims.

    A branch claims its selector while it exists, UNLESS its most recent PR was
    merged — then the code is in `main` and `main`'s container supersedes the dev
    channel. "Most recent" is the highest PR number, so reopening a new PR from a
    previously merged branch makes it live again instead of flapping.

    Branches that :data:`DEV_BRANCH_RE` rejects are skipped: they can never have
    published, so letting one mark a selector live would only shield an orphan.
    """
    live = set()
    for branch in api.branches():
        if not DEV_BRANCH_RE.match(branch):
            continue
        pr = api.latest_pr_for_branch(branch)
        if pr is not None and pr[1]:  # (number, merged)
            continue
        live.add(f"dev-{slugify(branch)}")
    return live


def sweep(api, discover_fn) -> Plan:
    """Compute what to delete, guarding against publishes racing the listing.

    ``discover_fn()`` returns ``(ecr_tags, {bucket: keys})`` for everything that
    currently exists. Listing ECR and two S3 buckets is not instantaneous, and
    `Publish Containers` can land a dev channel while it runs, so the live set is
    sampled BEFORE and AFTER discovery and their UNION is protected.

    Union rather than intersection, deliberately: keeping an artifact one run too
    long is reclaimed by the next scheduled sweep, while deleting a channel
    someone just published is not recoverable without a rebuild.
    """
    live_before = live_selectors(api)
    tags, s3_keys = discover_fn()
    live_after = live_selectors(api)
    protected = live_before | live_after

    only_after = live_after - live_before
    if only_after:
        print(f"Protecting selector(s) published during discovery: {', '.join(sorted(only_after))}")
    only_before = live_before - live_after
    if only_before:
        print(f"Deferring selector(s) that went away during discovery: {', '.join(sorted(only_before))}")

    return Plan(
        ecr_tags=orphan_ecr_tags(tags, protected),
        s3={bucket: orphan_s3_keys(keys, protected) for bucket, keys in s3_keys.items()},
    )


class CleanupApi(GitHubApi):
    """GitHubApi plus the branch/PR lookups the sweep needs.

    A subclass rather than new methods on GitHubApi so `authorize.py` — the
    module that gates publishing — is not modified by this change.
    """

    def _lines(self, path: str, jq: str) -> list[str]:
        """`gh api --paginate --jq` output as lines. Raises on failure.

        `--jq` is applied per page and the outputs concatenated, so this needs no
        cross-page JSON assembly. Unlike GitHubApi._get there is no 404 branch:
        every path here is a collection that exists, so any failure is real.
        """
        proc = subprocess.run(
            ["gh", "api", "--paginate", path, "--jq", jq],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            detail = proc.stderr.strip() or f"exit {proc.returncode}"
            raise RuntimeError(f"GitHub API request failed for {path}: {detail}")
        return [ln for ln in proc.stdout.splitlines() if ln.strip()]

    def branches(self) -> list[str]:
        return self._lines(f"repos/{self.repo}/branches", ".[].name")

    def latest_pr_for_branch(self, branch: str) -> tuple[int, bool] | None:
        """(number, merged) for the branch's highest-numbered PR, or None.

        Picking the max number rather than trusting a `sort=` query parameter
        keeps the "latest PR" definition in code where the tests can pin it.
        """
        owner = self.repo.split("/", 1)[0]
        rows = self._lines(
            f"repos/{self.repo}/pulls?head={owner}:{branch}&state=all&per_page=100",
            '.[] | "\\(.number)\\t\\(.merged_at != null)"',
        )
        best: tuple[int, bool] | None = None
        for row in rows:
            number_text, _, merged_text = row.partition("\t")
            number = int(number_text)
            if best is None or number > best[0]:
                best = (number, merged_text == "true")
        return best


def _run(cmd: list[str]) -> None:
    print("+ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def _aws_strings(cmd: list[str]) -> list[str]:
    """Run an `aws --output json` query that yields a list of strings.

    A `--query` that matches nothing prints `null` (or nothing), which becomes an
    empty list — an empty ECR repo or a bucket with no `dev/` objects is normal,
    not an error. Anything else that is not a list of strings is a malformed
    response and must not be silently treated as "nothing to delete".
    """
    proc = subprocess.run(cmd, check=True, capture_output=True, text=True)
    text = proc.stdout.strip()
    if not text:
        return []
    value = json.loads(text)
    if value is None:
        return []
    if not isinstance(value, list) or not all(isinstance(v, str) for v in value):
        raise RuntimeError(f"expected a list of strings from {' '.join(cmd)}, got {value!r}")
    return value


def _chunks(items: list, size: int):
    for i in range(0, len(items), size):
        yield items[i : i + size]


def discover(region: str, ecr_repo: str, buckets: list[str]):
    """Everything that currently exists: (ecr_tags, {bucket: keys}).

    No filtering — :func:`sweep` decides what is an orphan, so this stays a thin
    listing wrapper. The S3 listing is already narrowed to the `dev/` prefix
    because that is all the role may delete; ECR cannot be narrowed server-side,
    so `main-*` and release tags come back here and are dropped by the guard.
    """
    tags = _aws_strings([
        "aws", "ecr", "list-images", "--region", region,
        "--repository-name", ecr_repo, "--filter", "tagStatus=TAGGED",
        "--query", "imageIds[].imageTag", "--output", "json",
    ])

    s3: dict[str, list[str]] = {}
    for bucket in buckets:
        s3[bucket] = _aws_strings([
            "aws", "s3api", "list-objects-v2", "--region", region,
            "--bucket", bucket, "--prefix", "dev/",
            "--query", "Contents[].Key", "--output", "json",
        ])

    return tags, s3


def _describe(plan: Plan) -> str:
    lines = []
    for tag in plan.ecr_tags:
        lines.append(f"  ECR    {tag}")
    for bucket, keys in sorted(plan.s3.items()):
        for key in keys:
            lines.append(f"  S3     s3://{bucket}/{key}")
    return "\n".join(lines) if lines else "  (nothing)"


def assert_dev_only(plan: Plan) -> None:
    """Raise unless every artifact in `plan` is inside the dev namespace.

    The guards already ran in the orphan functions; re-checking at the boundary
    where deletes are actually issued means no non-dev artifact can be deleted
    even if a future caller builds a Plan some other way.
    """
    for tag in plan.ecr_tags:
        if ecr_tag_selector(tag) is None:
            raise ValueError(f"refusing to delete non-dev ECR tag '{tag}'")
    for bucket, keys in plan.s3.items():
        for key in keys:
            if s3_key_selector(key) is None:
                raise ValueError(f"refusing to delete non-dev S3 key '{bucket}/{key}'")


def execute(plan: Plan, region: str, ecr_repo: str) -> None:
    """Delete everything in `plan`, after re-asserting the guards."""
    assert_dev_only(plan)

    # batch-delete-image accepts at most 100 image ids per call.
    for chunk in _chunks(plan.ecr_tags, 100):
        _run([
            "aws", "ecr", "batch-delete-image", "--region", region,
            "--repository-name", ecr_repo,
            "--image-ids", *[f"imageTag={t}" for t in chunk],
        ])

    # delete-objects accepts at most 1000 keys per call.
    for bucket, keys in sorted(plan.s3.items()):
        for chunk in _chunks(keys, 1000):
            payload = json.dumps({"Objects": [{"Key": k} for k in chunk], "Quiet": True})
            _run([
                "aws", "s3api", "delete-objects", "--region", region,
                "--bucket", bucket, "--delete", payload,
            ])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="print the plan but issue no deletions",
    )
    args = parser.parse_args(argv)

    repo = os.environ["OWNER_REPO"]
    region = os.environ["AWS_REGION"]
    ecr_repo = os.environ["ECR_REPOSITORY"]
    buckets = [os.environ["S3_CONTAINERS_BUCKET"], os.environ["S3_SCHEMAS_BUCKET"]]

    plan = sweep(
        CleanupApi(repo),
        lambda: discover(region, ecr_repo, buckets),
    )

    # Always print the full plan, so every run leaves an auditable record of what
    # it decided even when it deletes nothing.
    print("Orphaned dev artifacts:")
    print(_describe(plan))

    if plan.is_empty():
        print("Nothing to delete.")
        return 0
    if args.dry_run:
        print("Dry run; no deletions issued.")
        return 0

    try:
        execute(plan, region, ecr_repo)
    except ValueError as exc:
        print(f"::error::{exc}", file=sys.stderr)
        return 1
    print(
        f"Deleted {len(plan.ecr_tags)} ECR tag(s) and "
        f"{sum(len(v) for v in plan.s3.values())} S3 object(s)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
