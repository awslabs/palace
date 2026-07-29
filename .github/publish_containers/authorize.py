#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Authorization decision for the trusted container publish workflow.

The publish workflow runs from the default branch (via ``workflow_run``), so
this code — not the PR-controlled build — decides whether a build may publish
and to which selector. The decision uses only trusted inputs: the event facts
GitHub reports for the build run, plus live lookups against the GitHub API. The
build passes no authorization data of its own.

Two entry points, one decision function:

- default: emit ``authorized=`` / ``selector=`` to ``$GITHUB_OUTPUT``.
- ``--reverify EXPECTED``: recompute the decision immediately before the writes
  (TOCTOU guard) and fail unless it still authorizes the same selector.

A build is authorized to publish only when, checked against the live GitHub
API:
- it did NOT come from a fork (head repo == base repo); and
- the triggering event is one of:
  - push to `main`, and `main` still points at the built commit
    -> publishes the `main` channel;
  - push of a `vX.Y.Z` release tag that still points at the built commit
    -> publishes the `X.Y.Z` channel;
  - manual `workflow_dispatch` on a branch that still points at the built
    commit -> publishes `dev-<branch>`;
  - a same-repo `pull_request` whose head (matched by sha AND ref) is an OPEN
    PR currently carrying the `push-containers` label -> publishes
    `dev-<branch>`.

For the two `dev-<branch>` paths the branch name must match `DEV_BRANCH_RE`.

The GitHub API is reached through a small seam (the `GitHubApi` class) so the
decision logic is unit-tested with a fake client and no network.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass

PUSH_CONTAINERS_LABEL = "push-containers"
RELEASE_TAG_RE = re.compile(r"^v[0-9]+\.[0-9]+\.[0-9]+\Z")  # \Z, not $: see DEV_BRANCH_RE

# Dev-channel (dispatch / labeled-PR) branch names must match this shape: a
# hyphen-free, slash-free prefix, exactly one "/", then a segment with no
# further "/" (e.g., `gbozzola/test-containers`). The reason for this is that
# ECR has restrictions on what characters are allowed and we sanitize them in
# slugify in such a way that could lead to collisions. We could remove this
# restriction by using hashes adding much more infrastructure, but, honestly,
# YAGNI.
#
# \Z (not $) so a trailing newline can't sneak past the anchor.
DEV_BRANCH_RE = re.compile(r"^[A-Za-z0-9_.]+/[A-Za-z0-9][A-Za-z0-9_.-]*\Z")


@dataclass(frozen=True)
class Facts:
    """Trusted facts about the build run, from the ``workflow_run`` payload."""

    repo: str  # owner/repo, e.g. "awslabs/palace"
    event: str  # "push" | "workflow_dispatch" | "pull_request" | ...
    head_sha: str
    head_branch: str  # branch or tag name (refs/heads/ and refs/tags/ stripped)
    head_repo: str  # head_repository.full_name
    base_repo: str  # repository.full_name


@dataclass(frozen=True)
class Decision:
    authorized: bool
    selector: str
    reason: str


def _reject(reason: str) -> Decision:
    return Decision(False, "", reason)


def _unsupported_branch_reason(branch: str) -> str:
    return (
        f"branch '{branch}' is not supported for dev publishing yet: the name "
        f"must look like 'prefix/name' (a hyphen-free prefix, one '/', then a "
        f"name of [A-Za-z0-9_.-] starting alphanumeric), e.g. 'user/my-feature'. "
        f"This restriction keeps the dev-<branch> selector collision-free. "
        f"main and release-tag publishing are unaffected."
    )


def slugify(branch: str) -> str:
    """Collapse a branch name to an ECR-tag / S3-prefix-safe slug.

    ECR tags allow only ``[A-Za-z0-9_.-]`` and ``/`` is illegal, so everything
    else becomes ``-``. ``dev-<branch>`` is intentionally a MOVING selector.
    """
    return re.sub(r"[^A-Za-z0-9_.-]", "-", branch)


class GitHubApi:
    """Live GitHub API access via the ``gh`` CLI (authenticated by GH_TOKEN)."""

    def __init__(self, repo: str) -> None:
        self.repo = repo

    def _get(self, path: str) -> object | None:
        """Decoded JSON body, or None if the resource does not exist (404).

        Only an explicit 404 means "absent". Every other failure — auth, rate
        limit, 5xx, transport, non-JSON body — raises, because silently reading
        it as "absent" would turn an outage into an unauthorized publish: a
        release could then finish green having published nothing.
        """
        # --include prints the status line and headers before the body, which is
        # the only way to tell 404 apart from the rest (gh exits 1 for all of
        # them).
        proc = subprocess.run(
            ["gh", "api", "--include", path],
            capture_output=True,
            text=True,
        )
        headers, _, body = proc.stdout.replace("\r\n", "\n").partition("\n\n")
        match = re.match(r"^HTTP/\S+ ([0-9]{3})(?:\s|\Z)", headers)
        if match is None:
            detail = proc.stderr.strip() or "missing HTTP status"
            raise RuntimeError(f"GitHub API request failed for {path}: {detail}")
        status = int(match.group(1))
        if status == 404:
            return None
        if proc.returncode != 0 or not 200 <= status < 300:
            raise RuntimeError(f"GitHub API request failed for {path}: HTTP {status}")
        try:
            return json.loads(body)
        except json.JSONDecodeError as exc:
            raise RuntimeError(f"GitHub API returned invalid JSON for {path}") from exc

    def ref_sha(self, kind: str, name: str) -> str | None:
        """Commit sha a ref points at, or None. Peels annotated tags.

        ``kind`` is "heads" or "tags". For an annotated tag, the ref object is
        the tag object, so we peel it to the commit it targets.
        """
        ref = self._get(f"repos/{self.repo}/git/ref/{kind}/{name}")
        if not isinstance(ref, dict):
            return None
        obj = ref.get("object") or {}
        sha, obj_type = obj.get("sha"), obj.get("type")
        if obj_type == "tag":
            tag = self._get(f"repos/{self.repo}/git/tags/{sha}")
            if not isinstance(tag, dict):
                return None
            sha = (tag.get("object") or {}).get("sha")
        return sha

    def pr_for_head(self, sha: str, ref: str) -> tuple[int, str, str] | None:
        """(PR number, head repo full_name, state) for the PR whose head is
        ``sha`` AND whose head branch is ``ref``.

        ``state`` is "open" or "closed".
        """
        pulls = self._get(f"repos/{self.repo}/commits/{sha}/pulls")
        if not isinstance(pulls, list):
            return None
        for pr in pulls:
            if not isinstance(pr, dict):
                continue
            head = pr.get("head") or {}
            # ``/commits/{sha}/pulls`` returns PRs that merely *contain* the
            # commit (filtered out by matching ``head.sha``), and multiple PRs
            # can share a head commit (e.g. two branches at the same sha), so
            # we ALSO match ``head.ref`` to pin the PR for the branch that
            # actually built, not some other PR that happens to sit on the same
            # commit.
            if head.get("sha") == sha and head.get("ref") == ref:
                repo = (head.get("repo") or {}).get("full_name") or ""
                state = pr.get("state") or ""
                number = pr.get("number")
                if number is None:
                    return None
                return int(number), repo, state
        return None

    def pr_has_label(self, number: int, label: str) -> bool:
        labels = self._get(f"repos/{self.repo}/issues/{number}/labels")
        if not isinstance(labels, list):
            return False
        return any(isinstance(lbl, dict) and lbl.get("name") == label for lbl in labels)


def decide(facts, api) -> Decision:
    """Decide whether ``facts`` authorizes a publish, and to which selector.

    ``api`` is anything providing ref_sha(kind, name), pr_for_head(sha, ref),
    and pr_has_label(number, label) — the live GitHubApi in production, a fake
    in tests. Every branch re-derives from trusted facts / live API state, so
    calling ``decide`` again just before the writes is a valid TOCTOU re-check.
    """
    # Forks never publish, regardless of event. A maintainer approving a fork
    # build authorizes EXECUTING its code on a runner, never publishing.
    if facts.head_repo != facts.base_repo:
        return _reject(f"head repo '{facts.head_repo}' != base repo '{facts.base_repo}' (fork)")

    if facts.event == "push":
        # head_branch strips both refs/heads/ and refs/tags/, so it is
        # ambiguous; resolve each ref kind explicitly and require it to point at
        # the built commit right now.
        if facts.head_branch == "main":
            if api.ref_sha("heads", "main") == facts.head_sha:
                return Decision(True, "main", f"push to main at {facts.head_sha}")
            return _reject("main does not currently point at the built sha")
        if RELEASE_TAG_RE.match(facts.head_branch):
            if api.ref_sha("tags", facts.head_branch) == facts.head_sha:
                selector = facts.head_branch[1:]  # drop leading "v"
                return Decision(True, selector, f"release tag {facts.head_branch} at {facts.head_sha}")
            return _reject(f"release tag {facts.head_branch} does not point at the built sha")
        return _reject(f"push ref '{facts.head_branch}' is neither main nor a vX.Y.Z release tag")

    if facts.event == "workflow_dispatch":
        # Dispatch = publish request for a branch. GitHub gates dispatch to
        # write-access users; require the ref to be a branch at the built sha.
        if api.ref_sha("heads", facts.head_branch) != facts.head_sha:
            return _reject(f"dispatch ref '{facts.head_branch}' is not a branch at the built sha")
        if not DEV_BRANCH_RE.match(facts.head_branch):
            return _reject(_unsupported_branch_reason(facts.head_branch))
        return Decision(True, f"dev-{slugify(facts.head_branch)}", f"dispatch on branch {facts.head_branch} at {facts.head_sha}")

    if facts.event == "pull_request":
        pr = api.pr_for_head(facts.head_sha, facts.head_branch)
        if pr is None:
            return _reject(f"no PR has head {facts.head_branch} at {facts.head_sha}")
        number, pr_head_repo, state = pr
        if pr_head_repo != facts.base_repo:
            return _reject(f"PR #{number} head repo '{pr_head_repo}' is a fork")
        if state != "open":
            return _reject(f"PR #{number} is not open (state '{state}')")
        if not api.pr_has_label(number, PUSH_CONTAINERS_LABEL):
            return _reject(f"PR #{number} does not currently carry the {PUSH_CONTAINERS_LABEL} label")
        if not DEV_BRANCH_RE.match(facts.head_branch):
            return _reject(_unsupported_branch_reason(facts.head_branch))
        return Decision(True, f"dev-{slugify(facts.head_branch)}", f"labeled same-repo PR #{number} at {facts.head_sha}")

    # Unknown / unsupported event → fail closed.
    return _reject(f"unsupported source event '{facts.event}'")


def facts_from_env() -> Facts:
    return Facts(
        repo=os.environ["OWNER_REPO"],
        event=os.environ["SOURCE_EVENT"],
        head_sha=os.environ["HEAD_SHA"],
        head_branch=os.environ["HEAD_BRANCH"],
        head_repo=os.environ["HEAD_REPO"],
        base_repo=os.environ["BASE_REPO"],
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reverify",
        metavar="EXPECTED_SELECTOR",
        help="TOCTOU re-check: recompute the decision and fail unless it still "
        "authorizes this exact selector.",
    )
    args = parser.parse_args(argv)

    facts = facts_from_env()
    decision = decide(facts, GitHubApi(facts.repo))

    if args.reverify is not None:
        expected = args.reverify
        if not decision.authorized or decision.selector != expected:
            print(
                f"::error::Re-verification failed for selector '{expected}' "
                f"({facts.event} at {facts.head_sha}): {decision.reason}; "
                f"aborting before any publish.",
                file=sys.stderr,
            )
            return 1
        print(f"Re-verified '{expected}' still current at {facts.head_sha}; proceeding to publish.")
        return 0

    authorized = "true" if decision.authorized else "false"
    print(f"Authorization: authorized={authorized} selector='{decision.selector}' ({decision.reason})")
    out_path = os.environ.get("GITHUB_OUTPUT")
    if out_path:
        with open(out_path, "a", encoding="utf-8") as fh:
            fh.write(f"authorized={authorized}\n")
            fh.write(f"selector={decision.selector}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
