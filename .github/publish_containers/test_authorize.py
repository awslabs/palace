#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests for the container-publish authorization decision.

Run: python3 -m unittest discover -s .github/publish_containers -p 'test_*.py'

The decision logic is exercised with a FakeApi returning canned responses, so
these run offline and deterministically. Cases mirror the real GitHub API
behaviors verified live during development (annotated-tag peeling, PRs that
merely contain a commit vs. head it, empty pull lists, missing refs).
"""

from __future__ import annotations

import unittest

from authorize import Decision, Facts, decide, slugify

SHA = "c45e4ab4f2356589ba40f4fde2403f5854282826"
OTHER_SHA = "deadbeefdeadbeefdeadbeefdeadbeefdeadbeef"
REPO = "awslabs/palace"


class FakeApi:
    """Stand-in for GitHubApi. All lookups come from constructor data."""

    def __init__(self, *, refs=None, pr_for_head=None, labels=None):
        # refs: {("heads"|"tags", name): commit_sha}
        self._refs = refs or {}
        # pr_for_head: {(head_sha, head_ref): (number, head_repo_full_name, state)}
        self._pr = pr_for_head or {}
        # labels: {pr_number: set(label names)}
        self._labels = labels or {}

    def ref_sha(self, kind, name):
        return self._refs.get((kind, name))

    def pr_for_head(self, sha, ref):
        return self._pr.get((sha, ref))

    def pr_has_label(self, number, label):
        return label in self._labels.get(number, set())


def facts(**kw):
    base = dict(
        repo=REPO,
        event="push",
        head_sha=SHA,
        head_branch="main",
        head_repo=REPO,
        base_repo=REPO,
    )
    base.update(kw)
    return Facts(**base)


class ForkRejection(unittest.TestCase):
    def test_fork_head_repo_rejected_regardless_of_event(self):
        for event in ("push", "workflow_dispatch", "pull_request"):
            d = decide(facts(event=event, head_repo="attacker/palace"), FakeApi())
            self.assertFalse(d.authorized)
            self.assertIn("fork", d.reason)


class PushMain(unittest.TestCase):
    def test_main_at_built_sha_authorizes(self):
        api = FakeApi(refs={("heads", "main"): SHA})
        d = decide(facts(head_branch="main", head_sha=SHA), api)
        self.assertEqual(d, Decision(True, "main", d.reason))

    def test_stale_main_rejected(self):
        api = FakeApi(refs={("heads", "main"): OTHER_SHA})
        d = decide(facts(head_branch="main", head_sha=SHA), api)
        self.assertFalse(d.authorized)
        self.assertEqual(d.selector, "")


class ReleaseTag(unittest.TestCase):
    def test_release_tag_authorizes_and_strips_v(self):
        api = FakeApi(refs={("tags", "v0.18.0"): SHA})
        d = decide(facts(head_branch="v0.18.0", head_sha=SHA), api)
        self.assertTrue(d.authorized)
        self.assertEqual(d.selector, "0.18.0")

    def test_release_tag_wrong_sha_rejected(self):
        api = FakeApi(refs={("tags", "v0.18.0"): OTHER_SHA})
        d = decide(facts(head_branch="v0.18.0", head_sha=SHA), api)
        self.assertFalse(d.authorized)

    def test_nonexistent_tag_rejected(self):
        d = decide(facts(head_branch="v9.9.9", head_sha=SHA), FakeApi())
        self.assertFalse(d.authorized)

    def test_bad_syntax_tag_rejected(self):
        # ref exists at the sha but name is not strict vX.Y.Z
        api = FakeApi(refs={("tags", "v1.2"): SHA})
        d = decide(facts(head_branch="v1.2", head_sha=SHA), api)
        self.assertFalse(d.authorized)

    def test_feature_branch_push_rejected(self):
        d = decide(facts(head_branch="some-feature", head_sha=SHA), FakeApi())
        self.assertFalse(d.authorized)


class Dispatch(unittest.TestCase):
    def test_branch_at_built_sha_authorizes_moving_dev_selector(self):
        api = FakeApi(refs={("heads", "my/feature"): SHA})
        d = decide(facts(event="workflow_dispatch", head_branch="my/feature", head_sha=SHA), api)
        self.assertTrue(d.authorized)
        self.assertEqual(d.selector, "dev-my-feature")  # slash slugified

    def test_dispatch_on_tag_ref_rejected(self):
        # A tag isn't a branch, so heads/<name> resolves to nothing.
        api = FakeApi(refs={("tags", "v0.18.0"): SHA})
        d = decide(facts(event="workflow_dispatch", head_branch="v0.18.0", head_sha=SHA), api)
        self.assertFalse(d.authorized)

    def test_dispatch_stale_branch_rejected(self):
        api = FakeApi(refs={("heads", "team/b"): OTHER_SHA})
        d = decide(facts(event="workflow_dispatch", head_branch="team/b", head_sha=SHA), api)
        self.assertFalse(d.authorized)

    def test_dispatch_unsupported_branch_shape_rejected(self):
        # Branch is real and at the built sha, but its name is not the
        # collision-free prefix/name shape -> unsupported for dev publishing.
        # Includes "team/b\n": $ would accept it (matching before the newline)
        # and it slug-collides with "team/b-"; \Z rejects it.
        for bad in ("myfeature", "feature-a/b", "a/b/c", "/leading", "team/-dash", "team/b\n"):
            api = FakeApi(refs={("heads", bad): SHA})
            d = decide(facts(event="workflow_dispatch", head_branch=bad, head_sha=SHA), api)
            self.assertFalse(d.authorized, f"{bad!r} should be rejected")
            self.assertIn("not supported", d.reason)


class PullRequest(unittest.TestCase):
    # FakeApi.pr_for_head is keyed on (head_sha, head_ref); facts' head_branch
    # is the ref, so the key must be (SHA, <branch>).
    def test_labeled_open_same_repo_pr_authorizes(self):
        api = FakeApi(pr_for_head={(SHA, "team/feat"): (840, REPO, "open")}, labels={840: {"push-containers"}})
        d = decide(facts(event="pull_request", head_branch="team/feat", head_sha=SHA), api)
        self.assertTrue(d.authorized)
        self.assertEqual(d.selector, "dev-team-feat")

    def test_labeled_pr_unsupported_branch_shape_rejected(self):
        # Passes fork/open/label, but the branch name is not the supported shape.
        api = FakeApi(pr_for_head={(SHA, "flat-name"): (843, REPO, "open")}, labels={843: {"push-containers"}})
        d = decide(facts(event="pull_request", head_branch="flat-name", head_sha=SHA), api)
        self.assertFalse(d.authorized)
        self.assertIn("not supported", d.reason)

    def test_unlabeled_pr_rejected(self):
        api = FakeApi(pr_for_head={(SHA, "team/x"): (841, REPO, "open")}, labels={841: {"no-long-tests"}})
        d = decide(facts(event="pull_request", head_branch="team/x", head_sha=SHA), api)
        self.assertFalse(d.authorized)
        self.assertIn("label", d.reason)

    def test_closed_but_labeled_pr_rejected(self):
        # A closed (or merged) PR must not publish even if it still carries the label.
        api = FakeApi(pr_for_head={(SHA, "team/x"): (842, REPO, "closed")}, labels={842: {"push-containers"}})
        d = decide(facts(event="pull_request", head_branch="team/x", head_sha=SHA), api)
        self.assertFalse(d.authorized)
        self.assertIn("not open", d.reason)

    def test_no_pr_heads_this_sha_rejected(self):
        # e.g. the sha is a main commit or only a member of a PR, not its head
        d = decide(facts(event="pull_request", head_branch="team/x", head_sha=SHA), FakeApi(pr_for_head={}))
        self.assertFalse(d.authorized)
        self.assertIn("no PR", d.reason)

    def test_fork_pr_via_head_repo_rejected(self):
        # Same-repo base check passed, but the PR's own head repo is a fork.
        api = FakeApi(pr_for_head={(SHA, "team/x"): (900, "attacker/palace", "open")}, labels={900: {"push-containers"}})
        d = decide(facts(event="pull_request", head_branch="team/x", head_sha=SHA), api)
        self.assertFalse(d.authorized)
        self.assertIn("fork", d.reason)

    def test_shared_commit_picks_pr_for_the_built_branch(self):
        # Two PRs share the built commit: an authorized labeled PR on the branch
        # that built (team/real), and an UNLABELED PR on another branch at the
        # same sha (team/evil). Matching on ref must pick the built branch's PR,
        # not inherit the other's authorization. Here the build is team/evil
        # (unlabeled) and must be rejected despite team/real carrying the label.
        api = FakeApi(
            pr_for_head={
                (SHA, "team/real"): (10, REPO, "open"),
                (SHA, "team/evil"): (11, REPO, "open"),
            },
            labels={10: {"push-containers"}},  # only the real branch's PR is labeled
        )
        d = decide(facts(event="pull_request", head_branch="team/evil", head_sha=SHA), api)
        self.assertFalse(d.authorized)
        self.assertIn("label", d.reason)  # checked #11's (absent) label, not #10's


class UnknownEvent(unittest.TestCase):
    def test_schedule_fails_closed(self):
        d = decide(facts(event="schedule"), FakeApi())
        self.assertFalse(d.authorized)
        self.assertIn("unsupported", d.reason)


class Slugify(unittest.TestCase):
    def test_slash_and_metachars_collapse(self):
        self.assertEqual(slugify("feature/a b"), "feature-a-b")
        self.assertEqual(slugify("ok.name_1-2"), "ok.name_1-2")  # allowed chars kept


if __name__ == "__main__":
    unittest.main()
