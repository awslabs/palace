#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests for the REAL GitHubApi methods (ref_sha / pr_for_head / pr_has_label).

test_authorize.py exercises decide() against a FakeApi, which necessarily
*pre-bakes* the sha+ref filtering these methods implement. So these tests cover
the actual logic — annotated-tag peeling, the head.sha+head.ref filter over the
"PRs containing this commit" list, and malformed-response handling — by mocking
only the subprocess boundary (GitHubApi._get) with canned API responses keyed on
the request path.

Run: python3 -m unittest discover -s .github/publish -p 'test_*.py'
"""

from __future__ import annotations

import unittest

from authorize import GitHubApi

REPO = "awslabs/palace"
COMMIT = "c45e4ab4f2356589ba40f4fde2403f5854282826"
TAG_OBJ = "b4523792f95f2e796762476255b6adefeb399914"


class StubApi(GitHubApi):
    """GitHubApi with the subprocess boundary (_get) replaced by canned data.

    ``responses`` maps an API path to the object _get would have returned (a
    dict/list on success, or None to simulate a 404 / non-JSON / non-zero exit).
    """

    def __init__(self, responses):
        super().__init__(REPO)
        self._responses = responses

    def _get(self, path):
        return self._responses.get(path)


class RefSha(unittest.TestCase):
    def test_lightweight_tag_or_branch_returns_commit(self):
        api = StubApi({
            f"repos/{REPO}/git/ref/heads/main": {"object": {"type": "commit", "sha": COMMIT}},
        })
        self.assertEqual(api.ref_sha("heads", "main"), COMMIT)

    def test_annotated_tag_is_peeled_to_commit(self):
        # Annotated tag: the ref object is the TAG object; must peel via git/tags.
        api = StubApi({
            f"repos/{REPO}/git/ref/tags/v0.18.0": {"object": {"type": "tag", "sha": TAG_OBJ}},
            f"repos/{REPO}/git/tags/{TAG_OBJ}": {"object": {"type": "commit", "sha": COMMIT}},
        })
        self.assertEqual(api.ref_sha("tags", "v0.18.0"), COMMIT)

    def test_missing_ref_returns_none(self):
        self.assertIsNone(StubApi({}).ref_sha("heads", "nope"))

    def test_malformed_ref_object_returns_none(self):
        api = StubApi({f"repos/{REPO}/git/ref/heads/x": {"unexpected": True}})
        self.assertIsNone(api.ref_sha("heads", "x"))

    def test_annotated_tag_peel_404_returns_none(self):
        api = StubApi({
            f"repos/{REPO}/git/ref/tags/v9.9.9": {"object": {"type": "tag", "sha": TAG_OBJ}},
            # git/tags/<obj> intentionally absent -> peel fails -> None
        })
        self.assertIsNone(api.ref_sha("tags", "v9.9.9"))


class PrForHead(unittest.TestCase):
    def _pulls(self, entries):
        return {f"repos/{REPO}/commits/{COMMIT}/pulls": entries}

    def test_matches_on_sha_and_ref(self):
        api = StubApi(self._pulls([
            {"number": 10, "state": "open",
             "head": {"sha": COMMIT, "ref": "team/real", "repo": {"full_name": REPO}}},
        ]))
        self.assertEqual(api.pr_for_head(COMMIT, "team/real"), (10, REPO, "open"))

    def test_ignores_pr_that_only_contains_commit_wrong_ref(self):
        # The API returns PRs merely CONTAINING the commit: same sha, different
        # ref must NOT match (this is the shared-commit / label-ride defense).
        api = StubApi(self._pulls([
            {"number": 11, "state": "open",
             "head": {"sha": COMMIT, "ref": "team/other", "repo": {"full_name": REPO}}},
        ]))
        self.assertIsNone(api.pr_for_head(COMMIT, "team/real"))

    def test_ignores_pr_whose_head_sha_differs(self):
        api = StubApi(self._pulls([
            {"number": 12, "state": "open",
             "head": {"sha": "differentsha", "ref": "team/real", "repo": {"full_name": REPO}}},
        ]))
        self.assertIsNone(api.pr_for_head(COMMIT, "team/real"))

    def test_picks_correct_pr_among_several_sharing_commit(self):
        api = StubApi(self._pulls([
            {"number": 11, "state": "open",
             "head": {"sha": COMMIT, "ref": "team/evil", "repo": {"full_name": REPO}}},
            {"number": 10, "state": "closed",
             "head": {"sha": COMMIT, "ref": "team/real", "repo": {"full_name": "fork/palace"}}},
        ]))
        self.assertEqual(api.pr_for_head(COMMIT, "team/real"), (10, "fork/palace", "closed"))

    def test_empty_or_missing_list_returns_none(self):
        self.assertIsNone(StubApi({}).pr_for_head(COMMIT, "team/real"))
        self.assertIsNone(StubApi(self._pulls([])).pr_for_head(COMMIT, "team/real"))

    def test_malformed_elements_do_not_crash(self):
        api = StubApi(self._pulls([
            "not-a-dict",
            {"head": {"sha": COMMIT, "ref": "team/real", "repo": {"full_name": REPO}}},  # no number
        ]))
        # First element is skipped (not a dict); second matches sha+ref but has
        # no number -> returns None rather than raising.
        self.assertIsNone(api.pr_for_head(COMMIT, "team/real"))


class PrHasLabel(unittest.TestCase):
    def _labels(self, num, entries):
        return {f"repos/{REPO}/issues/{num}/labels": entries}

    def test_present(self):
        api = StubApi(self._labels(10, [{"name": "push-containers"}, {"name": "x"}]))
        self.assertTrue(api.pr_has_label(10, "push-containers"))

    def test_absent(self):
        api = StubApi(self._labels(10, [{"name": "x"}]))
        self.assertFalse(api.pr_has_label(10, "push-containers"))

    def test_missing_list_is_false(self):
        self.assertFalse(StubApi({}).pr_has_label(10, "push-containers"))

    def test_malformed_element_does_not_crash(self):
        api = StubApi(self._labels(10, ["not-a-dict", {"name": "push-containers"}]))
        self.assertTrue(api.pr_has_label(10, "push-containers"))


if __name__ == "__main__":
    unittest.main()
