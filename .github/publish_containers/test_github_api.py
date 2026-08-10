#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Tests for the REAL GitHubApi methods (_get / ref_sha / pr_for_head /
pr_has_label).

test_authorize.py exercises decide() against a FakeApi, which necessarily
*pre-bakes* the sha+ref filtering these methods implement. So these tests cover
the actual logic — annotated-tag peeling, the head.sha+head.ref filter over the
"PRs containing this commit" list, and malformed-response handling — by mocking
only the subprocess boundary with canned API responses keyed on the request
path: StubApi replaces ``_get`` for the method tests, while ``gh_response``
replaces ``subprocess.run`` to cover ``_get`` itself.

Run: python3 -m unittest discover -s .github/publish_containers -p 'test_*.py'
"""

from __future__ import annotations

import subprocess
import unittest
from unittest import mock

from authorize import GitHubApi

REPO = "awslabs/palace"
COMMIT = "c45e4ab4f2356589ba40f4fde2403f5854282826"
TAG_OBJ = "b4523792f95f2e796762476255b6adefeb399914"


class StubApi(GitHubApi):
    """GitHubApi with the subprocess boundary (_get) replaced by canned data.

    ``responses`` maps an API path to the object _get would have returned (a
    dict/list on success, or None to simulate a 404).
    """

    def __init__(self, responses):
        super().__init__(REPO)
        self._responses = responses

    def _get(self, path):
        return self._responses.get(path)


def gh_response(status, body, *, returncode=None, stderr=""):
    """A `gh api --include` result: status line, headers, blank line, body.

    Headers use CRLF and the status line LF, matching real `gh` output. gh exits
    non-zero for any non-2xx, so returncode defaults to match ``status``.
    """
    if returncode is None:
        returncode = 0 if 200 <= status < 300 else 1
    stdout = (
        f"HTTP/2.0 {status} Reason\n"
        "Content-Type: application/json; charset=utf-8\r\n"
        "\r\n"
        f"{body}"
    )
    return subprocess.CompletedProcess(
        args=["gh"], returncode=returncode, stdout=stdout, stderr=stderr
    )


class Get(unittest.TestCase):
    """_get: only an explicit 404 is "absent"; every other failure raises."""

    def _get(self, completed):
        with mock.patch("subprocess.run", return_value=completed):
            return GitHubApi(REPO)._get("repos/x/y")

    def test_success_returns_decoded_body(self):
        self.assertEqual(self._get(gh_response(200, '{"ok": true}')), {"ok": True})

    def test_not_found_returns_none(self):
        self.assertIsNone(self._get(gh_response(404, '{"message": "Not Found"}')))

    def test_forbidden_raises(self):
        # A rate limit or a revoked token must fail the job, not read as absent:
        # otherwise a release publishes nothing and still exits green.
        for status in (401, 403, 429, 500, 502):
            with self.subTest(status=status):
                with self.assertRaises(RuntimeError):
                    self._get(gh_response(status, '{"message": "nope"}'))

    def test_transport_failure_raises(self):
        # gh could not reach the API, so it printed no status line at all.
        no_status = subprocess.CompletedProcess(
            args=["gh"], returncode=1, stdout="", stderr="dial tcp: i/o timeout"
        )
        with self.assertRaises(RuntimeError):
            self._get(no_status)

    def test_non_json_success_raises(self):
        with self.assertRaises(RuntimeError):
            self._get(gh_response(200, "<html>gateway</html>"))

    def test_success_status_but_nonzero_exit_raises(self):
        with self.assertRaises(RuntimeError):
            self._get(gh_response(200, '{"ok": true}', returncode=1))


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
