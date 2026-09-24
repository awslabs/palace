#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests for the dev-container cleanup sweep.

Run: python3 -m unittest discover -s .github/publish_containers -p 'test_*.py'

The sweep's delete decision is pure, so these run offline with a FakeApi.
"""

from __future__ import annotations

import json
import subprocess
import unittest
from unittest import mock

from authorize import DEV_BRANCH_RE, slugify
from cleanup import (
    Plan,
    assert_dev_only,
    ecr_tag_selector,
    execute,
    live_selectors,
    orphan_ecr_tags,
    orphan_s3_keys,
    s3_key_selector,
    sweep,
)


class FakeApi:
    """Stand-in for CleanupApi. All lookups come from constructor data."""

    def __init__(self, *, branches=None):
        # branches: list of branch names the repo currently has
        self._branches = list(branches or [])

    def branches(self):
        return list(self._branches)


class EcrTagParsing(unittest.TestCase):
    def test_arch_label_is_last_hyphen_component(self):
        # Arch labels are [a-z0-9_] (no hyphens), so the last hyphen splits
        # selector from arch even when the slug itself contains hyphens.
        self.assertEqual(ecr_tag_selector("dev-gbozzola-foo-x86_64_v3"), "dev-gbozzola-foo")

    def test_multi_underscore_arch_label_with_suffix(self):
        self.assertEqual(
            ecr_tag_selector("dev-team-my.feat-x86_64_v3_cuda"), "dev-team-my.feat"
        )

    def test_non_dev_tags_are_not_parsed(self):
        for tag in ("main-sapphirerapids", "0.16.1-x86_64_v3", "latest", "development-foo-arm64"):
            self.assertIsNone(ecr_tag_selector(tag), f"{tag!r} must not parse as a dev tag")


class S3KeyParsing(unittest.TestCase):
    def test_sif_and_schema_keys_map_to_selector(self):
        self.assertEqual(s3_key_selector("dev/gbozzola-foo/x86_64_v3.sif"), "dev-gbozzola-foo")
        self.assertEqual(s3_key_selector("dev/gbozzola-foo/schema.json"), "dev-gbozzola-foo")

    def test_keys_outside_dev_namespace_are_not_parsed(self):
        for key in (
            "main/x86_64_v3.sif",
            "0.16.1/x86_64_v3.sif",
            "development/foo/x.sif",
            "dev/x86_64_v3.sif",        # too shallow: no slug component
            "dev/foo/bar/x.sif",        # too deep
            "dev/../main/x86_64_v3.sif",  # traversal-shaped
            "dev/./x.sif",
            "dev/../x.sif",
        ):
            self.assertIsNone(s3_key_selector(key), f"{key!r} must not parse as a dev key")


class PrefixGuard(unittest.TestCase):
    """A non-dev artifact must never be emitted, for any `live` set."""

    def test_main_and_release_ecr_tags_never_emitted(self):
        tags = ["main-sapphirerapids", "0.16.1-x86_64_v3", "dev-gone-x86_64_v3"]
        for live in (set(), {"dev-gone"}, {"main", "0.16.1"}):
            got = orphan_ecr_tags(tags, live)
            self.assertNotIn("main-sapphirerapids", got)
            self.assertNotIn("0.16.1-x86_64_v3", got)

    def test_main_and_release_s3_keys_never_emitted(self):
        keys = [
            "main/x86_64_v3.sif",
            "0.16.1/x86_64_v3.sif",
            "dev/gone/x86_64_v3.sif",
        ]
        for live in (set(), {"dev-gone"}, {"main", "0.16.1"}):
            got = orphan_s3_keys(keys, live)
            self.assertNotIn("main/x86_64_v3.sif", got)
            self.assertNotIn("0.16.1/x86_64_v3.sif", got)

    def test_every_emitted_ecr_tag_starts_with_dev(self):
        tags = ["main-a", "0.16.1-b", "dev-x-c", "weird", "dev/slash-d"]
        for tag in orphan_ecr_tags(tags, set()):
            self.assertTrue(tag.startswith("dev-"), f"{tag!r} escaped the guard")

    def test_every_emitted_s3_key_starts_with_dev_slash(self):
        keys = ["main/a.sif", "dev/x/b.sif", "../dev/x/c.sif", "dev/../main/d.sif"]
        for key in orphan_s3_keys(keys, set()):
            self.assertTrue(key.startswith("dev/"), f"{key!r} escaped the guard")


class OrphanSelection(unittest.TestCase):
    def test_live_selector_retained_orphan_emitted(self):
        tags = ["dev-team-live-x86_64_v3", "dev-team-gone-x86_64_v3"]
        got = orphan_ecr_tags(tags, {"dev-team-live"})
        self.assertEqual(got, ["dev-team-gone-x86_64_v3"])

    def test_all_arches_of_an_orphan_are_emitted(self):
        tags = [
            "dev-team-gone-x86_64_v3",
            "dev-team-gone-x86_64_v3_cuda",
            "dev-team-gone-neoverse_v2",
        ]
        self.assertEqual(sorted(orphan_ecr_tags(tags, set())), sorted(tags))

    def test_both_sif_and_schema_of_an_orphan_are_emitted(self):
        keys = ["dev/team-gone/x86_64_v3.sif", "dev/team-gone/schema.json"]
        self.assertEqual(sorted(orphan_s3_keys(keys, set())), sorted(keys))


class LiveSelectors(unittest.TestCase):
    def test_existing_publishable_branch_is_live(self):
        # A branch is live purely because it exists and matches DEV_BRANCH_RE;
        # merging deletes the branch, so no PR-state check is needed.
        api = FakeApi(branches=["team/feat"])
        self.assertEqual(live_selectors(api), {"dev-team-feat"})

    def test_multiple_branches_map_to_their_selectors(self):
        api = FakeApi(branches=["team/feat", "user/other"])
        self.assertEqual(live_selectors(api), {"dev-team-feat", "dev-user-other"})

    def test_deleted_branch_contributes_nothing(self):
        api = FakeApi(branches=[])
        self.assertEqual(live_selectors(api), set())

    def test_unpublishable_branch_shapes_are_ignored(self):
        # Branches that DEV_BRANCH_RE rejects can never have published, so they
        # must not mark any selector live (which would protect an orphan).
        api = FakeApi(branches=["flat", "a/b/c", "feature-a/b", "main"])
        self.assertEqual(live_selectors(api), set())

    def test_main_never_marks_a_dev_selector_live(self):
        api = FakeApi(branches=["main"])
        self.assertNotIn("dev-main", live_selectors(api))


class Sweep(unittest.TestCase):
    """The sweep samples the live set once, before discovery, and protects it.

    A channel that lands while discovery runs comes from a branch that already
    existed when it was dispatched, so it is already in that sample. Nothing
    published from a branch created after the sample can have finished yet, so
    there is nothing of its to delete this run.
    """

    def _discover(self, tags, keys):
        return lambda: (tags, {"buck": keys})

    def test_stable_orphan_is_deleted(self):
        api = FakeApi(branches=["team/live"])
        plan = sweep(
            api,
            self._discover(
                ["dev-team-live-x86_64_v3", "dev-team-gone-x86_64_v3"],
                ["dev/team-live/x86_64_v3.sif", "dev/team-gone/schema.json"],
            ),
        )
        self.assertEqual(plan.ecr_tags, ["dev-team-gone-x86_64_v3"])
        self.assertEqual(plan.s3, {"buck": ["dev/team-gone/schema.json"]})

    def test_live_branch_artifacts_are_protected(self):
        api = FakeApi(branches=["team/live"])
        plan = sweep(
            api,
            self._discover(["dev-team-live-x86_64_v3"], ["dev/team-live/x86_64_v3.sif"]),
        )
        self.assertTrue(plan.is_empty(), f"live channel must survive, got {plan}")

    def test_sweep_never_emits_non_dev_artifacts(self):
        api = FakeApi(branches=[])
        plan = sweep(
            api,
            self._discover(["main-sapphirerapids", "0.16.1-x86_64_v3"], ["main/x86_64_v3.sif"]),
        )
        self.assertTrue(plan.is_empty())


class AssertDevOnly(unittest.TestCase):
    """Boundary guard: the last check before a delete is issued."""

    def test_dev_only_plan_is_accepted(self):
        plan = Plan(ecr_tags=["dev-a-x86_64_v3"], s3={"buck": ["dev/a/schema.json"]})
        assert_dev_only(plan)  # must not raise

    def test_non_dev_ecr_tag_raises(self):
        for tag in ("main-sapphirerapids", "0.16.1-x86_64_v3", "latest"):
            with self.assertRaises(ValueError, msg=tag):
                assert_dev_only(Plan(ecr_tags=[tag], s3={}))

    def test_non_dev_s3_key_raises(self):
        for key in ("main/x86_64_v3.sif", "0.16.1/x86_64_v3.sif", "dev/../main/x.sif"):
            with self.assertRaises(ValueError, msg=key):
                assert_dev_only(Plan(ecr_tags=[], s3={"buck": [key]}))

    def test_empty_plan_is_accepted(self):
        assert_dev_only(Plan())


class Execute(unittest.TestCase):
    """Both AWS delete calls exit 0 while reporting per-item failures in the
    body, so `execute` must inspect the body and raise on a real failure —
    otherwise a permission gap looks like a clean sweep and orphans persist.
    """

    def _fake_run(self, responses):
        """Return a `subprocess.run` stand-in that replies from `responses`.

        `responses` is a list of dicts; each call pops the next and returns it
        as the JSON stdout of a completed process.
        """
        queue = list(responses)

        def run(cmd, check=False, capture_output=False, text=False):
            body = queue.pop(0)
            return subprocess.CompletedProcess(cmd, 0, stdout=json.dumps(body), stderr="")

        return run

    def test_clean_delete_succeeds(self):
        plan = Plan(ecr_tags=["dev-a-x86_64_v3"], s3={"buck": ["dev/a/schema.json"]})
        with mock.patch(
            "cleanup.subprocess.run",
            side_effect=self._fake_run([{"failures": []}, {"Errors": []}]),
        ):
            execute(plan, "us-west-2", "palace")  # must not raise

    def test_ecr_failure_raises(self):
        plan = Plan(ecr_tags=["dev-a-x86_64_v3"], s3={})
        failure = {"failures": [{"imageId": {"imageTag": "dev-a-x86_64_v3"},
                                 "failureCode": "ServerException",
                                 "failureReason": "boom"}]}
        with mock.patch("cleanup.subprocess.run", side_effect=self._fake_run([failure])):
            with self.assertRaises(RuntimeError):
                execute(plan, "us-west-2", "palace")

    def test_ecr_image_not_found_is_success(self):
        # A deleted-then-gone tag reports ImageNotFound; reruns must stay idempotent.
        plan = Plan(ecr_tags=["dev-a-x86_64_v3"], s3={})
        failure = {"failures": [{"imageId": {"imageTag": "dev-a-x86_64_v3"},
                                 "failureCode": "ImageNotFound",
                                 "failureReason": "does not exist"}]}
        with mock.patch("cleanup.subprocess.run", side_effect=self._fake_run([failure])):
            execute(plan, "us-west-2", "palace")  # must not raise

    def test_s3_error_raises(self):
        plan = Plan(ecr_tags=[], s3={"buck": ["dev/a/schema.json"]})
        errors = {"Errors": [{"Key": "dev/a/schema.json", "Code": "AccessDenied",
                              "Message": "no s3:DeleteObject"}]}
        with mock.patch("cleanup.subprocess.run", side_effect=self._fake_run([errors])):
            with self.assertRaises(RuntimeError):
                execute(plan, "us-west-2", "palace")


class SlugifyInjectivity(unittest.TestCase):
    """The sweep reasons backwards from slug to branch; that needs injectivity.

    If DEV_BRANCH_RE is ever widened to allow hyphens in the prefix, two
    branches could share one selector and deleting one would delete the other's
    containers. This test pins the assumption.
    """

    def test_hyphenated_prefix_is_not_publishable(self):
        # Both of these slugify to "gbozzola-foo-bar"; only one is publishable.
        self.assertEqual(slugify("gbozzola/foo-bar"), "gbozzola-foo-bar")
        self.assertEqual(slugify("gbozzola-foo/bar"), "gbozzola-foo-bar")
        self.assertIsNotNone(DEV_BRANCH_RE.match("gbozzola/foo-bar"))
        self.assertIsNone(DEV_BRANCH_RE.match("gbozzola-foo/bar"))


if __name__ == "__main__":
    unittest.main()
