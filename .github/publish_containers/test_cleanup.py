#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Unit tests for the dev-container cleanup sweep.

Run: python3 -m unittest discover -s .github/publish_containers -p 'test_*.py'

The sweep's delete decision is pure, so these run offline with a FakeApi.
"""

from __future__ import annotations

import unittest

from authorize import DEV_BRANCH_RE, slugify
from cleanup import (
    Plan,
    assert_dev_only,
    ecr_tag_selector,
    live_selectors,
    orphan_ecr_tags,
    orphan_s3_keys,
    s3_key_selector,
    sweep,
)


class FakeApi:
    """Stand-in for GitHubApi. All lookups come from constructor data."""

    def __init__(self, *, branches=None, latest_pr=None):
        # branches: list of branch names
        self._branches = list(branches or [])
        # latest_pr: {branch: (pr_number, merged_bool)}; absent => no PR
        self._latest_pr = latest_pr or {}
        # Set to a list of branch lists to make successive branches() calls
        # return successive entries, simulating the repo changing mid-sweep. The
        # last entry repeats once exhausted.
        self.branches_sequence = None
        self._branches_calls = 0

    def branches(self):
        if self.branches_sequence is None:
            return list(self._branches)
        idx = min(self._branches_calls, len(self.branches_sequence) - 1)
        self._branches_calls += 1
        return list(self.branches_sequence[idx])

    def latest_pr_for_branch(self, branch):
        return self._latest_pr.get(branch)


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
    def test_branch_with_no_pr_is_live(self):
        # Dispatch-published branches have no PR at all.
        api = FakeApi(branches=["team/dispatch"])
        self.assertEqual(live_selectors(api), {"dev-team-dispatch"})

    def test_branch_with_open_latest_pr_is_live(self):
        api = FakeApi(branches=["team/feat"], latest_pr={"team/feat": (100, False)})
        self.assertEqual(live_selectors(api), {"dev-team-feat"})

    def test_branch_with_merged_latest_pr_is_not_live(self):
        api = FakeApi(branches=["team/feat"], latest_pr={"team/feat": (100, True)})
        self.assertEqual(live_selectors(api), set())

    def test_deleted_branch_contributes_nothing(self):
        api = FakeApi(branches=[])
        self.assertEqual(live_selectors(api), set())

    def test_reopened_pr_after_merge_makes_branch_live_again(self):
        # PR #100 merged, then #101 opened from the same branch. Keying on the
        # LATEST PR (highest number) is what stops the sweep from flapping.
        api = FakeApi(branches=["team/feat"], latest_pr={"team/feat": (101, False)})
        self.assertEqual(live_selectors(api), {"dev-team-feat"})

    def test_unpublishable_branch_shapes_are_ignored(self):
        # Branches that DEV_BRANCH_RE rejects can never have published, so they
        # must not mark any selector live (which would protect an orphan).
        api = FakeApi(branches=["flat", "a/b/c", "feature-a/b", "main"])
        self.assertEqual(live_selectors(api), set())

    def test_main_never_marks_a_dev_selector_live(self):
        api = FakeApi(branches=["main"])
        self.assertNotIn("dev-main", live_selectors(api))


class Sweep(unittest.TestCase):
    """Listing ECR/S3 takes time, and a publish can land while it runs.

    The sweep samples the live set before AND after discovery and protects the
    UNION, so a channel that was live at either instant survives. Union, not
    intersection: erring toward keeping is recoverable (the next run cleans it),
    erring toward deleting destroys a channel someone just published.
    """

    def _discover(self, tags, keys):
        return lambda: (tags, {"buck": keys})

    def test_selector_published_during_discovery_is_protected(self):
        # branches() returns nothing first, then the branch appears: a developer
        # pushed and published while the sweep was listing.
        api = FakeApi(branches=[])
        api.branches_sequence = [[], ["team/new"]]
        plan = sweep(api, self._discover(["dev-team-new-x86_64_v3"], ["dev/team-new/x86_64_v3.sif"]))
        self.assertTrue(plan.is_empty(), f"republished channel must survive, got {plan}")

    def test_selector_deleted_during_discovery_is_protected_this_run(self):
        # Live at the first sample, gone by the second. Keeping it is safe; the
        # next scheduled run reclaims it.
        api = FakeApi(branches=[])
        api.branches_sequence = [["team/going"], []]
        plan = sweep(api, self._discover(["dev-team-going-x86_64_v3"], []))
        self.assertTrue(plan.is_empty())

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
