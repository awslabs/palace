#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Publish a build run's container artifacts to ECR (OCI) and S3 (SIF).

Runs only after :mod:`authorize` has authorized the publish and produced the
selector. Splits cleanly into:

- :func:`plan` — pure logic: map each downloaded artifact leg to its ECR tag
  and S3 URI. Unit-tested, no I/O.
- :func:`main` — thin glue: ECR login, then run each planned copy via skopeo /
  aws with explicit argument lists (no shell, so no quoting pitfalls).

A native build produces one leg per arch; a release build one leg per
microarchitecture target. Each leg is a pair of downloaded artifacts:
``<image>-oci/<image>.tar`` and ``<image>-sif/<image>.sif`` where ``<image>`` is
``palace-<arch_label>-<shorthash>``.
"""

from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class PublishItem:
    image_name: str
    arch_label: str
    oci_tar: Path
    sif: Path
    ecr_tag: str  # full ECR ref, e.g. <registry>/palace:main-sapphirerapids
    s3_uri: str  # full s3:// destination for the SIF


IMAGE_NAME_RE = re.compile(r"^palace-(?P<label>[a-z0-9_]+)-(?P<hash>[0-9a-f]+)$")


def arch_label_from_image(image_name: str) -> str:
    """``palace-<arch_label>-<shorthash>`` -> ``<arch_label>``.

    Arch labels are ``[a-z0-9_]`` (no hyphens) and the short hash is hex, so the
    shape is unambiguous.
    """
    # We match that exact shape rather than split on the last hyphen: a loose
    # split would mis-read e.g. ``palace-foo-latest`` as label ``foo`` (hash
    # ``latest``) or absorb a hyphen into the label, and two names collapsing to
    # one label would then overwrite each other's tag/key.
    m = IMAGE_NAME_RE.match(image_name)
    if not m:
        raise ValueError(
            f"image name '{image_name}' is not palace-<label>-<hexhash> "
            f"(label [a-z0-9_], hash hex)"
        )
    return m.group("label")


def s3_prefix_for(selector: str) -> str:
    """S3 key prefix for a selector, preserving the dev/<branch> layout.

    ``dev-<branch>`` selectors publish under ``dev/<branch>/``; ``main`` and
    release selectors publish under ``<selector>/``.
    """
    if selector.startswith("dev-"):
        return f"dev/{selector[len('dev-'):]}"
    return selector


def plan(
    artifacts_dir: Path,
    selector: str,
    registry: str,
    ecr_repo: str,
    s3_bucket: str,
) -> list[PublishItem]:
    """Build the publish plan from the downloaded artifacts. Pure; no I/O beyond
    reading the artifact directory.

    Every leg must be complete (both the OCI tar and the SIF present) and its
    arch label unique, so an incomplete or mislabelled leg is a hard error
    rather than a silently partial publish. Completeness of the set of legs is
    guaranteed upstream: this workflow only runs when the whole `Containers`
    matrix succeeded (workflow_run.conclusion == 'success'), and the build's
    artifact uploads are unconditional, so a successful run has every leg.
    """
    items: list[PublishItem] = []
    seen_labels: set[str] = set()
    for oci_tar in sorted(artifacts_dir.glob("*-oci/*.tar")):
        image_name = oci_tar.parent.name[: -len("-oci")]  # strip "-oci" dir suffix
        arch_label = arch_label_from_image(image_name)
        if arch_label in seen_labels:
            raise ValueError(f"duplicate arch label '{arch_label}' among artifacts")
        seen_labels.add(arch_label)
        sif = artifacts_dir / f"{image_name}-sif" / f"{image_name}.sif"
        if not sif.is_file():
            raise ValueError(f"leg '{image_name}' has an OCI tar but no SIF at {sif}")
        ecr_tag = f"{registry}/{ecr_repo}:{selector}-{arch_label}"
        s3_uri = f"s3://{s3_bucket}/{s3_prefix_for(selector)}/{arch_label}.sif"
        items.append(PublishItem(image_name, arch_label, oci_tar, sif, ecr_tag, s3_uri))

    return items


def schema_s3_uri(selector: str, schema_bucket: str) -> str:
    """S3 destination for this build's config schema.

    The schema is arch-independent, so a single object per selector suffices; it
    mirrors the SIF prefix layout (``dev-<branch>`` -> ``dev/<branch>/``).
    """
    return f"s3://{schema_bucket}/{s3_prefix_for(selector)}/schema.json"


def find_schema(artifacts_dir: Path) -> Path:
    """Locate the bundled config schema among the downloaded artifacts.

    ``build-container`` uploads ``config-schema.json`` once per build leg as
    ``<image>-schema/config-schema.json``. The legs are byte-identical (the
    schema does not vary by architecture), so any one is authoritative — we take
    the first. Every build bundles it, so absence is a hard error (a
    silently-unpublished schema would leave downstream tools unable to validate
    against this build), mirroring the missing-SIF check in :func:`plan`.
    """
    matches = sorted(artifacts_dir.glob("*-schema/config-schema.json"))
    if not matches:
        raise ValueError("no config-schema.json artifact found among build artifacts")
    return matches[0]


def _run(cmd: list[str]) -> None:
    print("+ " + " ".join(cmd))
    subprocess.run(cmd, check=True)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--selector", required=True)
    parser.add_argument("--artifacts-dir", default="./artifacts")
    args = parser.parse_args(argv)

    region = os.environ["AWS_REGION"]
    ecr_repo = os.environ["ECR_REPOSITORY"]
    s3_bucket = os.environ["S3_CONTAINERS_BUCKET"]
    schema_bucket = os.environ["S3_SCHEMAS_BUCKET"]

    account_id = subprocess.run(
        ["aws", "sts", "get-caller-identity", "--query", "Account", "--output", "text"],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    print(f"::add-mask::{account_id}")
    registry = f"{account_id}.dkr.ecr.{region}.amazonaws.com"

    try:
        items = plan(Path(args.artifacts_dir), args.selector, registry, ecr_repo, s3_bucket)
        # Resolve the schema up front too, so a missing one fails before any
        # push rather than after the images are already live.
        schema = find_schema(Path(args.artifacts_dir))
    except ValueError as exc:
        print(f"::error::{exc}", file=sys.stderr)
        return 1
    if not items:
        print("::error::no OCI artifacts found to publish", file=sys.stderr)
        return 1

    # Log skopeo into ECR once (writes an auth file), so the registry password
    # is never placed on a command line / visible in the process list.
    ecr_pass = subprocess.run(
        ["aws", "ecr", "get-login-password", "--region", region],
        check=True, capture_output=True, text=True,
    ).stdout.strip()
    subprocess.run(
        ["skopeo", "login", "--username", "AWS", "--password-stdin", registry],
        input=ecr_pass, text=True, check=True,
    )

    for item in items:
        print(f"Publishing {item.image_name} -> {item.ecr_tag}")
        # skopeo is pre-installed on ubuntu-latest. The tar is
        # `podman save --format docker`, so use the docker-archive transport.
        _run(["skopeo", "copy", f"docker-archive:{item.oci_tar}", f"docker://{item.ecr_tag}"])
        _run(["aws", "s3", "cp", "--region", region, str(item.sif), item.s3_uri])

    print(f"Published {len(items)} image(s) for selector '{args.selector}'.")

    # Publish the config schema (resolved above) alongside the images so
    # downstream tools can fetch the exact contract this build was compiled
    # from. One arch-independent object per selector.
    schema_uri = schema_s3_uri(args.selector, schema_bucket)
    print(f"Publishing schema {schema} -> {schema_uri}")
    _run(["aws", "s3", "cp", "--region", region, str(schema), schema_uri])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
