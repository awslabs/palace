#!/usr/bin/env bash

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

readonly BREW_PACKAGES=(
  gcc
  autoconf
  automake
  cmake
  diffutils
  libtool
  m4
  perl
  pkgconf
  python
)
readonly SPACK_EXTERNALS=("${BREW_PACKAGES[@]:1}")

if [[ "$(uname -s)" != "Darwin" ]]; then
  echo "Error: this setup script only supports macOS." >&2
  exit 1
fi
if ! command -v brew >/dev/null; then
  echo "Error: Homebrew is not installed; see https://brew.sh/." >&2
  exit 1
fi
if ! command -v spack >/dev/null; then
  echo "Error: Spack is not available in the current shell." >&2
  exit 1
fi
if [[ ! -f spack.yaml ]]; then
  echo "Error: run this script from a directory containing spack.yaml." >&2
  exit 1
fi

# This is idempotent and upgrades installed formulae when they are outdated.
brew install "${BREW_PACKAGES[@]}"

spack -e . compiler find "$(brew --prefix gcc)/bin"

# Restrict discovery to Homebrew prefixes so runner-provided tools cannot
# accidentally differ from a developer's setup.
external_find_args=()
for tool in "${SPACK_EXTERNALS[@]}"; do
  external_find_args+=(--path "$(brew --prefix "${tool}")/bin")
done
spack -e . external find --not-buildable \
  "${external_find_args[@]}" "${SPACK_EXTERNALS[@]}"

# `external find` exits successfully when it finds nothing. Check now instead
# of silently falling back to a source build (which fails for Python with GCC).
packages_config="$(spack -e . config get packages)"
for tool in "${SPACK_EXTERNALS[@]}"; do
  if ! grep -q "^  ${tool}:" <<< "${packages_config}"; then
    echo "Error: expected external '${tool}' was not detected." >&2
    exit 1
  fi
done
