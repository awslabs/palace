#!/usr/bin/env bash

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail

tools=(
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

# This is idempotent and upgrades installed formulae when they are outdated.
brew install gcc "${tools[@]}"

spack -e . compiler find "$(brew --prefix gcc)/bin"

paths=()
for tool in "${tools[@]}"; do
  paths+=(--path "$(brew --prefix "${tool}")/bin")
done
spack -e . external find --not-buildable "${paths[@]}" "${tools[@]}"
