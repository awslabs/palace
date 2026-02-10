<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

# Spack Package Patches

This directory contains patches for upstream
[spack-packages](https://github.com/spack/spack-packages) required for Palace to
build correctly. These patches address upstream issues, apply urgent fixes
before PRs are merged, and ensure CI stability.

Generate patches using `git diff` with the `.diff` extension:

```bash
git diff > package-name-description.diff
```

Please remove outdated patches if you find any (e.g., when upstream PRs are
merged). To this end, we recommend including the PR number to the name of the
file, e.g., `pr2580.diff`.
