<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

# Changelog Entries

Each pull request must include a changelog entry file in this directory.

## File naming

Name the file to match the suffix of your branch name:

| Branch | File |
|--------|------|
| `hughcars/fix-port-excitation` | `changes/fix-port-excitation.md` |
| `schema-validation` | `changes/schema-validation.md` |

## Format

Use `### Category` headers with indented bullet descriptions. Valid categories:
`New Features`, `Bug Fixes`, `Interface Changes`, `Build system`.

```
### Bug Fixes

  - Fixed port excitation parsing to correctly handle signed integer values.
```

Multiple categories in one file are fine for multi-concern PRs.

**Do not** include the `[PR NNN]` link — it is added automatically at release time.
