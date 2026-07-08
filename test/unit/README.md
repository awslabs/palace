<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

# Palace unit tests

Build target: `palace-unit-tests`. Houses both fast unit tests and the
full-mesh regression cases (in `regression/cases.cpp`).

See [`docs/src/developer/testing.md`](../../docs/src/developer/testing.md)
for build instructions, how the `[Serial]` / `[Parallel]` / `[GPU]` /
`[Regression]` / `[Long]` tags interact, override flags, re-baselining,
and how to add a new regression case.
