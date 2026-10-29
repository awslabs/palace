<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

These meshes were made for the *Palace* unit tests. The tests read the MFEM example meshes
they use from the `data/` folder of the MFEM installation with `TestMeshPath` in
`test/unit/test-helpers.hpp`, which also generates `star-tri.mesh` and `fichera-tet.mesh`
from MFEM's `star.mesh` and `fichera.mesh` with `mfem::Mesh::MakeSimplicial`.
