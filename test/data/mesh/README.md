<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

Most of these meshes are copies of files in the `data/` folder of MFEM v4.10
([d964264](https://github.com/mfem/mfem/tree/d964264cdb9a13e94a201b6c236c7721e0c8765f/data)):

  - `fichera-mixed-p2.mesh`, `periodic-torus-sector.msh`, `star-mixed-p2.mesh`, and
    `tinyzoo-3d.mesh` are the MFEM files of the same name.
  - `fichera-hex.mesh` and `star-quad.mesh` are MFEM's `fichera.mesh` and `star.mesh`.

The header comment of `fichera-hex.mesh`, `fichera-mixed-p2.mesh`, `star-mixed-p2.mesh`,
and `star-quad.mesh` differs from the MFEM file; the rest of each file is the same.
`nc-boundary-submesh.mesh` was written for *Palace*; its header describes it.

The unit tests also use `star-tri.mesh` and `fichera-tet.mesh`, which are not stored here.
`TestMeshPath` in `test/unit/test-helpers.hpp` generates them from `star-quad.mesh` and
`fichera-hex.mesh` with `mfem::Mesh::MakeSimplicial`.
