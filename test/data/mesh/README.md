<!---
Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
SPDX-License-Identifier: Apache-2.0
--->

Most of these meshes come from the `data/` folder of the
[MFEM repository](https://github.com/mfem/mfem/tree/master/data). As of MFEM v4.10:

  - `fichera-mixed-p2.mesh`, `periodic-torus-sector.msh`, `star-mixed-p2.mesh`, and
    `tinyzoo-3d.mesh` are the MFEM meshes of the same name.
  - `fichera-hex.mesh` and `star-quad.mesh` are MFEM's `fichera.mesh` and `star.mesh`.
  - `fichera-tet.mesh` and `star-tri.mesh` are not in MFEM. They are `fichera.mesh` and
    `star.mesh` converted to simplices with `mfem::Mesh::MakeSimplicial` and written with
    `mfem::Mesh::Print`.

Some of these files were written with `mfem::Mesh::Print`, so their header comment differs
from the MFEM file. `nc-boundary-submesh.mesh` was written for *Palace*; its header
describes it.
