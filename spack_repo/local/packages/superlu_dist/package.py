# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.superlu_dist.package import SuperluDist as BuiltinSuperluDist

from spack.package import *


class SuperluDist(BuiltinSuperluDist):
    """SuperLU_DIST, with its METIS ordering also without ParMETIS."""

    # Without ParMETIS, SuperLU_DIST has no METIS ordering and aborts with "Invalid ISPEC"
    # when it is requested. The superbuild patch extern/patch/superlu_dist/patch_metis.diff
    # enables it with METIS alone; this is its version for 9.2.1.
    with when("@9.2.1~parmetis"):
        depends_on("metis@5:+int64", when="+int64")
        depends_on("metis@5:~int64", when="~int64")
        patch("metis-9.2.1.patch")
    conflicts(
        "~parmetis",
        when="@9.3:",
        msg="Add extern/patch/superlu_dist/patch_metis.diff for the METIS ordering",
    )

    def cmake_args(self):
        args = super().cmake_args()
        if self.spec.satisfies("@9.2.1~parmetis"):
            args.append(self.define("TPL_ENABLE_METISLIB", True))
            args.append(self.define("TPL_METIS_LIBRARIES", self.spec["metis"].libs.joined(";")))
        return args
