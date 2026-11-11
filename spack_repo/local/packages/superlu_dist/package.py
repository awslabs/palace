# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.superlu_dist.package import SuperluDist as BuiltinSuperluDist

from spack.package import *


class SuperluDist(BuiltinSuperluDist):
    """SuperLU_DIST, in this repository so that Palace can apply its patches to it, with the
    version of the superbuild, and with its METIS ordering also without ParMETIS."""

    version("9.3.1", tag="v9.3.1", commit="9f7c8196d8f692cbfc59d81b61a74ceb4b37789b")

    # Without ParMETIS, SuperLU_DIST has no METIS ordering and aborts with "Invalid ISPEC"
    # when it is requested. The superbuild patch extern/patch/superlu_dist/patch_metis.diff
    # (metis-9.3.1.patch) enables it with METIS alone; metis-9.2.1.patch is its version for
    # 9.2.1.
    with when("@9.2.1,9.3.1 ~parmetis"):
        depends_on("metis@5:+int64", when="+int64")
        depends_on("metis@5:~int64", when="~int64")
    patch("metis-9.2.1.patch", when="@9.2.1~parmetis")
    patch("metis-9.3.1.patch", when="@9.3.1~parmetis")
    conflicts(
        "~parmetis",
        when="@9.3.2:",
        msg="Add extern/patch/superlu_dist/patch_metis.diff for the METIS ordering",
    )

    def cmake_args(self):
        args = super().cmake_args()
        if self.spec.satisfies("@9.2.1,9.3.1 ~parmetis"):
            args.append(self.define("TPL_ENABLE_METISLIB", True))
            args.append(self.define("TPL_METIS_LIBRARIES", self.spec["metis"].libs.joined(";")))
        return args
