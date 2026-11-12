# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.superlu_dist.package import SuperluDist as BuiltinSuperluDist

from spack.package import *


class SuperluDist(BuiltinSuperluDist):
    """SuperLU_DIST, in this repository so that Palace can apply its patches to it, with the
    version of the superbuild."""

    version("9.3.1", tag="v9.3.1", commit="9f7c8196d8f692cbfc59d81b61a74ceb4b37789b")

    # Without ParMETIS, keep the METIS ordering, and its default
    # (https://github.com/xiaoyeli/superlu_dist/pull/258)
    with when("@9.3.1~parmetis"):
        depends_on("metis@5:+int64", when="+int64")
        depends_on("metis@5:~int64", when="~int64")
        patch(
            "https://github.com/xiaoyeli/superlu_dist/commit/"
            "a7dfeb58cccd9c8ea363e52bb5eead815f3c32f7.diff?full_index=1",
            sha256="534a69a4a9bd0b0464d8a46c367c8e92c78de69d78aaa94e542e411300f76f5c",
        )

    def cmake_args(self):
        args = super().cmake_args()
        if self.spec.satisfies("@9.3.1~parmetis"):
            args.append(self.define("TPL_ENABLE_METISLIB", True))
            args.append(self.define("TPL_METIS_LIBRARIES", self.spec["metis"].libs.joined(";")))
        return args
