# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.strumpack.package import Strumpack as BuiltinStrumpack

from spack.package import *


class Strumpack(BuiltinStrumpack):
    """STRUMPACK, in this repository so that Palace can apply its patches to it."""

    def cmake_args(self):
        args = super().cmake_args()
        # The PBLAS block size of the superbuild, which the patch of Palace for
        # https://github.com/pghysels/STRUMPACK/pull/169 makes configurable in STRUMPACK
        # 8.1.0 (STRUMPACK uses 32, or 512 with SLATE, and ignores the value without it)
        if self.spec.satisfies("@8.1.0 ~slate"):
            args.append(self.define("STRUMPACK_PBLAS_BLOCKSIZE", "128"))
        return args
