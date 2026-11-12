# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.strumpack.package import Strumpack as BuiltinStrumpack

from spack.package import *


class Strumpack(BuiltinStrumpack):
    """STRUMPACK, in this repository so that Palace can apply its patches to it."""

    def cmake_args(self):
        args = super().cmake_args()
        # The PBLAS block size of the superbuild, which the strumpack_pblas_blocksize.diff
        # patch of Palace for STRUMPACK 8.1.0 makes configurable (STRUMPACK uses 32, or 512
        # with SLATE, and ignores the value without the patch)
        if self.spec.satisfies("@8.1.0 ~slate"):
            args.append(self.define("STRUMPACK_PBLAS_BLOCKSIZE", "128"))
        return args
