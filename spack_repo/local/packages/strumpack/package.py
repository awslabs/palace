# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.strumpack.package import Strumpack as BuiltinStrumpack

from spack.package import *


class Strumpack(BuiltinStrumpack):
    """STRUMPACK, in this repository so that Palace can apply its patches to it."""

    # Without SLATE by default, as in the superbuild and the containers: with SLATE,
    # STRUMPACK factors and solves its distributed dense fronts with SLATE instead of
    # ScaLAPACK with the PBLAS block size set below, and SLATE requires OpenMP
    variant("slate", default=False, description="Build with SLATE support")

    def cmake_args(self):
        args = super().cmake_args()
        # The PBLAS block size of the superbuild, for builds without SLATE, which keep
        # the block size of 512 that STRUMPACK uses with it. The block size can be set
        # since https://github.com/pghysels/STRUMPACK/pull/169, which Palace applies to
        # STRUMPACK 8.1.0 as a patch (without it, 8.1.0 ignores the value and uses 32)
        if self.spec.satisfies("@8.1.0: ~slate"):
            args.append(self.define("STRUMPACK_PBLAS_BLOCKSIZE", "128"))
        return args
