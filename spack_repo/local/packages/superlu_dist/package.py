# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.superlu_dist.package import SuperluDist as BuiltinSuperluDist

from spack.package import *


class SuperluDist(BuiltinSuperluDist):
    """SuperLU_DIST, in this repository so that Palace can apply its patches to it, with the
    version of the superbuild."""

    version("9.3.1", tag="v9.3.1", commit="9f7c8196d8f692cbfc59d81b61a74ceb4b37789b")
