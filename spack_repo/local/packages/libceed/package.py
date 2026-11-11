# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.libceed.package import Libceed as BuiltinLibceed

from spack.package import *


class Libceed(BuiltinLibceed):
    """libCEED, in this repository so that Palace can apply its patches to it."""

    # The version of the superbuild (cmake/ExternalGitTags.cmake), which the patches of
    # Palace are made against.
    version("1.0.1", tag="v1.0.1", commit="cbcc9c934bce2efc324fbd2edf95032911109353")
