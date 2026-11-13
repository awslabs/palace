# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.libceed.package import Libceed as BuiltinLibceed

from spack.package import *


class Libceed(BuiltinLibceed):
    """libCEED, in this repository to build the linear-element-assembly branch of a fork,
    which adds CeedOperatorSetLinearElementAssembly."""

    version(
        "1.1.0.lea",
        git="https://github.com/dnpham23/libCEED.git",
        commit="2d55eaab7b9ce16969f7963d149e1329660d0bf0",
    )
