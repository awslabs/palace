# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.slepc.package import Slepc as BuiltinSlepc

from spack.package import *


class Slepc(BuiltinSlepc):
    """SLEPc, in this repository so that Palace can apply its patches to it."""
