# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.mfem.package import Mfem as BuiltinMfem

from spack.package import *


class Mfem(BuiltinMfem):
    """MFEM, in this repository so that Palace can apply its patches to it."""
