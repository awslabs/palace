# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.hypre.package import Hypre as BuiltinHypre

from spack.package import *


class Hypre(BuiltinHypre):
    """hypre, in this repository so that Palace can apply its patches to it."""
