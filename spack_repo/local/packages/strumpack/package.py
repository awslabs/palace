# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.strumpack.package import Strumpack as BuiltinStrumpack

from spack.package import *


class Strumpack(BuiltinStrumpack):
    """STRUMPACK, in this repository so that Palace can apply its patches to it."""
