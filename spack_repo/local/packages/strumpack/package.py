# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

from spack_repo.builtin.packages.strumpack.package import Strumpack as BuiltinStrumpack

from spack.package import *


class Strumpack(BuiltinStrumpack):
    """STRUMPACK, with the version used by the Palace superbuild."""

    version("8.1.0", tag="v8.1.0", commit="819ff1beaef0379e575fd77baaed8393948ca6d9")
