#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Test-only stand-in for the reviewed native MMG adapter command shape."""
from pathlib import Path
import shutil
import sys

arguments = list(sys.argv)
if "--required-tetrahedra" in arguments:
    position = arguments.index("--required-tetrahedra")
    if position + 1 >= len(arguments) or not Path(arguments[position + 1]).is_file():
        raise SystemExit(2)
    del arguments[position:position + 2]
else:
    raise SystemExit(2)
if len(arguments) < 9:
    raise SystemExit(2)
output = Path(arguments[4])
if output.exists():
    raise SystemExit(2)
shutil.copyfile(arguments[1], output)
