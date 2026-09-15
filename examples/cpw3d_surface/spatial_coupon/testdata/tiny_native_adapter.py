#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0
"""Test-only stand-in for the reviewed native MMG adapter command shape."""
from pathlib import Path
import shutil
import sys

if len(sys.argv) < 9:
    raise SystemExit(2)
output = Path(sys.argv[4])
if output.exists():
    raise SystemExit(2)
shutil.copyfile(sys.argv[1], output)
