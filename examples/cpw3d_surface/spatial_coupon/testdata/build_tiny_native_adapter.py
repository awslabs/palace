#!/usr/bin/env python3
# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Build the test-only native adapter and its fixture libmmg3d into a directory."""
from pathlib import Path
import shutil
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def compiler():
    return shutil.which("cc")


def build(directory):
    """Return (adapter, library) paths; the adapter links @rpath/$ORIGIN libmmg3d."""
    directory = Path(directory).resolve()
    directory.mkdir(parents=True, exist_ok=True)
    cc = compiler()
    if cc is None:
        raise RuntimeError("A C compiler (cc) is required to build the native fixture adapter")
    darwin = sys.platform == "darwin"
    library = directory / ("libmmg3d.dylib" if darwin else "libmmg3d.so")
    adapter = directory / "tiny_native_adapter"
    library_flags = (["-Wl,-install_name,@rpath/libmmg3d.dylib"] if darwin else
                     ["-fPIC", "-Wl,-soname,libmmg3d.so"])
    subprocess.run([cc, "-shared", *library_flags, "-o", str(library),
                    str(HERE / "tiny_mmg3d_fixture_library.c")], check=True)
    subprocess.run([cc, "-o", str(adapter), str(HERE / "tiny_native_adapter.c"),
                    f"-L{directory}", "-lmmg3d", f"-Wl,-rpath,{directory}"], check=True)
    return adapter, library


if __name__ == "__main__":
    print(*build(sys.argv[1]))
