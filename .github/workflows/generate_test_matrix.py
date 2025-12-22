#!/usr/bin/env python3

# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

"""Test Matrix Generator

Generates pairwise test combinations used in spack.yml using
[allpairspy](https://github.com/thombashi/allpairspy).

## Usage

```bash
pip install -U pyyaml allpairspy
python generate_test_matrix.py
```

It enforces certain constrains.

"""
from allpairspy import AllPairs
import yaml

parameters = [
    ["x86", "x86", "arm"], # Favor x86
    ["gcc", "llvm", "intel-oneapi-compilers"],
    ["openmpi", "mpich", "intel-oneapi-mpi"],
    ["openblas", "amdblis", "armpl-gcc", "intel-oneapi-mkl"],
    ["+shared", "~shared"],
    ["+int64", "~int64"],
    ["~openmp", "+openmp"],
    ["+arpack", "+slepc"],
    ["+mumps", "+superlu-dist", "+strumpack"],
    ["~cuda", "+cuda"],
]

def is_valid(combo):
    if len(combo) < 4:
        # AllPairs creates partial combinations
        return True
   
    arch, compiler, mpi, math_libs = combo[0], combo[1], combo[2], combo[3]

    # AllPairs creates partial combinations, so we don't always have 9 entries
    openmp = combo[6] if len(combo) > 6 else "~openmp"
    cuda = combo[9] if len(combo) > 9 else "~cuda"

    # Intel packages must all be together and only on x86

    # We also do not allow the CUDA+Intel combination because the "g" instances
    # are all AMD
    intel_packages = [compiler == "intel-oneapi-compilers", 
                      mpi == "intel-oneapi-mpi", 
                      math_libs == "intel-oneapi-mkl"]
    if any(intel_packages):
        if arch != "x86" or not all(intel_packages) or cuda == "+cuda":
            return False

    if cuda == "+cuda":
        # We only have GPUs on x86 machines
        if arch != "x86":
            return False
        # We don't want to mix GPU with OpenMP
        if openmp == "+openmp":
            return False

    # Use ARMPL only with GCC and ARM
    if math_libs == "armpl-gcc" and (compiler != "gcc" or arch != "arm"):
        return False

    # On ARM, use armpl-gcc or openblas
    if arch == "arm" and math_libs not in ("armpl-gcc", "openblas"):
        return False
    
    return True

matrix = []
for combo in AllPairs(parameters, filter_func=is_valid):
    matrix.append({
        "arch": combo[0],
        "compiler": combo[1],
        "mpi": combo[2],
        "math-libs": combo[3],
        "shared": combo[4],
        "int": combo[5],
        "openmp": combo[6],
        "eigensolver": combo[7],
        "solver": combo[8],
        "cuda": combo[9],
    })

# Add more CUDA test cases with different solver/eigensolver combinations
# (AllPairs is highly constrained with CUDA, so we manually add a couple of
# extra cases)
desired_cuda_cases = [
    {
        "arch": "x86", "compiler": "gcc", "mpi": "openmpi", "math-libs": "openblas",
        "shared": "+shared", "int": "+int64", "openmp": "~openmp", 
        "eigensolver": "+slepc", "solver": "+superlu-dist", "cuda": "+cuda"
    },
    {
        "arch": "x86", "compiler": "llvm", "mpi": "openmpi", "math-libs": "amdblis",
        "shared": "~shared", "int": "~int64", "openmp": "~openmp",
        "eigensolver": "+arpack", "solver": "+strumpack", "cuda": "+cuda"
    },

]

# Add one case where we turn a multiple solvers (to check that there's no
# problem with compiling multiple solvers)
for cuda in ("~cuda", "+cuda"):
    matrix.append({
        "arch": "x86",
        "compiler": "gcc",
        "mpi": "openmpi",
        "math-libs": "openblas",
        "shared": "+shared",
        "int": "~int64",
        "openmp": "~openmp",
        "eigensolver": "+slepc+arpack",
        "solver": "+superlu-dist+mumps+sundials+strumpack",
        "cuda": cuda
    })

for cuda_case in desired_cuda_cases:
    exists = any(
        all(existing.get(k) == v for k, v in cuda_case.items())
        for existing in matrix
    )
    if not exists:
        matrix.append(cuda_case)

        
# Custom list class for flow-style output, so that we can print it as [self-hosted, gpu]
# instead of
# - self-hosted
# - gpu
class FlowList(list):
    pass

def represent_flow_list(dumper, data):
    return dumper.represent_sequence('tag:yaml.org,2002:seq', data, flow_style=True)

yaml.add_representer(FlowList, represent_flow_list)

# Add runners
for entry in matrix:
    if entry["arch"] == "arm":
        entry["runner"] = FlowList(["self-hosted", "arm64"])
    elif entry["cuda"] == "+cuda":
        entry["runner"] = FlowList(["self-hosted", "gpu"])
    else:
        entry["runner"] = "palace_ubuntu-latest_16-core"

print(yaml.dump(matrix, default_flow_style=False, sort_keys=False))
