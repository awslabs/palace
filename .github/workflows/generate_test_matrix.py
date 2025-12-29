#!/usr/bin/env python3

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
    ["~cuda"],
]

def is_valid(combo):
    if len(combo) < 4:
        # AllPairs creates partial combinations
        return True
    
    arch, compiler, mpi, math_libs = combo[0], combo[1], combo[2], combo[3]
    
    # Intel packages must all be together and only on x86
    intel_packages = [compiler == "intel-oneapi-compilers", 
                      mpi == "intel-oneapi-mpi", 
                      math_libs == "intel-oneapi-mkl"]
    if any(intel_packages):
        if arch != "x86" or not all(intel_packages):
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

# Add one case where we turn a multiple solvers (to check that there's no
# problem with compiling multiple solvers)
matrix.append({
    "arch": "x86",
    "compiler": "gcc",
    "mpi": "openmpi",
    "math-libs": "openblas",
    "shared": "~shared",
    "int": "~int64",
    "openmp": "+openmp",
    "eigensolver": "+slepc+arpack",
    "solver": "+superlu-dist+mumps+sundials+strumpack",
    "cuda": "~cuda"
})

print(yaml.dump(matrix, default_flow_style=False, sort_keys=False))
