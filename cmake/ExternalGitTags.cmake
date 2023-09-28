# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Repository URLs and tags for external third-party dependencies
#

if(__extern_git_tags)
  return()
endif()
set(__extern_git_tags YES)

# ARPACK-NG
set(EXTERN_ARPACK_URL
  "https://github.com/opencollab/arpack-ng.git" CACHE STRING
  "URL for external ARPACK-NG build"
)
set(EXTERN_ARPACK_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external ARPACK-NG build"
)
set(EXTERN_ARPACK_GIT_TAG
  "e7cf106bc204ed7f7909d906fd9b1c3d18e86f76" CACHE STRING  # 08/14/2023
  "Git tag for external ARPACK-NG build"
)

# ButterflyPACK (for STRUMPACK)
set(EXTERN_BUTTERFLYPACK_URL
  "https://github.com/liuyangzhuan/ButterflyPACK.git" CACHE STRING
  "URL for external ButterflyPACK build"
)
set(EXTERN_BUTTERFLYPACK_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external ButterflyPACK build"
)
set(EXTERN_BUTTERFLYPACK_GIT_TAG
  "6b0afaae91be4cd5e5cbc85f51f9c959921c0059" CACHE STRING  # 07/19/2023
  "Git tag for external ButterflyPACK build"
)

# GKlib (for METIS and ParMETIS)
set(EXTERN_GKLIB_URL
  "https://github.com/KarypisLab/GKlib.git" CACHE STRING
  "URL for external GKlib build"
)
set(EXTERN_GKLIB_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external GKlib build"
)
set(EXTERN_GKLIB_GIT_TAG
  "8bd6bad750b2b0d90800c632cf18e8ee93ad72d7" CACHE STRING  # 03/26/2023
  "Git tag for external GKlib build"
)

# GSLIB
set(EXTERN_GSLIB_URL
  "https://github.com/Nek5000/gslib.git" CACHE STRING
  "URL for external GSLIB build"
)
set(EXTERN_GSLIB_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external GSLIB build"
)
set(EXTERN_GSLIB_GIT_TAG
  "39d1baae8f4bfebe3ebca6a234dcc8ba1ee5edc7" CACHE STRING  # 11/09/2022
  "Git tag for external GSLIB build"
)

# HYPRE (for MFEM)
set(EXTERN_HYPRE_URL
  "https://github.com/hypre-space/hypre.git" CACHE STRING
  "URL for external HYPRE build"
)
set(EXTERN_HYPRE_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external HYPRE build"
)
set(EXTERN_HYPRE_GIT_TAG
  "b372b31a116a25d1432c689cef75881aac145a8e" CACHE STRING  # 08/16/2023
  "Git tag for external HYPRE build"
)

# libCEED
set(EXTERN_LIBCEED_URL
  "https://github.com/CEED/libCEED.git" CACHE STRING
  "URL for external libCEED build"
)
set(EXTERN_LIBCEED_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external libCEED build"
)
set(EXTERN_LIBCEED_GIT_TAG
  "2bc6469880c9a85ea757b18808a9311af8c3b635" CACHE STRING  # main @ 10/06/2023
  "Git tag for external libCEED build"
)

# LIBXSMM (for libCEED)
set(EXTERN_LIBXSMM_URL
  "https://github.com/hfp/libxsmm.git" CACHE STRING
  "URL for external LIBXSMM build"
)
set(EXTERN_LIBXSMM_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external LIBXSMM build"
)
set(EXTERN_LIBXSMM_GIT_TAG
  "f4ba09230d84b78db0187320677d3c71d0d4854f" CACHE STRING  # 08/21/2023
  "Git tag for external LIBXSMM build"
)

# METIS
set(EXTERN_METIS_URL
  "https://github.com/KarypisLab/METIS.git" CACHE STRING
  "URL for external METIS build"
)
set(EXTERN_METIS_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external METIS build"
)
set(EXTERN_METIS_GIT_TAG
  "e0f1b88b8efcb24ffa0ec55eabb78fbe61e58ae7" CACHE STRING  # 04/02/2023
  "Git tag for external METIS build"
)

# MFEM
set(EXTERN_MFEM_URL
  "https://github.com/mfem/mfem.git" CACHE STRING
  "URL for external MFEM build"
)
set(EXTERN_MFEM_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external MFEM build"
)
set(EXTERN_MFEM_GIT_TAG
  "c0867ca009e9b5fafb5585318648f430b1d6eef0" CACHE STRING  # master @ 10/02/2023
  "Git tag for external MFEM build"
)

# MUMPS
set(EXTERN_MUMPS_URL
  "https://github.com/scivision/mumps.git" CACHE STRING
  "URL for external MUMPS build"
)
set(EXTERN_MUMPS_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external MUMPS build"
)
set(EXTERN_MUMPS_GIT_TAG
  "e7d12e84c35706433331c058c4bf00186d6b852d" CACHE STRING  # 08/23/2023
  "Git tag for external MUMPS build"
)

# ParMETIS
set(EXTERN_PARMETIS_URL
  "https://github.com/KarypisLab/ParMETIS.git" CACHE STRING
  "URL for external ParMETIS build"
)
set(EXTERN_PARMETIS_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external ParMETIS build"
)
set(EXTERN_PARMETIS_GIT_TAG
  "8ee6a372ca703836f593e3c450ca903f04be14df" CACHE STRING  # 03/26/2023
  "Git tag for external ParMETIS build"
)

# PETSc (for SLEPc)
set(EXTERN_PETSC_URL
  "https://gitlab.com/petsc/petsc.git" CACHE STRING
  "URL for external PETSc build"
)
set(EXTERN_PETSC_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external PETSc build"
)
set(EXTERN_PETSC_GIT_TAG
  "01a4aebd931b780f244a87406598026c571af919" CACHE STRING  # 08/27/2023
  "Git tag for external PETSc build"
)

# ScaLAPACK (for STRUMPACK and MUMPS)
set(EXTERN_SCALAPACK_URL
  "https://github.com/scivision/scalapack.git" CACHE STRING
  "URL for external ScaLAPACK build"
)
set(EXTERN_SCALAPACK_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external ScaLAPACK build"
)
set(EXTERN_SCALAPACK_GIT_TAG
  "c180ae76969a929a4ab7e95251cedd1f866edd48" CACHE STRING  # 07/24/2022
  "Git tag for external ScaLAPACK build"
)

# SLEPc
set(EXTERN_SLEPC_URL
  "https://gitlab.com/slepc/slepc.git" CACHE STRING
  "URL for external SLEPc build"
)
set(EXTERN_SLEPC_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external SLEPc build"
)
set(EXTERN_SLEPC_GIT_TAG
  "2f6668495ab0d84da740777486385283a890d952" CACHE STRING  # 08/27/2023
  "Git tag for external SLEPc build"
)

# STRUMPACK
set(EXTERN_STRUMPACK_URL
  "https://github.com/pghysels/STRUMPACK.git" CACHE STRING
  "URL for external STRUMPACK build"
)
set(EXTERN_STRUMPACK_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external STRUMPACK build"
)
set(EXTERN_STRUMPACK_GIT_TAG
  "59d267964850b3dc3a97b37f2dba9d1b8228100f" CACHE STRING  # 08/25/2023
  "Git tag for external STRUMPACK build"
)

# SuperLU_DIST
set(EXTERN_SUPERLU_URL
  "https://github.com/xiaoyeli/superlu_dist.git" CACHE STRING
  "URL for external SuperLU_DIST build"
)
set(EXTERN_SUPERLU_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external SuperLU_DIST build"
)
set(EXTERN_SUPERLU_GIT_TAG
  "7af577d8bfd6defb65a7e5c5d585214ced77959b" CACHE STRING  # 07/27/2023
  "Git tag for external SuperLU_DIST build"
)

# ZFP (for STRUMPACK)
set(EXTERN_ZFP_URL
  "https://github.com/LLNL/zfp.git" CACHE STRING
  "URL for external ZFP build"
)
set(EXTERN_ZFP_GIT_BRANCH
  "develop" CACHE STRING
  "Git branch for external ZFP build"
)
set(EXTERN_ZFP_GIT_TAG
  "bc5e24145ea6f168cf5dd672e0f11544013004da" CACHE STRING  # 08/24/2023
  "Git tag for external ZFP build"
)

# nlohmann/json
set(EXTERN_JSON_URL
  "https://github.com/nlohmann/json/releases/download/v3.11.2/json.tar.xz" CACHE STRING
  "URL for external nlohmann/json build"
)

# fmt
set(EXTERN_FMT_URL
  "https://github.com/fmtlib/fmt/releases/download/10.0.0/fmt-10.0.0.zip" CACHE STRING
  "URL for external fmt build"
)

# Eigen
set(EXTERN_EIGEN_URL
  "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz" CACHE STRING
  "URL for external Eigen build"
)
