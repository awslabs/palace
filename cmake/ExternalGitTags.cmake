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
set(EXTERN_ARPACK_GIT_TAG
  "6d625dbf34375a0b5823d6c03bd16044f405baae" CACHE STRING  # 05/27/2023
  "Git tag for external ARPACK-NG build"
)

# ButterflyPACK (for STRUMPACK)
set(EXTERN_BUTTERFLYPACK_URL
  "https://github.com/liuyangzhuan/ButterflyPACK.git" CACHE STRING
  "URL for external ButterflyPACK build"
)
set(EXTERN_BUTTERFLYPACK_GIT_TAG
  "5a92f43697699e962f694536e07be291d4e9a3a7" CACHE STRING  # 06/14/2023
  "Git tag for external ButterflyPACK build"
)

# GKlib (for METIS and ParMETIS)
set(EXTERN_GKLIB_URL
  "https://github.com/KarypisLab/GKlib.git" CACHE STRING
  "URL for external GKlib build"
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
set(EXTERN_GSLIB_GIT_TAG
  "39d1baae8f4bfebe3ebca6a234dcc8ba1ee5edc7" CACHE STRING  # 11/09/2022
  "Git tag for external GSLIB build"
)

# HYPRE (for MFEM)
set(EXTERN_HYPRE_URL
  "https://github.com/hypre-space/hypre.git" CACHE STRING
  "URL for external HYPRE build"
)
set(EXTERN_HYPRE_GIT_TAG
  "72f5f3e136e5adec25f7cacc4710fab0f47ad610" CACHE STRING  # 06/15/2023
  "Git tag for external HYPRE build"
)

# libCEED
set(EXTERN_LIBCEED_URL
  "https://github.com/CEED/libCEED.git" CACHE STRING
  "URL for external libCEED build"
)
set(EXTERN_LIBCEED_GIT_TAG
  "bd882c8a454763a096666645dc9a6229d5263694" CACHE STRING  # main @ 06/15/2023
  "Git tag for external libCEED build"
)

# LIBXSMM (for libCEED)
set(EXTERN_LIBXSMM_URL
  "https://github.com/hfp/libxsmm.git" CACHE STRING
  "URL for external LIBXSMM build"
)
set(EXTERN_LIBXSMM_GIT_TAG
  "7ae1b5b14d5b0cc0d0fe3abe58ed0bbb1b876d8b" CACHE STRING  # 06/07/2023
  "Git tag for external LIBXSMM build"
)

# METIS
set(EXTERN_METIS_URL
  "https://github.com/KarypisLab/METIS.git" CACHE STRING
  "URL for external METIS build"
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
set(EXTERN_MFEM_GIT_TAG
  "6470d3a7b2edf868aace2b9454d95d124ff98173" CACHE STRING  # master @ 06/26/2023
  "Git tag for external MFEM build"
)

# MUMPS
set(EXTERN_MUMPS_URL
  "https://github.com/scivision/mumps.git" CACHE STRING
  "URL for external MUMPS build"
)
set(EXTERN_MUMPS_GIT_TAG
  "dc37cf6e3413f75cb39e867c4b7d0ce09d02a4cd" CACHE STRING  # 05/04/2023
  "Git tag for external MUMPS build"
)

# ParMETIS
set(EXTERN_PARMETIS_URL
  "https://github.com/KarypisLab/ParMETIS.git" CACHE STRING
  "URL for external ParMETIS build"
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
set(EXTERN_PETSC_GIT_TAG
  "ed63cbed510a5c8206399273090c321f5182507d" CACHE STRING  # 06/15/2023
  "Git tag for external PETSc build"
)

# ScaLAPACK (for STRUMPACK and MUMPS)
set(EXTERN_SCALAPACK_URL
  "https://github.com/scivision/scalapack.git" CACHE STRING
  "URL for external ScaLAPACK build"
)
set(EXTERN_SCALAPACK_GIT_TAG
  "b24a040ce5d9f7d262cef223134bd12d372cd72f" CACHE STRING  # 02/22/2022
  "Git tag for external ScaLAPACK build"
)

# SLEPc
set(EXTERN_SLEPC_URL
  "https://gitlab.com/slepc/slepc.git" CACHE STRING
  "URL for external SLEPc build"
)
set(EXTERN_SLEPC_GIT_TAG
  "cb46717a0b8405ec0bc9db973003912fb563af7f" CACHE STRING  # 06/13/2023
  "Git tag for external SLEPc build"
)

# STRUMPACK
set(EXTERN_STRUMPACK_URL
  "https://github.com/pghysels/STRUMPACK.git" CACHE STRING
  "URL for external STRUMPACK build"
)
set(EXTERN_STRUMPACK_GIT_TAG
  "17de5e2c833c0c51888036f83abfa4a49fa2bc8d" CACHE STRING  # 06/15/2023
  "Git tag for external STRUMPACK build"
)

# SuperLU_DIST
set(EXTERN_SUPERLU_URL
  "https://github.com/xiaoyeli/superlu_dist.git" CACHE STRING
  "URL for external SuperLU_DIST build"
)
set(EXTERN_SUPERLU_GIT_TAG
  "02b7c0d71bc33e785d098b0f8e4c26414bb8e39a" CACHE STRING  # 05/08/2023
  "Git tag for external SuperLU_DIST build"
)

# ZFP (for STRUMPACK)
set(EXTERN_ZFP_URL
  "https://github.com/LLNL/zfp.git" CACHE STRING
  "URL for external ZFP build"
)
set(EXTERN_ZFP_GIT_TAG
  "300e77d12f25d3eaa4ba9461d937fb17a71d45f6" CACHE STRING  # 04/29/2023
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
