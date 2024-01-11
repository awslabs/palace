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
  "0eaf38bd51efd3046f0a7e61a08ddcd7ed1f1fcc" CACHE STRING  # 11/22/2023
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
  "e1ad6091e8dc2cb906ec426222f4acedb4eeeff2" CACHE STRING  # 01/10/2024
  "Git tag for external ButterflyPACK build"
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
  "57bfb26e268ddf003668c5d0b5938ae258922a83" CACHE STRING  # 01/08/2024
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
  "699fe8f87dbb33a262c5d9f777cbfdff72550eee" CACHE STRING  # main @ 01/11/2024
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
  "462cbc42f7c5423642a340a71f825731a013ecd6" CACHE STRING  # 01/11/2023
  "Git tag for external LIBXSMM build"
)

# MAGMA
set(EXTERN_MAGMA_URL
  "https://bitbucket.org/icl/magma.git" CACHE STRING
  "URL for external MAGMA build"
)
set(EXTERN_MAGMA_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external MAGMA build"
)
set(EXTERN_MAGMA_GIT_TAG
  "061a800467aab6d15f4353c4fc0770375863b43e" CACHE STRING  # 12/29/2023
  "Git tag for external MAGMA build"
)

# METIS
set(EXTERN_METIS_URL
  "https://bitbucket.org/petsc/pkg-metis.git" CACHE STRING
  "URL for external METIS build"
)
set(EXTERN_METIS_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external METIS build"
)
set(EXTERN_METIS_GIT_TAG
  "8b194fdf09661ac41b36fa16db0474d38f46f1ac" CACHE STRING  # 01/10/2023
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
  "2b27d4815f6e549ffb01065e9ad2d3549706ccec" CACHE STRING  # master @ 01/11/2024
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
  "89d12a54114636a040baf652a3b654530ef2590f" CACHE STRING  # 01/11/2023
  "Git tag for external MUMPS build"
)

# ParMETIS
set(EXTERN_PARMETIS_URL
  "https://bitbucket.org/petsc/pkg-parmetis.git" CACHE STRING
  "URL for external ParMETIS build"
)
set(EXTERN_PARMETIS_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external ParMETIS build"
)
set(EXTERN_PARMETIS_GIT_TAG
  "f5e3aab04fd5fe6e09fa02f885c1c29d349f9f8b" CACHE STRING  # 01/11/2023
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
  "e98d5aa1d9a8b761c58f921ec2d97e44eaf55cc3" CACHE STRING  # 01/12/2024
  "Git tag for external PETSc build"
)

# ScaLAPACK (for STRUMPACK and MUMPS)
set(EXTERN_SCALAPACK_URL
  "https://github.com/Reference-ScaLAPACK/scalapack.git" CACHE STRING
  "URL for external ScaLAPACK build"
)
set(EXTERN_SCALAPACK_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external ScaLAPACK build"
)
set(EXTERN_SCALAPACK_GIT_TAG
  "8435bf3bc435d9611d096abab46ebf6beb1d35ea" CACHE STRING  # 11/15/2023
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
  "9fccf47e48b0f46fc169e86349cb3ddd4f8dee5b" CACHE STRING  # 01/12/2024
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
  "90dcaaa3c7f06e69f94da43aa1043dd9e030cc2a" CACHE STRING  # 12/14/2023
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
  "b3eecd3eaac3a1332d0d2c5fc052d1af114df192" CACHE STRING  # 11/17/2023
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
  "f61aac5ca3445ac48bab94d0e3da62bb77bdd989" CACHE STRING  # 01/11/2024
  "Git tag for external ZFP build"
)

# nlohmann/json
set(EXTERN_JSON_URL
  "https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz" CACHE STRING
  "URL for external nlohmann/json build"
)

# fmt
set(EXTERN_FMT_URL
  "https://github.com/fmtlib/fmt/releases/download/10.1.1/fmt-10.1.1.zip" CACHE STRING
  "URL for external fmt build"
)

# Eigen
set(EXTERN_EIGEN_URL
  "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz" CACHE STRING
  "URL for external Eigen build"
)
