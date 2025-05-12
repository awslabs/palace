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
  "804fa3149a0f773064198a8e883bd021832157ca" CACHE STRING
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
  "cbbee2999ea8ec82a1f9d06320522db22d5d3ad2" CACHE STRING
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
  "95acf5b42301d6cb48fda88d662f1d784b863089" CACHE STRING
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
  "be52325a3ed8923fb93af348b1262ecfe44ab5d2" CACHE STRING
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
  "f8fffaefe2b0bf336e4f45347da405b61e8c02d8" CACHE STRING
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
  "6ea2c101ee6774a240ae6b4fcd8728d4a608e3cc" CACHE STRING
  "Git tag for external LIBXSMM build"
)

# MAGMA
set(EXTERN_MAGMA_URL
  "https://github.com/icl-utk-edu/magma.git" CACHE STRING
  "URL for external MAGMA build"
)
set(EXTERN_MAGMA_GIT_BRANCH
  "master" CACHE STRING
  "Git branch for external MAGMA build"
)
set(EXTERN_MAGMA_GIT_TAG
  "3b3e6058ed82ea2b05269275ede03586bc2a0229" CACHE STRING
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
  "08c3082720ff9114b8e3cbaa4484a26739cd7d2d" CACHE STRING
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
  "2b3c9ad1491e1bc907f5d9d5c26efadf130f0dca" CACHE STRING
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
  "1cfd19699702f9a64ff5d45827d6025ff5c3873a" CACHE STRING
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
  "53c9341b6c1ba876c97567cb52ddfc87c159dc36" CACHE STRING
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
  "022ab154a4ee4f77bad574e56fe9a469a926dbf2" CACHE STRING
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
  "0e8767285b7a201c7b1ff34d2c2bb009534145df" CACHE STRING
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
  "94cd514ee23a6958ceca39b49c8ef779eb18da6c" CACHE STRING
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
  "115b152be9a5d0d77846e3694f699c53c93fe394" CACHE STRING
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
  "461f3cf71c76da32d6bab8f44d1cf82e2a51a54a" CACHE STRING
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
  "7d1e3a21047a976599b562b3bbd53b1f34348f1a" CACHE STRING
  "Git tag for external ZFP build"
)

# nlohmann/json
set(EXTERN_JSON_URL
  "https://github.com/nlohmann/json/releases/download/v3.11.3/json.tar.xz" CACHE STRING
  "URL for external nlohmann/json build"
)

# fmt
set(EXTERN_FMT_URL
  "https://github.com/fmtlib/fmt/releases/download/10.2.1/fmt-10.2.1.zip" CACHE STRING
  "URL for external fmt build"
)

# Eigen
set(EXTERN_EIGEN_URL
  "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz" CACHE STRING
  "URL for external Eigen build"
)

# SUNDIALS
set(EXTERN_SUNDIALS_URL
  "https://github.com/LLNL/sundials.git" CACHE STRING
  "URL for external SUNDIALS build"
)
set(EXTERN_SUNDIALS_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external SUNDIALS build"
)
set(EXTERN_SUNDIALS_GIT_TAG
  "e941546af1a5b5e492dcac0a1872540e9961c556" CACHE STRING
  "Git tag for external SUNDIALS build"
)