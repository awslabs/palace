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
  "b7130fd234a23232620dfabc52fa1c7df6065842" CACHE STRING
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
  "v1.0.9" CACHE STRING
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
  "fca98afcd78a26b784dc24b0c3edb0e6140598af" CACHE STRING
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
  "95bd1e908b16e04a70015e3a9a7fddec5e9c3fc8" CACHE STRING
  "Git tag for external libCEED build"
)

# LIBXSMM (for libCEED)
set(EXTERN_LIBXSMM_URL
  "https://github.com/libxsmm/libxsmm.git" CACHE STRING
  "URL for external LIBXSMM build"
)
set(EXTERN_LIBXSMM_GIT_BRANCH
  "main" CACHE STRING
  "Git branch for external LIBXSMM build"
)
# libxsmm does not tag versions very often (last was Dec 2021)
set(EXTERN_LIBXSMM_GIT_TAG
  "ea0b20499a41377bab148257240adbbfe1b4a333" CACHE STRING
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
  "07b2b05635f0510ea4538f7ab68e50dcf0c0c815" CACHE STRING
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
  "v5.1.0-p13" CACHE STRING
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
  "bf1b25d82e232909e1d3aa1f416dbdae9931a550" CACHE STRING
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
  "v5.8.0.0" CACHE STRING
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
  "v4.0.3-p10" CACHE STRING
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
  "v3.24.3" CACHE STRING
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
  "6423f17933eb9a2522814b78ab3c0d6da25ee85a" CACHE STRING
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
  "v3.24.1" CACHE STRING
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
  "cfba574f60b792433ee87e7e2a874eeb96c74d1f" CACHE STRING
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
  "v9.2.1" CACHE STRING
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
# ZFP does not tag versions very often (last was Dec 2023)
set(EXTERN_ZFP_GIT_TAG
  "a888ab7606c58cc5ab5f34298fec25473e0d28ad" CACHE STRING
  "Git tag for external ZFP build"
)

# nlohmann/json
set(EXTERN_JSON_URL
  "https://github.com/nlohmann/json/releases/download/v3.12.0/json.tar.xz" CACHE STRING
  "URL for external nlohmann/json build"
)

# fmt
set(EXTERN_FMT_URL
  "https://github.com/fmtlib/fmt/releases/download/12.1.0/fmt-12.1.0.zip" CACHE STRING
  "URL for external fmt build"
)

# scn
set(EXTERN_SCN_URL
  "https://github.com/eliaskosunen/scnlib/archive/refs/tags/v4.0.1.zip" CACHE STRING
  "URL for external scn build"
)

# Eigen
set(EXTERN_EIGEN_URL
  "https://gitlab.com/libeigen/eigen/-/archive/5.0.0/eigen-5.0.0.tar.gz" CACHE STRING
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
  "v7.5.0" CACHE STRING
  "Git tag for external SUNDIALS build"
)

# Umpire (for HYPRE)
set(EXTERN_UMPIRE_URL
  "https://github.com/LLNL/Umpire" CACHE STRING
  "URL for external UMPIRE build"
)
set(EXTERN_UMPIRE_GIT_BRANCH
  "develop" CACHE STRING
  "Git branch for external UMPIRE build"
)
set(EXTERN_UMPIRE_GIT_TAG
  "ac034c775c7b4d4f4ff0f66774398373cf5fce73" CACHE STRING
  "Git tag for external UMPIRE build"
)
