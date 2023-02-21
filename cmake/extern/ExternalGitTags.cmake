# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Repository tags for external third-party dependencies
#

if(__extern_git_tags)
  return()
endif()
set(__extern_git_tags YES)

set(EXTERN_ARPACK_GIT_TAG
  "22172ff21222c92bb1360c39e8aa62747610a1dd" CACHE STRING  # 04/10/2023
  "Git tag for external ARPACK-NG build"
)
set(EXTERN_BUTTERFLYPACK_GIT_TAG
  "04d498f5ef170fade5fd9f3caed841354b215a5a" CACHE STRING  # 04/13/2023
  "Git tag for external ButterflyPACK build"
)
set(EXTERN_GKLIB_GIT_TAG
  "8bd6bad750b2b0d90800c632cf18e8ee93ad72d7" CACHE STRING  # 03/26/2023
  "Git tag for external GKlib build"
)
set(EXTERN_GSLIB_GIT_TAG
  "39d1baae8f4bfebe3ebca6a234dcc8ba1ee5edc7" CACHE STRING  # 11/09/2022
  "Git tag for external GSLIB build"
)
set(EXTERN_HYPRE_GIT_TAG
  "e351324df841fdb323e1791162915a55f1ba642a" CACHE STRING  # 05/02/2023
  "Git tag for external HYPRE build"
)
set(EXTERN_METIS_GIT_TAG
  "e0f1b88b8efcb24ffa0ec55eabb78fbe61e58ae7" CACHE STRING  # 04/03/2022
  "Git tag for external METIS build"
)
set(EXTERN_MUMPS_GIT_TAG
  "69b17761914d7dc71de6d67e4594c4cf55bceb5f" CACHE STRING  # 05/02/2023
  "Git tag for external MUMPS build"
)
set(EXTERN_PARMETIS_GIT_TAG
  "8ee6a372ca703836f593e3c450ca903f04be14df" CACHE STRING  # 03/26/2023
  "Git tag for external ParMETIS build"
)
set(EXTERN_PETSC_GIT_TAG
  "c377bc1b91c268822f21b2f974c0bfcf48ec1990" CACHE STRING  # 05/02/2023
  "Git tag for external PETSc build"
)
set(EXTERN_SLEPC_GIT_TAG
  "59bebb3495c16e27671f3799be79dd1fb541dd93" CACHE STRING  # 04/28/2023
  "Git tag for external SLEPc build"
)
set(EXTERN_SCALAPACK_GIT_TAG
  "b24a040ce5d9f7d262cef223134bd12d372cd72f" CACHE STRING  # 02/22/2022
  "Git tag for external ScaLAPACK build"
)
set(EXTERN_STRUMPACK_GIT_TAG
  "72d2236946ebb2e97e79faf6b90eaa1d7edd4769" CACHE STRING  # 05/02/2023
  "Git tag for external STRUMPACK build"
)
set(EXTERN_SUPERLU_GIT_TAG
  "09cb1430f7be288fd4d75b8ed461aa0b7e68fefe" CACHE STRING  # 05/02/2023
  "Git tag for external SuperLU_DIST build"
)
set(EXTERN_ZFP_GIT_TAG
  "300e77d12f25d3eaa4ba9461d937fb17a71d45f6" CACHE STRING  # 04/29/2023
  "Git tag for external ZFP build"
)
