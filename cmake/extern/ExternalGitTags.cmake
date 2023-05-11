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
  "3e3e12fef2c516fe21658d2e67dc8cbe86ec60bb" CACHE STRING  # 05/25/2023
  "Git tag for external ARPACK-NG build"
)
set(EXTERN_BUTTERFLYPACK_GIT_TAG
  "703d67848cc1901f3c827519e1780a1441fa49f6" CACHE STRING  # 05/09/2023
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
  "2b2e9d2eee238a72ee15ae68b1f3d7e22b6d495e" CACHE STRING  # 05/19/2023
  "Git tag for external HYPRE build"
)
set(EXTERN_METIS_GIT_TAG
  "e0f1b88b8efcb24ffa0ec55eabb78fbe61e58ae7" CACHE STRING  # 04/03/2023
  "Git tag for external METIS build"
)
set(EXTERN_MUMPS_GIT_TAG
  "dc37cf6e3413f75cb39e867c4b7d0ce09d02a4cd" CACHE STRING  # 05/04/2023
  "Git tag for external MUMPS build"
)
set(EXTERN_PARMETIS_GIT_TAG
  "8ee6a372ca703836f593e3c450ca903f04be14df" CACHE STRING  # 03/26/2023
  "Git tag for external ParMETIS build"
)
set(EXTERN_PETSC_GIT_TAG
  "6694d23023a85ac50ab9ebd4bb517715c165a429" CACHE STRING  # 05/25/2023
  "Git tag for external PETSc build"
)
set(EXTERN_SCALAPACK_GIT_TAG
  "b24a040ce5d9f7d262cef223134bd12d372cd72f" CACHE STRING  # 02/22/2022
  "Git tag for external ScaLAPACK build"
)
set(EXTERN_SLEPC_GIT_TAG
  "cc8b4002b6b025d3631b960c2f0fef89f469d199" CACHE STRING  # 05/23/2023
  "Git tag for external SLEPc build"
)
set(EXTERN_STRUMPACK_GIT_TAG
  "ce9d0a3dc5a37d08b1bf96e3ce023e80088a2576" CACHE STRING  # 05/25/2023
  "Git tag for external STRUMPACK build"
)
set(EXTERN_SUPERLU_GIT_TAG
  "02b7c0d71bc33e785d098b0f8e4c26414bb8e39a" CACHE STRING  # 05/08/2023
  "Git tag for external SuperLU_DIST build"
)
set(EXTERN_ZFP_GIT_TAG
  "300e77d12f25d3eaa4ba9461d937fb17a71d45f6" CACHE STRING  # 04/29/2023
  "Git tag for external ZFP build"
)
