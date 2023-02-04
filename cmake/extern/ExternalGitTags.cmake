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
  "f4d7302a49509ddaf9984653f0e4858948026d8b" CACHE STRING  # 01/19/2023
  "Git tag for external ARPACK-NG build"
)
set(EXTERN_BUTTERFLYPACK_GIT_TAG
  "1c552b273752bd724cf769310d08cefe56000623" CACHE STRING  # 02/02/2023
  "Git tag for external ButterflyPACK build"
)
set(EXTERN_GKLIB_GIT_TAG
  "08b72b2c41f0ac2a825438649ee7361bf0b488c3" CACHE STRING  # 11/27/2022
  "Git tag for external GKlib build"
)
set(EXTERN_GSLIB_GIT_TAG
  "39d1baae8f4bfebe3ebca6a234dcc8ba1ee5edc7" CACHE STRING  # 11/09/2022
  "Git tag for external GSLIB build"
)
set(EXTERN_HYPRE_GIT_TAG
  "832ad94de321b583424a958b4659a7bbe7204676" CACHE STRING  # 02/03/2023
  "Git tag for external HYPRE build"
)
set(EXTERN_METIS_GIT_TAG
  "f5ae915a84d3bbf1508b529af90292dd7085b9ec" CACHE STRING  # 12/05/2022
  "Git tag for external METIS build"
)
set(EXTERN_MUMPS_GIT_TAG
  "e478370809c46add2a8f148d83ab4881154c0809" CACHE STRING  # 01/26/2023
  "Git tag for external MUMPS build"
)
set(EXTERN_PARMETIS_GIT_TAG
  "44fadbf58c71a74b39abb110a7691355d2a760ca" CACHE STRING  # 01/09/2022
  "Git tag for external ParMETIS build"
)
set(EXTERN_PETSC_GIT_TAG
  "995ec06f924a86c4d28df68d1fdd6572768b0de1" CACHE STRING  # 02/03/2023
  "Git tag for external PETSc build"
)
set(EXTERN_SLEPC_GIT_TAG
  "a3f36a3399b1875a9b133f47dc6da13878af6437" CACHE STRING  # 02/03/2023
  "Git tag for external SLEPc build"
)
set(EXTERN_SCALAPACK_GIT_TAG
  "bb40089771fcbd84140ea35796e70094df0d77d5" CACHE STRING  # 11/29/2022
  "Git tag for external ScaLAPACK build"
)
set(EXTERN_STRUMPACK_GIT_TAG
  "ef8a48536c8b1d984dbb0afd6931a830fcc6a934" CACHE STRING  # 02/02/2023
  "Git tag for external STRUMPACK build"
)
set(EXTERN_SUPERLU_GIT_TAG
  "5517e3d3a2b09ec2c28ef70d61d2ffee9c34eae3" CACHE STRING  # 01/31/2023
  "Git tag for external SuperLU_DIST build"
)
