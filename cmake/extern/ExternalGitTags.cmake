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
  "06657372bb941d6e0bd099de5a12b076b3dfae27" CACHE STRING  # 01/19/2023
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
  "d3f6b03c87f12aabc94f8c2c606b6000d50786fc" CACHE STRING  # 01/23/2023
  "Git tag for external HYPRE build"
)
set(EXTERN_METIS_GIT_TAG
  "f5ae915a84d3bbf1508b529af90292dd7085b9ec" CACHE STRING  # 12/05/2022
  "Git tag for external METIS build"
)
set(EXTERN_MUMPS_GIT_TAG
  "38bf3b3f598d983a630cb7d62561c140c27e7019" CACHE STRING  # 01/26/2023
  "Git tag for external MUMPS build"
)
set(EXTERN_PARMETIS_GIT_TAG
  "44fadbf58c71a74b39abb110a7691355d2a760ca" CACHE STRING  # 01/09/2022
  "Git tag for external ParMETIS build"
)
set(EXTERN_PETSC_GIT_TAG
  "fc0b9b6be2232f473472b55ad745d8ca4e2668e2" CACHE STRING  # 01/26/2023
  "Git tag for external PETSc build"
)
set(EXTERN_SCALAPACK_GIT_TAG
  "bb40089771fcbd84140ea35796e70094df0d77d5" CACHE STRING  # 11/29/2022
  "Git tag for external ScaLAPACK build"
)
set(EXTERN_STRUMPACK_GIT_TAG
  "882f7611c1aab9f57cc2e08950530f0efeec6edf" CACHE STRING  # 01/24/2023
  "Git tag for external STRUMPACK build"
)
set(EXTERN_SUPERLU_GIT_TAG
  "1aa8e658586abed61bed519aeee2a20bec99e0d7" CACHE STRING  # 01/25/2023
  "Git tag for external SuperLU_DIST build"
)
