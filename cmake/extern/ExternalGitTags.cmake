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
  "db55a7ff4180bb440a980d618e302eded0db52f3" CACHE STRING  # 10/27/2022
  "Git tag for external ARPACK-NG build"
)
set(EXTERN_BUTTERFLYPACK_GIT_TAG
  "4f79c39f359454ccc091f36ebb803f68276f0645" CACHE STRING  # 10/21/2022
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
  "9bbdd9799fbb8def6f664eb52df74fa25b81380b" CACHE STRING  # 12/05/2022
  "Git tag for external HYPRE build"
)
set(EXTERN_METIS_GIT_TAG
  "f5ae915a84d3bbf1508b529af90292dd7085b9ec" CACHE STRING  # 12/05/2022
  "Git tag for external METIS build"
)
set(EXTERN_MUMPS_GIT_TAG
  "d293dfd680523d4d9b72c22df8fc5eddfd846703" CACHE STRING  # 12/01/2022
  "Git tag for external MUMPS build"
)
set(EXTERN_PARMETIS_GIT_TAG
  "44fadbf58c71a74b39abb110a7691355d2a760ca" CACHE STRING  # 01/09/2022
  "Git tag for external ParMETIS build"
)
set(EXTERN_PETSC_GIT_TAG
  "a44db9c6abf53727ac1441d1d480305cab36de07" CACHE STRING  # 12/06/2022
  "Git tag for external PETSc build"
)
set(EXTERN_SCALAPACK_GIT_TAG
  "bb40089771fcbd84140ea35796e70094df0d77d5" CACHE STRING  # 11/29/2022
  "Git tag for external ScaLAPACK build"
)
set(EXTERN_SLEPC_GIT_TAG
  "3ef5a0699975fbf8418b3a969196fe73958e7253" CACHE STRING  # 11/30/2022
  "Git tag for external SLEPc build"
)
set(EXTERN_STRUMPACK_GIT_TAG
  "e4b110b2d823c51a90575b77ec1531c699097a9f" CACHE STRING  # 12/01/2022
  "Git tag for external STRUMPACK build"
)
set(EXTERN_SUPERLU_GIT_TAG
  "72602c56e0fd332e150997f908c7edb314c5e915" CACHE STRING  # 12/05/2022
  "Git tag for external SuperLU_DIST build"
)
