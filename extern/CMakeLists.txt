# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Download and build third-party dependencies
#

# Only build if requested, user can alternatively provide dependencies directly
if(NOT PALACE_BUILD_EXTERNAL_DEPS)
  return()
endif()

# Avoid DOWNLOAD_EXTRACT_TIMESTAMP warning
if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.24.0")
  cmake_policy(SET CMP0135 NEW)
endif()

# Git tags for third-party dependency builds
include(ExternalGitTags)

# Add METIS/ParMETIS libraries
message(STATUS "================ Configuring METIS/ParMETIS dependencies ===============")
include(ExternalMETIS)

# Add SuperLU_DIST
if(PALACE_WITH_SUPERLU)
  message(STATUS "================= Configuring SuperLU_DIST dependency ==================")
  include(ExternalSuperLU_DIST)
endif()

# Add ScaLAPACK
if(PALACE_WITH_STRUMPACK OR PALACE_WITH_MUMPS)
  message(STATUS "=================== Configuring ScaLAPACK dependency ===================")
  include(ExternalScaLAPACK)
endif()

# Add STRUMPACK
if(PALACE_WITH_STRUMPACK)
  message(STATUS "=================== Configuring STRUMPACK dependency ===================")
  include(ExternalSTRUMPACK)
endif()

# Add MUMPS
if(PALACE_WITH_MUMPS)
  message(STATUS "===================== Configuring MUMPS dependency =====================")
  include(ExternalMUMPS)
endif()

# Add PETSc and SLEPc
if(PALACE_WITH_SLEPC)
  message(STATUS "=============== Configuring PETSc and SLEPc dependencies ===============")
  include(ExternalSLEPc)
endif()

# Add ARPACK
if(PALACE_WITH_ARPACK)
  message(STATUS "==================== Configuring ARPACK dependency =====================")
  include(ExternalARPACK)
endif()

# Add GSLIB
if(PALACE_WITH_GSLIB)
  message(STATUS "===================== Configuring GSLIB dependency =====================")
  include(ExternalGSLIB)
endif()

# Add HYPRE
message(STATUS "===================== Configuring HYPRE dependency =====================")
include(ExternalHYPRE)

# Add JSON
message(STATUS "================= Configuring nlohmann/json dependency =================")
include(ExternalJSON)

# Add fmt
message(STATUS "====================== Configuring fmt dependency ======================")
include(ExternalFmt)

# Add Eigen
message(STATUS "===================== Configuring Eigen dependency =====================")
include(ExternalEigen)

#XX TODO: libCEED with libxsmm (optional)
