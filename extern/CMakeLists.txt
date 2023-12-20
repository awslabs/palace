# Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
# SPDX-License-Identifier: Apache-2.0

#
# Download and build third-party dependencies
#

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

# Add SLATE (for STRUMPACK with GPU support)
if(PALACE_WITH_STRUMPACK AND (PALACE_WITH_CUDA OR PALACE_WITH_HIP))
  message(STATUS "===================== Configuring SLATE dependency =====================")
  include(ExternalSLATE)
endif()

# Add MAGMA (for STRUMPACK with GPU support and for libCEED)
if(PALACE_WITH_MAGMA)
  message(STATUS "===================== Configuring MAGMA dependency =====================")
  include(ExternalMAGMA)
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

# Add LIBXSMM (for libCEED)
if(PALACE_WITH_LIBXSMM)
  message(STATUS "==================== Configuring LIBXSMM dependency ====================")
  include(ExternalLIBXSMM)
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
