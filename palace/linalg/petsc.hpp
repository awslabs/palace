// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_PETSC_HPP
#define PALACE_LINALG_PETSC_HPP

#if defined(PALACE_WITH_SLEPC)

#include <petscsys.h>
#include <mfem.hpp>

#if !defined(PETSC_USE_REAL_DOUBLE)
#error "PETSc should be compiled with double precision!"
#endif
#if defined(PETSC_HAVE_HYPRE)
#error \
    "PETSc should be built without Hypre to avoid conflicts with MFEM's Hypre dependency!"
#endif
#if defined(PETSC_USE_64BIT_INDICES) && !(defined(HYPRE_BIGINT) || defined(HYPRE_MIXEDINT))
#warning "Mismatch between big HYPRE (32bit) and PETSc (64bit) integer types!"
#endif
#if !defined(PETSC_USE_64BIT_INDICES) && (defined(HYPRE_BIGINT) || defined(HYPRE_MIXEDINT))
#warning "Mismatch between big HYPRE (64bit) and PETSc (32bit) integer types!"
#endif
#if (defined(PETSC_HAVE_CUDA) && !defined(MFEM_USE_CUDA)) || \
    (!defined(PETSC_HAVE_CUDA) && defined(MFEM_USE_CUDA))
#error "Mismatch between MFEM and PETSc CUDA support!"
#endif
#if (defined(PETSC_HAVE_HIP) && !defined(MFEM_USE_HIP)) || \
    (!defined(PETSC_HAVE_HIP) && defined(MFEM_USE_HIP))
#error "Mismatch between MFEM and PETSc HIP support!"
#endif

// Forward declarations of PETSc objects.
typedef struct _p_Vec *Vec;
typedef struct _p_Mat *Mat;

// Error handling similar to Petsc's PetscCallAbort but always aborts on the global
// PETSC_COMM_WORLD communicator.
#define PalacePetscCall(...) PetscCallAbort(PETSC_COMM_WORLD, __VA_ARGS__)

#endif

#endif  // PALACE_LINALG_PETSC_HPP
