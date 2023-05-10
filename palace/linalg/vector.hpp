// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LINALG_VECTOR_HPP
#define PALACE_LINALG_VECTOR_HPP

#include <mpi.h>
#include <mfem.hpp>

namespace palace
{

class ComplexVector;
class ParOperator;

using Vector = mfem::Vector;

namespace linalg
{

//
// Basic functions for parallel vectors distributed across MPI processes.
//

// Returns the global vector size.
HYPRE_BigInt GlobalSize(MPI_Comm comm, const Vector &x);

// Sets all entries of the vector to random numbers sampled from the [-1, 1].
void SetRandom(MPI_Comm comm, Vector &x, int seed = 0);
void SetRandomSign(MPI_Comm comm, Vector &x, int seed = 0);

// Calculate the vector 2-norm.
double Norml2(MPI_Comm comm, const Vector &x);

// Calculate the vector infinity-norm.
double Normlinf(MPI_Comm comm, const Vector &x);

// Calculate the vector 1-norm.
double Norml1(MPI_Comm comm, const Vector &x);

// Normalize the vector, possibly with respect to an SPD matrix B.
double Normalize(MPI_Comm comm, Vector &x);
double Normalize(MPI_Comm comm, Vector &x, const ParOperator &B, Vector &Bx);
double Normalize(MPI_Comm comm, ComplexVector &x, const ParOperator &B, ComplexVector &Bx);

}  // namespace linalg

}  // namespace palace

#endif  // PALACE_LINALG_VECTOR_HPP
