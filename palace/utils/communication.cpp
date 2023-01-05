// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "communication.hpp"

namespace palace::mpi
{

template <>
const MPI_Datatype CommTrace<bool>::MPIType = MPI_C_BOOL;
template <>
const MPI_Datatype CommTrace<int>::MPIType = MPI_INTEGER;
template <>
const MPI_Datatype CommTrace<float>::MPIType = MPI_FLOAT;
template <>
const MPI_Datatype CommTrace<double>::MPIType = MPI_DOUBLE;
template <>
const MPI_Datatype CommTrace<std::complex<double>>::MPIType = MPI_DOUBLE;
template <>
const MPI_Datatype CommTrace<std::complex<float>>::MPIType = MPI_FLOAT;
#if defined(HYPRE_BIGINT) || defined(HYPRE_MIXEDINT)
template <>
const MPI_Datatype CommTrace<HYPRE_BigInt>::MPIType = HYPRE_MPI_BIG_INT;
#endif
template <>
const MPI_Datatype CommTrace<char>::MPIType = MPI_CHAR;

template <>
const int CommTrace<bool>::multiplicity = 1;
template <>
const int CommTrace<int>::multiplicity = 1;
template <>
const int CommTrace<float>::multiplicity = 1;
template <>
const int CommTrace<double>::multiplicity = 1;
template <>
const int CommTrace<std::complex<double>>::multiplicity = 2;
template <>
const int CommTrace<std::complex<float>>::multiplicity = 2;
#if defined(HYPRE_BIGINT) || defined(HYPRE_MIXEDINT)
template <>
const int CommTrace<HYPRE_BigInt>::multiplicity = 1;
#endif
template <>
const int CommTrace<char>::multiplicity = 1;

}  // namespace palace::mpi
