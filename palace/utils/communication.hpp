// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_COMMUNICATION_HPP
#define PALACE_COMMUNICATION_HPP

#include <complex>
#include <fmt/format.h>
#include <fmt/printf.h>
#include <fmt/ranges.h>
#include <mfem.hpp>

namespace palace
{

namespace mpi
{

template <typename T>
struct CommTrace
{
  static const MPI_Datatype MPIType;
  static const int multiplicity;
};

template <>
const MPI_Datatype CommTrace<bool>::MPIType;
template <>
const int CommTrace<bool>::multiplicity;
template <>
const MPI_Datatype CommTrace<int>::MPIType;
template <>
const int CommTrace<int>::multiplicity;
template <>
const MPI_Datatype CommTrace<float>::MPIType;
template <>
const int CommTrace<float>::multiplicity;
template <>
const MPI_Datatype CommTrace<double>::MPIType;
template <>
const int CommTrace<double>::multiplicity;
template <>
const MPI_Datatype CommTrace<std::complex<double>>::MPIType;
template <>
const int CommTrace<std::complex<double>>::multiplicity;
template <>
const MPI_Datatype CommTrace<std::complex<float>>::MPIType;
template <>
const int CommTrace<std::complex<float>>::multiplicity;
#if defined(HYPRE_BIGINT) || defined(HYPRE_MIXEDINT)
template <>
const MPI_Datatype CommTrace<HYPRE_BigInt>::MPIType;
template <>
const int CommTrace<HYPRE_BigInt>::multiplicity;
#endif
template <>
const MPI_Datatype CommTrace<char>::MPIType;
template <>
const int CommTrace<char>::multiplicity;

}  // namespace mpi

//
// A simple convenience class for easy access to some MPI functionality. This is similar to
// mfem::Mpi and ideally should inherit from it, but the constructor being private instead
// of protected doesn't allow for that.
//
class Mpi
{
public:
  Mpi(const Mpi &) = delete;
  Mpi(Mpi &&) = delete;
  Mpi &operator=(const Mpi &) = delete;
  Mpi &operator=(Mpi &&) = delete;

  // Singleton creation.
  static void Init() { Init(nullptr, nullptr); }
  static void Init(int &argc, char **&argv) { Init(&argc, &argv); }

  // Finalize MPI (if it has been initialized and not yet already finalized).
  static void Finalize()
  {
    if (IsInitialized() && !IsFinalized())
    {
      MPI_Finalize();
    }
  }

  // Return true if MPI has been initialized.
  static bool IsInitialized()
  {
    int is_init;
    int ierr = MPI_Initialized(&is_init);
    return ierr == MPI_SUCCESS && is_init;
  }

  // Return true if MPI has been finalized.
  static bool IsFinalized()
  {
    int is_finalized;
    int ierr = MPI_Finalized(&is_finalized);
    return ierr == MPI_SUCCESS && is_finalized;
  }

  // Call MPI_Abort with the given error code.
  static void Abort(int code, MPI_Comm comm = World()) { MPI_Abort(comm, code); }

  // Barrier on the communicator.
  static void Barrier(MPI_Comm comm = World()) { MPI_Barrier(comm); }

  // Return processor's rank in the communicator.
  static int Rank(MPI_Comm comm)
  {
    int rank;
    MPI_Comm_rank(comm, &rank);
    return rank;
  }

  // Return communicator size.
  static int Size(MPI_Comm comm)
  {
    int size;
    MPI_Comm_size(comm, &size);
    return size;
  }

  // Return communicator size.
  static bool Root(MPI_Comm comm) { return Rank(comm) == 0; }

  // Wrapper for MPI_AllReduce.
  template <typename T>
  static void GlobalOp(int len, T *buff, MPI_Op op, MPI_Comm comm)
  {
    MPI_Allreduce(MPI_IN_PLACE, buff, mpi::CommTrace<T>::multiplicity * len,
                  mpi::CommTrace<T>::MPIType, op, comm);
  }

  // Global minimum (in-place, result is broadcast to all processes).
  template <typename T>
  static void GlobalMin(int len, T *buff, MPI_Comm comm)
  {
    GlobalOp(len, buff, MPI_MIN, comm);
  }

  // Global maximum (in-place, result is broadcast to all processes).
  template <typename T>
  static void GlobalMax(int len, T *buff, MPI_Comm comm)
  {
    GlobalOp(len, buff, MPI_MAX, comm);
  }

  // Global sum (in-place, result is broadcast to all processes).
  template <typename T>
  static void GlobalSum(int len, T *buff, MPI_Comm comm)
  {
    GlobalOp(len, buff, MPI_SUM, comm);
  }

  // Global broadcast from root.
  template <typename T>
  static void Broadcast(int len, T *buff, int root, MPI_Comm comm)
  {
    MPI_Bcast(buff, mpi::CommTrace<T>::multiplicity * len, mpi::CommTrace<T>::MPIType, root,
              comm);
  }

  // Print methods only print on the root process of MPI_COMM_WORLD or a given MPI_Comm.
  template <typename... T>
  static void Print(MPI_Comm comm, fmt::format_string<T...> fmt, T &&...args)
  {
    if (Root(comm))
    {
      fmt::print(fmt, std::forward<T>(args)...);
    }
  }

  template <typename... T>
  static void Print(fmt::format_string<T...> fmt, T &&...args)
  {
    Print(World(), fmt, std::forward<T>(args)...);
  }

  template <typename... T>
  static void Printf(MPI_Comm comm, const char *format, T &&...args)
  {
    if (Root(comm))
    {
      fmt::printf(format, std::forward<T>(args)...);
    }
  }

  template <typename... T>
  static void Printf(const char *format, T &&...args)
  {
    Printf(World(), format, std::forward<T>(args)...);
  }

  template <typename... T>
  static void Warning(MPI_Comm comm, fmt::format_string<T...> fmt, T &&...args)
  {
    Print(comm, fmt::format("\nWarning:\n{}\n", fmt), std::forward<T>(args)...);
  }

  template <typename... T>
  static void Warning(fmt::format_string<T...> fmt, T &&...args)
  {
    Warning(World(), fmt, std::forward<T>(args)...);
  }

  // Return the global communicator.
  static MPI_Comm World() { return MPI_COMM_WORLD; }

private:
  Mpi() = default;
  ~Mpi() { Finalize(); }

  static void Init(int *argc, char ***argv)
  {
    // The Mpi object below needs to be created after MPI_Init() for some MPI
    // implementations.
    MFEM_VERIFY(!IsInitialized(), "MPI should not be initialized more than once!");
    MPI_Init(argc, argv);
    static Mpi mpi;
  }
};

}  // namespace palace

#endif  // PALACE_COMMUNICATION_HPP
