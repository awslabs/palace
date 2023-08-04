// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_COMMUNICATION_HPP
#define PALACE_UTILS_COMMUNICATION_HPP

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
inline MPI_Datatype DataType();

template <>
inline MPI_Datatype DataType<char>()
{
  return MPI_CHAR;
}

template <>
inline MPI_Datatype DataType<signed char>()
{
  return MPI_SIGNED_CHAR;
}

template <>
inline MPI_Datatype DataType<unsigned char>()
{
  return MPI_UNSIGNED_CHAR;
}

template <>
inline MPI_Datatype DataType<signed short>()
{
  return MPI_SHORT;
}

template <>
inline MPI_Datatype DataType<unsigned short>()
{
  return MPI_UNSIGNED_SHORT;
}

template <>
inline MPI_Datatype DataType<signed int>()
{
  return MPI_INT;
}

template <>
inline MPI_Datatype DataType<unsigned int>()
{
  return MPI_UNSIGNED;
}

template <>
inline MPI_Datatype DataType<signed long int>()
{
  return MPI_LONG;
}

template <>
inline MPI_Datatype DataType<unsigned long int>()
{
  return MPI_UNSIGNED_LONG;
}

template <>
inline MPI_Datatype DataType<signed long long int>()
{
  return MPI_LONG_LONG;
}

template <>
inline MPI_Datatype DataType<unsigned long long int>()
{
  return MPI_UNSIGNED_LONG_LONG;
}

template <>
inline MPI_Datatype DataType<float>()
{
  return MPI_FLOAT;
}

template <>
inline MPI_Datatype DataType<double>()
{
  return MPI_DOUBLE;
}

template <>
inline MPI_Datatype DataType<long double>()
{
  return MPI_LONG_DOUBLE;
}

template <>
inline MPI_Datatype DataType<std::complex<float>>()
{
  return MPI_C_COMPLEX;
}

template <>
inline MPI_Datatype DataType<std::complex<double>>()
{
  return MPI_C_DOUBLE_COMPLEX;
}

template <>
inline MPI_Datatype DataType<std::complex<long double>>()
{
  return MPI_C_LONG_DOUBLE_COMPLEX;
}

template <>
inline MPI_Datatype DataType<bool>()
{
  return MPI_C_BOOL;
}

// Helper for returning the rank with a MIN or MAX reduction.
template <typename T>
struct ValueAndRank
{
  // The reduced value from the MINLOC or MAXLOC operation.
  T value;
  // The processor rank that has the value/
  int rank;
  // Convenience decay operators to allow use as the underlying value.
  inline operator const T &() const { return value; }
  inline operator T &() { return value; }
};

template <>
inline MPI_Datatype DataType<ValueAndRank<double>>()
{
  return MPI_DOUBLE_INT;
}

template <>
inline MPI_Datatype DataType<ValueAndRank<int>>()
{
  return MPI_2INT;
}

}  // namespace mpi

//
// A simple convenience class for easy access to some MPI functionality. This is similar to
// mfem::Mpi and ideally should inherit from it, but the constructor being private instead
// of protected doesn't allow for that.
//
class Mpi
{
public:
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
    MPI_Allreduce(MPI_IN_PLACE, buff, len, mpi::DataType<T>(), op, comm);
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

  // Return the rank that has the minimum over the communicator
  template <typename T>
  static mpi::ValueAndRank<T> FindMin(const T val, MPI_Comm comm)
  {
    mpi::ValueAndRank<T> local{val, Mpi::Rank(comm)};
    GlobalOp(1, &local, MPI_MINLOC, comm);
    return local;
  }

  // Return the rank that has the maximum over the communicator
  template <typename T>
  static mpi::ValueAndRank<T> FindMax(const T &val, MPI_Comm comm)
  {
    mpi::ValueAndRank<T> local{val, Mpi::Rank(comm)};
    GlobalOp(1, &local, MPI_MAXLOC, comm);
    return local;
  }

  // Global broadcast from root.
  template <typename T>
  static void Broadcast(int len, T *buff, int root, MPI_Comm comm)
  {
    MPI_Bcast(buff, len, mpi::DataType<T>(), root, comm);
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
    Print(comm, "\nWarning!\n");
    Print(comm, fmt, std::forward<T>(args)...);
    Print(comm, "\n");
  }

  template <typename... T>
  static void Warning(fmt::format_string<T...> fmt, T &&...args)
  {
    Warning(World(), fmt, std::forward<T>(args)...);
  }

  // Return the global communicator.
  static MPI_Comm World() { return MPI_COMM_WORLD; }

  // Generate a custom MPI_Datatype for a given std::array<T,N>. The types are
  // deleted on Finalization of the singleton.
  template <typename T, std::size_t N>
  static MPI_Datatype GetArrayType()
  {
    // This only works if the array type is size equivalent to a C-array.
    static_assert(sizeof(std::array<T, N>) == N * sizeof(T));

    static MPI_Datatype type;
    static bool set = false;
    if (set)
    {
      return type;
    }
    else
    {
      int block_lengths[N];
      MPI_Aint displacements[N];
      MPI_Datatype types[N];

      for (int i = 0; i < N; ++i)
      {
        block_lengths[i] = 1;
        displacements[i] = i * sizeof(T);
        types[i] = mpi::DataType<T>();
      }

      MPI_Type_create_struct(N, block_lengths, displacements, types, &type);
      MPI_Type_commit(&type);

      Instance().custom_types.push_back(type);
      set = true;
    }
    return type;
  }

private:
  // Prevent direct construction of objects of this class.
  Mpi() = default;
  ~Mpi()
  {
    // Remove any custom types generated since initialization
    for (auto &dt : custom_types)
    {
      MPI_Type_free(&dt);
    }
    Finalize();
  }

  // Access the singleton instance.
  static Mpi &Instance()
  {
    static Mpi mpi;
    return mpi;
  }

  std::vector<MPI_Datatype> custom_types;

  static void Init(int *argc, char ***argv)
  {
    // The Mpi object below needs to be created after MPI_Init() for some MPI
    // implementations.
    MFEM_VERIFY(!IsInitialized(), "MPI should not be initialized more than once!");
#if defined(MFEM_USE_OPENMP)
    int provided, requested = MPI_THREAD_FUNNELED;  // MPI_THREAD_MULTIPLE
    MPI_Init_thread(argc, argv, requested, &provided);
    MFEM_VERIFY(provided >= requested,
                "MPI could not provide the requested level of thread support!");
#else
    MPI_Init(argc, argv);
#endif
    // Initialize the singleton Instance.
    Instance();
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_COMMUNICATION_HPP
