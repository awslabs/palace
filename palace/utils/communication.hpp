// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_COMMUNICATION_HPP
#define PALACE_UTILS_COMMUNICATION_HPP

#include <algorithm>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <vector>
#include <fmt/color.h>
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

template <typename T, typename U>
struct ValueAndLoc
{
  T val;
  U loc;
};

template <>
inline MPI_Datatype DataType<ValueAndLoc<float, signed int>>()
{
  return MPI_FLOAT_INT;
}

template <>
inline MPI_Datatype DataType<ValueAndLoc<double, signed int>>()
{
  return MPI_DOUBLE_INT;
}

template <>
inline MPI_Datatype DataType<ValueAndLoc<long double, signed int>>()
{
  return MPI_LONG_DOUBLE_INT;
}

template <>
inline MPI_Datatype DataType<ValueAndLoc<signed short, signed int>>()
{
  return MPI_SHORT_INT;
}

template <>
inline MPI_Datatype DataType<ValueAndLoc<signed int, signed int>>()
{
  return MPI_2INT;
}

template <>
inline MPI_Datatype DataType<ValueAndLoc<signed long int, signed int>>()
{
  return MPI_LONG_INT;
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
  static void Init(int requested = default_thread_required)
  {
    Init(nullptr, nullptr, requested);
  }
  static void Init(int &argc, char **&argv, int requested = default_thread_required)
  {
    Init(&argc, &argv, requested);
  }

  // Finalize MPI (if it has been initialized and not yet already finalized).
  static void Finalize()
  {
    if (IsInitialized() && !IsFinalized())
    {
      MPI_Finalize();
    }
  }

  // Record that the program is returning from main. Until then, the singleton's destructor
  // treats the process exit as abnormal: an exit() called from a library in the middle of
  // the run (gslib's die() on one rank, decision 346 (b)) would otherwise run ~Mpi and
  // enter MPI_Finalize, which blocks until every other rank finalizes — a silent hang
  // instead of a failure.
  static void MarkNormalShutdown() { normal_shutdown = true; }

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

  // Global minimum with index (in-place, result is broadcast to all processes).
  template <typename T, typename U>
  static void GlobalMinLoc(int len, T *val, U *loc, MPI_Comm comm)
  {
    std::vector<mpi::ValueAndLoc<T, U>> buffer(len);
    for (int i = 0; i < len; i++)
    {
      buffer[i].val = val[i];
      buffer[i].loc = loc[i];
    }
    GlobalOp(len, buffer.data(), MPI_MINLOC, comm);
    for (int i = 0; i < len; i++)
    {
      val[i] = buffer[i].val;
      loc[i] = buffer[i].loc;
    }
  }

  // Global maximum with index (in-place, result is broadcast to all processes).
  template <typename T, typename U>
  static void GlobalMaxLoc(int len, T *val, U *loc, MPI_Comm comm)
  {
    std::vector<mpi::ValueAndLoc<T, U>> buffer(len);
    for (int i = 0; i < len; i++)
    {
      buffer[i].val = val[i];
      buffer[i].loc = loc[i];
    }
    GlobalOp(len, buffer.data(), MPI_MAXLOC, comm);
    for (int i = 0; i < len; i++)
    {
      val[i] = buffer[i].val;
      loc[i] = buffer[i].loc;
    }
  }

  // Global logical or (in-place, result is broadcast to all processes).
  static void GlobalOr(int len, bool *buff, MPI_Comm comm)
  {
    GlobalOp(len, buff, MPI_LOR, comm);
  }

  // Global logical and (in-place, result is broadcast to all processes).
  static void GlobalAnd(int len, bool *buff, MPI_Comm comm)
  {
    GlobalOp(len, buff, MPI_LAND, comm);
  }

  // Global broadcast from root.
  template <typename T>
  static void Broadcast(int len, T *buff, int root, MPI_Comm comm)
  {
    MPI_Bcast(buff, len, mpi::DataType<T>(), root, comm);
  }

  // Global broadcast from root for buffers whose length may exceed INT_MAX (e.g. a
  // serialized multi-GB mesh). Chunks the transfer to respect the int count limit of the
  // classic MPI interface. The chunk length is exposed for testing the chunking logic
  // with small buffers.
  template <typename T>
  static void BroadcastLarge(std::int64_t len, T *buff, int root, MPI_Comm comm,
                             std::int64_t max_chunk = std::numeric_limits<int>::max())
  {
    MFEM_ASSERT(len >= 0 && max_chunk > 0, "Invalid BroadcastLarge lengths!");
    for (std::int64_t offset = 0; offset < len; offset += max_chunk)
    {
      const int n = static_cast<int>(std::min(max_chunk, len - offset));
      MPI_Bcast(buff + offset, n, mpi::DataType<T>(), root, comm);
    }
  }

  // OI-SCALE-GATHER (idfix DESIGN section 9 (ii); hardening-1 H16): the gathered perimeter
  // / facet / patch tables of a chip-scale run exceed the int count of the classic MPI
  // interface. The *Large forms below gather / scatter / reduce 64-bit-sized
  // variable-length arrays in rounds of at most LargeChunk() elements per rank and INT_MAX
  // / 2 per call, in rank order, so that every result is element-for-element the one-call
  // form's (the P == S identity the unit tests read at forced small chunks). The default
  // chunk is INT_MAX / 8; PALACE_LARGE_MESSAGE_CHUNK (a positive integer, read once)
  // overrides it for the exactness runs at forced small batches.
  static std::int64_t LargeChunk()
  {
    static const std::int64_t chunk = []()
    {
      const char *value = std::getenv("PALACE_LARGE_MESSAGE_CHUNK");
      if (value && *value)
      {
        const long long parsed = std::atoll(value);
        if (parsed > 0)
        {
          return static_cast<std::int64_t>(parsed);
        }
      }
      return static_cast<std::int64_t>(std::numeric_limits<int>::max() / 8);
    }();
    return chunk;
  }

  // The per-rank 64-bit counts and offsets of a variable-length gather (every rank).
  struct LargeGatherPlan
  {
    std::vector<std::int64_t> counts;
    std::vector<std::int64_t> offsets;
    std::int64_t total = 0;
  };
  static LargeGatherPlan PlanLargeGather(std::int64_t local_count, MPI_Comm comm)
  {
    LargeGatherPlan plan;
    plan.counts.resize(Size(comm));
    Allgather(1, &local_count, plan.counts.data(), comm);
    plan.offsets.resize(plan.counts.size());
    for (std::size_t rank = 0; rank < plan.counts.size(); rank++)
    {
      plan.offsets[rank] = plan.total;
      plan.total += plan.counts[rank];
    }
    return plan;
  }

  // The per-rank counts of one round: at most `max_chunk` of the remaining elements per
  // rank and INT_MAX / 2 in total.
  static bool NextLargeRound(const LargeGatherPlan &plan, std::vector<std::int64_t> &sent,
                             std::int64_t max_chunk, std::vector<int> &counts,
                             std::vector<int> &displs)
  {
    const std::int64_t per_rank = std::max<std::int64_t>(
        1,
        std::min(max_chunk, static_cast<std::int64_t>(std::numeric_limits<int>::max() / 2) /
                                static_cast<std::int64_t>(plan.counts.size())));
    counts.assign(plan.counts.size(), 0);
    displs.assign(plan.counts.size(), 0);
    int total = 0;
    bool any = false;
    for (std::size_t rank = 0; rank < plan.counts.size(); rank++)
    {
      const std::int64_t remaining = plan.counts[rank] - sent[rank];
      counts[rank] = static_cast<int>(std::min(per_rank, remaining));
      displs[rank] = total;
      total += counts[rank];
      any = any || remaining > 0;
    }
    return any;
  }

  // Allgather of a variable-length array per rank, 64-bit sized, in rounds; the result on
  // every rank is the concatenation in rank order. `plan` (optional) receives the per-rank
  // counts / offsets.
  template <typename T>
  static std::vector<T> AllgatherLarge(const std::vector<T> &local, MPI_Comm comm,
                                       LargeGatherPlan *plan_out = nullptr,
                                       std::int64_t max_chunk = LargeChunk())
  {
    const LargeGatherPlan plan =
        PlanLargeGather(static_cast<std::int64_t>(local.size()), comm);
    std::vector<T> result(static_cast<std::size_t>(plan.total));
    std::vector<std::int64_t> sent(plan.counts.size(), 0);
    std::vector<int> counts, displs;
    const int rank = Rank(comm);
    while (NextLargeRound(plan, sent, max_chunk, counts, displs))
    {
      std::vector<T> round(static_cast<std::size_t>(displs.back() + counts.back()));
      MPI_Allgatherv(local.data() + sent[rank], counts[rank], mpi::DataType<T>(),
                     round.data(), counts.data(), displs.data(), mpi::DataType<T>(), comm);
      for (std::size_t r = 0; r < counts.size(); r++)
      {
        std::copy_n(round.data() + displs[r], counts[r],
                    result.data() + plan.offsets[r] + sent[r]);
        sent[r] += counts[r];
      }
    }
    if (plan_out)
    {
      *plan_out = plan;
    }
    return result;
  }

  // Gather of a variable-length array per rank onto the root, 64-bit sized, in rounds; the
  // root's result is the concatenation in rank order, the other ranks' is empty.
  template <typename T>
  static std::vector<T> GatherLarge(const std::vector<T> &local, int root, MPI_Comm comm,
                                    LargeGatherPlan *plan_out = nullptr,
                                    std::int64_t max_chunk = LargeChunk())
  {
    const LargeGatherPlan plan =
        PlanLargeGather(static_cast<std::int64_t>(local.size()), comm);
    const int rank = Rank(comm);
    std::vector<T> result(rank == root ? static_cast<std::size_t>(plan.total) : 0);
    std::vector<std::int64_t> sent(plan.counts.size(), 0);
    std::vector<int> counts, displs;
    while (NextLargeRound(plan, sent, max_chunk, counts, displs))
    {
      std::vector<T> round(
          rank == root ? static_cast<std::size_t>(displs.back() + counts.back()) : 0);
      MPI_Gatherv(local.data() + sent[rank], counts[rank], mpi::DataType<T>(), round.data(),
                  counts.data(), displs.data(), mpi::DataType<T>(), root, comm);
      for (std::size_t r = 0; r < counts.size(); r++)
      {
        if (rank == root)
        {
          std::copy_n(round.data() + displs[r], counts[r],
                      result.data() + plan.offsets[r] + sent[r]);
        }
        sent[r] += counts[r];
      }
    }
    if (plan_out)
    {
      *plan_out = plan;
    }
    return result;
  }

  // Scatter of per-rank variable-length pieces from the root, 64-bit sized, in rounds: the
  // root's `data` holds the pieces in rank order with the 64-bit `counts` (every rank
  // passes its own expected count as counts[rank]; the counts are taken from the root).
  template <typename T>
  static std::vector<T> ScattervLarge(const std::vector<T> &data,
                                      std::vector<std::int64_t> counts, int root,
                                      MPI_Comm comm, std::int64_t max_chunk = LargeChunk())
  {
    const int rank = Rank(comm);
    counts.resize(Size(comm));
    Broadcast(static_cast<int>(counts.size()), counts.data(), root, comm);
    LargeGatherPlan plan;
    plan.counts = counts;
    plan.offsets.resize(counts.size());
    for (std::size_t r = 0; r < counts.size(); r++)
    {
      plan.offsets[r] = plan.total;
      plan.total += counts[r];
    }
    std::vector<T> result(static_cast<std::size_t>(plan.counts[rank]));
    std::vector<std::int64_t> sent(plan.counts.size(), 0);
    std::vector<int> round_counts, displs;
    while (NextLargeRound(plan, sent, max_chunk, round_counts, displs))
    {
      std::vector<T> round(
          rank == root ? static_cast<std::size_t>(displs.back() + round_counts.back()) : 0);
      if (rank == root)
      {
        for (std::size_t r = 0; r < round_counts.size(); r++)
        {
          std::copy_n(data.data() + plan.offsets[r] + sent[r], round_counts[r],
                      round.data() + displs[r]);
        }
      }
      MPI_Scatterv(round.data(), round_counts.data(), displs.data(), mpi::DataType<T>(),
                   result.data() + sent[rank], round_counts[rank], mpi::DataType<T>(), root,
                   comm);
      for (std::size_t r = 0; r < round_counts.size(); r++)
      {
        sent[r] += round_counts[r];
      }
    }
    return result;
  }

  // Element-wise MAX reduction of a 64-bit-sized array on every rank, in chunks.
  template <typename T>
  static void GlobalMaxLarge(std::vector<T> &values, MPI_Comm comm,
                             std::int64_t max_chunk = LargeChunk())
  {
    const auto len = static_cast<std::int64_t>(values.size());
    for (std::int64_t offset = 0; offset < len; offset += max_chunk)
    {
      const int n = static_cast<int>(std::min(max_chunk, len - offset));
      MPI_Allreduce(MPI_IN_PLACE, values.data() + offset, n, mpi::DataType<T>(), MPI_MAX,
                    comm);
    }
  }

  // Allgather a fixed count of items per rank; every rank ends up with the concatenated
  // result (sized sendcount * Size(comm)) in recvbuf.
  template <typename T>
  static void Allgather(int sendcount, const T *sendbuf, T *recvbuf, MPI_Comm comm)
  {
    MPI_Allgather(sendbuf, sendcount, mpi::DataType<T>(), recvbuf, sendcount,
                  mpi::DataType<T>(), comm);
  }

  // Allgather a variable number of items per rank; recvcounts/displs are per-rank arrays
  // sized Size(comm). Every rank ends up with the complete concatenated picture.
  template <typename T>
  static void Allgatherv(int sendcount, const T *sendbuf, T *recvbuf, const int *recvcounts,
                         const int *displs, MPI_Comm comm)
  {
    MPI_Allgatherv(sendbuf, sendcount, mpi::DataType<T>(), recvbuf, recvcounts, displs,
                   mpi::DataType<T>(), comm);
  }

  // All-to-all exchange of a fixed count of items with every rank.
  template <typename T>
  static void Alltoall(int sendcount, const T *sendbuf, T *recvbuf, MPI_Comm comm)
  {
    MPI_Alltoall(sendbuf, sendcount, mpi::DataType<T>(), recvbuf, sendcount,
                 mpi::DataType<T>(), comm);
  }

  // All-to-all exchange with variable per-rank counts and displacements.
  template <typename T>
  static void Alltoallv(const T *sendbuf, const int *sendcounts, const int *senddispls,
                        T *recvbuf, const int *recvcounts, const int *recvdispls,
                        MPI_Comm comm)
  {
    MPI_Alltoallv(sendbuf, sendcounts, senddispls, mpi::DataType<T>(), recvbuf, recvcounts,
                  recvdispls, mpi::DataType<T>(), comm);
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
    Print(comm, "\n{}\n", fmt::styled("--> Warning!", fmt::fg(fmt::color::yellow)));
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

  // Default level of threading used in MPI_Init_thread unless provided to Init.
#if defined(MFEM_USE_OPENMP)
  inline static int default_thread_required = MPI_THREAD_FUNNELED;
#else
  inline static int default_thread_required = MPI_THREAD_SINGLE;
#endif

private:
  inline static bool normal_shutdown = false;

  // Prevent direct construction of objects of this class.
  Mpi() = default;
  ~Mpi()
  {
    if (!normal_shutdown && IsInitialized() && !IsFinalized())
    {
      int rank = -1;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      std::fprintf(stderr,
                   "Abnormal exit on MPI rank %d before the end of the program: aborting "
                   "all ranks!\n",
                   rank);
      std::fflush(stderr);
      MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    Finalize();
  }

  // Access the singleton instance.
  static Mpi &Instance()
  {
    static Mpi mpi;
    return mpi;
  }

  static void Init(int *argc, char ***argv, int requested)
  {
    // The Mpi object below needs to be created after MPI_Init() for some MPI
    // implementations.
    MFEM_VERIFY(!IsInitialized(), "MPI should not be initialized more than once!");
    int provided;
    MPI_Init_thread(argc, argv, requested, &provided);
    MFEM_VERIFY(provided >= requested,
                "MPI could not provide the requested level of thread support!");
    // Initialize the singleton Instance.
    Instance();
  }
};

}  // namespace palace

#endif  // PALACE_UTILS_COMMUNICATION_HPP
