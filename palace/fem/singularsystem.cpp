// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "singularsystem.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include <mfem/general/forall.hpp>

#include "utils/communication.hpp"

namespace palace
{

namespace fem
{

namespace singular
{

namespace
{

bool SameCommunicator(MPI_Comm left, MPI_Comm right)
{
  int relation = MPI_UNEQUAL;
  return MPI_Comm_compare(left, right, &relation) == MPI_SUCCESS &&
         (relation == MPI_IDENT || relation == MPI_CONGRUENT);
}

std::array<HYPRE_BigInt, 2> LocalPartition(const mfem::HypreParMatrix &matrix, bool rows)
{
  const auto *partition = rows ? matrix.GetRowStarts() : matrix.GetColStarts();
  const int offset = HYPRE_AssumedPartitionCheck() ? 0 : Mpi::Rank(matrix.GetComm());
  return {partition[offset], partition[offset + 1]};
}

bool SameRowPartition(const mfem::HypreParMatrix &left, const mfem::HypreParMatrix &right)
{
  return left.GetGlobalNumRows() == right.GetGlobalNumRows() &&
         LocalPartition(left, true) == LocalPartition(right, true);
}

bool SameColumnPartition(const mfem::HypreParMatrix &left,
                         const mfem::HypreParMatrix &right)
{
  return left.GetGlobalNumCols() == right.GetGlobalNumCols() &&
         LocalPartition(left, false) == LocalPartition(right, false);
}

bool SameRowColumnPartition(const mfem::HypreParMatrix &matrix)
{
  return matrix.GetGlobalNumRows() == matrix.GetGlobalNumCols() &&
         LocalPartition(matrix, true) == LocalPartition(matrix, false);
}

void ValidateEnrichedOperatorBlocks(const mfem::HypreParMatrix &standard,
                                    const ParallelSparseOperatorBlocks &enrichment)
{
  const MPI_Comm comm = standard.GetComm();
  bool valid = enrichment.enrichment_enrichment && enrichment.standard_enrichment &&
               enrichment.enrichment_standard && standard.GetGlobalNumRows() > 0 &&
               SameRowColumnPartition(standard);
  if (valid)
  {
    const auto &ee = *enrichment.enrichment_enrichment;
    const auto &se = *enrichment.standard_enrichment;
    const auto &es = *enrichment.enrichment_standard;
    valid = SameCommunicator(comm, ee.GetComm()) && SameCommunicator(comm, se.GetComm()) &&
            SameCommunicator(comm, es.GetComm()) && ee.GetGlobalNumRows() > 0 &&
            SameRowColumnPartition(ee) && SameRowPartition(standard, se) &&
            SameColumnPartition(standard, es) && SameColumnPartition(ee, se) &&
            SameRowPartition(ee, es);
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular operator blocks have inconsistent true-DOF partitions!");
  }
}

bool ValidEssentialDofs(const mfem::Array<int> &dofs, int size)
{
  if (size < 0)
  {
    return false;
  }
  std::vector<bool> seen(size);
  for (const int dof : dofs)
  {
    if (dof < 0 || dof >= size || seen[dof])
    {
      return false;
    }
    seen[dof] = true;
  }
  return true;
}

template <std::size_t RecordSize>
std::vector<std::array<HYPRE_BigInt, RecordSize>> ExchangeBigIntRecords(
    MPI_Comm comm,
    const std::vector<std::vector<std::array<HYPRE_BigInt, RecordSize>>> &send_records)
{
  static_assert(RecordSize > 0);
  const int ranks = Mpi::Size(comm);
  bool valid = send_records.size() == static_cast<std::size_t>(ranks);
  std::vector<int> send_counts(ranks), receive_counts(ranks);
  std::int64_t send_total = 0;
  if (valid)
  {
    for (int destination = 0; destination < ranks; destination++)
    {
      if (send_records[destination].size() >
          static_cast<std::size_t>(std::numeric_limits<int>::max() / RecordSize))
      {
        valid = false;
        break;
      }
      send_counts[destination] =
          static_cast<int>(RecordSize * send_records[destination].size());
      send_total += send_counts[destination];
      if (send_total > std::numeric_limits<int>::max())
      {
        valid = false;
        break;
      }
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Parallel singular integer-record exchange exceeds MPI counts!");
  }

  MPI_Alltoall(send_counts.data(), 1, MPI_INT, receive_counts.data(), 1, MPI_INT, comm);
  std::vector<int> send_displacements(ranks), receive_displacements(ranks);
  std::int64_t receive_total = 0;
  int send_displacement = 0;
  for (int source = 0; source < ranks; source++)
  {
    send_displacements[source] = send_displacement;
    send_displacement += send_counts[source];
    if (receive_counts[source] < 0 ||
        receive_counts[source] % static_cast<int>(RecordSize) != 0 ||
        receive_total > std::numeric_limits<int>::max() - receive_counts[source])
    {
      valid = false;
      break;
    }
    receive_displacements[source] = static_cast<int>(receive_total);
    receive_total += receive_counts[source];
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Parallel singular integer-record exchange exceeds MPI counts!");
  }

  std::vector<HYPRE_BigInt> send_buffer(static_cast<std::size_t>(send_total));
  std::size_t next = 0;
  for (const auto &destination : send_records)
  {
    for (const auto &record : destination)
    {
      std::copy(record.begin(), record.end(), send_buffer.begin() + next);
      next += RecordSize;
    }
  }
  if (next != send_buffer.size())
  {
    throw std::logic_error(
        "Parallel singular integer-record exchange packed inconsistent dimensions!");
  }
  std::vector<HYPRE_BigInt> receive_buffer(static_cast<std::size_t>(receive_total));
  MPI_Alltoallv(send_buffer.data(), send_counts.data(), send_displacements.data(),
                mpi::DataType<HYPRE_BigInt>(), receive_buffer.data(), receive_counts.data(),
                receive_displacements.data(), mpi::DataType<HYPRE_BigInt>(), comm);

  std::vector<std::array<HYPRE_BigInt, RecordSize>> result(
      static_cast<std::size_t>(receive_total) / RecordSize);
  for (std::size_t record = 0; record < result.size(); record++)
  {
    std::copy(receive_buffer.begin() + record * RecordSize,
              receive_buffer.begin() + (record + 1) * RecordSize, result[record].begin());
  }
  return result;
}

std::map<HYPRE_BigInt, std::vector<std::size_t>>
BuildLocalFeatureMembership(MPI_Comm comm,
                            const std::vector<std::vector<std::size_t>> &local_membership,
                            const TrueDofMap &numbering, std::size_t number_features,
                            const std::vector<HYPRE_BigInt> &requested_true_dofs)
{
  const int ranks = Mpi::Size(comm);
  const int rank = Mpi::Rank(comm);
  bool valid = number_features > 0 &&
               number_features <=
                   static_cast<std::size_t>(std::numeric_limits<HYPRE_BigInt>::max()) &&
               local_membership.size() == numbering.local_to_true.size() &&
               local_membership.size() == numbering.owner.size() &&
               numbering.global_size > 0 && numbering.owned_offset >= 0 &&
               numbering.owned_size >= 0 &&
               numbering.owned_offset <= numbering.global_size - numbering.owned_size;
  std::map<HYPRE_BigInt, std::vector<std::size_t>> result;
  using MembershipRecord = std::array<HYPRE_BigInt, 3>;
  std::vector<std::vector<MembershipRecord>> send_records(ranks);
  if (valid)
  {
    for (std::size_t i = 0; i < local_membership.size(); i++)
    {
      const auto &membership = local_membership[i];
      if (membership.empty() || !std::is_sorted(membership.begin(), membership.end()) ||
          std::adjacent_find(membership.begin(), membership.end()) != membership.end() ||
          membership.back() >= number_features || numbering.local_to_true[i] < 0 ||
          numbering.local_to_true[i] >= numbering.global_size || numbering.owner[i] < 0 ||
          numbering.owner[i] >= ranks)
      {
        valid = false;
        break;
      }
      const auto [entry, inserted] = result.emplace(numbering.local_to_true[i], membership);
      if (!inserted && entry->second != membership)
      {
        valid = false;
        break;
      }
      for (std::size_t feature : membership)
      {
        send_records[numbering.owner[i]].push_back(
            {numbering.local_to_true[i], rank, static_cast<HYPRE_BigInt>(feature)});
      }
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular feature patches received invalid DOF membership!");
  }

  constexpr int record_size = std::tuple_size_v<MembershipRecord>;
  std::vector<int> send_counts(ranks), receive_counts(ranks);
  std::int64_t send_total = 0;
  for (int destination = 0; destination < ranks; destination++)
  {
    if (send_records[destination].size() >
        static_cast<std::size_t>(std::numeric_limits<int>::max() / record_size))
    {
      valid = false;
      break;
    }
    send_counts[destination] =
        static_cast<int>(record_size * send_records[destination].size());
    send_total += send_counts[destination];
    if (send_total > std::numeric_limits<int>::max())
    {
      valid = false;
      break;
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Parallel singular feature membership exceeds MPI integer counts!");
  }

  MPI_Alltoall(send_counts.data(), 1, MPI_INT, receive_counts.data(), 1, MPI_INT, comm);
  std::vector<int> send_displacements(ranks), receive_displacements(ranks);
  std::int64_t receive_total = 0;
  int send_displacement = 0;
  for (int source = 0; source < ranks; source++)
  {
    send_displacements[source] = send_displacement;
    send_displacement += send_counts[source];
    if (receive_counts[source] < 0 || receive_counts[source] % record_size != 0 ||
        receive_total > std::numeric_limits<int>::max() - receive_counts[source])
    {
      valid = false;
      break;
    }
    receive_displacements[source] = static_cast<int>(receive_total);
    receive_total += receive_counts[source];
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::overflow_error(
        "Parallel singular feature membership exceeds MPI integer counts!");
  }

  std::vector<HYPRE_BigInt> send_buffer(static_cast<std::size_t>(send_total));
  std::size_t next = 0;
  for (const auto &destination : send_records)
  {
    for (const auto &record : destination)
    {
      std::copy(record.begin(), record.end(), send_buffer.begin() + next);
      next += record_size;
    }
  }
  if (next != send_buffer.size())
  {
    throw std::logic_error(
        "Parallel singular feature membership packed inconsistent dimensions!");
  }
  std::vector<HYPRE_BigInt> receive_buffer(static_cast<std::size_t>(receive_total));
  MPI_Alltoallv(send_buffer.data(), send_counts.data(), send_displacements.data(),
                mpi::DataType<HYPRE_BigInt>(), receive_buffer.data(), receive_counts.data(),
                receive_displacements.data(), mpi::DataType<HYPRE_BigInt>(), comm);

  std::map<HYPRE_BigInt, std::map<int, std::set<std::size_t>>> received;
  for (std::size_t offset = 0; offset < receive_buffer.size(); offset += record_size)
  {
    const HYPRE_BigInt true_dof = receive_buffer[offset];
    const HYPRE_BigInt origin = receive_buffer[offset + 1];
    const HYPRE_BigInt feature = receive_buffer[offset + 2];
    if (true_dof < numbering.owned_offset ||
        true_dof >= numbering.owned_offset + numbering.owned_size || origin < 0 ||
        origin >= ranks || feature < 0 ||
        feature >= static_cast<HYPRE_BigInt>(number_features) ||
        !received[true_dof][static_cast<int>(origin)]
             .insert(static_cast<std::size_t>(feature))
             .second)
    {
      valid = false;
      break;
    }
  }
  if (valid)
  {
    for (HYPRE_BigInt local = 0; local < numbering.owned_size; local++)
    {
      const HYPRE_BigInt true_dof = numbering.owned_offset + local;
      const auto authoritative = result.find(true_dof);
      const auto occurrences = received.find(true_dof);
      if (authoritative == result.end() || occurrences == received.end())
      {
        valid = false;
        break;
      }
      const std::set<std::size_t> expected(authoritative->second.begin(),
                                           authoritative->second.end());
      for (const auto &[origin, membership] : occurrences->second)
      {
        (void)origin;
        if (membership != expected)
        {
          valid = false;
          break;
        }
      }
      if (!valid)
      {
        break;
      }
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Shared singular DOFs have inconsistent feature membership across ranks!");
  }

  std::vector<HYPRE_BigInt> owned_sizes(ranks);
  Mpi::Allgather(1, &numbering.owned_size, owned_sizes.data(), comm);
  std::vector<HYPRE_BigInt> owned_partition(ranks + 1);
  for (int owner = 0; owner < ranks; owner++)
  {
    if (owned_sizes[owner] < 0 ||
        owned_partition[owner] > numbering.global_size - owned_sizes[owner])
    {
      valid = false;
      break;
    }
    owned_partition[owner + 1] = owned_partition[owner] + owned_sizes[owner];
  }
  valid = valid && owned_partition.back() == numbering.global_size &&
          numbering.owned_offset == owned_partition[rank] &&
          numbering.owned_size == owned_sizes[rank];
  std::vector<HYPRE_BigInt> requests = requested_true_dofs;
  std::sort(requests.begin(), requests.end());
  requests.erase(std::unique(requests.begin(), requests.end()), requests.end());
  if (!requests.empty() &&
      (requests.front() < 0 || requests.back() >= numbering.global_size))
  {
    valid = false;
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular feature membership query has inconsistent ownership!");
  }

  using QueryRecord = std::array<HYPRE_BigInt, 2>;
  using ResponseRecord = std::array<HYPRE_BigInt, 2>;
  std::vector<std::vector<QueryRecord>> query_send(ranks);
  for (HYPRE_BigInt true_dof : requests)
  {
    if (result.find(true_dof) != result.end())
    {
      continue;
    }
    const auto upper =
        std::upper_bound(owned_partition.begin(), owned_partition.end(), true_dof);
    const int owner = static_cast<int>(std::distance(owned_partition.begin(), upper)) - 1;
    if (owner < 0 || owner >= ranks)
    {
      valid = false;
      break;
    }
    query_send[owner].push_back({true_dof, rank});
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular feature membership query has no valid owner!");
  }
  const auto received_queries = ExchangeBigIntRecords(comm, query_send);

  std::vector<std::vector<ResponseRecord>> response_send(ranks);
  std::set<std::pair<HYPRE_BigInt, int>> unique_queries;
  for (const auto &query : received_queries)
  {
    const HYPRE_BigInt true_dof = query[0];
    const HYPRE_BigInt requester = query[1];
    const auto membership = result.find(true_dof);
    if (true_dof < numbering.owned_offset ||
        true_dof >= numbering.owned_offset + numbering.owned_size || requester < 0 ||
        requester >= ranks || membership == result.end() ||
        !unique_queries.emplace(true_dof, static_cast<int>(requester)).second)
    {
      valid = false;
      break;
    }
    for (std::size_t feature : membership->second)
    {
      response_send[requester].push_back({true_dof, static_cast<HYPRE_BigInt>(feature)});
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular feature membership owner received an invalid query!");
  }
  const auto responses = ExchangeBigIntRecords(comm, response_send);

  std::map<HYPRE_BigInt, std::set<std::size_t>> fetched;
  for (const auto &response : responses)
  {
    const HYPRE_BigInt true_dof = response[0];
    const HYPRE_BigInt feature = response[1];
    if (!std::binary_search(requests.begin(), requests.end(), true_dof) || feature < 0 ||
        feature >= static_cast<HYPRE_BigInt>(number_features) ||
        !fetched[true_dof].insert(static_cast<std::size_t>(feature)).second)
    {
      valid = false;
      break;
    }
  }
  if (valid)
  {
    for (HYPRE_BigInt true_dof : requests)
    {
      if (result.find(true_dof) != result.end())
      {
        continue;
      }
      const auto membership = fetched.find(true_dof);
      if (membership == fetched.end() || membership->second.empty())
      {
        valid = false;
        break;
      }
      result.emplace(true_dof, std::vector<std::size_t>(membership->second.begin(),
                                                        membership->second.end()));
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular feature membership query received an invalid response!");
  }
  return result;
}

}  // namespace

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichedOperator(const mfem::HypreParMatrix &standard_standard,
                              const ParallelSparseOperatorBlocks &enrichment)
{
  ValidateEnrichedOperatorBlocks(standard_standard, enrichment);
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = &standard_standard;
  blocks(0, 1) = enrichment.standard_enrichment.get();
  blocks(1, 0) = enrichment.enrichment_standard.get();
  blocks(1, 1) = enrichment.enrichment_enrichment.get();
  return std::unique_ptr<mfem::HypreParMatrix>(mfem::HypreParMatrixFromBlocks(blocks));
}

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichedGradient(const mfem::HypreParMatrix &standard_gradient,
                              const mfem::HypreParMatrix &enrichment_gradient)
{
  const MPI_Comm comm = standard_gradient.GetComm();
  bool valid = SameCommunicator(comm, enrichment_gradient.GetComm()) &&
               standard_gradient.GetGlobalNumRows() > 0 &&
               standard_gradient.GetGlobalNumCols() > 0 &&
               enrichment_gradient.GetGlobalNumRows() > 0 &&
               enrichment_gradient.GetGlobalNumCols() > 0;
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular gradient blocks have inconsistent communicators or sizes!");
  }

  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks = nullptr;
  blocks(0, 0) = &standard_gradient;
  blocks(1, 1) = &enrichment_gradient;
  return std::unique_ptr<mfem::HypreParMatrix>(mfem::HypreParMatrixFromBlocks(blocks));
}

ParallelConstrainedOperatorBlocks BuildParallelConstrainedOperatorBlocks(
    const mfem::HypreParMatrix &standard_standard,
    const ParallelSparseOperatorBlocks &enrichment,
    const mfem::Array<int> &standard_essential_true_dofs,
    const mfem::Array<int> &enrichment_essential_true_dofs)
{
  ValidateEnrichedOperatorBlocks(standard_standard, enrichment);
  const int standard_local_size = standard_standard.Height();
  const int enrichment_local_size = enrichment.enrichment_enrichment->Height();
  bool valid = ValidEssentialDofs(standard_essential_true_dofs, standard_local_size) &&
               ValidEssentialDofs(enrichment_essential_true_dofs, enrichment_local_size);
  Mpi::GlobalAnd(1, &valid, standard_standard.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular essential true DOFs are inconsistent with the blocks!");
  }

  ParallelConstrainedOperatorBlocks constrained{
      std::make_unique<mfem::HypreParMatrix>(standard_standard),
      std::make_unique<mfem::HypreParMatrix>(*enrichment.standard_enrichment),
      std::make_unique<mfem::HypreParMatrix>(*enrichment.enrichment_standard),
      std::make_unique<mfem::HypreParMatrix>(*enrichment.enrichment_enrichment)};

  std::unique_ptr<mfem::HypreParMatrix> discarded(
      constrained.standard_standard->EliminateRowsCols(standard_essential_true_dofs));
  discarded.reset(
      constrained.enrichment_enrichment->EliminateRowsCols(enrichment_essential_true_dofs));

  constrained.standard_enrichment->EliminateRows(standard_essential_true_dofs);
  discarded.reset(
      constrained.standard_enrichment->EliminateCols(enrichment_essential_true_dofs));
  constrained.enrichment_standard->EliminateRows(enrichment_essential_true_dofs);
  discarded.reset(
      constrained.enrichment_standard->EliminateCols(standard_essential_true_dofs));
  return constrained;
}

ParallelCoupledPatch
BuildParallelCoupledPatch(const mfem::HypreParMatrix &constrained,
                          const mfem::HypreParMatrix &standard_enrichment,
                          int standard_local_size)
{
  const MPI_Comm comm = constrained.GetComm();
  const int enrichment_local_size = constrained.Height() - standard_local_size;
  bool valid =
      SameCommunicator(comm, standard_enrichment.GetComm()) &&
      SameRowColumnPartition(constrained) && standard_local_size >= 0 &&
      enrichment_local_size >= 0 && standard_enrichment.Height() == standard_local_size &&
      standard_enrichment.Width() == enrichment_local_size &&
      standard_enrichment.GetGlobalNumRows() + standard_enrichment.GetGlobalNumCols() ==
          constrained.GetGlobalNumRows();
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Cannot build a coupled singular patch from inconsistent operator blocks!");
  }

  mfem::SparseMatrix diagonal, off_diagonal;
  HYPRE_BigInt *off_diagonal_columns = nullptr;
  standard_enrichment.GetDiag(diagonal);
  standard_enrichment.GetOffd(off_diagonal, off_diagonal_columns);
  const auto *diagonal_offsets = diagonal.HostReadI();
  const auto *diagonal_values = diagonal.HostReadData();
  const auto *off_diagonal_offsets = off_diagonal.HostReadI();
  const auto *off_diagonal_values = off_diagonal.HostReadData();

  mfem::Array<int> patch_true_dofs;
  patch_true_dofs.Reserve(standard_local_size + enrichment_local_size);
  HYPRE_BigInt local_standard_dofs = 0;
  for (int i = 0; i < standard_local_size; i++)
  {
    bool coupled = false;
    for (int j = diagonal_offsets[i]; j < diagonal_offsets[i + 1]; j++)
    {
      coupled = coupled || diagonal_values[j] != 0.0;
    }
    for (int j = off_diagonal_offsets[i]; j < off_diagonal_offsets[i + 1]; j++)
    {
      coupled = coupled || off_diagonal_values[j] != 0.0;
    }
    if (coupled)
    {
      patch_true_dofs.Append(i);
      local_standard_dofs++;
    }
  }
  for (int i = 0; i < enrichment_local_size; i++)
  {
    patch_true_dofs.Append(standard_local_size + i);
  }

  HYPRE_BigInt global_standard_dofs = local_standard_dofs;
  HYPRE_BigInt global_enrichment_dofs = enrichment_local_size;
  Mpi::GlobalSum(1, &global_standard_dofs, comm);
  Mpi::GlobalSum(1, &global_enrichment_dofs, comm);
  if (global_enrichment_dofs == 0)
  {
    throw std::runtime_error("Coupled singular patch has no enrichment true DOFs!");
  }

#if MFEM_HYPRE_VERSION >= 21800
  std::unique_ptr<mfem::HypreParMatrix> patch(
      constrained.ExtractSubmatrix(patch_true_dofs));
#else
  std::unique_ptr<mfem::HypreParMatrix> patch;
#endif
  valid = patch && patch->Height() == patch_true_dofs.Size() &&
          patch->Width() == patch_true_dofs.Size() &&
          patch->GetGlobalNumRows() == global_standard_dofs + global_enrichment_dofs &&
          SameRowColumnPartition(*patch);
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::runtime_error(
        "Failed to extract the parallel coupled singular patch matrix!");
  }
  return {std::move(patch), std::move(patch_true_dofs), global_standard_dofs,
          global_enrichment_dofs};
}

ParallelFeaturePatches BuildParallelFeaturePatches(
    const mfem::HypreParMatrix &constrained,
    const mfem::HypreParMatrix &standard_enrichment, int standard_local_size,
    const std::vector<std::vector<std::size_t>> &local_enrichment_features,
    const TrueDofMap &enrichment_numbering, std::size_t number_features)
{
  const MPI_Comm comm = constrained.GetComm();
  const int enrichment_local_size = constrained.Height() - standard_local_size;
  bool valid =
      SameCommunicator(comm, standard_enrichment.GetComm()) &&
      SameRowColumnPartition(constrained) && standard_local_size >= 0 &&
      enrichment_local_size >= 0 && standard_enrichment.Height() == standard_local_size &&
      standard_enrichment.Width() == enrichment_local_size &&
      standard_enrichment.GetGlobalNumRows() + standard_enrichment.GetGlobalNumCols() ==
          constrained.GetGlobalNumRows() &&
      enrichment_numbering.owned_size == enrichment_local_size &&
      enrichment_numbering.global_size == standard_enrichment.GetGlobalNumCols();
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Cannot build singular feature patches from inconsistent operator blocks!");
  }

  mfem::SparseMatrix diagonal, off_diagonal;
  HYPRE_BigInt *off_diagonal_columns = nullptr;
  standard_enrichment.GetDiag(diagonal);
  standard_enrichment.GetOffd(off_diagonal, off_diagonal_columns);
  std::vector<HYPRE_BigInt> requested_membership;
  requested_membership.reserve(off_diagonal.Width());
  for (int column = 0; column < off_diagonal.Width(); column++)
  {
    requested_membership.push_back(off_diagonal_columns[column]);
  }
  const auto local_membership =
      BuildLocalFeatureMembership(comm, local_enrichment_features, enrichment_numbering,
                                  number_features, requested_membership);

  std::vector<std::vector<int>> local_enrichment(number_features);
  std::vector<int> enrichment_multiplicity(enrichment_local_size);
  const int rank = Mpi::Rank(comm);
  for (std::size_t i = 0; i < local_enrichment_features.size(); i++)
  {
    if (enrichment_numbering.owner[i] != rank)
    {
      continue;
    }
    const HYPRE_BigInt local =
        enrichment_numbering.local_to_true[i] - enrichment_numbering.owned_offset;
    if (local < 0 || local >= enrichment_numbering.owned_size)
    {
      throw std::invalid_argument(
          "An owned singular enrichment DOF is outside its parallel partition!");
    }
    for (std::size_t feature : local_enrichment_features[i])
    {
      local_enrichment[feature].push_back(static_cast<int>(local));
      enrichment_multiplicity[local]++;
    }
  }
  for (auto &indices : local_enrichment)
  {
    std::sort(indices.begin(), indices.end());
    if (std::adjacent_find(indices.begin(), indices.end()) != indices.end())
    {
      throw std::logic_error(
          "A singular enrichment DOF occurs twice in one feature patch!");
    }
  }

  int minimum_multiplicity = std::numeric_limits<int>::max();
  int maximum_multiplicity = 0;
  for (int multiplicity : enrichment_multiplicity)
  {
    minimum_multiplicity = std::min(minimum_multiplicity, multiplicity);
    maximum_multiplicity = std::max(maximum_multiplicity, multiplicity);
  }
  Mpi::GlobalMin(1, &minimum_multiplicity, comm);
  Mpi::GlobalMax(1, &maximum_multiplicity, comm);
  if (minimum_multiplicity < 1 || maximum_multiplicity < minimum_multiplicity)
  {
    throw std::invalid_argument(
        "Singular feature patches do not cover every enrichment true DOF!");
  }

  const auto *diagonal_offsets = diagonal.HostReadI();
  const auto *diagonal_columns = diagonal.HostReadJ();
  const auto *diagonal_values = diagonal.HostReadData();
  const auto *off_diagonal_offsets = off_diagonal.HostReadI();
  const auto *off_diagonal_column_indices = off_diagonal.HostReadJ();
  const auto *off_diagonal_values = off_diagonal.HostReadData();

  std::vector<std::vector<int>> local_standard(number_features);
  for (int row = 0; row < standard_local_size; row++)
  {
    std::set<std::size_t> row_features;
    for (int entry = diagonal_offsets[row]; entry < diagonal_offsets[row + 1]; entry++)
    {
      if (diagonal_values[entry] == 0.0)
      {
        continue;
      }
      const int column = diagonal_columns[entry];
      if (column < 0 || column >= enrichment_local_size)
      {
        throw std::runtime_error(
            "Singular standard-enrichment diagonal block has an invalid column!");
      }
      const HYPRE_BigInt global_column = enrichment_numbering.owned_offset + column;
      const auto membership = local_membership.find(global_column);
      if (membership == local_membership.end())
      {
        throw std::runtime_error(
            "A coupled singular enrichment DOF has no feature membership!");
      }
      row_features.insert(membership->second.begin(), membership->second.end());
    }
    for (int entry = off_diagonal_offsets[row]; entry < off_diagonal_offsets[row + 1];
         entry++)
    {
      if (off_diagonal_values[entry] == 0.0)
      {
        continue;
      }
      const int column = off_diagonal_column_indices[entry];
      if (column < 0 || column >= off_diagonal.Width())
      {
        throw std::runtime_error(
            "Singular standard-enrichment off-diagonal block has an invalid column!");
      }
      const auto membership = local_membership.find(off_diagonal_columns[column]);
      if (membership == local_membership.end())
      {
        throw std::runtime_error(
            "An off-rank singular enrichment DOF has no feature membership!");
      }
      row_features.insert(membership->second.begin(), membership->second.end());
    }
    for (std::size_t feature : row_features)
    {
      local_standard[feature].push_back(row);
    }
  }

  ParallelFeaturePatches result;
  result.sum_global_standard_dofs = 0;
  result.sum_global_enrichment_dofs = 0;
  result.maximum_global_standard_dofs = 0;
  result.maximum_global_enrichment_dofs = 0;
  result.minimum_enrichment_multiplicity = minimum_multiplicity;
  result.maximum_enrichment_multiplicity = maximum_multiplicity;
  result.patches.reserve(number_features);
  for (std::size_t feature = 0; feature < number_features; feature++)
  {
    mfem::Array<int> patch_true_dofs;
    patch_true_dofs.Reserve(local_standard[feature].size() +
                            local_enrichment[feature].size());
    for (int dof : local_standard[feature])
    {
      patch_true_dofs.Append(dof);
    }
    for (int dof : local_enrichment[feature])
    {
      patch_true_dofs.Append(standard_local_size + dof);
    }

    HYPRE_BigInt global_standard_dofs = local_standard[feature].size();
    HYPRE_BigInt global_enrichment_dofs = local_enrichment[feature].size();
    Mpi::GlobalSum(1, &global_standard_dofs, comm);
    Mpi::GlobalSum(1, &global_enrichment_dofs, comm);
    if (global_enrichment_dofs == 0)
    {
      throw std::runtime_error("A straight singular feature has no enrichment true DOFs!");
    }

#if MFEM_HYPRE_VERSION >= 21800
    std::unique_ptr<mfem::HypreParMatrix> patch(
        constrained.ExtractSubmatrix(patch_true_dofs));
#else
    std::unique_ptr<mfem::HypreParMatrix> patch;
#endif
    valid = patch && patch->Height() == patch_true_dofs.Size() &&
            patch->Width() == patch_true_dofs.Size() &&
            patch->GetGlobalNumRows() == global_standard_dofs + global_enrichment_dofs &&
            SameRowColumnPartition(*patch);
    Mpi::GlobalAnd(1, &valid, comm);
    if (!valid)
    {
      throw std::runtime_error(
          "Failed to extract a parallel singular feature patch matrix!");
    }

    result.sum_global_standard_dofs += global_standard_dofs;
    result.sum_global_enrichment_dofs += global_enrichment_dofs;
    result.maximum_global_standard_dofs =
        std::max(result.maximum_global_standard_dofs, global_standard_dofs);
    result.maximum_global_enrichment_dofs =
        std::max(result.maximum_global_enrichment_dofs, global_enrichment_dofs);
    result.patches.push_back({feature, std::move(patch), std::move(patch_true_dofs),
                              global_standard_dofs, global_enrichment_dofs});
  }
  return result;
}

SymmetricDiagonalScaling::SymmetricDiagonalScaling(const mfem::HypreParMatrix &unscaled)
  : matrix(std::make_unique<mfem::HypreParMatrix>(unscaled)),
    scaled_diagonal_minimum(std::numeric_limits<double>::infinity()),
    scaled_diagonal_maximum(0.0)
{
  const MPI_Comm comm = unscaled.GetComm();
  coordinate_scaling.SetSize(unscaled.Height());
  unscaled.GetDiag(coordinate_scaling);
  bool valid = unscaled.Height() == unscaled.Width() &&
               unscaled.GetGlobalNumRows() == unscaled.GetGlobalNumCols();
  if (valid)
  {
    auto *scaling = coordinate_scaling.HostReadWrite();
    for (int i = 0; i < coordinate_scaling.Size(); i++)
    {
      const double diagonal = scaling[i];
      if (!std::isfinite(diagonal) || !(diagonal > 0.0))
      {
        valid = false;
        break;
      }
      scaling[i] = 1.0 / std::sqrt(diagonal);
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Symmetric diagonal scaling requires a square positive-diagonal operator!");
  }

  mfem::HypreParVector scaling(comm, unscaled.GetGlobalNumRows(), coordinate_scaling, 0,
                               unscaled.GetRowStarts());
  scaling.HypreRead();
  matrix->HypreReadWrite();
  int error = HYPRE_ParCSRMatrixDiagScale(*matrix, scaling, scaling);
  Mpi::GlobalMax(1, &error, comm);
  if (error != 0)
  {
    throw std::runtime_error("HYPRE failed to symmetrically scale the singular system!");
  }

  Vector scaled_diagonal;
  matrix->GetDiag(scaled_diagonal);
  valid = scaled_diagonal.Size() == matrix->Height();
  if (valid)
  {
    const auto *diagonal = scaled_diagonal.HostRead();
    for (int i = 0; i < scaled_diagonal.Size(); i++)
    {
      if (!std::isfinite(diagonal[i]) || !(diagonal[i] > 0.0))
      {
        valid = false;
        break;
      }
      scaled_diagonal_minimum = std::min(scaled_diagonal_minimum, diagonal[i]);
      scaled_diagonal_maximum = std::max(scaled_diagonal_maximum, diagonal[i]);
    }
  }
  Mpi::GlobalAnd(1, &valid, comm);
  Mpi::GlobalMin(1, &scaled_diagonal_minimum, comm);
  Mpi::GlobalMax(1, &scaled_diagonal_maximum, comm);
  constexpr double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
  if (!valid || !std::isfinite(scaled_diagonal_minimum) ||
      !std::isfinite(scaled_diagonal_maximum) ||
      std::abs(scaled_diagonal_minimum - 1.0) > tolerance ||
      std::abs(scaled_diagonal_maximum - 1.0) > tolerance)
  {
    throw std::runtime_error(
        "Symmetrically scaled singular system does not have a unit diagonal!");
  }
}

void SymmetricDiagonalScaling::Apply(const Vector &input, Vector &output,
                                     bool inverse) const
{
  MFEM_VERIFY(input.Size() == coordinate_scaling.Size(),
              "Diagonal scaling received an inconsistent vector size!");
  output.SetSize(input.Size());
  const bool use_dev = input.UseDevice() || output.UseDevice();
  const auto *x = input.Read(use_dev);
  const auto *s = coordinate_scaling.Read(use_dev);
  auto *y = output.Write(use_dev);
  mfem::forall_switch(use_dev, input.Size(), [=] MFEM_HOST_DEVICE(int i)
                      { y[i] = inverse ? x[i] / s[i] : s[i] * x[i]; });
}

void SymmetricDiagonalScaling::ScaleRHS(const Vector &unscaled, Vector &scaled) const
{
  Apply(unscaled, scaled, false);
}

void SymmetricDiagonalScaling::ScaleInitialGuess(const Vector &unscaled,
                                                 Vector &scaled) const
{
  Apply(unscaled, scaled, true);
}

void SymmetricDiagonalScaling::RecoverSolution(const Vector &scaled, Vector &unscaled) const
{
  Apply(scaled, unscaled, false);
}

ParallelDirichletSystem
BuildParallelDirichletSystem(std::unique_ptr<mfem::HypreParMatrix> &&matrix,
                             int standard_local_size,
                             const mfem::Array<int> &standard_essential_true_dofs,
                             const mfem::Array<int> &enrichment_essential_true_dofs)
{
  if (!matrix)
  {
    throw std::invalid_argument(
        "Cannot eliminate essential DOFs from an empty singular operator!");
  }
  const MPI_Comm comm = matrix->GetComm();
  const int enrichment_local_size = matrix->Height() - standard_local_size;
  bool valid = matrix->Height() == matrix->Width() && standard_local_size >= 0 &&
               standard_local_size <= matrix->Height() &&
               ValidEssentialDofs(standard_essential_true_dofs, standard_local_size) &&
               ValidEssentialDofs(enrichment_essential_true_dofs, enrichment_local_size);
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular essential true DOFs are inconsistent with the operator!");
  }

  mfem::Array<int> essential;
  essential.Reserve(standard_essential_true_dofs.Size() +
                    enrichment_essential_true_dofs.Size());
  for (const int dof : standard_essential_true_dofs)
  {
    essential.Append(dof);
  }
  for (const int dof : enrichment_essential_true_dofs)
  {
    essential.Append(standard_local_size + dof);
  }
  essential.Sort();

  std::unique_ptr<mfem::HypreParMatrix> eliminated(matrix->EliminateRowsCols(essential));
  if (!eliminated)
  {
    throw std::runtime_error("Failed to retain the eliminated singular operator entries!");
  }
  return {std::move(matrix), std::move(eliminated), std::move(essential)};
}

void ParallelDirichletSystem::EliminateRHS(const mfem::Vector &x, mfem::Vector &b) const
{
  if (!constrained || !eliminated)
  {
    throw std::logic_error("The parallel singular Dirichlet system is incomplete!");
  }
  const MPI_Comm comm = constrained->GetComm();
  bool valid = SameCommunicator(comm, eliminated->GetComm()) &&
               constrained->Height() == constrained->Width() &&
               eliminated->Height() == constrained->Height() &&
               eliminated->Width() == constrained->Width() &&
               x.Size() == constrained->Width() && b.Size() == constrained->Height() &&
               ValidEssentialDofs(essential_true_dofs, constrained->Height());
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Cannot eliminate a singular-system RHS with inconsistent dimensions!");
  }
  constrained->EliminateBC(*eliminated, essential_true_dofs, x, b);
}

}  // namespace singular

}  // namespace fem

}  // namespace palace
