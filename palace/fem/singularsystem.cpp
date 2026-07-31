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
#include <Eigen/Eigenvalues>
#include <fmt/format.h>
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

std::vector<HYPRE_BigInt> BuildPartition(MPI_Comm comm, HYPRE_BigInt local_offset,
                                         HYPRE_BigInt local_size, HYPRE_BigInt global_size)
{
  const int ranks = Mpi::Size(comm);
  std::vector<HYPRE_BigInt> offsets(ranks), sizes(ranks), global_sizes(ranks);
  Mpi::Allgather(1, &local_offset, offsets.data(), comm);
  Mpi::Allgather(1, &local_size, sizes.data(), comm);
  Mpi::Allgather(1, &global_size, global_sizes.data(), comm);
  HYPRE_BigInt expected_offset = 0;
  bool valid = global_size > 0 && std::all_of(global_sizes.begin(), global_sizes.end(),
                                              [global_size](HYPRE_BigInt size)
                                              { return size == global_size; });
  for (int rank = 0; rank < ranks; rank++)
  {
    if (offsets[rank] != expected_offset || sizes[rank] < 0 ||
        sizes[rank] > std::numeric_limits<HYPRE_BigInt>::max() - expected_offset)
    {
      valid = false;
      break;
    }
    expected_offset += sizes[rank];
  }
  if (!valid || expected_offset != global_size)
  {
    throw std::invalid_argument("Parallel singular true-DOF partition is inconsistent!");
  }
  if (HYPRE_AssumedPartitionCheck())
  {
    return {local_offset, local_offset + local_size};
  }
  std::vector<HYPRE_BigInt> partition(ranks + 1);
  std::copy(offsets.begin(), offsets.end(), partition.begin());
  partition.back() = global_size;
  return partition;
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

int ParallelElementPatchInverse::DecodeSignedDof(int dof)
{
  return dof >= 0 ? dof : -1 - dof;
}

double ParallelElementPatchInverse::DecodeSignedDofSign(int dof)
{
  return dof >= 0 ? 1.0 : -1.0;
}

namespace
{

double FrobeniusNorm(const mfem::DenseMatrix &matrix)
{
  long double norm_squared = 0.0L;
  for (int i = 0; i < matrix.Height(); i++)
  {
    for (int j = 0; j < matrix.Width(); j++)
    {
      const long double value = matrix(i, j);
      norm_squared += value * value;
    }
  }
  return std::sqrt(static_cast<double>(norm_squared));
}

struct PatchSpectrum
{
  mfem::DenseMatrix inverse;
  int discarded_modes = 0;
  double minimum_eigenvalue = std::numeric_limits<double>::infinity();
  double maximum_eigenvalue = 0.0;
  double resolution = 0.0;
  bool materially_indefinite = false;
};

PatchSpectrum BuildPatchPseudoinverse(mfem::DenseMatrix &matrix,
                                      const mfem::DenseMatrix &estimated_absolute_error)
{
  const int size = matrix.Height();
  MFEM_VERIFY(size > 0 && matrix.Width() == size &&
                  estimated_absolute_error.Height() == size &&
                  estimated_absolute_error.Width() == size,
              "Invalid matrix for coupled singular element patch pseudoinverse!");

  const Eigen::Map<
      const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
      eigen_matrix(matrix.Data(), size, size);
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigensystem(eigen_matrix);
  if (eigensystem.info() != Eigen::Success)
  {
    throw std::runtime_error(
        "Failed to diagonalize a coupled singular element patch matrix!");
  }
  const auto &eigenvalues = eigensystem.eigenvalues();
  const auto &eigenvectors = eigensystem.eigenvectors();

  PatchSpectrum result;
  result.inverse.SetSize(size);
  result.inverse = 0.0;
  result.minimum_eigenvalue = eigenvalues[0];
  result.maximum_eigenvalue =
      std::max(std::abs(eigenvalues[0]), std::abs(eigenvalues[size - 1]));
  const double roundoff =
      256.0 * std::numeric_limits<double>::epsilon() * size *
      std::max(std::numeric_limits<double>::min(), result.maximum_eigenvalue);
  result.resolution = std::max(roundoff, 2.0 * FrobeniusNorm(estimated_absolute_error));
  result.materially_indefinite = result.minimum_eigenvalue < -result.resolution;

  for (int mode = 0; mode < size; mode++)
  {
    const double eigenvalue = eigenvalues[mode];
    if (eigenvalue <= result.resolution)
    {
      result.discarded_modes++;
      continue;
    }
    const double inverse_eigenvalue = 1.0 / eigenvalue;
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        result.inverse(i, j) +=
            inverse_eigenvalue * eigenvectors(i, mode) * eigenvectors(j, mode);
      }
    }
  }
  return result;
}

}  // namespace

void ParallelElementPatchInverse::ApplyTrueWeight(const Vector &input, Vector &output) const
{
  MFEM_VERIFY(input.Size() == height,
              "Coupled singular element patches received an invalid true vector!");
  output.SetSize(height);
  const auto *input_data = input.HostRead();
  const auto *weight_data = true_weight.HostRead();
  auto *output_data = output.HostWrite();
  for (int i = 0; i < height; i++)
  {
    output_data[i] = weight_data[i] * input_data[i];
  }
}

ParallelElementPatchInverse::ParallelElementPatchInverse(
    const mfem::ParFiniteElementSpace &standard_fespace,
    const TrueDofMap &enrichment_numbering,
    const std::vector<LocalNDElementPatchMatrices> &element_matrices,
    double stiffness_coefficient, double mass_coefficient,
    const mfem::Array<int> &standard_essential_true_dofs,
    const mfem::Array<int> &enrichment_essential_true_dofs)
  : Operator(standard_fespace.GetTrueVSize() +
             static_cast<int>(enrichment_numbering.owned_size)),
    standard_prolongation(standard_fespace.Dof_TrueDof_Matrix()),
    enrichment_prolongation(BuildParallelEnrichmentProlongation(standard_fespace.GetComm(),
                                                                enrichment_numbering)),
    standard_true_size(standard_fespace.GetTrueVSize()),
    standard_local_size(standard_fespace.GetVSize()),
    enrichment_local_size(static_cast<int>(enrichment_numbering.local_size))
{
  bool valid = standard_prolongation && !element_matrices.empty() &&
               std::isfinite(stiffness_coefficient) && std::isfinite(mass_coefficient) &&
               stiffness_coefficient >= 0.0 && mass_coefficient > 0.0 &&
               enrichment_numbering.owned_size <= std::numeric_limits<int>::max() &&
               enrichment_numbering.local_size <= std::numeric_limits<int>::max() &&
               standard_prolongation->Height() == standard_local_size &&
               standard_prolongation->Width() == standard_true_size &&
               enrichment_prolongation->Height() == enrichment_local_size &&
               enrichment_prolongation->Width() ==
                   static_cast<int>(enrichment_numbering.owned_size) &&
               ValidEssentialDofs(standard_essential_true_dofs, standard_true_size) &&
               ValidEssentialDofs(enrichment_essential_true_dofs,
                                  static_cast<int>(enrichment_numbering.owned_size));
  Mpi::GlobalAnd(1, &valid, standard_fespace.GetComm());
  if (!valid)
  {
    throw std::invalid_argument(
        "Coupled singular element patches received inconsistent spaces or coefficients!");
  }

  Vector standard_true_essential(standard_true_size);
  Vector enrichment_true_essential(static_cast<int>(enrichment_numbering.owned_size));
  standard_true_essential = 0.0;
  enrichment_true_essential = 0.0;
  for (int dof : standard_essential_true_dofs)
  {
    standard_true_essential[dof] = 1.0;
  }
  for (int dof : enrichment_essential_true_dofs)
  {
    enrichment_true_essential[dof] = 1.0;
  }
  Vector standard_local_essential, enrichment_local_essential;
  standard_local_essential.SetSize(standard_local_size);
  enrichment_local_essential.SetSize(enrichment_local_size);
  standard_prolongation->Mult(standard_true_essential, standard_local_essential);
  enrichment_prolongation->Mult(enrichment_true_essential, enrichment_local_essential);

  Vector standard_local_multiplicity(standard_local_size);
  Vector enrichment_local_multiplicity(enrichment_local_size);
  standard_local_multiplicity = 0.0;
  enrichment_local_multiplicity = 0.0;
  auto *standard_multiplicity = standard_local_multiplicity.HostReadWrite();
  auto *enrichment_multiplicity = enrichment_local_multiplicity.HostReadWrite();
  const auto *standard_essential = standard_local_essential.HostRead();
  const auto *enrichment_essential = enrichment_local_essential.HostRead();

  HYPRE_BigInt local_rank_deficient_patches = 0;
  HYPRE_BigInt local_discarded_modes = 0;
  HYPRE_BigInt local_materially_indefinite_patches = 0;
  double local_minimum_relative_eigenvalue = std::numeric_limits<double>::infinity();
  double local_maximum_relative_resolution = 0.0;
  patches.reserve(element_matrices.size());
  for (const auto &element : element_matrices)
  {
    const int standard_size = element.standard_dofs.Size();
    const int enrichment_size = element.enrichment_dofs.Size();
    const int patch_size = standard_size + enrichment_size;
    if (element.element < 0 || patch_size <= 0 || element.mass.Height() != patch_size ||
        element.mass.Width() != patch_size || element.curl_curl.Height() != patch_size ||
        element.curl_curl.Width() != patch_size ||
        element.mass_estimated_absolute_error.Height() != patch_size ||
        element.mass_estimated_absolute_error.Width() != patch_size ||
        element.curl_curl_estimated_absolute_error.Height() != patch_size ||
        element.curl_curl_estimated_absolute_error.Width() != patch_size)
    {
      throw std::invalid_argument(
          "Coupled singular element patch has inconsistent local matrices!");
    }

    Patch patch;
    patch.local_dofs.SetSize(patch_size);
    patch.signs.SetSize(patch_size);
    std::vector<char> essential(static_cast<std::size_t>(patch_size), false);
    for (int i = 0; i < standard_size; i++)
    {
      const int local = DecodeSignedDof(element.standard_dofs[i]);
      if (local < 0 || local >= standard_local_size)
      {
        throw std::invalid_argument(
            "Coupled singular element patch has an invalid standard local DOF!");
      }
      patch.local_dofs[i] = local;
      patch.signs[i] = DecodeSignedDofSign(element.standard_dofs[i]);
      essential[static_cast<std::size_t>(i)] =
          std::abs(standard_essential[local]) > 1.0e-12;
      standard_multiplicity[local] += 1.0;
    }
    for (int i = 0; i < enrichment_size; i++)
    {
      const int local = element.enrichment_dofs[i];
      if (local < 0 || local >= enrichment_local_size)
      {
        throw std::invalid_argument(
            "Coupled singular element patch has an invalid enrichment local DOF!");
      }
      patch.local_dofs[standard_size + i] = standard_local_size + local;
      patch.signs[standard_size + i] = 1.0;
      essential[static_cast<std::size_t>(standard_size + i)] =
          std::abs(enrichment_essential[local]) > 1.0e-12;
      enrichment_multiplicity[local] += 1.0;
    }

    mfem::DenseMatrix shifted(patch_size);
    mfem::DenseMatrix shifted_error(patch_size);
    for (int i = 0; i < patch_size; i++)
    {
      for (int j = 0; j < patch_size; j++)
      {
        shifted(i, j) = stiffness_coefficient * element.curl_curl(i, j) +
                        mass_coefficient * element.mass(i, j);
        shifted_error(i, j) =
            stiffness_coefficient * element.curl_curl_estimated_absolute_error(i, j) +
            mass_coefficient * element.mass_estimated_absolute_error(i, j);
      }
    }
    for (int i = 0; i < patch_size; i++)
    {
      for (int j = i + 1; j < patch_size; j++)
      {
        const double symmetric = 0.5 * (shifted(i, j) + shifted(j, i));
        shifted(i, j) = shifted(j, i) = symmetric;
      }
      if (essential[static_cast<std::size_t>(i)])
      {
        for (int j = 0; j < patch_size; j++)
        {
          shifted(i, j) = shifted(j, i) = 0.0;
          shifted_error(i, j) = shifted_error(j, i) = 0.0;
        }
        shifted(i, i) = 1.0;
      }
    }
    auto spectrum = BuildPatchPseudoinverse(shifted, shifted_error);
    if (spectrum.discarded_modes > 0)
    {
      local_rank_deficient_patches++;
      local_discarded_modes += spectrum.discarded_modes;
    }
    if (spectrum.materially_indefinite)
    {
      local_materially_indefinite_patches++;
    }
    const double spectral_scale =
        std::max(std::numeric_limits<double>::min(), spectrum.maximum_eigenvalue);
    local_minimum_relative_eigenvalue = std::min(
        local_minimum_relative_eigenvalue, spectrum.minimum_eigenvalue / spectral_scale);
    local_maximum_relative_resolution =
        std::max(local_maximum_relative_resolution, spectrum.resolution / spectral_scale);
    patch.inverse = std::move(spectrum.inverse);
    patches.push_back(std::move(patch));
  }

  HYPRE_BigInt global_patch_count = static_cast<HYPRE_BigInt>(patches.size());
  HYPRE_BigInt global_rank_deficient_patches = local_rank_deficient_patches;
  HYPRE_BigInt global_discarded_modes = local_discarded_modes;
  HYPRE_BigInt global_materially_indefinite_patches = local_materially_indefinite_patches;
  double global_minimum_relative_eigenvalue = local_minimum_relative_eigenvalue;
  double global_maximum_relative_resolution = local_maximum_relative_resolution;
  Mpi::GlobalSum(1, &global_patch_count, standard_fespace.GetComm());
  Mpi::GlobalSum(1, &global_rank_deficient_patches, standard_fespace.GetComm());
  Mpi::GlobalSum(1, &global_discarded_modes, standard_fespace.GetComm());
  Mpi::GlobalSum(1, &global_materially_indefinite_patches, standard_fespace.GetComm());
  Mpi::GlobalMin(1, &global_minimum_relative_eigenvalue, standard_fespace.GetComm());
  Mpi::GlobalMax(1, &global_maximum_relative_resolution, standard_fespace.GetComm());
  Mpi::Print(" Singular coupled element-patch spectra: {:d}/{:d} rank-deficient patches, "
             "{:d} discarded modes, min. relative eigenvalue {:.3e}, max. relative "
             "resolution {:.3e}\n",
             global_rank_deficient_patches, global_patch_count, global_discarded_modes,
             global_minimum_relative_eigenvalue, global_maximum_relative_resolution);
  if (global_materially_indefinite_patches > 0)
  {
    throw std::runtime_error(fmt::format(
        "Found {} materially indefinite coupled singular element patches; the most "
        "negative eigenvalue exceeds the certified quadrature and roundoff resolution!",
        global_materially_indefinite_patches));
  }

  Vector standard_true_multiplicity, enrichment_true_multiplicity;
  standard_true_multiplicity.SetSize(standard_true_size);
  enrichment_true_multiplicity.SetSize(static_cast<int>(enrichment_numbering.owned_size));
  standard_prolongation->AbsMultTranspose(standard_local_multiplicity,
                                          standard_true_multiplicity);
  enrichment_prolongation->AbsMultTranspose(enrichment_local_multiplicity,
                                            enrichment_true_multiplicity);
  true_weight.SetSize(height);
  auto *weight = true_weight.HostWrite();
  const auto *standard_count = standard_true_multiplicity.HostRead();
  const auto *enrichment_count = enrichment_true_multiplicity.HostRead();
  for (int i = 0; i < standard_true_size; i++)
  {
    weight[i] = standard_count[i] > 0.0 ? 1.0 / std::sqrt(standard_count[i]) : 0.0;
  }
  for (int i = 0; i < enrichment_true_multiplicity.Size(); i++)
  {
    weight[standard_true_size + i] =
        enrichment_count[i] > 0.0 ? 1.0 / std::sqrt(enrichment_count[i]) : 0.0;
  }
  for (int dof : standard_essential_true_dofs)
  {
    weight[dof] = 0.0;
  }
  for (int dof : enrichment_essential_true_dofs)
  {
    weight[standard_true_size + dof] = 0.0;
  }
}

void ParallelElementPatchInverse::Mult(const Vector &input, Vector &output) const
{
  MFEM_VERIFY(input.Size() == width,
              "Coupled singular element patches received an invalid input vector!");
  ApplyTrueWeight(input, scaled_input);
  Vector scaled_standard, scaled_enrichment;
  scaled_standard.MakeRef(scaled_input, 0, standard_true_size);
  scaled_enrichment.MakeRef(scaled_input, standard_true_size, height - standard_true_size);
  standard_local_rhs.SetSize(standard_local_size);
  enrichment_local_rhs.SetSize(enrichment_local_size);
  standard_prolongation->Mult(scaled_standard, standard_local_rhs);
  enrichment_prolongation->Mult(scaled_enrichment, enrichment_local_rhs);

  standard_local_correction.SetSize(standard_local_size);
  enrichment_local_correction.SetSize(enrichment_local_size);
  standard_local_correction = 0.0;
  enrichment_local_correction = 0.0;
  const auto *standard_rhs = standard_local_rhs.HostRead();
  const auto *enrichment_rhs = enrichment_local_rhs.HostRead();
  auto *standard_correction = standard_local_correction.HostReadWrite();
  auto *enrichment_correction = enrichment_local_correction.HostReadWrite();
  for (auto &patch : patches)
  {
    patch.rhs.SetSize(patch.local_dofs.Size());
    auto *patch_rhs = patch.rhs.HostWrite();
    const auto *sign = patch.signs.HostRead();
    for (int i = 0; i < patch.local_dofs.Size(); i++)
    {
      const int local = patch.local_dofs[i];
      patch_rhs[i] = sign[i] * (local < standard_local_size
                                    ? standard_rhs[local]
                                    : enrichment_rhs[local - standard_local_size]);
    }
    patch.correction.SetSize(patch.local_dofs.Size());
    patch.inverse.Mult(patch.rhs, patch.correction);
    const auto *patch_correction = patch.correction.HostRead();
    for (int i = 0; i < patch.local_dofs.Size(); i++)
    {
      const int local = patch.local_dofs[i];
      if (local < standard_local_size)
      {
        standard_correction[local] += sign[i] * patch_correction[i];
      }
      else
      {
        enrichment_correction[local - standard_local_size] += sign[i] * patch_correction[i];
      }
    }
  }

  output.SetSize(height);
  Vector output_standard, output_enrichment;
  output_standard.MakeRef(output, 0, standard_true_size);
  output_enrichment.MakeRef(output, standard_true_size, height - standard_true_size);
  standard_prolongation->MultTranspose(standard_local_correction, output_standard);
  enrichment_prolongation->MultTranspose(enrichment_local_correction, output_enrichment);
  auto *output_data = output.HostReadWrite();
  const auto *weight = true_weight.HostRead();
  for (int i = 0; i < height; i++)
  {
    output_data[i] *= weight[i];
  }
}

namespace
{

ParallelSparseOperatorBlocks
CopyParallelSparseOperatorBlocks(const ParallelSparseOperatorBlocks &source)
{
  ParallelSparseOperatorBlocks copy;
  const auto copy_matrix = [](const std::unique_ptr<mfem::HypreParMatrix> &matrix)
  { return matrix ? std::make_unique<mfem::HypreParMatrix>(*matrix) : nullptr; };
  copy.enrichment_enrichment = copy_matrix(source.enrichment_enrichment);
  copy.standard_enrichment = copy_matrix(source.standard_enrichment);
  copy.enrichment_standard = copy_matrix(source.enrichment_standard);
  copy.transformed_enrichment_diagonal =
      source.transformed_enrichment_diagonal
          ? std::make_unique<Vector>(*source.transformed_enrichment_diagonal)
          : nullptr;
  return copy;
}

}  // namespace

ParallelHybridEnrichedOperator::ParallelHybridEnrichedOperator(
    std::unique_ptr<Operator> &&standard_standard,
    const ParallelSparseOperatorBlocks &enrichment,
    const mfem::Array<int> &standard_essential_true_dofs,
    const mfem::Array<int> &enrichment_essential_true_dofs,
    Operator::DiagonalPolicy diagonal_policy,
    std::shared_ptr<const Operator> coupled_patch_inverse)
  : ParallelHybridEnrichedOperator(
        std::move(standard_standard), CopyParallelSparseOperatorBlocks(enrichment),
        standard_essential_true_dofs, enrichment_essential_true_dofs, diagonal_policy,
        std::move(coupled_patch_inverse))
{
}

ParallelHybridEnrichedOperator::ParallelHybridEnrichedOperator(
    std::unique_ptr<Operator> &&standard_standard,
    ParallelSparseOperatorBlocks &&enrichment,
    const mfem::Array<int> &standard_essential_true_dofs,
    const mfem::Array<int> &enrichment_essential_true_dofs,
    Operator::DiagonalPolicy diagonal_policy,
    std::shared_ptr<const Operator> coupled_patch_inverse)
  : Operator(0), standard_standard(std::move(standard_standard)),
    enrichment(std::move(enrichment)),
    coupled_patch_inverse(std::move(coupled_patch_inverse)),
    standard_size(this->standard_standard ? this->standard_standard->Height() : 0)
{
  if (!this->standard_standard || !this->enrichment.standard_enrichment ||
      !this->enrichment.enrichment_standard || !this->enrichment.enrichment_enrichment)
  {
    throw std::invalid_argument(
        "A hybrid singular operator requires all four operator blocks!");
  }
  const auto &se = *this->enrichment.standard_enrichment;
  const auto &es = *this->enrichment.enrichment_standard;
  const auto &ee = *this->enrichment.enrichment_enrichment;
  const MPI_Comm comm = ee.GetComm();
  bool valid =
      diagonal_policy == Operator::DIAG_ONE || diagonal_policy == Operator::DIAG_ZERO;
  valid = valid && this->standard_standard->Height() == this->standard_standard->Width() &&
          standard_size == se.Height() && standard_size == es.Width() &&
          ee.Height() == ee.Width() && ee.Height() == se.Width() &&
          ee.Height() == es.Height() && SameCommunicator(comm, se.GetComm()) &&
          SameCommunicator(comm, es.GetComm()) && SameRowColumnPartition(ee) &&
          SameColumnPartition(ee, se) && SameRowPartition(ee, es) &&
          se.GetGlobalNumRows() == es.GetGlobalNumCols() &&
          ValidEssentialDofs(standard_essential_true_dofs, standard_size) &&
          ValidEssentialDofs(enrichment_essential_true_dofs, ee.Height());
  valid = valid && (!this->coupled_patch_inverse ||
                    (this->coupled_patch_inverse->Height() == standard_size + ee.Height() &&
                     this->coupled_patch_inverse->Width() == standard_size + ee.Height()));
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "A hybrid singular operator received inconsistent blocks or essential DOFs!");
  }

  std::unique_ptr<mfem::HypreParMatrix> discarded;
  this->enrichment.standard_enrichment->EliminateRows(standard_essential_true_dofs);
  discarded.reset(
      this->enrichment.enrichment_standard->EliminateCols(standard_essential_true_dofs));
  this->enrichment.enrichment_standard->EliminateRows(enrichment_essential_true_dofs);
  discarded.reset(
      this->enrichment.standard_enrichment->EliminateCols(enrichment_essential_true_dofs));
  this->enrichment.enrichment_enrichment->EliminateBC(enrichment_essential_true_dofs,
                                                      diagonal_policy);
  if (this->enrichment.transformed_enrichment_diagonal)
  {
    MFEM_VERIFY(this->enrichment.transformed_enrichment_diagonal->Size() == ee.Height(),
                "A transformed singular diagonal has inconsistent dimensions!");
    linalg::SetSubVector(*this->enrichment.transformed_enrichment_diagonal,
                         enrichment_essential_true_dofs,
                         diagonal_policy == Operator::DIAG_ONE ? 1.0 : 0.0);
  }
  RemoveExplicitZeros(*this->enrichment.standard_enrichment);
  RemoveExplicitZeros(*this->enrichment.enrichment_standard);
  RemoveExplicitZeros(*this->enrichment.enrichment_enrichment);

  height = width = standard_size + ee.Height();
}

void ParallelHybridEnrichedOperator::MakeInputBlocks(const Vector &input, Vector &standard,
                                                     Vector &enriched) const
{
  MFEM_ASSERT(input.Size() == width,
              "Incompatible input dimensions for a hybrid singular operator!");
  auto &mutable_input = const_cast<Vector &>(input);
  standard.MakeRef(mutable_input, 0, standard_size);
  enriched.MakeRef(mutable_input, standard_size, width - standard_size);
}

void ParallelHybridEnrichedOperator::MakeOutputBlocks(Vector &output, Vector &standard,
                                                      Vector &enriched) const
{
  MFEM_ASSERT(output.Size() == height,
              "Incompatible output dimensions for a hybrid singular operator!");
  standard.MakeRef(output, 0, standard_size);
  enriched.MakeRef(output, standard_size, height - standard_size);
}

void ParallelHybridEnrichedOperator::AssembleDiagonal(Vector &diagonal) const
{
  diagonal.SetSize(height);
  diagonal.UseDevice(true);
  Vector standard_diagonal(standard_size);
  Vector enrichment_diagonal(height - standard_size);
  standard_standard->AssembleDiagonal(standard_diagonal);
  enrichment.enrichment_enrichment->AssembleDiagonal(enrichment_diagonal);
  linalg::SetSubVector(diagonal, 0, standard_diagonal);
  linalg::SetSubVector(diagonal, standard_size, enrichment_diagonal);
}

void ParallelHybridEnrichedOperator::Mult(const Vector &input, Vector &output) const
{
  Vector input_standard, input_enriched, output_standard, output_enriched;
  MakeInputBlocks(input, input_standard, input_enriched);
  MakeOutputBlocks(output, output_standard, output_enriched);
  standard_standard->Mult(input_standard, output_standard);
  enrichment.standard_enrichment->AddMult(input_enriched, output_standard);
  enrichment.enrichment_standard->Mult(input_standard, output_enriched);
  enrichment.enrichment_enrichment->AddMult(input_enriched, output_enriched);
}

void ParallelHybridEnrichedOperator::MultTranspose(const Vector &input,
                                                   Vector &output) const
{
  Vector input_standard, input_enriched, output_standard, output_enriched;
  MakeInputBlocks(input, input_standard, input_enriched);
  MakeOutputBlocks(output, output_standard, output_enriched);
  standard_standard->MultTranspose(input_standard, output_standard);
  enrichment.enrichment_standard->AddMultTranspose(input_enriched, output_standard);
  enrichment.standard_enrichment->MultTranspose(input_standard, output_enriched);
  enrichment.enrichment_enrichment->AddMultTranspose(input_enriched, output_enriched);
}

void ParallelHybridEnrichedOperator::AddMult(const Vector &input, Vector &output,
                                             double coefficient) const
{
  Vector input_standard, input_enriched, output_standard, output_enriched;
  MakeInputBlocks(input, input_standard, input_enriched);
  MakeOutputBlocks(output, output_standard, output_enriched);
  standard_standard->AddMult(input_standard, output_standard, coefficient);
  enrichment.standard_enrichment->AddMult(input_enriched, output_standard, coefficient);
  enrichment.enrichment_standard->AddMult(input_standard, output_enriched, coefficient);
  enrichment.enrichment_enrichment->AddMult(input_enriched, output_enriched, coefficient);
}

void ParallelHybridEnrichedOperator::AddMultTranspose(const Vector &input, Vector &output,
                                                      double coefficient) const
{
  Vector input_standard, input_enriched, output_standard, output_enriched;
  MakeInputBlocks(input, input_standard, input_enriched);
  MakeOutputBlocks(output, output_standard, output_enriched);
  standard_standard->AddMultTranspose(input_standard, output_standard, coefficient);
  enrichment.enrichment_standard->AddMultTranspose(input_enriched, output_standard,
                                                   coefficient);
  enrichment.standard_enrichment->AddMultTranspose(input_standard, output_enriched,
                                                   coefficient);
  enrichment.enrichment_enrichment->AddMultTranspose(input_enriched, output_enriched,
                                                     coefficient);
}

ParallelTransformedHybridEnrichedOperator::ParallelTransformedHybridEnrichedOperator(
    std::unique_ptr<ParallelHybridEnrichedOperator> &&untransformed,
    const mfem::HypreParMatrix &coordinate_shift)
  : untransformed(std::move(untransformed)), coordinate_shift(&coordinate_shift),
    standard_size(this->untransformed ? this->untransformed->GetStandardSize() : 0)
{
  MFEM_VERIFY(this->untransformed && coordinate_shift.Height() == standard_size &&
                  coordinate_shift.Width() == this->untransformed->Height() - standard_size,
              "Stabilized singular coordinate shift has inconsistent dimensions!");
  height = width = this->untransformed->Height();
}

void ParallelTransformedHybridEnrichedOperator::ApplyCoordinateTransform(
    const Vector &input, Vector &output) const
{
  MFEM_VERIFY(input.Size() == width,
              "Stabilized singular coordinate transform received an invalid vector!");
  output.SetSize(width);
  output.UseDevice(input.UseDevice());
  output = input;
  auto &mutable_input = const_cast<Vector &>(input);
  Vector input_enrichment, output_standard;
  input_enrichment.MakeRef(mutable_input, standard_size, width - standard_size);
  output_standard.MakeRef(output, 0, standard_size);
  standard_work.SetSize(standard_size);
  coordinate_shift->Mult(input_enrichment, standard_work);
  output_standard -= standard_work;
}

void ParallelTransformedHybridEnrichedOperator::ApplyCoordinateTransformTranspose(
    const Vector &input, Vector &output) const
{
  MFEM_VERIFY(input.Size() == height,
              "Stabilized singular transpose transform received an invalid vector!");
  output.SetSize(height);
  output.UseDevice(input.UseDevice());
  output = input;
  auto &mutable_input = const_cast<Vector &>(input);
  Vector input_standard, output_enrichment;
  input_standard.MakeRef(mutable_input, 0, standard_size);
  output_enrichment.MakeRef(output, standard_size, height - standard_size);
  enrichment_work.SetSize(height - standard_size);
  coordinate_shift->MultTranspose(input_standard, enrichment_work);
  output_enrichment -= enrichment_work;
}

void ParallelTransformedHybridEnrichedOperator::Mult(const Vector &input,
                                                     Vector &output) const
{
  ApplyCoordinateTransform(input, transformed_input);
  transformed_action.SetSize(height);
  untransformed->Mult(transformed_input, transformed_action);
  ApplyCoordinateTransformTranspose(transformed_action, output);
}

void ParallelTransformedHybridEnrichedOperator::MultTranspose(const Vector &input,
                                                              Vector &output) const
{
  ApplyCoordinateTransform(input, transformed_input);
  transformed_action.SetSize(height);
  untransformed->MultTranspose(transformed_input, transformed_action);
  ApplyCoordinateTransformTranspose(transformed_action, output);
}

void ParallelTransformedHybridEnrichedOperator::AddMult(const Vector &input, Vector &output,
                                                        double coefficient) const
{
  Mult(input, transformed_action);
  output.Add(coefficient, transformed_action);
}

void ParallelTransformedHybridEnrichedOperator::AddMultTranspose(const Vector &input,
                                                                 Vector &output,
                                                                 double coefficient) const
{
  MultTranspose(input, transformed_action);
  output.Add(coefficient, transformed_action);
}

ParallelTransformedEnrichmentOperator::ParallelTransformedEnrichmentOperator(
    const ParallelTransformedHybridEnrichedOperator &transformed)
  : transformed(&transformed), standard_size(transformed.GetStandardSize())
{
  MFEM_VERIFY(standard_size > 0 && transformed.Height() >= standard_size &&
                  transformed.GetCoordinateShift().GetGlobalNumCols() > 0,
              "A transformed enrichment block requires nonempty standard and "
              "enrichment spaces!");
  height = width = transformed.Height() - standard_size;
}

void ParallelTransformedEnrichmentOperator::AssembleDiagonal(Vector &diagonal) const
{
  const auto &blocks = transformed->GetUntransformedOperator().GetEnrichmentBlocks();
  if (blocks.transformed_enrichment_diagonal)
  {
    diagonal = *blocks.transformed_enrichment_diagonal;
  }
  else
  {
    blocks.enrichment_enrichment->AssembleDiagonal(diagonal);
  }
}

void ParallelTransformedEnrichmentOperator::Mult(const Vector &input, Vector &output) const
{
  MFEM_VERIFY(input.Size() == width,
              "A transformed enrichment block received an invalid input vector!");
  combined_input.SetSize(standard_size + width);
  combined_input.UseDevice(input.UseDevice());
  combined_input = 0.0;
  linalg::SetSubVector(combined_input, standard_size, input);
  transformed->Mult(combined_input, combined_action);
  output.SetSize(height);
  output.UseDevice(combined_action.UseDevice());
  Vector enrichment_action;
  enrichment_action.MakeRef(combined_action, standard_size, height);
  output = enrichment_action;
}

void ParallelTransformedEnrichmentOperator::MultTranspose(const Vector &input,
                                                          Vector &output) const
{
  MFEM_VERIFY(input.Size() == width,
              "A transformed enrichment block received an invalid transpose input!");
  combined_input.SetSize(standard_size + width);
  combined_input.UseDevice(input.UseDevice());
  combined_input = 0.0;
  linalg::SetSubVector(combined_input, standard_size, input);
  transformed->MultTranspose(combined_input, combined_action);
  output.SetSize(height);
  output.UseDevice(combined_action.UseDevice());
  Vector enrichment_action;
  enrichment_action.MakeRef(combined_action, standard_size, height);
  output = enrichment_action;
}

HYPRE_BigInt RemoveExplicitZeros(mfem::HypreParMatrix &matrix)
{
  auto *parallel = static_cast<hypre_ParCSRMatrix *>(matrix);
  if (hypre_ParCSRMatrixMemoryLocation(parallel) != HYPRE_MEMORY_HOST)
  {
    throw std::invalid_argument(
        "Singular sparse compaction currently requires a CPU HypreParMatrix!");
  }

  HYPRE_BigInt local_removed = 0;
  const auto compact = [&local_removed](hypre_CSRMatrix *block, bool diagonal)
  {
    auto *offsets = hypre_CSRMatrixI(block);
    auto *columns = hypre_CSRMatrixJ(block);
    auto *values = hypre_CSRMatrixData(block);
    const HYPRE_Int rows = hypre_CSRMatrixNumRows(block);
    HYPRE_Int read_begin = offsets[0];
    HYPRE_Int write = 0;
    offsets[0] = 0;
    for (HYPRE_Int row = 0; row < rows; row++)
    {
      const HYPRE_Int read_end = offsets[row + 1];
      for (HYPRE_Int entry = read_begin; entry < read_end; entry++)
      {
        if (values[entry] != 0.0 || (diagonal && columns[entry] == row))
        {
          columns[write] = columns[entry];
          values[write] = values[entry];
          write++;
        }
        else
        {
          local_removed++;
        }
      }
      offsets[row + 1] = write;
      read_begin = read_end;
    }
    hypre_CSRMatrixNumNonzeros(block) = write;
  };
  compact(hypre_ParCSRMatrixDiag(parallel), true);
  compact(hypre_ParCSRMatrixOffd(parallel), false);

  if (hypre_ParCSRMatrixSetNumNonzeros(parallel) != 0 ||
      hypre_ParCSRMatrixSetDNumNonzeros(parallel) != 0)
  {
    throw std::runtime_error(
        "Hypre failed to update singular sparse-matrix nonzero counts!");
  }
  HYPRE_BigInt global_removed = local_removed;
  Mpi::GlobalSum(1, &global_removed, matrix.GetComm());
  return global_removed;
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

std::unique_ptr<mfem::HypreParMatrix>
BuildParallelEnrichedProlongation(const mfem::HypreParMatrix &standard_prolongation,
                                  const TrueDofMap &enrichment_numbering)
{
  const MPI_Comm comm = standard_prolongation.GetComm();
  bool valid = standard_prolongation.GetGlobalNumRows() > 0 &&
               standard_prolongation.GetGlobalNumCols() > 0 &&
               enrichment_numbering.global_size > 0 &&
               enrichment_numbering.owned_offset >= 0 &&
               enrichment_numbering.owned_size >= 0 &&
               enrichment_numbering.owned_offset <= enrichment_numbering.global_size &&
               enrichment_numbering.owned_size <=
                   enrichment_numbering.global_size - enrichment_numbering.owned_offset &&
               enrichment_numbering.owned_size <= std::numeric_limits<int>::max();
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Parallel singular prolongation received inconsistent standard or enrichment "
        "dimensions!");
  }

  const int local_size = static_cast<int>(enrichment_numbering.owned_size);
  std::vector<int> rows(local_size + 1);
  std::vector<HYPRE_BigInt> columns(local_size);
  std::vector<double> values(local_size, 1.0);
  for (int i = 0; i < local_size; i++)
  {
    rows[i] = i;
    columns[i] = enrichment_numbering.owned_offset + i;
  }
  rows[local_size] = local_size;
  const auto partition =
      BuildPartition(comm, enrichment_numbering.owned_offset,
                     enrichment_numbering.owned_size, enrichment_numbering.global_size);
  auto identity = std::make_unique<mfem::HypreParMatrix>(
      comm, local_size, enrichment_numbering.global_size, enrichment_numbering.global_size,
      rows.data(), columns.data(), values.data(), partition.data(), partition.data());

  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks = nullptr;
  blocks(0, 0) = &standard_prolongation;
  blocks(1, 1) = identity.get();
  return std::unique_ptr<mfem::HypreParMatrix>(mfem::HypreParMatrixFromBlocks(blocks));
}

std::unique_ptr<mfem::HypreParMatrix> BuildParallelEnrichedProlongation(
    const mfem::HypreParMatrix &standard_prolongation,
    const mfem::HypreParMatrix &standard_enrichment_correction,
    const TrueDofMap &enrichment_numbering)
{
  const MPI_Comm comm = standard_prolongation.GetComm();
  bool valid = standard_prolongation.GetGlobalNumRows() > 0 &&
               standard_prolongation.GetGlobalNumCols() > 0 &&
               standard_enrichment_correction.GetComm() == comm &&
               standard_enrichment_correction.Height() == standard_prolongation.Height() &&
               standard_enrichment_correction.GetGlobalNumRows() ==
                   standard_prolongation.GetGlobalNumRows() &&
               standard_enrichment_correction.GetGlobalNumCols() ==
                   enrichment_numbering.global_size &&
               enrichment_numbering.global_size > 0 &&
               enrichment_numbering.owned_offset >= 0 &&
               enrichment_numbering.owned_size >= 0 &&
               enrichment_numbering.owned_offset <= enrichment_numbering.global_size &&
               enrichment_numbering.owned_size <=
                   enrichment_numbering.global_size - enrichment_numbering.owned_offset &&
               enrichment_numbering.owned_size <= std::numeric_limits<int>::max();
  Mpi::GlobalAnd(1, &valid, comm);
  if (!valid)
  {
    throw std::invalid_argument(
        "Corrected parallel singular prolongation received inconsistent dimensions!");
  }

  const int local_size = static_cast<int>(enrichment_numbering.owned_size);
  std::vector<int> rows(local_size + 1);
  std::vector<HYPRE_BigInt> columns(local_size);
  std::vector<double> values(local_size, 1.0);
  for (int i = 0; i < local_size; i++)
  {
    rows[i] = i;
    columns[i] = enrichment_numbering.owned_offset + i;
  }
  rows[local_size] = local_size;
  const auto partition =
      BuildPartition(comm, enrichment_numbering.owned_offset,
                     enrichment_numbering.owned_size, enrichment_numbering.global_size);
  auto identity = std::make_unique<mfem::HypreParMatrix>(
      comm, local_size, enrichment_numbering.global_size, enrichment_numbering.global_size,
      rows.data(), columns.data(), values.data(), partition.data(), partition.data());

  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks = nullptr;
  blocks(0, 0) = &standard_prolongation;
  blocks(0, 1) = &standard_enrichment_correction;
  blocks(1, 1) = identity.get();
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
  RemoveExplicitZeros(*matrix);
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
