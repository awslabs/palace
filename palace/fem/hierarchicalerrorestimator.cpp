// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hierarchicalerrorestimator.hpp"

#include <algorithm>
#include <map>

namespace palace::fem::hierarchical
{

namespace
{

int UnsignedDof(int dof)
{
  return dof >= 0 ? dof : -1 - dof;
}

double DofSign(int dof)
{
  return dof >= 0 ? 1.0 : -1.0;
}

}  // namespace

mfem::Vector AssembleResidual(int combined_size,
                              const std::vector<LocalOperatorContribution> &contributions,
                              const mfem::Vector &injected,
                              const std::vector<bool> &essential)
{
  MFEM_VERIFY(injected.Size() == combined_size &&
                  essential.size() == static_cast<std::size_t>(combined_size),
              "Hierarchical residual received inconsistent combined vectors!");
  mfem::Vector residual(combined_size);
  residual = 0.0;
  for (const auto &data : contributions)
  {
    MFEM_VERIFY(data.matrix.Height() == static_cast<int>(data.dofs.size()) &&
                    data.matrix.Width() == static_cast<int>(data.dofs.size()) &&
                    data.rhs.Size() == static_cast<int>(data.dofs.size()),
                "Hierarchical residual received an inconsistent local contribution!");
    mfem::Vector local(static_cast<int>(data.dofs.size()));
    for (int i = 0; i < local.Size(); i++)
    {
      MFEM_VERIFY(data.dofs[i] >= 0 && data.dofs[i] < combined_size,
                  "Hierarchical residual received an invalid combined DOF!");
      local(i) = injected(data.dofs[i]);
    }
    mfem::Vector action(local.Size());
    data.matrix.Mult(local, action);
    for (int i = 0; i < local.Size(); i++)
    {
      residual(data.dofs[i]) += data.rhs(i) - action(i);
    }
  }
  for (int dof = 0; dof < combined_size; dof++)
  {
    if (essential[dof])
    {
      residual(dof) = 0.0;
    }
  }
  return residual;
}

long long
AssembleRestrictedOperator(const std::vector<LocalOperatorContribution> &contributions,
                           const std::set<int> &support_elements,
                           const std::vector<SparseColumn> &basis,
                           const mfem::Vector &residual, mfem::DenseMatrix &restricted,
                           mfem::Vector &restricted_rhs)
{
  const int dimension = static_cast<int>(basis.size());
  restricted.SetSize(dimension);
  restricted = 0.0;
  restricted_rhs.SetSize(dimension);
  restricted_rhs = 0.0;
  std::map<int, std::vector<std::pair<int, double>>> row_basis;
  for (int column = 0; column < dimension; column++)
  {
    MFEM_VERIFY(basis[column].dofs.size() == basis[column].values.size(),
                "Hierarchical basis column has inconsistent sparse storage!");
    for (std::size_t entry = 0; entry < basis[column].dofs.size(); entry++)
    {
      const int dof = basis[column].dofs[entry];
      MFEM_VERIFY(dof >= 0 && dof < residual.Size(),
                  "Hierarchical basis column contains an invalid DOF!");
      row_basis[dof].push_back({column, basis[column].values[entry]});
      restricted_rhs(column) += basis[column].values[entry] * residual(dof);
    }
  }

  long long touched = 0;
  for (const auto &data : contributions)
  {
    if (data.support_element >= 0 &&
        support_elements.find(data.support_element) == support_elements.end())
    {
      continue;
    }
    for (int i = 0; i < static_cast<int>(data.dofs.size()); i++)
    {
      const auto row = row_basis.find(data.dofs[i]);
      if (row == row_basis.end())
      {
        continue;
      }
      for (int j = 0; j < static_cast<int>(data.dofs.size()); j++)
      {
        const auto column = row_basis.find(data.dofs[j]);
        if (column == row_basis.end())
        {
          continue;
        }
        touched++;
        for (const auto &[row_basis_index, row_value] : row->second)
        {
          for (const auto &[column_basis_index, column_value] : column->second)
          {
            restricted(row_basis_index, column_basis_index) +=
                row_value * data.matrix(i, j) * column_value;
          }
        }
      }
    }
  }
  return touched;
}

double Energy(const std::vector<LocalOperatorContribution> &contributions,
              const mfem::Vector &vector)
{
  double energy = 0.0;
  for (const auto &data : contributions)
  {
    mfem::Vector local(static_cast<int>(data.dofs.size()));
    for (int i = 0; i < local.Size(); i++)
    {
      MFEM_VERIFY(data.dofs[i] >= 0 && data.dofs[i] < vector.Size(),
                  "Hierarchical energy received an invalid combined DOF!");
      local(i) = vector(data.dofs[i]);
    }
    mfem::Vector action(local.Size());
    data.matrix.Mult(local, action);
    energy += mfem::InnerProduct(local, action);
  }
  return energy;
}

SparseInjection BuildSparsePInjection(mfem::Mesh &mesh,
                                      mfem::FiniteElementSpace &coarse_space,
                                      mfem::FiniteElementSpace &fine_space)
{
  SparseInjection injection;
  std::vector<std::map<int, double>> entries(coarse_space.GetVSize());
  mfem::DenseMatrix local_transfer;
  mfem::Array<int> coarse_dofs, fine_dofs;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mfem::IsoparametricTransformation identity;
    identity.SetIdentityTransformation(mesh.GetElementGeometry(element));
    mfem::DofTransformation coarse_transformation, fine_transformation;
    coarse_space.GetElementVDofs(element, coarse_dofs, coarse_transformation);
    fine_space.GetElementVDofs(element, fine_dofs, fine_transformation);
    injection.signed_coarse_dofs += static_cast<int>(std::count_if(
        coarse_dofs.begin(), coarse_dofs.end(), [](int dof) { return dof < 0; }));
    injection.signed_fine_dofs += static_cast<int>(
        std::count_if(fine_dofs.begin(), fine_dofs.end(), [](int dof) { return dof < 0; }));
    injection.nonidentity_transformations += !coarse_transformation.IsIdentity();
    injection.nonidentity_transformations += !fine_transformation.IsIdentity();
    fine_space.GetFE(element)->GetTransferMatrix(*coarse_space.GetFE(element), identity,
                                                 local_transfer);
    std::set<int> element_coarse;
    for (int dof : coarse_dofs)
    {
      element_coarse.insert(UnsignedDof(dof));
    }
    for (int global_coarse : element_coarse)
    {
      mfem::Vector local_coarse(coarse_dofs.Size());
      local_coarse = 0.0;
      for (int i = 0; i < coarse_dofs.Size(); i++)
      {
        if (UnsignedDof(coarse_dofs[i]) == global_coarse)
        {
          local_coarse(i) = DofSign(coarse_dofs[i]);
        }
      }
      coarse_transformation.InvTransformPrimal(local_coarse);
      mfem::Vector local_fine(fine_dofs.Size());
      local_transfer.Mult(local_coarse, local_fine);
      fine_transformation.TransformPrimal(local_fine);
      auto &column = entries[global_coarse];
      for (int i = 0; i < fine_dofs.Size(); i++)
      {
        const double value = DofSign(fine_dofs[i]) * local_fine(i);
        if (std::abs(value) <= 1.0e-13)
        {
          continue;
        }
        const int row = UnsignedDof(fine_dofs[i]);
        const auto [found, inserted] = column.emplace(row, value);
        if (!inserted)
        {
          injection.consistency_error =
              std::max(injection.consistency_error, std::abs(found->second - value));
        }
      }
    }
  }
  injection.columns.resize(entries.size());
  for (std::size_t column = 0; column < entries.size(); column++)
  {
    for (const auto &[row, value] : entries[column])
    {
      injection.columns[column].dofs.push_back(row);
      injection.columns[column].values.push_back(value);
    }
  }
  return injection;
}

PatchLiftingReport EstimateByPatchLifting(
    mfem::Mesh &mesh, mfem::FiniteElementSpace &coarse_space,
    mfem::FiniteElementSpace &fine_space, const SparseInjection &injection,
    const std::vector<LocalOperatorContribution> &residual_contributions,
    const std::vector<LocalOperatorContribution> &metric_contributions,
    const std::vector<bool> &fine_essential, const std::vector<bool> &coarse_essential,
    const mfem::Vector &coarse_combined,
    const std::vector<std::vector<int>> &element_enrichment_guests,
    const PatchLiftingOptions &options)
{
  PatchLiftingReport report;
  const int coarse_standard = coarse_space.GetVSize();
  const int fine_standard = fine_space.GetVSize();
  const int enrichment = static_cast<int>(fine_essential.size()) - fine_standard;
  const int combined_size = fine_standard + enrichment;
  MFEM_VERIFY(enrichment >= 0 && coarse_combined.Size() == coarse_standard + enrichment &&
                  coarse_essential.size() ==
                      static_cast<std::size_t>(coarse_standard + enrichment) &&
                  static_cast<int>(injection.columns.size()) == coarse_standard &&
                  static_cast<int>(element_enrichment_guests.size()) == mesh.GetNE(),
              "Hierarchical patch lifting received inconsistent combined layouts!");

  // Inject the coarse combined solution into the fine combined layout and assemble the
  // uneliminated residual there.
  mfem::Vector injected(combined_size);
  injected = 0.0;
  for (int column = 0; column < coarse_standard; column++)
  {
    for (std::size_t entry = 0; entry < injection.columns[column].dofs.size(); entry++)
    {
      injected(injection.columns[column].dofs[entry]) +=
          coarse_combined(column) * injection.columns[column].values[entry];
    }
  }
  for (int dof = 0; dof < enrichment; dof++)
  {
    injected(fine_standard + dof) = coarse_combined(coarse_standard + dof);
  }
  mfem::Vector residual =
      AssembleResidual(combined_size, residual_contributions, injected, fine_essential);

  std::vector<std::set<int>> dof_elements(combined_size);
  for (const auto &data : metric_contributions)
  {
    for (int dof : data.dofs)
    {
      dof_elements[dof].insert(data.support_element);
    }
  }

  // Conforming complement of the coarse entity trace inside the fine entity trace. The
  // result spans the fine directions on this entity that the injected coarse space cannot
  // represent; each such direction is owned by exactly one patch.
  const auto entity_complement =
      [&](mfem::Array<int> fine_entity, mfem::Array<int> coarse_entity)
  {
    mfem::FiniteElementSpace::AdjustVDofs(fine_entity);
    mfem::FiniteElementSpace::AdjustVDofs(coarse_entity);
    std::vector<int> rows, columns;
    for (int dof : fine_entity)
    {
      if (!fine_essential[dof])
      {
        rows.push_back(dof);
      }
    }
    for (int dof : coarse_entity)
    {
      if (!coarse_essential[dof])
      {
        columns.push_back(dof);
      }
    }
    std::vector<mfem::Vector> range;
    for (int column : columns)
    {
      mfem::Vector vector(static_cast<int>(rows.size()));
      vector = 0.0;
      for (int i = 0; i < static_cast<int>(rows.size()); i++)
      {
        for (std::size_t entry = 0; entry < injection.columns[column].dofs.size(); entry++)
        {
          if (injection.columns[column].dofs[entry] == rows[i])
          {
            vector(i) = injection.columns[column].values[entry];
          }
        }
      }
      for (int repeat = 0; repeat < 2; repeat++)
      {
        for (const auto &basis : range)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
      }
      const double norm = vector.Norml2();
      if (norm > 1.0e-12)
      {
        vector /= norm;
        range.push_back(vector);
      }
    }
    std::vector<mfem::Vector> complement;
    for (int direction = 0; direction < static_cast<int>(rows.size()); direction++)
    {
      mfem::Vector vector(static_cast<int>(rows.size()));
      vector = 0.0;
      vector(direction) = 1.0;
      for (int repeat = 0; repeat < 2; repeat++)
      {
        for (const auto &basis : range)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
        for (const auto &basis : complement)
        {
          vector.Add(-mfem::InnerProduct(vector, basis), basis);
        }
      }
      const double norm = vector.Norml2();
      if (norm > 1.0e-10)
      {
        vector /= norm;
        complement.push_back(vector);
      }
    }
    std::vector<SparseColumn> result;
    for (const auto &local : complement)
    {
      SparseColumn column;
      for (int i = 0; i < static_cast<int>(rows.size()); i++)
      {
        if (std::abs(local(i)) > options.drop_tolerance)
        {
          column.dofs.push_back(rows[i]);
          column.values.push_back(local(i));
        }
      }
      result.push_back(std::move(column));
    }
    return result;
  };

  struct Patch
  {
    int owned = 0;
    std::vector<SparseColumn> basis;
    std::vector<int> coarse_guests;
    std::vector<int> enrichment_guests;
    mfem::Vector coefficients;
    mfem::DenseMatrix restricted;
    std::set<int> support;
  };
  std::vector<Patch> patches;
  std::vector<int> overlap(mesh.GetNE(), 0);
  mfem::Array<int> element_dofs;
  enum class PatchEntity
  {
    EDGE,
    FACE,
    INTERIOR
  };
  const auto add_patch = [&](std::vector<SparseColumn> owned,
                             const std::vector<int> &owner_elements, PatchEntity entity)
  {
    if (owned.empty())
    {
      return;
    }
    Patch patch;
    patch.owned = static_cast<int>(owned.size());
    patch.basis = std::move(owned);
    report.owned_modes += patch.owned;
    report.edge_patches += entity == PatchEntity::EDGE;
    report.face_patches += entity == PatchEntity::FACE;
    report.interior_patches += entity == PatchEntity::INTERIOR;
    std::set<int> coarse_guests, enrichment_guests;
    for (int element : owner_elements)
    {
      mfem::DofTransformation transformation;
      coarse_space.GetElementVDofs(element, element_dofs, transformation);
      for (int dof : element_dofs)
      {
        dof = UnsignedDof(dof);
        if (!coarse_essential[dof])
        {
          coarse_guests.insert(dof);
        }
      }
      for (int dof : element_enrichment_guests[element])
      {
        MFEM_VERIFY(dof >= 0 && dof < enrichment,
                    "Hierarchical patch lifting received an invalid enrichment guest!");
        if (!fine_essential[fine_standard + dof])
        {
          enrichment_guests.insert(dof);
        }
      }
    }
    for (int dof : coarse_guests)
    {
      patch.coarse_guests.push_back(dof);
      patch.basis.push_back(injection.columns[dof]);
    }
    for (int dof : enrichment_guests)
    {
      SparseColumn column;
      column.dofs.push_back(fine_standard + dof);
      column.values.push_back(1.0);
      patch.enrichment_guests.push_back(dof);
      patch.basis.push_back(column);
    }
    const int dimension = static_cast<int>(patch.basis.size());
    report.maximum_patch_dimension = std::max(report.maximum_patch_dimension, dimension);
    bool touches_essential = false;
    for (const auto &column : patch.basis)
    {
      for (int dof : column.dofs)
      {
        touches_essential = touches_essential || fine_essential[dof];
        patch.support.insert(dof_elements[dof].begin(), dof_elements[dof].end());
      }
    }
    MFEM_VERIFY(!touches_essential,
                "Hierarchical patch basis columns must avoid essential fine equations!");
    report.maximum_support_elements =
        std::max(report.maximum_support_elements, static_cast<int>(patch.support.size()));
    for (int element : patch.support)
    {
      overlap[element]++;
    }
    mfem::DenseMatrix restricted;
    mfem::Vector patch_rhs;
    AssembleRestrictedOperator(metric_contributions, patch.support, patch.basis, residual,
                               restricted, patch_rhs);
    patch.restricted = restricted;
    mfem::DenseMatrixInverse inverse(restricted, true);
    mfem::DenseMatrix inverse_matrix;
    inverse.GetInverseMatrix(inverse_matrix);
    report.maximum_patch_condition = std::max(report.maximum_patch_condition,
                                              restricted.FNorm() * inverse_matrix.FNorm());
    patch.coefficients.SetSize(dimension);
    inverse.Mult(patch_rhs, patch.coefficients);
    mfem::Vector solve_residual(dimension);
    restricted.Mult(patch.coefficients, solve_residual);
    solve_residual -= patch_rhs;
    report.maximum_patch_residual =
        std::max(report.maximum_patch_residual,
                 solve_residual.Norml2() / std::max(patch_rhs.Norml2(), 1.0e-30));
    patches.push_back(std::move(patch));
  };

  std::vector<std::vector<int>> edge_elements(mesh.GetNEdges());
  mfem::Array<int> edges, orientations;
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mesh.GetElementEdges(element, edges, orientations);
    for (int edge : edges)
    {
      edge_elements[edge].push_back(element);
    }
  }
  for (int edge = 0; edge < mesh.GetNEdges(); edge++)
  {
    mfem::Array<int> coarse_entity, fine_entity;
    coarse_space.GetEdgeDofs(edge, coarse_entity);
    fine_space.GetEdgeDofs(edge, fine_entity);
    add_patch(entity_complement(fine_entity, coarse_entity), edge_elements[edge],
              PatchEntity::EDGE);
  }
  if (mesh.Dimension() == 3)
  {
    for (int face = 0; face < mesh.GetNFaces(); face++)
    {
      mfem::Array<int> coarse_entity, fine_entity;
      coarse_space.GetFaceDofs(face, coarse_entity);
      fine_space.GetFaceDofs(face, fine_entity);
      int first, second;
      mesh.GetFaceElements(face, &first, &second);
      std::vector<int> owner_elements{first};
      if (second >= 0)
      {
        owner_elements.push_back(second);
      }
      add_patch(entity_complement(fine_entity, coarse_entity), owner_elements,
                PatchEntity::FACE);
    }
  }
  for (int element = 0; element < mesh.GetNE(); element++)
  {
    mfem::Array<int> coarse_entity, fine_entity;
    coarse_space.GetElementInteriorDofs(element, coarse_entity);
    fine_space.GetElementInteriorDofs(element, fine_entity);
    add_patch(entity_complement(fine_entity, coarse_entity), {element},
              PatchEntity::INTERIOR);
  }
  int coarse_free = 0, fine_free = 0;
  for (int dof = 0; dof < coarse_standard; dof++)
  {
    coarse_free += !coarse_essential[dof];
  }
  for (int dof = 0; dof < fine_standard; dof++)
  {
    fine_free += !fine_essential[dof];
  }
  MFEM_VERIFY(report.owned_modes == fine_free - coarse_free,
              "Hierarchical patch ownership must cover every conforming complement "
              "direction exactly once!");
  report.maximum_element_overlap =
      overlap.empty() ? 0 : *std::max_element(overlap.begin(), overlap.end());

  // Partition-of-unity reconstruction: owned modes are inserted once; coarse and
  // enrichment guests are averaged by stable DOF identity across their hosting patches.
  report.coarse_guest_counts.assign(coarse_standard, 0);
  report.enrichment_guest_counts.assign(enrichment, 0);
  mfem::Vector raw(combined_size);
  raw = 0.0;
  std::vector<double> coarse_sum(coarse_standard, 0.0), enrichment_sum(enrichment, 0.0);
  std::vector<int> coarse_count(coarse_standard, 0), enrichment_count(enrichment, 0);
  const auto reconstruct = [&](mfem::Vector &target)
  {
    target = 0.0;
    std::fill(coarse_sum.begin(), coarse_sum.end(), 0.0);
    std::fill(enrichment_sum.begin(), enrichment_sum.end(), 0.0);
    std::fill(coarse_count.begin(), coarse_count.end(), 0);
    std::fill(enrichment_count.begin(), enrichment_count.end(), 0);
    for (const auto &patch : patches)
    {
      for (int basis = 0; basis < patch.owned; basis++)
      {
        for (std::size_t entry = 0; entry < patch.basis[basis].dofs.size(); entry++)
        {
          target(patch.basis[basis].dofs[entry]) +=
              patch.coefficients(basis) * patch.basis[basis].values[entry];
        }
      }
      int coefficient = patch.owned;
      for (int dof : patch.coarse_guests)
      {
        coarse_sum[dof] += patch.coefficients(coefficient++);
        coarse_count[dof]++;
      }
      for (int dof : patch.enrichment_guests)
      {
        enrichment_sum[dof] += patch.coefficients(coefficient++);
        enrichment_count[dof]++;
      }
    }
    for (int dof = 0; dof < coarse_standard; dof++)
    {
      if (coarse_count[dof] > 0)
      {
        const double coefficient = coarse_sum[dof] / coarse_count[dof];
        for (std::size_t entry = 0; entry < injection.columns[dof].dofs.size(); entry++)
        {
          target(injection.columns[dof].dofs[entry]) +=
              coefficient * injection.columns[dof].values[entry];
        }
      }
    }
    for (int dof = 0; dof < enrichment; dof++)
    {
      if (enrichment_count[dof] > 0)
      {
        target(fine_standard + dof) += enrichment_sum[dof] / enrichment_count[dof];
      }
    }
  };
  reconstruct(raw);
  report.coarse_guest_counts = coarse_count;
  report.enrichment_guest_counts = enrichment_count;
  double raw_energy = Energy(metric_contributions, raw);
  double raw_work = mfem::InnerProduct(raw, residual);
  const double alpha = raw_energy > 0.0 ? raw_work / raw_energy : 0.0;
  raw *= alpha;

  // Additional additive-Schwarz defect-correction sweeps relaxing the coercive metric
  // equation whose right-hand side is the (possibly indefinite) combined residual.
  for (int sweep = 1; sweep < options.sweeps; sweep++)
  {
    mfem::Vector current_residual(residual);
    for (const auto &data : metric_contributions)
    {
      mfem::Vector local(static_cast<int>(data.dofs.size()));
      for (int i = 0; i < local.Size(); i++)
      {
        local(i) = raw(data.dofs[i]);
      }
      mfem::Vector action(local.Size());
      data.matrix.Mult(local, action);
      for (int i = 0; i < local.Size(); i++)
      {
        current_residual(data.dofs[i]) -= action(i);
      }
    }
    for (auto &patch : patches)
    {
      mfem::Vector patch_rhs(static_cast<int>(patch.basis.size()));
      for (int basis = 0; basis < static_cast<int>(patch.basis.size()); basis++)
      {
        patch_rhs(basis) = 0.0;
        for (std::size_t entry = 0; entry < patch.basis[basis].dofs.size(); entry++)
        {
          patch_rhs(basis) += patch.basis[basis].values[entry] *
                              current_residual(patch.basis[basis].dofs[entry]);
        }
      }
      mfem::DenseMatrixInverse inverse(patch.restricted, true);
      inverse.Mult(patch_rhs, patch.coefficients);
    }
    mfem::Vector delta(combined_size);
    reconstruct(delta);
    const double delta_energy = Energy(metric_contributions, delta);
    const double delta_work = mfem::InnerProduct(delta, current_residual);
    if (delta_energy > 0.0)
    {
      raw.Add(delta_work / delta_energy, delta);
    }
  }

  // A final scalar projection restores the energy/work identity without changing ranking.
  raw_energy = Energy(metric_contributions, raw);
  raw_work = mfem::InnerProduct(raw, residual);
  const double final_scale = raw_energy > 0.0 ? raw_work / raw_energy : 0.0;
  raw *= final_scale;
  report.energy = final_scale * final_scale * raw_energy;
  report.work = final_scale * raw_work;
  report.indicator.assign(mesh.GetNE(), 0.0);
  for (const auto &data : metric_contributions)
  {
    mfem::Vector local(static_cast<int>(data.dofs.size()));
    for (int i = 0; i < local.Size(); i++)
    {
      local(i) = raw(data.dofs[i]);
    }
    mfem::Vector action(local.Size());
    data.matrix.Mult(local, action);
    report.indicator[data.support_element] += mfem::InnerProduct(local, action);
  }
  return report;
}

}  // namespace palace::fem::hierarchical
