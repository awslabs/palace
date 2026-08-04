// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "hierarchicalerrorestimator.hpp"

#include <map>

namespace palace::fem::hierarchical
{

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

}  // namespace palace::fem::hierarchical
