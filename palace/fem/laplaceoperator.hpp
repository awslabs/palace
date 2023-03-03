// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LAPLACE_OPERATOR_HPP
#define PALACE_LAPLACE_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/materialoperator.hpp"

namespace palace
{

class IoData;

//
// A class handling discretization of Laplace problems for electrostatics.
//
class LaplaceOperator
{
private:
  // Essential boundary condition markers.
  mfem::Array<int> dbc_marker, dbc_tdof_list;

  // Options for system matrix assembly.
  const int skip_zeros;  // Whether to skip the zeros during assembly of operators
  const bool pc_gmg;     // Whether to use geometric multigrid in preconditioning

  // Helper variable and function for log file printing.
  bool print_hdr;
  void PrintHeader();

  // Objects defining the finite element spaces for the electrostatic potential (H1) and
  // electric field (Nedelec) on the given mesh.
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_fecs;
  mfem::ND_FECollection nd_fec;
  mfem::ParFiniteElementSpaceHierarchy h1_fespaces;
  mfem::ParFiniteElementSpace nd_fespace;

  // Operator for domain material properties.
  MaterialOperator mat_op;

  // Boundary attributes for each terminal index.
  std::map<int, mfem::Array<int>> source_attr_lists;

public:
  LaplaceOperator(const IoData &iodata,
                  const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh);

  // Returns array marking Dirichlet BC attributes and local subdomain vdofs.
  const auto &GetDbcMarker() const { return dbc_marker; }
  const auto &GetDbcTDofList() const { return dbc_tdof_list; }

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access source attribute lists.
  const auto &GetSources() const { return source_attr_lists; }

  // Return the parallel finite element space objects.
  auto &GetH1Spaces() { return h1_fespaces; }
  auto &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  auto &GetNDSpace() { return nd_fespace; }

  // Return the number of true (conforming) dofs on the finest H1 space.
  auto GetNDof() { return GetH1Space().GetNConformingDofs(); }

  // Construct and return system matrix representing discretized Laplace operator for
  // Gauss's law.
  void GetStiffnessMatrix(std::vector<std::unique_ptr<mfem::Operator>> &K,
                          std::vector<std::unique_ptr<mfem::Operator>> &Ke);

  // Construct and return the discrete negative gradient matrix.
  std::unique_ptr<mfem::Operator> GetNegGradMatrix();

  // Assemble the solution boundary conditions and right-hand side vector for a nonzero
  // prescribed voltage on the specified surface index.
  void GetExcitationVector(int idx, const mfem::Operator &K, const mfem::Operator &Ke,
                           mfem::Vector &X, mfem::Vector &RHS);
};

}  // namespace palace

#endif  // PALACE_LAPLACE_OPERATOR_HPP
