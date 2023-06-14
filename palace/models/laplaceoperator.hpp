// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_LAPLACE_OPERATOR_HPP
#define PALACE_MODELS_LAPLACE_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

class IoData;

//
// A class handling discretization of Laplace problems for electrostatics.
//
class LaplaceOperator
{
private:
  const mfem::AssemblyLevel assembly_level;  // Use full or partial assembly for operators
  const int skip_zeros;                      // Skip zeros during full assembly of operators
  const bool pc_mg;                          // Use geometric multigrid in preconditioning

  // Helper variable for log file printing.
  bool print_hdr;

  // Essential boundary condition markers.
  mfem::Array<int> dbc_marker;
  std::vector<mfem::Array<int>> dbc_tdof_lists;

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

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access source attribute lists.
  const auto &GetSources() const { return source_attr_lists; }

  // Return the parallel finite element space objects.
  auto &GetH1Spaces() { return h1_fespaces; }
  auto &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  auto &GetNDSpace() { return nd_fespace; }

  // Construct and return system matrix representing discretized Laplace operator for
  // Gauss's law.
  void GetStiffnessMatrix(std::vector<std::unique_ptr<ParOperator>> &K);

  // Construct and return the discrete gradient matrix.
  std::unique_ptr<ParOperator> GetGradMatrix();

  // Assemble the solution boundary conditions and right-hand side vector for a nonzero
  // prescribed voltage on the specified surface index.
  void GetExcitationVector(int idx, const ParOperator &K, Vector &X, Vector &RHS);
};

}  // namespace palace

#endif  // PALACE_MODELS_LAPLACE_OPERATOR_HPP
