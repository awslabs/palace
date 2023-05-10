// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_CURL_CURL_OPERATOR_HPP
#define PALACE_MODELS_CURL_CURL_OPERATOR_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
#include "models/surfacecurrentoperator.hpp"

namespace palace
{

class IoData;

//
// A class handling discretization of curl-curl problems for magnetostatics.
//
class CurlCurlOperator
{
private:
  const mfem::AssemblyLevel assembly_level;  // Use full or partial assembly for operators
  const int skip_zeros;                      // Skip zeros during full assembly of operators
  const bool pc_gmg;                         // Use geometric multigrid in preconditioning

  // Helper variable for log file printing.
  bool print_hdr;

  // Essential boundary condition markers.
  mfem::Array<int> dbc_marker;
  std::vector<mfem::Array<int>> dbc_tdof_lists;
  void CheckBoundaryProperties();

  // Objects defining the finite element spaces for the magnetic vector potential
  // (Nedelec) and magnetic flux density (Raviart-Thomas) on the given mesh. The H1 spaces
  // are used for various purposes throughout the code including postprocessing.
  std::vector<std::unique_ptr<mfem::ND_FECollection>> nd_fecs;
  mfem::H1_FECollection h1_fec;
  mfem::RT_FECollection rt_fec;
  mfem::ParFiniteElementSpaceHierarchy nd_fespaces;
  mfem::ParFiniteElementSpace h1_fespace, rt_fespace;

  // Operator for domain material properties.
  MaterialOperator mat_op;

  // Operator for source current excitation.
  SurfaceCurrentOperator surf_j_op;

public:
  CurlCurlOperator(const IoData &iodata,
                   const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh);

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access to underlying BC operator objects for postprocessing.
  const auto &GetSurfaceCurrentOp() const { return surf_j_op; }

  // Return the parallel finite element space objects.
  auto &GetNDSpaces() { return nd_fespaces; }
  auto &GetNDSpace() { return nd_fespaces.GetFinestFESpace(); }
  auto &GetH1Space() { return h1_fespace; }
  auto &GetRTSpace() { return rt_fespace; }

  // Construct and return system matrix representing discretized curl-curl operator for
  // Ampere's law.
  void GetStiffnessMatrix(std::vector<std::unique_ptr<ParOperator>> &K);

  // Construct and return the discrete curl matrix.
  std::unique_ptr<ParOperator> GetCurlMatrix();

  // Assemble the right-hand side source term vector for a current source applied on
  // specified excited boundaries.
  void GetExcitationVector(int idx, Vector &RHS);
};

}  // namespace palace

#endif  // PALACE_MODELS_CURL_CURL_OPERATOR_HPP
