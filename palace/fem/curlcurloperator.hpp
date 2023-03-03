// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_CURL_CURL_OPERATOR_HPP
#define PALACE_CURL_CURL_OPERATOR_HPP

#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/materialoperator.hpp"
#include "fem/surfacecurrentoperator.hpp"

namespace palace
{

class IoData;

//
// A class handling discretization of curl-curl problems for magnetostatics.
//
class CurlCurlOperator
{
private:
  // Essential boundary condition markers.
  mfem::Array<int> dbc_marker, dbc_tdof_list;
  void CheckBoundaryProperties();

  // Options for system matrix assembly.
  const int skip_zeros;  // Whether to skip the zeros during assembly of operators
  const bool pc_gmg;     // Whether to use geometric multigrid in preconditioning

  // Helper variable and function for log file printing.
  bool print_hdr;
  void PrintHeader();

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

  // Returns array marking Dirichlet BC attributes and local subdomain vdofs.
  const auto &GetDbcMarker() const { return dbc_marker; }
  const auto &GetDbcTDofList() const { return dbc_tdof_list; }

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access to underlying BC operator objects for postprocessing.
  const auto &GetSurfaceCurrentOp() const { return surf_j_op; }

  // Return the parallel finite element space objects.
  auto &GetNDSpaces() { return nd_fespaces; }
  auto &GetNDSpace() { return nd_fespaces.GetFinestFESpace(); }
  auto &GetH1Space() { return h1_fespace; }
  auto &GetRTSpace() { return rt_fespace; }

  // Return the number of true (conforming) dofs on the finest ND space.
  auto GetNDof() { return GetNDSpace().GetNConformingDofs(); }

  // Construct and return system matrix representing discretized curl-curl operator for
  // Ampere's law.
  void GetStiffnessMatrix(std::vector<std::unique_ptr<mfem::Operator>> &K);

  // Construct and return the discrete curl matrix.
  std::unique_ptr<mfem::Operator> GetCurlMatrix();

  // Assemble the right-hand side source term vector for a current source applied on
  // specified excited boundaries.
  void GetExcitationVector(int idx, mfem::Vector &RHS);
};

}  // namespace palace

#endif  // PALACE_CURL_CURL_OPERATOR_HPP
