// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_LAPLACE_OPERATOR_HPP
#define PALACE_MODELS_LAPLACE_OPERATOR_HPP

#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"

namespace palace
{

class IoData;
class Mesh;

//
// A class handling discretization of Laplace problems for electrostatics.
//
class LaplaceOperator
{
private:
  // Helper variable for log file printing.
  bool print_hdr;

  // Essential boundary condition markers.
  mfem::Array<int> dbc_attr;
  std::vector<mfem::Array<int>> dbc_tdof_lists;

  // Objects defining the finite element spaces for the electrostatic potential (H1) and
  // electric field (Nedelec) on the given mesh. The RT spaces are used for error
  // estimation.
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_fecs;
  std::unique_ptr<mfem::ND_FECollection> nd_fec;
  std::vector<std::unique_ptr<mfem::RT_FECollection>> rt_fecs;
  FiniteElementSpaceHierarchy h1_fespaces;
  FiniteElementSpace nd_fespace;
  FiniteElementSpaceHierarchy rt_fespaces;

  // Operator for domain material properties.
  MaterialOperator mat_op;

  // Boundary attributes for each terminal index.
  std::map<int, mfem::Array<int>> source_attr_lists;

  mfem::Array<int> SetUpBoundaryProperties(const IoData &iodata, const mfem::ParMesh &mesh);
  std::map<int, mfem::Array<int>> ConstructSources(const IoData &iodata);

public:
  LaplaceOperator(const IoData &iodata, const std::vector<std::unique_ptr<Mesh>> &mesh);

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access source attribute lists.
  const auto &GetSources() const { return source_attr_lists; }

  // Return the parallel finite element space objects.
  auto &GetH1Spaces() { return h1_fespaces; }
  const auto &GetH1Spaces() const { return h1_fespaces; }
  auto &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  const auto &GetH1Space() const { return h1_fespaces.GetFinestFESpace(); }
  auto &GetNDSpace() { return nd_fespace; }
  const auto &GetNDSpace() const { return nd_fespace; }
  auto &GetRTSpaces() { return rt_fespaces; }
  const auto &GetRTSpaces() const { return rt_fespaces; }
  auto &GetRTSpace() { return rt_fespaces.GetFinestFESpace(); }
  const auto &GetRTSpace() const { return rt_fespaces.GetFinestFESpace(); }

  // Access the underlying mesh object.
  const auto &GetMesh() const { return GetH1Space().GetMesh(); }

  // Return the number of true (conforming) dofs on the finest H1 space.
  auto GlobalTrueVSize() const { return GetH1Space().GlobalTrueVSize(); }

  // Construct and return system matrix representing discretized Laplace operator for
  // Gauss's law.
  std::unique_ptr<Operator> GetStiffnessMatrix();

  // Construct and return the discrete gradient matrix.
  const Operator &GetGradMatrix() const
  {
    return GetNDSpace().GetDiscreteInterpolator(GetH1Space());
  }

  // Assemble the solution boundary conditions and right-hand side vector for a nonzero
  // prescribed voltage on the specified surface index.
  void GetExcitationVector(int idx, const Operator &K, Vector &X, Vector &RHS);

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return GetH1Space().GetComm(); }
};

}  // namespace palace

#endif  // PALACE_MODELS_LAPLACE_OPERATOR_HPP
