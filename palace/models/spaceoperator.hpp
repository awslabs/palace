// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_SPACE_OPERATOR_HPP
#define PALACE_MODELS_SPACE_OPERATOR_HPP

#include <complex>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/coefficient.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/lumpedportoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfacecurrentoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/waveportoperator.hpp"

namespace palace
{

class IoData;

//
// A class handling spatial discretization of the governing equations.
//
class SpaceOperator
{
private:
  const int pa_order_threshold;  // Order above which to use partial assembly vs. full
  const int skip_zeros;          // Skip zeros during full assembly of operators
  const bool pc_mat_real;        // Use real-valued matrix for preconditioner
  const bool pc_mat_shifted;     // Use shifted mass matrix for preconditioner
  const bool pc_mat_lor;         // Use low-order refined (LOR) space for preconditioner

  // Helper variables for log file printing.
  bool print_hdr, print_prec_hdr;

  // Perfect electrical conductor essential boundary condition markers.
  mfem::Array<int> dbc_marker, aux_bdr_marker;
  std::vector<mfem::Array<int>> nd_dbc_tdof_lists, h1_dbc_tdof_lists, aux_bdr_tdof_lists;
  void CheckBoundaryProperties();

  // Objects defining the finite element spaces for the electric field (Nedelec) and
  // magnetic flux density (Raviart-Thomas) on the given mesh. The H1 spaces are used for
  // various purposes throughout the code including postprocessing.
  std::vector<std::unique_ptr<mfem::ND_FECollection>> nd_fecs;
  std::vector<std::unique_ptr<mfem::H1_FECollection>> h1_fecs;
  mfem::RT_FECollection rt_fec;
  mfem::ParFiniteElementSpaceHierarchy nd_fespaces, h1_fespaces;
  mfem::ParFiniteElementSpace rt_fespace;

  // Operator for domain material properties.
  MaterialOperator mat_op;

  // Operators for boundary conditions and source excitations.
  FarfieldBoundaryOperator farfield_op;
  SurfaceConductivityOperator surf_sigma_op;
  SurfaceImpedanceOperator surf_z_op;
  LumpedPortOperator lumped_port_op;
  WavePortOperator wave_port_op;
  SurfaceCurrentOperator surf_j_op;

  // Helper functions for building the bilinear forms corresponding to the discretized
  // operators in Maxwell's equations.
  void AddStiffnessCoefficients(double coef, SumMatrixCoefficient &df,
                                SumMatrixCoefficient &f);
  void AddStiffnessBdrCoefficients(double coef, SumMatrixCoefficient &fb);
  void AddDampingCoefficients(double coef, SumMatrixCoefficient &f);
  void AddDampingBdrCoefficients(double coef, SumMatrixCoefficient &fb);
  void AddRealMassCoefficients(double coef, SumMatrixCoefficient &f);
  void AddRealMassBdrCoefficients(double coef, SumMatrixCoefficient &fb);
  void AddImagMassCoefficients(double coef, SumMatrixCoefficient &f);
  void AddAbsMassCoefficients(double coef, SumMatrixCoefficient &f);
  void AddExtraSystemBdrCoefficients(double omega, SumCoefficient &dfbr,
                                     SumCoefficient &dfbi, SumMatrixCoefficient &fbr,
                                     SumMatrixCoefficient &fbi);

  // Helper functions for excitation vector assembly.
  bool AddExcitationVector1Internal(Vector &RHS);
  bool AddExcitationVector2Internal(double omega, ComplexVector &RHS);

public:
  SpaceOperator(const IoData &iodata,
                const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh);

  // Return list of all PEC boundary true dofs for all finite element space levels.
  const std::vector<mfem::Array<int>> &GetNDDbcTDofLists() const
  {
    return nd_dbc_tdof_lists;
  }
  const std::vector<mfem::Array<int>> &GetH1DbcTDofLists() const
  {
    return h1_dbc_tdof_lists;
  }

  // Returns lists of all boundary condition true dofs, PEC included, for the auxiliary
  // H1 space hierarchy. These are all boundaries which affect the stiffness and damping
  // (K and C) matrices, used for nullspace corrections.
  const std::vector<mfem::Array<int>> &GetAuxBdrTDofLists() const
  {
    return aux_bdr_tdof_lists;
  }

  // Return material operator for postprocessing.
  const MaterialOperator &GetMaterialOp() const { return mat_op; }

  // Access to underlying BC operator objects for postprocessing.
  auto &GetLumpedPortOp() { return lumped_port_op; }
  auto &GetWavePortOp() { return wave_port_op; }
  auto &GetSurfaceCurrentOp() { return surf_j_op; }
  const auto &GetLumpedPortOp() const { return lumped_port_op; }
  const auto &GetWavePortOp() const { return wave_port_op; }
  const auto &GetSurfaceCurrentOp() const { return surf_j_op; }

  // Return the parallel finite element space objects.
  auto &GetNDSpaces() { return nd_fespaces; }
  auto &GetNDSpace() { return nd_fespaces.GetFinestFESpace(); }
  const auto &GetNDSpace() const { return nd_fespaces.GetFinestFESpace(); }
  auto &GetH1Spaces() { return h1_fespaces; }
  auto &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  const auto &GetH1Space() const { return h1_fespaces.GetFinestFESpace(); }
  auto &GetRTSpace() { return rt_fespace; }
  const auto &GetRTSpace() const { return rt_fespace; }

  // Return the number of true (conforming) dofs on the finest ND space.
  auto GlobalTrueVSize() { return GetNDSpace().GlobalTrueVSize(); }

  // Construct any part of the frequency-dependent complex linear system matrix:
  //                     A = K + iω C - ω² (Mr + i Mi) + A2(ω) .
  // For time domain problems, any one of K, C, or M = Mr can be constructed. The argument
  // ω is required only for the constructing the "extra" matrix A2(ω).
  template <typename OperType>
  std::unique_ptr<OperType> GetStiffnessMatrix(Operator::DiagonalPolicy diag_policy);
  template <typename OperType>
  std::unique_ptr<OperType> GetDampingMatrix(Operator::DiagonalPolicy diag_policy);
  template <typename OperType>
  std::unique_ptr<OperType> GetMassMatrix(Operator::DiagonalPolicy diag_policy);
  template <typename OperType>
  std::unique_ptr<OperType> GetExtraSystemMatrix(double omega,
                                                 Operator::DiagonalPolicy diag_policy);

  // Construct the complete frequency or time domain system matrix using the provided
  // stiffness, damping, mass, and extra matrices:
  //                     A = a0 K + a1 C + a2 (Mr + i Mi) + A2 .
  // It is assumed that the inputs have been constructed using previous calls to
  // GetSystemMatrix() and the returned operator does not inherit ownership of any of them.
  template <typename OperType, typename ScalarType>
  std::unique_ptr<OperType>
  GetSystemMatrix(ScalarType a0, ScalarType a1, ScalarType a2, const OperType *K,
                  const OperType *C, const OperType *M, const OperType *A2 = nullptr);

  // Construct the real, SPD matrix for weighted L2 or H(curl) inner products:
  //                           B = a0 Kr + a2 Mr .
  // It is assumed that the inputs have been constructed using previous calls to
  // GetSystemMatrix() and the returned operator does not inherit ownership of any of them.
  // If K or M have eliminated boundary conditions, they are not eliminated from the
  // returned operator.
  std::unique_ptr<Operator> GetInnerProductMatrix(double a0, double a2,
                                                  const ComplexOperator *K,
                                                  const ComplexOperator *M);

  // Construct the real, optionally SPD matrix for frequency or time domain linear system
  // preconditioning (Mr > 0, Mi < 0, |Mr + i Mi| is done on the material property
  // coefficient, not the matrix entries themselves):
  //             B = a0 K + a1 C -/+ a2 |Mr + i Mi| + A2r(a3) + A2i(a3) .
  template <typename OperType>
  std::unique_ptr<OperType> GetPreconditionerMatrix(double a0, double a1, double a2,
                                                    double a3);

  // Construct and return the discrete curl or gradient matrices. The complex variants
  // return a matrix suitable for applying to complex-valued vectors.
  template <typename OperType>
  std::unique_ptr<OperType> GetCurlMatrix();
  template <typename OperType>
  std::unique_ptr<OperType> GetGradMatrix();

  // Assemble the right-hand side source term vector for an incident field or current source
  // applied on specified excited boundaries. The return value indicates whether or not the
  // excitation is nonzero (and thus is true most of the time).
  bool GetExcitationVector(Vector &RHS);
  bool GetExcitationVector(double omega, ComplexVector &RHS);

  // Separate out RHS vector as RHS = iω RHS1 + RHS2(ω). The return value indicates whether
  // or not the excitation is nonzero (and thus is true most of the time).
  bool GetExcitationVector1(ComplexVector &RHS1);
  bool GetExcitationVector2(double omega, ComplexVector &RHS2);

  // Construct a constant or randomly initialized vector which satisfies the PEC essential
  // boundary conditions.
  void GetRandomInitialVector(ComplexVector &v);
  void GetConstantInitialVector(ComplexVector &v);

  // Get the associated MPI communicator.
  MPI_Comm GetComm() const { return GetNDSpace().GetComm(); }
};

}  // namespace palace

#endif  // PALACE_MODELS_SPACE_OPERATOR_HPP
