// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_SPACE_OPERATOR_HPP
#define PALACE_SPACE_OPERATOR_HPP

#include <functional>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/farfieldboundaryoperator.hpp"
#include "fem/lumpedportoperator.hpp"
#include "fem/materialoperator.hpp"
#include "fem/surfaceconductivityoperator.hpp"
#include "fem/surfacecurrentoperator.hpp"
#include "fem/surfaceimpedanceoperator.hpp"
#include "fem/waveportoperator.hpp"

namespace palace
{

class IoData;
class SumCoefficient;
class SumMatrixCoefficient;

namespace petsc
{

class PetscParMatrix;
class PetscParVector;

}  // namespace petsc

//
// A class handling spatial discretization of the governing equations.
//
class SpaceOperator
{
private:
  // Perfect electrical conductor essential boundary condition markers.
  mfem::Array<int> dbc_marker, dbc_tdof_list, aux_bdr_marker;
  void CheckBoundaryProperties();

  // Options for system matrix assembly.
  const int skip_zeros;   // Whether to skip the zeros during assembly of operators
  const bool pc_gmg;      // Whether to use geometric multigrid in preconditioning
  const bool pc_lor;      // Whether to use low-order refined (LOR) preconditioner
  const bool pc_shifted;  // Whether the preconditioner uses the shifted mass matrix

  // Helper variable and function for log file printing.
  bool print_hdr;
  void PrintHeader();

  // Objects defining the finite element spaces for the electric field(Nedelec) and magnetic
  // flux density (Raviart-Thomas) on the given mesh. The H1 spaces are used for various
  // purposes throughout the code including postprocessing.
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

  // Helper function to assemble preconditioner matrix data structures.
  void GetPreconditionerInternal(
      const std::function<void(SumMatrixCoefficient &, SumMatrixCoefficient &,
                               SumCoefficient &, SumMatrixCoefficient &)> &AddCoefficients,
      std::vector<std::unique_ptr<mfem::Operator>> &B,
      std::vector<std::unique_ptr<mfem::Operator>> &AuxB, bool print);

  // Helper functions for building the bilinear forms corresponding to the discretized
  // operators in Maxwell's equations.
  void AddStiffnessCoefficients(double coef, SumMatrixCoefficient &df,
                                SumMatrixCoefficient &f, SumMatrixCoefficient &fb);
  void AddRealMassCoefficients(double coef, bool abs_coef, SumMatrixCoefficient &f,
                               SumMatrixCoefficient &fb);
  void AddImagMassCoefficients(double coef, SumMatrixCoefficient &f,
                               SumMatrixCoefficient &fb);
  void AddDampingCoefficients(double coef, SumMatrixCoefficient &f,
                              SumMatrixCoefficient &fb);
  void AddExtraSystemBdrCoefficients(double omega, SumCoefficient &dfbr,
                                     SumCoefficient &dfbi, SumMatrixCoefficient &fbr,
                                     SumMatrixCoefficient &fbi);

  // Helper functions for excitation vector assembly.
  bool GetExcitationVector1Internal(mfem::Vector &RHS);
  bool GetExcitationVector2Internal(double omega, mfem::Vector &RHSr, mfem::Vector &RHSi);

public:
  SpaceOperator(const IoData &iodata,
                const std::vector<std::unique_ptr<mfem::ParMesh>> &mesh);

  // Returns array marking Dirichlet BC (PEC) attributes and list of local true dofs.
  const mfem::Array<int> &GetDbcMarker() const { return dbc_marker; }
  const mfem::Array<int> &GetDbcTDofList() const { return dbc_tdof_list; }

  // Returns array marking all boundary condition attributes, PEC included. These are all
  // boundaries which affect the stiffness and damping (K and C) matrices, used for
  // nullspace corrections.
  const mfem::Array<int> &GetAuxBdrMarker() const { return aux_bdr_marker; }

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
  auto &GetH1Spaces() { return h1_fespaces; }
  auto &GetH1Space() { return h1_fespaces.GetFinestFESpace(); }
  auto &GetRTSpace() { return rt_fespace; }

  // Construct the frequency-dependent complex linear system matrix:
  //                 A = K + iω C - ω² (Mr + i Mi) + A2(ω)
  // or any one of its terms. The type parameter controls which terms of the above
  // formulation to include in the resulting matrix. The argument ω is only used for
  // the "complete" or "extra" system matrix options, all others come unscaled.
  enum class OperatorType
  {
    COMPLETE,
    STIFFNESS,
    MASS,
    DAMPING,
    EXTRA
  };
  std::unique_ptr<petsc::PetscParMatrix>
  GetSystemMatrixPetsc(OperatorType type, double omega,
                       mfem::Operator::DiagonalPolicy ess_diag, bool print = true);
  std::unique_ptr<petsc::PetscParMatrix>
  GetSystemMatrixPetsc(OperatorType type, mfem::Operator::DiagonalPolicy ess_diag,
                       bool print = true)
  {
    return GetSystemMatrixPetsc(type, 0.0, ess_diag, print);
  }
  std::unique_ptr<mfem::Operator> GetSystemMatrix(OperatorType type, double omega,
                                                  mfem::Operator::DiagonalPolicy ess_diag,
                                                  bool print = true);
  std::unique_ptr<mfem::Operator> GetSystemMatrix(OperatorType type,
                                                  mfem::Operator::DiagonalPolicy ess_diag,
                                                  bool print = true)
  {
    return GetSystemMatrix(type, 0.0, ess_diag, print);
  }

  // Construct the real, optionally SPD matrix for frequency or time domain preconditioning
  // (Mr > 0, Mi < 0):
  //              B =    K +  ω C + ω² (-/+ Mr - Mi) , or
  //              B = a0 K + a1 C +         Mr .
  void GetPreconditionerMatrix(double omega,
                               std::vector<std::unique_ptr<mfem::Operator>> &B,
                               std::vector<std::unique_ptr<mfem::Operator>> &AuxB,
                               bool print = true);
  void GetPreconditionerMatrix(double a0, double a1,
                               std::vector<std::unique_ptr<mfem::Operator>> &B,
                               std::vector<std::unique_ptr<mfem::Operator>> &AuxB,
                               bool print = true);

  // Construct and return the discrete negative curl or gradient matrices.
  std::unique_ptr<mfem::Operator> GetNegCurlMatrix();
  std::unique_ptr<petsc::PetscParMatrix> GetNegCurlMatrixPetsc();
  std::unique_ptr<mfem::Operator> GetGradMatrix();
  std::unique_ptr<petsc::PetscParMatrix> GetGradMatrixPetsc();

  // Assemble the right-hand side source term vector for an incident field or current source
  // applied on specified excited boundaries.
  bool GetTimeDomainExcitationVector(mfem::Vector &RHS);
  bool GetFreqDomainExcitationVector(double omega, petsc::PetscParVector &RHS);

  // Separate out RHS vector as RHS = iω RHS1 + RHS2(ω).
  bool GetFreqDomainExcitationVector1(petsc::PetscParVector &RHS1);
  bool GetFreqDomainExcitationVector2(double omega, petsc::PetscParVector &RHS2);

  // Compute elemental error estimates given a provided electric field, E
  std::vector<double> GetErrorEstimates(const mfem::ParComplexGridFunction &E);
};

}  // namespace palace

#endif  // PALACE_SPACE_OPERATOR_HPP
