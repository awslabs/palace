// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_WAVE_PORT_OPERATOR_HPP
#define PALACE_MODELS_WAVE_PORT_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <mfem.hpp>
#include "linalg/eigen.hpp"
#include "linalg/ksp.hpp"
#include "linalg/petsc.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class SumMatrixCoefficient;
class SumVectorCoefficient;

namespace config
{

struct WavePortData;

}  // namespace config

//
// Helper class for wave port boundaries in a model.
//
class WavePortData
{
private:
  bool excitation;
  int mode_idx;
  double d_offset;

  // Marker for all boundary attributes making up this port boundary. Mutable because
  // some MFEM API calls are not const correct.
  mutable mfem::Array<int> attr_marker;

  // Lists of non-essential true degrees of freedom associated with the port boundary.
  mfem::Array<int> nd_attr_tdof_list, h1_attr_tdof_list;

  // Operator storage for repeated boundary mode eigenvalue problem solves.
  std::unique_ptr<petsc::PetscParMatrix> A, B, A1, A2, B3, B4;
  std::unique_ptr<petsc::PetscParVector> e, e0, y0;
  std::unique_ptr<petsc::PetscScatter> scatter;
  double muepsmax;

  // Grid functions storing the last computed electric field mode on the port and the
  // associated propagation constant.
  std::unique_ptr<mfem::ParComplexGridFunction> E0t, E0n;
  std::complex<double> kn0;
  double omega0;

  // Coefficients storing the incident port mode (n x H_inc) and linear forms for
  // postprocessing integrated quantities on the port.
  std::unique_ptr<mfem::VectorCoefficient> nxH0r_func, nxH0i_func;
  std::unique_ptr<mfem::ParLinearForm> sr, si;

  // Eigenvalue solver for boundary modes.
  std::unique_ptr<EigenSolverBase> eigen;
  std::unique_ptr<KspSolver> ksp;

  // Helper function to get true degrees of freedom on the port.
  void GetTrueDofs(const mfem::Array<int> &dbc_marker,
                   mfem::ParFiniteElementSpace &nd_fespace,
                   mfem::ParFiniteElementSpace &h1_fespace, mfem::Array<int> &nd_tdof_list,
                   mfem::Array<int> &h1_tdof_list);

  // Configure and solve the linear eigenvalue problem for the boundary mode.
  void GetInitialSpace(int nt, int nn, petsc::PetscParVector &y0);
  std::complex<double> Solve(petsc::PetscParVector &y0, petsc::PetscParVector &e0,
                             petsc::PetscParVector &e, petsc::PetscScatter &scatter);

public:
  WavePortData(const config::WavePortData &data, const MaterialOperator &mat_op,
               const mfem::Array<int> &dbc_marker, mfem::ParFiniteElementSpace &nd_fespace,
               mfem::ParFiniteElementSpace &h1_fespace);

  const mfem::Array<int> &GetMarker() const { return attr_marker; }
  mfem::Array<int> &GetMarker() { return attr_marker; }

  void Initialize(double omega);

  const petsc::PetscParMatrix *GetA() const { return A.get(); }
  const petsc::PetscParMatrix *GetB() const { return B.get(); }

  std::complex<double> GetPropagationConstant() const { return kn0; }
  double GetOperatingFrequency() const { return omega0; }

  bool IsExcited() const { return excitation; }
  int GetModeIndex() const { return mode_idx; }
  double GetOffsetDistance() const { return d_offset; }

  const std::unique_ptr<mfem::VectorCoefficient> &GetModeCoefficientReal() const
  {
    return nxH0r_func;
  }
  const std::unique_ptr<mfem::VectorCoefficient> &GetModeCoefficientImag() const
  {
    return nxH0i_func;
  }

  std::complex<double> GetCharacteristicImpedance() const
  {
    MFEM_ABORT("GetImpedance is not yet implemented for wave port boundaries!");
    return 0.0;
  }

  double GetExcitationPower() const;
  std::complex<double> GetExcitationVoltage() const
  {
    MFEM_ABORT("GetExcitationVoltage is not yet implemented for wave port boundaries!");
    return 0.0;
  }

  std::complex<double> GetSParameter(mfem::ParComplexGridFunction &E) const;
  std::complex<double> GetPower(mfem::ParComplexGridFunction &E,
                                mfem::ParComplexGridFunction &B,
                                const MaterialOperator &mat_op,
                                const std::map<int, int> &local_to_shared) const;
  std::complex<double> GetVoltage(mfem::ParComplexGridFunction &E) const
  {
    MFEM_ABORT("GetVoltage is not yet implemented for wave port boundaries!");
    return 0.0;
  }
};

//
// A class handling wave port boundaries and their postprocessing.
//
class WavePortOperator
{
private:
  // References to configuration file and material property data (not owned).
  const IoData &iodata;
  const MaterialOperator &mat_op;

  // Flag which forces no printing during WavePortData::Print().
  bool suppress_output;

  // Mapping from port index to data structure containing port information.
  std::map<int, WavePortData> ports;
  mfem::Array<int> port_marker;
  void SetUpBoundaryProperties(const IoData &iodata, const MaterialOperator &mat_op,
                               mfem::ParFiniteElementSpace &nd_fespace,
                               mfem::ParFiniteElementSpace &h1_fespace);
  void PrintBoundaryInfo(const IoData &iodata, mfem::ParMesh &mesh);

  // Compute boundary modes for all wave port boundaries at the specified frequency.
  void Initialize(double omega);

public:
  WavePortOperator(const IoData &iod, const MaterialOperator &mat,
                   mfem::ParFiniteElementSpace &nd_fespace,
                   mfem::ParFiniteElementSpace &h1_fespace);

  // Access data structures for the wave port with the given index.
  const WavePortData &GetPort(int idx) const;
  auto begin() const { return ports.begin(); }
  auto end() const { return ports.end(); }
  auto rbegin() const { return ports.rbegin(); }
  auto rend() const { return ports.rend(); }
  auto Size() const { return ports.size(); }

  // Enable or suppress all outputs (log printing and fields to disk).
  void SetSuppressOutput(bool suppress) { suppress_output = suppress; }

  // Returns array marking wave port attributes.
  const mfem::Array<int> &GetMarker() const { return port_marker; }

  // Add contributions to system matrix from wave ports.
  void AddExtraSystemBdrCoefficients(double omega, SumMatrixCoefficient &fbr,
                                     SumMatrixCoefficient &fbi);

  // Add contributions to the right-hand side source term vector for an incident field at
  // excited port boundaries.
  void AddExcitationBdrCoefficients(double omega, SumVectorCoefficient &fbr,
                                    SumVectorCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_MODELS_WAVE_PORT_OPERATOR_HPP