// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_WAVE_PORT_OPERATOR_HPP
#define PALACE_MODELS_WAVE_PORT_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <mfem.hpp>
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

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

  // Attribute list and marker for all boundary attributes making up this port boundary.
  // Mutable because some MFEM API calls are not const correct.
  mfem::Array<int> attr_list;
  mutable mfem::Array<int> attr_marker;

  // SubMesh data structures to define finite element spaces and grid functions on the
  // SubMesh corresponding to this port boundary.
  std::unique_ptr<mfem::ParSubMesh> port_mesh;
  std::unique_ptr<mfem::FiniteElementCollection> port_nd_fec, port_h1_fec;
  std::unique_ptr<mfem::ParFiniteElementSpace> port_nd_fespace, port_h1_fespace;
  std::unique_ptr<mfem::ParTransferMap> port_nd_transfer, port_h1_transfer;

  // Operator storage for repeated boundary mode eigenvalue problem solves.
  double mu_eps_max;
  std::unique_ptr<mfem::HypreParMatrix> A2r, A2i, B3;
  std::unique_ptr<ComplexOperator> A, B, P;
  ComplexVector v0, e0, e0t, e0n;

  // Eigenvalue solver for boundary modes.
  MPI_Comm port_comm;
  int port_root;
  std::unique_ptr<EigenvalueSolver> eigen;
  std::unique_ptr<ComplexKspSolver> ksp;

  // Grid functions storing the last computed electric field mode on the port and the
  // associated propagation constant. Also the coefficient for the incident port mode
  // (n x H_inc) computed from the electric field mode.
  std::unique_ptr<mfem::ParComplexGridFunction> port_E0t, port_E0n;
  std::unique_ptr<mfem::VectorCoefficient> port_nxH0r_func, port_nxH0i_func;
  std::unique_ptr<mfem::LinearForm> port_sr, port_si;
  std::unique_ptr<mfem::ParGridFunction> port_S0t;
  std::complex<double> kn0;
  double omega0;

public:
  WavePortData(const config::WavePortData &data, const MaterialOperator &mat_op,
               const mfem::ParFiniteElementSpace &nd_fespace,
               const mfem::ParFiniteElementSpace &h1_fespace,
               const mfem::Array<int> &dbc_marker);
  ~WavePortData();

  const mfem::Array<int> &GetMarker() const { return attr_marker; }
  mfem::Array<int> &GetMarker() { return attr_marker; }

  void Initialize(double omega);

  HYPRE_BigInt GlobalTrueNDSize() const { return port_nd_fespace->GlobalTrueVSize(); }
  HYPRE_BigInt GlobalTrueH1Size() const { return port_h1_fespace->GlobalTrueVSize(); }

  std::complex<double> GetPropagationConstant() const { return kn0; }
  double GetOperatingFrequency() const { return omega0; }

  bool IsExcited() const { return excitation; }
  int GetModeIndex() const { return mode_idx; }
  double GetOffsetDistance() const { return d_offset; }

  const mfem::VectorCoefficient &GetModeCoefficientReal() const { return *port_nxH0r_func; }
  mfem::VectorCoefficient &GetModeCoefficientReal() { return *port_nxH0r_func; }
  const mfem::VectorCoefficient &GetModeCoefficientImag() const { return *port_nxH0i_func; }
  mfem::VectorCoefficient &GetModeCoefficientImag() { return *port_nxH0i_func; }

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
                                const MaterialOperator &mat_op) const;
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
                               const mfem::ParFiniteElementSpace &nd_fespace,
                               const mfem::ParFiniteElementSpace &h1_fespace);
  void PrintBoundaryInfo(const IoData &iodata, mfem::ParMesh &mesh);

  // Compute boundary modes for all wave port boundaries at the specified frequency.
  void Initialize(double omega);

public:
  WavePortOperator(const IoData &iod, const MaterialOperator &mat,
                   const mfem::ParFiniteElementSpace &nd_fespace,
                   const mfem::ParFiniteElementSpace &h1_fespace);

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
