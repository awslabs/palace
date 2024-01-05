// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_WAVE_PORT_OPERATOR_HPP
#define PALACE_MODELS_WAVE_PORT_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <unordered_map>
#include <mfem.hpp>
#include "fem/fespace.hpp"
#include "linalg/eps.hpp"
#include "linalg/ksp.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
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
public:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // Wave port properties.
  int mode_idx;
  double d_offset;
  bool excitation;
  std::complex<double> kn0;
  double omega0;

private:
  // SubMesh data structures to define finite element spaces and grid functions on the
  // SubMesh corresponding to this port boundary.
  std::unique_ptr<mfem::ParSubMesh> port_mesh;
  std::unique_ptr<mfem::FiniteElementCollection> port_nd_fec, port_h1_fec;
  std::unique_ptr<FiniteElementSpace> port_nd_fespace, port_h1_fespace;
  std::unique_ptr<mfem::ParTransferMap> port_nd_transfer, port_h1_transfer;
  std::unordered_map<int, int> submesh_parent_elems;

  // List of all boundary attributes making up this port boundary.
  mfem::Array<int> attr_list;

  // Operator storage for repeated boundary mode eigenvalue problem solves.
  double mu_eps_max;
  std::unique_ptr<mfem::HypreParMatrix> A1, A2r, A2i, B1r, B1i, B3;
  ComplexVector v0, e0, e0t, e0n;

  // Eigenvalue solver for boundary modes.
  MPI_Comm port_comm;
  int port_root;
  std::unique_ptr<EigenvalueSolver> eigen;
  std::unique_ptr<ComplexKspSolver> ksp;

  // Stored objects for computing functions of the port modes for use as an excitation or
  // in postprocessing.
  std::unique_ptr<mfem::ParGridFunction> port_S0t;
  std::unique_ptr<mfem::LinearForm> port_sr, port_si;

  // Grid functions storing the last computed electric field mode on the port.
  std::unique_ptr<mfem::ParComplexGridFunction> port_E0t, port_E0n;

public:
  WavePortData(const config::WavePortData &data, const MaterialOperator &mat_op,
               mfem::ParFiniteElementSpace &nd_fespace,
               mfem::ParFiniteElementSpace &h1_fespace, const mfem::Array<int> &dbc_attr);
  ~WavePortData();

  const auto &GetAttrList() const { return attr_list; }

  void Initialize(double omega);

  HYPRE_BigInt GlobalTrueNDSize() const { return port_nd_fespace->GlobalTrueVSize(); }
  HYPRE_BigInt GlobalTrueH1Size() const { return port_h1_fespace->GlobalTrueVSize(); }

  std::unique_ptr<mfem::VectorCoefficient> GetModeExcitationCoefficientReal() const;
  std::unique_ptr<mfem::VectorCoefficient> GetModeExcitationCoefficientImag() const;

  std::unique_ptr<mfem::VectorCoefficient> GetModeFieldCoefficientReal() const;
  std::unique_ptr<mfem::VectorCoefficient> GetModeFieldCoefficientImag() const;

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
                                mfem::ParComplexGridFunction &B) const;
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
  // Mapping from port index to data structure containing port information.
  std::map<int, WavePortData> ports;

  // Flag which forces no printing during WavePortData::Print().
  bool suppress_output;
  double fc, kc;

  void SetUpBoundaryProperties(const IoData &iodata, const MaterialOperator &mat_op,
                               mfem::ParFiniteElementSpace &nd_fespace,
                               mfem::ParFiniteElementSpace &h1_fespace);
  void PrintBoundaryInfo(const IoData &iodata, const mfem::ParMesh &mesh);

  // Compute boundary modes for all wave port boundaries at the specified frequency.
  void Initialize(double omega);

public:
  WavePortOperator(const IoData &iodata, const MaterialOperator &mat_op,
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

  // Returns array of wave port attributes.
  mfem::Array<int> GetAttrList() const;

  // Add contributions to system matrix from wave ports.
  void AddExtraSystemBdrCoefficients(double omega, MaterialPropertyCoefficient &fbr,
                                     MaterialPropertyCoefficient &fbi);

  // Add contributions to the right-hand side source term vector for an incident field at
  // excited port boundaries.
  void AddExcitationBdrCoefficients(double omega, SumVectorCoefficient &fbr,
                                    SumVectorCoefficient &fbi);
};

}  // namespace palace

#endif  // PALACE_MODELS_WAVE_PORT_OPERATOR_HPP
