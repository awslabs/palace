// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_LUMPED_PORT_OPERATOR_HPP
#define PALACE_MODELS_LUMPED_PORT_OPERATOR_HPP

#include <complex>
#include <map>
#include <memory>
#include <vector>
#include <mfem.hpp>
#include "fem/lumpedelement.hpp"

namespace palace
{

class IoData;
class MaterialOperator;
class SumMatrixCoefficient;
class SumVectorCoefficient;

namespace config
{

struct LumpedPortData;

}  // namespace config

//
// Helper class for lumped ports in a model.
//
class LumpedPortData
{
private:
  bool excitation;
  double R, L, C;

  // To accomodate multielement lumped ports, a port may be made up of elements with
  // different attributes and directions which add in parallel.
  std::vector<std::unique_ptr<LumpedElementData>> elems;

  // Linear forms for postprocessing integrated quantities on the port.
  mutable std::unique_ptr<mfem::ParLinearForm> s, v;

public:
  LumpedPortData(const config::LumpedPortData &data,
                 mfem::ParFiniteElementSpace &h1_fespace);

  const std::vector<std::unique_ptr<LumpedElementData>> &GetElements() const
  {
    return elems;
  }

  double GetToSquare(const LumpedElementData &elem) const
  {
    return elem.GetGeometryWidth() / elem.GetGeometryLength() * elems.size();
  }

  bool IsExcited() const { return excitation; }
  double GetR() const { return R; }
  double GetL() const { return L; }
  double GetC() const { return C; }

  std::complex<double> GetCharacteristicImpedance(double omega = 0.0) const;

  double GetExcitationPower() const;
  double GetExcitationVoltage() const;

  std::complex<double> GetSParameter(mfem::ParComplexGridFunction &E) const;
  std::complex<double> GetPower(mfem::ParComplexGridFunction &E,
                                mfem::ParComplexGridFunction &B,
                                const MaterialOperator &mat_op) const;
  double GetPower(mfem::ParGridFunction &E, mfem::ParGridFunction &B,
                  const MaterialOperator &mat_op) const;
  std::complex<double> GetVoltage(mfem::ParComplexGridFunction &E) const;
  double GetVoltage(mfem::ParGridFunction &E) const;
};

//
// A class handling lumped port boundaries and their postprocessing.
//
class LumpedPortOperator
{
private:
  // Mapping from port index to data structure containing port information and methods to
  // calculate circuit properties like voltage and current on lumped or multielement lumped
  // ports.
  std::map<int, LumpedPortData> ports;
  mfem::Array<int> port_marker, port_Rs_marker, port_Ls_marker, port_Cs_marker;
  void SetUpBoundaryProperties(const IoData &iodata,
                               mfem::ParFiniteElementSpace &h1_fespace);
  void PrintBoundaryInfo(const IoData &iodata, mfem::ParMesh &mesh);

public:
  LumpedPortOperator(const IoData &iodata, mfem::ParFiniteElementSpace &h1_fespace);

  // Access data structures for the lumped port with the given index.
  const LumpedPortData &GetPort(int idx) const;
  auto begin() const { return ports.begin(); }
  auto end() const { return ports.end(); }
  auto rbegin() const { return ports.rbegin(); }
  auto rend() const { return ports.rend(); }
  auto Size() const { return ports.size(); }

  // Returns array marking lumped port attributes.
  const mfem::Array<int> &GetMarker() const { return port_marker; }
  const mfem::Array<int> &GetRsMarker() const { return port_Rs_marker; }
  const mfem::Array<int> &GetLsMarker() const { return port_Ls_marker; }
  const mfem::Array<int> &GetCsMarker() const { return port_Cs_marker; }

  // Add contributions to system matrices from lumped elements with nonzero inductance,
  // capacitance, and/or resistance.
  void AddStiffnessBdrCoefficients(double coef, SumMatrixCoefficient &fb);
  void AddMassBdrCoefficients(double coef, SumMatrixCoefficient &fb);
  void AddDampingBdrCoefficients(double coef, SumMatrixCoefficient &fb);

  // Add contributions to the right-hand side source term vector for an incident field at
  // excited port boundaries, -U_inc/(iω) for the real version (versus the full -U_inc for
  // the complex one).
  void AddExcitationBdrCoefficients(SumVectorCoefficient &fb);
};

}  // namespace palace

#endif  // PALACE_MODELS_LUMPED_PORT_OPERATOR_HPP
