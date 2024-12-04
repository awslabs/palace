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

class GridFunction;
class IoData;
class MaterialOperator;
class MaterialPropertyCoefficient;
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
public:
  // Reference to material property data (not owned).
  const MaterialOperator &mat_op;

  // To accomodate multielement lumped ports, a port may be made up of elements with
  // different attributes and directions which add in parallel.
  std::vector<std::unique_ptr<LumpedElementData>> elems;

  // Lumped port properties.
  double R, L, C;
  bool active;
  int excitation;

private:
  // Linear forms for postprocessing integrated quantities on the port.
  mutable std::unique_ptr<mfem::LinearForm> s, v;

  void InitializeLinearForms(mfem::ParFiniteElementSpace &nd_fespace) const;

public:
  LumpedPortData(const config::LumpedPortData &data, const MaterialOperator &mat_op,
                 const mfem::ParMesh &mesh);

  double GetToSquare(const LumpedElementData &elem) const
  {
    return elem.GetGeometryWidth() / elem.GetGeometryLength() * elems.size();
  }

  enum class Branch
  {
    TOTAL,
    R,
    L,
    C
  };
  std::complex<double> GetCharacteristicImpedance(double omega = 0.0,
                                                  Branch branch = Branch::TOTAL) const;

  double GetExcitationPower() const;
  double GetExcitationVoltage() const;

  std::complex<double> GetPower(GridFunction &E, GridFunction &B) const;
  std::complex<double> GetSParameter(GridFunction &E) const;
  std::complex<double> GetVoltage(GridFunction &E) const;
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

  void SetUpBoundaryProperties(const IoData &iodata, const MaterialOperator &mat_op,
                               const mfem::ParMesh &mesh);
  void PrintBoundaryInfo(const IoData &iodata, const mfem::ParMesh &mesh);

public:
  LumpedPortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                     const mfem::ParMesh &mesh);

  // Access data structures for the lumped port with the given index.
  const LumpedPortData &GetPort(int idx) const;
  auto begin() const { return ports.begin(); }
  auto end() const { return ports.end(); }
  auto rbegin() const { return ports.rbegin(); }
  auto rend() const { return ports.rend(); }
  auto Size() const { return ports.size(); }

  // Returns array of lumped port attributes.
  mfem::Array<int> GetAttrList() const;
  mfem::Array<int> GetRsAttrList() const;
  mfem::Array<int> GetLsAttrList() const;
  mfem::Array<int> GetCsAttrList() const;

  // Add contributions to system matrices from lumped elements with nonzero inductance,
  // resistance, and/or capacitance.
  void AddStiffnessBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddDampingBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);
  void AddMassBdrCoefficients(double coeff, MaterialPropertyCoefficient &fb);

  // Add contributions to the right-hand side source term vector for an incident field at
  // excited port boundaries, -U_inc/(iω) for the real version (versus the full -U_inc for
  // the complex one).
  void AddExcitationBdrCoefficients(SumVectorCoefficient &fb);
};

}  // namespace palace

#endif  // PALACE_MODELS_LUMPED_PORT_OPERATOR_HPP
