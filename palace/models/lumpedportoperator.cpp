// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "lumpedportoperator.hpp"

#include "fem/coefficient.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

LumpedPortData::LumpedPortData(const config::LumpedPortData &data,
                               const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op), active(data.active), excitation(data.excitation)
{
  // Check inputs. Only one of the circuit or per square properties should be specified
  // for the port boundary.
  bool has_circ = (std::abs(data.R) + std::abs(data.L) + std::abs(data.C) > 0.0);
  bool has_surf = (std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) > 0.0);
  MFEM_VERIFY(has_circ || has_surf,
              "Lumped port boundary has no R/L/C or Rs/Ls/Cs defined, needs "
              "at least one!");
  MFEM_VERIFY(!(has_circ && has_surf),
              "Lumped port boundary has both R/L/C and Rs/Ls/Cs defined, "
              "should only use one!");

  if (excitation != 0)
  {
    if (has_circ)
    {
      MFEM_VERIFY(data.R > 0.0, "Excited lumped port must have nonzero resistance!");
      MFEM_VERIFY(data.C == 0.0 && data.L == 0.0,
                  "Lumped port excitations do not support nonzero reactance!");
    }
    else
    {
      MFEM_VERIFY(data.Rs > 0.0, "Excited lumped port must have nonzero resistance!");
      MFEM_VERIFY(data.Cs == 0.0 && data.Ls == 0.0,
                  "Lumped port excitations do not support nonzero reactance!");
    }
  }

  // Construct the port elements allowing for a possible multielement lumped port.
  for (const auto &elem : data.elements)
  {
    mfem::Array<int> attr_list;
    attr_list.Append(elem.attributes.data(), elem.attributes.size());
    switch (elem.coordinate_system)
    {
      case config::internal::ElementData::CoordinateSystem::CYLINDRICAL:
        elems.push_back(
            std::make_unique<CoaxialElementData>(elem.direction, attr_list, mesh));
        break;
      case config::internal::ElementData::CoordinateSystem::CARTESIAN:
        elems.push_back(
            std::make_unique<UniformElementData>(elem.direction, attr_list, mesh));
        break;
    }
  }

  // Populate the property data for the lumped port.
  if (std::abs(data.Rs) + std::abs(data.Ls) + std::abs(data.Cs) == 0.0)
  {
    R = data.R;
    L = data.L;
    C = data.C;
  }
  else
  {
    // If defined by surface properties, need to compute circuit properties for the
    // multielement port.
    double ooR = 0.0, ooL = 0.0;
    R = L = C = 0.0;
    for (const auto &elem : elems)
    {
      const double sq = elem->GetGeometryWidth() / elem->GetGeometryLength();
      if (std::abs(data.Rs) > 0.0)
      {
        ooR += sq / data.Rs;
      }
      if (std::abs(data.Ls) > 0.0)
      {
        ooL += sq / data.Ls;
      }
      if (std::abs(data.Cs) > 0.0)
      {
        C += sq * data.Cs;
      }
    }
    if (std::abs(ooR) > 0.0)
    {
      R = 1.0 / ooR;
    }
    if (std::abs(ooL) > 0.0)
    {
      L = 1.0 / ooL;
    }
  }
}

std::complex<double>
LumpedPortData::GetCharacteristicImpedance(double omega,
                                           LumpedPortData::Branch branch) const
{
  MFEM_VERIFY((L == 0.0 && C == 0.0) || branch == Branch::R || omega > 0.0,
              "Lumped port with nonzero reactance requires frequency in order to define "
              "characteristic impedance!");
  std::complex<double> Y = 0.0;
  if (std::abs(R) > 0.0 && (branch == Branch::TOTAL || branch == Branch::R))
  {
    Y += 1.0 / R;
  }
  if (std::abs(L) > 0.0 && (branch == Branch::TOTAL || branch == Branch::L))
  {
    Y += 1.0 / (1i * omega * L);
  }
  if (std::abs(C) > 0.0 && (branch == Branch::TOTAL || branch == Branch::C))
  {
    Y += 1i * omega * C;
  }
  MFEM_VERIFY(std::abs(Y) > 0.0,
              "Characteristic impedance requested for lumped port with zero admittance!")
  return 1.0 / Y;
}

double LumpedPortData::GetExcitationPower() const
{
  // The lumped port excitation is normalized such that the power integrated over the port
  // is 1: ∫ (E_inc x H_inc) ⋅ n dS = 1.
  return excitation ? 1.0 : 0.0;
}

double LumpedPortData::GetExcitationVoltage() const
{
  // Incident voltage should be the same across all elements of an excited lumped port.
  if (excitation)
  {
    double V_inc = 0.0;
    for (const auto &elem : elems)
    {
      const double Rs = R * GetToSquare(*elem);
      const double E_inc = std::sqrt(
          Rs / (elem->GetGeometryWidth() * elem->GetGeometryLength() * elems.size()));
      V_inc += E_inc * elem->GetGeometryLength() / elems.size();
    }
    return V_inc;
  }
  else
  {
    return 0.0;
  }
}

void LumpedPortData::InitializeLinearForms(mfem::ParFiniteElementSpace &nd_fespace) const
{
  const auto &mesh = *nd_fespace.GetParMesh();
  mfem::Array<int> attr_marker;
  if (!s || !v)
  {
    mfem::Array<int> attr_list;
    for (const auto &elem : elems)
    {
      attr_list.Append(elem->GetAttrList());
    }
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mesh::AttrToMarker(bdr_attr_max, attr_list, attr_marker);
  }

  // The port S-parameter, or the projection of the field onto the port mode, is computed
  // as: (E x H_inc) ⋅ n = E ⋅ (E_inc / Z_s), integrated over the port surface.
  if (!s)
  {
    SumVectorCoefficient fb(mesh.SpaceDimension());
    for (const auto &elem : elems)
    {
      const double Rs = R * GetToSquare(*elem);
      const double Hinc = (std::abs(Rs) > 0.0)
                              ? 1.0 / std::sqrt(Rs * elem->GetGeometryWidth() *
                                                elem->GetGeometryLength() * elems.size())
                              : 0.0;
      fb.AddCoefficient(elem->GetModeCoefficient(Hinc));
    }
    s = std::make_unique<mfem::LinearForm>(&nd_fespace);
    s->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb), attr_marker);
    s->UseFastAssembly(false);
    s->UseDevice(false);
    s->Assemble();
    s->UseDevice(true);
  }

  // The voltage across a port is computed using the electric field solution.
  // We have:
  //             V = ∫ E ⋅ l̂ dl = 1/w ∫ E ⋅ l̂ dS  (for rectangular ports)
  // or,
  //             V = 1/(2π) ∫ E ⋅ r̂ / r dS        (for coaxial ports).
  // We compute the surface integral via an inner product between the linear form with the
  // averaging function as a vector coefficient and the solution expansion coefficients.
  if (!v)
  {
    SumVectorCoefficient fb(mesh.SpaceDimension());
    for (const auto &elem : elems)
    {
      fb.AddCoefficient(
          elem->GetModeCoefficient(1.0 / (elem->GetGeometryWidth() * elems.size())));
    }
    v = std::make_unique<mfem::LinearForm>(&nd_fespace);
    v->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb), attr_marker);
    v->UseFastAssembly(false);
    v->UseDevice(false);
    v->Assemble();
    v->UseDevice(true);
  }
}

std::complex<double> LumpedPortData::GetPower(GridFunction &E, GridFunction &B) const
{
  // Compute port power, (E x H) ⋅ n = E ⋅ (-n x H), integrated over the port surface using
  // the computed E and H = μ⁻¹ B fields, where +n is the direction of propagation (into the
  // domain). The BdrSurfaceCurrentVectorCoefficient computes -n x H for an outward normal,
  // so we multiply by -1. The linear form is reconstructed from scratch each time due to
  // changing H.
  MFEM_VERIFY((E.HasImag() && B.HasImag()) || (!E.HasImag() && !B.HasImag()),
              "Mismatch between real- and complex-valued E and B fields in port power "
              "calculation!");
  const bool has_imag = E.HasImag();
  auto &nd_fespace = *E.ParFESpace();
  const auto &mesh = *nd_fespace.GetParMesh();
  SumVectorCoefficient fbr(mesh.SpaceDimension()), fbi(mesh.SpaceDimension());
  mfem::Array<int> attr_list;
  for (const auto &elem : elems)
  {
    fbr.AddCoefficient(
        std::make_unique<RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient>>(
            elem->GetAttrList(), B.Real(), mat_op));
    if (has_imag)
    {
      fbi.AddCoefficient(
          std::make_unique<RestrictedVectorCoefficient<BdrSurfaceCurrentVectorCoefficient>>(
              elem->GetAttrList(), B.Imag(), mat_op));
    }
    attr_list.Append(elem->GetAttrList());
  }
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  std::complex<double> dot;
  {
    mfem::LinearForm pr(&nd_fespace);
    pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbr), attr_marker);
    pr.UseFastAssembly(false);
    pr.UseDevice(false);
    pr.Assemble();
    pr.UseDevice(true);
    dot = -(pr * E.Real()) + (has_imag ? -1i * (pr * E.Imag()) : 0.0);
  }
  if (has_imag)
  {
    mfem::LinearForm pi(&nd_fespace);
    pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fbi), attr_marker);
    pi.UseFastAssembly(false);
    pi.UseDevice(false);
    pi.Assemble();
    pi.UseDevice(true);
    dot += -(pi * E.Imag()) + 1i * (pi * E.Real());
    Mpi::GlobalSum(1, &dot, E.ParFESpace()->GetComm());
    return dot;
  }
  else
  {
    double rdot = dot.real();
    Mpi::GlobalSum(1, &rdot, E.ParFESpace()->GetComm());
    return rdot;
  }
}

std::complex<double> LumpedPortData::GetSParameter(GridFunction &E) const
{
  // Compute port S-parameter, or the projection of the field onto the port mode.
  InitializeLinearForms(*E.ParFESpace());
  std::complex<double> dot((*s) * E.Real(), 0.0);
  if (E.HasImag())
  {
    dot.imag((*s) * E.Imag());
  }
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

std::complex<double> LumpedPortData::GetVoltage(GridFunction &E) const
{
  // Compute the average voltage across the port.
  InitializeLinearForms(*E.ParFESpace());
  std::complex<double> dot((*v) * E.Real(), 0.0);
  if (E.HasImag())
  {
    dot.imag((*v) * E.Imag());
  }
  Mpi::GlobalSum(1, &dot, E.GetComm());
  return dot;
}

LumpedPortOperator::LumpedPortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                       const mfem::ParMesh &mesh)
{
  // Set up lumped port boundary conditions.
  SetUpBoundaryProperties(iodata, mat_op, mesh);
  PrintBoundaryInfo(iodata, mesh);
}

void LumpedPortOperator::SetUpBoundaryProperties(const IoData &iodata,
                                                 const MaterialOperator &mat_op,
                                                 const mfem::ParMesh &mesh)
{
  // Check that lumped port boundary attributes have been specified correctly.
  if (!iodata.boundaries.lumpedport.empty())
  {
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), port_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    port_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : iodata.boundaries.lumpedport)
    {
      for (const auto &elem : data.elements)
      {
        for (auto attr : elem.attributes)
        {
          MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                      "Port boundary attribute tags must be non-negative and correspond to "
                      "boundaries in the mesh!");
          MFEM_VERIFY(bdr_attr_marker[attr - 1],
                      "Unknown port boundary attribute " << attr << "!");
          MFEM_VERIFY(!data.active || !port_marker[attr - 1],
                      "Boundary attribute is assigned to more than one lumped port!");
          port_marker[attr - 1] = 1;
        }
      }
    }
  }

  // Set up lumped port data structures.
  for (const auto &[idx, data] : iodata.boundaries.lumpedport)
  {
    ports.try_emplace(idx, data, mat_op, mesh);
  }
}

void LumpedPortOperator::PrintBoundaryInfo(const IoData &iodata, const mfem::ParMesh &mesh)
{
  if (ports.empty())
  {
    return;
  }

  fmt::memory_buffer buf{};  // Output buffer & buffer append lambda for cleaner code
  auto to = [](auto &buf, auto fmt, auto &&...args)
  { fmt::format_to(std::back_inserter(buf), fmt, std::forward<decltype(args)>(args)...); };
  using VT = IoData::ValueType;

  // Print out BC info for all port attributes, for both active and inactive ports.
  to(buf, "\nConfiguring Robin impedance BC for lumped ports at attributes:\n");
  for (const auto &[idx, data] : ports)
  {
    for (const auto &elem : data.elems)
    {
      for (auto attr : elem->GetAttrList())
      {
        to(buf, " {:d}:", attr);
        if (std::abs(data.R) > 0.0)
        {
          double Rs = data.R * data.GetToSquare(*elem);
          to(buf, " Rs = {:.3e} Ω/sq,", iodata.DimensionalizeValue(VT::IMPEDANCE, Rs));
        }
        if (std::abs(data.L) > 0.0)
        {
          double Ls = data.L * data.GetToSquare(*elem);
          to(buf, " Ls = {:.3e} H/sq,", iodata.DimensionalizeValue(VT::INDUCTANCE, Ls));
        }
        if (std::abs(data.C) > 0.0)
        {
          double Cs = data.C / data.GetToSquare(*elem);
          to(buf, " Cs = {:.3e} F/sq,", iodata.DimensionalizeValue(VT::CAPACITANCE, Cs));
        }
        to(buf, " n = ({:+.1f})\n", fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
      }
    }
  }

  // Print out port info for all active ports.
  fmt::memory_buffer buf_a{};
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    to(buf_a, " Index = {:d}: ", idx);
    if (std::abs(data.R) > 0.0)
    {
      to(buf_a, "R = {:.3e} Ω,", iodata.DimensionalizeValue(VT::IMPEDANCE, data.R));
    }
    if (std::abs(data.L) > 0.0)
    {
      to(buf_a, "L = {:.3e} H,", iodata.DimensionalizeValue(VT::INDUCTANCE, data.L));
    }
    if (std::abs(data.C) > 0.0)
    {
      to(buf_a, "C = {:.3e} F,", iodata.DimensionalizeValue(VT::CAPACITANCE, data.C));
    }
    buf_a.resize(buf_a.size() - 1);  // Remove last ","
    to(buf_a, "\n");
  }
  if (buf_a.size() > 0)
  {
    to(buf, "\nConfiguring lumped port circuit properties:\n");
    buf.append(buf_a);
    buf_a.clear();
  }

  // Print some information for excited lumped ports.
  for (const auto &[idx, data] : ports)
  {
    if (data.excitation == 0)
    {
      continue;
    }

    for (const auto &elem : data.elems)
    {
      for (auto attr : elem->GetAttrList())
      {
        to(buf_a, " {:d}: Index = {:d}\n", attr, idx);
      }
    }
  }
  if (buf_a.size() > 0)
  {
    to(buf, "\nConfiguring lumped port excitation source term at attributes:\n");
    buf.append(buf_a);
  }

  Mpi::Print("{}", fmt::to_string(buf));
}

const LumpedPortData &LumpedPortOperator::GetPort(int idx) const
{
  auto it = ports.find(idx);
  MFEM_VERIFY(it != ports.end(), "Unknown lumped port index requested!");
  return it->second;
}

mfem::Array<int> LumpedPortOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    for (const auto &elem : data.elems)
    {
      attr_list.Append(elem->GetAttrList());
    }
  }
  return attr_list;
}

mfem::Array<int> LumpedPortOperator::GetRsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.R) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        attr_list.Append(elem->GetAttrList());
      }
    }
  }
  return attr_list;
}

mfem::Array<int> LumpedPortOperator::GetLsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.L) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        attr_list.Append(elem->GetAttrList());
      }
    }
  }
  return attr_list;
}

mfem::Array<int> LumpedPortOperator::GetCsAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.C) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        attr_list.Append(elem->GetAttrList());
      }
    }
  }
  return attr_list;
}

void LumpedPortOperator::AddStiffnessBdrCoefficients(double coeff,
                                                     MaterialPropertyCoefficient &fb)
{
  // Add lumped inductor boundaries to the bilinear form.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.L) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        const double Ls = data.L * data.GetToSquare(*elem);
        fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(elem->GetAttrList()),
                               coeff / Ls);
      }
    }
  }
}

void LumpedPortOperator::AddDampingBdrCoefficients(double coeff,
                                                   MaterialPropertyCoefficient &fb)
{
  // Add lumped resistor boundaries to the bilinear form.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.R) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        const double Rs = data.R * data.GetToSquare(*elem);
        fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(elem->GetAttrList()),
                               coeff / Rs);
      }
    }
  }
}

void LumpedPortOperator::AddMassBdrCoefficients(double coeff,
                                                MaterialPropertyCoefficient &fb)
{
  // Add lumped capacitance boundaries to the bilinear form.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (std::abs(data.C) > 0.0)
    {
      for (const auto &elem : data.elems)
      {
        const double Cs = data.C / data.GetToSquare(*elem);
        fb.AddMaterialProperty(data.mat_op.GetCeedBdrAttributes(elem->GetAttrList()),
                               coeff * Cs);
      }
    }
  }
}

void LumpedPortOperator::AddExcitationBdrCoefficients(int excitation_idx,
                                                      SumVectorCoefficient &fb)
{
  // Construct the RHS source term for lumped port boundaries, which looks like -U_inc =
  // +2 iω/Z_s E_inc for a port boundary with an incident field E_inc. The chosen incident
  // field magnitude corresponds to a unit incident power over the full port boundary. See
  // p. 49 and p. 82 of the COMSOL RF Module manual for more detail.
  // Note: The real RHS returned here does not yet have the factor of (iω) included, so
  // works for time domain simulations requiring RHS -U_inc(t).
  for (const auto &[idx, data] : ports)
  {
    if (data.excitation != excitation_idx)
    {
      continue;
    }
    MFEM_VERIFY(std::abs(data.R) > 0.0,
                "Unexpected zero resistance in excited lumped port!");
    for (const auto &elem : data.elems)
    {
      const double Rs = data.R * data.GetToSquare(*elem);
      const double Hinc = 1.0 / std::sqrt(Rs * elem->GetGeometryWidth() *
                                          elem->GetGeometryLength() * data.elems.size());
      fb.AddCoefficient(elem->GetModeCoefficient(2.0 * Hinc));
    }
  }
}

}  // namespace palace
