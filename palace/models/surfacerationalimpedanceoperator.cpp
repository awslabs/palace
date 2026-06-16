// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "surfacerationalimpedanceoperator.hpp"

#include <complex>
#include <set>
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

namespace
{

// Horner evaluation of a polynomial with real coefficients stored highest-degree-first,
// i.e. coeffs = {c_n, ..., c_1, c_0} represents c_n s^n + ... + c_1 s + c_0.
std::complex<double> EvalPoly(const std::vector<double> &coeffs, std::complex<double> s)
{
  std::complex<double> val = 0.0;
  for (double c : coeffs)
  {
    val = val * s + c;
  }
  return val;
}

}  // namespace

SurfaceRationalImpedanceOperator::SurfaceRationalImpedanceOperator(
    const std::vector<config::RationalImpedanceData> &impedance,
    const std::unordered_set<int> &cracked_attributes, ProblemType problem_type,
    const Units &units, const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  SetUpBoundaryProperties(impedance, cracked_attributes, problem_type, mesh);
  PrintBoundaryInfo(units, mesh);
}

SurfaceRationalImpedanceOperator::SurfaceRationalImpedanceOperator(
    const IoData &iodata, const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : SurfaceRationalImpedanceOperator(iodata.boundaries.rational_impedance,
                                     iodata.boundaries.cracked_attributes,
                                     iodata.problem.type, iodata.units, mat_op, mesh)
{
}

void SurfaceRationalImpedanceOperator::SetUpBoundaryProperties(
    const std::vector<config::RationalImpedanceData> &impedance,
    const std::unordered_set<int> &cracked_attributes, ProblemType problem_type,
    const mfem::ParMesh &mesh)
{
  // Check that rational impedance boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!impedance.empty())
  {
    mfem::Array<int> impedance_marker(bdr_attr_max);
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    impedance_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (const auto &data : impedance)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(!impedance_marker[attr - 1],
                    "Multiple definitions of rational impedance boundary properties for "
                    "boundary attribute "
                        << attr << "!");
        impedance_marker[attr - 1] = 1;
        if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
        {
          bdr_warn_list.insert(attr);
        }
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning("Unknown rational impedance boundary attributes!\nSolver will just "
                   "ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
    }
  }

  // A rational surface impedance contributes only to the frequency-dependent system matrix
  // A2(ω); it has no time-domain (transient) or static realization here.
  MFEM_VERIFY(impedance.empty() || problem_type == ProblemType::DRIVEN ||
                  problem_type == ProblemType::EIGENMODE ||
                  problem_type == ProblemType::BOUNDARYMODE,
              "Rational impedance boundaries are only available for frequency-domain "
              "simulations (driven, eigenmode, boundary mode)!");

  boundaries.reserve(impedance.size());
  for (const auto &data : impedance)
  {
    MFEM_VERIFY(!data.num.empty() && !data.den.empty(),
                "Rational impedance boundary requires nonempty \"Numerator\" and "
                "\"Denominator\" coefficient lists!");
    auto &bdr = boundaries.emplace_back();
    bdr.num = data.num;
    bdr.den = data.den;
    bdr.attr_list.Reserve(static_cast<int>(data.attributes.size()));
    for (auto attr : data.attributes)
    {
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        continue;  // Can just ignore if wrong
      }
      bdr.attr_list.Append(attr);
      // Per-attribute scaling to account for increased area when using mesh cracking.
      bdr.attr_scaling[attr] =
          (cracked_attributes.find(attr) != cracked_attributes.end()) ? 2.0 : 1.0;
    }
  }
}

void SurfaceRationalImpedanceOperator::PrintBoundaryInfo(const Units &units,
                                                         const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin rational impedance BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      Mpi::Print(" {:d}: Zs(iω) = N(iω)/D(iω), N deg = {:d}, D deg = {:d}, n = ({:+.1f})\n",
                 attr, static_cast<int>(bdr.num.size()) - 1,
                 static_cast<int>(bdr.den.size()) - 1,
                 fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
    }
  }
}

mfem::Array<int> SurfaceRationalImpedanceOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    attr_list.Append(bdr.attr_list);
  }
  return attr_list;
}

void SurfaceRationalImpedanceOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi)
{
  // Rational surface impedance boundaries Zs(s) = N(s)/D(s) with s = iω. The Robin BC adds
  // the term iω / Zs(iω) per square to the frequency-dependent system matrix A2(ω), exactly
  // as the parallel admittance Ys = 1/Rs + 1/(iωLs) + iωCs contributes iω·Ys to
  // K + iωC - ω²M. Coefficients are stored nondimensionalized and highest-degree-first.
  if (omega == 0.0)
  {
    return;  // iω/Zs -> 0 for a finite-impedance boundary; frequency-domain drivers never
             // evaluate A2 at ω = 0.
  }
  const std::complex<double> s(0.0, omega);  // s = iω (nondimensional)
  for (const auto &bdr : boundaries)
  {
    const std::complex<double> N = EvalPoly(bdr.num, s);
    const std::complex<double> D = EvalPoly(bdr.den, s);
    MFEM_VERIFY(std::abs(N) > 0.0,
                "Rational impedance boundary has a transmission zero (Zs = 0) at the "
                "evaluation frequency; the admittance iω/Zs is singular!");
    const std::complex<double> Z = N / D;
    for (auto attr : bdr.attr_list)
    {
      const double sc = bdr.attr_scaling.at(attr);
      const std::complex<double> coef = s / (Z * sc);  // iω·Ys per square
      fbr.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coef.real());
      fbi.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coef.imag());
    }
  }
}

}  // namespace palace
