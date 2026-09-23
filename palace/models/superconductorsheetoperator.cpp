// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superconductorsheetoperator.hpp"

#include <array>
#include <cmath>
#include <map>
#include <set>
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

SuperconductorSheetOperator::SuperconductorSheetOperator(
    const std::vector<config::SuperconductorData> &superconductor,
    const std::unordered_set<int> &cracked_attributes, const Units &units,
    const MaterialOperator &mat_op, const mfem::ParMesh &mesh)
  : mat_op(mat_op)
{
  SetUpBoundaryProperties(superconductor, cracked_attributes, mesh);
  PrintBoundaryInfo(units, mesh);
}

SuperconductorSheetOperator::SuperconductorSheetOperator(const IoData &iodata,
                                                         const MaterialOperator &mat_op,
                                                         const mfem::ParMesh &mesh)
  : SuperconductorSheetOperator(iodata.boundaries.superconductor,
                                iodata.boundaries.cracked_attributes, iodata.units, mat_op,
                                mesh)
{
}

void SuperconductorSheetOperator::SetUpBoundaryProperties(
    const std::vector<config::SuperconductorData> &superconductor,
    const std::unordered_set<int> &cracked_attributes, const mfem::ParMesh &mesh)
{
  // Check that superconductor sheet boundary attributes have been specified correctly.
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> bdr_attr_marker;
  if (!superconductor.empty())
  {
    mfem::Array<int> superconductor_marker(bdr_attr_max);
    bdr_attr_marker.SetSize(bdr_attr_max);
    bdr_attr_marker = 0;
    superconductor_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    std::set<int> bdr_warn_list;
    for (const auto &data : superconductor)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(!superconductor_marker[attr - 1],
                    "Multiple definitions of superconductor sheet boundary properties for "
                    "boundary attribute "
                        << attr << "!");
        superconductor_marker[attr - 1] = 1;
        if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
        {
          bdr_warn_list.insert(attr);
        }
      }
    }
    if (!bdr_warn_list.empty())
    {
      Mpi::Print("\n");
      Mpi::Warning("Unknown superconductor sheet boundary attributes!\nSolver will just "
                   "ignore them!");
      utils::PrettyPrint(bdr_warn_list, "Boundary attribute list:");
      Mpi::Print("\n");
    }
  }

  // Kinetic sheet inductance L_ksq [H/sq]: supplied directly, or from (lambda, d) via the
  // finite-thickness London sheet L_ksq = lambda*coth(d/lambda) (nondimensional, mu0
  // absorbed). Reduces to the thin-film Pearl limit lambda^2/d for d << lambda and
  // saturates at lambda for d >> lambda, capturing the through-thickness current profile
  // without resolving the film thickness.
  boundaries.reserve(superconductor.size());
  for (const auto &data : superconductor)
  {
    // Two-sided (two-port) sheet: the diagonal self-term inductance is
    // lambda*tanh(d/lambda), so 1/Ls = coth(d/lambda)/lambda = cosh/(lambda*sinh); the
    // off-diagonal cross-face term is added separately via BuildTwoPortCoupling. Otherwise
    // the single-sheet L_ksq applies.
    double Ls;
    if (data.two_sided)
    {
      Ls = data.lambda_L * std::tanh(data.thickness / data.lambda_L);
      has_two_port_ = true;
    }
    else
    {
      Ls =
          (data.Ls > 0.0) ? data.Ls : KineticSheetInductance(data.lambda_L, data.thickness);
    }
    MFEM_VERIFY(Ls > 0.0,
                "Superconductor sheet has non-positive kinetic sheet inductance!");
    auto &bdr = boundaries.emplace_back();
    bdr.Ls = Ls;
    bdr.lambda_L = data.lambda_L;
    bdr.thickness = data.thickness;
    bdr.two_sided = data.two_sided;
    bdr.attr_list.Reserve(static_cast<int>(data.attributes.size()));
    for (auto attr : data.attributes)
    {
      if (attr <= 0 || attr > bdr_attr_max || !bdr_attr_marker[attr - 1])
      {
        continue;  // Can just ignore if wrong
      }
      bdr.attr_list.Append(attr);
      // Per-attribute scaling for the doubled area of a cracked interface. A two-sided
      // sheet keeps each face at its full self-term (the two faces are physically
      // distinct), so it is not halved even though it cracks.
      bdr.attr_scaling[attr] =
          (!data.two_sided && cracked_attributes.find(attr) != cracked_attributes.end())
              ? 2.0
              : 1.0;
    }
  }
}

void SuperconductorSheetOperator::PrintBoundaryInfo(const Units &units,
                                                    const mfem::ParMesh &mesh)
{
  if (boundaries.empty())
  {
    return;
  }

  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};
  using VT = Units::ValueType;

  fmt::format_to(out, "\nConfiguring superconductor sheet BC at attributes:\n");
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      fmt::format_to(out, " {:d}: L_ksq = {:.3e} H/sq, n = ({:+.1f})\n", attr,
                     units.Dimensionalize<VT::INDUCTANCE>(bdr.Ls),
                     fmt::join(mesh::GetSurfaceNormal(mesh, attr), ","));
    }
  }
  Mpi::Print("{}", fmt::to_string(buffer));
}

double SuperconductorSheetOperator::KineticSheetInductance(double lambda_L,
                                                           double thickness)
{
  // Finite-thickness London sheet inductance L_ksq = lambda * coth(d/lambda), with
  // coth(x) = 1/tanh(x). Reduces to lambda^2/d for d << lambda; saturates at lambda for
  // d >> lambda.
  return lambda_L / std::tanh(thickness / lambda_L);
}

std::unique_ptr<mfem::SparseMatrix> SuperconductorSheetOperator::BuildTwoPortCoupling(
    mfem::ParFiniteElementSpace &nd_fespace) const
{
  if (!has_two_port_)
  {
    return nullptr;
  }
  mfem::ParMesh &pmesh = *nd_fespace.GetParMesh();
  const int sd = pmesh.SpaceDimension();
  auto C =
      std::make_unique<mfem::SparseMatrix>(nd_fespace.GetVSize(), nd_fespace.GetVSize());

  for (const auto &bdr : boundaries)
  {
    if (!bdr.two_sided)
    {
      continue;
    }
    // Cross-face coupling coefficient -1/(lambda*sinh(d/lambda)).
    const double off = -1.0 / (bdr.lambda_L * std::sinh(bdr.thickness / bdr.lambda_L));
    for (auto attr : bdr.attr_list)
    {
      // Group this attribute's boundary elements by centroid: the two coincident cracked
      // faces.
      std::map<std::array<long, 3>, std::vector<int>> by_centroid;
      for (int be = 0; be < pmesh.GetNBE(); be++)
      {
        if (pmesh.GetBdrAttribute(be) != attr)
        {
          continue;
        }
        mfem::Array<int> v;
        pmesh.GetBdrElementVertices(be, v);
        double cc[3] = {0.0, 0.0, 0.0};
        for (int k = 0; k < v.Size(); k++)
        {
          const double *x = pmesh.GetVertex(v[k]);
          cc[0] += x[0];
          cc[1] += x[1];
          cc[2] += x[2];
        }
        const double inv = 1.0 / v.Size(), q = 1e6;
        by_centroid[{std::lround(cc[0] * inv * q), std::lround(cc[1] * inv * q),
                     std::lround(cc[2] * inv * q)}]
            .push_back(be);
      }
      mfem::DenseMatrix vshp, vshm, cross;
      for (auto &[key, bes] : by_centroid)
      {
        if (bes.size() != 2)
        {
          continue;  // Only genuine coincident pairs couple.
        }
        const int bp = bes[0], bm = bes[1];
        const mfem::FiniteElement *fep = nd_fespace.GetBE(bp), *fem = nd_fespace.GetBE(bm);
        mfem::ElementTransformation *Tp = pmesh.GetBdrElementTransformation(bp);
        mfem::ElementTransformation *Tm = pmesh.GetBdrElementTransformation(bm);
        const int nd = fep->GetDof();
        mfem::Array<int> dp, dm;
        nd_fespace.GetBdrElementDofs(bp, dp);
        nd_fespace.GetBdrElementDofs(bm, dm);
        cross.SetSize(nd);
        cross = 0.0;
        vshp.SetSize(nd, sd);
        vshm.SetSize(nd, sd);
        // Map face-p reference points to face-m by the triangle vertex permutation (the
        // faces are geometrically coincident, but cracking may reorder their local
        // vertices). This is exact and avoids ElementTransformation::TransformBack, which
        // is unreliable for a surface element (2D reference in 3D space).
        mfem::Array<int> vp, vm;
        pmesh.GetBdrElementVertices(bp, vp);
        pmesh.GetBdrElementVertices(bm, vm);
        int perm[3] = {0, 1, 2};
        for (int i = 0; i < vp.Size(); i++)
        {
          const double *xi = pmesh.GetVertex(vp[i]);
          for (int j = 0; j < vm.Size(); j++)
          {
            const double *xj = pmesh.GetVertex(vm[j]);
            double d2 = 0.0;
            for (int dd = 0; dd < sd; dd++)
            {
              d2 += (xi[dd] - xj[dd]) * (xi[dd] - xj[dd]);
            }
            if (d2 < 1.0e-18)
            {
              perm[i] = j;
              break;
            }
          }
        }
        const mfem::IntegrationRule &ir =
            mfem::IntRules.Get(fep->GetGeomType(), 2 * fep->GetOrder() + 2);
        for (int iq = 0; iq < ir.GetNPoints(); iq++)
        {
          const mfem::IntegrationPoint &ip = ir.IntPoint(iq);
          Tp->SetIntPoint(&ip);
          const double w = ip.weight * Tp->Weight();
          fep->CalcVShape(*Tp, vshp);
          // Permute barycentric coordinates (v0,v1,v2) -> (perm) to land on the same
          // physical point in face m.
          double bary[3] = {1.0 - ip.x - ip.y, ip.x, ip.y}, barym[3] = {0.0, 0.0, 0.0};
          barym[perm[0]] = bary[0];
          barym[perm[1]] = bary[1];
          barym[perm[2]] = bary[2];
          mfem::IntegrationPoint ipm;
          ipm.Set2(barym[1], barym[2]);
          Tm->SetIntPoint(&ipm);
          fem->CalcVShape(*Tm, vshm);
          for (int a = 0; a < nd; a++)
          {
            for (int b = 0; b < nd; b++)
            {
              double dot = 0.0;
              for (int dd = 0; dd < sd; dd++)
              {
                dot += vshp(a, dd) * vshm(b, dd);
              }
              cross(a, b) += w * dot;
            }
          }
        }
        for (int a = 0; a < nd; a++)
        {
          int ga = dp[a];
          double sa = 1.0;
          if (ga < 0)
          {
            ga = -1 - ga;
            sa = -1.0;
          }
          for (int b = 0; b < nd; b++)
          {
            int gb = dm[b];
            double sb = 1.0;
            if (gb < 0)
            {
              gb = -1 - gb;
              sb = -1.0;
            }
            const double val = off * sa * sb * cross(a, b);
            C->Add(ga, gb, val);
            C->Add(gb, ga, val);  // symmetric partner <phi_b^-, phi_a^+>
          }
        }
      }
    }
  }
  C->Finalize();
  return C;
}

mfem::Array<int> SuperconductorSheetOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &bdr : boundaries)
  {
    attr_list.Append(bdr.attr_list);
  }
  return attr_list;
}

void SuperconductorSheetOperator::AddStiffnessBdrCoefficients(
    double coeff, MaterialPropertyCoefficient &fb) const
{
  // Kinetic sheet inductance boundaries: add (coeff / L_ksq) as a tangential surface mass.
  for (const auto &bdr : boundaries)
  {
    for (auto attr : bdr.attr_list)
    {
      const double s = bdr.attr_scaling.at(attr);
      fb.AddMaterialProperty(mat_op.GetCeedBdrAttributes(attr), coeff / (bdr.Ls * s));
    }
  }
}

}  // namespace palace
