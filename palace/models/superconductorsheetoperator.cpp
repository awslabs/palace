// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "superconductorsheetoperator.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <set>
#include <vector>
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

std::unique_ptr<mfem::HypreParMatrix> SuperconductorSheetOperator::BuildTwoPortCoupling(
    mfem::ParFiniteElementSpace &nd_fespace) const
{
  if (!has_two_port_)
  {
    return nullptr;
  }
  mfem::ParMesh &pmesh = *nd_fespace.GetParMesh();
  MPI_Comm comm = nd_fespace.GetComm();
  const int sd = pmesh.SpaceDimension();

  // Cross-face coupling coefficient -1/(lambda*sinh(d/lambda)) per two-sided attribute.
  std::unordered_map<int, double> attr_off;
  for (const auto &bdr : boundaries)
  {
    if (!bdr.two_sided)
    {
      continue;
    }
    const double off = -1.0 / (bdr.lambda_L * std::sinh(bdr.thickness / bdr.lambda_L));
    for (auto attr : bdr.attr_list)
    {
      attr_off[attr] = off;
    }
  }

  // Reference triangle FE and quadrature, consistent across ranks (needed even where the
  // film is absent). Only triangular film faces are supported.
  const mfem::FiniteElement *tri_fe =
      nd_fespace.FEColl()->FiniteElementForGeometry(mfem::Geometry::TRIANGLE);
  MFEM_VERIFY(tri_fe, "Two-sided superconductor sheets require a triangular surface mesh!");
  const int nd = tri_fe->GetDof();
  const mfem::IntegrationRule &ir =
      mfem::IntRules.Get(mfem::Geometry::TRIANGLE, 2 * tri_fe->GetOrder() + 2);
  const int nqp = ir.GetNPoints();

  // Fixed-size per-face record so records concatenate for an Allgatherv. The two coincident
  // faces may be on different ranks, so we gather every film face and each rank assembles
  // the matrix rows it owns. Layout: off | 3*sd canonical vertex coords (sorted, identical
  // for the two faces) | nd global true-dofs | nd signs | nqp*nd*sd Piola basis | nqp
  // weights.
  const int o_vert = 1, o_gtd = o_vert + 3 * sd, o_sgn = o_gtd + nd, o_vsh = o_sgn + nd,
            o_w = o_vsh + nqp * nd * sd, rec = o_w + nqp;

  std::vector<double> local;
  mfem::DenseMatrix vsh(nd, sd);
  mfem::Array<int> vids, dofs;
  for (int be = 0; be < pmesh.GetNBE(); be++)
  {
    auto it = attr_off.find(pmesh.GetBdrAttribute(be));
    if (it == attr_off.end())
    {
      continue;
    }
    pmesh.GetBdrElementVertices(be, vids);
    MFEM_VERIFY(vids.Size() == 3, "Two-sided sheet faces must be triangles!");
    // Canonical vertex order (sorted lexicographically), identical for the two coincident
    // faces regardless of how cracking ordered their local vertices.
    int canon[3] = {0, 1, 2};
    auto less = [&](int p, int q)
    {
      const double *a = pmesh.GetVertex(vids[p]), *b = pmesh.GetVertex(vids[q]);
      for (int d = 0; d < sd; d++)
      {
        if (a[d] != b[d])
        {
          return a[d] < b[d];
        }
      }
      return false;
    };
    std::sort(canon, canon + 3, less);

    const mfem::FiniteElement *fe = nd_fespace.GetBE(be);
    mfem::ElementTransformation *T = pmesh.GetBdrElementTransformation(be);
    nd_fespace.GetBdrElementDofs(be, dofs);

    const std::size_t base = local.size();
    local.resize(base + rec, 0.0);
    local[base] = it->second;  // off
    for (int k = 0; k < 3; k++)
    {
      const double *x = pmesh.GetVertex(vids[canon[k]]);
      for (int d = 0; d < sd; d++)
      {
        local[base + o_vert + k * sd + d] = x[d];
      }
    }
    for (int a = 0; a < nd; a++)
    {
      int ld = dofs[a];
      double sgn = 1.0;
      if (ld < 0)
      {
        ld = -1 - ld;
        sgn = -1.0;
      }
      local[base + o_gtd + a] = static_cast<double>(nd_fespace.GetGlobalTDofNumber(ld));
      local[base + o_sgn + a] = sgn;
    }
    for (int q = 0; q < nqp; q++)
    {
      // Barycentric (b0,b1,b2) on canonical vertices -> element reference coords. The
      // physical quadrature point is then identical for the two coincident faces.
      const mfem::IntegrationPoint &ip = ir.IntPoint(q);
      const double bc[3] = {1.0 - ip.x - ip.y, ip.x, ip.y};
      double eb[3];
      for (int k = 0; k < 3; k++)
      {
        eb[canon[k]] = bc[k];
      }
      mfem::IntegrationPoint ipe;
      ipe.Set2(eb[1], eb[2]);
      T->SetIntPoint(&ipe);
      fe->CalcVShape(*T, vsh);
      local[base + o_w + q] = ip.weight * T->Weight();
      for (int a = 0; a < nd; a++)
      {
        for (int d = 0; d < sd; d++)
        {
          local[base + o_vsh + (q * nd + a) * sd + d] = vsh(a, d);
        }
      }
    }
  }

  // Gather all film-face records across ranks.
  int sendcount = static_cast<int>(local.size()), nranks;
  MPI_Comm_size(comm, &nranks);
  std::vector<int> counts(nranks), displs(nranks);
  Mpi::Allgather(1, &sendcount, counts.data(), comm);
  int total = 0;
  for (int r = 0; r < nranks; r++)
  {
    displs[r] = total;
    total += counts[r];
  }
  std::vector<double> all(total);
  Mpi::Allgatherv(sendcount, local.data(), all.data(), counts.data(), displs.data(), comm);

  // Group faces by canonical vertex key; each physical film triangle appears exactly twice.
  std::map<std::array<long, 9>, std::vector<int>> by_key;  // key -> record start offsets
  for (int off = 0; off + rec <= total; off += rec)
  {
    std::array<long, 9> key{};
    for (int i = 0; i < 3 * sd && i < 9; i++)
    {
      key[i] = std::lround(all[off + o_vert + i] * 1e6);
    }
    by_key[key].push_back(off);
  }

  // Assemble C(row, col) = off * sign_p[a] * sign_m[b] * <phi_a^p, phi_b^m> into the rows
  // this rank owns; columns use global true-dof numbers.
  const HYPRE_BigInt my_off = nd_fespace.GetMyTDofOffset();
  const int lt = nd_fespace.GetTrueVSize();
  const HYPRE_BigInt glob = nd_fespace.GlobalTrueVSize();
  mfem::SparseMatrix Cloc(lt, static_cast<int>(glob));
  int n_unpaired = 0;
  for (auto &[key, recs] : by_key)
  {
    if (recs.size() != 2)
    {
      n_unpaired += static_cast<int>(recs.size());
      continue;
    }
    const int P = recs[0], M = recs[1];
    const double off = all[P];
    auto add_block = [&](int R, int S)  // rows from face R, columns from face S
    {
      for (int a = 0; a < nd; a++)
      {
        const HYPRE_BigInt row = static_cast<HYPRE_BigInt>(all[R + o_gtd + a]);
        if (row < my_off || row >= my_off + lt)
        {
          continue;
        }
        const double sa = all[R + o_sgn + a];
        for (int b = 0; b < nd; b++)
        {
          double crossab = 0.0;
          for (int q = 0; q < nqp; q++)
          {
            double dot = 0.0;
            for (int d = 0; d < sd; d++)
            {
              dot += all[R + o_vsh + (q * nd + a) * sd + d] *
                     all[S + o_vsh + (q * nd + b) * sd + d];
            }
            crossab += all[R + o_w + q] * dot;
          }
          const double val = off * sa * all[S + o_sgn + b] * crossab;
          Cloc.Add(static_cast<int>(row - my_off), static_cast<int>(all[S + o_gtd + b]),
                   val);
        }
      }
    };
    add_block(P, M);
    add_block(M, P);
  }
  Cloc.Finalize();
  Mpi::GlobalSum(1, &n_unpaired, comm);
  MFEM_VERIFY(
      n_unpaired == 0,
      "Two-sided (two-port) sheet has "
          << n_unpaired
          << " unpaired face(s): each cracked film face must have a coincident twin "
             "(the film must be fully interior).");

  // Build the parallel matrix from local rows with global column indices.
  std::vector<HYPRE_BigInt> J(Cloc.NumNonZeroElems());
  for (int k = 0; k < Cloc.NumNonZeroElems(); k++)
  {
    J[k] = Cloc.GetJ()[k];
  }
  HYPRE_BigInt *rows = nd_fespace.GetTrueDofOffsets();
  auto C = std::make_unique<mfem::HypreParMatrix>(comm, lt, glob, glob, Cloc.GetI(),
                                                  J.data(), Cloc.GetData(), rows, rows);
  C->CopyRowStarts();
  C->CopyColStarts();
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
