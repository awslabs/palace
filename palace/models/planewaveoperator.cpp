// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "planewaveoperator.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <span>
#include <string>
#include <fmt/ranges.h>
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

namespace
{

// The local dof of a signed vector dof index.
inline int LocalDof(int d)
{
  return (d >= 0) ? d : -1 - d;
}

}  // namespace

PlaneWaveOperator::PlaneWaveOperator(const std::map<int, config::PlaneWaveData> &plane_wave,
                                     const Units &units, const mfem::ParMesh &mesh,
                                     bool print)
  : print(print)
{
  MFEM_VERIFY(plane_wave.empty() || mesh.SpaceDimension() == 3,
              "Plane wave excitations require a 3D mesh!");
  mesh_scale = units.GetMeshLengthRelativeScale();
  unit_GHz =
      units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) / (2.0 * std::numbers::pi);
  const int attr_max = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
  for (const auto &[idx, data] : plane_wave)
  {
    PlaneWaveSource src;
    src.direction.SetSize(3);
    src.polarization.SetSize(3);
    src.origin.SetSize(3);
    for (int d = 0; d < 3; d++)
    {
      src.direction(d) = data.direction[d];
      src.polarization(d) = data.polarization[d];
      src.origin(d) = data.origin[d];
    }
    const double kn = src.direction.Norml2();
    MFEM_VERIFY(kn > 0.0, "Plane wave direction must be nonzero!");
    src.direction /= kn;
    src.polarization.Add(-(src.polarization * src.direction), src.direction);
    const double pn = src.polarization.Norml2();
    MFEM_VERIFY(pn > 1.0e-12,
                "Plane wave polarization must not be parallel to the direction!");
    src.polarization /= pn;
    src.amplitude = units.Nondimensionalize<Units::ValueType::FIELD_E>(data.amplitude);
    src.attr_marker.SetSize(attr_max);
    src.attr_marker = 0;
    for (auto attr : data.attributes)
    {
      MFEM_VERIFY(attr > 0 && attr <= attr_max,
                  "Plane wave total-field attribute " << attr << " is not in the mesh!");
      src.attr_marker[attr - 1] = 1;
    }
    sources.emplace(idx, std::move(src));
  }

  if (sources.empty() || !print)
  {
    return;
  }
  Mpi::Print("\nConfiguring incident plane wave (total-field/scattered-field) sources:\n");
  for (const auto &[idx, src] : sources)
  {
    Mpi::Print(" Plane wave {:d}: amplitude = {:.3e} V/m\n"
               " \tDirection = ({:.3e})\n"
               " \tPolarization = ({:.3e})\n"
               " \tTotal-field attributes: {}\n",
               idx, plane_wave.at(idx).amplitude,
               fmt::join(std::span(src.direction.GetData(), 3), ", "),
               fmt::join(std::span(src.polarization.GetData(), 3), ", "),
               fmt::join(plane_wave.at(idx).attributes, " "));
  }
}

void PlaneWaveOperator::GetLocalTotalFieldMarker(int idx,
                                                 const mfem::ParFiniteElementSpace &fes,
                                                 mfem::Array<int> &lmark) const
{
  const auto &src = sources.at(idx);
  const mfem::ParMesh &mesh = *fes.GetParMesh();
  lmark.SetSize(fes.GetVSize());
  lmark = 0;
  mfem::Array<int> vdofs;
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    if (!src.attr_marker[mesh.GetAttribute(e) - 1])
    {
      continue;
    }
    fes.GetElementVDofs(e, vdofs);
    for (int d : vdofs)
    {
      lmark[d >= 0 ? d : -1 - d] = 1;
    }
  }
  fes.Synchronize(lmark);
}

const Vector &
PlaneWaveOperator::GetTotalFieldMarker(int idx, const FiniteElementSpace &nd_fespace) const
{
  if (auto it = tf_markers.find(idx);
      it != tf_markers.end() && it->second.Size() == nd_fespace.GetTrueVSize())
  {
    return it->second;
  }
  // Plane waves with the same total-field region share the marker.
  for (const auto &[other, tmark] : tf_markers)
  {
    if (other != idx && tmark.Size() == nd_fespace.GetTrueVSize() &&
        sources.at(other).attr_marker == sources.at(idx).attr_marker)
    {
      return tmark;
    }
  }
  // Mark the local dofs of the elements in the total-field region, then combine the marks
  // of the processes sharing each dof with a boolean OR and restrict to the true dofs, as
  // for essential dofs: a true dof belongs to the closure of the region if any element of
  // the region contains it. P^T is not used, since for Nédélec spaces with face dofs on
  // shared triangles it contains the face orientation transformations, so that it does
  // not map a 0/1 marker to a 0/1 marker.
  const mfem::ParFiniteElementSpace &fes = nd_fespace.Get();
  mfem::Array<int> lmark, tmark_bool;
  GetLocalTotalFieldMarker(idx, fes, lmark);
  fes.GetRestrictionMatrix()->BooleanMult(lmark, tmark_bool);
  Vector &tmark = tf_markers[idx];
  tmark.SetSize(fes.GetTrueVSize());
  const int *tb = tmark_bool.HostRead();
  auto *m = tmark.HostWrite();
  int n_tf = 0;
  for (int i = 0; i < tmark.Size(); i++)
  {
    m[i] = tb[i] ? 1.0 : 0.0;
    n_tf += (tb[i] != 0);
  }
  Mpi::GlobalSum(1, &n_tf, fes.GetComm());
  if (print)
  {
    Mpi::Print(" Plane wave {:d}: {:d} total-field unknowns (of {:d})\n", idx, n_tf,
               fes.GlobalTrueVSize());
  }
  return tmark;
}

bool PlaneWaveOperator::GetIncidentField(int idx, double omega,
                                         const FiniteElementSpace &nd_fespace,
                                         ComplexVector &g) const
{
  auto it = sources.find(idx);
  if (it == sources.end())
  {
    return false;
  }
  const auto &src = it->second;

  // Interpolate E_inc = E0 p exp(-i k khat . (x - x0)), with k = omega / c in the
  // background (nondimensional units), times exp(+i k_F . x) for the periodic envelope
  // which is the unknown with a Floquet wave vector k_F.
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(nd_fespace.Get());
  auto &mesh = *fes.GetParMesh();
  mfem::Vector kf(3);
  kf = 0.0;
  if (src.kf.Size() == 3)
  {
    kf = src.kf;
    if (src.kf_scaled)
    {
      kf *= omega;
    }
  }
  const int n = fes.GetTrueVSize();
  g.SetSize(n);
  g.UseDevice(true);
  double jump = 0.0, scale = 0.0;
  mfem::Array<int> vdofs;
  mfem::DofTransformation doftrans;
  mfem::Vector elvect;
  for (int part = 0; part < 2; part++)
  {
    mfem::VectorFunctionCoefficient coef(
        3,
        [&src, &kf, omega, part](const mfem::Vector &x, mfem::Vector &E)
        {
          double phase = 0.0, phase_f = 0.0;
          for (int d = 0; d < 3; d++)
          {
            phase += src.direction(d) * (x(d) - src.origin(d));
            phase_f += kf(d) * x(d);
          }
          phase = phase * omega / src.c - phase_f;
          const double f =
              src.amplitude * ((part == 0) ? std::cos(phase) : -std::sin(phase));
          E.SetSize(3);
          for (int d = 0; d < 3; d++)
          {
            E(d) = f * src.polarization(d);
          }
        });
    mfem::ParGridFunction gf(&fes);
    gf.ProjectCoefficient(coef);
    Vector &gt = (part == 0) ? g.Real() : g.Imag();
    gf.ParallelProject(gt);
    if (!check_interpolated)
    {
      continue;
    }

    // The interpolated field must be the same from every element sharing a dof: it is,
    // up to rounding, except across the periodic boundaries when the envelope is not
    // periodic, that is when the incident field is not quasi-periodic with k_F. Compare
    // the values of each element with the stored ones (the last element wins) on each
    // process, and with those of the owning process on the shared dofs.
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      fes.GetElementVDofs(e, vdofs, doftrans);
      elvect.SetSize(vdofs.Size());
      fes.GetFE(e)->Project(coef, *fes.GetElementTransformation(e), elvect);
      doftrans.TransformPrimal(elvect);
      for (int i = 0; i < vdofs.Size(); i++)
      {
        const int d = (vdofs[i] >= 0) ? vdofs[i] : -1 - vdofs[i];
        const double v = (vdofs[i] >= 0) ? elvect(i) : -elvect(i);
        jump = std::max(jump, std::abs(v - gf(d)));
        scale = std::max(scale, std::abs(v));
      }
    }
    mfem::ParGridFunction owner(&fes);
    owner.Distribute(gt);
    for (int d = 0; d < gf.Size(); d++)
    {
      jump = std::max(jump, std::abs(owner(d) - gf(d)));
    }
  }
  if (!check_interpolated)
  {
    return true;
  }
  Mpi::GlobalMax(1, &jump, fes.GetComm());
  Mpi::GlobalMax(1, &scale, fes.GetComm());
  const double rel = (scale > 0.0) ? jump / scale : 0.0;
  if (rel > 1.0e-6)
  {
    // Wave vectors in radians per mesh length unit, as the Floquet wave vector in the
    // configuration file.
    mfem::Vector k(src.direction);
    k *= omega / src.c / mesh_scale;
    kf *= 1.0 / mesh_scale;
    const auto msg = fmt::format(
        "Plane wave {:d}: the incident field is not quasi-periodic with the Floquet wave "
        "vector of the periodic boundaries (relative jump {:.2e} across them). The "
        "Floquet wave vector [{:.6e}, {:.6e}, {:.6e}] must be the part of the wave vector "
        "of the plane wave [{:.6e}, {:.6e}, {:.6e}] along the periodic directions, up to a "
        "reciprocal lattice vector",
        idx, rel, kf(0), kf(1), kf(2), k(0), k(1), k(2));
    MFEM_VERIFY(rel <= 1.0e-2, msg << "!");
    if (warned_interpolated.insert(idx).second)
    {
      Mpi::Warning(fes.GetComm(),
                   "{}; the source leaks into the scattered field at the level "
                   "of this jump (warned once, at the first frequency where it occurs).",
                   msg);
    }
  }
  return true;
}

void PlaneWaveOperator::GetIncidentField(int idx, double omega,
                                         const FiniteElementSpace &nd_fespace,
                                         const mfem::Array<int> &elements,
                                         ComplexVector &g) const
{
  const auto &src = sources.at(idx);
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(nd_fespace.Get());
  mfem::Vector kf(3);
  kf = 0.0;
  if (src.kf.Size() == 3)
  {
    kf = src.kf;
    if (src.kf_scaled)
    {
      kf *= omega;
    }
  }
  g.SetSize(fes.GetTrueVSize());
  g.UseDevice(true);
  mfem::Array<int> vdofs;
  mfem::DofTransformation doftrans;
  mfem::Vector elvect;
  mfem::ParGridFunction gf(&fes);
  for (int part = 0; part < 2; part++)
  {
    // The same interpolation as ParGridFunction::ProjectCoefficient, on the given elements.
    mfem::VectorFunctionCoefficient coef(
        3,
        [&src, &kf, omega, part](const mfem::Vector &x, mfem::Vector &E)
        {
          double phase = 0.0, phase_f = 0.0;
          for (int d = 0; d < 3; d++)
          {
            phase += src.direction(d) * (x(d) - src.origin(d));
            phase_f += kf(d) * x(d);
          }
          phase = phase * omega / src.c - phase_f;
          const double f =
              src.amplitude * ((part == 0) ? std::cos(phase) : -std::sin(phase));
          E.SetSize(3);
          for (int d = 0; d < 3; d++)
          {
            E(d) = f * src.polarization(d);
          }
        });
    gf = 0.0;
    for (int e : elements)
    {
      fes.GetElementVDofs(e, vdofs, doftrans);
      elvect.SetSize(vdofs.Size());
      fes.GetFE(e)->Project(coef, *fes.GetElementTransformation(e), elvect);
      doftrans.TransformPrimal(elvect);
      gf.SetSubVector(vdofs, elvect);
    }
    gf.ParallelProject((part == 0) ? g.Real() : g.Imag());
  }
}

void PlaneWaveOperator::GetCoupling(int idx, const FiniteElementSpace &nd_fespace,
                                    mfem::Array<int> &l_c,
                                    mfem::Array<int> *coupling_elements) const
{
  // The system matrix couples the true dofs of the closure of the total-field region (the
  // marker m) and the others through the elements whose local dofs depend on both, with the
  // prolongation P from the true dofs to the local dofs: on a nonconforming mesh, also
  // through the constraints of hanging faces and edges. The local dofs which depend on
  // marked or on unmarked true dofs are |P| m and |P| (1 - m).
  const auto &src = sources.at(idx);
  const mfem::ParFiniteElementSpace &fes = nd_fespace.Get();
  const mfem::ParMesh &mesh = *fes.GetParMesh();
  mfem::HypreParMatrix &P = *fes.Dof_TrueDof_Matrix();
  const Vector &tmark = GetTotalFieldMarker(idx, nd_fespace);
  const int nt = fes.GetTrueVSize(), nl = fes.GetVSize();
  mfem::Array<int> m_t(nt), m_o(nt), l_t(nl), l_o(nl), vdofs;
  {
    const auto *m = tmark.HostRead();
    for (int i = 0; i < nt; i++)
    {
      m_t[i] = (m[i] > 0.0);
      m_o[i] = !m_t[i];
    }
  }
  P.BooleanMult(1, m_t.HostRead(), 0, l_t.HostWrite());
  P.BooleanMult(1, m_o.HostRead(), 0, l_o.HostWrite());
  l_c.SetSize(nl);
  l_c = 0;
  if (coupling_elements)
  {
    coupling_elements->SetSize(0);
  }
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    fes.GetElementVDofs(e, vdofs);
    bool has_t = false, has_o = false;
    for (int d : vdofs)
    {
      has_t = has_t || l_t[LocalDof(d)];
      has_o = has_o || l_o[LocalDof(d)];
    }
    if (has_t && has_o)
    {
      for (int d : vdofs)
      {
        l_c[LocalDof(d)] = 1;
      }
      if (coupling_elements && !src.attr_marker[mesh.GetAttribute(e) - 1])
      {
        coupling_elements->Append(e);
      }
    }
  }
}

void PlaneWaveOperator::GetLayer(int idx, const FiniteElementSpace &nd_fespace,
                                 mfem::Array<int> &tdofs, mfem::Array<int> &elements,
                                 mfem::Array<int> *coupling_elements) const
{
  // The layer true dofs are those which the local dofs of the coupling elements depend on,
  // |P|ᵀ, and the incident field is interpolated on the elements with a local dof which
  // depends on a layer true dof.
  const mfem::ParFiniteElementSpace &fes = nd_fespace.Get();
  const mfem::ParMesh &mesh = *fes.GetParMesh();
  mfem::HypreParMatrix &P = *fes.Dof_TrueDof_Matrix();
  const int nt = fes.GetTrueVSize(), nl = fes.GetVSize();
  mfem::Array<int> l_c, t_c(nt), l_i(nl), vdofs;
  GetCoupling(idx, nd_fespace, l_c, coupling_elements);
  P.BooleanMultTranspose(1, l_c.HostRead(), 0, t_c.HostWrite());
  P.BooleanMult(1, t_c.HostRead(), 0, l_i.HostWrite());
  tdofs.SetSize(0);
  for (int i = 0; i < nt; i++)
  {
    if (t_c[i])
    {
      tdofs.Append(i);
    }
  }
  elements.SetSize(0);
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    fes.GetElementVDofs(e, vdofs);
    for (int d : vdofs)
    {
      if (l_i[LocalDof(d)])
      {
        elements.Append(e);
        break;
      }
    }
  }
}

void PlaneWaveOperator::GetVertexNeighbors(int idx, mfem::ParMesh &mesh,
                                           mfem::Array<int> &elements,
                                           mfem::Array<int> &bdr_elements,
                                           bool inside) const
{
  // Mark the vertices of the elements of the total-field region (or of the rest of the
  // domain, inside = true) on an order-1 H1 space, and extend the marker through the
  // prolongation to the vertices shared with other processes and to the hanging vertices
  // which depend on them: |P| |P|ᵀ.
  const auto &src = sources.at(idx);
  mfem::H1_FECollection fec(1, mesh.Dimension());
  mfem::ParFiniteElementSpace fes(&mesh, &fec);
  mfem::HypreParMatrix &P = *fes.Dof_TrueDof_Matrix();
  mfem::Array<int> l_v(fes.GetVSize()), t_v(fes.GetTrueVSize()), dofs;
  l_v = 0;
  auto InSource = [&](int e)
  { return src.attr_marker[mesh.GetAttribute(e) - 1] != inside; };
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    if (InSource(e))
    {
      fes.GetElementDofs(e, dofs);
      for (int d : dofs)
      {
        l_v[d] = 1;
      }
    }
  }
  P.BooleanMultTranspose(1, l_v.HostRead(), 0, t_v.HostWrite());
  P.BooleanMult(1, t_v.HostRead(), 0, l_v.HostWrite());
  auto Touches = [&l_v](const mfem::Array<int> &dofs)
  { return std::ranges::any_of(dofs, [&l_v](int d) { return l_v[d] != 0; }); };
  elements.SetSize(0);
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    if (!InSource(e))
    {
      fes.GetElementDofs(e, dofs);
      if (Touches(dofs))
      {
        elements.Append(e);
      }
    }
  }
  bdr_elements.SetSize(0);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    fes.GetBdrElementDofs(be, dofs);
    if (Touches(dofs))
    {
      bdr_elements.Append(be);
    }
  }
}

void PlaneWaveOperator::AddLayerToMarkedElements(Mesh &mesh,
                                                 mfem::Array<int> &marked_elements) const
{
  // A nonconforming refinement of part of the layer of elements around the total-field
  // region of a plane wave leaves hanging entities inside the layer, where the source is
  // then not consistent with the incident field, and the solution degrades. Refining the
  // whole layer whenever part of it is marked avoids this, at the cost of its elements.
  // The layer here takes the elements on both sides of the boundary S of the region which
  // touch it, also only at a vertex, so that S stays conforming: hanging entities on S
  // would constrain the incident field of the refined side to its trace on the coarse
  // side, which is no longer the interpolant of a plane wave in the refined elements, and
  // this mismatch radiates into the scattered field. The hanging entities of the refined
  // layer are then inside the region, or outside of the layer, where the source does not
  // enter.
  mfem::ND_FECollection nd_fec(1, mesh.Get().Dimension());
  FiniteElementSpace nd_fespace(mesh, &nd_fec);
  mfem::Array<int> marked(mesh.Get().GetNE()), tdofs, elements, coupling, layer,
      bdr_neighbors;
  marked = 0;
  for (int e : marked_elements)
  {
    marked[e] = 1;
  }
  for (const auto &[idx, src] : sources)
  {
    GetLayer(idx, nd_fespace, tdofs, elements, &coupling);
    {
      mfem::Array<int> inner, in_layer(mesh.Get().GetNE());
      GetVertexNeighbors(idx, mesh.Get(), layer, bdr_neighbors);
      GetVertexNeighbors(idx, mesh.Get(), inner, bdr_neighbors, true);
      in_layer = 0;
      for (int e : layer)
      {
        in_layer[e] = 1;
      }
      for (const auto *list : {&coupling, &inner})
      {
        for (int e : *list)
        {
          if (!in_layer[e])
          {
            in_layer[e] = 1;
            layer.Append(e);
          }
        }
      }
    }
    int counts[2] = {0, layer.Size()};  // Marked and all elements of the layer
    for (int e : layer)
    {
      counts[0] += marked[e];
    }
    Mpi::GlobalSum(2, counts, mesh.Get().GetComm());
    if (counts[0] == 0 || counts[0] == counts[1])
    {
      continue;
    }
    for (int e : layer)
    {
      if (!marked[e])
      {
        marked[e] = 1;
        marked_elements.Append(e);
      }
    }
    Mpi::Print(" Plane wave {:d}: marked the whole layer of {:d} elements around the "
               "boundary of the total-field region ({:d} of them were marked)\n",
               idx, counts[1], counts[0]);
  }
}

void PlaneWaveOperator::CheckBoundaries(const config::BoundaryData &boundaries,
                                        const FiniteElementSpace &nd_fespace) const
{
  // The incident field does not satisfy the boundary conditions of the domain, so the
  // source leaks into the scattered field where a boundary with a condition touches the
  // layer of elements around the total-field region, crossing its boundary or touching it
  // from inside or outside. PEC and PMC boundaries are accepted where the incident field
  // satisfies them, as symmetry planes along the direction of propagation: with the
  // polarization (PEC) or the magnetic field (PMC) along the normal. Boundaries inside the
  // region which do not reach the layer are part of the scatterer, except absorbing
  // boundaries and Floquet ports: they truncate the scattered field, so the total-field
  // region must be separated from them by background elements. Where they touch the
  // region, the total field is driven only by the residual of the incident field on the
  // truncation outside of it, and it vanishes when the region encloses the rest.
  if (sources.empty())
  {
    return;
  }
  const mfem::ParFiniteElementSpace &fes = nd_fespace.Get();
  mfem::ParMesh &mesh = *fes.GetParMesh();
  const int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> checked(bdr_attr_max), pec(bdr_attr_max), pmc(bdr_attr_max);
  checked = 0;
  pec = 0;
  pmc = 0;
  for (auto attr : boundaries.attributes)
  {
    if (attr > 0 && attr <= bdr_attr_max)
    {
      checked[attr - 1] = 1;
    }
  }
  for (const auto &pair : boundaries.periodic.boundary_pairs)
  {
    for (const auto *attrs : {&pair.donor_attributes, &pair.receiver_attributes})
    {
      for (auto attr : *attrs)
      {
        if (attr > 0 && attr <= bdr_attr_max)
        {
          checked[attr - 1] = 0;
        }
      }
    }
  }
  for (auto attr : boundaries.pec.attributes)
  {
    if (attr > 0 && attr <= bdr_attr_max)
    {
      pec[attr - 1] = 1;
    }
  }
  for (auto attr : boundaries.pmc.attributes)
  {
    if (attr > 0 && attr <= bdr_attr_max)
    {
      pmc[attr - 1] = 1;
    }
  }
  mfem::Array<int> truncation(bdr_attr_max);
  truncation = 0;
  auto MarkTruncation = [&](const std::vector<int> &attrs)
  {
    for (auto attr : attrs)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        truncation[attr - 1] = 1;
      }
    }
  };
  MarkTruncation(boundaries.farfield.attributes);
  for (const auto &[port_idx, data] : boundaries.floquetport)
  {
    MarkTruncation(data.attributes);
  }
  mfem::Array<int> l_c, vdofs, neighbors, bdr_neighbors;
  mfem::Vector normal(3), h(3);
  for (const auto &[idx, src] : sources)
  {
    GetCoupling(idx, nd_fespace, l_c);
    std::vector<int> bad(bdr_attr_max, 0);
    // The magnetic field of the plane wave is along khat x p.
    h(0) = src.direction(1) * src.polarization(2) - src.direction(2) * src.polarization(1);
    h(1) = src.direction(2) * src.polarization(0) - src.direction(0) * src.polarization(2);
    h(2) = src.direction(0) * src.polarization(1) - src.direction(1) * src.polarization(0);
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      const int attr = mesh.GetBdrAttribute(be);
      if (attr <= 0 || attr > bdr_attr_max || !checked[attr - 1])
      {
        continue;
      }
      fes.GetBdrElementVDofs(be, vdofs);
      if (std::ranges::none_of(vdofs, [&l_c](int d) { return l_c[LocalDof(d)] != 0; }))
      {
        continue;
      }
      if (pec[attr - 1] || pmc[attr - 1])
      {
        auto &T = *mesh.GetBdrElementTransformation(be);
        T.SetIntPoint(&mfem::Geometries.GetCenter(T.GetGeometryType()));
        mfem::CalcOrtho(T.Jacobian(), normal);
        normal /= normal.Norml2();
        const mfem::Vector &field = pec[attr - 1] ? src.polarization : h;
        const double along = std::abs(normal * field) / field.Norml2();
        if (std::abs(normal * src.direction) <= 1.0e-6 && along >= 1.0 - 1.0e-6)
        {
          continue;
        }
      }
      bad[attr - 1] = 1;
    }
    Mpi::GlobalMax(bdr_attr_max, bad.data(), mesh.GetComm());
    std::vector<int> bad_attrs;
    for (int i = 0; i < bdr_attr_max; i++)
    {
      if (bad[i])
      {
        bad_attrs.push_back(i + 1);
      }
    }
    MFEM_VERIFY(bad_attrs.empty(),
                "Plane wave " << idx << ": the boundaries with attributes "
                              << fmt::format("{}", fmt::join(bad_attrs, ", "))
                              << " touch the layer of elements around the total-field "
                                 "region, where the incident field does not satisfy their "
                                 "boundary condition: the source would leak into the "
                                 "scattered field! The total-field region and its layer "
                                 "must be separated from the boundaries, except for PEC "
                                 "or PMC symmetry planes along the direction of "
                                 "propagation which the incident field satisfies.");

    // Absorbing boundaries and Floquet ports must not touch the total-field region at all.
    GetVertexNeighbors(idx, mesh, neighbors, bdr_neighbors);
    std::fill(bad.begin(), bad.end(), 0);
    for (int be : bdr_neighbors)
    {
      const int attr = mesh.GetBdrAttribute(be);
      if (attr > 0 && attr <= bdr_attr_max && truncation[attr - 1])
      {
        bad[attr - 1] = 1;
      }
    }
    Mpi::GlobalMax(bdr_attr_max, bad.data(), mesh.GetComm());
    bad_attrs.clear();
    for (int i = 0; i < bdr_attr_max; i++)
    {
      if (bad[i])
      {
        bad_attrs.push_back(i + 1);
      }
    }
    MFEM_VERIFY(bad_attrs.empty(),
                "Plane wave " << idx << ": the absorbing boundaries or Floquet ports with "
                              << "attributes "
                              << fmt::format("{}", fmt::join(bad_attrs, ", "))
                              << " touch the total-field region: they truncate the "
                                 "scattered field, so the total-field region must be "
                                 "separated from them by a layer of background elements!");
  }
}

const ComplexVector &
PlaneWaveOperator::GetLastIncidentField(int idx, double omega,
                                        const FiniteElementSpace &nd_fespace) const
{
  auto &[last_idx, last_omega, g] = last_incident;
  if (last_idx != idx || last_omega != omega || g.Size() != nd_fespace.GetTrueVSize())
  {
    g.SetSize(nd_fespace.GetTrueVSize());
    g.UseDevice(true);
    GetIncidentField(idx, omega, nd_fespace, g);
    last_idx = idx;
    last_omega = omega;
  }
  return g;
}

bool PlaneWaveOperator::AddExcitationVector(int idx, double omega, const ComplexOperator &A,
                                            const FiniteElementSpace &nd_fespace,
                                            ComplexVector &RHS) const
{
  if (sources.find(idx) == sources.end())
  {
    return false;
  }
  const Vector &tmark = GetTotalFieldMarker(idx, nd_fespace);

  const int n = nd_fespace.GetTrueVSize();
  const ComplexVector &g = GetLastIncidentField(idx, omega, nd_fespace);
  ComplexVector y(n);
  y.UseDevice(true);

  // b_T = -(A g_O)_T.
  ComplexVector gs(g);
  {
    const auto *m = tmark.HostRead();
    auto *gr = gs.Real().HostReadWrite();
    auto *gi = gs.Imag().HostReadWrite();
    for (int i = 0; i < n; i++)
    {
      gr[i] *= (1.0 - m[i]);
      gi[i] *= (1.0 - m[i]);
    }
  }
  A.Mult(gs, y);
  {
    const auto *m = tmark.HostRead();
    const auto *yr = y.Real().HostRead();
    const auto *yi = y.Imag().HostRead();
    auto *br = RHS.Real().HostReadWrite();
    auto *bi = RHS.Imag().HostReadWrite();
    for (int i = 0; i < n; i++)
    {
      br[i] -= m[i] * yr[i];
      bi[i] -= m[i] * yi[i];
    }
  }

  // b_O = (A g_T)_O.
  gs = g;
  {
    const auto *m = tmark.HostRead();
    auto *gr = gs.Real().HostReadWrite();
    auto *gi = gs.Imag().HostReadWrite();
    for (int i = 0; i < n; i++)
    {
      gr[i] *= m[i];
      gi[i] *= m[i];
    }
  }
  A.Mult(gs, y);
  {
    const auto *m = tmark.HostRead();
    const auto *yr = y.Real().HostRead();
    const auto *yi = y.Imag().HostRead();
    auto *br = RHS.Real().HostReadWrite();
    auto *bi = RHS.Imag().HostReadWrite();
    for (int i = 0; i < n; i++)
    {
      br[i] += (1.0 - m[i]) * yr[i];
      bi[i] += (1.0 - m[i]) * yi[i];
    }
  }
  return true;
}

void PlaneWaveOperator::SetUpBackground(const MaterialOperator &mat_op,
                                        const FiniteElementSpace &nd_fespace)
{
  // The incident field only enters the right-hand side through the elements outside of the
  // total-field region which are coupled to it, and it is the plane wave of a homogeneous
  // medium: these elements must have one isotropic, lossless, non-dispersive material. Each
  // process checks the attributes of its own elements, since the material operator only
  // knows the local attributes, and the results are reduced before any process aborts.
  constexpr double tol = 1.0e-12;
  auto IsScalar = [](const mfem::DenseMatrix &M, double diag)
  {
    for (int i = 0; i < M.Height(); i++)
    {
      for (int j = 0; j < M.Width(); j++)
      {
        if (std::abs(M(i, j) - ((i == j) ? diag : 0.0)) > tol * std::max(1.0, diag))
        {
          return false;
        }
      }
    }
    return true;
  };
  const mfem::ParFiniteElementSpace &fes = nd_fespace.Get();
  const mfem::ParMesh &mesh = *fes.GetParMesh();
  mfem::Array<int> l_c, layer;
  for (auto &[idx, src] : sources)
  {
    // The elements of the layer follow the prolongation (see GetCoupling): on a
    // nonconforming mesh, a refined element outside of the region shares no local dof with
    // the coarse elements of the region, and is coupled to them through the constraints of
    // its hanging faces and edges.
    GetCoupling(idx, nd_fespace, l_c, &layer);
    int bad_attr = 0, n_layer = 0;
    double eps_min = mfem::infinity(), eps_max = -mfem::infinity(),
           muinv_min = mfem::infinity(), muinv_max = -mfem::infinity();
    for (int e : layer)
    {
      const int attr = mesh.GetAttribute(e);
      n_layer++;
      const double eps = mat_op.GetPermittivityReal(attr)(0, 0),
                   muinv = mat_op.GetInvPermeability(attr)(0, 0);
      if (!IsScalar(mat_op.GetPermittivityReal(attr), eps) ||
          !IsScalar(mat_op.GetInvPermeability(attr), muinv) ||
          !IsScalar(mat_op.GetPermittivityImag(attr), 0.0) ||
          !IsScalar(mat_op.GetConductivity(attr), 0.0) ||
          !IsScalar(mat_op.GetInvLondonDepth(attr), 0.0) ||
          mat_op.GetFrequencyDependentPermittivityModelMaterial(attr) >= 0)
      {
        bad_attr = std::max(bad_attr, attr);
      }
      eps_min = std::min(eps_min, eps);
      eps_max = std::max(eps_max, eps);
      muinv_min = std::min(muinv_min, muinv);
      muinv_max = std::max(muinv_max, muinv);
    }
    Mpi::GlobalMax(1, &bad_attr, fes.GetComm());
    Mpi::GlobalSum(1, &n_layer, fes.GetComm());
    Mpi::GlobalMin(1, &eps_min, fes.GetComm());
    Mpi::GlobalMax(1, &eps_max, fes.GetComm());
    Mpi::GlobalMin(1, &muinv_min, fes.GetComm());
    Mpi::GlobalMax(1, &muinv_max, fes.GetComm());
    MFEM_VERIFY(
        n_layer > 0,
        "Plane wave " << idx << ": the total-field region must not fill the whole domain!");
    MFEM_VERIFY(bad_attr == 0,
                "Plane wave " << idx
                              << ": the boundary of the total-field region must be "
                                 "surrounded by an isotropic, lossless, non-dispersive "
                                 "material, but it "
                                 "touches domain attribute "
                              << bad_attr << " outside of the region!");
    MFEM_VERIFY(eps_max - eps_min <= tol * eps_max &&
                    muinv_max - muinv_min <= tol * muinv_max,
                "Plane wave " << idx
                              << ": the boundary of the total-field region must be "
                                 "surrounded by a single homogeneous material, but the "
                                 "elements outside of the region which touch it have "
                                 "different permittivities or permeabilities!");
    // PML regions truncate the scattered field, as absorbing boundaries do, and their
    // attributes have the material of the background, so they would pass the check above:
    // they must not overlap the total-field region nor touch its closure, even at a vertex.
    if (mat_op.HasPML())
    {
      const auto pml_attr = mat_op.GetPMLAttributes();
      auto IsPML = [&pml_attr](int attr)
      { return std::ranges::binary_search(pml_attr, attr); };
      mfem::Array<int> neighbors, bdr_neighbors;
      GetVertexNeighbors(idx, const_cast<mfem::ParMesh &>(mesh), neighbors, bdr_neighbors);
      int touching = 0;
      for (int e : neighbors)
      {
        touching =
            std::max(touching, IsPML(mesh.GetAttribute(e)) ? mesh.GetAttribute(e) : 0);
      }
      for (int e = 0; e < mesh.GetNE(); e++)
      {
        const int attr = mesh.GetAttribute(e);
        if (src.attr_marker[attr - 1] && IsPML(attr))
        {
          touching = std::max(touching, attr);
        }
      }
      Mpi::GlobalMax(1, &touching, fes.GetComm());
      MFEM_VERIFY(touching == 0,
                  "Plane wave " << idx << ": the PML region with domain attribute "
                                << touching
                                << " touches the total-field region: it truncates the "
                                   "scattered field, so the total-field region must be "
                                   "separated from it by a layer of background elements!");
    }

    src.c = std::sqrt(muinv_max / eps_max);
    if (mat_op.HasWaveVector())
    {
      src.kf = mat_op.GetWaveVector();
      src.kf_scaled = mat_op.HasFloquetFrequencyScaling();
    }
    Mpi::Print(" Plane wave {:d}: background relative permittivity = {:.3e}, relative "
               "permeability = {:.3e}\n",
               idx, eps_max, 1.0 / muinv_max);
  }
}

namespace
{

// The lattice translation of a periodic boundary pair without an explicit translation (as
// detected from the mesh), for an axis-aligned cell as for Floquet ports: the extent of the
// mesh along the common normal of the faces of the pair, if they are all normal to the same
// coordinate axis.
std::optional<mfem::Vector> AxisAlignedTranslation(mfem::ParMesh &mesh,
                                                   const config::PeriodicData &pair)
{
  int count[4] = {0, 0, 0, 0};  // Faces normal to x, y, z, and to none of them
  mfem::Vector normal(3);
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    const int attr = mesh.GetBdrAttribute(be);
    if (std::ranges::find(pair.donor_attributes, attr) == pair.donor_attributes.end() &&
        std::ranges::find(pair.receiver_attributes, attr) == pair.receiver_attributes.end())
    {
      continue;
    }
    auto &T = *mesh.GetBdrElementTransformation(be);
    T.SetIntPoint(&mfem::Geometries.GetCenter(T.GetGeometryType()));
    mfem::CalcOrtho(T.Jacobian(), normal);
    normal /= normal.Norml2();
    int d = 0;
    for (int i = 1; i < 3; i++)
    {
      d = (std::abs(normal(i)) > std::abs(normal(d))) ? i : d;
    }
    count[(std::abs(normal(d)) >= 1.0 - 1.0e-6) ? d : 3]++;
  }
  Mpi::GlobalSum(4, count, mesh.GetComm());
  const int n_axes = (count[0] > 0) + (count[1] > 0) + (count[2] > 0);
  if (count[3] > 0 || n_axes != 1)
  {
    return std::nullopt;
  }
  const int d = (count[0] > 0) ? 0 : ((count[1] > 0) ? 1 : 2);
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(
      mesh, mesh::AttrToMarker(mesh.attributes.Max(), mesh.attributes), false, bbmin,
      bbmax);
  mfem::Vector t(3);
  t = 0.0;
  t(d) = bbmax(d) - bbmin(d);
  return t;
}

}  // namespace

void PlaneWaveOperator::CheckPeriodicity(const config::PeriodicBoundaryData &periodic,
                                         const std::vector<double> &omega,
                                         mfem::ParMesh &mesh)
{
  periodic_cell = !periodic.boundary_pairs.empty();
  if (sources.empty() || periodic.boundary_pairs.empty())
  {
    return;
  }

  // Lattice translations of the periodic boundary pairs (nondimensional).
  std::vector<mfem::Vector> a;
  for (const auto &pair : periodic.boundary_pairs)
  {
    const auto &T = pair.affine_transform;
    if (std::ranges::none_of(T, [](double t) { return t != 0.0; }))
    {
      // Detected from the mesh: the translation of an axis-aligned cell is its extent, and
      // otherwise the interpolated field is checked at each frequency instead.
      if (auto t = AxisAlignedTranslation(mesh, pair))
      {
        a.push_back(std::move(*t));
      }
      else
      {
        check_interpolated = true;
      }
      continue;
    }
    for (int i = 0; i < 3; i++)
    {
      for (int j = 0; j < 3; j++)
      {
        MFEM_VERIFY(std::abs(T[i * 4 + j] - ((i == j) ? 1.0 : 0.0)) <= 1.0e-12,
                    "Plane wave excitations require periodic boundary pairs which are "
                    "translations, without rotation!");
      }
    }
    mfem::Vector t(3);
    for (int i = 0; i < 3; i++)
    {
      t(i) = T[i * 4 + 3] / mesh_scale;
    }
    a.push_back(std::move(t));
  }

  for (const auto &[idx, src] : sources)
  {
    // The largest mismatch over the frequencies of the sweep, and its message.
    double worst = 0.0;
    std::string worst_msg;
    for (const double w : omega)
    {
      mfem::Vector kf(3);
      kf = 0.0;
      if (src.kf.Size() == 3)
      {
        kf = src.kf;
        if (src.kf_scaled)
        {
          kf *= w;
        }
      }
      mfem::Vector dk(src.direction);
      dk *= w / src.c;
      dk -= kf;
      double mismatch = 0.0;
      for (const auto &t : a)
      {
        // A fixed k_F may differ from the tangential wave vector by a reciprocal lattice
        // vector, while a k_F which scales with the frequency must match it.
        const double phase = dk * t;
        mismatch = std::max(mismatch,
                            std::abs(src.kf_scaled
                                         ? phase
                                         : std::remainder(phase, 2.0 * std::numbers::pi)));
      }
      if (mismatch > 1.0e-6)
      {
        mfem::Vector k(src.direction);
        k *= w / src.c / mesh_scale;
        kf *= 1.0 / mesh_scale;
        const auto msg = fmt::format(
            "Plane wave {:d}: the incident field is not quasi-periodic with the Floquet "
            "wave vector across the periodic boundaries at {:.6e} GHz (phase mismatch "
            "{:.2e} rad over a lattice translation). The "
            "Floquet wave vector [{:.6e}, {:.6e}, {:.6e}] (zero if not given) must be the "
            "part of the wave vector of the plane wave [{:.6e}, {:.6e}, {:.6e}] along the "
            "periodic directions{}",
            idx, w * unit_GHz, mismatch, kf(0), kf(1), kf(2), k(0), k(1), k(2),
            src.kf_scaled ? ", with \"FloquetReferenceFrequency\""
                          : ", up to a reciprocal lattice vector");
        MFEM_VERIFY(mismatch <= 1.0e-2, msg << "!");
        if (mismatch > worst)
        {
          worst = mismatch;
          worst_msg = msg;
        }
      }
    }
    if (worst > 0.0)
    {
      Mpi::Warning("{}; the source leaks into the scattered field at the level of this "
                   "mismatch (the largest over the frequencies of the sweep).",
                   worst_msg);
    }
  }
}

bool PlaneWaveOperator::GetTotalAndScatteredField(int idx, double omega,
                                                  const FiniteElementSpace &nd_fespace,
                                                  const ComplexVector &x,
                                                  ComplexVector &E_t,
                                                  ComplexVector &E_s) const
{
  if (sources.find(idx) == sources.end())
  {
    return false;
  }
  const Vector &tmark = GetTotalFieldMarker(idx, nd_fespace);
  const int n = nd_fespace.GetTrueVSize();
  const ComplexVector &g = GetLastIncidentField(idx, omega, nd_fespace);
  E_t.SetSize(n);
  E_s.SetSize(n);
  const auto *m = tmark.HostRead();
  const auto *xr = x.Real().HostRead();
  const auto *xi = x.Imag().HostRead();
  const auto *gr = g.Real().HostRead();
  const auto *gi = g.Imag().HostRead();
  auto *tr = E_t.Real().HostWrite();
  auto *ti = E_t.Imag().HostWrite();
  auto *sr = E_s.Real().HostWrite();
  auto *si = E_s.Imag().HostWrite();
  for (int i = 0; i < n; i++)
  {
    tr[i] = xr[i] + (1.0 - m[i]) * gr[i];
    ti[i] = xi[i] + (1.0 - m[i]) * gi[i];
    sr[i] = xr[i] - m[i] * gr[i];
    si[i] = xi[i] - m[i] * gi[i];
  }
  return true;
}

}  // namespace palace
