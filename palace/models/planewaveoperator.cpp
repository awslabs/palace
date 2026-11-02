// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "planewaveoperator.hpp"

#include <cmath>
#include <span>
#include <fmt/ranges.h>
#include "fem/fespace.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

namespace palace
{

PlaneWaveOperator::PlaneWaveOperator(const std::map<int, config::PlaneWaveData> &plane_wave,
                                     const Units &units, const mfem::ParMesh &mesh)
{
  MFEM_VERIFY(plane_wave.empty() || mesh.SpaceDimension() == 3,
              "Plane wave excitations require a 3D mesh!");
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

  if (sources.empty())
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
  Mpi::Print(" Plane wave {:d}: {:d} total-field unknowns (of {:d})\n", idx, n_tf,
             fes.GlobalTrueVSize());
  return tmark;
}

bool PlaneWaveOperator::AddExcitationVector(int idx, double omega, const ComplexOperator &A,
                                            const FiniteElementSpace &nd_fespace,
                                            ComplexVector &RHS) const
{
  auto it = sources.find(idx);
  if (it == sources.end())
  {
    return false;
  }
  const auto &src = it->second;
  const Vector &tmark = GetTotalFieldMarker(idx, nd_fespace);

  // Interpolate the incident field E_inc = E0 p exp(-i k khat . (x - x0)), with k = omega
  // in the (vacuum) background in nondimensional units.
  auto &fes = const_cast<mfem::ParFiniteElementSpace &>(nd_fespace.Get());
  const int n = fes.GetTrueVSize();
  ComplexVector g(n), y(n);
  g.UseDevice(true);
  y.UseDevice(true);
  for (int part = 0; part < 2; part++)
  {
    mfem::VectorFunctionCoefficient coef(
        3,
        [&src, omega, part](const mfem::Vector &x, mfem::Vector &E)
        {
          double phase = 0.0;
          for (int d = 0; d < 3; d++)
          {
            phase += src.direction(d) * (x(d) - src.origin(d));
          }
          phase *= omega;
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
    gf.ParallelProject((part == 0) ? g.Real() : g.Imag());
  }

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

bool PlaneWaveOperator::TouchesTotalField(int idx, const FiniteElementSpace &nd_fespace,
                                          const mfem::Array<int> &bdr_attr_marker) const
{
  const mfem::ParFiniteElementSpace &fes = nd_fespace.Get();
  const mfem::ParMesh &mesh = *fes.GetParMesh();
  mfem::Array<int> lmark, vdofs;
  GetLocalTotalFieldMarker(idx, fes, lmark);
  bool touches = false;
  for (int i = 0; i < mesh.GetNBE() && !touches; i++)
  {
    const int attr = mesh.GetBdrAttribute(i);
    if (attr > bdr_attr_marker.Size() || !bdr_attr_marker[attr - 1])
    {
      continue;
    }
    fes.GetBdrElementVDofs(i, vdofs);
    for (int d : vdofs)
    {
      if (lmark[d >= 0 ? d : -1 - d])
      {
        touches = true;
        break;
      }
    }
  }
  Mpi::GlobalOr(1, &touches, fes.GetComm());
  return touches;
}

}  // namespace palace
