// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "condensedenvironment.hpp"

#include <algorithm>
#include <string>
#include <fmt/core.h>
#include <fmt/ranges.h>
#include "fem/substructure.hpp"
#include "models/drivensubstructure.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

std::vector<int> EnvironmentPorts(const IoData &iodata, const DrivenSubstructure &ds,
                                  std::map<int, std::vector<char>> &port_dofs)
{
  std::vector<int> env_ports;
  for (const auto &[idx, data] : iodata.boundaries.lumpedport)
  {
    std::vector<int> attrs;
    for (const auto &elem : data.elements)
    {
      attrs.insert(attrs.end(), elem.attributes.begin(), elem.attributes.end());
    }
    auto mark = ds.BoundaryTrueDofs(attrs);
    double n[3] = {0.0, 0.0, 0.0};  // environment interior, interface, region interior
    for (std::size_t i = 0; i < mark.size(); i++)
    {
      n[0] += mark[i] && ds.EnvironmentInterior()[i];
      n[1] += mark[i] && ds.Interface()[i];
      n[2] += mark[i] && ds.RegionInterior()[i];
    }
    Mpi::GlobalSum(3, n, ds.GetComm());
    MFEM_VERIFY(n[1] == 0.0 && (n[0] == 0.0 || n[2] == 0.0),
                "Lumped port " << idx
                               << " touches the substructuring interface: a saved driven "
                                  "model needs the lumped ports away from it!");
    if (n[0] > 0.0)
    {
      env_ports.push_back(idx);
      port_dofs[idx] = std::move(mark);
    }
  }
  for (int idx : ds.EnvironmentWavePorts())
  {
    env_ports.push_back(-idx);
  }
  return env_ports;
}

std::vector<double> EnvironmentPortData(const IoData &iodata,
                                        const std::vector<int> &env_ports)
{
  std::vector<double> data;
  for (int p : env_ports)
  {
    if (p > 0)
    {
      const auto &lp = iodata.boundaries.lumpedport.at(p);
      data.insert(data.end(),
                  {lp.R, lp.L, lp.C, lp.Rs, lp.Ls, lp.Cs, lp.active ? 1.0 : 0.0});
    }
    else
    {
      const auto &wp = iodata.boundaries.waveport.at(-p);
      data.insert(data.end(), {static_cast<double>(wp.mode_idx), wp.d_offset,
                               wp.active ? 1.0 : 0.0, 0.0, 0.0, 0.0, 0.0});
    }
  }
  return data;
}

std::vector<std::complex<double>> FullSchur(const std::complex<double> *lower, int n,
                                            const SignatureMap *map)
{
  auto at = [&](int a, int b)
  {
    if (a < b)
    {
      std::swap(a, b);
    }
    return lower[static_cast<std::size_t>(b) * n -
                 static_cast<std::size_t>(b) * (b - 1) / 2 + (a - b)];
  };
  std::vector<std::complex<double>> S(static_cast<std::size_t>(n) * n);
  for (int j = 0; j < n; j++)
  {
    for (int i = j; i < n; i++)
    {
      std::complex<double> v = 0.0;
      if (map)
      {
        for (const auto &[a, ca] : map->rows[i])
        {
          for (const auto &[b, cb] : map->rows[j])
          {
            v += ca * cb * at(a, b);
          }
        }
      }
      else
      {
        v = at(i, j);
      }
      S[static_cast<std::size_t>(j) * n + i] = S[static_cast<std::size_t>(i) * n + j] = v;
    }
  }
  return S;
}

int FindValue(const std::vector<double> &y, double x, double tol)
{
  for (std::size_t j = 0; j < y.size(); j++)
  {
    if (std::abs(y[j] - x) <= tol * std::max(std::abs(x), std::abs(y[j])))
    {
      return static_cast<int>(j);
    }
  }
  return -1;
}

void RestrictPostprocessingToRegion(const mfem::ParMesh &mesh,
                                    const config::SubstructuringData &sub,
                                    const DrivenSubstructure &ds,
                                    const DrivenSubstructureModel &model, bool env_ports,
                                    const char *mode, config::DomainData &domains,
                                    config::BoundaryData &boundaries)
{
  domains.postpro.partial = true;
  std::vector<int> skipped;
  for (auto it = domains.postpro.energy.begin(); it != domains.postpro.energy.end();)
  {
    const bool region =
        std::ranges::all_of(it->second.attributes,
                            [&](int a)
                            {
                              return std::ranges::find(sub.region_attributes, a) !=
                                     sub.region_attributes.end();
                            });
    if (!region)
    {
      skipped.push_back(it->first);
      it = domains.postpro.energy.erase(it);
    }
    else
    {
      ++it;
    }
  }
  // On the mesh of the region, which lacks the environment's part of a surface, the
  // boundaries in the region are those of the model too.
  auto contains = [](const std::vector<int> &v, int a)
  { return std::ranges::find(v, a) != v.end(); };
  auto bdr_in_region = BoundaryAttributesInRegion(mesh, sub.region_attributes);
  if (ds.RegionOnly())
  {
    std::erase_if(bdr_in_region,
                  [&](int a) { return !contains(model.region_bdr_attributes, a); });
  }
  std::vector<int> skipped_flux, skipped_probe;
  for (auto it = boundaries.postpro.flux.begin(); it != boundaries.postpro.flux.end();)
  {
    const auto &attrs = it->second.attributes;
    if (attrs.empty() ||
        !std::ranges::all_of(attrs, [&](int a) { return contains(bdr_in_region, a); }))
    {
      skipped_flux.push_back(it->first);
      it = boundaries.postpro.flux.erase(it);
    }
    else
    {
      ++it;
    }
  }
  for (auto it = domains.postpro.probe.begin(); it != domains.postpro.probe.end();)
  {
    if (!PointInRegion(mesh, it->second.center, sub.region_attributes))
    {
      skipped_probe.push_back(it->first);
      it = domains.postpro.probe.erase(it);
    }
    else
    {
      ++it;
    }
  }
  auto none = [](const std::vector<int> &idx, const char *what)
  {
    return idx.empty() ? std::string()
                       : fmt::format(", no {} ({})", what, fmt::join(idx, ", "));
  };
  Mpi::Warning("{} driven substructuring writes no total energies or participation "
               "ratios{}{}{}{}!\n",
               mode, none(skipped, "energies of the domain postprocessing outside the region"),
               none(skipped_flux, "surface fluxes outside the region"),
               none(skipped_probe, "probes outside the region"),
               !env_ports ? std::string()
                          : ", and no power (or wave port voltage) of the "
                            "environment's ports");
}

}  // namespace palace
