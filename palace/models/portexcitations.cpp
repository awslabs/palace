// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "portexcitations.hpp"

#include "currentdipoleoperator.hpp"
#include "lumpedportoperator.hpp"
#include "surfacecurrentoperator.hpp"
#include "waveportoperator.hpp"

#include <fmt/format.h>
#include <fmt/ranges.h>
#include <nlohmann/json.hpp>

namespace palace
{

[[nodiscard]] std::string PortExcitations::FmtLog() const
{
  fmt::memory_buffer buf{};
  auto to = [&buf](auto f, auto &&...a)  // mini-lambda for cleaner code
  { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

  int i = 1;
  for (const auto &[idx, ex] : excitations)
  {
    to("Excitation{} with index {:d} has contributions from:\n",
       (Size() > 1) ? fmt::format(" {:d}/{:d}", i, Size()) : "", idx);
    if (!ex.lumped_port.empty())
    {
      to(" Lumped port{} {:2d}\n", (ex.lumped_port.size() > 1) ? "s" : "",
         fmt::join(ex.lumped_port, " "));
    }
    if (!ex.wave_port.empty())
    {
      to(" Wave port{} {:2d}\n", (ex.wave_port.size() > 1) ? "s" : "",
         fmt::join(ex.wave_port, " "));
    }
    if (!ex.current_port.empty())
    {
      to(" Surface current port{} {:2d}\n", (ex.current_port.size() > 1) ? "s" : "",
         fmt::join(ex.current_port, " "));
    }
    if (!ex.current_dipole.empty())
    {
      to(" Current dipole{} {:2d}\n", (ex.current_dipole.size() > 1) ? "s" : "",
         fmt::join(ex.current_dipole, " "));
    }
    i++;
  }
  return fmt::to_string(buf);
}

void to_json(nlohmann::json &j, const PortExcitations::SingleExcitationSpec &p)
{
  j = nlohmann::json{{"LumpedPort", p.lumped_port},
                     {"WavePort", p.wave_port},
                     {"SurfaceCurrent", p.current_port},
                     {"CurrentDipole", p.current_dipole}};
}

void from_json(const nlohmann::json &j, PortExcitations::SingleExcitationSpec &p)
{
  j.at("LumpedPort").get_to(p.lumped_port);
  j.at("WavePort").get_to(p.wave_port);
  j.at("SurfaceCurrent").get_to(p.current_port);
  j.at("CurrentDipole").get_to(p.current_dipole);
}

void to_json(nlohmann::json &j, const PortExcitations &p)
{
  j = nlohmann::json{p.excitations};
}

void from_json(const nlohmann::json &j, PortExcitations &p)
{
  j.get_to(p.excitations);
}

PortExcitations::PortExcitations(const LumpedPortOperator &lumped_port_op,
                                 const WavePortOperator &wave_port_op,
                                 const SurfaceCurrentOperator &surf_j_op,
                                 const CurrentDipoleOperator &dipole_op)
{
  for (const auto &[idx, port] : lumped_port_op)
  {
    if (!port.HasExcitation())
    {
      continue;
    }
    excitations.try_emplace(port.excitation, SingleExcitationSpec{});  // If not present
    excitations.at(port.excitation).lumped_port.push_back(idx);
  }
  for (const auto &[idx, port] : wave_port_op)
  {
    if (!port.HasExcitation())
    {
      continue;
    }
    excitations.try_emplace(port.excitation, SingleExcitationSpec{});
    excitations.at(port.excitation).wave_port.push_back(idx);
  }

  // Surface currents are always excited. Add them to all single existing excitations.
  // TODO: Add excitation 1 if not present already?
  std::vector<int> current_port_idx;
  for (const auto &[idx, port] : surf_j_op)
  {
    current_port_idx.push_back(idx);
  }
  if (!current_port_idx.empty())
  {
    for (auto &[ex_idx, ex_spec] : excitations)
    {
      ex_spec.current_port = current_port_idx;
    }
  }

  // Current dipoles are always excited. Add them to all existing excitations.
  std::vector<int> current_dipole_idx;
  for (const auto &[idx, dipole] : dipole_op)
  {
    current_dipole_idx.push_back(idx);
  }
  if (!current_dipole_idx.empty())
  {
    if (excitations.empty())
    {
      excitations.try_emplace(1, SingleExcitationSpec{});
    }
    for (auto &[ex_idx, ex_spec] : excitations)
    {
      ex_spec.current_dipole = current_dipole_idx;
    }
  }
};

}  // namespace palace
