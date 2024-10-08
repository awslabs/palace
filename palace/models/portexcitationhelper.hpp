// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP
#define PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP

#include <map>
#include <fmt/format.h>

namespace palace
{

// Small helper class to collect data of what (lumped / wave / surface) ports are
// excited in driven and transient simulation, as stored in space_op;
// Manges indices.
class PortExcitationHelper
{
public:
  struct SingleExcitationSpec
  {
    std::vector<int> lumped_port = {};
    std::vector<int> wave_port = {};
    std::vector<int> current_port = {};
  };

  std::map<int, SingleExcitationSpec> excitations = {};

  template <typename LumpedportIt, typename WaveportIt, typename CurrentIt>
  PortExcitationHelper(const LumpedportIt &lumped_port_op, const WaveportIt &wave_port_op,
                        const CurrentIt &surf_j_op)
  {
    for (const auto &[idx, port] : lumped_port_op)
    {
      if (port.excitation == 0)  // Pre-ranges filter
      {
        continue;
      }
      excitations.try_emplace(port.excitation, SingleExcitationSpec{});  // If not present
      excitations.at(port.excitation).lumped_port.push_back(idx);
    }
    for (const auto &[idx, port] : wave_port_op)
    {
      if (port.excitation == 0)
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
  };

  [[nodiscard]] int MaxIdx() const
  {
    if (excitations.empty())
    {
      return 0;
    }
    // Map is stored order by key so max key is last item
    return std::next(std::rend(excitations))->first;
  }
  [[nodiscard]] auto Size() const { return excitations.size(); }
  [[nodiscard]] auto Empty() const { return excitations.empty(); }

  [[nodiscard]] std::string FmtLog() const
  {
    fmt::memory_buffer buf{};
    auto to = [&buf](auto f, auto &&...a)  // mini-lambda for cleaner code
    { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

    int i = 0;
    for (const auto &[idx, ex] : excitations)
    {
      to("Excitation{} with index {:2d} has contribution from:\n",
         (Size() > 1) ? fmt::format(" {:2d}/{:2d}", i, Size()) : "", idx);
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
    }
    return fmt::to_string(buf);
  }
};

}  // namespace palace

#endif  // PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP
