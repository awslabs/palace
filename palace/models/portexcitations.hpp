// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP
#define PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP

#include <map>
#include <vector>
#include <nlohmann/json_fwd.hpp>

namespace palace
{
class LumpedPortOperator;
class WavePortOperator;
class SurfaceCurrentOperator;

// Small helper class to collect data of what (lumped / wave / surface) ports are
// excited in driven and transient simulation, as stored in space_op;
// Manages indices.

enum class PortType : std::uint8_t
{
  LumpedPort = 0,
  WavePort = 1,
  CurrentPort = 2,
  Undefined = 3
};

class PortExcitations
{
public:
  struct SingleExcitationSpec
  {
    std::vector<int> lumped_port = {};
    std::vector<int> wave_port = {};
    std::vector<int> current_port = {};

    // TODO: C++20 to replace this with iterator over joined range.
    auto FlatternPortIndices() const
    {
      std::vector<int> out;
      out.insert(out.end(), lumped_port.cbegin(), lumped_port.cend());
      out.insert(out.end(), wave_port.cbegin(), wave_port.cend());
      out.insert(out.end(), current_port.cbegin(), current_port.cend());
      return out;
    }

    // Only a single port is excited.
    std::tuple<bool, PortType, int> IsSimple() const
    {
      auto n_lumped = lumped_port.size();
      auto n_wave = wave_port.size();
      auto n_current = current_port.size();

      if (n_lumped == 1 && n_wave == 0 && n_current == 0)
      {
        return std::make_tuple(true, PortType::LumpedPort, lumped_port.at(0));
      }
      else if (n_lumped == 0 && n_wave == 1 && n_current == 0)
      {
        return std::make_tuple(true, PortType::WavePort, wave_port.at(0));
      }
      else if (n_lumped == 0 && n_wave == 0 && n_current == 1)
      {
        return std::make_tuple(true, PortType::CurrentPort, current_port.at(0));
      }
      else
      {
        return std::make_tuple(false, PortType::Undefined, 0);
      }
    }
  };

  std::map<int, SingleExcitationSpec> excitations = {};

  auto begin() { return excitations.begin(); }
  auto end() { return excitations.end(); }
  auto begin() const { return excitations.begin(); }
  auto end() const { return excitations.end(); }

  PortExcitations(const LumpedPortOperator &lumped_port_op,
                  const WavePortOperator &wave_port_op,
                  const SurfaceCurrentOperator &surf_j_op);

  [[nodiscard]] int MaxIdx() const
  {
    // Map is stored order by key so max key is last item.
    return excitations.empty() ? 0 : std::next(std::rend(excitations))->first;
  }
  [[nodiscard]] auto Size() const { return excitations.size(); }
  [[nodiscard]] auto Empty() const { return excitations.empty(); }

  [[nodiscard]] std::string FmtLog() const;

  // Single Simple (only 1 port per excitation) Excitation.
  [[nodiscard]] std::tuple<bool, int, PortType, int> IsSingleSimple() const
  {
    if (Size() == 1)
    {
      const auto &[ex_idx, ex_spec] = *excitations.begin();
      const auto [is_simple, port_type, port_idx] = ex_spec.IsSimple();
      if (is_simple)
      {
        return std::make_tuple(true, ex_idx, port_type, port_idx);
      }
    }
    return std::make_tuple(false, 0, PortType::Undefined, 0);
  }

  // Multiple Simple (only 1 port per excitation) Excitation.
  [[nodiscard]] bool IsMultipleSimple() const
  {
    return std::all_of(excitations.begin(), excitations.end(),
                       [](const auto &ex) { return std::get<0>(ex.second.IsSimple()); });
  }
};

void to_json(nlohmann::json &j, const PortExcitations::SingleExcitationSpec &p);

void from_json(const nlohmann::json &j, PortExcitations::SingleExcitationSpec &p);

void to_json(nlohmann::json &j, const PortExcitations &p);

void from_json(const nlohmann::json &j, PortExcitations &p);

}  // namespace palace

#endif  // PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP
