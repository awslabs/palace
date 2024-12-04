// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP
#define PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP

#include <map>
#include <vector>
#include <nlohmann/json_fwd.hpp>
#include "utils/strongtype.hpp"

namespace palace
{

// Small helper class to collect data of what (lumped / wave / surface) ports are
// excited in driven and transient simulation, as stored in space_op;
// Manges indices.

enum class PortType : std::uint8_t
{
  LumpedPort = 0,
  WavePort = 1,
  CurrentPort = 2,
  Undefined = 3
};

class PortExcitationHelper
{
public:
  struct SingleExcitationSpec
  {
    std::vector<int> lumped_port = {};
    std::vector<int> wave_port = {};
    std::vector<int> current_port = {};

    // TODO: C++20 to replace this with iterator over joined range
    auto flatten_port_indices() const
    {
      std::vector<int> out;
      out.insert(out.end(), lumped_port.cbegin(), lumped_port.cend());
      out.insert(out.end(), wave_port.cbegin(), wave_port.cend());
      out.insert(out.end(), current_port.cbegin(), current_port.cend());
      return out;
    }

    // Only a single port is excitated in the excitation
    std::tuple<bool, PortType, int> is_simple() const
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

  std::map<ExcitationIdx, SingleExcitationSpec> excitations = {};

  auto begin() { return excitations.begin(); }
  auto end() { return excitations.end(); }
  auto begin() const { return excitations.begin(); }
  auto end() const { return excitations.end(); }

  template <typename LumpedportIt, typename WaveportIt, typename CurrentIt>
  PortExcitationHelper(const LumpedportIt &lumped_port_op, const WaveportIt &wave_port_op,
                       const CurrentIt &surf_j_op)
  {
    for (const auto &[idx, port] : lumped_port_op)
    {
      if (port.excitation ==
          ExcitationIdx(0))  // LumpedPortData does not have HasExcitaiton
      {
        continue;
      }
      excitations.try_emplace(port.excitation, SingleExcitationSpec{});  // If not present
      excitations.at(port.excitation).lumped_port.push_back(idx);
    }
    for (const auto &[idx, port] : wave_port_op)
    {
      if (port.excitation == ExcitationIdx(0))  // WavePortData does not have HasExcitaiton
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

  [[nodiscard]] ExcitationIdx MaxIdx() const
  {
    if (excitations.empty())
    {
      return ExcitationIdx(0);
    }
    // Map is stored order by key so max key is last item
    return std::next(std::rend(excitations))->first;
  }
  [[nodiscard]] auto Size() const { return excitations.size(); }
  [[nodiscard]] auto Empty() const { return excitations.empty(); }

  [[nodiscard]] std::string FmtLog() const;

  // Single Simple (only 1 port per excitation) Excitation 
  [[nodiscard]] std::tuple<bool, ExcitationIdx, PortType, int> IsSingleSimple() const
  {
    if (Size() == 1)
    {
      const auto &[ex_idx, ex_spec] = *excitations.begin();
      const auto [is_simple, port_type, port_idx] = ex_spec.is_simple();
      if (is_simple)
      {
        return std::make_tuple(true, ex_idx, port_type, port_idx);
      }
    }
    return std::make_tuple(false, ExcitationIdx(0), PortType::Undefined, 0);
  }

  // Multiple Simple (only 1 port per excitation) Excitation 
  [[nodiscard]] bool IsMultipleSimple() const
  {
    return std::all_of(excitations.begin(), excitations.end(),
                       [](const auto &ex) { return std::get<0>(ex.second.is_simple()); });
  }
};

void to_json(nlohmann::json &j, const PortExcitationHelper::SingleExcitationSpec &p);

void from_json(const nlohmann::json &j, PortExcitationHelper::SingleExcitationSpec &p);

void to_json(nlohmann::json &j, const PortExcitationHelper &p);

void from_json(const nlohmann::json &j, PortExcitationHelper &p);

}  // namespace palace

#endif  // PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP
