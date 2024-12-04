// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP
#define PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP

#include <map>
#include <vector>
#include <nlohmann/json_fwd.hpp>

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

    // TODO: C++20 to replace this with iterator over joined range
    auto flatten_port_indices() const
    {
      std::vector<int> out;
      out.insert(out.end(), lumped_port.cbegin(), lumped_port.cend());
      out.insert(out.end(), wave_port.cbegin(), wave_port.cend());
      out.insert(out.end(), current_port.cbegin(), current_port.cend());
      return out;
    }
  };

  std::map<int, SingleExcitationSpec> excitations = {};

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

  [[nodiscard]] std::string FmtLog() const;
};

void to_json(nlohmann::json &j, const PortExcitationHelper::SingleExcitationSpec &p);

void from_json(const nlohmann::json &j, PortExcitationHelper::SingleExcitationSpec &p);

void to_json(nlohmann::json &j, const PortExcitationHelper &p);

void from_json(const nlohmann::json &j, PortExcitationHelper &p);

}  // namespace palace

#endif  // PALACE_MODELS_PORT_EXCITATION_MANAGER_HPP
