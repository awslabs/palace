// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "portexcitationhelper.hpp"
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <nlohmann/json.hpp>

namespace palace
{

[[nodiscard]] std::string PortExcitationHelper::FmtLog() const
{
  fmt::memory_buffer buf{};
  auto to = [&buf](auto f, auto &&...a)  // mini-lambda for cleaner code
  { fmt::format_to(std::back_inserter(buf), f, std::forward<decltype(a)>(a)...); };

  int i = 1;
  for (const auto &[idx, ex] : excitations)
  {
    to("Excitation{} with index {:d} has contribution from:\n",
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
    i++;
  }
  return fmt::to_string(buf);
}

void to_json(nlohmann::json &j, const PortExcitationHelper::SingleExcitationSpec &p)
{
  j = nlohmann::json{{"LumpedPort", p.lumped_port},
                     {"WavePort", p.wave_port},
                     {"SurfaceCurrent", p.current_port}};
}

void from_json(const nlohmann::json &j, PortExcitationHelper::SingleExcitationSpec &p)
{
  j.at("LumpedPort").get_to(p.lumped_port);
  j.at("WavePort").get_to(p.wave_port);
  j.at("SurfaceCurrent").get_to(p.current_port);
}

void to_json(nlohmann::json &j, const PortExcitationHelper &p)
{
  j = nlohmann::json{p.excitations};
}

void from_json(const nlohmann::json &j, PortExcitationHelper &p)
{
  j.get_to(p.excitations);
}

}  // namespace palace
