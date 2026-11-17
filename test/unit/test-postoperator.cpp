// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <algorithm>
#include <cmath>
#include <complex>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <numbers>
#include <vector>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "drivers/fluxerrorestimator.hpp"
#include "fem/errorindicator.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
#include "fem/mesh.hpp"
#include "fixtures.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/laplaceoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/spaceoperator.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
// Helpers

class DrivenPostOperatorTest : public PostOperator<ProblemType::DRIVEN>
{
public:
  using PostOperator<ProblemType::DRIVEN>::PostOperator;

  Measurement::PortPostData MeasureLumpedPortForTest(int idx, const ComplexVector &e,
                                                     const ComplexVector &b,
                                                     std::complex<double> omega)
  {
    measurement_cache.freq = omega;
    SetEGridFunction(e);
    SetBGridFunction(b);
    MeasureLumpedPorts();
    return measurement_cache.lumped_port_vi.at(idx);
  }
};

// random integer between 0 and n as a double.
auto randd(int n)
{
  return double(rand() % n);
};

// Random measurement sample over ndomains.
auto RandomMeasurement(int ndomain = 5)
{
  Measurement cache;
  cache.domain_E_field_energy_all = 0.0;
  cache.domain_H_field_energy_all = 0.0;
  for (int i = 0; i < ndomain; i++)
  {
    cache.domain_H_field_energy_i.emplace_back(
        Measurement::DomainData{i, 1 + randd(100), 0});
    cache.domain_H_field_energy_all += cache.domain_H_field_energy_i.back().energy;
    cache.domain_E_field_energy_i.emplace_back(
        Measurement::DomainData{i, 1 + randd(100), 0});
    cache.domain_E_field_energy_all += cache.domain_H_field_energy_i.back().energy;
  }
  for (auto &H : cache.domain_H_field_energy_i)
  {
    H.participation_ratio = H.energy / cache.domain_H_field_energy_all;
  }
  for (auto &E : cache.domain_E_field_energy_i)
  {
    E.participation_ratio = E.energy / cache.domain_E_field_energy_all;
  }

  cache.freq = std::complex(1 + randd(99), randd(100));
  cache.eigenmode_Q = 5e6 * (randd(99) + 1);
  cache.error_abs = 1e-9 * (randd(99) + 1);
  cache.error_bkwd = 1e-9 * (randd(99) + 1);

  int nport = 2 + (rand() % 10);
  for (int i = 1; i <= nport - 2; i++)
  {
    if (rand() % 2 == 0)
    {
      cache.lumped_port_vi[i];
    }
    else
    {
      cache.wave_port_vi[i];
    }
  }
  // always make last two one of each. Others are randomly allocated.
  cache.lumped_port_vi[nport - 1];
  cache.wave_port_vi[nport];
  cache.lumped_port_capacitor_energy = 0;
  cache.lumped_port_inductor_energy = 0;
  using Branch = LumpedPortData::Branch;
  auto impedance = [](std::complex<double> omega, Branch b, double R, double L, double C)
  {
    std::complex<double> Y;
    if (std::abs(R) > 0.0 && (b == Branch::TOTAL || b == Branch::R))
    {
      Y += 1.0 / R;
    }
    if (std::abs(L) > 0.0 && (b == Branch::TOTAL || b == Branch::L))
    {
      Y += 1.0 / (std::complex(0.0, 1.0) * omega * L);
    }
    if (std::abs(C) > 0.0 && (b == Branch::TOTAL || b == Branch::C))
    {
      Y += std::complex(0.0, 1.0) * omega * C;
    }
    MFEM_VERIFY(std::abs(Y) > 0.0,
                "Characteristic impedance requested for lumped port with zero admittance!")
    return 1.0 / Y;
  };
  for (auto &[i, l] : cache.lumped_port_vi)
  {
    double R = 1 + randd(100), L = 1 + randd(100), C = 1 + randd(100);
    l.V = {1 + randd(100), -50 + randd(100)};
    l.I_RLC[0] = l.V / impedance(cache.freq, Branch::R, R, L, C);
    l.I_RLC[1] = l.V / impedance(cache.freq, Branch::L, R, L, C);
    l.I_RLC[2] = l.V / impedance(cache.freq, Branch::C, R, L, C);
    l.I = std::accumulate(l.I_RLC.begin(), l.I_RLC.end(), std::complex(0.0, 0.0));
    l.P = l.V * std::conj(l.I);
    // Random sub-unit magnitude and random phase
    l.S = randd(1000) / 1000 * std::exp(std::complex(randd(100) / 100, randd(100) / 100));

    l.inductor_energy = (1 + randd(100));
    l.capacitor_energy = (1 + randd(100));
    l.inductive_energy_participation = (1 + randd(9999)) / 10000;
    l.mode_port_kappa = (1 + randd(100));
    l.quality_factor = 1e9 / (1 + randd(9999));
  }

  for (auto &[i, l] : cache.wave_port_vi)
  {
    l.P = {1 + randd(100), -50 + randd(100)};
    // Random sub-unit magnitude and random phase
    l.S = randd(1000) / 1000 * std::exp(std::complex(randd(100) / 100, randd(100) / 100));
  }

  int nprobe = 5;
  for (int i = 0; i < nprobe; i++)
  {
    cache.probe_E_field.push_back(std::complex(randd(100) / 100, randd(100) / 100));
    cache.probe_B_field.push_back(std::complex(randd(100) / 100, randd(100) / 100));
  }

  int nsurf = 5;
  for (int i = 0; i < nsurf; i++)
  {
    cache.surface_flux_i.emplace_back(
        Measurement::FluxData{i, std::complex((1 + randd(99)) / 100, randd(100) / 100),
                              static_cast<SurfaceFlux>(rand() % 3)});
  }
  int ninterface = 5;
  for (int i = 0; i < ninterface; i++)
  {
    cache.interface_eps_i.emplace_back(
        Measurement::InterfaceData{i, 1 + randd(100), (1 + randd(9999)) / 10000,
                                   (1 + randd(9999)) / 10000, 1e9 / (1 + randd(9999))});
    cache.interface_edge_i.emplace_back(Measurement::InterfaceEdgeData{
        i, 1 + randd(100), 1 + randd(100), (1 + randd(9999)) / 10000, 1 + randd(100),
        (1 + randd(9999)) / 10000});
    cache.interface_local_edge_i.emplace_back(Measurement::InterfaceLocalEdgeData{
        i,
        i + 1,
        true,
        i + 2,
        i + 3,
        {1, 2},
        {0.0, 0.0, 1.0},
        {1 + randd(100), 1 + randd(100), 1 + randd(100)},
        {1 + randd(100), 1 + randd(100), 1 + randd(100)},
        1 + randd(100),
        1 + randd(100),
        1 + randd(100),
        (1 + randd(9999)) / 10000,
        {1 + randd(100), 1 + randd(100)},
        {(1 + randd(9999)) / 10000, (1 + randd(9999)) / 10000},
        1 + randd(100),
        (1 + randd(9999)) / 10000,
        {1 + randd(100), 1 + randd(100)},
        {(1 + randd(9999)) / 10000, (1 + randd(9999)) / 10000},
        1 + randd(100),
        (1 + randd(9999)) / 10000,
        {1 + randd(100), 1 + randd(100)},
        {(1 + randd(9999)) / 10000, (1 + randd(9999)) / 10000},
        1 + randd(100),
        (1 + randd(9999)) / 10000,
        1 + randd(100),
        (1 + randd(9999)) / 10000,
        1 + randd(100),
        (1 + randd(9999)) / 10000,
        {1 + randd(100), 1 + randd(100), 1 + randd(100), 1 + randd(100), 1 + randd(100),
         1 + randd(100)},
        {(1 + randd(9999)) / 10000, (1 + randd(9999)) / 10000, (1 + randd(9999)) / 10000,
         (1 + randd(9999)) / 10000, (1 + randd(9999)) / 10000, (1 + randd(9999)) / 10000}});
  }

  int nfarfield = 5;
  for (int i = 0; i < nfarfield; i++)
  {
    cache.farfield.thetaphis.emplace_back(randd(180), randd(360));
    cache.farfield.E_field.push_back(
        {std::complex((1 + randd(99)) / 100, randd(100) / 100),
         std::complex((1 + randd(99)) / 100, randd(100) / 100),
         std::complex((1 + randd(99)) / 100, randd(100) / 100)});
  }

  cache.mode_data.kn = std::complex((1 + randd(99)) / 100, randd(100) / 100);
  cache.mode_data.n_eff = std::complex(1 + randd(100) / 100, randd(100) / 1000);
  int nmode_postpro = 2;
  for (int i = 0; i < nmode_postpro; i++)
  {
    // The first mode has no V/I impedance and the second no Z0, which must be left
    // unscaled.
    cache.mode_data.impedance[i] = {std::complex((1 + randd(99)) / 100, -randd(100) / 100),
                                    std::complex((1 + randd(99)) / 100, -randd(100) / 100),
                                    i == 0, i == 1};
    cache.mode_data.voltage[i] = {std::complex((1 + randd(99)) / 100, randd(100) / 100)};
  }

  return cache;
}

TEST_CASE("PostOperator", "[idempotent][Serial]")
{
  auto cache = RandomMeasurement();
  Units units(1e-6, 152 * 1e-6);  // μm, 152μm
  auto dim_cache = Measurement::Dimensionalize(units, cache);
  auto non_dim_cache = Measurement::Nondimensionalize(units, dim_cache);
  CHECK_THAT(cache.domain_E_field_energy_all,
             Catch::Matchers::WithinRel(non_dim_cache.domain_E_field_energy_all));
  CHECK_THAT(cache.domain_E_field_energy_all,
             !Catch::Matchers::WithinRel(dim_cache.domain_E_field_energy_all));
  CHECK_THAT(cache.domain_H_field_energy_all,
             Catch::Matchers::WithinRel(non_dim_cache.domain_H_field_energy_all));
  CHECK_THAT(cache.domain_H_field_energy_all,
             !Catch::Matchers::WithinRel(dim_cache.domain_H_field_energy_all));

  for (std::size_t i = 0; i < cache.domain_E_field_energy_i.size(); i++)
  {
    {
      const auto &c = cache.domain_E_field_energy_i[i];
      const auto &ndc = non_dim_cache.domain_E_field_energy_i[i];
      CHECK(c.idx == ndc.idx);
      CHECK_THAT(c.energy, Catch::Matchers::WithinRel(ndc.energy));
      CHECK_THAT(c.participation_ratio,
                 Catch::Matchers::WithinRel(ndc.participation_ratio));

      const auto &dc = dim_cache.domain_E_field_energy_i[i];
      CHECK(c.idx == dc.idx);
      CHECK_THAT(c.energy, !Catch::Matchers::WithinRel(dc.energy));
      CHECK_THAT(c.participation_ratio, Catch::Matchers::WithinRel(dc.participation_ratio));
    }
    {
      const auto &c = cache.domain_H_field_energy_i[i];
      const auto &ndc = non_dim_cache.domain_H_field_energy_i[i];
      CHECK(c.idx == ndc.idx);
      CHECK_THAT(c.energy, Catch::Matchers::WithinRel(ndc.energy));
      CHECK_THAT(c.participation_ratio,
                 Catch::Matchers::WithinRel(ndc.participation_ratio));

      const auto &dc = dim_cache.domain_H_field_energy_i[i];
      CHECK(c.idx == dc.idx);
      CHECK_THAT(c.energy, !Catch::Matchers::WithinRel(dc.energy));
      CHECK_THAT(c.participation_ratio, Catch::Matchers::WithinRel(dc.participation_ratio));
    }
  }
  auto check_port_data = [&](const auto &c, const auto &dc, const auto &ndc)
  {
    CAPTURE(c, dc, ndc);
    // Nondimensional
    CHECK_THAT(std::abs(c.P), Catch::Matchers::WithinRel(std::abs(ndc.P)));
    CHECK_THAT(std::arg(c.P), Catch::Matchers::WithinRel(std::arg(ndc.P)));
    CHECK_THAT(std::abs(c.S), Catch::Matchers::WithinRel(std::abs(ndc.S)));
    CHECK_THAT(std::arg(c.S), Catch::Matchers::WithinRel(std::arg(ndc.S)));
    CHECK_THAT(std::abs(c.V), Catch::Matchers::WithinRel(std::abs(ndc.V)));
    CHECK_THAT(std::arg(c.V), Catch::Matchers::WithinRel(std::arg(ndc.V)));
    CHECK_THAT(std::abs(c.I), Catch::Matchers::WithinRel(std::abs(ndc.I)));
    CHECK_THAT(std::arg(c.I), Catch::Matchers::WithinRel(std::arg(ndc.I)));
    for (auto i : {0, 1, 2})
    {
      CHECK_THAT(std::abs(c.I_RLC[i]), Catch::Matchers::WithinRel(std::abs(ndc.I_RLC[i])));
      CHECK_THAT(std::arg(c.I_RLC[i]), Catch::Matchers::WithinRel(std::arg(ndc.I_RLC[i])));
    }
    auto sum_cI = std::accumulate(c.I_RLC.begin(), c.I_RLC.end(), std::complex{0.0, 0.0});
    CHECK_THAT(std::abs(c.I), Catch::Matchers::WithinRel(std::abs(sum_cI)));
    CHECK_THAT(std::arg(c.I), Catch::Matchers::WithinRel(std::arg(sum_cI)));
    auto sum_ndcI =
        std::accumulate(ndc.I_RLC.begin(), ndc.I_RLC.end(), std::complex{0.0, 0.0});
    CHECK_THAT(std::abs(ndc.I), Catch::Matchers::WithinRel(std::abs(sum_ndcI)));
    CHECK_THAT(std::arg(ndc.I), Catch::Matchers::WithinRel(std::arg(sum_ndcI)));

    CHECK_THAT(c.inductor_energy, Catch::Matchers::WithinRel(ndc.inductor_energy));
    CHECK_THAT(c.capacitor_energy, Catch::Matchers::WithinRel(ndc.capacitor_energy));
    CHECK_THAT(c.mode_port_kappa, Catch::Matchers::WithinRel(ndc.mode_port_kappa));
    CHECK_THAT(c.inductive_energy_participation,
               Catch::Matchers::WithinRel(ndc.inductive_energy_participation));

    // Dimensional
    REQUIRE_THAT(units.GetScaleFactor<Units::ValueType::POWER>(),
                 Catch::Matchers::WithinRel(1.0));
    CHECK_THAT(std::abs(c.P),
               Catch::Matchers::WithinRel(std::abs(dc.P)));  // Power always unit normalized
    CHECK_THAT(std::arg(c.P), Catch::Matchers::WithinRel(
                                  std::arg(dc.P)));  // Phase unchanged by normalization
    CHECK_THAT(std::abs(c.S),
               Catch::Matchers::WithinRel(std::abs(dc.S)));  // Scattering always non-dim
    CHECK_THAT(std::arg(c.S), Catch::Matchers::WithinRel(
                                  std::arg(dc.S)));  // Phase unchanged by normalization
    if (std::abs(c.V) > 0)                           // lumped
    {
      CHECK_THAT(std::abs(c.V), !Catch::Matchers::WithinRel(std::abs(dc.V)));
      CHECK_THAT(std::arg(c.V), Catch::Matchers::WithinRel(std::arg(dc.V)));
      CHECK_THAT(std::abs(c.I), !Catch::Matchers::WithinRel(std::abs(dc.I)));
      CHECK_THAT(std::arg(c.I), Catch::Matchers::WithinRel(std::arg(dc.I)));
      for (auto i : {0, 1, 2})
      {
        CHECK_THAT(std::abs(c.I_RLC[i]),
                   !Catch::Matchers::WithinRel(std::abs(dc.I_RLC[i])));
        CHECK_THAT(std::arg(c.I_RLC[i]), Catch::Matchers::WithinRel(std::arg(dc.I_RLC[i])));
      }
      auto sum_dcI =
          std::accumulate(dc.I_RLC.begin(), dc.I_RLC.end(), std::complex{0.0, 0.0});
      CHECK_THAT(std::abs(dc.I), Catch::Matchers::WithinRel(std::abs(sum_dcI)));
      CHECK_THAT(std::arg(dc.I), Catch::Matchers::WithinRel(std::arg(sum_dcI)));

      // Power voltage current relations are stable through unit conversion
      auto ndcP = ndc.V * std::conj(ndc.I);
      CHECK_THAT(std::abs(ndc.P), Catch::Matchers::WithinRel(std::abs(ndcP)));
      CHECK_THAT(std::arg(ndc.P), Catch::Matchers::WithinRel(std::arg(ndcP)));

      auto dcP = dc.V * std::conj(dc.I);
      CHECK_THAT(std::abs(dc.P), Catch::Matchers::WithinRel(std::abs(dcP)));
      CHECK_THAT(std::arg(dc.P), Catch::Matchers::WithinRel(std::arg(dcP)));

      CHECK_THAT(c.inductor_energy, !Catch::Matchers::WithinRel(dc.inductor_energy));
      CHECK_THAT(c.capacitor_energy, !Catch::Matchers::WithinRel(dc.capacitor_energy));
      CHECK_THAT(c.mode_port_kappa, !Catch::Matchers::WithinRel(dc.mode_port_kappa));
      CHECK_THAT(c.inductive_energy_participation,
                 Catch::Matchers::WithinRel(dc.inductive_energy_participation));
    }
  };
  for (const auto &[k, v] : cache.lumped_port_vi)
  {
    check_port_data(cache.lumped_port_vi[k], dim_cache.lumped_port_vi[k],
                    non_dim_cache.lumped_port_vi[k]);
  }
  for (const auto &[k, v] : cache.wave_port_vi)
  {
    check_port_data(cache.wave_port_vi[k], dim_cache.wave_port_vi[k],
                    non_dim_cache.wave_port_vi[k]);
  }

  for (std::size_t i = 0; i < cache.probe_E_field.size(); i++)
  {
    const auto &cE = cache.probe_E_field[i];
    const auto &ndcE = non_dim_cache.probe_E_field[i];
    CHECK_THAT(std::abs(cE), Catch::Matchers::WithinRel(std::abs(ndcE)));
    CHECK_THAT(std::arg(cE), Catch::Matchers::WithinRel(std::arg(ndcE)));

    const auto &dcE = dim_cache.probe_E_field[i];
    CHECK_THAT(std::abs(cE), !Catch::Matchers::WithinRel(std::abs(dcE)));
    CHECK_THAT(std::arg(cE), Catch::Matchers::WithinRel(std::arg(dcE)));

    const auto &cB = cache.probe_B_field[i];
    const auto &ndcB = non_dim_cache.probe_B_field[i];
    CHECK_THAT(std::abs(cB), Catch::Matchers::WithinRel(std::abs(ndcB)));
    CHECK_THAT(std::arg(cB), Catch::Matchers::WithinRel(std::arg(ndcB)));

    const auto &dcB = dim_cache.probe_B_field[i];
    CHECK_THAT(std::abs(cB), !Catch::Matchers::WithinRel(std::abs(dcB)));
    CHECK_THAT(std::arg(cB), Catch::Matchers::WithinRel(std::arg(dcB)));

    // Compute Poynting vector H x E^ᴴ, W/m^2
    auto cP = cE * std::conj(cB);
    auto dcP = dcE * std::conj(dcB) / electromagnetics::mu0_;
    auto l_c = units.GetScaleFactor<Units::ValueType::LENGTH>();
    CHECK_THAT(std::abs(cP), Catch::Matchers::WithinRel(std::abs(dcP) * l_c * l_c));
    CHECK_THAT(std::arg(cP), Catch::Matchers::WithinRel(std::arg(dcP)));
  }

  for (std::size_t i = 0; i < cache.surface_flux_i.size(); i++)
  {
    auto &c = cache.surface_flux_i[i];
    auto &ndc = non_dim_cache.surface_flux_i[i];
    CHECK(c.idx == ndc.idx);
    CHECK(c.type == ndc.type);
    CHECK_THAT(std::abs(c.Phi), Catch::Matchers::WithinRel(std::abs(ndc.Phi)));
    CHECK_THAT(std::arg(c.Phi), Catch::Matchers::WithinRel(std::arg(ndc.Phi)));

    auto &dc = dim_cache.surface_flux_i[i];
    CHECK(c.idx == dc.idx);
    CHECK(c.type == dc.type);
    if (c.type == SurfaceFlux::POWER)
    {
      CHECK_THAT(std::abs(c.Phi), Catch::Matchers::WithinRel(std::abs(dc.Phi)));
    }
    else
    {
      CHECK_THAT(std::abs(c.Phi), !Catch::Matchers::WithinRel(std::abs(dc.Phi)));
    }
    CHECK_THAT(std::arg(c.Phi), Catch::Matchers::WithinRel(std::arg(dc.Phi)));
  }

  for (std::size_t i = 0; i < cache.interface_eps_i.size(); i++)
  {
    auto &c = cache.interface_eps_i[i];
    auto &ndc = non_dim_cache.interface_eps_i[i];
    CHECK(c.idx == ndc.idx);
    CHECK_THAT(c.energy, Catch::Matchers::WithinRel(ndc.energy));
    CHECK_THAT(c.tandelta, Catch::Matchers::WithinRel(ndc.tandelta));
    CHECK_THAT(c.energy_participation,
               Catch::Matchers::WithinRel(ndc.energy_participation));
    CHECK_THAT(c.quality_factor, Catch::Matchers::WithinRel(ndc.quality_factor));

    auto &dc = dim_cache.interface_eps_i[i];
    CHECK(c.idx == dc.idx);
    CHECK_THAT(c.energy, !Catch::Matchers::WithinRel(dc.energy));
    CHECK_THAT(c.tandelta, Catch::Matchers::WithinRel(dc.tandelta));
    CHECK_THAT(c.energy_participation, Catch::Matchers::WithinRel(dc.energy_participation));
    CHECK_THAT(c.quality_factor, Catch::Matchers::WithinRel(dc.quality_factor));
  }

  for (std::size_t i = 0; i < cache.interface_edge_i.size(); i++)
  {
    auto &c = cache.interface_edge_i[i];
    auto &ndc = non_dim_cache.interface_edge_i[i];
    CHECK(c.idx == ndc.idx);
    CHECK_THAT(c.distance, Catch::Matchers::WithinRel(ndc.distance));
    CHECK_THAT(c.energy_outside, Catch::Matchers::WithinRel(ndc.energy_outside));
    CHECK_THAT(c.participation_outside,
               Catch::Matchers::WithinRel(ndc.participation_outside));
    CHECK_THAT(c.energy_annulus, Catch::Matchers::WithinRel(ndc.energy_annulus));
    CHECK_THAT(c.participation_annulus,
               Catch::Matchers::WithinRel(ndc.participation_annulus));

    auto &dc = dim_cache.interface_edge_i[i];
    CHECK(c.idx == dc.idx);
    CHECK_THAT(c.distance, !Catch::Matchers::WithinRel(dc.distance));
    CHECK_THAT(c.energy_outside, !Catch::Matchers::WithinRel(dc.energy_outside));
    CHECK_THAT(c.participation_outside,
               Catch::Matchers::WithinRel(dc.participation_outside));
    CHECK_THAT(c.energy_annulus, !Catch::Matchers::WithinRel(dc.energy_annulus));
    CHECK_THAT(c.participation_annulus,
               Catch::Matchers::WithinRel(dc.participation_annulus));
  }

  for (std::size_t i = 0; i < cache.interface_local_edge_i.size(); i++)
  {
    auto &c = cache.interface_local_edge_i[i];
    auto &ndc = non_dim_cache.interface_local_edge_i[i];
    CHECK(c.idx == ndc.idx);
    CHECK(c.edge == ndc.edge);
    CHECK(c.automatic == ndc.automatic);
    CHECK(c.component == ndc.component);
    CHECK(c.chain == ndc.chain);
    CHECK(c.vertex_types == ndc.vertex_types);
    CHECK(c.process_normal == ndc.process_normal);
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(c.p0[d], Catch::Matchers::WithinRel(ndc.p0[d]));
      CHECK_THAT(c.p1[d], Catch::Matchers::WithinRel(ndc.p1[d]));
    }
    CHECK_THAT(c.length, Catch::Matchers::WithinRel(ndc.length));
    CHECK_THAT(c.distance, Catch::Matchers::WithinRel(ndc.distance));
    CHECK_THAT(c.energy_total, Catch::Matchers::WithinRel(ndc.energy_total));
    CHECK_THAT(c.participation_total, Catch::Matchers::WithinRel(ndc.participation_total));
    CHECK_THAT(c.energy_inside, Catch::Matchers::WithinRel(ndc.energy_inside));
    CHECK_THAT(c.participation_inside,
               Catch::Matchers::WithinRel(ndc.participation_inside));
    CHECK_THAT(c.energy_annulus, Catch::Matchers::WithinRel(ndc.energy_annulus));
    CHECK_THAT(c.participation_annulus,
               Catch::Matchers::WithinRel(ndc.participation_annulus));
    CHECK_THAT(c.energy_vertex_inside,
               Catch::Matchers::WithinRel(ndc.energy_vertex_inside));
    CHECK_THAT(c.participation_vertex_inside,
               Catch::Matchers::WithinRel(ndc.participation_vertex_inside));
    CHECK_THAT(c.energy_volume_annulus,
               Catch::Matchers::WithinRel(ndc.energy_volume_annulus));
    CHECK_THAT(c.participation_volume_annulus,
               Catch::Matchers::WithinRel(ndc.participation_volume_annulus));
    CHECK_THAT(c.energy_volume_vertex_annulus,
               Catch::Matchers::WithinRel(ndc.energy_volume_vertex_annulus));
    CHECK_THAT(c.participation_volume_vertex_annulus,
               Catch::Matchers::WithinRel(ndc.participation_volume_vertex_annulus));
    for (std::size_t component = 0; component < 2; component++)
    {
      CHECK_THAT(c.energy_total_polarized[component],
                 Catch::Matchers::WithinRel(ndc.energy_total_polarized[component]));
      CHECK_THAT(c.participation_total_polarized[component],
                 Catch::Matchers::WithinRel(ndc.participation_total_polarized[component]));
      CHECK_THAT(c.energy_inside_polarized[component],
                 Catch::Matchers::WithinRel(ndc.energy_inside_polarized[component]));
      CHECK_THAT(c.participation_inside_polarized[component],
                 Catch::Matchers::WithinRel(ndc.participation_inside_polarized[component]));
      CHECK_THAT(c.energy_annulus_polarized[component],
                 Catch::Matchers::WithinRel(ndc.energy_annulus_polarized[component]));
      CHECK_THAT(
          c.participation_annulus_polarized[component],
          Catch::Matchers::WithinRel(ndc.participation_annulus_polarized[component]));
    }
    for (std::size_t component = 0; component < 6; component++)
    {
      CHECK_THAT(
          c.energy_volume_annulus_polarized[component],
          Catch::Matchers::WithinRel(ndc.energy_volume_annulus_polarized[component]));
      CHECK_THAT(c.participation_volume_annulus_polarized[component],
                 Catch::Matchers::WithinRel(
                     ndc.participation_volume_annulus_polarized[component]));
    }

    auto &dc = dim_cache.interface_local_edge_i[i];
    CHECK(c.idx == dc.idx);
    CHECK(c.edge == dc.edge);
    CHECK(c.automatic == dc.automatic);
    CHECK(c.component == dc.component);
    CHECK(c.chain == dc.chain);
    CHECK(c.vertex_types == dc.vertex_types);
    CHECK(c.process_normal == dc.process_normal);
    for (int d = 0; d < 3; d++)
    {
      CHECK_THAT(c.p0[d], !Catch::Matchers::WithinRel(dc.p0[d]));
      CHECK_THAT(c.p1[d], !Catch::Matchers::WithinRel(dc.p1[d]));
    }
    CHECK_THAT(c.length, !Catch::Matchers::WithinRel(dc.length));
    CHECK_THAT(c.distance, !Catch::Matchers::WithinRel(dc.distance));
    CHECK_THAT(c.energy_total, !Catch::Matchers::WithinRel(dc.energy_total));
    CHECK_THAT(c.participation_total, Catch::Matchers::WithinRel(dc.participation_total));
    CHECK_THAT(c.energy_inside, !Catch::Matchers::WithinRel(dc.energy_inside));
    CHECK_THAT(c.participation_inside, Catch::Matchers::WithinRel(dc.participation_inside));
    CHECK_THAT(c.energy_annulus, !Catch::Matchers::WithinRel(dc.energy_annulus));
    CHECK_THAT(c.participation_annulus,
               Catch::Matchers::WithinRel(dc.participation_annulus));
    CHECK_THAT(c.energy_vertex_inside,
               !Catch::Matchers::WithinRel(dc.energy_vertex_inside));
    CHECK_THAT(c.participation_vertex_inside,
               Catch::Matchers::WithinRel(dc.participation_vertex_inside));
    CHECK_THAT(c.energy_volume_annulus,
               !Catch::Matchers::WithinRel(dc.energy_volume_annulus));
    CHECK_THAT(c.participation_volume_annulus,
               Catch::Matchers::WithinRel(dc.participation_volume_annulus));
    CHECK_THAT(c.energy_volume_vertex_annulus,
               !Catch::Matchers::WithinRel(dc.energy_volume_vertex_annulus));
    CHECK_THAT(c.participation_volume_vertex_annulus,
               Catch::Matchers::WithinRel(dc.participation_volume_vertex_annulus));
    for (std::size_t component = 0; component < 2; component++)
    {
      CHECK_THAT(c.energy_total_polarized[component],
                 !Catch::Matchers::WithinRel(dc.energy_total_polarized[component]));
      CHECK_THAT(c.participation_total_polarized[component],
                 Catch::Matchers::WithinRel(dc.participation_total_polarized[component]));
      CHECK_THAT(c.energy_inside_polarized[component],
                 !Catch::Matchers::WithinRel(dc.energy_inside_polarized[component]));
      CHECK_THAT(c.participation_inside_polarized[component],
                 Catch::Matchers::WithinRel(dc.participation_inside_polarized[component]));
      CHECK_THAT(c.energy_annulus_polarized[component],
                 !Catch::Matchers::WithinRel(dc.energy_annulus_polarized[component]));
      CHECK_THAT(c.participation_annulus_polarized[component],
                 Catch::Matchers::WithinRel(dc.participation_annulus_polarized[component]));
    }
    for (std::size_t component = 0; component < 6; component++)
    {
      CHECK_THAT(
          c.energy_volume_annulus_polarized[component],
          !Catch::Matchers::WithinRel(dc.energy_volume_annulus_polarized[component]));
      CHECK_THAT(
          c.participation_volume_annulus_polarized[component],
          Catch::Matchers::WithinRel(dc.participation_volume_annulus_polarized[component]));
    }
  }

  REQUIRE(dim_cache.farfield.E_field.size() == cache.farfield.E_field.size());
  REQUIRE(non_dim_cache.farfield.E_field.size() == cache.farfield.E_field.size());
  for (std::size_t i = 0; i < cache.farfield.E_field.size(); i++)
  {
    CHECK(cache.farfield.thetaphis[i] == dim_cache.farfield.thetaphis[i]);
    CHECK(cache.farfield.thetaphis[i] == non_dim_cache.farfield.thetaphis[i]);
    for (std::size_t d = 0; d < 3; d++)
    {
      const auto &c = cache.farfield.E_field[i][d];
      const auto &ndc = non_dim_cache.farfield.E_field[i][d];
      CHECK_THAT(std::abs(c), Catch::Matchers::WithinRel(std::abs(ndc)));
      CHECK_THAT(std::arg(c), Catch::Matchers::WithinRel(std::arg(ndc)));

      // The far field r E_∞ is a field times a distance: it must dimensionalize like the
      // probe E field (V/m) times the length scale (m), i.e. as a voltage.
      const auto &dc = dim_cache.farfield.E_field[i][d];
      const auto e_scale = units.GetScaleFactor<Units::ValueType::FIELD_E>();
      const auto l_scale = units.GetScaleFactor<Units::ValueType::LENGTH>();
      CHECK_THAT(std::abs(dc), Catch::Matchers::WithinRel(std::abs(c) * e_scale * l_scale));
      CHECK_THAT(std::arg(dc), Catch::Matchers::WithinRel(std::arg(c)));
    }
  }

  // Mode analysis data: kn [1/m], V [V], and Z [Ω] dimensionalize; n_eff is dimensionless.
  const auto l_scale = units.GetScaleFactor<Units::ValueType::LENGTH>();
  const auto v_scale = units.GetScaleFactor<Units::ValueType::VOLTAGE>();
  const auto z_scale = units.GetScaleFactor<Units::ValueType::IMPEDANCE>();
  const auto &m = cache.mode_data;
  const auto &dm = dim_cache.mode_data;
  const auto &ndm = non_dim_cache.mode_data;
  CHECK_THAT(std::abs(m.kn), Catch::Matchers::WithinRel(std::abs(ndm.kn)));
  CHECK_THAT(std::arg(m.kn), Catch::Matchers::WithinRel(std::arg(ndm.kn)));
  CHECK_THAT(std::abs(dm.kn), Catch::Matchers::WithinRel(std::abs(m.kn) / l_scale));
  CHECK_THAT(std::arg(dm.kn), Catch::Matchers::WithinRel(std::arg(m.kn)));
  CHECK(m.n_eff == ndm.n_eff);
  CHECK(m.n_eff == dm.n_eff);
  REQUIRE(ndm.voltage.size() == m.voltage.size());
  for (const auto &[idx, v] : m.voltage)
  {
    CHECK_THAT(std::abs(v.V), Catch::Matchers::WithinRel(std::abs(ndm.voltage.at(idx).V)));
    CHECK_THAT(std::arg(v.V), Catch::Matchers::WithinRel(std::arg(ndm.voltage.at(idx).V)));
    CHECK_THAT(std::abs(dm.voltage.at(idx).V),
               Catch::Matchers::WithinRel(std::abs(v.V) * v_scale));
    CHECK_THAT(std::arg(dm.voltage.at(idx).V), Catch::Matchers::WithinRel(std::arg(v.V)));
  }
  REQUIRE(ndm.impedance.size() == m.impedance.size());
  for (const auto &[idx, z] : m.impedance)
  {
    CHECK_THAT(std::abs(ndm.impedance.at(idx).Z0 - z.Z0),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12 * std::abs(z.Z0)));
    CHECK_THAT(std::abs(ndm.impedance.at(idx).Z_VI - z.Z_VI),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12 * std::abs(z.Z_VI)));
    const auto Z0 = z.has_impedance ? z.Z0 * z_scale : z.Z0;
    const auto Z_VI = z.has_vi_impedance ? z.Z_VI * z_scale : z.Z_VI;
    CHECK_THAT(std::abs(dm.impedance.at(idx).Z0 - Z0),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12 * std::abs(Z0)));
    CHECK_THAT(std::abs(dm.impedance.at(idx).Z_VI - Z_VI),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12 * std::abs(Z_VI)));
  }
}

TEST_CASE_METHOD(test::SharedTempDir, "Field export",
                 "[gridfunction][paraview][Serial][Parallel]")
{
  // Create iodata.
  Units units(0.496, 1.453);
  IoData iodata = IoData(units);

  iodata.domains.materials.emplace_back().attributes = {1};
  iodata.boundaries.pec.attributes = {1};
  iodata.problem.output_formats.gridfunction = true;
  // Direct output into a shared temporary directory. The default Problem.Output
  // ("postpro") is relative to the working directory, which is not guaranteed to be
  // writable during the test run.
  iodata.problem.output = temp_dir.string();
  iodata.CheckConfiguration();  // initializes quadrature

  // Setup lumped port boundary data for driven and transient.
  json boundaries = {{"LumpedPort",
                      {{{"Attributes", {2}},
                        {"Index", 1},
                        {"R", 50.0},
                        {"Direction", "+X"},
                        {"Excitation", true}}}}};
  config::BoundaryData boundary_port(boundaries);

  // Create serial mesh.
  int resolution = 3;
  std::unique_ptr<mfem::Mesh> serial_mesh =
      std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
          resolution, resolution, resolution, mfem::Element::TETRAHEDRON));

  // Create parallel mesh.
  auto comm = Mpi::World();
  int size = Mpi::Size(comm);
  iodata.model.Lc = mesh::ComputeReferenceLength(serial_mesh, comm);
  iodata.NondimensionalizeInputs(serial_mesh);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));

  auto check_files = [&](const std::string &subdir, int step, int pad_digits,
                         const std::vector<std::string> &fields)
  {
    // Root-owned metadata and every rank's payload must be complete before any rank
    // inspects the shared output tree.
    Mpi::Barrier(comm);
    for (int i = 0; i < size; i++)
    {
      for (const auto &field : fields)
      {
        auto path =
            fs::path(iodata.problem.output) / "gridfunction" / subdir /
            fmt::format("{}_{:0{}d}.gf.{:0{}d}", field, step, pad_digits, i, pad_digits);
        CHECK(fs::exists(path));
        CHECK(!fs::is_empty(path));
      }
    }
  };
  auto check_paraview_fields =
      [&](const std::string &subdir, const std::vector<std::string> &fields)
  {
    std::string metadata;
    const auto paraview_dir = fs::path(iodata.problem.output) / "paraview" / subdir;
    for (const auto &entry : fs::recursive_directory_iterator(paraview_dir))
    {
      const auto extension = entry.path().extension();
      if (!entry.is_regular_file() ||
          (extension != ".pvd" && extension != ".pvtu" && extension != ".vtu"))
      {
        continue;
      }
      std::ifstream file(entry.path(), std::ios::binary);
      metadata.append(std::istreambuf_iterator<char>(file), {});
    }
    for (const auto &field : fields)
    {
      CHECK(metadata.find(field) != std::string::npos);
    }
  };

  SECTION("Electrostatic")
  {
    // Create operator.
    iodata.problem.type = ProblemType::ELECTROSTATIC;
    iodata.solver.electrostatic.n_post = 1;
    LaplaceOperator laplace_op(iodata, mesh);
    PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);

    // Write fields.
    const int pad_digits = post_op.GetPadDigitsDefault();
    const auto &Grad = laplace_op.GetGradMatrix();
    Vector V(Grad.Width()), E(Grad.Height());
    V = 0.0;
    E = 0.0;
    post_op.MeasureAndPrintAll(0, V, E, 0);
    check_files("electrostatic", 0, post_op.GetPadDigitsDefault(), {"E", "V", "U_e"});
  }

  SECTION("Magnetostatic")
  {
    // Create operator.
    iodata.problem.type = ProblemType::MAGNETOSTATIC;
    iodata.solver.magnetostatic.n_post = 1;
    CurlCurlOperator curlcurl_op(iodata, mesh);
    PostOperator<ProblemType::MAGNETOSTATIC> post_op(iodata, curlcurl_op);

    // Write fields.
    const int pad_digits = post_op.GetPadDigitsDefault();
    const auto &Curl = curlcurl_op.GetCurlMatrix();
    Vector A(Curl.Width()), B(Curl.Height());
    A = 0.0;
    B = 0.0;
    post_op.MeasureAndPrintAll(0, A, B, 0);
    check_files("magnetostatic", 0, post_op.GetPadDigitsDefault(), {"B", "A", "U_m"});
  }

  SECTION("Transient")
  {
    // Create operator.
    iodata.problem.type = ProblemType::TRANSIENT;
    iodata.solver.transient.delta_post = 1;
    iodata.boundaries.lumpedport = boundary_port.lumpedport;
    SpaceOperator space_op(iodata, mesh);
    PostOperator<ProblemType::TRANSIENT> post_op(iodata, space_op);

    // Write fields.
    const int pad_digits = post_op.GetPadDigitsDefault();
    Vector E(space_op.GetNDSpace().GetTrueVSize()), B(space_op.GetRTSpace().GetTrueVSize());
    E = 0.0;
    B = 0.0;
    std::function<double(double)> J_coef = [](double x) -> double { return x; };
    post_op.MeasureAndPrintAll(0, E, B, 0.0, J_coef(0.0));
    check_files("transient", 0, post_op.GetPadDigitsDefault(),
                {"E", "B", "S", "U_e", "U_m"});
  }

  SECTION("Driven")
  {
    // Create operator.
    iodata.problem.type = ProblemType::DRIVEN;
    iodata.solver.driven.save_indices = {0};
    iodata.solver.driven.sample_f = {1.0};
    iodata.boundaries.lumpedport = boundary_port.lumpedport;
    SpaceOperator space_op(iodata, mesh);
    PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);

    // Write fields.
    const int pad_digits = post_op.GetPadDigitsDefault();
    ComplexVector E(space_op.GetNDSpace().GetTrueVSize()),
        B(space_op.GetRTSpace().GetTrueVSize()),
        E_corrected(space_op.GetNDSpace().GetTrueVSize()),
        D_corrected(space_op.GetRTSpace().GetTrueVSize());
    E = 0.0;
    B = 0.0;
    E_corrected = 1.0;
    D_corrected = 2.0;
    post_op.SetSurfaceResponseCorrectedField(E_corrected, &D_corrected);
    post_op.MeasureAndPrintAll(1, 0, E, B, 1.0);
    check_files("driven", 1, post_op.GetPadDigitsDefault(),
                {"E_real", "E_imag", "B_real", "B_imag", "E_corrected_real",
                 "E_corrected_imag", "D_corrected_real", "D_corrected_imag", "S", "U_e",
                 "U_m"});
    check_paraview_fields("driven", {"E_corrected_real", "E_corrected_imag",
                                     "D_corrected_real", "D_corrected_imag"});

    // The ParaView "Time" field for a driven solve is the physical frequency in
    // GHz, not the nondimensional angular frequency omega passed above.
    std::ifstream pvd(fs::path(iodata.problem.output) / "paraview" / "driven" /
                      "driven.pvd");
    std::string contents(std::istreambuf_iterator<char>(pvd), {});
    auto pos = contents.find("timestep=\"") + std::string("timestep=\"").size();
    double pvd_time = std::stod(contents.substr(pos, contents.find('"', pos) - pos));
    const double expected_freq =
        iodata.units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0) /
        (2.0 * std::numbers::pi);
    CHECK_THAT(pvd_time, Catch::Matchers::WithinRel(expected_freq, 1.0e-5));
  }

  SECTION("Driven multiple excitations")
  {
    // Create operator with two excitations.
    iodata.problem.type = ProblemType::DRIVEN;
    iodata.solver.driven.save_indices = {0};
    iodata.solver.driven.sample_f = {1.0};
    json multi_port_boundaries = {{"LumpedPort",
                                   {{{"Attributes", {2}},
                                     {"Index", 1},
                                     {"R", 50.0},
                                     {"Direction", "+X"},
                                     {"Excitation", 1}},
                                    {{"Attributes", {3}},
                                     {"Index", 2},
                                     {"R", 50.0},
                                     {"Direction", "+Y"},
                                     {"Excitation", 2}}}}};
    config::BoundaryData multi_port(multi_port_boundaries);
    iodata.boundaries.lumpedport = multi_port.lumpedport;
    SpaceOperator space_op(iodata, mesh);
    PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);

    // Use different fields for the two excitations, to tell their output apart.
    ComplexVector E(space_op.GetNDSpace().GetTrueVSize()),
        B(space_op.GetRTSpace().GetTrueVSize());
    E = 0.0;
    B = 0.0;
    post_op.InitializeParaviewDataCollection(1);
    post_op.MeasureAndPrintAll(1, 0, E, B, 1.0);
    E = 1.0;
    post_op.InitializeParaviewDataCollection(2);
    post_op.MeasureAndPrintAll(2, 0, E, B, 1.0);
    post_op.MeasureFinalize(ErrorIndicator());

    // Switching excitations must preserve both collections.
    auto root = fs::path(iodata.problem.output) / "paraview";
    for (int excitation : {1, 2})
    {
      auto folder = fmt::format("excitation_{}", excitation);
      CHECK(fs::is_regular_file(root / "driven" / folder / (folder + ".pvd")));
      CHECK(fs::is_regular_file(root / "driven_boundary" / folder / (folder + ".pvd")));
    }

    // Grid functions use the same excitation subfolders, while the mesh and the rank
    // partitioning are written once for all excitations.
    const int pad_digits = post_op.GetPadDigitsDefault();
    const std::vector<std::string> fields = {"E_real", "E_imag", "B_real", "B_imag",
                                             "S",      "U_e",    "U_m"};
    check_files("driven/excitation_1", 1, pad_digits, fields);
    check_files("driven/excitation_2", 1, pad_digits, fields);
    auto gf_root = fs::path(iodata.problem.output) / "gridfunction" / "driven";
    auto read_file = [](const fs::path &path)
    {
      std::ifstream file(path);
      return std::string(std::istreambuf_iterator<char>(file), {});
    };
    for (int i = 0; i < size; i++)
    {
      const auto rank = fmt::format("{:0{}d}", i, pad_digits);
      const auto E_real = fmt::format("E_real_{:0{}d}.gf.{}", 1, pad_digits, rank);
      CHECK(read_file(gf_root / "excitation_1" / E_real) !=
            read_file(gf_root / "excitation_2" / E_real));
      CHECK_FALSE(fs::exists(gf_root / E_real));
      CHECK(fs::is_regular_file(gf_root / ("mesh." + rank)));
      CHECK(fs::is_regular_file(gf_root / ("rank.gf." + rank)));
      CHECK_FALSE(fs::exists(gf_root / "excitation_1" / ("mesh." + rank)));
    }
  }

  SECTION("Driven reactive lumped port measurement")
  {
    // Purely reactive ports use the unit real reference resistance for excitation and
    // S-parameter normalization. Exercise the batched PostOperator measurement path,
    // which must agree with LumpedPortData::GetSParameter rather than special-case R = 0.
    json reactive_boundaries = {{"LumpedPort",
                                 {{{"Attributes", {2}},
                                   {"Index", 1},
                                   {"L", 1.0},
                                   {"Direction", "+X"},
                                   {"Excitation", true}}}}};
    config::BoundaryData reactive_port(reactive_boundaries);

    iodata.problem.type = ProblemType::DRIVEN;
    iodata.solver.driven.sample_f = {1.0};
    iodata.boundaries.lumpedport = reactive_port.lumpedport;
    SpaceOperator space_op(iodata, mesh);
    DrivenPostOperatorTest post_op(iodata, space_op);

    ComplexVector E(space_op.GetNDSpace().GetTrueVSize()),
        B(space_op.GetRTSpace().GetTrueVSize());
    E = 1.0;
    B = 0.0;
    const auto vi = post_op.MeasureLumpedPortForTest(1, E, B, 1.0);
    const auto &port = space_op.GetLumpedPortOp().GetPort(1);
    REQUIRE(std::abs(vi.V) > 0.0);
    CHECK_THAT(std::abs(vi.S - vi.V / std::sqrt(port.GetExcitationRefResistance())),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12 * std::abs(vi.V)));

    GridFunction E_gf(space_op.GetNDSpace(), true);
    E_gf.Real().SetFromTrueDofs(E.Real());
    E_gf.Imag().SetFromTrueDofs(E.Imag());
    const auto S_legacy = port.GetSParameter(E_gf);
    CHECK_THAT(std::abs(vi.S - S_legacy),
               Catch::Matchers::WithinAbs(0.0, 1.0e-12 * std::abs(S_legacy)));
  }

  SECTION("Eigenmode")
  {
    // Create operator.
    iodata.problem.type = ProblemType::EIGENMODE;
    iodata.solver.eigenmode.n_post = 1;
    SpaceOperator space_op(iodata, mesh);
    PostOperator<ProblemType::EIGENMODE> post_op(iodata, space_op);

    // Write fields.
    const int pad_digits = post_op.GetPadDigitsDefault();
    ComplexVector E(space_op.GetNDSpace().GetTrueVSize()),
        B(space_op.GetRTSpace().GetTrueVSize());
    E = 0.0;
    B = 0.0;
    constexpr int num_conv = 3;
    for (int i = 0; i < num_conv; i++)
    {
      post_op.MeasureAndPrintAll(i, E, B, 1.0 + i, 0.0, 0.0, num_conv);
    }
    check_files("eigenmode", 1, post_op.GetPadDigitsDefault(),
                {"E_real", "E_imag", "B_real", "B_imag", "S", "U_e", "U_m"});

    // Every converged mode is reported, while Save = 1 limits field data to the first
    // requested mode.
    std::ifstream eig_csv(fs::path(iodata.problem.output) / "eig.csv");
    REQUIRE(eig_csv.good());
    const int num_lines = static_cast<int>(std::count(
        std::istreambuf_iterator<char>(eig_csv), std::istreambuf_iterator<char>(), '\n'));
    CHECK(num_lines == num_conv + 1);  // Header plus one row per converged mode.
    for (int rank = 0; rank < size; rank++)
    {
      auto extra_field =
          fs::path(iodata.problem.output) / "gridfunction" / "eigenmode" /
          fmt::format("E_real_{:0{}d}.gf.{:0{}d}", 2, pad_digits, rank, pad_digits);
      CHECK_FALSE(fs::exists(extra_field));
    }
  }
}

TEST_CASE("Dimensional field output", "[postoperator][Serial][Parallel]")
{
  // The output of the fields redimensionalizes the mesh, undoes the effect of the mesh
  // scaling on the fields mapped from the reference element (ScaleGridFunctions), and
  // scales the fields to SI units (DimensionalizeGridFunctions). The field values at fixed
  // reference points, including in the face-neighbor elements owned by other processes,
  // must then be the nondimensional ones times the unit scale factors, independent of the
  // mesh length unit. The spaces are those of the PostOperator: E and A in H(curl), V in
  // H1, and B in H(div) in 3D or the scalar out-of-plane field in L2 with INTEGRAL map type
  // in 2D. The energy densities and Poynting vector projected from the dimensional
  // coefficients are checked the same way.
  constexpr double rtol = 1.0e-12;
  constexpr int order = 2;
  const int dim = GENERATE(2, 3);
  const bool simplex = GENERATE(true, false);
  CAPTURE(dim, simplex);

  // Create iodata.
  Units units(0.496, 1.453);
  IoData iodata = IoData(units);
  iodata.domains.materials.emplace_back().attributes = {1};

  // Create serial mesh.
  constexpr int resolution = 3;
  constexpr double length = 12.2;
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      (dim == 2) ? mfem::Mesh::MakeCartesian2D(resolution, resolution,
                                               simplex ? mfem::Element::TRIANGLE
                                                       : mfem::Element::QUADRILATERAL,
                                               true, length, length)
                 : mfem::Mesh::MakeCartesian3D(resolution, resolution, resolution,
                                               simplex ? mfem::Element::TETRAHEDRON
                                                       : mfem::Element::HEXAHEDRON,
                                               length, length, length));
  serial_mesh->EnsureNodes();

  // Create parallel mesh.
  auto comm = Mpi::World();
  if (Mpi::Size(comm) > serial_mesh->GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  iodata.model.Lc = mesh::ComputeReferenceLength(serial_mesh, comm);
  iodata.NondimensionalizeInputs(serial_mesh);
  Mesh mesh(std::make_unique<mfem::ParMesh>(comm, *serial_mesh));
  auto &pmesh = mesh.Get();
  const auto mesh_Lc0 = iodata.units.GetMeshLengthRelativeScale();
  CAPTURE(mesh_Lc0);
  REQUIRE(std::abs(mesh_Lc0 - 1.0) > 0.1);
  if (Mpi::Size(comm) > 1)
  {
    // Each process has elements on other processes neighboring its own.
    int num_face_nbr = pmesh.GetNFaceNeighborElements();
    Mpi::GlobalMin(1, &num_face_nbr, comm);
    REQUIRE(num_face_nbr > 0);
  }

  // Create gridfunctions.
  mfem::H1_FECollection h1_fec(order, dim);
  mfem::ND_FECollection nd_fec(order, dim);
  mfem::RT_FECollection rt_fec(order - 1, dim);
  std::unique_ptr<mfem::FiniteElementCollection> b_fec;
  if (dim == 2)
  {
    b_fec = std::make_unique<mfem::L2_FECollection>(
        order - 1, dim, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL);
  }
  else
  {
    b_fec = std::make_unique<mfem::RT_FECollection>(order - 1, dim);
  }
  mfem::ParFiniteElementSpace h1_fespace(&pmesh, &h1_fec), nd_fespace(&pmesh, &nd_fec),
      rt_fespace(&pmesh, &rt_fec), b_fespace(&pmesh, b_fec.get());
  mfem::ParGridFunction Ue_gf(&h1_fespace), Um_gf(&h1_fespace), S_gf(&rt_fespace);

  std::unique_ptr<GridFunction> E, B, V, A;
  E = std::make_unique<GridFunction>(nd_fespace, true);
  B = std::make_unique<GridFunction>(b_fespace, true);
  V = std::make_unique<GridFunction>(h1_fespace);
  A = std::make_unique<GridFunction>(nd_fespace);

  // Initialize the fields to random values, and make them ready for evaluation on shared
  // faces as in PostOperator.
  int seed = 1;
  for (auto *gf : {&E->Real(), &E->Imag(), &B->Real(), &B->Imag(), &V->Real(), &A->Real()})
  {
    Vector x(gf->ParFESpace()->GetTrueVSize());
    linalg::SetRandom(comm, x, seed++);
    gf->SetFromTrueDofs(x);
    gf->ExchangeFaceNbrData();
  }

  // Create some coefficients.
  MaterialOperator mat_op(iodata, mesh);
  std::unique_ptr<mfem::VectorCoefficient> S_nondim, S_dim;
  std::unique_ptr<mfem::Coefficient> Ue_nondim, Ue_dim, Um_nondim, Um_dim;

  // Omit scaling argument (e.g. scaling = 1.0) for non-dimensional coefficients.
  Ue_nondim =
      std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC>>(*E, mat_op);
  Um_nondim =
      std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC>>(*B, mat_op);
  S_nondim = std::make_unique<PoyntingVectorCoefficient>(*E, *B, mat_op);

  // Use same scaling as in postoperator.cpp for dimensional coefficients.
  auto scaling = iodata.units.Dimensionalize<Units::ValueType::FIELD_D>(1.0) /
                 iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
  Ue_dim = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::ELECTRIC>>(
      *E, mat_op, scaling);
  scaling = iodata.units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) /
            iodata.units.Dimensionalize<Units::ValueType::FIELD_B>(1.0);
  Um_dim = std::make_unique<EnergyDensityCoefficient<EnergyDensityType::MAGNETIC>>(
      *B, mat_op, scaling);
  S_dim = std::make_unique<PoyntingVectorCoefficient>(*E, *B, mat_op, scaling);

  auto update_coeff_gridfunctions = [&](bool dimensional)
  {
    Ue_gf.ProjectCoefficient(dimensional ? *Ue_dim : *Ue_nondim);
    Um_gf.ProjectCoefficient(dimensional ? *Um_dim : *Um_nondim);
    S_gf.ProjectCoefficient(dimensional ? *S_dim : *S_nondim);
  };

  // Fields with their units scale factors. The coefficient grid functions are only
  // evaluated in local elements, since their face-neighbor data is not exchanged.
  struct Field
  {
    std::string name;
    const mfem::ParGridFunction *gf;
    double scale;
    bool face_nbr;
  };
  const auto scale_E = iodata.units.GetScaleFactor<Units::ValueType::FIELD_E>();
  const auto scale_B = iodata.units.GetScaleFactor<Units::ValueType::FIELD_B>();
  const auto scale_D = iodata.units.GetScaleFactor<Units::ValueType::FIELD_D>();
  const auto scale_H = iodata.units.GetScaleFactor<Units::ValueType::FIELD_H>();
  const std::vector<Field> fields = {
      {"E_real", &E->Real(), scale_E, true},
      {"E_imag", &E->Imag(), scale_E, true},
      {"B_real", &B->Real(), scale_B, true},
      {"B_imag", &B->Imag(), scale_B, true},
      {"V", &V->Real(), iodata.units.GetScaleFactor<Units::ValueType::VOLTAGE>(), true},
      {"A", &A->Real(), iodata.units.GetScaleFactor<Units::ValueType::CURRENT>(), true},
      {"U_e", &Ue_gf, scale_D * scale_E, false},
      {"U_m", &Um_gf, scale_H * scale_B, false},
      {"S", &S_gf, scale_E * scale_H, false}};

  // Evaluate the fields at a set of reference points in every local and face-neighbor
  // element.
  auto Sample = [&]()
  {
    std::vector<std::vector<double>> values(fields.size());
    mfem::IsoparametricTransformation T_nbr;
    mfem::Vector val;
    auto SampleElement =
        [&](mfem::ElementTransformation &T, mfem::Geometry::Type geom, bool face_nbr)
    {
      const auto &ir = mfem::IntRules.Get(geom, 2 * order);
      for (int q = 0; q < ir.GetNPoints(); q++)
      {
        const auto &ip = ir.IntPoint(q);
        T.SetIntPoint(&ip);
        for (std::size_t i = 0; i < fields.size(); i++)
        {
          if (face_nbr && !fields[i].face_nbr)
          {
            continue;
          }
          fields[i].gf->GetVectorValue(T, ip, val);
          values[i].insert(values[i].end(), val.begin(), val.end());
        }
      }
    };
    for (int e = 0; e < pmesh.GetNE(); e++)
    {
      SampleElement(*pmesh.GetElementTransformation(e), pmesh.GetElementGeometry(e), false);
    }
    for (int e = 0; e < pmesh.GetNFaceNeighborElements(); e++)
    {
      pmesh.GetFaceNbrElementTransformation(e, T_nbr);
      SampleElement(T_nbr, pmesh.face_nbr_elements[e]->GetGeometryType(), true);
    }
    return values;
  };

  // Check that values = scale * ref for each field, with scale its units scale factor if
  // dimensional and one otherwise.
  auto Compare = [&](const std::vector<std::vector<double>> &ref,
                     const std::vector<std::vector<double>> &values, bool dimensional)
  {
    for (std::size_t i = 0; i < fields.size(); i++)
    {
      REQUIRE(ref[i].size() == values[i].size());
      const double scale = dimensional ? fields[i].scale : 1.0;
      double max_ref = 0.0, max_diff = 0.0;
      for (std::size_t j = 0; j < ref[i].size(); j++)
      {
        max_ref = std::max(max_ref, std::abs(scale * ref[i][j]));
        max_diff = std::max(max_diff, std::abs(values[i][j] - scale * ref[i][j]));
      }
      Mpi::GlobalMax(1, &max_ref, comm);
      Mpi::GlobalMax(1, &max_diff, comm);
      CAPTURE(fields[i].name, dimensional, max_ref, max_diff);
      CHECK(max_ref > 0.0);
      CHECK(max_diff <= rtol * max_ref);
    }
  };

  // Sample the fields before dimensionalizing.
  update_coeff_gridfunctions(false);
  const auto values_before = Sample();

  // Dimensionalize.
  mesh::DimensionalizeMesh(pmesh, mesh_Lc0);
  ScaleGridFunctions(mesh_Lc0, {E.get(), B.get(), V.get(), A.get()});
  DimensionalizeGridFunctions(iodata.units, E, B, V, A);
  update_coeff_gridfunctions(true);
  Compare(values_before, Sample(), true);

  // Nondimensionalize.
  mesh::NondimensionalizeMesh(pmesh, mesh_Lc0);
  ScaleGridFunctions(1.0 / mesh_Lc0, {E.get(), B.get(), V.get(), A.get()});
  NondimensionalizeGridFunctions(iodata.units, E, B, V, A);
  update_coeff_gridfunctions(false);
  Compare(values_before, Sample(), false);
}

TEST_CASE("GridFunction scaling keeps face-neighbor data consistent",
          "[postoperator][Parallel]")
{
  // In-place scaling of a GridFunction must scale its cached face-neighbor data together
  // with the field, so the cache stays consistent with the field until the next
  // ExchangeFaceNbrData. Check this for a real- and a complex-valued field, and for the
  // Units in-place (non)dimensionalization round trip.
  auto comm = Mpi::World();
  Units units(1.0e-3, 2.0e-3);

  // Mesh partitioned across ranks, so each rank caches face-neighbor data.
  constexpr int resolution = 3;
  auto serial_mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
      resolution, resolution, resolution, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0));
  mfem::ParMesh mesh(comm, *serial_mesh);
  mfem::ND_FECollection nd_fec(1, mesh.Dimension());
  mfem::ParFiniteElementSpace nd_fespace(&mesh, &nd_fec);

  // Check that values = scale * baseline, component-wise within relative tolerance.
  constexpr double rtol = 1.0e-13;
  auto CheckScaled =
      [](const Vector &baseline, const Vector &values, double scale, double tolerance)
  {
    REQUIRE(baseline.Size() == values.Size());
    double max_diff = 0.0, max_val = 0.0;
    for (int i = 0; i < baseline.Size(); i++)
    {
      max_diff = std::max(max_diff, std::abs(scale * baseline(i) - values(i)));
      max_val = std::max(max_val, std::abs(values(i)));
    }
    CHECK(max_diff <= tolerance * std::max(max_val, 1.0));
  };

  // Check that the cached face-neighbor data equals a fresh exchange of the current
  // field values on the neighboring ranks.
  auto CheckFaceNbrData = [&](GridFunction &gf)
  {
    Vector nbr_r = gf.Real().FaceNbrData();
    Vector nbr_i;
    if (gf.HasImag())
    {
      nbr_i = gf.Imag().FaceNbrData();
    }
    gf.Real().ExchangeFaceNbrData();
    if (gf.HasImag())
    {
      gf.Imag().ExchangeFaceNbrData();
    }
    CheckScaled(nbr_r, gf.Real().FaceNbrData(), 1.0, rtol);
    if (gf.HasImag())
    {
      CheckScaled(nbr_i, gf.Imag().FaceNbrData(), 1.0, rtol);
    }
  };

  auto CheckScaling = [&](GridFunction &gf, double scale)
  {
    // Populate the field through the true-dof distribution, then cache the face-neighbor
    // data and keep a baseline of the local field values.
    Vector x(nd_fespace.GetTrueVSize());
    linalg::SetRandom(comm, x, 1);
    gf.Real().SetFromTrueDofs(x);
    if (gf.HasImag())
    {
      linalg::SetRandom(comm, x, 2);
      gf.Imag().SetFromTrueDofs(x);
    }
    gf.Real().ExchangeFaceNbrData();
    if (gf.HasImag())
    {
      gf.Imag().ExchangeFaceNbrData();
    }
    REQUIRE(gf.Real().FaceNbrData().Size() > 0);
    Vector base_r = gf.Real();
    Vector base_i;
    if (gf.HasImag())
    {
      REQUIRE(gf.Imag().FaceNbrData().Size() > 0);
      base_i = gf.Imag();
    }

    // Scaling in place scales the cached face-neighbor data together with the field.
    gf *= scale;
    CheckScaled(base_r, gf.Real(), scale, rtol);
    if (gf.HasImag())
    {
      CheckScaled(base_i, gf.Imag(), scale, rtol);
    }
    CheckFaceNbrData(gf);

    // The Units in-place scaling round trip returns to the scaled baseline and also
    // keeps the cached face-neighbor data consistent.
    units.DimensionalizeInPlace<Units::ValueType::FIELD_E>(gf);
    units.NondimensionalizeInPlace<Units::ValueType::FIELD_E>(gf);
    CheckScaled(base_r, gf.Real(), scale, rtol);
    if (gf.HasImag())
    {
      CheckScaled(base_i, gf.Imag(), scale, rtol);
    }
    CheckFaceNbrData(gf);
  };

  constexpr double scale = 2.5;
  GridFunction gf_real(nd_fespace);
  CheckScaling(gf_real, scale);

  GridFunction gf_complex(nd_fespace, true);
  CheckScaling(gf_complex, scale);
}

namespace
{

// The FluxRecovery cases below (decision 1578; SD-605 / SD-606 / SD-508 / SD-504): the unit
// cube with substrate (element attribute 1, z < 0.5) under vacuum (2, z > 0.5), the
// exterior boundary attributes 1-6 of mfem::Mesh::MakeCartesian3D (1 = the face z = 0),
// optionally a sheet of attribute 7 on the interior plane x = 0.5 (not cracked: its faces
// are shared by the elements on both sides), and the interior substrate-air interface of
// attribute 8 on the plane z = 0.5 (no boundary condition there: never a sheet).
constexpr int flux_sheet_attr = 7, flux_interface_attr = 8;

std::unique_ptr<mfem::Mesh> MakeFluxRecoveryCubeMesh(bool interior_sheet)
{
  auto base =
      mfem::Mesh::MakeCartesian3D(4, 4, 4, mfem::Element::TETRAHEDRON, 1.0, 1.0, 1.0);
  std::vector<std::pair<int, int>> interior_faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    if (e2 < 0)
    {
      continue;
    }
    base.GetFaceVertices(f, fv);
    auto OnPlane = [&](int d)
    {
      return std::all_of(fv.begin(), fv.end(), [&](int v)
                         { return std::abs(base.GetVertex(v)[d] - 0.5) < 1.0e-12; });
    };
    if (interior_sheet && OnPlane(0))
    {
      interior_faces.emplace_back(f, flux_sheet_attr);
    }
    else if (OnPlane(2))
    {
      interior_faces.emplace_back(f, flux_interface_attr);
    }
  }
  auto mesh =
      std::make_unique<mfem::Mesh>(3, base.GetNV(), base.GetNE(),
                                   base.GetNBE() + static_cast<int>(interior_faces.size()));
  for (int v = 0; v < base.GetNV(); v++)
  {
    mesh->AddVertex(base.GetVertex(v));
  }
  mfem::Vector center(3);
  for (int e = 0; e < base.GetNE(); e++)
  {
    base.GetElementCenter(e, center);
    auto *el = base.GetElement(e)->Duplicate(mesh.get());
    el->SetAttribute((center(2) < 0.5) ? 1 : 2);
    mesh->AddElement(el);
  }
  for (int be = 0; be < base.GetNBE(); be++)
  {
    mesh->AddBdrElement(base.GetBdrElement(be)->Duplicate(mesh.get()));
  }
  for (const auto &[f, attr] : interior_faces)
  {
    auto *el = base.GetFace(f)->Duplicate(mesh.get());
    el->SetAttribute(attr);
    mesh->AddBdrElement(el);
  }
  mesh->FinalizeTopology();
  mesh->Finalize();
  mesh->SetAttributes();
  return mesh;
}

// The two-dimensional counterpart for boundary mode: the unit square with substrate
// (y < 0.5) under vacuum, the exterior boundary attributes 1-4 of MakeCartesian2D, an
// optional sheet of attribute 5 on the interior line x = 0.5 and the substrate-air
// interface of attribute 6 on the line y = 0.5.
std::unique_ptr<mfem::Mesh> MakeFluxRecoverySquareMesh(bool interior_sheet)
{
  auto mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian2D(4, 4, mfem::Element::TRIANGLE, false, 1.0, 1.0));
  mfem::Vector center(2);
  for (int e = 0; e < mesh->GetNE(); e++)
  {
    mesh->GetElementCenter(e, center);
    mesh->SetAttribute(e, (center(1) < 0.5) ? 1 : 2);
  }
  mfem::Array<int> fv;
  for (int f = 0; f < mesh->GetNumFaces(); f++)
  {
    int e1, e2;
    mesh->GetFaceElements(f, &e1, &e2);
    if (e2 < 0)
    {
      continue;
    }
    mesh->GetFaceVertices(f, fv);
    auto OnLine = [&](int d)
    {
      return std::all_of(fv.begin(), fv.end(), [&](int v)
                         { return std::abs(mesh->GetVertex(v)[d] - 0.5) < 1.0e-12; });
    };
    const int attr = (interior_sheet && OnLine(0)) ? 5 : (OnLine(1) ? 6 : 0);
    if (attr > 0)
    {
      auto *el = mesh->GetFace(f)->Duplicate(mesh.get());
      el->SetAttribute(attr);
      mesh->AddBdrElement(el);
    }
  }
  mesh->FinalizeTopology();
  mesh->Finalize();
  mesh->SetAttributes();
  return mesh;
}

// A configuration of the cube or square above with ONE interface dielectric entry (index 1)
// and the PEC boundaries given (+ the interior sheet's attribute when present).
IoData MakeFluxRecoveryIoData(ProblemType type, const std::vector<int> &pec,
                              InterfaceDielectric dielectric_type, int attr,
                              bool flux_recovery, const std::string &output)
{
  IoData iodata(Units(1.0, 1.0));
  iodata.problem.type = type;
  iodata.problem.verbose = 0;
  iodata.problem.output = output;
  iodata.problem.output_formats.paraview = false;
  iodata.problem.output_formats.gridfunction = false;
  iodata.model.L0 = 1.0;
  iodata.model.Lc = 1.0;
  iodata.solver.order = 1;
  iodata.solver.linear.tol = 1.0e-8;
  iodata.solver.linear.max_it = 100;
  iodata.solver.linear.estimator_tol = 1.0e-12;
  iodata.solver.linear.estimator_max_it = 1000;
  iodata.solver.boundary_mode.freq = 1.0;
  iodata.solver.boundary_mode.n = 1;
  iodata.solver.boundary_mode.n_post = 1;
  iodata.solver.transient.delta_post = 1;
  config::MaterialData substrate, vacuum;
  substrate.attributes = {1};
  substrate.epsilon_r.s = {4.0, 4.0, 4.0};
  vacuum.attributes = {2};
  iodata.domains.materials = {substrate, vacuum};
  iodata.boundaries.pec.attributes = pec;
  config::InterfaceDielectricData interface;
  interface.type = dielectric_type;
  interface.t = 2.0e-3;
  interface.epsilon_r = 10.0;
  interface.tandelta = 1.0e-3;
  interface.attributes = {attr};
  interface.flux_recovery = flux_recovery;
  iodata.boundaries.postpro.dielectric.emplace(1, interface);
  fem::DefaultIntegrationOrder::p_trial = 1;
  fem::DefaultIntegrationOrder::q_order_jac = true;
  fem::DefaultIntegrationOrder::q_order_extra_pk = 0;
  fem::DefaultIntegrationOrder::q_order_extra_qk = 0;
  return iodata;
}

std::vector<std::unique_ptr<Mesh>>
MakeFluxRecoveryParMesh(MPI_Comm comm, std::unique_ptr<mfem::Mesh> smesh)
{
  if (Mpi::Size(comm) > smesh->GetNE())
  {
    SKIP("Mesh has fewer elements than processes");
  }
  smesh->EnsureNodes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *smesh)));
  return mesh;
}

// The FluxRecovery entries of an electrostatic postprocessing on the cube of
// MakeFluxRecoveryCubeMesh without its sheet (no interior sheet face, so the recovery space
// is not broken, as on the cracked device meshes of the legacy key): 1 = SA on the
// substrate-air interface z = 0.5 with FluxRecovery, 2 = the same without (the raw field),
// 3 = MS on the exterior PEC face z = 0, 4 = MA on the exterior face z = 1 (vacuum side).
IoData FluxRecoveryNoSheetIoData(const std::string &output)
{
  IoData iodata(Units(1.0, 1.0));
  iodata.problem.type = ProblemType::ELECTROSTATIC;
  iodata.problem.verbose = 0;
  iodata.problem.output = output;
  iodata.problem.output_formats.paraview = false;
  iodata.problem.output_formats.gridfunction = false;
  iodata.model.L0 = 1.0;
  iodata.model.Lc = 1.0;
  iodata.solver.order = 2;
  iodata.solver.linear.estimator_tol = 1.0e-14;
  iodata.solver.linear.estimator_max_it = 10000;
  config::MaterialData substrate, vacuum;
  substrate.attributes = {1};
  substrate.epsilon_r.s = {4.0, 4.0, 4.0};
  vacuum.attributes = {2};
  iodata.domains.materials = {substrate, vacuum};
  iodata.boundaries.pec.attributes = {1};
  auto AddEntry = [&](int idx, InterfaceDielectric type, int attr, bool flux_recovery)
  {
    config::InterfaceDielectricData data;
    data.type = type;
    data.t = 2.0e-3;
    data.epsilon_r = 10.0;
    data.tandelta = 1.0e-3;
    data.attributes = {attr};
    data.flux_recovery = flux_recovery;
    iodata.boundaries.postpro.dielectric.emplace(idx, data);
  };
  AddEntry(1, InterfaceDielectric::SA, 8, true);
  AddEntry(2, InterfaceDielectric::SA, 8, false);
  AddEntry(3, InterfaceDielectric::MS, 1, true);
  AddEntry(4, InterfaceDielectric::MA, 6, true);
  fem::DefaultIntegrationOrder::p_trial = iodata.solver.order;
  fem::DefaultIntegrationOrder::q_order_jac = iodata.solver.q_order_jac;
  fem::DefaultIntegrationOrder::q_order_extra_pk = iodata.solver.q_order_extra;
  fem::DefaultIntegrationOrder::q_order_extra_qk = iodata.solver.q_order_extra;
  return iodata;
}

// The domain energy and the energies of the entries 1-4 for V = 1 - z + 0.2 x y and
// E = -grad V, with D recovered from E by the electrostatic error estimator as the
// electrostatic driver builds it.
std::vector<double> FluxRecoveryNoSheetEnergies(MPI_Comm comm, const std::string &output)
{
  IoData iodata = FluxRecoveryNoSheetIoData(output);
  auto smesh = MakeFluxRecoveryCubeMesh(false);
  smesh->EnsureNodes();
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::make_unique<mfem::ParMesh>(comm, *smesh)));
  LaplaceOperator laplace_op(iodata, mesh);
  PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);
  REQUIRE(post_op.NeedsRecoveredElectricFlux());
  mfem::FunctionCoefficient V_coeff([](const mfem::Vector &x)
                                    { return 1.0 - x(2) + 0.2 * x(0) * x(1); });
  mfem::VectorFunctionCoefficient E_coeff(3,
                                          [](const mfem::Vector &x, mfem::Vector &v)
                                          {
                                            v(0) = -0.2 * x(1);
                                            v(1) = -0.2 * x(0);
                                            v(2) = 1.0;
                                          });
  GridFunction V(laplace_op.GetH1Space()), E(laplace_op.GetNDSpace());
  V.Real().ProjectCoefficient(V_coeff);
  E.Real().ProjectCoefficient(E_coeff);
  Vector v(laplace_op.GetH1Space().GetTrueVSize()),
      e(laplace_op.GetNDSpace().GetTrueVSize());
  V.Real().GetTrueDofs(v);
  E.Real().GetTrueDofs(e);
  const auto &linear = iodata.solver.linear;
  GradFluxErrorEstimator<Vector> estimator(
      laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
      linear.estimator_tol, linear.estimator_max_it, 0, linear.estimator_mg,
      iodata.boundaries.GetSheetAttributes());
  REQUIRE(!estimator.GetFluxSpace().IsBroken());
  Vector d(estimator.GetFluxSpace().GetTrueVSize());
  estimator.RecoverFlux(e, d);
  const auto energies = post_op.GetElectrostaticEnergies(v, e, &d);
  std::vector<double> values = {energies.domain};
  for (int idx : {1, 2, 3, 4})
  {
    values.push_back(energies.interfaces.at(idx).energy);
  }
  return values;
}

// Measurement access for the FluxRecovery cases.
template <ProblemType solver_t>
class FluxRecoveryPostOperator : public PostOperator<solver_t>
{
public:
  using PostOperator<solver_t>::PostOperator;
  const Measurement &Cache() const { return this->measurement_cache; }
  const GridFunction &CorrectedFluxOutput() const
  {
    return *this->surface_response_corrected_flux;
  }
};

}  // namespace

TEST_CASE_METHOD(test::SharedTempDir,
                 "FluxRecovery refused before the solve on an interior sheet face",
                 "[postoperator][Serial][Parallel]")
{
  // An interface dielectric "FluxRecovery" is refused at the construction of the
  // postprocessing, before any solve, on an interior sheet face, cracked or not, by the
  // SurfacePostOperator (our method rule: no controlled normal trace of the volume
  // projection on a zero-thickness PEC sheet; SD-606 (1); kept by ruling 1600). Elsewhere
  // on a mesh with interior sheet faces it is accepted: the recovered flux lives on a
  // recovery space broken across them, which the PostOperator holds per side (P3-14,
  // SD-3202; the refusal of decision 1578,
  // OI-UPSYNC-FLUXRECOVERY-REFUSED-ON-BROKEN-RECOVERY-SPACE, is lifted). An exterior sheet
  // boundary keeps FluxRecovery, and without FluxRecovery nothing is refused. The sheet
  // faces cross the process boundaries in parallel (collective checks).
  auto comm = Mpi::World();
  const std::string sheet_refusal = "is not supported on internal boundary sheets";
  auto Electrostatic =
      [&](bool interior_sheet, InterfaceDielectric type, int attr, bool flux_recovery)
  {
    std::vector<int> pec = {1};
    if (interior_sheet)
    {
      pec.push_back(flux_sheet_attr);
    }
    IoData iodata = MakeFluxRecoveryIoData(ProblemType::ELECTROSTATIC, pec, type, attr,
                                           flux_recovery, temp_dir.string());
    auto mesh = MakeFluxRecoveryParMesh(comm, MakeFluxRecoveryCubeMesh(interior_sheet));
    LaplaceOperator laplace_op(iodata, mesh);
    PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);
    return post_op.NeedsRecoveredElectricFlux();
  };
  auto BoundaryMode = [&](bool interior_sheet, InterfaceDielectric type, int attr)
  {
    std::vector<int> pec = {1, 2, 3, 4};
    if (interior_sheet)
    {
      pec.push_back(5);
    }
    IoData iodata = MakeFluxRecoveryIoData(ProblemType::BOUNDARYMODE, pec, type, attr, true,
                                           temp_dir.string());
    auto mesh = MakeFluxRecoveryParMesh(comm, MakeFluxRecoverySquareMesh(interior_sheet));
    MaterialOperator mat_op(iodata, *mesh.front());
    BoundaryModeOperator mode_op(iodata, mesh, mat_op);
    PostOperator<ProblemType::BOUNDARYMODE> post_op(iodata, mode_op);
    return post_op.NeedsRecoveredElectricFlux();
  };

  SECTION("On the interior sheet")
  {
    CHECK_THROWS_WITH(Electrostatic(true, InterfaceDielectric::MA, flux_sheet_attr, true),
                      ContainsSubstring(sheet_refusal));
  }
  SECTION("Elsewhere on a mesh with an interior sheet")
  {
    bool needs_flux = false;
    CHECK_NOTHROW(needs_flux = Electrostatic(true, InterfaceDielectric::SA,
                                             flux_interface_attr, true));
    CHECK(needs_flux);
  }
  SECTION("Exterior sheet boundary")
  {
    bool needs_flux = false;
    CHECK_NOTHROW(needs_flux = Electrostatic(false, InterfaceDielectric::MS, 1, true));
    CHECK(needs_flux);
  }
  SECTION("Without FluxRecovery")
  {
    bool needs_flux = true;
    CHECK_NOTHROW(needs_flux =
                      Electrostatic(true, InterfaceDielectric::MA, flux_sheet_attr, false));
    CHECK(!needs_flux);
  }
  SECTION("Boundary mode")
  {
    CHECK_THROWS_WITH(BoundaryMode(true, InterfaceDielectric::MA, 5),
                      ContainsSubstring(sheet_refusal));
    bool needs_flux = false;
    CHECK_NOTHROW(needs_flux = BoundaryMode(true, InterfaceDielectric::SA, 6));
    CHECK(needs_flux);
    needs_flux = false;
    CHECK_NOTHROW(needs_flux = BoundaryMode(false, InterfaceDielectric::SA, 6));
    CHECK(needs_flux);
  }
}

TEST_CASE_METHOD(test::SharedTempDir, "Recovered flux size mismatch aborts",
                 "[postoperator][Serial][Parallel]")
{
  // The recovered electric flux D handed to the postprocessing must live on the space of
  // the estimators' flux recovery; a vector of another size is refused before it is read
  // (SD-502 / SD-605 backstop), for FluxRecovery and for the self-consistent corrected
  // flux.
  auto comm = Mpi::World();
  const std::string size_refusal = "does not match the flux recovery space";
  SECTION("Electrostatic FluxRecovery")
  {
    IoData iodata =
        MakeFluxRecoveryIoData(ProblemType::ELECTROSTATIC, {1}, InterfaceDielectric::MS, 1,
                               true, temp_dir.string());
    auto mesh = MakeFluxRecoveryParMesh(comm, MakeFluxRecoveryCubeMesh(false));
    LaplaceOperator laplace_op(iodata, mesh);
    PostOperator<ProblemType::ELECTROSTATIC> post_op(iodata, laplace_op);
    REQUIRE(post_op.NeedsRecoveredElectricFlux());
    GradFluxErrorEstimator<Vector> estimator(
        laplace_op.GetMaterialOp(), laplace_op.GetNDSpace(), laplace_op.GetRTSpaces(),
        iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0, false,
        iodata.boundaries.GetSheetAttributes());
    const int size = estimator.GetFluxSpace().GetTrueVSize();
    Vector d(size), d_wrong(size + 1);
    d = 1.0;
    d_wrong = 1.0;
    CHECK_THROWS_WITH(post_op.SetRecoveredElectricFlux(d_wrong),
                      ContainsSubstring(size_refusal));
    CHECK_NOTHROW(post_op.SetRecoveredElectricFlux(d));
  }
  SECTION("Driven corrected flux")
  {
    IoData iodata = MakeFluxRecoveryIoData(
        ProblemType::DRIVEN, {1}, InterfaceDielectric::MS, 1, false, temp_dir.string());
    iodata.solver.driven.sample_f = {1.0};
    config::BoundaryData port(json{{"LumpedPort",
                                    {{{"Attributes", {2}},
                                      {"Index", 1},
                                      {"R", 50.0},
                                      {"Direction", "+X"},
                                      {"Excitation", true}}}}});
    iodata.boundaries.lumpedport = port.lumpedport;
    auto mesh = MakeFluxRecoveryParMesh(comm, MakeFluxRecoveryCubeMesh(false));
    SpaceOperator space_op(iodata, mesh);
    PostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);
    FluxErrorEstimator<ComplexVector> estimator(space_op, iodata);
    const int size = estimator.GetElectricFluxSpace().GetTrueVSize();
    ComplexVector e(space_op.GetNDSpace().GetTrueVSize()), d(size), d_wrong(size + 1);
    e = 1.0;
    d = 2.0;
    d_wrong = 2.0;
    CHECK_THROWS_WITH(post_op.SetSurfaceResponseCorrectedField(e, &d_wrong),
                      ContainsSubstring(size_refusal));
    CHECK_NOTHROW(post_op.SetSurfaceResponseCorrectedField(e, &d));
  }
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "FluxRecovery without interior sheet faces matches pre-integration values",
                 "[postoperator][Serial]")
{
  // Where the recovery space is not broken (no interior sheet face: a mesh without sheets,
  // or the cracked device meshes of the legacy key) FluxRecovery keeps its pre-integration
  // values. The reference values are the pre-integration tree's (9269317f46, the tree of
  // palace-5c6ca28) for the same computation (P3-17; the probe commit 8534f9063f of branch
  // upsync/P3B/probe-preinteg, built as the dev build p3b-probe-9269317, np 1), within the
  // flux recovery's tolerance (EstimatorTol 1e-14; MEASURED: every value within 2e-15
  // relative, the raw entry 2 included, which takes upstream's surface functional).
  // EXCEPT the three entries that read the recovered D (1, 3, 4): upstream #1072
  // (0233a2f299, the energy-norm flux projection in FluxProjector) changed the recovered
  // flux, so they are re-pinned to the merged head's values (P4 SY1; ruling 2036 (1);
  // MEASURED: with that projector change reversed the old values return to <= 3e-16):
  //   1: 8.92520238653244898e-04 -> 5.2421642661993142e-04
  //   3: 1.59993232981066910e-03 -> 1.5998888424706354e-03
  //   4: 1.00016919572039689e-04 -> 1.0000662242343952e-04
  // The accuracy of FluxRecovery readings under the new projector is not re-established
  // by this test.
  const auto values = FluxRecoveryNoSheetEnergies(Mpi::World(), temp_dir.string());
  REQUIRE(values.size() == 5);
  const std::vector<double> reference = {1.28333333333333122e+00, 5.2421642661993142e-04,
                                         3.66666666666666833e-04, 1.5998888424706354e-03,
                                         1.0000662242343952e-04};
  for (std::size_t i = 0; i < values.size(); i++)
  {
    INFO("value " << i << " = " << std::setprecision(17) << values[i] << ", reference "
                  << reference[i]);
    CHECK_THAT(values[i], WithinRel(reference[i], 1.0e-12));
  }
  // The FluxRecovery entry reads D, not the raw field.
  CHECK(std::abs(values[1] - values[2]) > 1.0e-3 * values[2]);
}

TEST_CASE_METHOD(test::SharedTempDir, "Transient FluxRecovery sets the recovered D",
                 "[postoperator][Serial][Parallel]")
{
  // The transient driver recovers D with the estimator and sets it before each measurement
  // when FluxRecovery is requested (ruling 1574; SD-504): the measured interface energy is
  // then the one of that D (the substrate-air normal part reads D, the tangential part E),
  // and a measurement without a D for its time step is refused (P3-20).
  auto comm = Mpi::World();
  IoData iodata =
      MakeFluxRecoveryIoData(ProblemType::TRANSIENT, {1}, InterfaceDielectric::SA,
                             flux_interface_attr, true, temp_dir.string());
  config::BoundaryData port(json{{"LumpedPort",
                                  {{{"Attributes", {2}},
                                    {"Index", 1},
                                    {"R", 50.0},
                                    {"Direction", "+X"},
                                    {"Excitation", true}}}}});
  iodata.boundaries.lumpedport = port.lumpedport;
  auto mesh = MakeFluxRecoveryParMesh(comm, MakeFluxRecoveryCubeMesh(false));
  SpaceOperator space_op(iodata, mesh);
  FluxRecoveryPostOperator<ProblemType::TRANSIENT> post_op(iodata, space_op);
  REQUIRE(post_op.NeedsRecoveredElectricFlux());

  // A smooth electric field with normal and tangential parts on the interface z = 0.5.
  GridFunction E(space_op.GetNDSpace());
  mfem::VectorFunctionCoefficient f(3,
                                    [](const mfem::Vector &x, mfem::Vector &v)
                                    {
                                      v(0) = 1.0 + 0.5 * x(1);
                                      v(1) = 0.3 * x(0) * x(2);
                                      v(2) = 2.0 - x(0) + 0.2 * x(1);
                                    });
  E.Real().ProjectCoefficient(f);
  Vector e(space_op.GetNDSpace().GetTrueVSize()), b(space_op.GetRTSpace().GetTrueVSize());
  E.Real().GetTrueDofs(e);
  E.Real().ExchangeFaceNbrData();
  b = 0.0;
  FluxErrorEstimator<Vector> estimator(space_op, iodata);
  Vector d(estimator.GetElectricFluxSpace().GetTrueVSize());
  estimator.RecoverElectricFlux(e, d);

  auto Measured = [&](const Vector &flux, int step)
  {
    post_op.SetRecoveredElectricFlux(flux);
    post_op.MeasureAndPrintAll(step, e, b, 0.1 * step, 0.0);
    REQUIRE(post_op.Cache().interface_eps_i.size() == 1);
    return post_op.Cache().interface_eps_i.front().energy;
  };
  auto Expected = [&](const Vector &flux)
  {
    GridFunction D(space_op.GetRTSpace());
    D.Real().SetFromTrueDofs(flux);
    D.Real().ExchangeFaceNbrData();
    return post_op.GetSurfacePostOp().GetInterfaceElectricFieldEnergy(1, E, &D);
  };
  const double energy = Measured(d, 0);
  CHECK(energy > 0.0);
  CHECK_THAT(energy, WithinRel(Expected(d), 1.0e-12));

  // A different D changes the measured energy (the normal part reads D, not E).
  Vector d2(d);
  d2 *= 2.0;
  const double energy2 = Measured(d2, 1);
  CHECK_THAT(energy2, WithinRel(Expected(d2), 1.0e-12));
  CHECK(std::abs(energy2 - energy) > 1.0e-6 * energy);

  // Each D serves one measurement: the next time step without a new D is refused.
  CHECK_THROWS_WITH(
      post_op.MeasureAndPrintAll(2, e, b, 0.2, 0.0),
      ContainsSubstring("requires the recovered electric flux of the current time step"));
}

TEST_CASE_METHOD(test::SharedTempDir,
                 "FluxRecovery next to an interior lumped port on the broken flux space",
                 "[postoperator][Serial][Parallel]")
{
  // P3-14 (SD-3202 (c)): an interior lumped port is never cracked, so under the legacy key
  // too the drivers' recovered D lives on a flux space broken across it (refused on M1,
  // decision 1578). The postprocessing accepts it with that space's size, holds it per side
  // for the energies, and writes the corrected flux (ParaView, grid functions) as its FIRST
  // version on the RT space (upstream's rule for a broken field); a D of the RT size is
  // refused.
  auto comm = Mpi::World();
  IoData iodata = MakeFluxRecoveryIoData(ProblemType::DRIVEN, {1}, InterfaceDielectric::SA,
                                         flux_interface_attr, true, temp_dir.string());
  iodata.solver.driven.sample_f = {1.0};
  config::BoundaryData port(json{{"LumpedPort",
                                  {{{"Attributes", {flux_sheet_attr}},
                                    {"Index", 1},
                                    {"R", 50.0},
                                    {"Direction", "+Y"},
                                    {"Excitation", true}}}}});
  iodata.boundaries.lumpedport = port.lumpedport;
  auto mesh = MakeFluxRecoveryParMesh(comm, MakeFluxRecoveryCubeMesh(true));
  SpaceOperator space_op(iodata, mesh);
  FluxRecoveryPostOperator<ProblemType::DRIVEN> post_op(iodata, space_op);
  REQUIRE(post_op.NeedsRecoveredElectricFlux());
  FluxErrorEstimator<ComplexVector> estimator(space_op, iodata);
  const auto &flux_space = estimator.GetElectricFluxSpace();
  REQUIRE(flux_space.IsBroken());
  CHECK(post_op.GetRecoveredElectricFluxTrueVSize() == flux_space.GetTrueVSize());
  CHECK(flux_space.GlobalTrueVSize() > space_op.GetRTSpace().GlobalTrueVSize());

  // A smooth complex field and its recovered flux on the broken space.
  GridFunction E(space_op.GetNDSpace(), true);
  mfem::VectorFunctionCoefficient fr(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = 1.0 + 0.5 * x(1);
                                       v(1) = 0.3 * x(0) * x(2);
                                       v(2) = 2.0 - x(0) + 0.2 * x(1);
                                     });
  mfem::VectorFunctionCoefficient fi(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = -0.2 * x(2);
                                       v(1) = 0.4 + 0.1 * x(0);
                                       v(2) = 0.3 * x(1);
                                     });
  E.Real().ProjectCoefficient(fr);
  E.Imag().ProjectCoefficient(fi);
  ComplexVector e(space_op.GetNDSpace().GetTrueVSize()), d(flux_space.GetTrueVSize());
  E.Real().GetTrueDofs(e.Real());
  E.Imag().GetTrueDofs(e.Imag());
  estimator.RecoverElectricFlux(e, d);
  CHECK_NOTHROW(post_op.SetRecoveredElectricFlux(d));
  ComplexVector d_rt(space_op.GetRTSpace().GetTrueVSize());
  d_rt = 0.0;
  CHECK_THROWS_WITH(post_op.SetRecoveredElectricFlux(d_rt),
                    ContainsSubstring("does not match the flux recovery space"));

  // The corrected flux output = the FIRST version of D on the RT space.
  post_op.SetSurfaceResponseCorrectedField(e, &d);
  const auto first = flux_space.BuildVersionMap(FiniteElementSpace::VersionMap::FIRST);
  for (const bool imag : {false, true})
  {
    Vector expected(space_op.GetRTSpace().GetTrueVSize()),
        actual(space_op.GetRTSpace().GetTrueVSize());
    first->Mult(imag ? d.Imag() : d.Real(), expected);
    const auto &output = post_op.CorrectedFluxOutput();
    (imag ? output.Imag() : output.Real()).GetTrueDofs(actual);
    double max_expected = expected.Size() ? expected.Normlinf() : 0.0, max_diff = 0.0;
    for (int i = 0; i < expected.Size(); i++)
    {
      max_diff = std::max(max_diff, std::abs(actual(i) - expected(i)));
    }
    Mpi::GlobalMax(1, &max_expected, comm);
    Mpi::GlobalMax(1, &max_diff, comm);
    INFO((imag ? "imaginary" : "real")
         << " part: max. difference " << max_diff << " (max. " << max_expected << ")");
    CHECK(max_expected > 0.0);
    CHECK(max_diff <= 1.0e-14 * max_expected);
  }
}
