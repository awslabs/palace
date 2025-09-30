// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "models/postoperator.hpp"

using namespace palace;

// Helpers

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
                                   (1 + randd(9999) / 10000), 1e9 / (1 + randd(9999))});
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
}
