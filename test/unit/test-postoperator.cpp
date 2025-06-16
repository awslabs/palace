#include <complex>
#include <iterator>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "models/postoperator.hpp"

using namespace palace;

auto RandomMeasurement(int ndomain = 5)
{
  Measurement cache;
  auto randd = [](int n) { return double(rand() % n); };
  cache.domain_E_field_energy_all = 0.0;
  cache.domain_H_field_energy_all = 0.0;
  for (int i = 0; i < ndomain; i++)
  {
    cache.domain_H_field_energy_i.emplace_back(Measurement::DomainData{i, randd(100), 0});
    cache.domain_H_field_energy_all += cache.domain_H_field_energy_i.back().energy;
    cache.domain_E_field_energy_i.emplace_back(Measurement::DomainData{i, randd(100), 0});
    cache.domain_H_field_energy_all += cache.domain_H_field_energy_i.back().energy;
  }
  for (auto &H : cache.domain_H_field_energy_i)
  {
    H.participation_ratio = H.energy / cache.domain_H_field_energy_all;
  }
  for (auto &E : cache.domain_E_field_energy_i)
  {
    E.participation_ratio = E.energy / cache.domain_E_field_energy_all;
  }

  cache.freq = std::complex(randd(100), randd(100));
  cache.eigenmode_Q = 5e6 * (randd(99) + 1);
  cache.error_abs = 1e-9 * randd(100);
  cache.error_bkwd = 1e-9 * randd(100);

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
    double R = randd(100), L = randd(100), C = randd(100);
    l.P = {1 + randd(100), -50 + randd(100)};
    l.V = {randd(100), randd(100)};
    // l.I = std::conj(l.P / l.V);
    l.I_RLC[0] = l.V / impedance(cache.freq, Branch::R, R, L, C);
    l.I_RLC[1] = l.V / impedance(cache.freq, Branch::L, R, L, C);
    l.I_RLC[2] = l.V / impedance(cache.freq, Branch::C, R, L, C);
    l.I = std::accumulate(l.I_RLC.begin(), l.I_RLC.end(), std::complex(0.0, 0.0));
    // Random sub-unit magnitude and random phase
    l.S = randd(1000) / 1000 * std::exp(std::complex(randd(100) / 100, randd(100) / 100));
    l.abs_S_ij = 20.0 * std::log10(std::abs(l.S));  // TODO: eliminate these
    l.arg_S_ij = std::arg(l.S) * 180 / M_PI;
  }

  for (auto &[i, l] : cache.wave_port_vi)
  {
    l.P = {1 + randd(100), -50 + randd(100)};
    // Random sub-unit magnitude and random phase
    l.S = randd(1000) / 1000 * std::exp(std::complex(randd(100) / 100, randd(100) / 100));
    l.abs_S_ij = 20.0 * std::log10(std::abs(l.S));  // TODO: eliminate these
    l.arg_S_ij = std::arg(l.S) * 180 / M_PI;
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
        Measurement::FluxData{i, std::complex(randd(100) / 100, randd(100) / 100),
                              static_cast<SurfaceFlux>(rand() % 3)});
  }
  int ninterface = 5;
  for (int i = 0; i < ninterface; i++)
  {
    cache.interface_eps_i.emplace_back(Measurement::InterfaceData{
        i, randd(100), randd(10000) / 10000, randd(10000) / 10000, 1e9 / randd(10000)});
  }

  return cache;
}

TEST_CASE("PostOperator", "[idempotent]")
{
  auto cache = RandomMeasurement();
}