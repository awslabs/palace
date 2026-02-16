// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <complex>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "fem/integrator.hpp"
#include "models/postoperator.hpp"
#include "utils/iodata.hpp"
#include "utils/units.hpp"

using namespace palace;
using namespace Catch::Matchers;
using json = nlohmann::json;
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

TEST_CASE("GridFunction export", "[gridfunction][Serial][Parallel]")
{
  // Create iodata.
  Units units(0.496, 1.453);
  IoData iodata = IoData(units);

  iodata.domains.materials.emplace_back().attributes = {1};
  iodata.boundaries.pec.attributes = {1};
  iodata.problem.output_formats.gridfunction = true;
  iodata.CheckConfiguration();  // initializes quadrature

  // Setup lumped port boundary data for driven and transient.
  auto filename = fmt::format("{}/{}", PALACE_TEST_DIR, "config/boundary_configs.json");
  auto jsonstream = PreprocessFile(filename.c_str());  // Apply custom palace json
  auto config = json::parse(jsonstream);
  auto entry = config.find("boundaries_lumped_port_X_2")->find("Boundaries");
  config::BoundaryData boundary_port(*entry);

  // Create serial mesh.
  int resolution = 3;
  std::unique_ptr<mfem::Mesh> serial_mesh =
      std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
          resolution, resolution, resolution, mfem::Element::TETRAHEDRON));

  // Create parallel mesh.
  auto comm = Mpi::World();
  int size = Mpi::Size(comm);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.NondimensionalizeInputs(*par_mesh);
  std::vector<std::unique_ptr<Mesh>> mesh;
  mesh.push_back(std::make_unique<Mesh>(std::move(par_mesh)));

  auto check_files = [&](const std::string &subdir, int step, int pad_digits,
                         const std::vector<std::string> &fields)
  {
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
        B(space_op.GetRTSpace().GetTrueVSize());
    E = 0.0;
    B = 0.0;
    post_op.MeasureAndPrintAll(1, 0, E, B, 1.0);
    check_files("driven", 1, post_op.GetPadDigitsDefault(),
                {"E_real", "E_imag", "B_real", "B_imag", "S", "U_e", "U_m"});
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
    post_op.MeasureAndPrintAll(0, E, B, 1.0, 0.0, 0.0, 1);
    check_files("eigenmode", 1, post_op.GetPadDigitsDefault(),
                {"E_real", "E_imag", "B_real", "B_imag", "S", "U_e", "U_m"});
  }
}

auto CompareGridFunctions(std::vector<mfem::ParGridFunction> &gf1,
                          std::vector<mfem::ParGridFunction> &gf2,
                          const std::vector<double> &expected_ratios, double rtol)
{
  CHECK(gf1.size() == gf2.size());
  CHECK(gf1.size() == expected_ratios.size());
  for (auto i = 0; i < gf1.size(); i++)
  {
    // Get the gridfunction data.
    Vector v1 = gf1[i].GetTrueVector();
    Vector v2 = gf2[i].GetTrueVector();
    v2 /= v1;
    // Get the min/max of the ratio.
    double min = v2.Min();
    double max = v2.Max();
    // Check that min == max == expected_ratio.
    CHECK_THAT(min, Catch::Matchers::WithinRel(max, rtol));
    CHECK_THAT(min, Catch::Matchers::WithinRel(expected_ratios[i], rtol));
  }
}

TEST_CASE("Dimensional field output", "[postoperator][Serial][Parallel]")
{
  constexpr double rtol = 1.0e-10;

  // Create iodata.
  Units units(0.496, 1.453);
  IoData iodata = IoData(units);
  iodata.domains.materials.emplace_back().attributes = {1};

  // Create serial mesh.
  int resolution = 3;
  double length = 12.2;
  std::unique_ptr<mfem::Mesh> serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(resolution, resolution, resolution,
                                  mfem::Element::TETRAHEDRON, length, length, length));

  // Create parallel mesh.
  auto comm = Mpi::World();
  int size = Mpi::Size(comm);
  auto par_mesh = std::make_unique<mfem::ParMesh>(comm, *serial_mesh);
  iodata.NondimensionalizeInputs(*par_mesh);
  Mesh mesh(std::move(par_mesh));

  // Create gridfunctions.
  int order = 1;
  mfem::H1_FECollection h1_fec(order, mesh.Dimension());
  mfem::ParFiniteElementSpace h1_fespace(&mesh.Get(), &h1_fec);
  mfem::ParGridFunction Ue_gf(&h1_fespace), Um_gf(&h1_fespace);

  mfem::ND_FECollection nd_fec(order, mesh.Dimension());
  mfem::ParFiniteElementSpace nd_fespace(&mesh.Get(), &nd_fec);

  mfem::RT_FECollection rt_fec(order, mesh.Dimension());
  mfem::ParFiniteElementSpace rt_fespace(&mesh.Get(), &rt_fec);
  mfem::ParGridFunction S_gf(&rt_fespace);

  std::unique_ptr<GridFunction> E, B, V, A;
  E = std::make_unique<GridFunction>(nd_fespace, true);
  B = std::make_unique<GridFunction>(rt_fespace, true);
  A = std::make_unique<GridFunction>(nd_fespace);
  V = std::make_unique<GridFunction>(h1_fespace);

  // Create fields, initialize them to random values, and set the corresponding
  // gridfunctions.
  ComplexVector e(nd_fespace.GetTrueVSize()), b(rt_fespace.GetTrueVSize());
  Vector a(nd_fespace.GetTrueVSize()), v(h1_fespace.GetTrueVSize());
  linalg::SetRandom(comm, e);
  linalg::SetRandom(comm, b);
  linalg::SetRandom(comm, a);
  linalg::SetRandom(comm, v);

  E->Real().SetFromTrueDofs(e.Real());
  E->Imag().SetFromTrueDofs(e.Imag());
  B->Real().SetFromTrueDofs(b.Real());
  B->Imag().SetFromTrueDofs(b.Imag());
  A->Real().SetFromTrueDofs(a);
  V->Real().SetFromTrueDofs(v);

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

  // Expected scalings (including both units dimensionalization and mesh scaling Lc0 =
  // Lc/L0). See comment in ScaleGridFunctions() in postoperator.cpp for mesh scaling
  // factor.
  auto mesh_Lc0 = iodata.units.GetMeshLengthRelativeScale();
  // E [V/m], H(curl) requires * L.
  auto scaling_E = iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0) * mesh_Lc0;
  // B [Wb/m²], H(div) requires * L^2.
  auto scaling_B =
      iodata.units.Dimensionalize<Units::ValueType::FIELD_B>(1.0) * mesh_Lc0 * mesh_Lc0;
  // V [V], H1 requires no mesh scaling.
  auto scaling_V = iodata.units.Dimensionalize<Units::ValueType::VOLTAGE>(1.0);
  // A [A/m], H(curl) requires * L.
  auto scaling_A = iodata.units.Dimensionalize<Units::ValueType::CURRENT>(1.0) * mesh_Lc0;
  // U_e = 1/2 Dᴴ E, H1 requires no mesh scaling.
  auto scaling_Ue = iodata.units.Dimensionalize<Units::ValueType::FIELD_D>(1.0) *
                    iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0);
  // U_m = 1/2 Hᴴ B, H1 requires no mesh scaling.
  auto scaling_Um = iodata.units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) *
                    iodata.units.Dimensionalize<Units::ValueType::FIELD_B>(1.0);
  // S = Re{E x H⋆}, H(div) requires * L^2.
  auto scaling_S = iodata.units.Dimensionalize<Units::ValueType::FIELD_E>(1.0) *
                   iodata.units.Dimensionalize<Units::ValueType::FIELD_H>(1.0) * mesh_Lc0 *
                   mesh_Lc0;
  std::vector<double> expected_scaling = {scaling_E,  scaling_E,  scaling_B,
                                          scaling_B,  scaling_V,  scaling_A,
                                          scaling_Ue, scaling_Um, scaling_S};

  auto update_coeff_gridfunctions = [&](bool dimensional)
  {
    Ue_gf = 0.0;
    Um_gf = 0.0;
    S_gf = 0.0;
    if (dimensional)
    {
      Ue_gf.ProjectCoefficient(*Ue_dim.get());
      Um_gf.ProjectCoefficient(*Um_dim.get());
      S_gf.ProjectCoefficient(*S_dim.get());
    }
    else
    {
      Ue_gf.ProjectCoefficient(*Ue_nondim.get());
      Um_gf.ProjectCoefficient(*Um_nondim.get());
      S_gf.ProjectCoefficient(*S_nondim.get());
    }
  };

  auto copy_gridfunctions = [&]()
  {
    mfem::ParGridFunction Er_copy(E->Real()), Ei_copy(E->Imag()), Br_copy(B->Real()),
        Bi_copy(B->Real()), V_copy(V->Real()), A_copy(A->Real()), Ue_copy(Ue_gf),
        Um_copy(Um_gf), S_copy(S_gf);
    return std::vector<mfem::ParGridFunction>(
        {Er_copy, Ei_copy, Br_copy, Bi_copy, V_copy, A_copy, Ue_copy, Um_copy, S_copy});
  };

  // Copy gridfunctions before dimensionalizing.
  update_coeff_gridfunctions(false);
  std::vector<mfem::ParGridFunction> gf_before = copy_gridfunctions();

  // Dimensionalize.
  mesh::DimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(mesh_Lc0, mesh.Dimension(), E, B, V, A);
  DimensionalizeGridFunctions(iodata.units, E, B, V, A);
  update_coeff_gridfunctions(true);
  std::vector<mfem::ParGridFunction> gf_dim = copy_gridfunctions();
  CompareGridFunctions(gf_before, gf_dim, expected_scaling, rtol);

  // Nondimensionalize.
  mesh::NondimensionalizeMesh(mesh, mesh_Lc0);
  ScaleGridFunctions(1.0 / mesh_Lc0, mesh.Dimension(), E, B, V, A);
  NondimensionalizeGridFunctions(iodata.units, E, B, V, A);
  update_coeff_gridfunctions(false);
  std::vector<mfem::ParGridFunction> gf_after = copy_gridfunctions();
  CompareGridFunctions(gf_before, gf_after, std::vector<double>(gf_after.size(), 1.0),
                       rtol);
}