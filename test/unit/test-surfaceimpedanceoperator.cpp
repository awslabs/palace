// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <optional>
#include <unordered_set>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/units.hpp"

namespace palace
{
using namespace Catch::Matchers;

namespace
{

// Returns the scalar coefficient value assigned to bdr_attr, or std::nullopt if the
// attribute is not present on this rank.
std::optional<double> GetBdrCoeffValue(const MaterialPropertyCoefficient &fb,
                                       const Mesh &mesh, int bdr_attr)
{
  auto ceed_attrs = mesh.GetCeedBdrAttributes(bdr_attr);
  if (ceed_attrs.Size() == 0)
  {
    return std::nullopt;
  }
  int ceed_attr = ceed_attrs[0];
  const auto &attr_mat = fb.GetAttributeToMaterial();
  int mat_idx = attr_mat[ceed_attr - 1];
  if (mat_idx < 0)
  {
    return std::nullopt;
  }
  return fb.GetMaterialProperties()(0, 0, mat_idx);
}

}  // namespace

TEST_CASE("SurfaceImpedanceOperator", "[surfaceimpedanceoperator][Serial][Parallel]")
{
  auto serial_mesh = std::make_unique<mfem::Mesh>(
      mfem::Mesh::MakeCartesian3D(1, 1, 1, mfem::Element::TETRAHEDRON));
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), *serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));

  // Sanity check: ParMesh's METIS partitioning must give every rank at least one element,
  // otherwise this rank contributes nothing and the [Parallel] run is misleading.
  REQUIRE(palace_mesh.GetNE() > 0);

  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);

  Units units(1.0, 1.0);
  const double coeff = 1.0;
  const double Rs = 50.0, Ls = 1e-9, Cs = 1e-12;

  // Globally REQUIRE that boundary attributes 1 and 2 were each visible on at least one
  // rank, so under [Parallel] no SECTION silently passes with zero CHECKs everywhere.
  auto require_global_coverage = [](bool local_v1, bool local_v2)
  {
    bool flags[2] = {local_v1, local_v2};
    Mpi::GlobalOr(2, flags, Mpi::World());
    REQUIRE(flags[0]);
    REQUIRE(flags[1]);
  };

  SECTION("Per-attribute scaling, all uncracked")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> cracked = {};

    SurfaceImpedanceOperator op({imp}, cracked, units, mat_op, palace_mesh);

    MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
    op.AddStiffnessBdrCoefficients(coeff, fb);

    auto v1 = GetBdrCoeffValue(fb, palace_mesh, 1);
    auto v2 = GetBdrCoeffValue(fb, palace_mesh, 2);
    if (v1)
    {
      CHECK_THAT(*v1, WithinRel(coeff / Ls, 1e-12));
    }
    if (v2)
    {
      CHECK_THAT(*v2, WithinRel(coeff / Ls, 1e-12));
    }
    require_global_coverage(v1.has_value(), v2.has_value());
  }

  SECTION("Per-attribute scaling, all cracked")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> cracked = {1, 2};

    SurfaceImpedanceOperator op({imp}, cracked, units, mat_op, palace_mesh);

    MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
    op.AddStiffnessBdrCoefficients(coeff, fb);

    auto v1 = GetBdrCoeffValue(fb, palace_mesh, 1);
    auto v2 = GetBdrCoeffValue(fb, palace_mesh, 2);
    if (v1)
    {
      CHECK_THAT(*v1, WithinRel(coeff / (Ls * 2.0), 1e-12));
    }
    if (v2)
    {
      CHECK_THAT(*v2, WithinRel(coeff / (Ls * 2.0), 1e-12));
    }
    require_global_coverage(v1.has_value(), v2.has_value());
  }

  SECTION("Per-attribute scaling, mixed cracked/uncracked - stiffness")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> cracked = {2};

    SurfaceImpedanceOperator op({imp}, cracked, units, mat_op, palace_mesh);

    MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
    op.AddStiffnessBdrCoefficients(coeff, fb);

    auto v1 = GetBdrCoeffValue(fb, palace_mesh, 1);
    auto v2 = GetBdrCoeffValue(fb, palace_mesh, 2);
    if (v1)
    {
      CHECK_THAT(*v1, WithinRel(coeff / (Ls * 1.0), 1e-12));
    }
    if (v2)
    {
      CHECK_THAT(*v2, WithinRel(coeff / (Ls * 2.0), 1e-12));
    }
    require_global_coverage(v1.has_value(), v2.has_value());
  }

  SECTION("Per-attribute scaling, mixed cracked/uncracked - damping")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> cracked = {2};

    SurfaceImpedanceOperator op({imp}, cracked, units, mat_op, palace_mesh);

    MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
    op.AddDampingBdrCoefficients(coeff, fb);

    auto v1 = GetBdrCoeffValue(fb, palace_mesh, 1);
    auto v2 = GetBdrCoeffValue(fb, palace_mesh, 2);
    if (v1)
    {
      CHECK_THAT(*v1, WithinRel(coeff / (Rs * 1.0), 1e-12));
    }
    if (v2)
    {
      CHECK_THAT(*v2, WithinRel(coeff / (Rs * 2.0), 1e-12));
    }
    require_global_coverage(v1.has_value(), v2.has_value());
  }

  SECTION("Per-attribute scaling, mixed cracked/uncracked - mass")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> cracked = {2};

    SurfaceImpedanceOperator op({imp}, cracked, units, mat_op, palace_mesh);

    MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
    op.AddMassBdrCoefficients(coeff, fb);

    auto v1 = GetBdrCoeffValue(fb, palace_mesh, 1);
    auto v2 = GetBdrCoeffValue(fb, palace_mesh, 2);
    if (v1)
    {
      CHECK_THAT(*v1, WithinRel(coeff * Cs / 1.0, 1e-12));
    }
    if (v2)
    {
      CHECK_THAT(*v2, WithinRel(coeff * Cs / 2.0, 1e-12));
    }
    require_global_coverage(v1.has_value(), v2.has_value());
  }
}

}  // namespace palace
