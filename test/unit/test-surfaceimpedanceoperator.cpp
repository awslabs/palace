// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <optional>
#include <unordered_set>
#include <vector>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
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

// Unit cube mesh with interior boundary elements of attribute 7 on the plane x = 0.5.
mfem::Mesh MakeInteriorBoundaryMesh()
{
  auto base = mfem::Mesh::MakeCartesian3D(2, 2, 2, mfem::Element::TETRAHEDRON);
  std::vector<int> faces;
  mfem::Array<int> fv;
  for (int f = 0; f < base.GetNumFaces(); f++)
  {
    int e1, e2;
    base.GetFaceElements(f, &e1, &e2);
    base.GetFaceVertices(f, fv);
    bool on_plane = (e2 >= 0);
    for (auto v : fv)
    {
      on_plane = on_plane && std::abs(base.GetVertex(v)[0] - 0.5) < 1.0e-12;
    }
    if (on_plane)
    {
      faces.push_back(f);
    }
  }
  mfem::Mesh mesh(3, base.GetNV(), base.GetNE(),
                  base.GetNBE() + static_cast<int>(faces.size()));
  for (int v = 0; v < base.GetNV(); v++)
  {
    mesh.AddVertex(base.GetVertex(v));
  }
  for (int e = 0; e < base.GetNE(); e++)
  {
    mesh.AddElement(base.GetElement(e)->Duplicate(&mesh));
  }
  for (int be = 0; be < base.GetNBE(); be++)
  {
    mesh.AddBdrElement(base.GetBdrElement(be)->Duplicate(&mesh));
  }
  for (auto f : faces)
  {
    auto *el = base.GetFace(f)->Duplicate(&mesh);
    el->SetAttribute(7);
    mesh.AddBdrElement(el);
  }
  mesh.FinalizeTopology();
  mesh.Finalize();
  mesh.SetAttributes();
  return mesh;
}

}  // namespace

TEST_CASE("SurfaceConductivityOperator interior sheets",
          "[surfaceimpedanceoperator][Serial][Parallel]")
{
  // An interior conductivity boundary (not cracked) is a conducting sheet with two
  // conductor surfaces, so its admittance is twice that of the exterior boundary. The
  // "External" flag only affects the thickness correction, of the exterior attributes (it
  // does not apply to interior sheets, where it is ignored with a warning).
  auto serial_mesh = MakeInteriorBoundaryMesh();
  auto par_mesh = std::make_unique<mfem::ParMesh>(Mpi::World(), serial_mesh);
  Mesh palace_mesh(std::move(par_mesh));
  REQUIRE(palace_mesh.GetNE() > 0);
  config::MaterialData material;
  material.attributes = {1};
  config::PeriodicBoundaryData periodic;
  MaterialOperator mat_op({material}, periodic, ProblemType::DRIVEN, palace_mesh);
  Units units(1.0, 1.0);
  const bool external = GENERATE(false, true);
  const double h = GENERATE(0.0, 0.1);

  // Coefficients of the extra system matrix (real and imaginary parts) and of the boundary
  // mass summed over the groups, for the exterior attribute 1 and the interior attribute 7.
  auto Values = [&](const config::ConductivityData &cond)
  {
    SurfaceConductivityOperator op({cond}, ProblemType::DRIVEN, units, mat_op,
                                   palace_mesh.Get());
    MaterialPropertyCoefficient fbr(mat_op.MaxCeedBdrAttribute()),
        fbi(mat_op.MaxCeedBdrAttribute()), fb(mat_op.MaxCeedBdrAttribute());
    op.AddExtraSystemBdrCoefficients(2.0, fbr, fbi);
    for (std::size_t g = 0; g < op.Size(); g++)
    {
      op.AddBoundaryMassBdrCoefficients(g, fb);
    }
    std::array<double, 6> vals = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    for (int k = 0; k < 3; k++)
    {
      const auto &f = (k == 0) ? fbr : ((k == 1) ? fbi : fb);
      auto v1 = GetBdrCoeffValue(f, palace_mesh, 1);
      auto v7 = GetBdrCoeffValue(f, palace_mesh, 7);
      vals[2 * k] = v1 ? *v1 : 0.0;
      vals[2 * k + 1] = v7 ? *v7 : 0.0;
    }
    // Attributes are not present on all processes (the values are positive).
    Mpi::GlobalMax(6, vals.data(), Mpi::World());
    return std::make_pair(vals, op.Size());
  };

  config::ConductivityData cond;
  cond.sigma = 1.0e3;
  cond.h = h;
  cond.external = external;
  cond.attributes = {1, 7};
  const auto [vals, num_groups] = Values(cond);

  // References: the exterior attribute alone, without "External", with the thickness h
  // and 2h.
  config::ConductivityData ref = cond;
  ref.external = false;
  ref.attributes = {1};
  const auto ref_h = Values(ref).first;
  ref.h = 2.0 * h;
  const auto ref_2h = Values(ref).first;
  const auto &ref_ext = external ? ref_2h : ref_h;

  CAPTURE(external, h);
  REQUIRE(ref_h[0] > 0.0);
  REQUIRE(ref_h[2] > 0.0);
  CHECK(num_groups == ((external && h > 0.0) ? 2 : 1));
  CHECK_THAT(vals[0], WithinRel(ref_ext[0], 1e-12));
  CHECK_THAT(vals[2], WithinRel(ref_ext[2], 1e-12));
  CHECK_THAT(vals[1], WithinRel(2.0 * ref_h[0], 1e-12));
  CHECK_THAT(vals[3], WithinRel(2.0 * ref_h[2], 1e-12));
  CHECK(vals[4] == 1.0);
  CHECK(vals[5] == 2.0);
  if (external && h > 0.0)
  {
    // The thickness correction makes a difference.
    CHECK(std::abs(ref_2h[0] - ref_h[0]) > 1.0e-6 * ref_h[0]);
  }
}

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

  SECTION("Per-attribute scaling, no split impedance")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> split = {};

    SurfaceImpedanceOperator op({imp}, split, units, mat_op, palace_mesh);

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

  SECTION("Per-attribute scaling, all split (legacy cracking)")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> split = {1, 2};

    SurfaceImpedanceOperator op({imp}, split, units, mat_op, palace_mesh);

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

  SECTION("Per-attribute scaling, mixed split - stiffness")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> split = {2};

    SurfaceImpedanceOperator op({imp}, split, units, mat_op, palace_mesh);

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

  SECTION("Per-attribute scaling, mixed split - damping")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> split = {2};

    SurfaceImpedanceOperator op({imp}, split, units, mat_op, palace_mesh);

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

  SECTION("Per-attribute scaling, mixed split - mass")
  {
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls;
    imp.Cs = Cs;
    imp.attributes = {1, 2};
    std::unordered_set<int> split = {2};

    SurfaceImpedanceOperator op({imp}, split, units, mat_op, palace_mesh);

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

  SECTION("Combined terms are added once per attribute")
  {
    // Both attributes get the same stiffness term first, so they share one entry of the
    // coefficient; the per-attribute mass term must then be added once to each of them
    // (regression: it was added once per attribute to the shared entry). Values of similar
    // magnitude so that the sum is resolved.
    const double Ls_c = 0.5, Cs_c = 0.25;
    config::ImpedanceData imp;
    imp.Rs = Rs;
    imp.Ls = Ls_c;
    imp.Cs = Cs_c;
    imp.attributes = {1, 2};
    std::unordered_set<int> split = {};

    SurfaceImpedanceOperator op({imp}, split, units, mat_op, palace_mesh);

    MaterialPropertyCoefficient fb(mat_op.MaxCeedBdrAttribute());
    op.AddStiffnessBdrCoefficients(coeff, fb);
    op.AddMassBdrCoefficients(coeff, fb);

    auto v1 = GetBdrCoeffValue(fb, palace_mesh, 1);
    auto v2 = GetBdrCoeffValue(fb, palace_mesh, 2);
    if (v1)
    {
      CHECK_THAT(*v1, WithinRel(coeff / Ls_c + coeff * Cs_c, 1e-12));
    }
    if (v2)
    {
      CHECK_THAT(*v2, WithinRel(coeff / Ls_c + coeff * Cs_c, 1e-12));
    }
    require_global_coverage(v1.has_value(), v2.has_value());
  }
}

}  // namespace palace
