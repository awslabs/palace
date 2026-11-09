// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include <array>
#include <cmath>
#include <complex>
#include <numbers>
#include <Eigen/Dense>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "drivers/drivensolver.hpp"
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/mesh.hpp"
#include "models/materialoperator.hpp"
#include "models/portexcitations.hpp"
#include "models/postoperator.hpp"
#include "models/postoperatorcsv.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/filesystem.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

using namespace palace;
using namespace nlohmann;
using namespace Catch::Matchers;

class LumpedPortDataTest : public LumpedPortData
{
public:
  void InitializeLinearFormsForTest(mfem::ParFiniteElementSpace &nd_fespace) const
  {
    InitializeLinearForms(nd_fespace);
  }
  mfem::LinearForm *GetLinearFormS() const { return s.get(); }
  mfem::LinearForm *GetLinearFormV() const { return v.get(); }
};

auto LoadScaleParMesh(IoData &iodata, MPI_Comm world_comm)
{
  // Load Mesh — copy from main.cpp
  std::vector<std::unique_ptr<Mesh>> mesh_;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    auto smesh = mesh::Load(iodata, world_comm);
    int ne = smesh ? smesh->GetNE() : 0;
    Mpi::GlobalMax(1, &ne, world_comm);
    if (Mpi::Size(world_comm) > ne)
    {
      SKIP("Mesh has fewer elements than processes");
    }
    if (iodata.model.Lc <= 0.0)
    {
      iodata.model.Lc = mesh::ComputeReferenceLength(smesh, world_comm);
    }
    iodata.NondimensionalizeInputs(smesh);
    mfem_mesh.push_back(mesh::Partition(iodata, std::move(smesh), world_comm));
    mesh::RefineMesh(iodata, mfem_mesh);
    for (auto &m : mfem_mesh)
    {
      mesh_.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }
  return mesh_;
}

TEST_CASE("LumpedPort batched voltages", "[lumped_port][Serial][Parallel][GPU]")
{
  MPI_Comm comm = Mpi::World();
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_3_2_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_3_2_1_tet.msh"));
  auto order = GENERATE(1, 2);
  CAPTURE(mesh_is_hex, order);

  IoData iodata{Units(1.0e-6, 7.0)};
  iodata.model.mesh = mesh_path;
  iodata.model.L0 = 1.0e-6;
  iodata.model.Lc = 7.0;
  iodata.model.crack_bdr_elements = false;
  iodata.solver.order = order;
  json domains_json = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  iodata.domains = config::DomainData(domains_json);
  json boundaries_json = {
      {"LumpedPort",
       json::array(
           {json::object({{"Index", 1},
                          {"R", 50.0},
                          {"Excitation", uint(1)},
                          {"Attributes", json::array({1})},
                          {"Direction", "+Y"}}),
            json::object(
                {{"Index", 2},
                 {"R", 75.0},
                 {"Elements",
                  json::array({json::object({{"Attributes", json::array({9, 11})},
                                             {"Direction", "+Y"}}),
                               json::object({{"Attributes", json::array({2, 6, 10})},
                                             {"Direction", "+X"}})})}})})}};
  iodata.boundaries = config::BoundaryData(boundaries_json);
  iodata.CheckConfiguration();

  auto mesh = LoadScaleParMesh(iodata, comm);
  SpaceOperator space_op(iodata, mesh);
  GridFunction E(space_op.GetNDSpace(), true);
  mfem::VectorFunctionCoefficient Er(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = 1.0 + x(1);
                                       v(1) = 0.5 + 2.0 * x(0);
                                       v(2) = x(0) - x(1);
                                     });
  mfem::VectorFunctionCoefficient Ei(3,
                                     [](const mfem::Vector &x, mfem::Vector &v)
                                     {
                                       v(0) = x(0) + x(2);
                                       v(1) = 1.0 - x(1);
                                       v(2) = 0.25 + x(0) * x(1);
                                     });
  E.Real().ProjectCoefficient(Er);
  E.Imag().ProjectCoefficient(Ei);

  auto CheckBatched = [&]()
  {
    std::map<int, std::complex<double>> expected;
    for (const auto &[idx, data] : space_op.GetLumpedPortOp())
    {
      expected.emplace(idx, data.GetVoltage(E));
    }
    const auto values = space_op.GetLumpedPortOp().GetVoltages(E);
    REQUIRE(values.size() == expected.size());
    for (const auto &[idx, ref] : expected)
    {
      CAPTURE(idx, ref.real(), ref.imag(), values.at(idx).real(), values.at(idx).imag());
      const double tol = 1.0e-11 * std::max(std::abs(ref), 1.0);
      CHECK_THAT(values.at(idx).real(), WithinAbs(ref.real(), tol));
      CHECK_THAT(values.at(idx).imag(), WithinAbs(ref.imag(), tol));
    }
  };

  CheckBatched();
  E *= 1.75;
  CheckBatched();
}

TEST_CASE("LumpedPort_BasicTests_1ElementPort_Cube321", "[lumped_port][Serial][Parallel]")
{
  // This is a test of a square lumped port with R=50 Ohm. It tests the geometry,
  // integrators, and port vectors associated with lumped port properties.
  //
  // The domain is a rectangle of dimension (dx,dy,dz) = (3,2,1). We are going to paint the
  // lumped port on the boundary square with diagonal (0,0,0) to (1,1,0). The port direction
  // is "+Y" (0,1,0). The port is made from a single attribute label in the mesh: for a hex
  // mesh this label is a single square, for the tet mesh this is two triangles. The bulk
  // domain is a trivial material (ε=1.0, µ=1.0, no loss).
  //
  // This test iterates: (a) without any PEC boundary, (b) with a metal (PEC) boundary
  // neighbor next to the port (SIDE), and (c) with a metal boundary neighbor at the end of
  // the port (END). In (b), we set the boundary square with diagonal (1,0,0)-(2,1,0) to be
  // PEC. This will cause the electric field of port on the shared edged to be set to zero
  // and change integrals and normalizations as the effective width of the port is changed.
  // In (c), we set the boundary square with diagonal (0, 1, 0)-(1, 2, 0) to be PEC. This
  // should not alter the fields which are perpendicularly incident on the metal so should
  // give the same result as (a). However, this is not true for tet meshes
  //
  // All other boundaries are left unspecified (so PMC).
  //
  // This tests also iterates over solver orders 1-3.

  using VT = palace::Units::ValueType;
  MPI_Comm world_comm = Mpi::World();

  // Generate 3 x 2 x 3 different test configuration.
  auto solver_order = GENERATE(1, 2, 3);
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_3_2_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_3_2_1_tet.msh"));

  // Cases: no PEC, one where element neighboring port
  enum class PEC : std::uint8_t
  {
    NONE,
    SIDE,
    ENDS
  };
  const auto &[boundary_pec_type, boundary_pec_attr] =
      GENERATE(std::make_tuple(PEC::NONE, json::array({})),
               std::make_tuple(PEC::SIDE, json::array({5})),
               std::make_tuple(PEC::ENDS, json::array({3})));

  // Construct IoData from scratch.
  double L0 = 1.0e-6;
  double Lc = 7.0;

  IoData iodata{Units(L0, Lc)};
  iodata.model.mesh = mesh_path;
  iodata.model.L0 = L0;
  iodata.model.Lc = Lc;
  iodata.model.crack_bdr_elements = false;

  iodata.solver.order = solver_order;

  json domains_json = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  iodata.domains = config::DomainData(domains_json);

  // Case 1: Put in a single port with single attribute.
  json boundary_json_1 = {
      {"PEC", json::object({{"Attributes", boundary_pec_attr}})},
      {"LumpedPort", json::array({json::object({{"Index", 1},
                                                {"R", 50.0},
                                                {"Excitation", uint(1)},
                                                {"Attributes", json::array({1})},
                                                {"Direction", "+Y"}})})}};
  iodata.boundaries = config::BoundaryData(boundary_json_1);
  iodata.CheckConfiguration();

  auto mesh_io = LoadScaleParMesh(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  const auto &mesh = space_op.GetNDSpace().GetParMesh();
  CHECK(mesh.SpaceDimension() == 3);

  const auto &port_1 = space_op.GetLumpedPortOp().GetPort(1);

  // LumpedPortData Basics.
  CHECK(port_1.HasExcitation());
  CHECK(port_1.elems.size() == 1);

  // Power normalization and corresponding excitation voltage.
  // These analytic values are actually not valid for PEC::SIDE see below.
  CHECK_THAT(port_1.GetExcitationPower(),
             WithinRel(iodata.units.Nondimensionalize<VT::POWER>(1.0)));
  CHECK_THAT(port_1.GetExcitationVoltage(),
             WithinRel(iodata.units.Nondimensionalize<VT::VOLTAGE>(std::sqrt(50.))));

  // Properties of single rectangular element in port.
  const UniformElementData *el_ptr =
      dynamic_cast<const UniformElementData *>(port_1.elems.begin()->get());
  REQUIRE(el_ptr != nullptr);

  double length_dx_m = iodata.units.Nondimensionalize<VT::LENGTH>(1.0 * L0);  // in [m]

  CHECK_THAT(el_ptr->GetGeometryLength(), WithinRel(length_dx_m));
  CHECK_THAT(el_ptr->GetGeometryWidth(), WithinRel(length_dx_m));

  // Scale factor between element wave impedance and port impedance: (W / L) * n_elems
  CHECK_THAT(port_1.GetToSquare(*el_ptr), WithinRel(1.0));

  // Validate against mesh properties at the attribute level.
  // TODO: Factor this out as a separate test for geometries of boundaries.
  // TODO: Need test for non-axis aligned port.
  {
    CHECK(el_ptr->GetAttrList().Size() == 1);
    auto el_attr = el_ptr->GetAttrList()[0];
    auto bbox = mesh::GetBoundingBox(mesh, el_attr, true);

    CHECK(bbox.planar);
    CHECK_THAT(bbox.Area(), WithinRel(length_dx_m * length_dx_m));
    CHECK(bbox.Volume() == 0.0);

    // Can do exact equals below as axis alignment of port means that there should be
    // strictly no double rounding.

    auto bbox_lengths = bbox.Lengths();
    CHECK(bbox_lengths.Size() == 3);
    CHECK(bbox_lengths(0) == length_dx_m);
    CHECK(bbox_lengths(1) == length_dx_m);
    CHECK(bbox_lengths(2) == 0.0);
  }

  mfem::Array<int> attr_marker_loc;
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mesh::AttrToMarker(bdr_attr_max, el_ptr->GetAttrList(), attr_marker_loc);

  // Geometry Tests: not via point-cloud construction as for bbox above
  CHECK_THAT(mesh::GetSurfaceArea(mesh, attr_marker_loc),
             WithinRel(length_dx_m * length_dx_m));

  // Orientation of normal hard to extract here, check magnitude
  auto normal = mesh::GetSurfaceNormal(mesh, attr_marker_loc);
  CHECK(std::abs(normal[0]) == 0.0);
  CHECK(std::abs(normal[1]) == 0.0);
  CHECK(std::abs(normal[2]) == 1.0);

  // Make Port Bilinear forms:

  // - Identity coefficient restricted to the boundary which is ~ \delta_{i,j}
  mfem::ParBilinearForm g_boundary_id(&space_op.GetNDSpace().Get());
  SumCoefficient fb_id{};

  // - Mode coefficient on boundary, here ~ \hat{y}_i \hat{y}_j.
  mfem::ParBilinearForm g_boundary_mode(&space_op.GetNDSpace().Get());
  SumVectorCoefficient fb_mode(mesh.SpaceDimension());

  // - Filter coefficient with in-plane vector polarization perpendicular to port, here ~
  //   \hat{x}_i \hat{x}_j. Only zeros out if mode is perfectly along \hat{y}, which is
  //   *not* the case for neighbouring PEC boundary conditions.
  mfem::ParBilinearForm g_boundary_perp(&space_op.GetNDSpace().Get());
  SumVectorCoefficient fb_perp(mesh.SpaceDimension());

  fb_id.AddCoefficient(std::make_unique<RestrictedCoefficient<mfem::ConstantCoefficient>>(
      el_ptr->GetAttrList(), 1.0));
  g_boundary_id.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_id));
  g_boundary_id.Assemble(true);

  fb_mode.AddCoefficient(el_ptr->GetModeCoefficient(1.0));
  g_boundary_mode.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_mode));
  g_boundary_mode.Assemble(true);

  // Port mode is in y direction, so filter in x.
  StaticVector<3> x_dir;
  x_dir.UseDevice(true);
  x_dir = 0.0;
  x_dir(0) = 1.0;
  fb_perp.AddCoefficient(
      std::make_unique<RestrictedVectorCoefficient<mfem::VectorConstantCoefficient>>(
          el_ptr->GetAttrList(), x_dir));
  g_boundary_perp.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_perp));
  g_boundary_perp.Assemble(true);

  // Assembly primary vector h_t x n = e_t / eta of port; since TEM this is purely
  // tangential. Field normalization in GetLumpedPortExcitationVectorPrimaryHt is such that
  // reference impedance is 1.0 in internal units, so that it is Z_0 in physical units.
  ComplexVector port_primary_ht_cn;
  port_primary_ht_cn.UseDevice(true);
  space_op.GetLumpedPortExcitationVectorPrimaryHtcn(1, port_primary_ht_cn);
  // As palace GridFunction: Undo restriction
  GridFunction port_primary_gf_ht_cn(space_op.GetNDSpace());
  port_primary_gf_ht_cn = 0.0;
  port_primary_gf_ht_cn.Real().SetFromTrueDofs(port_primary_ht_cn.Real());

  ComplexVector port_primary_et;
  port_primary_et.UseDevice(true);
  space_op.GetLumpedPortExcitationVectorPrimaryEt(1, port_primary_et);
  GridFunction port_primary_gf_et(space_op.GetNDSpace());
  port_primary_gf_et = 0.0;
  port_primary_gf_et.Real().SetFromTrueDofs(port_primary_et.Real());

  // Total power for the (real) fields is (e_t / eta) \cdot e_t = 1.
  //
  // Note: If PEC::SIDE then the power normalization is wrong, since we have just removed
  // some of the electric field without properly compensating in the integral. This defines
  // a new alpha normalization of the field that we want to divide out to return to the
  // original p_0 = 1.
  double port_power = g_boundary_id.ParInnerProduct(port_primary_gf_ht_cn.Real(),
                                                    port_primary_gf_et.Real());
  double port_normalization_alpha = std::sqrt(port_power);  // Marks-Williams alpha
  if (boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(port_power, WithinRel(1.0));
  }

  // Plain normalization of et
  double norm_id_et =
      g_boundary_id.ParInnerProduct(port_primary_gf_et.Real(), port_primary_gf_et.Real());

  // The normalization of e_t is Z_R \sum_e  W_e / L_e but with Z_R = 1.0 in the
  // internal units.
  double et_norm_expected = port_1.GetExcitationFieldEtNormSqWithUnityZR();

  if (boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(norm_id_et, WithinRel(et_norm_expected, 1e-13));
  }
  else
  {
    // Need to add power correction factor, since the power != 1.
    // Can only fix this here since ther is a single element, so we can compensate error.
    CHECK_THAT(norm_id_et / (port_normalization_alpha * port_normalization_alpha),
               WithinRel(et_norm_expected));
  }

  // Plain normalization of h_t x n = e_t / eta
  double norm_id_ht_cn = g_boundary_id.ParInnerProduct(port_primary_gf_ht_cn.Real(),
                                                       port_primary_gf_ht_cn.Real());

  // The normalization of e_t / eta is 1 / (Z_R n_el^2) \sum_e L_e / W_e but with Z_R = 1.0
  // in the internal units, since we took out the 1 / sqrt(\vert Z_0 \vert) factor out of
  // e_t / eta.
  double ht_cn_norm_expected = port_1.GetExcitationFieldHtNormSqWithUnityZR();

  if (boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(norm_id_ht_cn, WithinRel(ht_cn_norm_expected));
  }
  else
  {
    // Need to add power correction factor, since the power != 1.
    // Can only fix this here since ther is a single element, so we can compensate error.
    CHECK_THAT(norm_id_ht_cn / (port_normalization_alpha * port_normalization_alpha),
               WithinRel(ht_cn_norm_expected));
  }

  double norm_mode_ht_cn = g_boundary_mode.ParInnerProduct(port_primary_gf_ht_cn.Real(),
                                                           port_primary_gf_ht_cn.Real());
  double norm_perp_ht_cn = g_boundary_perp.ParInnerProduct(port_primary_gf_ht_cn.Real(),
                                                           port_primary_gf_ht_cn.Real());

  // Validate port integration normalisation. This is always true for the single hex mesh.
  // or if there is no proximate metal. Tet mesh mixes directions.
  if (mesh_is_hex || boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(norm_id_ht_cn, WithinRel(norm_mode_ht_cn));
    CHECK_THAT(norm_perp_ht_cn, WithinAbs(0.0, 1e-12));
  }
  else
  {
    // If PEC side, then e_t is a non-uniform field pattern that depends on the details of
    // Nédélec tet elements. Here only check trivial statement g_id = g_mode + g_perp.
    CHECK_THAT(norm_id_ht_cn, WithinRel(norm_mode_ht_cn + norm_perp_ht_cn));
  }

  // Measure voltage.
  //
  // Note: On paper, we have normalized $\vert v_0 \vert^2 = \vert Z_R \vert$ on ports, with
  // the assumption of unit power.
  //
  // For the case with side metal, we have \alpha != 1. We need to divide out a factor alpha
  // to normalize e_t to the original p_0 = 1 field conventions.
  //
  // However, this is not sufficient! Because the e_t vector is now no longer equal the to
  // vector field in the LinearForm used in GetVoltage (which does not have PEC zeroed).
  // There are two options here. (1) Either redefine the meaning of voltage as an integral
  // with a e_t field with metal zeroed in, which is strange but internally consistent. (2)
  // Keep the voltage definition as is but this introduces another factor, that effectively
  // is a Marks-Williams \beta factor.
  //
  // Need to correct for different mode. However, this is not just the port power, since the
  // voltage integral along \hat{y} is no longer in the same direction everywhere as $e_t$.
  auto v_in_measure = port_1.GetVoltage(port_primary_gf_et);
  CHECK(v_in_measure.imag() == 0.0);
  double port_normalization_beta =
      iodata.units.Nondimensionalize<palace::Units::ValueType::VOLTAGE>(
          std::sqrt(electromagnetics::Z0_)) /
      (v_in_measure.real());
  if (boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(port_normalization_beta, WithinRel(1.0));
  }

  // Now "s" form in LumpedPortData is just e_t / eta (with Z_R = R). So check inner
  // product.
  //
  // Note: s should not change during change of field normalization alpha, but formula
  // assumes 1.0 and does not normalized. Also, the integral seems wrong with metal, as the
  // SParameter Linear form should be the field form e_t / eta.
  if (boundary_pec_type != PEC::SIDE)
  {
    std::complex<double> s_param = port_1.GetSParameter(port_primary_gf_ht_cn);
    CHECK_THAT(s_param.real(),
               WithinRel(ht_cn_norm_expected /
                         iodata.units.Nondimensionalize<VT::VOLTAGE>(std::sqrt(50.))));
    CHECK_THAT(s_param.imag(), WithinAbs(0.0, 1e-12));
  }

  // Extract s and v linear forms to check equality. Use reinterpret cast to identical
  // layout Test structure from original structure. Done not to populate code with friend
  // functions, but still access private members for testing.
  const auto &port_1_test_cast = reinterpret_cast<const LumpedPortDataTest &>(port_1);
  port_1_test_cast.InitializeLinearFormsForTest(space_op.GetNDSpace().Get());
  auto *form_s = port_1_test_cast.GetLinearFormS();
  ComplexVector VecFormS;
  VecFormS.SetSize(space_op.GetNDSpace().GetTrueVSize());
  VecFormS.UseDevice(true);
  space_op.GetNDSpace().GetProlongationMatrix()->MultTranspose(*form_s, VecFormS.Real());

  auto *form_v = port_1_test_cast.GetLinearFormV();
  ComplexVector VecFormV;
  VecFormV.SetSize(space_op.GetNDSpace().GetTrueVSize());
  VecFormV.UseDevice(true);
  space_op.GetNDSpace().GetProlongationMatrix()->MultTranspose(*form_v, VecFormV.Real());

  // Now GetExcitationVector1 in space op, returns the dual of 2 e_t / eta (with Z_R = R).
  ComplexVector RHS;
  RHS.UseDevice(true);
  space_op.GetExcitationVector1(1, RHS);
  RHS *= 0.5;

  // In the case where there is a neighbouring PEC condition, RHS and VecFormS are no
  // longer equal — the DoF on the shared edge with the PEC neighbours are zeroed.
  const auto &nd_dbc_tdof = space_op.GetNDDbcTDofLists().back();
  if (boundary_pec_type != PEC::NONE)
  {
    // If more than 1 MPI rank only know that total tdof is non-zero.
    auto nd_dbc_tdof_size = nd_dbc_tdof.Size();
    Mpi::GlobalSum(1, &nd_dbc_tdof_size, world_comm);
    CHECK(nd_dbc_tdof_size > 0);
  }

  REQUIRE(RHS.Real().Size() == VecFormS.Real().Size());
  // This should true up to normalization, but is not currently due to non-subtracted metal.
  for (int i = 0; i < RHS.Real().Size(); i++)
  {
    // If subtracted metal, check RHS is zero and continue. Tet mesh also subtracts active
    // dofs if PEC is at the end of the port.
    if (boundary_pec_type == PEC::SIDE || (!mesh_is_hex && boundary_pec_type == PEC::ENDS))
    {
      const auto *find_pec = std::ranges::find(nd_dbc_tdof, i);
      if (find_pec != nd_dbc_tdof.end())
      {
        CHECK(RHS.Real()[i] == 0.0);
        continue;
      }
    }
    CHECK_THAT(RHS.Real()[i], WithinAbs(VecFormS.Real()[i], 1e-12));
  }

  // Check that form_s and form_v are the same apart from the factor 1 / sqrt(R). In actual
  // fact, this should only true in the ideal case, without any metal subtraction. Here it
  // is always true because of the form definitions.
  REQUIRE(VecFormS.Real().Size() == VecFormV.Real().Size());
  for (int i = 0; i < VecFormS.Real().Size(); i++)
  {
    CHECK_THAT(VecFormV.Real()[i],
               WithinAbs(VecFormS.Real()[i] *
                             iodata.units.Nondimensionalize<VT::VOLTAGE>(std::sqrt(50.)),
                         1e-12));
  }
}

TEST_CASE("LumpedPort_BasicTests_3ElementPort_Cube321", "[lumped_port][Serial][Parallel]")
{
  // Similar to LumpedPort_BasicTests_1ElementPort_Cube321 test above,
  // but now the single lump port consists of three disjoint elements to test

  using VT = palace::Units::ValueType;
  MPI_Comm world_comm = Mpi::World();

  // Generate 3 x 2 x 3 different test configuration.
  auto solver_order = GENERATE(1, 2, 3);
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_3_2_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_3_2_1_tet.msh"));

  // Cases: no PEC, one where element neighboring port
  enum class PEC : std::uint8_t
  {
    NONE,
    SIDE,
    ENDS
  };
  const auto &[boundary_pec_type, boundary_pec_attr] =
      GENERATE(std::make_tuple(PEC::NONE, json::array({})),
               std::make_tuple(PEC::SIDE, json::array({4, 5, 8})),
               std::make_tuple(PEC::ENDS, json::array({3, 18})));

  // Construct IoData from scratch.
  double L0 = 1.0e-6;
  double Lc = 7.0;

  IoData iodata{Units(L0, Lc)};
  iodata.model.mesh = mesh_path;
  iodata.model.L0 = L0;
  iodata.model.Lc = Lc;
  iodata.model.crack_bdr_elements = false;

  iodata.solver.order = solver_order;

  json domains_json = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  iodata.domains = config::DomainData(domains_json);

  json boundary_json_1 = {
      {"PEC", json::object({{"Attributes", boundary_pec_attr}})},
      {"LumpedPort",
       json::array({json::object(
           {{"Index", 1},
            {"R", 50.0},
            {"Excitation", uint(1)},
            {"Elements",
             json::array(
                 {json::object({{"Attributes", json::array({1})},  // dx,dy=1,1
                                {"Direction", "+Y"}}),
                  json::object({{"Attributes", json::array({9, 11})},  // dx,dy=1,2
                                {"Direction", "+Y"}}),
                  json::object({{"Attributes", json::array({2, 6, 10})},  // dx,dy = 3,1
                                {"Direction", "+X"}})})}})})}};
  iodata.boundaries = config::BoundaryData(boundary_json_1);
  iodata.CheckConfiguration();

  auto mesh_io = LoadScaleParMesh(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);

  const auto &mesh = space_op.GetNDSpace().GetParMesh();
  CHECK(mesh.SpaceDimension() == 3);

  const auto &port_1 = space_op.GetLumpedPortOp().GetPort(1);

  // LumpedPortData Basics.
  CHECK(port_1.HasExcitation());
  CHECK(port_1.elems.size() == 3);

  // Power normalization and corresponding excitation voltage.
  // These analytic values are actually not valid for PEC::SIDE see below.
  CHECK_THAT(port_1.GetExcitationPower(),
             WithinRel(iodata.units.Nondimensionalize<VT::POWER>(1.0)));
  CHECK_THAT(port_1.GetExcitationVoltage(),
             WithinRel(iodata.units.Nondimensionalize<VT::VOLTAGE>(std::sqrt(50.))));

  double length_dx_m = iodata.units.Nondimensionalize<VT::LENGTH>(1.0 * L0);  // in [m]

  std::vector<std::array<double, 2>> expected_length_width = {{1., 1.}, {2., 1.}, {3., 1.}};
  // Properties of single rectangular element in port.
  for (std::size_t i = 0; i < port_1.elems.size(); i++)
  {
    const auto [dl, dw] = expected_length_width.at(i);
    const auto *el_ptr = port_1.elems.at(i).get();
    REQUIRE(el_ptr != nullptr);

    CHECK_THAT(el_ptr->GetGeometryLength(), WithinRel(dl * length_dx_m));
    CHECK_THAT(el_ptr->GetGeometryWidth(), WithinRel(dw * length_dx_m));

    // Scale factor between element wave impedance and port impedance: (W / L) * n_elems
    CHECK_THAT(port_1.GetToSquare(*el_ptr), WithinRel(3 * dw / dl));
  }

  // Make Port Bilinear forms:

  // - Identity coefficient restricted to the boundary which is ~ \delta_{i,j}
  mfem::ParBilinearForm g_boundary_id(&space_op.GetNDSpace().Get());
  SumCoefficient fb_id{};

  // - Mode coefficient on boundary, here ~ \hat{y}_i \hat{y}_j.
  mfem::ParBilinearForm g_boundary_mode(&space_op.GetNDSpace().Get());
  SumVectorCoefficient fb_mode(mesh.SpaceDimension());

  // - Filter coefficient with in-plane vector polarization perpendicular to port, here ~
  //   \hat{x}_i \hat{x}_j. Only zeros out if mode is perfectly along \hat{y}, which is
  //   *not* the case for neighbouring PEC boundary conditions.
  mfem::ParBilinearForm g_boundary_perp(&space_op.GetNDSpace().Get());
  SumVectorCoefficient fb_perp(mesh.SpaceDimension());

  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;

  std::vector<int> perp_vec_dir_index = {0, 0, 1};

  for (std::size_t i = 0; i < port_1.elems.size(); i++)
  {
    const auto *el_ptr = port_1.elems.at(i).get();
    REQUIRE(el_ptr != nullptr);

    mfem::Array<int> attr_marker_loc;
    mesh::AttrToMarker(bdr_attr_max, el_ptr->GetAttrList(), attr_marker_loc);

    fb_id.AddCoefficient(std::make_unique<RestrictedCoefficient<mfem::ConstantCoefficient>>(
        el_ptr->GetAttrList(), 1.0));

    fb_mode.AddCoefficient(el_ptr->GetModeCoefficient(1.0));

    // Port mode is in y direction, so filter in x.
    StaticVector<3> perp_dir;
    perp_dir = 0.0;
    perp_dir(perp_vec_dir_index.at(i)) = 1.0;
    fb_perp.AddCoefficient(
        std::make_unique<RestrictedVectorCoefficient<mfem::VectorConstantCoefficient>>(
            el_ptr->GetAttrList(), perp_dir));
  }

  g_boundary_id.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_id));
  g_boundary_id.Assemble(true);

  g_boundary_mode.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_mode));
  g_boundary_mode.Assemble(true);

  g_boundary_perp.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_perp));
  g_boundary_perp.Assemble(true);

  // Assembly primary vector h_t x n = e_t / eta of port; since TEM this is purely
  // tangential. Field normalization in GetLumpedPortExcitationVectorPrimaryHt is such that
  // reference impedance is 1.0 in internal units, so that it is Z_0 in physical units.
  ComplexVector port_primary_ht_cn;
  space_op.GetLumpedPortExcitationVectorPrimaryHtcn(1, port_primary_ht_cn);
  GridFunction port_primary_gf_ht_cn(space_op.GetNDSpace());  // As palace GridFunction
  port_primary_gf_ht_cn = 0.0;
  port_primary_gf_ht_cn.Real().SetFromTrueDofs(port_primary_ht_cn.Real());

  ComplexVector port_primary_et;
  space_op.GetLumpedPortExcitationVectorPrimaryEt(1, port_primary_et);
  GridFunction port_primary_gf_et(space_op.GetNDSpace());  // As palace GridFunction
  port_primary_gf_et = 0.0;
  port_primary_gf_et.Real().SetFromTrueDofs(port_primary_et.Real());

  // Total power for the (real) fields is (e_t / eta) \cdot e_t = 1.
  //
  // Note: If PEC::SIDE then the power normalization is wrong, since we have just removed
  // some of the electric field without properly compensating in the integral. This defines
  // a new alpha normalization of the field that we want to divide out to return to the
  // original p_0 = 1.
  double port_power = g_boundary_id.ParInnerProduct(port_primary_gf_ht_cn.Real(),
                                                    port_primary_gf_et.Real());
  // double port_normalization_alpha = std::sqrt(port_power);  // Marks-Williams alpha
  if (boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(port_power, WithinRel(1.0));
  }

  // Plain normalization of et
  double norm_id_et =
      g_boundary_id.ParInnerProduct(port_primary_gf_et.Real(), port_primary_gf_et.Real());

  // The normalization of e_t is Z_R \sum_e  W_e / L_e but with Z_R = 1.0 in the
  // internal units.
  double et_norm_expected = port_1.GetExcitationFieldEtNormSqWithUnityZR();

  if (boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(norm_id_et, WithinRel(et_norm_expected));
  }
  else
  {
    // Need to add power correction factor, since the power != 1. However, we can't
    // compensate for this here since it is a multi-port situation, so expression below is
    // not quite accurate.
    //
    // CHECK_THAT(norm_id_et / (port_normalization_alpha * port_normalization_alpha),
    //            WithinRel(et_norm_expected));
  }

  // Plain normalization of h_t x n = e_t / eta
  double norm_id_ht_cn = g_boundary_id.ParInnerProduct(port_primary_gf_ht_cn.Real(),
                                                       port_primary_gf_ht_cn.Real());

  // The normalization of e_t / eta is 1 / (Z_R n_el^2) \sum_e L_e / W_e but with Z_R = 1.0
  // in the internal units, since we took out the 1 / sqrt(\vert Z_0 \vert) factor out of
  // e_t / eta.
  double ht_cn_norm_expected = port_1.GetExcitationFieldHtNormSqWithUnityZR();

  if (boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(norm_id_ht_cn, WithinRel(ht_cn_norm_expected));
  }
  else
  {
    // Need to add power correction factor, since the power != 1. However, we can't
    // compensate for this here since it is a multi-port situation, so expression below is
    // not quite accurate.
    //
    // CHECK_THAT(norm_id_ht_cn / (port_normalization_alpha * port_normalization_alpha),
    //            WithinRel(ht_cn_norm_expected));
  }

  double norm_mode_ht_cn = g_boundary_mode.ParInnerProduct(port_primary_gf_ht_cn.Real(),
                                                           port_primary_gf_ht_cn.Real());
  double norm_perp_ht_cn = g_boundary_perp.ParInnerProduct(port_primary_gf_ht_cn.Real(),
                                                           port_primary_gf_ht_cn.Real());

  // Validate port integration normalisation. This is always true for the single hex mesh.
  // or if there is no proximate metal. Tet mesh mixes directions.
  if (mesh_is_hex || boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(norm_id_ht_cn, WithinRel(norm_mode_ht_cn));
    CHECK_THAT(norm_perp_ht_cn, WithinAbs(0.0, 1e-12));
  }
  else
  {
    // If PEC side, then e_t is a non-uniform field pattern that depends on the details of
    // Nédélec tet elements. Here only check trivial statement g_id = g_mode + g_perp.
    CHECK_THAT(norm_id_ht_cn, WithinRel(norm_mode_ht_cn + norm_perp_ht_cn));
  }

  // Measure voltage.
  //
  // Note: On paper, we have normalized $\vert v_0 \vert^2 = \vert Z_R \vert$ on ports, with
  // the assumption of unit power.
  //
  // For the case with side metal, we have \alpha != 1. We need to divide out a factor alpha
  // to normalize e_t to the original p_0 = 1 field conventions.
  //
  // However, this is not sufficient! Because the e_t vector is now no longer equal the to
  // vector field in the LinearForm used in GetVoltage (which does not have PEC zeroed).
  // There are two options here. (1) Either redefine the meaning of voltage as an integral
  // with a e_t field with metal zeroed in, which is strange but internally consistent. (2)
  // Keep the voltage definition as is but this introduces another factor, that effectively
  // is a Marks-Williams \beta factor.
  //
  // Need to correct for different mode. However, this is not just the port power, since the
  // voltage integral along \hat{y} is no longer in the same direction everywhere as $e_t$.
  auto v_in_measure = port_1.GetVoltage(port_primary_gf_et);
  CHECK(v_in_measure.imag() == 0.0);
  double port_normalization_beta =
      iodata.units.Nondimensionalize<palace::Units::ValueType::VOLTAGE>(
          std::sqrt(electromagnetics::Z0_)) /
      (v_in_measure.real());
  if (boundary_pec_type != PEC::SIDE)
  {
    CHECK_THAT(port_normalization_beta, WithinRel(1.0));
  }

  // Now "s" form in LumpedPortData is just e_t / eta (with Z_R = R). So check inner
  // product.
  //
  // Note: s should not change during change of field normalization alpha, but formula
  // assumes 1.0 and does not normalized. Also, the integral seems wrong with metal, as the
  // SParameter Linear form should be the field form e_t / eta.
  if (boundary_pec_type != PEC::SIDE)
  {
    std::complex<double> s_param = port_1.GetSParameter(port_primary_gf_ht_cn);
    CHECK_THAT(s_param.real(),
               WithinRel(ht_cn_norm_expected /
                         iodata.units.Nondimensionalize<VT::VOLTAGE>(std::sqrt(50.))));
    CHECK_THAT(s_param.imag(), WithinAbs(0.0, 1e-12));
  }

  // Extract s and v linear forms to check equality. Use reinterpret cast to identical
  // layout Test structure from original structure. Done not to populate code with friend
  // functions, but still access private members for testing.
  const auto &port_1_test_cast = reinterpret_cast<const LumpedPortDataTest &>(port_1);
  port_1_test_cast.InitializeLinearFormsForTest(space_op.GetNDSpace().Get());
  auto *form_s = port_1_test_cast.GetLinearFormS();
  ComplexVector VecFormS;
  VecFormS.SetSize(space_op.GetNDSpace().GetTrueVSize());
  VecFormS.UseDevice(true);
  space_op.GetNDSpace().GetProlongationMatrix()->MultTranspose(*form_s, VecFormS.Real());

  auto *form_v = port_1_test_cast.GetLinearFormV();
  ComplexVector VecFormV;
  VecFormV.SetSize(space_op.GetNDSpace().GetTrueVSize());
  VecFormV.UseDevice(true);
  space_op.GetNDSpace().GetProlongationMatrix()->MultTranspose(*form_v, VecFormV.Real());

  // Now GetExcitationVector1 in space op, returns the dual of 2 e_t / eta (with Z_R = R).
  ComplexVector RHS;
  space_op.GetExcitationVector1(1, RHS);
  RHS *= 0.5;

  // In the case where there is a neighbouring PEC condition, RHS and VecFormS are no
  // longer equal — the DoF on the shared edge with the PEC neighbours are zeroed.
  const auto &nd_dbc_tdof = space_op.GetNDDbcTDofLists().back();
  if (boundary_pec_type != PEC::NONE)
  {
    auto nd_dbc_tdof_size = nd_dbc_tdof.Size();
    Mpi::GlobalSum(1, &nd_dbc_tdof_size, world_comm);
    CHECK(nd_dbc_tdof_size > 0);
  }

  REQUIRE(RHS.Real().Size() == VecFormS.Real().Size());
  // This should true up to normalization, but is not currently due to non-subtracted metal.
  for (int i = 0; i < RHS.Real().Size(); i++)
  {
    // If subtracted metal, check RHS is zero and continue. Tet mesh also subtracts active
    // dofs if PEC is at the end of the port.
    if (boundary_pec_type == PEC::SIDE || (!mesh_is_hex && boundary_pec_type == PEC::ENDS))
    {
      const auto *find_pec = std::ranges::find(nd_dbc_tdof, i);
      if (find_pec != nd_dbc_tdof.end())
      {
        CHECK(RHS.Real()[i] == 0.0);
        continue;
      }
    }
    CHECK_THAT(RHS.Real()[i], WithinAbs(VecFormS.Real()[i], 1e-12));
  }

  // Check that form_s and form_v are the same apart from the factor 1 / sqrt(R). In actual
  // fact, this should only true in the ideal case, without any metal subtraction. Here it
  // is always true because of the form definitions.
  REQUIRE(VecFormS.Real().Size() == VecFormV.Real().Size());
  for (int i = 0; i < VecFormS.Real().Size(); i++)
  {
    CHECK_THAT(VecFormV.Real()[i],
               WithinAbs(VecFormS.Real()[i] *
                             iodata.units.Nondimensionalize<VT::VOLTAGE>(std::sqrt(50.)),
                         1e-12));
  }
}

TEST_CASE("LumpedPort_ReactiveExcitation_Cube321", "[lumped_port][Serial][Parallel]")
{
  // Exercise a reactive EXCITED lumped port — a port driven through a lumped inductor /
  // capacitor (e.g. a Josephson-junction port). Historically an excited port had to be
  // purely resistive; it may now carry reactance (including a purely reactive R = 0 port).
  // The reactance enters the system matrix as a physical termination; the incident-drive
  // normalization references a real resistance (the port R, or the unit internal reference
  // when R = 0) so a purely reactive drive does not divide by zero. This test verifies:
  //   (1) a reactive excited port constructs without tripping the guardrail,
  //   (2) GetExcitationRefResistance() returns R when R > 0, else the unit reference,
  //   (3) GetCharacteristicImpedance() returns the correct complex Z_ref(w) =
  //   R||iwL||1/iwC, (4) the excitation RHS is finite (no NaN/Inf) for a purely reactive R
  //   = 0 drive.
  using VT = palace::Units::ValueType;
  MPI_Comm world_comm = Mpi::World();

  const auto &[case_R, case_L, case_C] =
      GENERATE(std::make_tuple(50.0, 1.0e-9, 0.0),   // R + L (well-defined S)
               std::make_tuple(0.0, 1.0e-9, 0.0),    // pure inductor (R = 0)
               std::make_tuple(0.0, 0.0, 1.0e-12));  // pure capacitor (R = 0)

  double L0 = 1.0e-6;
  double Lc = 7.0;
  IoData iodata{Units(L0, Lc)};
  iodata.model.mesh =
      fs::path(PALACE_TEST_DATA_DIR) / "lumpedport_mesh/cube_mesh_3_2_1_hex.msh";
  iodata.model.L0 = L0;
  iodata.model.Lc = Lc;
  iodata.model.crack_bdr_elements = false;
  iodata.solver.order = 2;

  json domains_json = {
      {"Materials",
       json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                  {"Permeability", 1.0},
                                  {"Permittivity", 1.0},
                                  {"LossTan", 0.0}})})}};
  iodata.domains = config::DomainData(domains_json);

  // Build the excited reactive port. Only include nonzero circuit elements.
  json port = {{"Index", 1},
               {"Excitation", uint(1)},
               {"Attributes", json::array({1})},
               {"Direction", "+Y"}};
  if (case_R > 0.0)
  {
    port["R"] = case_R;
  }
  if (case_L > 0.0)
  {
    port["L"] = case_L;
  }
  if (case_C > 0.0)
  {
    port["C"] = case_C;
  }
  json boundary_json = {{"LumpedPort", json::array({port})}};
  // (1) Construction must succeed — the R >= 0 guardrail admits reactance.
  iodata.boundaries = config::BoundaryData(boundary_json);
  iodata.CheckConfiguration();

  auto mesh_io = LoadScaleParMesh(iodata, world_comm);
  SpaceOperator space_op(iodata, mesh_io);
  const auto &port_1 = space_op.GetLumpedPortOp().GetPort(1);

  CHECK(port_1.HasExcitation());

  // (2) Reference resistance: R when R > 0, else the unit internal reference (= 1 in
  // internal units, i.e. Z_freespace dimensionally).
  const double R_nondim = iodata.units.Nondimensionalize<VT::IMPEDANCE>(case_R);
  if (case_R > 0.0)
  {
    CHECK_THAT(port_1.GetExcitationRefResistance(), WithinRel(R_nondim));
  }
  else
  {
    CHECK_THAT(port_1.GetExcitationRefResistance(), WithinRel(1.0));
  }

  // (3) Characteristic impedance Z_ref(w) = 1 / (1/R + 1/(iwL) + iwC). Check against a
  // direct evaluation at a representative (nondimensional) frequency.
  const double f_GHz = 10.0;
  const double omega =
      2.0 * std::numbers::pi * iodata.units.Nondimensionalize<VT::FREQUENCY>(f_GHz);
  const std::complex<double> imag_unit{0.0, 1.0};
  std::complex<double> Y = 0.0;
  if (case_R > 0.0)
  {
    Y += 1.0 / R_nondim;
  }
  if (case_L > 0.0)
  {
    Y += 1.0 / (imag_unit * omega * iodata.units.Nondimensionalize<VT::INDUCTANCE>(case_L));
  }
  if (case_C > 0.0)
  {
    Y += imag_unit * omega * iodata.units.Nondimensionalize<VT::CAPACITANCE>(case_C);
  }
  const std::complex<double> Z_expect = 1.0 / Y;
  const std::complex<double> Z_ref =
      port_1.GetCharacteristicImpedance(omega, LumpedPortData::Branch::TOTAL);
  CHECK_THAT(std::abs(Z_ref - Z_expect), WithinAbs(0.0, 1e-9 * std::abs(Z_expect)));

  // (4) The excitation RHS must be finite even for a purely reactive R = 0 drive (the key
  // fix: the incident-field normalization references GetExcitationRefResistance(), not R,
  // so it does not divide by zero).
  ComplexVector RHS;
  space_op.GetExcitationVector1(1, RHS);
  double max_abs = 0.0;
  for (int i = 0; i < RHS.Real().Size(); i++)
  {
    CHECK(std::isfinite(RHS.Real()[i]));
    max_abs = std::max(max_abs, std::abs(RHS.Real()[i]));
  }
  Mpi::GlobalMax(1, &max_abs, world_comm);
  // The drive is nonzero (a genuine incident field was assembled).
  CHECK(max_abs > 0.0);

  // (5) Output paths reference the same real resistance. GetExcitationVoltage() must be
  // finite and nonzero for R = 0 (it references GetExcitationRefResistance(), not R), and
  // the S-parameter projection linear form must produce a finite, nonzero projection for
  // a nonzero field (its normalization also references the real reference resistance).
  const double V_inc = port_1.GetExcitationVoltage();
  CHECK(std::isfinite(V_inc));
  CHECK(V_inc > 0.0);

  GridFunction E_test(space_op.GetNDSpace(), true);
  E_test.Real() = 1.0;  // Uniform nonzero field: projection onto the port mode is nonzero
  E_test.Imag() = 0.0;
  const std::complex<double> S_proj = port_1.GetSParameter(E_test);
  CHECK(std::isfinite(S_proj.real()));
  CHECK(std::isfinite(S_proj.imag()));
  CHECK(std::abs(S_proj) > 0.0);
}

TEST_CASE("LumpedPort excitation weight from amplitude and phase", "[lumped_port][Serial]")
{
  using C = std::complex<double>;

  // Multiples of 90 degrees are exact, so a weight which is real (or imaginary) in exact
  // arithmetic is also real (or imaginary) in floating point.
  CHECK(ExcitationWeight(2.0, 0.0) == C(2.0, 0.0));
  CHECK(ExcitationWeight(2.0, 90.0) == C(0.0, 2.0));
  CHECK(ExcitationWeight(2.0, 180.0) == C(-2.0, 0.0));
  CHECK(ExcitationWeight(2.0, 270.0) == C(0.0, -2.0));
  CHECK(ExcitationWeight(2.0, -90.0) == C(0.0, -2.0));
  CHECK(ExcitationWeight(2.0, -180.0) == C(-2.0, 0.0));
  CHECK(ExcitationWeight(2.0, 360.0) == C(2.0, 0.0));
  CHECK(ExcitationWeight(2.0, -360.0) == C(2.0, 0.0));
  CHECK(ExcitationWeight(2.0, 810.0) == C(0.0, 2.0));
  CHECK(ExcitationWeight(2.0, -1.0e-20) == C(2.0, 0.0));
  CHECK(ExcitationWeight(2.0, 180.0).imag() == 0.0);
  CHECK(ExcitationWeight(3.0, 90.0).real() == 0.0);

  // Other phases, reduced modulo 360 degrees.
  for (const double phase : {30.0, -110.0, 75.0, 359.0})
  {
    CAPTURE(phase);
    const C expected = std::polar(0.5, phase * std::numbers::pi / 180.0);
    CHECK(std::abs(ExcitationWeight(0.5, phase) - expected) <= 1.0e-15);
    CHECK(std::abs(ExcitationWeight(0.5, phase + 360.0) - expected) <= 1.0e-14);
    CHECK(std::abs(ExcitationWeight(0.5, phase - 720.0) - expected) <= 1.0e-14);
  }
}

TEST_CASE("LumpedPort weighted excitation vector", "[lumped_port][Serial][Parallel]")
{
  MPI_Comm comm = Mpi::World();
  auto [mesh_is_hex, mesh_path] =
      GENERATE(std::make_tuple(true, fs::path(PALACE_TEST_DATA_DIR) /
                                         "lumpedport_mesh/cube_mesh_3_2_1_hex.msh"),
               std::make_tuple(false, fs::path(PALACE_TEST_DATA_DIR) /
                                          "lumpedport_mesh/cube_mesh_3_2_1_tet.msh"));
  auto order = GENERATE(1, 2);
  CAPTURE(mesh_is_hex, order);
  const fs::path mesh_file = mesh_path;
  const int solver_order = order;

  // Port 1 is a single element port on attribute 1 and port 2 is a multielement port, as
  // in the other lumped port tests. The extra JSON entries override the defaults and add
  // the excitation options.
  auto MakePort1 = [](const json &extra)
  {
    json port = {{"Index", 1},
                 {"R", 50.0},
                 {"Attributes", json::array({1})},
                 {"Direction", "+Y"}};
    port.update(extra);
    return port;
  };
  auto MakePort2 = [](const json &extra)
  {
    json port = {
        {"Index", 2},
        {"R", 75.0},
        {"Elements",
         json::array({json::object({{"Attributes", json::array({9, 11})},
                                    {"Direction", "+Y"}}),
                      json::object({{"Attributes", json::array({2, 6, 10})},
                                    {"Direction", "+X"}})})}};
    port.update(extra);
    return port;
  };

  // Excitation vectors and incident quantities of the ports of one configuration.
  struct Drive
  {
    ComplexVector rhs1;
    Vector rhs_time;
    std::array<std::complex<double>, 2> weight, v_inc, i_inc;
    std::array<double, 2> v_exc, p_exc;
  };
  auto Assemble = [&](const json &port1, const json &port2, int excitation_idx,
                      bool time_domain)
  {
    IoData iodata{Units(1.0e-6, 7.0)};
    iodata.model.mesh = mesh_file;
    iodata.model.L0 = 1.0e-6;
    iodata.model.Lc = 7.0;
    iodata.model.crack_bdr_elements = false;
    iodata.solver.order = solver_order;
    json domains_json = {
        {"Materials",
         json::array({json::object({{"Attributes", json::array({1, 2, 3, 4, 5, 6})},
                                    {"Permeability", 1.0},
                                    {"Permittivity", 1.0},
                                    {"LossTan", 0.0}})})}};
    iodata.domains = config::DomainData(domains_json);
    json boundaries_json = {
        {"LumpedPort", json::array({MakePort1(port1), MakePort2(port2)})}};
    iodata.boundaries = config::BoundaryData(boundaries_json);
    iodata.CheckConfiguration();

    auto mesh = LoadScaleParMesh(iodata, comm);
    SpaceOperator space_op(iodata, mesh);

    Drive out;
    out.rhs1.UseDevice(true);
    REQUIRE(space_op.GetExcitationVector1(excitation_idx, out.rhs1));
    if (time_domain)
    {
      REQUIRE(space_op.GetExcitationVector(excitation_idx, out.rhs_time));
    }
    for (int i = 0; i < 2; i++)
    {
      const auto &port = space_op.GetLumpedPortOp().GetPort(i + 1);
      out.weight[i] = port.GetExcitationWeight();
      out.v_inc[i] = port.GetIncidentVoltage();
      out.i_inc[i] = port.GetIncidentCurrent();
      out.v_exc[i] = port.GetExcitationVoltage();
      out.p_exc[i] = port.GetExcitationPower();
    }
    return out;
  };

  const std::complex<double> a1 = std::polar(2.0, 30.0 * std::numbers::pi / 180.0);
  const std::complex<double> a2 = std::polar(0.5, -110.0 * std::numbers::pi / 180.0);
  const json no_port = json::object();

  // Each port driven alone with the default weight, as references.
  const auto drive_1 = Assemble(json{{"Excitation", 1}}, no_port, 1, true);
  const auto drive_2 = Assemble(no_port, json{{"Excitation", 2}}, 2, false);
  const double norm_1 = linalg::Norml2(comm, drive_1.rhs1);
  const double norm_2 = linalg::Norml2(comm, drive_2.rhs1);
  REQUIRE(norm_1 > 0.0);
  REQUIRE(norm_2 > 0.0);
  CHECK(linalg::Norml2(comm, drive_1.rhs1.Imag()) == 0.0);
  CHECK(linalg::Norml2(comm, drive_2.rhs1.Imag()) == 0.0);

  // The weights default to one for the excited port and zero for the other.
  CHECK(drive_1.weight[0] == std::complex<double>(1.0, 0.0));
  CHECK(drive_1.weight[1] == std::complex<double>(0.0, 0.0));
  CHECK(drive_1.v_inc[1] == std::complex<double>(0.0, 0.0));
  CHECK(drive_1.i_inc[1] == std::complex<double>(0.0, 0.0));
  CHECK(drive_1.p_exc[1] == 0.0);

  SECTION("Both ports in one excitation add with their complex weights")
  {
    const auto drive_w = Assemble(
        json{{"Excitation", 1}, {"ExcitationAmplitude", 2.0}, {"ExcitationPhase", 30.0}},
        json{{"Excitation", 1}, {"ExcitationAmplitude", 0.5}, {"ExcitationPhase", -110.0}},
        1, false);
    const double norm_w = linalg::Norml2(comm, drive_w.rhs1);
    REQUIRE(norm_w > 0.0);
    CHECK(linalg::Norml2(comm, drive_w.rhs1.Imag()) > 0.0);

    // RHS1 = a1 RHS1(port 1 alone) + a2 RHS1(port 2 alone).
    ComplexVector diff(drive_w.rhs1);
    diff.AXPY(-a1, drive_1.rhs1);
    diff.AXPY(-a2, drive_2.rhs1);
    CHECK(linalg::Norml2(comm, diff) <= 1.0e-12 * norm_w);

    // The weights are applied to the excitation and the incident quantities, but not to
    // the unit-normalized ones.
    CHECK(std::abs(drive_w.weight[0] - a1) <= 1.0e-15);
    CHECK(std::abs(drive_w.weight[1] - a2) <= 1.0e-15);
    CHECK(drive_w.v_exc[0] == drive_1.v_exc[0]);
    CHECK(drive_w.v_exc[1] == drive_2.v_exc[1]);
    CHECK(drive_w.p_exc[0] == 1.0);
    CHECK(drive_w.p_exc[1] == 1.0);
    const std::array<std::complex<double>, 2> a = {a1, a2};
    for (int i = 0; i < 2; i++)
    {
      CAPTURE(i);
      REQUIRE(drive_w.v_exc[i] > 0.0);
      CHECK(std::abs(drive_w.v_inc[i] - a[i] * drive_w.v_exc[i]) <=
            1.0e-14 * std::abs(a[i]) * drive_w.v_exc[i]);
      CHECK(std::abs(drive_w.i_inc[i] - a[i] / drive_w.v_exc[i]) <=
            1.0e-14 * std::abs(a[i]) / drive_w.v_exc[i]);
    }
  }

  SECTION("A phase of 180 degrees is a flipped direction")
  {
    const auto drive_180 = Assemble(
        json{{"Excitation", 1}, {"ExcitationPhase", 180.0}}, no_port, 1, false);
    const auto drive_flip =
        Assemble(json{{"Excitation", 1}, {"Direction", "-Y"}}, no_port, 1, false);

    // The weight is exactly real, so the imaginary part is not assembled.
    CHECK(drive_180.weight[0] == std::complex<double>(-1.0, 0.0));
    CHECK(linalg::Norml2(comm, drive_180.rhs1.Imag()) == 0.0);

    ComplexVector diff(drive_180.rhs1);
    diff -= drive_flip.rhs1;
    CHECK(linalg::Norml2(comm, diff) <= 1.0e-12 * norm_1);

    // Also the negative of the unweighted vector.
    diff.Set(drive_180.rhs1);
    diff += drive_1.rhs1;
    CHECK(linalg::Norml2(comm, diff) <= 1.0e-12 * norm_1);
  }

  SECTION("The amplitude scales the real and time domain excitation")
  {
    const auto drive_3 = Assemble(json{{"Excitation", 1}, {"ExcitationAmplitude", 3.0}},
                                  no_port, 1, true);
    CHECK(drive_3.weight[0] == std::complex<double>(3.0, 0.0));
    CHECK(linalg::Norml2(comm, drive_3.rhs1.Imag()) == 0.0);
    CHECK(drive_3.p_exc[0] == 1.0);
    CHECK(drive_3.v_exc[0] == drive_1.v_exc[0]);
    CHECK(std::abs(drive_3.v_inc[0] - 3.0 * drive_1.v_exc[0]) <=
          1.0e-14 * 3.0 * drive_1.v_exc[0]);

    // Frequency domain RHS1.
    ComplexVector diff(drive_3.rhs1);
    diff.AXPY(-3.0, drive_1.rhs1);
    CHECK(linalg::Norml2(comm, diff) <= 1.0e-12 * norm_1);

    // Time domain vector, which also equals the real part of RHS1.
    REQUIRE(drive_3.rhs_time.Size() == drive_1.rhs_time.Size());
    Vector diff_time(drive_3.rhs_time);
    diff_time.Add(-3.0, drive_1.rhs_time);
    CHECK(linalg::Norml2(comm, diff_time) <= 1.0e-12 * norm_1);
    Vector diff_real(drive_3.rhs_time);
    diff_real -= drive_3.rhs1.Real();
    CHECK(linalg::Norml2(comm, diff_real) <= 1.0e-12 * norm_1);
  }
}
