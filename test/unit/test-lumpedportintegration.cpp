// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

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
  mfem::LinearForm *GetLinearFormS() const { return s.get(); }
  mfem::LinearForm *GetLinearFormV() const { return v.get(); }
};

auto LoadScaleParMesh(IoData &iodata, MPI_Comm world_comm)
{
  // Load Mesh — copy from main.cpp
  std::vector<std::unique_ptr<Mesh>> mesh_;
  {
    std::vector<std::unique_ptr<mfem::ParMesh>> mfem_mesh;
    mfem_mesh.push_back(mesh::ReadMesh(iodata, world_comm));
    iodata.NondimensionalizeInputs(*mfem_mesh[0]);
    mesh::RefineMesh(iodata, mfem_mesh);
    for (auto &m : mfem_mesh)
    {
      mesh_.push_back(std::make_unique<Mesh>(std::move(m)));
    }
  }
  return mesh_;
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

    CHECK(bbox.Lengths().size() == 3);
    std::array<double, 3> bbox_lengths_out = {length_dx_m, length_dx_m, 0.0};
    CHECK(bbox.Lengths() == bbox_lengths_out);
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
  space_op.GetLumpedPortExcitationVectorPrimaryHtcn(1, port_primary_ht_cn, true);
  // As palace GridFunction: Undo restriction
  GridFunction port_primary_gf_ht_cn(space_op.GetNDSpace());
  port_primary_gf_ht_cn = 0.0;
  port_primary_gf_ht_cn.Real().SetFromTrueDofs(port_primary_ht_cn.Real());

  ComplexVector port_primary_et;
  port_primary_et.UseDevice(true);
  space_op.GetLumpedPortExcitationVectorPrimaryEt(1, port_primary_et, true);
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
      const auto *find_pec = std::find(nd_dbc_tdof.begin(), nd_dbc_tdof.end(), i);
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
  space_op.GetLumpedPortExcitationVectorPrimaryHtcn(1, port_primary_ht_cn, true);
  GridFunction port_primary_gf_ht_cn(space_op.GetNDSpace());  // As palace GridFunction
  port_primary_gf_ht_cn = 0.0;
  port_primary_gf_ht_cn.Real().SetFromTrueDofs(port_primary_ht_cn.Real());

  ComplexVector port_primary_et;
  space_op.GetLumpedPortExcitationVectorPrimaryEt(1, port_primary_et, true);
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
      const auto *find_pec = std::find(nd_dbc_tdof.begin(), nd_dbc_tdof.end(), i);
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
