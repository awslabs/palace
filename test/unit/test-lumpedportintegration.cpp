#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/src/Core/IO.h>
#include <fmt/format.h>
#include <catch2/catch_approx.hpp>
#include <catch2/catch_session.hpp>
#include <catch2/catch_test_macros.hpp>
#include <nlohmann/json.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <catch2/matchers/catch_matchers_vector.hpp>
#include "catch2/matchers/catch_matchers_floating_point.hpp"
#include "drivers/drivensolver.hpp"
#include "fem/bilinearform.hpp"
#include "fem/fespace.hpp"
#include "fem/gridfunction.hpp"
#include "fem/integrator.hpp"
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
#include "utils/meshio.hpp"
#include "utils/omp.hpp"
#include "utils/outputdir.hpp"

using namespace palace;

void PrintViaDenseEigen(const mfem::BilinearForm &m_ij, double round)
{
  // Debug Printing Only
  if constexpr (true)
  {
    return;
  }
  auto n = m_ij.NumRows();
  auto m = m_ij.NumCols();

  Eigen::MatrixXd eigen_m_ij = Eigen::MatrixXd::Zero(n, m);

  Vector x(n), y(m);
  x = 0.0;
  y = 0.0;

  for (int i = 0; i < n; i++)
  {
    x(i) = 1.0;
    m_ij.Mult(x, y);
    for (int j = 0; j < m; j++)
    {
      if (std::abs(y(j)) > round)
      {
        eigen_m_ij(i, j) = y(j);
      }
    }
    x(i) = 0.0;
  }

  Eigen::IOFormat HeavyFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "[", "]", "[", "]");
  std::cout << eigen_m_ij.format(HeavyFmt) << std::endl;
}

void run_lumped_port_case(std::string_view name, IoData iodata, MPI_Comm world_comm)
{
  MakeOutputFolder(iodata, world_comm);
  std::string path_loc = fmt::format("{}/phd_lumpedport_mesh/", PALACE_TEST_DIR);

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

  SpaceOperator space_op{iodata, mesh_};
  auto &lumped_port_op = space_op.GetLumpedPortOp();
  const int port_idx = 1;
  const auto &port_1 = lumped_port_op.GetPort(port_idx);

  mfem::Array<int> port_attr_list;

  auto &mesh = space_op.GetNDSpace().GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;

  SumCoefficient fb_id{};
  SumVectorCoefficient fb_x(mesh.SpaceDimension());
  SumVectorCoefficient fb_y(mesh.SpaceDimension());

  SumVectorCoefficient fb_x_port(mesh.SpaceDimension());

  for (std::size_t elem_i = 0; elem_i < port_1.elems.size(); elem_i++)
  {
    const auto &elem = port_1.elems[elem_i];

    // Print Port Information:
    Mpi::Print("Port[{}], sub-element[{}] information:", port_idx, elem_i);
    Mpi::Print("..Width {}\n.. Length {}, Nr {}\n", elem->GetGeometryWidth(),
               elem->GetGeometryLength(), port_1.elems.size());
    Mpi::Print("..AttrList {}\n", elem->GetAttrList());

    auto bbox = mesh::GetBoundingBox(mesh, elem->GetAttrList()[0], true);
    Mpi::Print("Attr bounding boxL\n .. (center): {}\n", bbox.center);
    Mpi::Print(".. (axes): {}n\n", bbox.axes);
    Mpi::Print(".. (planar): {}\n", bbox.planar);
    Mpi::Print(".. (area): {}\n", bbox.Area());
    Mpi::Print(".. (volume): {}\n", bbox.Volume());
    Mpi::Print(".. (lengths): {}\n", bbox.Lengths());
    Mpi::Print(".. (Normals): {}\n", bbox.Normals());
    Mpi::Print(".. (Deviations): {}", bbox.Deviations());

    auto port_element_scale_factor =
        1.0 / std::sqrt(elem->GetGeometryWidth() * elem->GetGeometryLength() *
                        port_1.elems.size());
    Mpi::Print(".. port_element_scale_factor {}\n", port_element_scale_factor);

    port_attr_list.Append(elem->GetAttrList());

    mfem::Array<int> attr_marker_loc;
    mesh::AttrToMarker(bdr_attr_max, elem->GetAttrList(), attr_marker_loc);
    Mpi::Print("..AttrMarker {}\n", attr_marker_loc);

    auto surface_area = mesh::GetSurfaceArea(mesh, attr_marker_loc);
    Mpi::Print("Mesh GetSurfaceArea: {}\n", surface_area);

    // Make identity coefficient restricted to the boundary which is ~ \delta_{i,j}
    fb_id.AddCoefficient(std::make_unique<RestrictedCoefficient<mfem::ConstantCoefficient>>(
        elem->GetAttrList(), 1.0));

    // Make identity coefficient restricted to the boundary which is ~ \hat{x}_i \hat{x}_j
    StaticVector<3> x_dir;
    x_dir = 0.0;
    x_dir(0) = 1.0;
    fb_x.AddCoefficient(
        std::make_unique<RestrictedVectorCoefficient<mfem::VectorConstantCoefficient>>(
            port_attr_list, x_dir));

    // Make identity coefficient restricted to the boundary which is ~ \hat{y}_i \hat{y}_j
    StaticVector<3> y_dir;
    y_dir = 0.0;
    y_dir(1) = 1.0;
    fb_y.AddCoefficient(
        std::make_unique<RestrictedVectorCoefficient<mfem::VectorConstantCoefficient>>(
            port_attr_list, y_dir));

    // Make scaled coefficient: this should be the same as fb_x but normalizes by
    // port_element_scale_factor (which should work).
    fb_x_port.AddCoefficient(elem->GetModeCoefficient(port_element_scale_factor));
  }

  // Total Port Attr Marker
  mfem::Array<int> port_attr_marker;
  mesh::AttrToMarker(bdr_attr_max, port_attr_list, port_attr_marker);
  Mpi::Print("Combined port attr markers: {}\n", port_attr_marker);

  Mpi::Print("{:*^30}\n", "");

  // GridFunction: contravariant (primary) vector epb^i. Project on tangential boundary.
  mfem::GridFunction e_primary_b_x(&space_op.GetNDSpace().Get());
  e_primary_b_x = 0.0;
  e_primary_b_x.ProjectBdrCoefficientTangent(fb_x_port, port_attr_marker);

  // LinearForm: covariant (dual) vector edb_i from boundary integrator.
  mfem::LinearForm e_dual_b_x(&space_op.GetNDSpace().Get());
  e_dual_b_x.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(fb_x_port),
                                   port_attr_marker);
  e_dual_b_x.Assemble();

  // Print
  {
    Mpi::Print("GridFunction e_primary_b_x: \n");
    e_primary_b_x.Print(std::cout, 64);
    Mpi::Print("LinearForm e_dual_b_x: \n");
    e_dual_b_x.Print(std::cout, 64);
  }

  // By construction primary filed and dual field should be same on inner product space of
  // boundary edb_i epb^i = 1. (Up to the port_element_scale_factor which is 1 here).
  {
    double dot(e_dual_b_x * e_primary_b_x);
    Mpi::GlobalSum(1, &dot, world_comm);
    Mpi::Print("e_dual_b_x * e_primary_b_x: {:e}\n", dot);
    CHECK_THAT(e_dual_b_x * e_primary_b_x, Catch::Matchers::WithinULP(1.0, 4));
  }

  Mpi::Print("{:*^30}\n", "");

  GridFunction palace_gridfunciton_e_primary_b_x(space_op.GetNDSpace());
  palace_gridfunciton_e_primary_b_x.Real() = e_primary_b_x;
  palace_gridfunciton_e_primary_b_x.Imag() = 0.0;

  auto V_out = port_1.GetVoltage(palace_gridfunciton_e_primary_b_x);
  Mpi::Print("Lumped Port GetVoltage on palace_gridfunciton_e_primary_b_x: {:e}, {:e}\n",
             V_out.real(), V_out.imag());

  Mpi::Print("{:*^30}\n", "");

  // Now let us assemble various bilinear forms:

  // Domain mass matrix g_domain_id_IJ
  mfem::BilinearForm g_domain_id(&space_op.GetNDSpace().Get());
  g_domain_id.AddDomainIntegrator(new mfem::VectorFEMassIntegrator());
  g_domain_id.Assemble(true);

  {
    std::ofstream myfile;
    myfile.open(fmt::format("{}/{}_g_domain_id.txt", path_loc, name));
    g_domain_id.PrintMatlab(myfile);
    myfile.close();
    Mpi::Print("g_domain_id:\n");
    PrintViaDenseEigen(g_domain_id, 1.0e-12);
    // Minor Q: Why are there so many tiny rounding values even thought this is a trivial
    // cube element and basically the reference element?
  }

  Mpi::Print("{:*^30}\n", "");

  // // Domain mass matrix g_domain_id_IJ restricted to the boundary port attr only
  // // Presumably, this does not work since port marker is in boundary marker list
  // // not bulk marker list.
  // mfem::BilinearForm g_domain_id_port_restricted(&space_op.GetNDSpace().Get());
  // g_domain_id_port_restricted.AddDomainIntegrator(new mfem::VectorFEMassIntegrator(),
  //                                                 port_attr_marker);
  // g_domain_id_port_restricted.Assemble();
  // {
  //   std::ofstream myfile;
  //   myfile.open(fmt::format("{}/{}_g_bulk_port_id_attr.txt", path_loc, name));
  //   g_domain_id_port_restricted.PrintMatlab(myfile);
  //   myfile.close();
  //   Mpi::Print("g_domain_id_port_restricted:\n");
  //   PrintViaDenseEigen(g_domain_id_port_restricted, 1.0e-12);
  // }

  // Boundary mass matrix g_boundary_id_ij with delta_ij coefficient
  mfem::BilinearForm g_boundary_id(&space_op.GetNDSpace().Get());
  g_boundary_id.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_id));
  g_boundary_id.Assemble(true);

  {
    std::ofstream myfile;
    myfile.open(fmt::format("{}/{}_g_boundary_id.txt", path_loc, name));
    g_boundary_id.PrintMatlab(myfile);
    myfile.close();
    Mpi::Print("g_boundary_id:\n");
    PrintViaDenseEigen(g_boundary_id, 1.0e-12);
  }

  Mpi::Print("{:*^30}\n", "");

  // Boundary mass matrix g_boundary_x_ij with vector x coefficient
  mfem::BilinearForm g_boundary_x(&space_op.GetNDSpace().Get());
  g_boundary_x.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_x));
  g_boundary_x.Assemble(true);

  {
    std::ofstream myfile;
    myfile.open(fmt::format("{}/{}_g_boundary_x.txt", path_loc, name));
    g_boundary_x.PrintMatlab(myfile);
    myfile.close();
    Mpi::Print("g_boundary_x:\n");
    PrintViaDenseEigen(g_boundary_x, 1.0e-12);
  }

  Mpi::Print("{:*^30}\n", "");

  // Boundary mass matrix g_boundary_y_ij with vector y coefficient
  mfem::BilinearForm g_boundary_y(&space_op.GetNDSpace().Get());
  g_boundary_y.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(fb_y));
  g_boundary_y.Assemble(true);

  {
    std::ofstream myfile;
    myfile.open(fmt::format("{}/{}_g_boundary_y.txt", path_loc, name));
    g_boundary_y.PrintMatlab(myfile);
    myfile.close();
    Mpi::Print("g_boundary_y:\n");
    PrintViaDenseEigen(g_boundary_y, 1.0e-12);
  }

  Mpi::Print("{:*^30}\n", "");

  // Now start taking inner products with bi-linear forms with primary field e_primary_b_x
  // Domain
  {
    double dot_domain_id = g_domain_id.InnerProduct(e_primary_b_x, e_primary_b_x);
    Mpi::Print("Inner product with domain identity:\n");
    Mpi::Print(".. e_primary_b_x * g_domain_id * e_primary_b_x = {:e}\n", dot_domain_id);
  }
  // Boundary
  {
    double dot_boundary_id = g_boundary_id.InnerProduct(e_primary_b_x, e_primary_b_x);
    Mpi::Print("Inner product with boundary identity:\n");
    Mpi::Print(".. e_primary_b_x * g_boundary_id * e_primary_b_x = {:e}\n",
               dot_boundary_id);
    CHECK_THAT(dot_boundary_id, Catch::Matchers::WithinULP(1.0, 4));

    double dot_boundary_x = g_boundary_x.InnerProduct(e_primary_b_x, e_primary_b_x);
    Mpi::Print("Inner product with boundary x polarized:\n");
    Mpi::Print(".. e_primary_b_x * g_boundary_x * e_primary_b_x = {:e}\n", dot_boundary_x);
    CHECK_THAT(dot_boundary_x, Catch::Matchers::WithinULP(1.0, 4));
    CHECK_THAT(dot_boundary_x, Catch::Matchers::WithinULP(dot_boundary_id, 4));

    double dot_boundary_y = g_boundary_y.InnerProduct(e_primary_b_x, e_primary_b_x);
    Mpi::Print("Inner product with boundary y polarized:\n");
    Mpi::Print(".. e_primary_b_x * g_boundary_y * e_primary_b_x = {:e}\n", dot_boundary_y);
    CHECK_THAT(dot_boundary_y, Catch::Matchers::WithinAbs(0.0, 1e-15));
  }

  // The fact that the domain and boundary inner product are different is due to the
  // presence of volume elements "broadening of the delta function". This is the effect we
  // have to deal with.

  Mpi::Print("{:*^30}\n", "");

  // Validate Linear Form as equivalent to g_ij e^j
  {
    const Vector &v_primary_b_x = e_primary_b_x.GetTrueVector();
    Vector v_form(v_primary_b_x.Size());
    v_form = 0.0;
    space_op.GetNDSpace().GetProlongationMatrix()->AddMultTranspose(e_dual_b_x, v_form);

    Vector v_form_test(v_primary_b_x.Size());
    g_boundary_id.Mult(v_primary_b_x, v_form_test);

    Mpi::Print("LinearForm e_dual_b_x as vector: \n");
    v_form.Print(std::cout, 64);
    Mpi::Print("Test bilinear-primary as matrix-vector g_boundary_id * e_primary_b_x :\n");
    v_form_test.Print(std::cout, 64);

    g_boundary_x.Mult(v_primary_b_x, v_form_test);
    Mpi::Print("Test bilinear-primary as matrix-vector g_boundary_x * e_primary_b_x :\n");
    v_form_test.Print(std::cout, 64);
  }

  Mpi::Print("{:*^30}\n", "");

  // Here were are going to "lift" by boundary-vector into the bulk, by performing
  // [g_bulk_id]^{-1} * g_boundary_id * e_primary_b_x with a linear solve.

  mfem::GridFunction e_primary_b_x_lifted(&space_op.GetNDSpace().Get());
  {
    e_primary_b_x_lifted = 0.0;

    // Form matrix linear system with MFEM machinery and solve. Don't impose essential BC.
    mfem::SparseMatrix A;
    mfem::Vector B, X;
    mfem::Array<int> ess_dof;
    g_domain_id.FormLinearSystem(ess_dof, e_primary_b_x_lifted, e_dual_b_x, A, X, B);
    mfem::CGSolver cg(MPI_COMM_WORLD);
    cg.SetRelTol(1e-14);
    cg.SetMaxIter(2000);
    cg.SetPrintLevel(1);
    cg.SetOperator(A);
    cg.Mult(B, X);
    g_domain_id.RecoverFEMSolution(X, e_primary_b_x, e_primary_b_x_lifted);

    Mpi::Print("\"Lifted\" primary boundary vector  [g_bulk_id]^(-1) * g_boundary_id * "
               "e_primary_b_x:\n");
    e_primary_b_x_lifted.Print(std::cout, 64);
  }

  // Inner products with lifted form
  {
    double dot_domain_id =
        g_domain_id.InnerProduct(e_primary_b_x_lifted, e_primary_b_x_lifted);
    Mpi::Print("Inner product with domain identity:\n");
    Mpi::Print(".. e_primary_b_x_lifted * g_domain_id * e_primary_b_x_lifted = {:e}\n",
               dot_domain_id);

    double dot_boundary_id =
        g_boundary_id.InnerProduct(e_primary_b_x_lifted, e_primary_b_x_lifted);
    Mpi::Print("Inner product with boundary identity:\n");
    Mpi::Print(".. e_primary_b_x_lifted * g_boundary_id * e_primary_b_x_lifted = {:e}\n",
               dot_boundary_id);

    double dot_boundary_x =
        g_boundary_x.InnerProduct(e_primary_b_x_lifted, e_primary_b_x_lifted);
    Mpi::Print("Inner product with boundary x polarized:\n");
    Mpi::Print(".. e_primary_b_x_lifted * g_boundary_x * e_primary_b_x_lifted = {:e}\n",
               dot_boundary_x);

    CHECK(dot_boundary_x == dot_boundary_id);

    // Believe this should be true?
    CHECK_THAT(dot_domain_id * dot_domain_id,
               Catch::Matchers::WithinRel(dot_boundary_x, 1e-12));

    double dot_boundary_y =
        g_boundary_y.InnerProduct(e_primary_b_x_lifted, e_primary_b_x_lifted);
    Mpi::Print("Inner product with boundary y polarized:\n");
    Mpi::Print(".. e_primary_b_x_lifted * g_boundary_y * e_primary_b_x_lifted = {:e}\n",
               dot_boundary_y);

    CHECK_THAT(dot_boundary_y, Catch::Matchers::WithinAbs(0.0, 1e-15));
  }

  // TODO:
  // Zero Metal Factor
  //   if (zero_metal)
  // {
  //   linalg::SetSubVector(RHS.Real(), nd_dbc_tdof_lists.back(), 0.0);
  // }
  // {
  //   std::complex<double> dot3(RHS.Real() * E_port_vec.Real());
  //   Mpi::GlobalSum(1, &dot3, GetComm());
  //   Mpi::Print("GetLumpedPortExcitationVector MetalZeroFactor: {:e}, {:e}\n",
  //   dot3.real(),
  //              dot3.imag());
  // }
}

// Run in serial only as not yet MPI safe
TEST_CASE("LumpedPortBasicTest_Hex111", "[driven_solver][Serial]")
{
  MPI_Comm world_comm = Mpi::World();
  auto palace_dir = fs::path(PALACE_TEST_DIR);

  std::string filename = palace_dir / "phd_lumpedport_mesh/cube_mesh_1_1_1.json";
  IoData iodata{filename.c_str(), false};

  iodata.model.mesh = palace_dir / "phd_lumpedport_mesh/cube_mesh_1_1_1_hex.msh";

  run_lumped_port_case("x_port", iodata, world_comm);
}

TEST_CASE("LumpedPortBasicTest_Hex321", "[driven_solver][Serial]")
{
  MPI_Comm world_comm = Mpi::World();
  auto palace_dir = fs::path(PALACE_TEST_DIR);

  std::string filename = palace_dir / "phd_lumpedport_mesh/cube_mesh_3_2_1.json";
  IoData iodata{filename.c_str(), false};

  iodata.model.mesh = palace_dir / "phd_lumpedport_mesh/cube_mesh_3_2_1_hex.msh";

  run_lumped_port_case("x_port", iodata, world_comm);
}

TEST_CASE("LumpedPortBasicTest_Tet111", "[driven_solver][Serial]")
{
  MPI_Comm world_comm = Mpi::World();
  auto palace_dir = fs::path(PALACE_TEST_DIR);

  std::string filename = palace_dir / "phd_lumpedport_mesh/cube_mesh_1_1_1.json";
  IoData iodata{filename.c_str(), false};

  // iodata.solver.order = 2;
  iodata.model.mesh = palace_dir / "phd_lumpedport_mesh/cube_mesh_1_1_1_tet.msh";

  run_lumped_port_case("x_port", iodata, world_comm);
}

TEST_CASE("LumpedPortBasicTest_Tet321", "[driven_solver][Serial]")
{
  MPI_Comm world_comm = Mpi::World();
  auto palace_dir = fs::path(PALACE_TEST_DIR);

  std::string filename = palace_dir / "phd_lumpedport_mesh/cube_mesh_1_1_1.json";
  IoData iodata{filename.c_str(), false};

  // iodata.solver.order = 2;
  iodata.model.mesh = palace_dir / "phd_lumpedport_mesh/cube_mesh_3_2_1_tet.msh";

  run_lumped_port_case("x_port", iodata, world_comm);
}

TEST_CASE("LumpedPortBasicTest_RefTet", "[driven_solver][Serial]")
{
  MPI_Comm world_comm = Mpi::World();
  auto palace_dir = fs::path(PALACE_TEST_DIR);

  std::string filename = palace_dir / "phd_lumpedport_mesh/ref-tetrahedron.json";
  IoData iodata{filename.c_str(), false};

  // iodata.solver.order = 2;
  iodata.model.mesh = palace_dir / "phd_lumpedport_mesh/ref-tetrahedron.mesh";

  run_lumped_port_case("x_port", iodata, world_comm);
}

// TODO Idea: wafer layer.
