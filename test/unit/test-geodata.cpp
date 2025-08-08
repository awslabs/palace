#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>

#include <vector>

#include "utils/geodata.hpp"
#include "utils/geodata_impl.hpp"

#include "fem/interpolator.hpp"
#include "utils/communication.hpp"
#include "utils/filesystem.hpp"

namespace palace
{
using namespace Catch;
namespace fs = std::filesystem;

TEST_CASE("TwoDimensionalDiagonalSquarePort", "[geodata]")
{
  std::vector<Eigen::Vector3d> vertices{
      {-0.19942181818181828, -0.5838274545454543, 0},
      {-0.19926108900667502, -0.5836667480978528, 0},
      {-0.19926061925486355, -0.5838279242308415, -8.56131891164518e-17},
      {-0.19926014950305207, -0.5839891003638303, 0},
      {-0.19910035983153176, -0.5835060416502512, 0},
      {-0.1990998900797203, -0.5836672177832399, -8.56131891164518e-17},
      {-0.1990994203279088, -0.5838283939162288, -1.3698110258632289e-16},
      {-0.19909895057609733, -0.5839895700492176, -8.56131891164518e-17},
      {-0.19909848082428586, -0.5841507461822064, 0},
      {-0.1989396306563885, -0.5833453352026496, 0},
      {-0.19893916090457703, -0.5835065113356385, -8.56131891164518e-17},
      {-0.19893869115276555, -0.5836676874686273, -1.3698110258632289e-16},
      {-0.19893822140095407, -0.583828863601616, -1.5410374040961326e-16},
      {-0.1989377516491426, -0.583990039734605, -1.3698110258632289e-16},
      {-0.19893728189733112, -0.5841512158675937, -8.56131891164518e-17},
      {-0.19893681214551964, -0.5843123920005826, 0},
      {-0.19877890148124525, -0.5831846287550482, 0},
      {-0.19877843172943377, -0.5833458048880369, -8.56131891164518e-17},
      {-0.1987779619776223, -0.5835069810210257, -1.3698110258632289e-16},
      {-0.19877749222581081, -0.5836681571540145, -1.5410374040961326e-16},
      {-0.19877702247399934, -0.5838293332870034, -1.3698110258632289e-16},
      {-0.19877655272218786, -0.5839905094199922, -1.5410374040961326e-16},
      {-0.19877608297037638, -0.584151685552981, -1.3698110258632289e-16},
      {-0.1987756132185649, -0.5843128616859699, -8.56131891164518e-17},
      {-0.19877514346675343, -0.5844740378189587, 0},
      {-0.19861817611093396, -0.5830239261117406, 0},
      {-0.1986153576000651, -0.5846338010914915, 0},
      {-0.19845746797462832, -0.583506036157253, -8.56131891164518e-17},
      {-0.19845745074062265, -0.582863223468433, 0},
      {-0.19845699822281687, -0.5836672122902418, -1.3698110258632289e-16},
      {-0.19845652969381722, -0.5831855725288573, -8.56131891164518e-17},
      {-0.19845652847100537, -0.5838283884232307, -1.5410374040961326e-16},
      {-0.1984560587191939, -0.5839895645562194, -1.3698110258632289e-16},
      {-0.1984555889673824, -0.5841507406892084, -1.5410374040961326e-16},
      {-0.19845557173337677, -0.5847935643640245, 0},
      {-0.19845511921557096, -0.5843119168221971, -1.3698110258632289e-16},
      {-0.19845464946375949, -0.5844730929551859, -8.56131891164518e-17},
      {-0.19829672537031137, -0.5827025208251255, 0},
      {-0.19829580432350594, -0.5830248698855497, -8.56131891164518e-17},
      {-0.19829578586668845, -0.5849533276365575, 0},
      {-0.19829486359707116, -0.5846328562277189, -8.56131891164518e-17},
      {-0.1981360344680114, -0.583827443559458, -1.3698110258632289e-16},
      {-0.1981360000000001, -0.5851130909090906, 0},
      {-0.1981360000000001, -0.5825418181818179, 0},
      {-0.19813556471619992, -0.5839886196924468, -1.5410374040961324e-16},
      {-0.19813509618720032, -0.5835069799310624, -1.3698110258632289e-16},
      {-0.19813509496438844, -0.5841497958254356, -1.3698110258632289e-16},
      {-0.19813507895319465, -0.5828641672422421, -8.56131891164518e-17},
      {-0.19813507773038283, -0.5847926195002519, -8.56131891164518e-17},
      {-0.198134625212577, -0.5843109719584244, -1.5410374040961324e-16},
      {-0.19813415790638922, -0.5831865163026665, -1.3698110258632289e-16},
      {-0.19813415546076552, -0.5844721480914132, -1.3698110258632289e-16},
      {-0.19797529186369447, -0.5849523827727849, -8.56131891164518e-17},
      {-0.1979743695940772, -0.5846319113639462, -1.3698110258632289e-16},
      {-0.19797435358288337, -0.5827034645989346, -8.56131891164518e-17},
      {-0.19797343253607796, -0.5830258136593589, -1.3698110258632289e-16},
      {-0.19781460096139447, -0.584148850961663, -1.5410374040961326e-16},
      {-0.19781458372738883, -0.5847916746364793, -1.3698110258632289e-16},
      {-0.197814131209583, -0.5843100270946517, -1.3698110258632289e-16},
      {-0.19781366268058337, -0.5838283873332673, -1.5410374040961326e-16},
      {-0.19781366145777152, -0.5844712032276406, -1.5410374040961326e-16},
      {-0.19781272439977232, -0.5835079237048715, -1.5410374040961324e-16},
      {-0.19781270716576665, -0.5828651110160513, -1.3698110258632289e-16},
      {-0.19781178611896122, -0.5831874600764757, -1.5410374040961326e-16},
      {-0.1976538755910832, -0.5846309665001737, -1.5410374040961326e-16},
      {-0.1976510607486499, -0.583026757433168, -1.5410374040961326e-16},
      {-0.19749316745477752, -0.584470258363868, -1.3698110258632289e-16},
      {-0.19749222917396644, -0.5841497947354722, -1.3698110258632289e-16},
      {-0.19749129089315537, -0.5838293311070765, -1.3698110258632289e-16},
      {-0.19749035261234427, -0.5835088674786807, -1.3698110258632289e-16},
      {-0.1974894143315332, -0.5831884038502848, -1.3698110258632289e-16},
      {-0.1973324665001741, -0.5843095574092645, -1.3698110258632289e-16},
      {-0.19733152821936303, -0.5839890937808687, -1.3698110258632289e-16},
      {-0.19733058993855196, -0.583668630152473, 1.3698110258632289e-16},
      {-0.19732965165774086, -0.5833481665240772, -1.3698110258632289e-16},
      {-0.19717176554557067, -0.5841488564546611, -1.3698110258632289e-16},
      {-0.1971708272647596, -0.5838283928262653, -1.3698110258632289e-16},
      {-0.1971698889839485, -0.5835079291978695, -1.3698110258632289e-16},
      {-0.1970110645909672, -0.5839881555000577, -1.3698110258632289e-16},
      {-0.19701012631015613, -0.5836676918716619, -1.3698110258632289e-16},
      {-0.19685036363636374, -0.5838274545454543, -1.3698110258632289e-16}};

  auto comm = Mpi::World();
  auto box = mesh::BoundingBoxFromPointCloud(comm, vertices, 0);

  // True box at 45 degrees
  auto invsqrt2 = 1.0 / std::sqrt(2);
  std::array<double, 3> ax0{invsqrt2, -invsqrt2, 0.0}, ax1{invsqrt2, invsqrt2, 0.0};

  // Find the bounding points from knowing its at 45.
  auto inf = std::numeric_limits<double>::infinity();
  double min_x = inf, min_y = inf, max_x = -inf, max_y = -inf;
  for (const auto &v : vertices)
  {
    min_x = min_x > v(0) ? v(0) : min_x;
    min_y = min_y > v(1) ? v(1) : min_y;
    max_x = max_x < v(0) ? v(0) : max_x;
    max_y = max_y < v(1) ? v(1) : max_y;
  }
  auto length_x = (max_x - min_x) * invsqrt2;
  auto length_y = (max_y - min_y) * invsqrt2;

  CHECK(length_x == Approx(length_y).margin(1e-6));

  auto length = (length_x + length_y) / 2;
  auto lengths = box.Lengths();
  CHECK(lengths[0] == Approx(length).margin(1e-6));
  CHECK(lengths[1] == Approx(length).margin(1e-6));
  CHECK(lengths[2] == Approx(0.0));

  auto normals = box.Normals();
  CHECK(normals[0][0] == Approx(ax0[0]).margin(1e-4));
  CHECK(normals[0][1] == Approx(ax0[1]).margin(1e-4));
  CHECK(normals[0][2] == Approx(ax0[2]).margin(1e-4));
  CHECK(normals[1][0] == Approx(ax1[0]).margin(1e-4));
  CHECK(normals[1][1] == Approx(ax1[1]).margin(1e-4));
  CHECK(normals[1][2] == Approx(ax1[2]).margin(1e-4));
  CHECK(box.planar);
}

TEST_CASE("TetToHex", "[geodata]")
{
  // Pull from the mfem externals data folder.
  auto ref_tet_path = fs::path(MFEM_DATA_PATH) / "ref-tetrahedron.mesh";
  mfem::Mesh single_tet(ref_tet_path.string());

  SECTION("Linear")
  {
    int order = 1;
    single_tet.EnsureNodes();
    single_tet.SetCurvature(order);
    auto four_hex = mesh::MeshTetToHex(single_tet);
    CHECK(four_hex.GetNE() == 4);

    // DOFs are added in vert -> edge -> face order, based on local ordering in vertex,
    // which has dofs (3,2,1,0).
    const std::vector<std::array<double, 3>> global_dof_vals{
        {0.0, 0.0, 0.0},              // in elem 3
        {1.0, 0.0, 0.0},              // in elem 2
        {0.0, 1.0, 0.0},              // in elem 1
        {0.0, 0.0, 1.0},              // in elem 0
        {0.0, 0.5, 0.5},              // 3 -> 2
        {0.5, 0.0, 0.5},              // 3 -> 1
        {0.0, 0.0, 0.5},              // 3 -> 0
        {0.5, 0.5, 0.0},              // 2 -> 1
        {0.0, 0.5, 0.0},              // 2 -> 0
        {0.5, 0.0, 0.0},              // 1 -> 0
        {1.0 / 3, 1.0 / 3, 0.0},      // opp 3
        {1.0 / 3, 0.0, 1.0 / 3},      // opp 2
        {0.0, 1.0 / 3, 1.0 / 3},      // opp 1
        {1.0 / 3, 1.0 / 3, 1.0 / 3},  // opp 0
        {0.25, 0.25, 0.25}};
    // From drawing out dof diagram.
    const std::vector<std::array<int, 8>> elem_dofs{{3, 4, 13, 5, 6, 12, 14, 11},
                                                    {2, 7, 13, 4, 8, 10, 14, 12},
                                                    {1, 5, 13, 7, 9, 11, 14, 10},
                                                    {0, 9, 10, 8, 6, 11, 14, 12}};

    for (int i = 0; i < 14; i++)
      for (int j = 0; j < 3; j++)
      {
        // margin(1e-12) for comparing zeros.
        CHECK((*four_hex.GetNodes())(j + 3 * i) ==
              Approx(global_dof_vals[i][j]).margin(1e-12));
      }

    mfem::Vector vdof_vals, col;
    four_hex.GetNodes()->GetElementDofValues(0, vdof_vals);
    REQUIRE(four_hex.GetNodes()->FESpace()->GetOrdering() == mfem::Ordering::byVDIM);
    REQUIRE(vdof_vals.Size() == 3 * 8);
    mfem::DenseMatrix vdof_vals_mat(vdof_vals.GetData(), 8, 3);

    auto check_mat =
        [&col, &global_dof_vals, &vdof_vals_mat](const std::array<int, 8> &verts)
    {
      for (std::size_t i = 0; i < 3; i++)
      {
        vdof_vals_mat.GetColumn(i, col);
        for (int j = 0; j < col.Size(); j++)
        {
          CAPTURE(i, j, global_dof_vals[verts[j]][i], col(j));
          // margin(1e-12) for comparing zeros.
          CHECK(col(j) == Approx(global_dof_vals[verts[j]][i]).margin(1e-12));
        }
      }
    };

    for (int i = 0; i < 4; i++)
    {
      four_hex.GetNodes()->GetElementDofValues(i, vdof_vals);
      check_mat(elem_dofs[i]);
    }
  }
#if defined(PALACE_WITH_GSLIB)
  SECTION("UniformSampler")
  {
    // Use GSLIB to find a mapping for the single tet mesh, with some randomly perturbed
    // position data.
    int order = GENERATE(2, 3);
    const int sdim = 3;
    // Create linear meshes, to copy the nodes to.
    mfem::Mesh linear_single_tet(single_tet);
    auto linear_four_hex = mesh::MeshTetToHex(linear_single_tet);
    linear_single_tet.EnsureNodes();
    linear_four_hex.EnsureNodes();
    REQUIRE(linear_single_tet.GetNodes());
    REQUIRE(linear_single_tet.GetNodes()->FESpace()->GetMaxElementOrder() == 1);
    REQUIRE(linear_four_hex.GetNodes());
    REQUIRE(linear_four_hex.GetNodes()->FESpace()->GetMaxElementOrder() == 1);

    single_tet.EnsureNodes();
    single_tet.SetCurvature(order);

    // Randomly perturb the non-vertex data with positive values (ensures non-zeros in later
    // comparison).
    for (int i = 0; i < single_tet.GetNodes()->Size(); i++)
    {
      (*single_tet.GetNodes())(i) += 0.05 * (1.0 + (double)rand() / RAND_MAX);
    }

    auto four_hex = mesh::MeshTetToHex(single_tet);
    REQUIRE(four_hex.GetNE() == 4);

    // Helper to generate a set of locations to sample tet at.
    // n_sample is number of points along an edge, thus >= 2.
    // In byVDIM ordering [x1,y1,z1,x2,y2,z2,...]
    auto gen_samples = [](int n_sample)
    {
      REQUIRE(n_sample >= 2);
      mfem::Vector xyz_samples(3 * (n_sample * (n_sample + 1) * (n_sample + 2) / 6));
      int o = 0;
      for (double k = 0; k < n_sample; k++)
        for (double j = 0; j < n_sample - k; j++)
          for (double i = 0; i < n_sample - j - k; i++)
          {
            xyz_samples(o) = i / (n_sample - 1);
            xyz_samples(o + 1) = j / (n_sample - 1);
            xyz_samples(o + 2) = k / (n_sample - 1);
            o++;
          }
      return xyz_samples;
    };

    // Uniform sampling over the tet as "xyz" coords of the original reference tet.
    auto xyz_samples = gen_samples(order + 2);

    // Create FiniteElementSpace on the linear meshes, with the same dofs from the higher
    // mesh nodes, then sample the node functions using coordinates from the linear meshes.
    // These should be equal to each other, as the sample points correspond to the original
    // reference space on the tet.
    const auto &tet_FESpace = single_tet.GetNodes()->FESpace();
    const auto &hex_FESpace = four_hex.GetNodes()->FESpace();
    mfem::FiniteElementSpace linear_tet_FESpace(&linear_single_tet, tet_FESpace->FEColl(),
                                                sdim, tet_FESpace->GetOrdering());
    mfem::FiniteElementSpace linear_hex_FESpace(&linear_four_hex, hex_FESpace->FEColl(),
                                                sdim, hex_FESpace->GetOrdering());
    mfem::GridFunction tet_nodes_on_linear_tet(&linear_tet_FESpace);
    mfem::GridFunction hex_nodes_on_linear_hex(&linear_hex_FESpace);
    REQUIRE(tet_nodes_on_linear_tet.Size() == single_tet.GetNodes()->Size());
    REQUIRE(hex_nodes_on_linear_hex.Size() == four_hex.GetNodes()->Size());
    tet_nodes_on_linear_tet = *single_tet.GetNodes();
    hex_nodes_on_linear_hex = *four_hex.GetNodes();

    mfem::Vector tet_vals(xyz_samples.Size()), hex_vals(xyz_samples.Size());
    fem::InterpolateFunction(xyz_samples, tet_nodes_on_linear_tet, tet_vals,
                             tet_FESpace->GetOrdering());
    fem::InterpolateFunction(xyz_samples, hex_nodes_on_linear_hex, hex_vals,
                             hex_FESpace->GetOrdering());
    for (int i = 0; i < tet_vals.Size(); i++)
    {
      CHECK(tet_vals(i) == Approx(hex_vals(i)).margin(1e-9));
    }
  }
#endif
}

}  // namespace palace