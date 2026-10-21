// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <numbers>
#include <unordered_set>
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/modeeigensolver.hpp"
#include "models/postoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

namespace
{

// Tangent frame of the extracted 2D submesh in the parent 3D coordinates. Used by
// Preprocess to rotate material tensors and project path coordinates into the local
// frame.
struct SubmeshFrame
{
  mfem::Vector centroid{0.0, 0.0, 0.0};
  mfem::Vector e1{0.0, 0.0, 0.0};
  mfem::Vector e2{0.0, 0.0, 0.0};
  mfem::Vector normal{0.0, 0.0, 0.0};
};

// Full 3D-boundary → 2D-submesh pipeline on the pre-partitioned serial mesh:
// CreateFromBoundary → attribute remap → add internal-boundary edges → 3D→2D
// projection. Runs on ranks that hold the serial mesh; the returned 2D mesh has its
// own node coordinates and no longer references the parent.
std::unique_ptr<mfem::Mesh>
ExtractBoundary2DSubmesh(mfem::Mesh &parent, const mfem::Array<int> &surface_attrs,
                         const std::vector<int> &internal_bdr_attrs, SubmeshFrame &frame)
{
  auto sub = std::make_unique<mfem::SubMesh>(
      mfem::SubMesh::CreateFromBoundary(parent, surface_attrs));
  // Project to 2D before remapping element attributes: ProjectSubmeshTo2D's call to
  // GetSurfaceNormal walks the unique-attribute list, which the remap rewrites from the
  // surface tags to the neighboring domain materials.
  frame.normal = mesh::ProjectSubmeshTo2D(*sub, frame.centroid, frame.e1, frame.e2);
  mesh::RemapSubMeshAttributes(*sub);
  mesh::RemapSubMeshBdrAttributes(*sub, surface_attrs);
  mesh::AddSubMeshInternalBoundaryElements(*sub, surface_attrs, internal_bdr_attrs);
  return sub;
}

// Area of an element or face, integrated exactly for the polynomial transformations of
// elements in a plane, so that the faces of the 3D mesh and the elements of the 2D
// cross-section (with the same geometry) give the same value.
double GetArea(mfem::ElementTransformation &T)
{
  const auto &ir = mfem::IntRules.Get(T.GetGeometryType(), 2 * T.OrderJ());
  double area = 0.0;
  for (int q = 0; q < ir.GetNPoints(); q++)
  {
    const auto &ip = ir.IntPoint(q);
    T.SetIntPoint(&ip);
    area += ip.weight * T.Weight();
  }
  return area;
}

// Center and area of a face of the cross-section.
struct CrossSectionFace
{
  mfem::Vector x;
  double area;
};

// The faces of the cross-section of a 3D mesh: its boundary elements with the given
// attributes.
std::vector<CrossSectionFace> GetCrossSectionFaces(const mfem::Mesh &mesh,
                                                   const mfem::Array<int> &surface_attrs)
{
  std::vector<CrossSectionFace> faces;
  mfem::IsoparametricTransformation T;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (surface_attrs.Find(mesh.GetBdrAttribute(be)) < 0)
    {
      continue;
    }
    mesh.GetBdrElementTransformation(be, &T);
    auto &face = faces.emplace_back();
    T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), face.x);
    face.area = GetArea(T);
  }
  return faces;
}

// Derefine a nonconforming 3D mesh until it has no hanging entities, one level of its
// refinement hierarchy at a time as in mesh::ReconstructCrackSides, and return the number
// of levels derefined. The mesh is not modified if its refinements cannot be derefined
// (pyramids, anisotropic refinements), and only partially derefined if the coarsest mesh
// of the hierarchy has hanging entities.
int DerefineToConforming(mfem::Mesh &mesh)
{
  auto Hanging = [&mesh]()
  {
    auto &ncmesh = *mesh.ncmesh;
    return ncmesh.GetEdgeList().slaves.Size() > 0 || ncmesh.GetFaceList().slaves.Size() > 0;
  };
  for (int e = 0; e < mesh.GetNE(); e++)
  {
    const int depth = mesh.ncmesh->GetElementDepth(e);
    if (mesh.GetElementGeometry(e) == mfem::Geometry::PYRAMID || depth > 10 ||
        mesh.ncmesh->GetElementSizeReduction(e) != (1 << (3 * depth)))
    {
      return 0;
    }
  }
  int levels = 0;
  while (Hanging())
  {
    mfem::Vector zero(mesh.GetNE());
    zero = 0.0;
    if (!mesh.DerefineByError(zero, 1.0))
    {
      break;
    }
    levels++;
  }
  return levels;
}

// Refine a 2D mesh, in at most the given number of levels, until its elements are the
// given faces (with centers in the plane of the mesh), which must refine its elements
// isotropically as the faces of the boundary of an isotropically refined 3D mesh do.
// Returns false if the faces are not recovered, for example on a cross-section which is not
// planar.
bool RefineToFaces(mfem::Mesh &mesh, const std::vector<CrossSectionFace> &faces,
                   int max_levels)
{
  if (faces.empty())
  {
    return false;
  }
  mesh.EnsureNCMesh(true);

  // Look the faces up by the cell of their center in a grid with the size of the smallest
  // face. Faces match an element with the same center and area (a refined triangle has the
  // center of its middle child).
  constexpr double tol = 1.0e-6;
  double h = mfem::infinity();
  for (const auto &face : faces)
  {
    h = std::min(h, std::sqrt(face.area));
  }
  auto Cell = [h](const mfem::Vector &x)
  { return std::array<long long, 2>{std::llround(x(0) / h), std::llround(x(1) / h)}; };
  std::map<std::array<long long, 2>, std::vector<int>> grid;
  for (int i = 0; i < static_cast<int>(faces.size()); i++)
  {
    grid[Cell(faces[i].x)].push_back(i);
  }
  auto IsFace = [&](const mfem::Vector &x, double area)
  {
    const auto cell = Cell(x);
    for (long long i = cell[0] - 1; i <= cell[0] + 1; i++)
    {
      for (long long j = cell[1] - 1; j <= cell[1] + 1; j++)
      {
        auto it = grid.find({i, j});
        if (it == grid.end())
        {
          continue;
        }
        for (int f : it->second)
        {
          if (std::abs(area - faces[f].area) < tol * faces[f].area &&
              x.DistanceTo(faces[f].x) < tol * std::sqrt(faces[f].area))
          {
            return true;
          }
        }
      }
    }
    return false;
  };

  mfem::IsoparametricTransformation T;
  mfem::Vector x;
  for (int level = 0;; level++)
  {
    mfem::Array<int> marked;
    for (int e = 0; e < mesh.GetNE(); e++)
    {
      mesh.GetElementTransformation(e, &T);
      T.Transform(mfem::Geometries.GetCenter(T.GetGeometryType()), x);
      if (!IsFace(x, GetArea(T)))
      {
        marked.Append(e);
      }
    }
    if (marked.Size() == 0)
    {
      break;
    }
    if (level == max_levels)
    {
      return false;
    }
    mesh.GeneralRefinement(marked, 1, 0);
  }
  return mesh.GetNE() == static_cast<int>(faces.size());
}

// Project each 3D point in `path` onto the 2D local frame (centroid, e1, e2). Leaves
// entries that are already 2D untouched. `path` is the iodata-side vector-of-vectors
// representation; we write back 2-element entries so the downstream PostOperator reads
// them as 2D coordinates.
void ProjectPathTo2D(std::vector<std::vector<double>> &path, const mfem::Vector &centroid,
                     const mfem::Vector &e1, const mfem::Vector &e2)
{
  for (auto &p : path)
  {
    if (p.size() != 3)
    {
      continue;
    }
    mfem::Vector p3d(p.data(), 3);
    mfem::Vector p2d = mesh::Project3Dto2D(p3d, centroid, e1, e2);
    p.assign({p2d(0), p2d(1)});
  }
}

}  // namespace

BoundaryModeSolver::BoundaryModeSolver(const IoData &iodata, bool root, int size,
                                       int num_thread, const char *git_tag)
  : BaseSolver(iodata, root, size, num_thread, git_tag)
{
}

void BoundaryModeSolver::Preprocess(IoData &iodata, std::unique_ptr<mfem::Mesh> &smesh,
                                    MPI_Comm comm) const
{
  const auto &bm = iodata.solver.boundary_mode;
  if (bm.attributes.empty())
  {
    // Direct 2D path: caller-supplied mesh is already the solve mesh.
    BaseSolver::Preprocess(iodata, smesh, comm);
    return;
  }

  SubmeshFrame frame;

  // Ranks holding the serial mesh run the extraction locally via mfem::SubMesh;
  // non-holders receive the frame via broadcast below.
  if (smesh)
  {
    MFEM_VERIFY(smesh->Dimension() == 3,
                "BoundaryMode with \"Attributes\" requires a 3D input mesh!");
    Mpi::Print(" Extracting 2D submesh from 3D boundary attributes...\n");

    mfem::Array<int> attr_list;
    attr_list.Append(bm.attributes.data(), bm.attributes.size());

    // Parent boundary face types whose intersection with the cross-section should
    // become 2D boundary edges: BCs that pin the tangential E-field (PEC, auxpec,
    // impedance, conductivity, rational impedance, farfield) plus other-waveports
    // (relabeled to PEC below).
    std::vector<int> internal_bdr_attrs;
    const auto &bdr = iodata.boundaries;
    auto append = [&](const auto &src)
    { internal_bdr_attrs.insert(internal_bdr_attrs.end(), src.begin(), src.end()); };
    append(bdr.pec.attributes);
    append(bdr.auxpec.attributes);
    for (const auto &d : bdr.impedance)
    {
      append(d.attributes);
    }
    for (const auto &d : bdr.conductivity)
    {
      append(d.attributes);
    }
    for (const auto &d : bdr.rational_impedance)
    {
      append(d.attributes);
    }
    append(bdr.farfield.attributes);

    std::unordered_set<int> other_waveport_attrs;
    for (const auto &[idx, data] : bdr.waveport)
    {
      for (auto a : data.attributes)
      {
        if (std::ranges::find(bm.attributes, a) != bm.attributes.end())
        {
          continue;
        }
        other_waveport_attrs.insert(a);
        internal_bdr_attrs.push_back(a);
      }
    }

    // Extract the cross-section of a 3D mesh, relabelling the other-waveport edges as PEC:
    // for the BoundaryMode solve they act as conducting boundaries. The boundary elements
    // of the submesh of a nonconforming mesh, including those added or relabelled here, are
    // not in its mfem::NCMesh, and would be lost when the mesh is distributed, so it is
    // rebuilt from its elements and hanging vertices.
    auto Extract = [&](mfem::Mesh &parent, SubmeshFrame &f, int &relabelled)
    {
      std::unique_ptr<mfem::Mesh> sub =
          ExtractBoundary2DSubmesh(parent, attr_list, internal_bdr_attrs, f);
      relabelled = 0;
      if (!other_waveport_attrs.empty())
      {
        MFEM_VERIFY(!bdr.pec.attributes.empty(),
                    "BoundaryMode submesh extraction found other-waveport edges on the "
                    "cross-section. Define at least one PEC boundary attribute to "
                    "relabel them to.");
        const int pec_attr = *std::ranges::min_element(bdr.pec.attributes);
        for (int sbe = 0; sbe < sub->GetNBE(); sbe++)
        {
          if (other_waveport_attrs.count(sub->GetBdrAttribute(sbe)))
          {
            sub->SetBdrAttribute(sbe, pec_attr);
            relabelled++;
          }
        }
      }
      if (sub->Nonconforming())
      {
        sub = mesh::RebuildNonconformingMesh(*sub);
      }
      return sub;
    };
    int relabelled = 0;
    auto extracted = Extract(*smesh, frame, relabelled);
    if (relabelled > 0)
    {
      Mpi::Print(" Relabelled {:d} other-waveport boundary edge(s) as PEC (attr {:d})\n",
                 relabelled, *std::ranges::min_element(bdr.pec.attributes));
    }

    // A rebuilt cross-section with hanging vertices has its leaf elements as root elements,
    // so the sides of its interior boundaries cannot be identified for error estimation.
    // Extract the cross-section instead from the coarsest conforming mesh of the refinement
    // hierarchy of the 3D mesh, and refine it into the faces of the cross-section of the
    // given mesh, as for the refined cross-section of a conforming mesh. The rebuilt
    // cross-section is kept if the 3D mesh cannot be derefined or the faces are not
    // recovered.
    if (extracted->Nonconforming())
    {
      auto faces = GetCrossSectionFaces(*smesh, attr_list);
      if (const int levels = DerefineToConforming(*smesh); levels > 0)
      {
        SubmeshFrame coarse_frame;
        int coarse_relabelled;
        auto coarse = Extract(*smesh, coarse_frame, coarse_relabelled);
        for (auto &face : faces)
        {
          face.x = mesh::Project3Dto2D(face.x, coarse_frame.centroid, coarse_frame.e1,
                                       coarse_frame.e2);
        }
        if (RefineToFaces(*coarse, faces, levels))
        {
          extracted = std::move(coarse);
          frame = coarse_frame;
        }
      }
    }

    smesh = std::move(extracted);
    Mpi::Print(" Surface normal = ({:+.3e}, {:+.3e}, {:+.3e})\n", frame.normal(0),
               frame.normal(1), frame.normal(2));
  }

  // Broadcast user-units frame to non-holding ranks.
  Mpi::Broadcast(3, frame.centroid.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame.e1.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame.e2.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame.normal.HostReadWrite(), 0, comm);

  // Bake the tangent frame into iodata (pre-nondim, user units throughout):
  //   - material tensors rotated into the local frame so MaterialOperator constructs
  //     natively-2D downstream;
  //   - impedance / voltage 3D path coordinates projected to 2D local coords so
  //     PostOperator parses them as 2D paths directly.
  RotateMaterialDefinitions(iodata.domains.materials, frame.e1, frame.e2, frame.normal);
  for (auto &[idx, imp] : iodata.boundaries.postpro.impedance)
  {
    ProjectPathTo2D(imp.voltage_path, frame.centroid, frame.e1, frame.e2);
    ProjectPathTo2D(imp.current_path, frame.centroid, frame.e1, frame.e2);
  }
  for (auto &[idx, vol] : iodata.boundaries.postpro.voltage)
  {
    ProjectPathTo2D(vol.voltage_path, frame.centroid, frame.e1, frame.e2);
  }

  // Resolve Lc from the 2D solve mesh and nondimensionalize in one pass.
  BaseSolver::Preprocess(iodata, smesh, comm);
}

std::pair<ErrorIndicator, long long int>
BoundaryModeSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &bm = iodata.solver.boundary_mode;
  const int num_modes = bm.n;
  const double freq_GHz = bm.freq;
  const double omega =
      2.0 * std::numbers::pi *
      iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide mode analysis at f = {:.3e} GHz "
             "(omega = {:.6e})\n",
             freq_GHz, omega);

  BlockTimer bt0(Timer::CONSTRUCT);
  MaterialOperator mat_op(iodata, *mesh.back());

  BoundaryModeOperator mode_op(iodata, mesh, mat_op);

  // Combined block system DBC tdof list: ND indices, then H1 indices offset by nd_size.
  const int nd_size = mode_op.GetNDTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(mode_op.GetNDDbcTDofLists().back());
  for (int i = 0; i < mode_op.GetH1DbcTDofLists().back().Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + mode_op.GetH1DbcTDofLists().back()[i]);
  }

  if (const auto n_levels = mode_op.GetNDSpaceHierarchy().GetNumLevels(); n_levels > 1)
  {
    Mpi::Print(" Using p-multigrid preconditioning with {:d} levels\n",
               static_cast<int>(n_levels));
  }
  const auto which_eig = (bm.target > 0.0) ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                                           : EigenvalueSolver::WhichType::LARGEST_REAL;
  ModeEigenSolver eig(mode_op, dbc_tdof_list, num_modes, bm.max_size, bm.tol, which_eig,
                      iodata.solver.linear, bm.type, iodata.problem.verbose);

  BoundaryModeFluxErrorEstimator<ComplexVector> estimator(
      mode_op.GetMaterialOp(), mode_op.GetNDSpaceHierarchy(), mode_op.GetRTSpaceHierarchy(),
      mode_op.GetCurlSpace(), mode_op.GetH1SpaceHierarchy(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg, iodata.boundaries.GetSheetAttributes());

  PostOperator<ProblemType::BOUNDARYMODE> post_op(iodata, mode_op);

  ErrorIndicator indicator;

  // Determine kn_target.
  double kn_target;
  if (bm.target > 0.0)
  {
    kn_target = bm.target * omega;
    Mpi::Print(" Target n_eff = {:.6e}, kn_target = {:.6e}\n", bm.target, kn_target);
  }
  else
  {
    const double mu_eps_max = mode_op.GetMaterialOp().GetMaxMuEpsilon();
    kn_target = omega * std::sqrt(1.1 * mu_eps_max);
    Mpi::Print(" Auto kn_target = {:.6e} (from max(mu_r) * max(epsilon_r) = {:.6e})\n",
               kn_target, mu_eps_max);
  }

  // Solve the eigenvalue problem.
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving GEP for {:d} propagation mode(s)...\n", num_modes);
  const double sigma = -kn_target * kn_target;
  auto result = eig.Solve(omega, sigma);
  int num_conv = result.num_converged;
  Mpi::Print(" Found {:d} converged eigenvalue{} (sigma = {:.6e})\n", num_conv,
             (num_conv != 1) ? "s" : "", result.sigma);

  for (int i = 0; i < num_conv; i++)
  {
    auto kn = eig.GetPropagationConstant(i);
    Mpi::Print(" eig {:d}: kn = {:.6e}{:+.6e}i, n_eff = {:.6e}{:+.6e}i\n", i, kn.real(),
               kn.imag(), kn.real() / omega, kn.imag() / omega);
  }

  // Postprocessing.
  BlockTimer bt2(Timer::POSTPRO);
  if (const auto *ksp = eig.GetLinearSolver())
  {
    SaveMetadata(*ksp);
  }
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n\n");

  const int h1_size = mode_op.GetH1TrueVSize();
  const int l2_size = mode_op.GetCurlSpace().GetTrueVSize();

  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    // Load eigenvector i, apply the shared VD back-transform so en holds physical En.
    ComplexVector e0(nd_size + h1_size);
    e0.UseDevice(true);
    eig.GetEigenvector(i, e0);
    ComplexVector et, en;
    const std::complex<double> kn = eig.GetPropagationConstant(i);
    mode_op.ApplyVDBackTransform(e0, kn, et, en);

    // Power-normalize the mode to |P| = 1. PostOperator re-evaluates P on the
    // normalized mode for its impedance postprocessing.
    const std::complex<double> P_initial = mode_op.ComputePoyntingPower(omega, kn, et, en);
    if (std::abs(P_initial) > 0.0)
    {
      e0 *= 1.0 / std::sqrt(std::abs(P_initial));
    }

    const double error_bkwd = eig.GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    const double error_abs = eig.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, et, en, kn, omega, error_abs, error_bkwd, n_print);

    if (i < num_modes && ModeEigenSolver::IsPropagating(kn))
    {
      // Bz = curl(Et) / (i·omega) = (Im(curl Et) - i·Re(curl Et)) / omega.
      ComplexVector bz(l2_size);
      bz.UseDevice(true);
      const auto &CurlOp =
          mode_op.GetCurlSpace().GetDiscreteInterpolator(mode_op.GetNDSpace());
      Vector curl_etr(l2_size), curl_eti(l2_size);
      curl_etr.UseDevice(true);
      curl_eti.UseDevice(true);
      CurlOp.Mult(et.Real(), curl_etr);
      CurlOp.Mult(et.Imag(), curl_eti);
      bz.Real() = curl_eti;
      bz.Real() *= 1.0 / omega;
      bz.Imag() = curl_etr;
      bz.Imag() *= -1.0 / omega;
      estimator.AddErrorIndicator(et, bz, total_domain_energy, indicator);
    }
  }
  Mpi::Print("\n");

  post_op.MeasureFinalize(indicator);

  return {indicator,
          mode_op.GetNDSpace().GlobalTrueVSize() + mode_op.GetH1Space().GlobalTrueVSize()};
}

}  // namespace palace
