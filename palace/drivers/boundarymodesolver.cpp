// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

#include <algorithm>
#include <complex>
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
  // GetSurfaceNormal walks the mesh's unique-attribute list, which RemapSubMeshAttributes
  // leaves stale relative to the per-element attributes it just rewrote.
  frame.normal = mesh::ProjectSubmeshTo2D(*sub, frame.centroid, frame.e1, frame.e2);
  mesh::RemapSubMeshAttributes(*sub);
  mesh::RemapSubMeshBdrAttributes(*sub, surface_attrs);
  mesh::AddSubMeshInternalBoundaryElements(*sub, surface_attrs, internal_bdr_attrs);
  return sub;
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
    // impedance, conductivity, farfield) plus other-waveports (relabeled to PEC below).
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
    append(bdr.farfield.attributes);

    std::unordered_set<int> other_waveport_attrs;
    for (const auto &[idx, data] : bdr.waveport)
    {
      for (auto a : data.attributes)
      {
        if (std::find(bm.attributes.begin(), bm.attributes.end(), a) != bm.attributes.end())
        {
          continue;
        }
        other_waveport_attrs.insert(a);
        internal_bdr_attrs.push_back(a);
      }
    }

    auto extracted = ExtractBoundary2DSubmesh(*smesh, attr_list, internal_bdr_attrs, frame);

    // Relabel other-waveport edges as PEC on this cross-section: for the
    // BoundaryMode solve they act as conducting boundaries.
    if (!other_waveport_attrs.empty())
    {
      MFEM_VERIFY(!bdr.pec.attributes.empty(),
                  "BoundaryMode submesh extraction found other-waveport edges on the "
                  "cross-section. Define at least one PEC boundary attribute to "
                  "relabel them to.");
      const int pec_attr =
          *std::min_element(bdr.pec.attributes.begin(), bdr.pec.attributes.end());
      int relabelled = 0;
      for (int sbe = 0; sbe < extracted->GetNBE(); sbe++)
      {
        if (other_waveport_attrs.count(extracted->GetBdrAttribute(sbe)))
        {
          extracted->SetBdrAttribute(sbe, pec_attr);
          relabelled++;
        }
      }
      if (relabelled > 0)
      {
        Mpi::Print(" Relabelled {:d} other-waveport boundary edge(s) as PEC (attr {:d})\n",
                   relabelled, pec_attr);
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
      iodata.solver.linear.estimator_mg);

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
    double c_min = mode_op.GetMaterialOp().GetLightSpeedMax().Min();
    Mpi::GlobalMin(1, &c_min, mode_op.GetComm());
    MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
                "Invalid material speed of light!");
    kn_target = omega / c_min * std::sqrt(1.1);
    Mpi::Print(" Auto kn_target = {:.6e} (from c_min = {:.6e})\n", kn_target, c_min);
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
