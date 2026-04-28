// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

#include <algorithm>
#include <complex>
#include <unordered_set>
#include "linalg/errorestimator.hpp"
#include "linalg/operator.hpp"
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

BoundaryModeSolver::BoundaryModeSolver(const IoData &iodata, bool root, int size,
                                       int num_thread, const char *git_tag)
  : BaseSolver(iodata, root, size, num_thread, git_tag)
{
}

double BoundaryModeSolver::PreprocessMesh(std::unique_ptr<mfem::Mesh> &smesh,
                                          MPI_Comm comm) const
{
  const auto &bm_data = iodata.solver.boundary_mode;
  if (bm_data.attributes.empty())
  {
    // Direct 2D path: caller-supplied mesh is already the solve mesh.
    return BaseSolver::PreprocessMesh(smesh, comm);
  }

  frame_ = std::make_unique<mesh::SubmeshFrame>();
  frame_->centroid.SetSize(3);
  frame_->centroid = 0.0;
  frame_->e1.SetSize(3);
  frame_->e1 = 0.0;
  frame_->e2.SetSize(3);
  frame_->e2 = 0.0;
  frame_->normal.SetSize(3);
  frame_->normal = 0.0;

  // Ranks holding the serial mesh (root + node-roots in the byte-string path) run the
  // extraction locally via mfem::SubMesh; non-holders receive the frame via broadcast.
  if (smesh)
  {
    MFEM_VERIFY(smesh->Dimension() == 3,
                "BoundaryMode with \"Attributes\" requires a 3D input mesh!");
    Mpi::Print(" Extracting 2D submesh from 3D boundary attributes...\n");

    mfem::Array<int> attr_list;
    attr_list.Append(bm_data.attributes.data(), bm_data.attributes.size());

    // Parent boundary face types whose intersection with the cross-section should
    // become 2D boundary edges: BCs that pin the tangential E-field (PEC, auxpec,
    // impedance, conductivity, farfield) plus other-waveports (relabeled to PEC below).
    std::vector<int> internal_bdr_attrs;
    const auto &bdr = iodata.boundaries;
    auto append = [&](const auto &src) {
      internal_bdr_attrs.insert(internal_bdr_attrs.end(), src.begin(), src.end());
    };
    append(bdr.pec.attributes);
    append(bdr.auxpec.attributes);
    for (const auto &d : bdr.impedance) { append(d.attributes); }
    for (const auto &d : bdr.conductivity) { append(d.attributes); }
    append(bdr.farfield.attributes);

    std::unordered_set<int> other_waveport_attrs;
    for (const auto &[idx, data] : bdr.waveport)
    {
      for (auto a : data.attributes)
      {
        if (std::find(bm_data.attributes.begin(), bm_data.attributes.end(), a) !=
            bm_data.attributes.end())
        {
          continue;
        }
        other_waveport_attrs.insert(a);
        internal_bdr_attrs.push_back(a);
      }
    }

    auto extr = mesh::ExtractBoundary2DSubmesh(*smesh, attr_list, internal_bdr_attrs);

    // Bake "other waveports act as PEC on this cross-section" into the mesh so the
    // operator sees plain PEC edges and has no frame-dependent branch. Per-waveport
    // frames compose this way without mutating iodata.
    if (!other_waveport_attrs.empty())
    {
      MFEM_VERIFY(!bdr.pec.attributes.empty(),
                  "BoundaryMode submesh extraction found other-waveport edges on the "
                  "cross-section but no PEC attribute is defined in the config to "
                  "relabel them to \u2014 define at least one PEC boundary attribute.");
      const int pec_attr = *std::min_element(bdr.pec.attributes.begin(),
                                             bdr.pec.attributes.end());
      int relabelled = 0;
      for (int sbe = 0; sbe < extr.mesh->GetNBE(); sbe++)
      {
        if (other_waveport_attrs.count(extr.mesh->GetBdrAttribute(sbe)))
        {
          extr.mesh->SetBdrAttribute(sbe, pec_attr);
          relabelled++;
        }
      }
      if (relabelled > 0)
      {
        Mpi::Print(" Relabelled {:d} other-waveport boundary edge(s) as PEC (attr {:d})\n",
                   relabelled, pec_attr);
      }
    }

    *frame_ = extr.frame;
    smesh = std::move(extr.mesh);
    Mpi::Print(" Surface normal = ({:+.3e}, {:+.3e}, {:+.3e})\n", frame_->normal(0),
               frame_->normal(1), frame_->normal(2));
  }

  // Broadcast user-units frame to non-holding ranks.
  Mpi::Broadcast(3, frame_->centroid.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame_->e1.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame_->e2.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame_->normal.HostReadWrite(), 0, comm);

  // Resolve Lc from the 2D solve mesh, then rescale the centroid into nondim units so
  // the frame is ready for Solve to consume without further rescaling (e1/e2/normal are
  // unit vectors).
  const double Lc = BaseSolver::PreprocessMesh(smesh, comm);
  MFEM_VERIFY(Lc > 0.0, "BoundaryMode: failed to resolve a positive reference length "
                        "from the extracted 2D submesh!");
  frame_->centroid /= Lc;
  return Lc;
}

std::pair<ErrorIndicator, long long int>
BoundaryModeSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &bm_data = iodata.solver.boundary_mode;
  const int num_modes = bm_data.n;
  const double freq_GHz = bm_data.freq;
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide mode analysis at f = {:.3e} GHz "
             "(omega = {:.6e})\n",
             freq_GHz, omega);

  // Construct MaterialOperator on the driver side: when the solve mesh came from a 3D
  // parent (frame_ non-null), tensors need rotating into the tangent plane and the frame
  // lives here. Injecting it into BoundaryModeOperator keeps the operator frame-agnostic.
  BlockTimer bt0(Timer::CONSTRUCT);
  auto mat_op = std::make_unique<MaterialOperator>(iodata, *mesh.back());
  if (frame_)
  {
    mat_op->RotateMaterialTensors(iodata, frame_->e1, frame_->e2, frame_->normal);
  }

  BoundaryModeOperator mode_op(iodata, mesh, std::move(mat_op));

  // Combined block system DBC tdof list: ND indices, then H1 indices offset by nd_size.
  const int nd_size = mode_op.GetNDTrueVSize();
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(mode_op.GetNDDbcTDofLists().back());
  for (int i = 0; i < mode_op.GetH1DbcTDofLists().back().Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + mode_op.GetH1DbcTDofLists().back()[i]);
  }

  const bool have_mg = (mode_op.GetNDSpaceHierarchy().GetNumLevels() > 1);
  if (have_mg)
  {
    Mpi::Print(" Using p-multigrid preconditioning with {:d} levels\n",
               static_cast<int>(mode_op.GetNDSpaceHierarchy().GetNumLevels()));
  }
  const auto which_eig = (bm_data.target > 0.0)
                             ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                             : EigenvalueSolver::WhichType::LARGEST_REAL;
  const mfem::Vector *submesh_normal = frame_ ? &frame_->normal : nullptr;
  ModeEigenSolver eig(mode_op.GetMaterialOp(), submesh_normal, mode_op.GetSurfZOp(),
                      mode_op.GetFarfieldOp(), mode_op.GetSurfSigmaOp(),
                      mode_op.GetNDSpace(), mode_op.GetH1Space(), dbc_tdof_list, num_modes,
                      bm_data.max_size, bm_data.tol, which_eig, iodata.solver.linear,
                      bm_data.type, iodata.problem.verbose, mode_op.GetComm(),
                      have_mg ? &mode_op : nullptr);

  BoundaryModeFluxErrorEstimator<ComplexVector> estimator(
      mode_op.GetMaterialOp(), mode_op.GetNDSpaceHierarchy(),
      mode_op.GetRTSpaceHierarchy(), mode_op.GetCurlSpace(), mode_op.GetH1SpaceHierarchy(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);

  // Project iodata 3D path coordinates onto the 2D local frame when we came from a 3D
  // parent. Frame is already in nondim units (see PreprocessMesh).
  PostOperator<ProblemType::BOUNDARYMODE> post_op(iodata, mode_op);
  if (frame_)
  {
    post_op.ProjectImpedancePaths(frame_->centroid, frame_->e1, frame_->e2);
    post_op.ProjectVoltagePaths(frame_->centroid, frame_->e1, frame_->e2);
  }

  ErrorIndicator indicator;

  // Determine kn_target.
  double kn_target;
  if (bm_data.target > 0.0)
  {
    kn_target = bm_data.target * omega;
    Mpi::Print(" Target n_eff = {:.6e}, kn_target = {:.6e}\n", bm_data.target, kn_target);
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
    eig.ApplyVDBackTransform(e0, kn, et, en);

    // Power-normalize to |P| = 1. Recompute P on the normalized mode so its phase
    // (used downstream for impedance postprocessing) reflects the final field.
    auto P = eig.ComputePoyntingPower(omega, kn, et, en);
    if (std::abs(P) > 0.0)
    {
      e0 *= 1.0 / std::sqrt(std::abs(P));
      P = eig.ComputePoyntingPower(omega, kn, et, en);
    }

    const double error_bkwd = eig.GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    const double error_abs = eig.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

    auto total_domain_energy = post_op.MeasureAndPrintAll(i, et, en, kn, P, omega,
                                                          error_abs, error_bkwd, n_print);

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
