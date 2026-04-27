// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

#include <algorithm>
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
    // Direct 2D path: caller-supplied mesh is already the solve mesh. Let the base
    // implementation resolve Lc from it.
    return BaseSolver::PreprocessMesh(smesh, comm);
  }

  frame_ = std::make_unique<SubmeshFrame>();
  frame_->centroid.SetSize(3);
  frame_->e1.SetSize(3);
  frame_->e2.SetSize(3);
  frame_->normal.SetSize(3);
  frame_->centroid = 0.0;
  frame_->e1 = 0.0;
  frame_->e2 = 0.0;
  frame_->normal = 0.0;

  // Ranks that hold a copy of the serial mesh (root always; one-per-node roots in the
  // byte-string distribution path) perform the extraction locally. The mesh is still
  // serial at this point (Partition runs next in main.cpp) so all submesh primitives run
  // on mfem::SubMesh directly, with no MPI collectives and no ParMesh(COMM_SELF) wrapper.
  if (smesh)
  {
    MFEM_VERIFY(smesh->Dimension() == 3,
                "BoundaryMode with \"Attributes\" requires a 3D input mesh!");
    Mpi::Print(" Extracting 2D submesh from 3D boundary attributes...\n");

    mfem::Array<int> attr_list;
    attr_list.Append(bm_data.attributes.data(), bm_data.attributes.size());

    // Internal boundary attributes: parent boundary faces whose intersection with the
    // cross-section should become 2D boundary edges. This is every BC type that pins
    // the tangential E-field on the 2D eigenmode problem — PEC, auxpec, impedance,
    // conductivity, farfield (absorbing) — plus any *other* waveports so their
    // intersection edges appear as boundary elements and can be relabeled to PEC below.
    std::vector<int> internal_bdr_attrs;
    const auto &bdr = iodata.boundaries;
    internal_bdr_attrs.insert(internal_bdr_attrs.end(), bdr.pec.attributes.begin(),
                              bdr.pec.attributes.end());
    internal_bdr_attrs.insert(internal_bdr_attrs.end(), bdr.auxpec.attributes.begin(),
                              bdr.auxpec.attributes.end());
    for (const auto &d : bdr.impedance)
    {
      internal_bdr_attrs.insert(internal_bdr_attrs.end(), d.attributes.begin(),
                                d.attributes.end());
    }
    for (const auto &d : bdr.conductivity)
    {
      internal_bdr_attrs.insert(internal_bdr_attrs.end(), d.attributes.begin(),
                                d.attributes.end());
    }
    internal_bdr_attrs.insert(internal_bdr_attrs.end(), bdr.farfield.attributes.begin(),
                              bdr.farfield.attributes.end());

    // Collect other-waveport attributes: waveports that are not part of the surface
    // currently being solved. Edges where the cross-section intersects such a waveport
    // behave as PEC for this mode problem (this cross-section doesn't see the other
    // waveport as a port, only as a bounding conductor).
    std::unordered_set<int> other_waveport_attrs;
    {
      const auto &bm_attrs = bm_data.attributes;
      for (const auto &[idx, data] : bdr.waveport)
      {
        for (auto a : data.attributes)
        {
          if (std::find(bm_attrs.begin(), bm_attrs.end(), a) != bm_attrs.end())
          {
            continue;
          }
          other_waveport_attrs.insert(a);
          internal_bdr_attrs.push_back(a);
        }
      }
    }

    // Step 1: extract the 2D SubMesh from the 3D boundary. SubMesh inherits from Mesh,
    // so we hold it through a base-class pointer; downstream code never needs the
    // SubMesh-specific API once extraction is done.
    auto sub = std::make_unique<mfem::SubMesh>(
        mfem::SubMesh::CreateFromBoundary(*smesh, attr_list));

    // Step 2: remap domain attributes (from neighboring 3D elements) and boundary
    // attributes (from parent boundary faces), then add internal boundary edges at
    // surface/internal-boundary intersections. Must run while the parent is alive.
    mesh::RemapSubMeshAttributes(*sub);
    mesh::RemapSubMeshBdrAttributes(*sub, attr_list);
    mesh::AddSubMeshInternalBoundaryElements(*sub, attr_list, internal_bdr_attrs);

    // Step 3: relabel other-waveport attributes on the 2D boundary to the first PEC
    // attribute. This bakes the "other waveports act as PEC on this cross-section"
    // rule into the mesh itself, so BoundaryModeOperator sees plain PEC edges and has
    // no parent/frame-dependent branch. Per-waveport frames therefore compose without
    // mutating iodata.
    if (!other_waveport_attrs.empty())
    {
      MFEM_VERIFY(!bdr.pec.attributes.empty(),
                  "BoundaryMode submesh extraction found other-waveport edges on the "
                  "cross-section but no PEC attribute is defined in the config to relabel "
                  "them to \u2014 define at least one PEC boundary attribute.");
      const int pec_attr = *std::min_element(bdr.pec.attributes.begin(),
                                             bdr.pec.attributes.end());
      int relabelled = 0;
      for (int sbe = 0; sbe < sub->GetNBE(); sbe++)
      {
        const int a = sub->GetBdrAttribute(sbe);
        if (other_waveport_attrs.count(a))
        {
          sub->SetBdrAttribute(sbe, pec_attr);
          relabelled++;
        }
      }
      if (relabelled > 0)
      {
        Mpi::Print(" Relabelled {:d} other-waveport boundary edge(s) as PEC (attr {:d})\n",
                   relabelled, pec_attr);
      }
    }

    // Step 4: project from 3D ambient to true 2D coordinates. Captures the surface
    // normal and tangent frame (centroid, e1, e2) used by Solve for material tensor
    // rotation and 3D-to-2D path projection.
    frame_->normal = mesh::ProjectSubmeshTo2D(*sub, &frame_->centroid, &frame_->e1,
                                              &frame_->e2);
    Mpi::Print(" Surface normal = ({:+.3e}, {:+.3e}, {:+.3e})\n", frame_->normal(0),
               frame_->normal(1), frame_->normal(2));

    // Hand ownership over as a plain Mesh pointer. The SubMesh's parent_ pointer is now
    // dangling but no downstream code (nondim, Partition, Refine, FE-space construction)
    // dereferences it; only the Mesh base interface is used from here on.
    smesh = std::move(sub);
  }

  // Broadcast the frame from world root so all ranks have it for downstream material-
  // tensor rotation and path projection. The frame is captured from the pre-nondim
  // serial mesh — Solve scales the centroid into the solve-mesh frame when projecting
  // path points (iodata paths are scaled by NondimensionalizeInputs after this hook).
  Mpi::Broadcast(3, frame_->centroid.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame_->e1.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame_->e2.HostReadWrite(), 0, comm);
  Mpi::Broadcast(3, frame_->normal.HostReadWrite(), 0, comm);

  // Delegate to the base to resolve Lc from the now-reshaped 2D solve mesh (which may
  // have a bbox very different from the 3D parent's). Running Lc on the submesh is the
  // whole reason this hook runs before nondimensionalization.
  return BaseSolver::PreprocessMesh(smesh, comm);
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

  // Construct the MaterialOperator here rather than inside BoundaryModeOperator: when the
  // solve mesh was extracted from a 3D parent (frame_ non-null) the tensors have to be
  // rotated into the local tangent plane, and the frame lives on the driver. Keeping this
  // out of BoundaryModeOperator lets the operator remain frame-agnostic (no knowledge of
  // whether it came from a submesh) and makes per-waveport extractions compose cleanly
  // without mutating iodata.
  BlockTimer bt0(Timer::CONSTRUCT);
  auto mat_op = std::make_unique<MaterialOperator>(iodata, *mesh.back());
  if (frame_)
  {
    mat_op->RotateMaterialTensors(iodata, frame_->e1, frame_->e2, frame_->normal);
  }

  BoundaryModeOperator mode_op(iodata, mesh, std::move(mat_op));

  // Construct the eigenvalue solver on top of the operator's FE context.
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

  // Construct the error estimator. BoundaryModeOperator owns the RT hierarchy for
  // flux recovery (analogous to SpaceOperator::rt_fespaces).
  BoundaryModeFluxErrorEstimator<ComplexVector> estimator(
      mode_op.GetMaterialOp(), mode_op.GetNDSpaceHierarchy(),
      mode_op.GetRTSpaceHierarchy(), mode_op.GetCurlSpace(), mode_op.GetH1SpaceHierarchy(),
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);

  // Construct PostOperator. When the solve mesh was extracted from a 3D parent, project
  // iodata 3D path coordinates onto the 2D local frame. The frame's centroid was captured
  // in user units before NondimensionalizeInputs ran; scale it to match the now-nondim
  // iodata path coordinates before projection.
  PostOperator<ProblemType::BOUNDARYMODE> post_op(iodata, mode_op);
  if (frame_)
  {
    mfem::Vector centroid_nd(frame_->centroid);
    centroid_nd /= iodata.units.GetMeshLengthRelativeScale();
    post_op.ProjectImpedancePaths(centroid_nd, frame_->e1, frame_->e2);
    post_op.ProjectVoltagePaths(centroid_nd, frame_->e1, frame_->e2);
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
    ComplexVector e0(nd_size + h1_size);
    e0.UseDevice(true);
    ComplexVector et, en;
    std::complex<double> kn = eig.GetPhysicalMode(i, omega, e0, et, en);
    double error_bkwd = eig.GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = eig.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

    // Poynting power P for impedance postprocessing (|P| ≈ 1 after normalization).
    std::complex<double> P = eig.ComputePoyntingPower(omega, kn, et, en);

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
