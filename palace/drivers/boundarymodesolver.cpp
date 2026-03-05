// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodesolver.hpp"

#include <algorithm>
#include <complex>
#include <vector>

#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "fem/multigrid.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/postoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

std::pair<ErrorIndicator, long long int>
BoundaryModeSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &bm_data = iodata.solver.boundary_mode;
  const bool use_submesh = !bm_data.attributes.empty();

  if (use_submesh)
  {
    MFEM_VERIFY(mesh.back()->Dimension() == 3,
                "BoundaryMode with \"Attributes\" requires a 3D mesh!");
  }
  else
  {
    MFEM_VERIFY(mesh.back()->Dimension() == 2,
                "BoundaryMode solver requires a 2D mesh (waveguide cross-section), "
                "or a 3D mesh with \"Attributes\" specifying the cross-section boundary!");
  }

  const double freq_GHz = bm_data.freq;
  const int num_modes = bm_data.n;
  const double tol = bm_data.tol;
  const double omega =
      2.0 * M_PI * iodata.units.Nondimensionalize<Units::ValueType::FREQUENCY>(freq_GHz);

  Mpi::Print("\nConfiguring 2D waveguide mode analysis at f = {:.3e} GHz "
             "(omega = {:.6e})\n",
             freq_GHz, omega);

  BlockTimer bt0(Timer::CONSTRUCT);

  std::unique_ptr<Mesh> submesh_holder;
  std::unique_ptr<MaterialOperator> owned_mat_op;
  mfem::Vector submesh_centroid, submesh_e1, submesh_e2;
  Mesh *solve_mesh;

  if (use_submesh)
  {
    Mpi::Print(" Extracting 2D submesh from 3D boundary attributes...\n");
    const auto &parent_mesh = mesh.back()->Get();
    MPI_Comm comm = parent_mesh.GetComm();
    mfem::Array<int> attr_list;
    attr_list.Append(bm_data.attributes.data(), bm_data.attributes.size());

    // Collect all boundary attributes that need internal boundary elements on the
    // cross-section. This includes PEC, AuxPEC, impedance, conductivity, and absorbing
    // BCs — any 3D boundary face whose edges should become boundary elements in the 2D
    // mesh.
    std::vector<int> internal_bdr_attrs;
    for (auto a : iodata.boundaries.pec.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }
    for (auto a : iodata.boundaries.auxpec.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }
    for (const auto &d : iodata.boundaries.impedance)
    {
      for (auto a : d.attributes)
      {
        internal_bdr_attrs.push_back(a);
      }
    }
    for (const auto &d : iodata.boundaries.conductivity)
    {
      for (auto a : d.attributes)
      {
        internal_bdr_attrs.push_back(a);
      }
    }
    for (auto a : iodata.boundaries.farfield.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }

    // Extract a standalone 2D serial mesh from the 3D boundary, with proper domain
    // attributes, boundary attributes, internal boundary edges, and 2D projection.
    mfem::Vector surface_normal;
    auto serial_mesh = mesh::ExtractStandalone2DSubmesh(
        parent_mesh, attr_list, internal_bdr_attrs, surface_normal, submesh_centroid,
        submesh_e1, submesh_e2);

    // Repartition across all MPI ranks using the MeshPartitioner-based distribution
    // pipeline, which correctly handles shared entity topology and edge orientations.
    if (!Mpi::Root(comm))
    {
      serial_mesh.reset();
    }
    auto mesh_2d = mesh::DistributeSerialMesh(comm, serial_mesh);

    submesh_holder = std::make_unique<Mesh>(std::move(mesh_2d));
    solve_mesh = submesh_holder.get();

    // Construct MaterialOperator from the standalone 2D mesh.
    owned_mat_op = std::make_unique<MaterialOperator>(iodata, *solve_mesh);

    // Rotate material tensors from global frame to the local tangent frame.
    owned_mat_op->RotateMaterialTensors(iodata, submesh_e1, submesh_e2, surface_normal);
  }
  else
  {
    solve_mesh = mesh.back().get();
    owned_mat_op = std::make_unique<MaterialOperator>(iodata, *mesh.back());
  }
  MaterialOperator &mat_op = *owned_mat_op;

  // Collect Dirichlet boundary attributes (PEC + AuxPEC, plus other wave port boundaries
  // on the cross-section when using submesh extraction).
  mfem::Array<int> dbc_bcs;
  {
    const auto &pmesh = solve_mesh->Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    for (auto attr : iodata.boundaries.pec.attributes)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        dbc_bcs.Append(attr);
      }
    }
    for (auto attr : iodata.boundaries.auxpec.attributes)
    {
      if (attr > 0 && attr <= bdr_attr_max)
      {
        dbc_bcs.Append(attr);
      }
    }
    if (use_submesh)
    {
      for (const auto &[idx, data] : iodata.boundaries.waveport)
      {
        for (auto attr : data.attributes)
        {
          if (std::find(bm_data.attributes.begin(), bm_data.attributes.end(), attr) !=
              bm_data.attributes.end())
          {
            continue;
          }
          if (attr > 0 && attr <= bdr_attr_max)
          {
            dbc_bcs.Append(attr);
          }
        }
      }
    }
    dbc_bcs.Sort();
    dbc_bcs.Unique();
  }

  // Construct FE space hierarchies for p-multigrid. For solver order p, this builds spaces
  // at p, p-1, ..., 1 (linear coarsening) or logarithmic coarsening. When mg_max_levels <=
  // 1, each hierarchy has a single level (no multigrid — falls back to sparse direct).
  // Wrap the solve mesh in a single-element vector for
  // ConstructFiniteElementSpaceHierarchy. For the submesh case, temporarily move ownership;
  // for the non-submesh case, use a non-owning raw pointer wrapper (released after
  // hierarchy construction).
  std::vector<std::unique_ptr<Mesh>> solve_mesh_vec;
  if (use_submesh)
  {
    solve_mesh_vec.push_back(std::move(submesh_holder));
  }
  else
  {
    solve_mesh_vec.emplace_back(solve_mesh);  // Non-owning
  }

  const auto &mg = iodata.solver.linear;
  auto nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(
      iodata.solver.order, solve_mesh->Dimension(), mg.mg_max_levels, mg.mg_coarsening,
      false);
  auto h1_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      iodata.solver.order, solve_mesh->Dimension(), mg.mg_max_levels, mg.mg_coarsening,
      false);

  std::vector<mfem::Array<int>> nd_dbc_tdof_lists, h1_dbc_tdof_lists;
  auto nd_fespaces = fem::ConstructFiniteElementSpaceHierarchy(
      mg.mg_max_levels, solve_mesh_vec, nd_fecs, &dbc_bcs, &nd_dbc_tdof_lists);
  auto h1_fespaces = fem::ConstructFiniteElementSpaceHierarchy(
      mg.mg_max_levels, solve_mesh_vec, h1_fecs, &dbc_bcs, &h1_dbc_tdof_lists);

  // H1 auxiliary space hierarchy for Hiptmair distributive relaxation in the ND block.
  auto h1_aux_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      iodata.solver.order, solve_mesh->Dimension(), mg.mg_max_levels, mg.mg_coarsening,
      false);
  std::vector<mfem::Array<int>> h1_aux_dbc_tdof_lists;
  auto h1_aux_fespaces = fem::ConstructFiniteElementSpaceHierarchy(
      mg.mg_max_levels, solve_mesh_vec, h1_aux_fecs, &dbc_bcs, &h1_aux_dbc_tdof_lists);

  // Restore mesh ownership.
  if (use_submesh)
  {
    submesh_holder = std::move(solve_mesh_vec[0]);
  }
  else
  {
    solve_mesh_vec[0].release();  // Release the non-owning wrapper.
  }

  // Get finest-level spaces (used by the rest of the code).
  auto &nd_fespace = nd_fespaces.GetFinestFESpace();
  auto &h1_fespace = h1_fespaces.GetFinestFESpace();

  // Get finest-level DBC tdof lists.
  const auto &nd_dbc_tdof_list = nd_dbc_tdof_lists.back();
  const auto &h1_dbc_tdof_list = h1_dbc_tdof_lists.back();

  const int nd_size = nd_fespace.GetTrueVSize();
  const int h1_size = h1_fespace.GetTrueVSize();
  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             nd_fespace.GlobalTrueVSize(), h1_fespace.GlobalTrueVSize(),
             nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize());

  // Material operator and PostOperator.
  BoundaryModeFemOp mode_op(mat_op, nd_fespace, h1_fespace, *solve_mesh,
                            iodata.solver.order);
  PostOperator<ProblemType::BOUNDARYMODE> post_op(iodata, mode_op);

  // Project impedance/voltage path coordinates from 3D to the 2D local frame.
  if (use_submesh)
  {
    post_op.ProjectImpedancePaths(submesh_centroid, submesh_e1, submesh_e2);
    post_op.ProjectVoltagePaths(submesh_centroid, submesh_e1, submesh_e2);
  }

  // Error estimator setup. These are single-level hierarchies (no p-multigrid) wrapping the
  // finest-level FE spaces for the error estimation flux recovery.
  auto rt_fec = std::make_unique<mfem::RT_FECollection>(iodata.solver.order - 1,
                                                        solve_mesh->Dimension());
  FiniteElementSpaceHierarchy nd_fespaces_est(
      std::make_unique<FiniteElementSpace>(*solve_mesh, nd_fecs.back().get()));
  FiniteElementSpaceHierarchy rt_fespaces_est(
      std::make_unique<FiniteElementSpace>(*solve_mesh, rt_fec.get()));
  FiniteElementSpaceHierarchy h1_fespaces_est(
      std::make_unique<FiniteElementSpace>(*solve_mesh, h1_fecs.back().get()));
  TimeDependentFluxErrorEstimator<ComplexVector> estimator(
      mat_op, nd_fespaces_est, rt_fespaces_est, iodata.solver.linear.estimator_tol,
      iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg,
      &mode_op.GetCurlSpace(), &h1_fespaces_est);
  ErrorIndicator indicator;

  // Determine kn_target from user-specified n_eff target or material properties.
  double kn_target;
  if (bm_data.target > 0.0)
  {
    kn_target = bm_data.target * omega;
    Mpi::Print(" Target n_eff = {:.6e}, kn_target = {:.6e}\n", bm_data.target, kn_target);
  }
  else
  {
    double c_min = mat_op.GetLightSpeedMax().Min();
    Mpi::GlobalMin(1, &c_min, nd_fespace.GetComm());
    MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
                "Invalid material speed of light!");
    kn_target = omega / c_min * std::sqrt(1.1);
    Mpi::Print(" Auto kn_target = {:.6e} (from c_min = {:.6e})\n", kn_target, c_min);
  }

  // Construct boundary operators for impedance, absorbing, and conductivity BCs.
  SurfaceImpedanceOperator surf_z_op(iodata, mat_op, *solve_mesh);
  FarfieldBoundaryOperator farfield_op(iodata, mat_op, *solve_mesh);
  SurfaceConductivityOperator surf_sigma_op(iodata, mat_op, *solve_mesh);

  // Configure solver.
  BoundaryModeOperatorConfig config;
  config.attr_to_material = &mat_op.GetAttributeToMaterial();
  config.inv_permeability = &mat_op.GetInvPermeability();
  config.curlcurl_inv_permeability = &mat_op.GetCurlCurlInvPermeability();
  config.permittivity_real = &mat_op.GetPermittivityReal();
  config.permittivity_scalar = &mat_op.GetPermittivityScalar();
  config.normal = nullptr;  // 2D domain mesh (projected for submesh case)
  config.permittivity_imag =
      mat_op.HasLossTangent() ? &mat_op.GetPermittivityImag() : nullptr;
  config.permittivity_imag_scalar =
      mat_op.HasLossTangent() ? &mat_op.GetPermittivityImagScalar() : nullptr;
  config.has_loss_tangent = mat_op.HasLossTangent();
  config.conductivity = mat_op.HasConductivity() ? &mat_op.GetConductivity() : nullptr;
  config.has_conductivity = mat_op.HasConductivity();
  config.inv_london_depth = mat_op.HasLondonDepth() ? &mat_op.GetInvLondonDepth() : nullptr;
  config.inv_london_depth_scalar =
      mat_op.HasLondonDepth() ? &mat_op.GetInvLondonDepthScalar() : nullptr;
  config.has_london_depth = mat_op.HasLondonDepth();
  config.mat_op = &mat_op;
  config.surf_z_op = &surf_z_op;
  config.farfield_op = &farfield_op;
  config.surf_sigma_op = &surf_sigma_op;
  config.num_modes = num_modes;
  config.num_vec = bm_data.max_size;
  config.eig_tol = tol;
  config.which_eig = (bm_data.target > 0.0) ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                                            : EigenvalueSolver::WhichType::LARGEST_REAL;
  config.linear = &iodata.solver.linear;
  config.eigen_backend = bm_data.type;
  config.verbose = iodata.problem.verbose;

  // Build combined dbc_tdof_list for the block system (PEC only on both ND and H1).
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  // Construct the boundary mode solver. When the FE space hierarchy has multiple levels,
  // configure p-multigrid preconditioning with block-diagonal GMG (ND + H1).
  std::unique_ptr<BoundaryModeMultigridConfig> mg_config_ptr;
  if (nd_fespaces.GetNumLevels() > 1)
  {
    mg_config_ptr = std::make_unique<BoundaryModeMultigridConfig>();
    mg_config_ptr->nd_fespaces = &nd_fespaces;
    mg_config_ptr->h1_fespaces = &h1_fespaces;
    mg_config_ptr->h1_aux_fespaces = &h1_aux_fespaces;
    mg_config_ptr->nd_dbc_tdof_lists = &nd_dbc_tdof_lists;
    mg_config_ptr->h1_dbc_tdof_lists = &h1_dbc_tdof_lists;
    mg_config_ptr->h1_aux_dbc_tdof_lists = &h1_aux_dbc_tdof_lists;
    for (std::size_t l = 0; l < nd_fespaces.GetNumLevels(); l++)
    {
      int nd_dbc = nd_dbc_tdof_lists[l].Size();
      int h1_dbc = h1_dbc_tdof_lists[l].Size();
      Mpi::Print(" MG level {:d}: ND DOFs={:d} (DBC={:d}), H1 DOFs={:d} (DBC={:d})\n",
                 static_cast<int>(l), nd_fespaces.GetFESpaceAtLevel(l).GlobalTrueVSize(),
                 nd_dbc, h1_fespaces.GetFESpaceAtLevel(l).GlobalTrueVSize(), h1_dbc);
    }
    Mpi::Print(" Using p-multigrid preconditioning with {:d} levels\n",
               static_cast<int>(nd_fespaces.GetNumLevels()));
  }
  BoundaryModeOperator mode_solver(config, nd_fespace, h1_fespace, dbc_tdof_list,
                                   solve_mesh->GetComm(), mg_config_ptr.get());

  // Store Btt for impedance postprocessing.
  mode_op.SetBttMatrix(std::make_unique<mfem::HypreParMatrix>(*mode_solver.GetBtt()));

  // Solve the GEP. The eigenvalue is lambda = -kn^2 (shift-and-invert with
  // sigma = -kn_target^2).
  BlockTimer bt1(Timer::EPS);
  Mpi::Print("\nSolving GEP for {:d} propagation mode(s)...\n", num_modes);

  double sigma = -kn_target * kn_target;
  auto result = mode_solver.Solve(omega, sigma);
  int num_conv = result.num_converged;
  Mpi::Print(" Found {:d} converged eigenvalue{} (sigma = {:.6e})\n", num_conv,
             (num_conv != 1) ? "s" : "", sigma);

  // Print ALL converged eigenvalues. The shift-and-invert eigenvalue lambda is related
  // to the original eigenvalue -kn^2 by: kn^2 = -sigma - 1/lambda.
  for (int i = 0; i < num_conv; i++)
  {
    std::complex<double> lambda = mode_solver.GetEigenvalue(i);
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);
    Mpi::Print(" eig {:d}: kn = {:.6e}{:+.6e}i, n_eff = {:.6e}{:+.6e}i\n", i, kn.real(),
               kn.imag(), kn.real() / omega, kn.imag() / omega);
  }

  // Postprocessing.
  BlockTimer bt2(Timer::POSTPRO);
  if (const auto *ksp = mode_solver.GetLinearSolver())
  {
    SaveMetadata(*ksp);
  }
  Mpi::Print("\nComputing solution error estimates and performing postprocessing\n\n");

  auto &l2_curl_fespace = mode_op.GetCurlSpace();
  const auto &CurlOp = l2_curl_fespace.GetDiscreteInterpolator(nd_fespace);
  const int l2_size = l2_curl_fespace.GetTrueVSize();

  const int n_print = std::min(num_conv, num_modes);
  for (int i = 0; i < n_print; i++)
  {
    // Recover kn from eigenvalue (shift-and-invert recovery).
    std::complex<double> lambda = mode_solver.GetEigenvalue(i);
    std::complex<double> kn = std::sqrt(-sigma - 1.0 / lambda);
    double error_bkwd = mode_solver.GetError(i, EigenvalueSolver::ErrorType::BACKWARD);
    double error_abs = mode_solver.GetError(i, EigenvalueSolver::ErrorType::ABSOLUTE);

    // Extract et (ND) and en_tilde (H1) as views into the block eigenvector.
    ComplexVector e0(nd_size + h1_size);
    e0.UseDevice(true);
    mode_solver.GetEigenvector(i, e0);
    ComplexVector et, en;
    et.Real().MakeRef(e0.Real(), 0, nd_size);
    et.Imag().MakeRef(e0.Imag(), 0, nd_size);
    en.Real().MakeRef(e0.Real(), nd_size, h1_size);
    en.Imag().MakeRef(e0.Imag(), nd_size, h1_size);

    // Power-normalize eigenvector.
    {
      const auto *Btt = mode_op.GetBtt();
      if (Btt)
      {
        Vector Btt_etr(nd_size), Btt_eti(nd_size);
        Btt_etr.UseDevice(true);
        Btt_eti.UseDevice(true);
        Btt->Mult(et.Real(), Btt_etr);
        Btt->Mult(et.Imag(), Btt_eti);
        double p_rr = mfem::InnerProduct(nd_fespace.GetComm(), et.Real(), Btt_etr);
        double p_ii = mfem::InnerProduct(nd_fespace.GetComm(), et.Imag(), Btt_eti);
        double p_ri = mfem::InnerProduct(nd_fespace.GetComm(), et.Real(), Btt_eti);
        double p_ir = mfem::InnerProduct(nd_fespace.GetComm(), et.Imag(), Btt_etr);
        std::complex<double> etH_Btt_et(p_rr + p_ii, p_ri - p_ir);
        std::complex<double> P = 0.5 * std::conj(kn) / omega * etH_Btt_et;
        double P_abs = std::abs(P);
        if (P_abs > 0.0)
        {
          double scale = 1.0 / std::sqrt(P_abs);
          et.Real() *= scale;
          et.Imag() *= scale;
          en.Real() *= scale;
          en.Imag() *= scale;
        }
        else
        {
          double norm2 = p_rr + p_ii;
          if (norm2 > 0.0)
          {
            double scale = 1.0 / std::sqrt(norm2);
            et.Real() *= scale;
            et.Imag() *= scale;
            en.Real() *= scale;
            en.Imag() *= scale;
          }
        }
      }
    }

    auto total_domain_energy =
        post_op.MeasureAndPrintAll(i, et, en, kn, omega, error_abs, error_bkwd, n_print);

    const bool is_propagating =
        std::abs(kn.imag()) < 0.1 * std::abs(kn.real()) && std::abs(kn.real()) > 0.0;
    if (i < num_modes && is_propagating)
    {
      ComplexVector bz(l2_size);
      bz.UseDevice(true);
      {
        Vector curl_etr(l2_size), curl_eti(l2_size);
        curl_etr.UseDevice(true);
        curl_eti.UseDevice(true);
        CurlOp.Mult(et.Real(), curl_etr);
        CurlOp.Mult(et.Imag(), curl_eti);
        bz.Real() = curl_eti;
        bz.Real() *= 1.0 / omega;
        bz.Imag() = curl_etr;
        bz.Imag() *= -1.0 / omega;
      }
      estimator.AddErrorIndicator(et, bz, total_domain_energy, indicator);
    }
  }
  Mpi::Print("\n");

  post_op.MeasureFinalize(indicator);
  return {indicator, nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize()};
}

}  // namespace palace
