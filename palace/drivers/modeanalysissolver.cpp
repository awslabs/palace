// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "modeanalysissolver.hpp"

#include <algorithm>
#include <complex>
#include <vector>

#include "fem/errorindicator.hpp"
#include "fem/fespace.hpp"
#include "fem/mesh.hpp"
#include "linalg/errorestimator.hpp"
#include "linalg/vector.hpp"
#include "models/boundarymodeoperator.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/modeanalysisoperator.hpp"
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
ModeAnalysisSolver::Solve(const std::vector<std::unique_ptr<Mesh>> &mesh) const
{
  const auto &ma_data = iodata.solver.mode_analysis;
  const bool use_submesh = !ma_data.attributes.empty();

  if (use_submesh)
  {
    MFEM_VERIFY(mesh.back()->Dimension() == 3,
                "ModeAnalysis with \"Attributes\" requires a 3D mesh!");
  }
  else
  {
    MFEM_VERIFY(mesh.back()->Dimension() == 2,
                "ModeAnalysis solver requires a 2D mesh (waveguide cross-section), "
                "or a 3D mesh with \"Attributes\" specifying the cross-section boundary!");
  }

  const double freq_GHz = ma_data.freq;
  const int num_modes = ma_data.n;
  const double tol = ma_data.tol;
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
    attr_list.Append(ma_data.attributes.data(), ma_data.attributes.size());

    // Collect PEC-related boundary attributes from config.
    std::vector<int> pec_bdr_attrs;
    for (auto a : iodata.boundaries.pec.attributes)
    {
      pec_bdr_attrs.push_back(a);
    }
    for (auto a : iodata.boundaries.auxpec.attributes)
    {
      pec_bdr_attrs.push_back(a);
    }
    for (const auto &d : iodata.boundaries.conductivity)
    {
      for (auto a : d.attributes)
      {
        pec_bdr_attrs.push_back(a);
      }
    }

    // Extract a standalone 2D serial mesh from the 3D boundary, with proper domain
    // attributes, boundary attributes, PEC internal edges, and 2D projection.
    mfem::Vector surface_normal;
    auto serial_mesh = mesh::ExtractStandalone2DSubmesh(
        parent_mesh, attr_list, pec_bdr_attrs, surface_normal, submesh_centroid,
        submesh_e1, submesh_e2);

    // Repartition across all MPI ranks and construct a distributed ParMesh.
    int nprocs = Mpi::Size(comm);
    auto *partitioning = serial_mesh->GeneratePartitioning(nprocs);
    auto mesh_2d = std::make_unique<mfem::ParMesh>(comm, *serial_mesh, partitioning);
    delete[] partitioning;

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

  // Construct FE spaces: ND for tangential E, H1 for normal (out-of-plane) E.
  auto nd_fec = std::make_unique<mfem::ND_FECollection>(iodata.solver.order,
                                                        solve_mesh->Dimension());
  auto h1_fec = std::make_unique<mfem::H1_FECollection>(iodata.solver.order,
                                                        solve_mesh->Dimension());
  FiniteElementSpace nd_fespace(*solve_mesh, nd_fec.get());
  FiniteElementSpace h1_fespace(*solve_mesh, h1_fec.get());

  // Essential (Dirichlet) BCs: PEC + AuxPEC + conductivity boundaries. For submesh mode
  // analysis, other wave port boundaries on the cross-section are also treated as PEC.
  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  {
    const auto &pmesh = solve_mesh->Get();
    int bdr_attr_max = pmesh.bdr_attributes.Size() ? pmesh.bdr_attributes.Max() : 0;
    mfem::Array<int> dbc_bcs;
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
    for (const auto &data : iodata.boundaries.conductivity)
    {
      for (auto attr : data.attributes)
      {
        if (attr > 0 && attr <= bdr_attr_max)
        {
          dbc_bcs.Append(attr);
        }
      }
    }
    if (use_submesh)
    {
      // Other wave port boundaries become PEC on this cross-section.
      for (const auto &[idx, data] : iodata.boundaries.waveport)
      {
        for (auto attr : data.attributes)
        {
          if (std::find(ma_data.attributes.begin(), ma_data.attributes.end(), attr) !=
              ma_data.attributes.end())
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

    auto dbc_marker = mesh::AttrToMarker(bdr_attr_max, dbc_bcs);
    nd_fespace.Get().GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
    h1_fespace.Get().GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  }

  const int nd_size = nd_fespace.GetTrueVSize();
  const int h1_size = h1_fespace.GetTrueVSize();
  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             nd_fespace.GlobalTrueVSize(), h1_fespace.GlobalTrueVSize(),
             nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize());

  // Material operator and PostOperator.
  ModeAnalysisOperator mode_op(mat_op, nd_fespace, h1_fespace, *solve_mesh,
                               iodata.solver.order);
  PostOperator<ProblemType::MODEANALYSIS> post_op(iodata, mode_op);

  // Project impedance path coordinates from 3D to the 2D local frame.
  if (use_submesh)
  {
    post_op.ProjectImpedancePaths(submesh_centroid, submesh_e1, submesh_e2);
  }

  // Error estimator setup (skip for submesh for now).
  std::unique_ptr<mfem::RT_FECollection> rt_fec;
  std::unique_ptr<FiniteElementSpaceHierarchy> nd_fespaces_h, rt_fespaces_h,
      h1_fespaces_est;
  std::unique_ptr<TimeDependentFluxErrorEstimator<ComplexVector>> estimator;
  if (!use_submesh)
  {
    rt_fec = std::make_unique<mfem::RT_FECollection>(iodata.solver.order - 1,
                                                     solve_mesh->Dimension());
    nd_fespaces_h = std::make_unique<FiniteElementSpaceHierarchy>(
        std::make_unique<FiniteElementSpace>(*solve_mesh, nd_fec.get()));
    rt_fespaces_h = std::make_unique<FiniteElementSpaceHierarchy>(
        std::make_unique<FiniteElementSpace>(*solve_mesh, rt_fec.get()));
    h1_fespaces_est = std::make_unique<FiniteElementSpaceHierarchy>(
        std::make_unique<FiniteElementSpace>(*solve_mesh, h1_fec.get()));
    estimator = std::make_unique<TimeDependentFluxErrorEstimator<ComplexVector>>(
        mat_op, *nd_fespaces_h, *rt_fespaces_h, iodata.solver.linear.estimator_tol,
        iodata.solver.linear.estimator_max_it, 0, iodata.solver.linear.estimator_mg,
        &mode_op.GetCurlSpace(), h1_fespaces_est.get());
  }
  ErrorIndicator indicator;

  // Determine kn_target from user-specified n_eff target or material properties.
  double kn_target;
  if (ma_data.target > 0.0)
  {
    kn_target = ma_data.target * omega;
    Mpi::Print(" Target n_eff = {:.6e}, kn_target = {:.6e}\n", ma_data.target, kn_target);
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
  // For submesh mode analysis, these operators can't be constructed on 1D boundary edges
  // with 2D material tensors. The BoundaryModeOperator handles impedance/absorbing BCs
  // internally through material coefficients instead.
  std::unique_ptr<SurfaceImpedanceOperator> surf_z_op;
  std::unique_ptr<FarfieldBoundaryOperator> farfield_op;
  std::unique_ptr<SurfaceConductivityOperator> surf_sigma_op;
  if (!use_submesh)
  {
    surf_z_op = std::make_unique<SurfaceImpedanceOperator>(iodata, mat_op, *solve_mesh);
    farfield_op = std::make_unique<FarfieldBoundaryOperator>(iodata, mat_op, *solve_mesh);
    surf_sigma_op =
        std::make_unique<SurfaceConductivityOperator>(iodata, mat_op, *solve_mesh);
  }

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
  config.inv_london_depth =
      mat_op.HasLondonDepth() ? &mat_op.GetInvLondonDepth() : nullptr;
  config.has_london_depth = mat_op.HasLondonDepth();
  config.mat_op = &mat_op;
  config.surf_z_op = surf_z_op.get();
  config.farfield_op = farfield_op.get();
  config.surf_sigma_op = surf_sigma_op.get();
  config.num_modes = num_modes;
  config.num_vec = -1;
  config.eig_tol = tol;
  config.which_eig = (ma_data.target > 0.0)
                         ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                         : EigenvalueSolver::WhichType::LARGEST_REAL;
  config.linear = &iodata.solver.linear;
  config.eigen_backend = ma_data.type;
  config.verbose = iodata.problem.verbose;

  // Build combined dbc_tdof_list for the block system (PEC only on both ND and H1).
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_list);
  for (int i = 0; i < h1_dbc_tdof_list.Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_list[i]);
  }

  // Construct the boundary mode solver.
  BoundaryModeOperator mode_solver(config, nd_fespace, h1_fespace, dbc_tdof_list,
                                   solve_mesh->GetComm());

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
    Mpi::Print(" eig {:d}: kn = {:.6e}{:+.6e}i, n_eff = {:.6e}{:+.6e}i\n", i,
               kn.real(), kn.imag(), kn.real() / omega, kn.imag() / omega);
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

    // Extract et (ND) and en_tilde (H1) from eigenvector.
    ComplexVector et(nd_size), en(h1_size);
    {
      ComplexVector e0(nd_size + h1_size);
      mode_solver.GetEigenvector(i, e0);
      std::copy_n(e0.Real().begin(), nd_size, et.Real().begin());
      std::copy_n(e0.Imag().begin(), nd_size, et.Imag().begin());
      std::copy_n(e0.Real().begin() + nd_size, h1_size, en.Real().begin());
      std::copy_n(e0.Imag().begin() + nd_size, h1_size, en.Imag().begin());
    }

    // Power-normalize eigenvector.
    {
      const auto *Btt = mode_op.GetBtt();
      if (Btt)
      {
        Vector Btt_etr(nd_size), Btt_eti(nd_size);
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
    if (i < num_modes && is_propagating && estimator)
    {
      ComplexVector bz(l2_size);
      {
        Vector curl_etr(l2_size), curl_eti(l2_size);
        CurlOp.Mult(et.Real(), curl_etr);
        CurlOp.Mult(et.Imag(), curl_eti);
        bz.Real() = curl_eti;
        bz.Real() *= 1.0 / omega;
        bz.Imag() = curl_etr;
        bz.Imag() *= -1.0 / omega;
      }
      estimator->AddErrorIndicator(et, bz, total_domain_energy, indicator);
    }
  }
  Mpi::Print("\n");

  post_op.MeasureFinalize(indicator);
  return {indicator, nd_fespace.GlobalTrueVSize() + h1_fespace.GlobalTrueVSize()};
}

}  // namespace palace
