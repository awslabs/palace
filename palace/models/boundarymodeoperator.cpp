// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "boundarymodeoperator.hpp"

#include <algorithm>
#include "fem/multigrid.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/waveportoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

BoundaryModeOperator::BoundaryModeOperator(const IoData &iodata_,
                                           const std::vector<std::unique_ptr<Mesh>> &mesh)
  : iodata(iodata_), solver_order(iodata_.solver.order),
    use_submesh(!iodata_.solver.boundary_mode.attributes.empty()), solve_mesh(nullptr)
{
  // Phase 1: Set up mesh (extract submesh from 3D or use direct 2D).
  SetUpMesh(mesh);

  // Phase 2: Material operator.
  mat_op = std::make_unique<MaterialOperator>(iodata, *solve_mesh);
  if (use_submesh)
  {
    mat_op->RotateMaterialTensors(iodata, submesh_e1, submesh_e2, submesh_normal);
  }

  // Phase 3: FE spaces.
  SetUpFESpaces(mesh);

  // Phase 4: Boundary operators.
  auto &pmesh = solve_mesh->Get();
  surf_z_op = std::make_unique<SurfaceImpedanceOperator>(iodata, *mat_op, pmesh);
  farfield_op = std::make_unique<FarfieldBoundaryOperator>(iodata, *mat_op, pmesh);
  surf_sigma_op = std::make_unique<SurfaceConductivityOperator>(iodata, *mat_op, pmesh);

  // Phase 5: Eigenvalue solver.
  SetUpEigenSolver();

  // Phase 6: Error estimator.
  const int dim = solve_mesh->Dimension();
  rt_fec_est = std::make_unique<mfem::RT_FECollection>(solver_order - 1, dim);
  nd_fespaces_est = std::make_unique<FiniteElementSpaceHierarchy>(
      std::make_unique<FiniteElementSpace>(*solve_mesh, nd_fecs.back().get()));
  rt_fespaces_est = std::make_unique<FiniteElementSpaceHierarchy>(
      std::make_unique<FiniteElementSpace>(*solve_mesh, rt_fec_est.get()));
  h1_fespaces_est = std::make_unique<FiniteElementSpaceHierarchy>(
      std::make_unique<FiniteElementSpace>(*solve_mesh, h1_fecs.back().get()));
  estimator = std::make_unique<BoundaryModeFluxErrorEstimator<ComplexVector>>(
      *mat_op, *nd_fespaces_est, *rt_fespaces_est, *l2_curl_fespace, *h1_fespaces_est,
      iodata.solver.linear.estimator_tol, iodata.solver.linear.estimator_max_it, 0,
      iodata.solver.linear.estimator_mg);

  Mpi::Print(" ND space: {:d} DOFs, H1 space: {:d} DOFs, total: {:d}\n",
             GetNDSpace().GlobalTrueVSize(), GetH1Space().GlobalTrueVSize(),
             GetNDSpace().GlobalTrueVSize() + GetH1Space().GlobalTrueVSize());
}

void BoundaryModeOperator::SetUpMesh(const std::vector<std::unique_ptr<Mesh>> &mesh)
{
  const auto &bm_data = iodata.solver.boundary_mode;

  if (use_submesh)
  {
    MFEM_VERIFY(mesh.back()->Dimension() == 3,
                "BoundaryMode with \"Attributes\" requires a 3D mesh!");
    Mpi::Print(" Extracting 2D submesh from 3D boundary attributes...\n");
    const auto &parent_mesh = mesh.back()->Get();
    MPI_Comm comm = parent_mesh.GetComm();
    mfem::Array<int> attr_list;
    attr_list.Append(bm_data.attributes.data(), bm_data.attributes.size());

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

    // Extract a standalone 2D serial mesh from the 3D parallel boundary, then redistribute.
    // The serial roundtrip (parallel → serial → parallel) is needed because building a
    // standalone 2D mesh requires global topology: domain attribute remapping from parent
    // volume elements, internal boundary edges at BC intersections, and consistent 2D
    // coordinate projection. ParSubMesh::CreateFromBoundary doesn't provide these.
    // Note: does not support nonconforming (NCMesh) meshes.
    auto serial_mesh = mesh::ExtractStandalone2DSubmesh(
        parent_mesh, attr_list, internal_bdr_attrs, submesh_normal, submesh_centroid,
        submesh_e1, submesh_e2);
    if (!Mpi::Root(comm))
    {
      serial_mesh.reset();
    }
    auto mesh_2d = mesh::DistributeSerialMesh(comm, serial_mesh);
    owned_mesh = std::make_unique<Mesh>(std::move(mesh_2d));
    solve_mesh = owned_mesh.get();
  }
  else
  {
    MFEM_VERIFY(mesh.back()->Dimension() == 2,
                "BoundaryMode solver requires a 2D mesh (waveguide cross-section), "
                "or a 3D mesh with \"Attributes\" specifying the cross-section boundary!");
    solve_mesh = mesh.back().get();
  }
}

void BoundaryModeOperator::SetUpFESpaces(const std::vector<std::unique_ptr<Mesh>> &mesh)
{
  const auto &mg = iodata.solver.linear;
  const int dim = solve_mesh->Dimension();

  // Collect Dirichlet boundary attributes.
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
        const auto &bm_attrs = iodata.solver.boundary_mode.attributes;
        for (auto attr : data.attributes)
        {
          if (std::find(bm_attrs.begin(), bm_attrs.end(), attr) != bm_attrs.end())
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

  // FE collections.
  nd_fecs = fem::ConstructFECollections<mfem::ND_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);
  h1_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);
  h1_aux_fecs = fem::ConstructFECollections<mfem::H1_FECollection>(
      solver_order, dim, mg.mg_max_levels, mg.mg_coarsening, false);

  // For submesh: temporarily move owned_mesh into a vector for hierarchy construction.
  // For direct 2D: use the caller's mesh vector directly.
  std::vector<std::unique_ptr<Mesh>> submesh_vec;
  if (use_submesh)
  {
    submesh_vec.push_back(std::move(owned_mesh));
  }
  const auto &fespace_mesh = use_submesh ? submesh_vec : mesh;

  nd_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::ND_FECollection>(
      mg.mg_max_levels, fespace_mesh, nd_fecs, &dbc_bcs, &nd_dbc_tdof_lists);
  h1_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
      mg.mg_max_levels, fespace_mesh, h1_fecs, &dbc_bcs, &h1_dbc_tdof_lists);
  h1_aux_fespaces = fem::ConstructFiniteElementSpaceHierarchy<mfem::H1_FECollection>(
      mg.mg_max_levels, fespace_mesh, h1_aux_fecs, &dbc_bcs, &h1_aux_dbc_tdof_lists);

  if (use_submesh)
  {
    owned_mesh = std::move(submesh_vec[0]);
  }

  // L2 curl space for 2D B-field.
  l2_curl_fec = std::make_unique<mfem::L2_FECollection>(
      solver_order - 1, dim, mfem::BasisType::GaussLegendre, mfem::FiniteElement::INTEGRAL);
  l2_curl_fespace = std::make_unique<FiniteElementSpace>(*solve_mesh, l2_curl_fec.get());
}

void BoundaryModeOperator::SetUpEigenSolver()
{
  const auto &bm_data = iodata.solver.boundary_mode;
  auto &nd_fespace = GetNDSpace();
  auto &h1_fespace = GetH1Space();
  const int nd_size = nd_fespace.GetTrueVSize();

  const auto which_eig = (bm_data.target > 0.0)
                             ? EigenvalueSolver::WhichType::LARGEST_MAGNITUDE
                             : EigenvalueSolver::WhichType::LARGEST_REAL;

  // Combined DBC tdof list for the block system.
  mfem::Array<int> dbc_tdof_list;
  dbc_tdof_list.Append(nd_dbc_tdof_lists.back());
  for (int i = 0; i < h1_dbc_tdof_lists.back().Size(); i++)
  {
    dbc_tdof_list.Append(nd_size + h1_dbc_tdof_lists.back()[i]);
  }

  // Multigrid config.
  std::unique_ptr<ModeEigenSolverMultigridConfig> mg_config_ptr;
  if (nd_fespaces.GetNumLevels() > 1)
  {
    mg_config_ptr = std::make_unique<ModeEigenSolverMultigridConfig>();
    mg_config_ptr->nd_fespaces = &nd_fespaces;
    mg_config_ptr->h1_fespaces = &h1_fespaces;
    mg_config_ptr->h1_aux_fespaces = &h1_aux_fespaces;
    mg_config_ptr->nd_dbc_tdof_lists = &nd_dbc_tdof_lists;
    mg_config_ptr->h1_dbc_tdof_lists = &h1_dbc_tdof_lists;
    mg_config_ptr->h1_aux_dbc_tdof_lists = &h1_aux_dbc_tdof_lists;
    Mpi::Print(" Using p-multigrid preconditioning with {:d} levels\n",
               static_cast<int>(nd_fespaces.GetNumLevels()));
  }

  port_data = std::make_unique<WavePortData>(
      *mat_op, nullptr, surf_z_op.get(), farfield_op.get(), surf_sigma_op.get(), nd_fespace,
      h1_fespace, dbc_tdof_list, bm_data.n, bm_data.max_size, bm_data.tol, which_eig,
      iodata.solver.linear, bm_data.type, iodata.problem.verbose, solve_mesh->GetComm(),
      mg_config_ptr.get());
}

BoundaryModeOperator::SolveResult BoundaryModeOperator::Solve(double omega,
                                                              double kn_target)
{
  double sigma = -kn_target * kn_target;
  auto result = port_data->GetEigenSolver().Solve(omega, sigma);
  return {result.num_converged, sigma};
}

std::complex<double> BoundaryModeOperator::GetEigenvalue(int i) const
{
  return port_data->GetEigenSolver().GetEigenvalue(i);
}

void BoundaryModeOperator::GetEigenvector(int i, ComplexVector &x) const
{
  port_data->GetEigenSolver().GetEigenvector(i, x);
}

double BoundaryModeOperator::GetError(int i, EigenvalueSolver::ErrorType type) const
{
  return port_data->GetEigenSolver().GetError(i, type);
}

const mfem::HypreParMatrix *BoundaryModeOperator::GetBtt() const
{
  return port_data->GetEigenSolver().GetBtt();
}

const mfem::HypreParMatrix *BoundaryModeOperator::GetAtnr() const
{
  return port_data->GetEigenSolver().GetAtnr();
}

const mfem::HypreParMatrix *BoundaryModeOperator::GetAtni() const
{
  return port_data->GetEigenSolver().GetAtni();
}

const ComplexKspSolver *BoundaryModeOperator::GetLinearSolver() const
{
  return port_data->GetEigenSolver().GetLinearSolver();
}

void BoundaryModeOperator::AddErrorIndicator(const ComplexVector &et,
                                             const ComplexVector &bz,
                                             double total_domain_energy,
                                             ErrorIndicator &indicator)
{
  estimator->AddErrorIndicator(et, bz, total_domain_energy, indicator);
}

}  // namespace palace
