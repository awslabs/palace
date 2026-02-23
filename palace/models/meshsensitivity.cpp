// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "meshsensitivity.hpp"

#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <set>
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/spaceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/iodata.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// Compute v^H * (A * u) for a possibly-null operator A. Returns 0 if A is null.
std::complex<double> SesquilinearForm(MPI_Comm comm, const ComplexOperator *A,
                                      const ComplexVector &v, const ComplexVector &u,
                                      ComplexVector &Au)
{
  if (!A)
  {
    return 0.0;
  }
  A->Mult(u, Au);
  return linalg::Dot(comm, Au, v);
}

// Get the VDof index for a scalar DOF and component in a vector-valued FE space.
int NodeDofIndex(const mfem::FiniteElementSpace &fes, int scalar_dof, int component)
{
  return fes.DofToVDof(scalar_dof, component);
}

// Core FD loop: perturb each boundary vertex in each spatial direction, evaluate the
// system form, and compute the FD derivative.
MeshSensitivity::SensitivityResult
ComputeFDSensitivities(const IoData &iodata, SpaceOperator &space_op,
                       const ComplexVector &u, const ComplexVector &v,
                       std::complex<double> omega, bool has_damping, double perturbation,
                       const char *label)
{
  MeshSensitivity::SensitivityResult result;
  MPI_Comm comm = space_op.GetComm();
  mfem::ParMesh &mesh = space_op.GetNDSpace().GetMesh().Get();
  const int sdim = mesh.SpaceDimension();

  // Determine FD step size.
  double h = perturbation;
  if (h <= 0.0)
  {
    h = MeshSensitivity::DefaultPerturbation(mesh, comm);
  }

  // Get boundary vertex indices (local to this process).
  auto bdr_verts = MeshSensitivity::GetBoundaryVertexIndices(mesh);
  result.num_bdr_nodes = static_cast<int>(bdr_verts.size());
  result.bdr_node_indices = bdr_verts;
  result.dJ_dx.SetSize(sdim * result.num_bdr_nodes);
  result.dJ_dx = 0.0;

  Mpi::Print(
      "\nComputing mesh sensitivities via finite differences ({})\n"
      " Boundary nodes on this process: {:d}, FD step size: {:.6e}\n",
      label, result.num_bdr_nodes, h);

  // Compute baseline: v^H * A(x) * u.
  std::complex<double> val_0 =
      MeshSensitivity::EvalSystemForm(space_op, u, v, omega, has_damping, comm);
  Mpi::Print(" Baseline v^H A u = ({:.6e}, {:.6e})\n", val_0.real(), val_0.imag());

  // Access mesh node coordinates. For isoparametric elements, mesh nodes are stored in a
  // GridFunction on the mesh's nodal FE space.
  mfem::GridFunction *nodes = mesh.GetNodes();
  MFEM_VERIFY(nodes, "Mesh has no nodal grid function for mesh sensitivity computation!");
  const mfem::FiniteElementSpace &node_fes = *nodes->FESpace();
  MFEM_VERIFY(node_fes.GetVDim() == sdim,
              "Nodal FE space vdim must match space dimension for mesh sensitivity!");

  // Perturb each boundary vertex in each spatial direction.
  // NOTE: In parallel, vertices on processor boundaries are only perturbed locally. This
  // may cause inconsistencies for shared boundary nodes. For the FD prototype this is
  // acceptable; the future analytic implementation will handle this correctly.
  for (int ni = 0; ni < result.num_bdr_nodes; ni++)
  {
    const int vert_idx = bdr_verts[ni];

    for (int d = 0; d < sdim; d++)
    {
      // Get the VDof index in the nodal GridFunction for this vertex and spatial component.
      // In MFEM's H1 nodal FE space, vertex v has scalar DOF index v (vertices are the
      // first NV DOFs). DofToVDof maps (scalar_dof, component) -> vector_dof_index.
      const int vdof = NodeDofIndex(node_fes, vert_idx, d);

      // Perturb the node coordinate.
      const double orig_val = (*nodes)(vdof);
      (*nodes)(vdof) = orig_val + h;

      // Evaluate perturbed system form.
      std::complex<double> val_h =
          MeshSensitivity::EvalSystemForm(space_op, u, v, omega, has_damping, comm);

      // Restore the node coordinate.
      (*nodes)(vdof) = orig_val;

      // Forward FD approximation: dJ/dx_{i,d} = (J(x+h*e) - J(x)) / h.
      std::complex<double> dval = (val_h - val_0) / h;
      result.dJ_dx(ni * sdim + d) = dval.real();
    }

    if ((ni + 1) % 100 == 0 || ni + 1 == result.num_bdr_nodes)
    {
      Mpi::Print(" Processed {:d}/{:d} boundary nodes\n", ni + 1, result.num_bdr_nodes);
    }
  }

  // Restore geometry cache to baseline state.
  space_op.GetNDSpace().GetMesh().ResetCeedObjects();

  return result;
}

}  // namespace

std::vector<int> MeshSensitivity::GetBoundaryVertexIndices(const mfem::ParMesh &mesh)
{
  // Collect vertex indices from all boundary elements on this process.
  std::set<int> bdr_verts;
  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    mfem::Array<int> verts;
    mesh.GetBdrElementVertices(i, verts);
    for (int j = 0; j < verts.Size(); j++)
    {
      bdr_verts.insert(verts[j]);
    }
  }
  return {bdr_verts.begin(), bdr_verts.end()};
}

double MeshSensitivity::DefaultPerturbation(const mfem::ParMesh &mesh, MPI_Comm comm)
{
  // Compute characteristic length from vertex coordinate ranges.
  const int sdim = mesh.SpaceDimension();
  const int nv = mesh.GetNV();
  mfem::Vector bb_min(sdim), bb_max(sdim);
  bb_min = std::numeric_limits<double>::max();
  bb_max = std::numeric_limits<double>::lowest();
  for (int i = 0; i < nv; i++)
  {
    const double *coord = mesh.GetVertex(i);
    for (int d = 0; d < sdim; d++)
    {
      bb_min(d) = std::min(bb_min(d), coord[d]);
      bb_max(d) = std::max(bb_max(d), coord[d]);
    }
  }
  // Reduce across processes to get the global bounding box.
  Mpi::GlobalMin(sdim, bb_min.GetData(), comm);
  Mpi::GlobalMax(sdim, bb_max.GetData(), comm);

  double diag_sq = 0.0;
  for (int d = 0; d < sdim; d++)
  {
    double diff = bb_max(d) - bb_min(d);
    diag_sq += diff * diff;
  }
  double diag = std::sqrt(diag_sq);
  // Use 1e-7 * characteristic_length as the default FD step.
  return 1.0e-7 * diag;
}

std::complex<double> MeshSensitivity::EvalSystemForm(SpaceOperator &space_op,
                                                     const ComplexVector &u,
                                                     const ComplexVector &v,
                                                     std::complex<double> omega,
                                                     bool has_damping, MPI_Comm comm)
{
  // Reset cached geometry data so operators are assembled with current mesh node positions.
  space_op.GetNDSpace().GetMesh().ResetCeedObjects();

  // Assemble K, C, M from current mesh configuration.
  auto K = space_op.GetStiffnessMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
  auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);

  // Compute v^H * A * u where A = K + i*omega*C - omega^2*M.
  const int n = u.Size();
  ComplexVector Au(n);
  Au.UseDevice(true);

  std::complex<double> val = 0.0;
  val += SesquilinearForm(comm, K.get(), v, u, Au);
  if (has_damping && C)
  {
    val += (1i * omega) * SesquilinearForm(comm, C.get(), v, u, Au);
  }
  if (M)
  {
    val += (-omega * omega) * SesquilinearForm(comm, M.get(), v, u, Au);
  }

  return val;
}

MeshSensitivity::SensitivityResult
MeshSensitivity::ComputeEigen(const IoData &iodata, SpaceOperator &space_op,
                               const ComplexVector &u, std::complex<double> omega,
                               bool has_damping, double perturbation)
{
  // For the eigenvalue problem, the adjoint vector is the eigenvector itself.
  // The FD loop computes u^H (dA/dx) u = numerator of dmu/dx or dlambda/dx.
  // We must convert to df/dx by applying the denominator and frequency conversion.
  auto result = ComputeFDSensitivities(iodata, space_op, u, u, omega, has_damping,
                                       perturbation, "eigenvalue");

  // Compute the denominator for the eigenvalue sensitivity formula.
  MPI_Comm comm = space_op.GetComm();
  const int n = u.Size();
  ComplexVector Au(n);
  Au.UseDevice(true);

  if (!has_damping)
  {
    // Linear EVP: mu = omega^2. dmu/dx = numerator / (u^H M u).
    // Then df/dx = dmu/dx / (2*omega) / (2*pi).
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    std::complex<double> uMu = 1.0;
    if (M)
    {
      M->Mult(u, Au);
      uMu = linalg::Dot(comm, Au, u);
    }
    // Conversion factor: raw -> df in GHz.
    // raw = u^H dA/dx u (eigenvalue-scale per mesh-length)
    // dmu/dx = raw / uMu
    // domega/dx = dmu/dx / (2*omega)
    // df/dx = domega/dx / (2*pi)
    std::complex<double> factor = 1.0 / (uMu * 2.0 * omega * 2.0 * M_PI);
    Mpi::Print(" Eigenvalue conversion: u^H M u = ({:.6e}, {:.6e}), omega = ({:.6e}, "
               "{:.6e})\n",
               uMu.real(), uMu.imag(), omega.real(), omega.imag());
    const int sdim = space_op.GetNDSpace().GetMesh().Get().SpaceDimension();
    for (int i = 0; i < result.num_bdr_nodes * sdim; i++)
    {
      result.dJ_dx(i) *= factor.real();
    }
  }
  else
  {
    // Quadratic EVP: lambda = i*omega. dlambda/dx = -numerator / (u^H(C+2*lambda*M)u).
    // Then df/dx = (-i/(2*pi)) * dlambda/dx.
    std::complex<double> lambda = std::complex<double>(0.0, 1.0) * omega;
    auto C = space_op.GetDampingMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    std::complex<double> denom = 0.0;
    if (C)
    {
      C->Mult(u, Au);
      denom += linalg::Dot(comm, Au, u);
    }
    if (M)
    {
      M->Mult(u, Au);
      denom += 2.0 * lambda * linalg::Dot(comm, Au, u);
    }
    // factor = -1/denom * (-i/(2*pi)) = i / (denom * 2*pi)
    std::complex<double> factor =
        std::complex<double>(0.0, 1.0) / (denom * 2.0 * M_PI);
    const int sdim = space_op.GetNDSpace().GetMesh().Get().SpaceDimension();
    for (int i = 0; i < result.num_bdr_nodes * sdim; i++)
    {
      result.dJ_dx(i) *= factor.real();
    }
  }

  return result;
}

MeshSensitivity::SensitivityResult
MeshSensitivity::ComputeDriven(const IoData &iodata, SpaceOperator &space_op,
                                const ComplexVector &u, const ComplexVector &v,
                                double omega, double perturbation)
{
  // For the driven problem, use A = K + i*omega*C - omega^2*M with real omega.
  const std::complex<double> omega_c(omega, 0.0);
  constexpr bool has_damping = true;
  return ComputeFDSensitivities(iodata, space_op, u, v, omega_c, has_damping, perturbation,
                                "driven");
}

void MeshSensitivity::WriteSensitivitySolb(const std::string &filename,
                                            const SensitivityResult &result,
                                            mfem::ParMesh &mesh, MPI_Comm comm)
{
  // Gather the full mesh to rank 0 for serial vertex count and coordinate mapping.
  // The serial mesh has a consistent global vertex ordering.
  mfem::Mesh smesh = mesh.GetSerialMesh(0);

  if (!Mpi::Root(comm))
  {
    return;
  }

  const int sdim = smesh.SpaceDimension();
  const int nv = smesh.GetNV();

  // Build the full-mesh sensitivity vector: zero everywhere, then scatter boundary values.
  // The SensitivityResult stores boundary-only sensitivities indexed by the local boundary
  // vertex list. We expand to all nv vertices so the .solb has the same vertex count as
  // the companion .meshb file.
  //
  // For single-process runs, the boundary vertex indices in result.bdr_node_indices are
  // directly usable as global indices. For multi-process runs, the serial mesh gathered
  // by GetSerialMesh(0) has a vertex ordering that matches the original serial mesh —
  // which is the same ordering that metric::WriteMeshMeshb produces. Since the FD
  // computation currently runs independently on each process with local indices, the
  // sensitivities should be gathered and mapped to serial indices in a production
  // implementation. For this initial version (designed for single-process verification),
  // we use the boundary indices directly.
  mfem::Vector full_sens(nv * sdim);
  full_sens = 0.0;

  for (int i = 0; i < result.num_bdr_nodes; i++)
  {
    const int vi = result.bdr_node_indices[i];
    if (vi < nv)
    {
      for (int d = 0; d < sdim; d++)
      {
        full_sens(vi * sdim + d) = result.dJ_dx(i * sdim + d);
      }
    }
  }

  // Write the .solb file with a vector field (field type 2 = GmfVec).
  std::ofstream fo(filename);
  MFEM_VERIFY(fo.is_open(),
              "Unable to open sensitivity output file \"" << filename << "\"!");

  fo << "MeshVersionFormatted 2\n";
  fo << "Dimension " << sdim << "\n\n";
  fo << "SolAtVertices\n";
  // Header: num_vertices num_fields field_type
  fo << nv << " 1 2\n";

  fo << std::scientific;
  fo.precision(15);

  for (int i = 0; i < nv; i++)
  {
    for (int d = 0; d < sdim; d++)
    {
      if (d > 0)
      {
        fo << ' ';
      }
      fo << full_sens(i * sdim + d);
    }
    fo << '\n';
  }

  fo << "\nEnd\n";
  fo.close();

  Mpi::Print(" Wrote mesh sensitivity dJ/dx to {}\n"
             " Total vertices: {:d}, boundary vertices: {:d}\n",
             filename, nv, result.num_bdr_nodes);
}

}  // namespace palace
