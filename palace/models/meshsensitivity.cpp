// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "meshsensitivity.hpp"

#include <cmath>
#include <complex>
#include <fstream>
#include <limits>
#include <set>
#include "fem/integ/meshsensitivity.hpp"
#include "linalg/operator.hpp"
#include "linalg/vector.hpp"
#include "models/materialoperator.hpp"
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
  MPI_Comm comm = space_op.GetComm();
  const int n = u.Size();

  const auto &mat = space_op.GetMaterialOp();
  if (!has_damping && !mat.HasLondonDepth())
  {
    // Linear EVP (lossless, no London superconductors): analytic sensitivity via a single
    // CeedOperator pass.
    // The sensitivity QFunction differentiates the integrand f = f_K - lambda * f_M
    // w.r.t. the Jacobian J at each quadrature point, using Enzyme AD or dual numbers.
    //
    // The raw output is df/dx_nodes (L-vector on the mesh nodal FE space). The
    // CeedOperator's element restriction transpose handles DOF assembly (summing
    // contributions from elements sharing a node). We index into the L-vector directly
    // via DofToVDof, then apply the eigenvalue-to-frequency conversion factor.

    const double lambda_val = (omega * omega).real();
    auto &mesh = space_op.GetNDSpace().GetMesh();
    const auto &nd_fespace = space_op.GetNDSpace();
    const mfem::ParMesh &pmesh = mesh.Get();
    const mfem::GridFunction &mesh_nodes = *pmesh.GetNodes();
    const int sdim = pmesh.SpaceDimension();

    // The analytic sensitivity path currently only supports single-process runs.
    // In parallel, shared boundary vertices need MPI reduction and the local-to-serial
    // vertex mapping in WriteSensitivitySolb needs proper handling.
    MFEM_VERIFY(Mpi::Size(comm) == 1,
                "Analytic mesh sensitivity currently requires single-process execution. "
                "Use the FD fallback for multi-process runs.");

    // For a lossless eigenvalue problem, the eigenvector is real up to a global phase
    // factor: u = e^{iθ} * u_real. Rotate the eigenvector to maximize the real part by
    // finding θ = arg(u^T u) and multiplying by e^{-iθ/2}. This makes u effectively real
    // so that the real-part-only QFunction captures the full sensitivity.
    ComplexVector u_rotated(u);
    {
      // Compute the "phase angle" via u^T u (not u^H u, which is always real).
      // u^T u = |u_r|^2 - |u_i|^2 + 2i*(u_r · u_i). The angle of this complex number
      // is 2θ where θ is the global phase. Rotating by -θ aligns u with the real axis.
      double ur_dot_ur = 0.0, ui_dot_ui = 0.0, ur_dot_ui = 0.0;
      const mfem::Vector &ur = u.Real(), &ui = u.Imag();
      for (int i = 0; i < ur.Size(); i++)
      {
        ur_dot_ur += ur(i) * ur(i);
        ui_dot_ui += ui(i) * ui(i);
        ur_dot_ui += ur(i) * ui(i);
      }
      double two_theta = std::atan2(2.0 * ur_dot_ui, ur_dot_ur - ui_dot_ui);
      double theta = two_theta / 2.0;
      double cos_t = std::cos(theta), sin_t = std::sin(theta);

      // u_rotated = u * e^{-iθ} = (u_r + i*u_i)(cos(-θ) + i*sin(-θ))
      //           = (u_r*cosθ + u_i*sinθ) + i*(-u_r*sinθ + u_i*cosθ)
      for (int i = 0; i < ur.Size(); i++)
      {
        u_rotated.Real()(i) = ur(i) * cos_t + ui(i) * sin_t;
        u_rotated.Imag()(i) = -ur(i) * sin_t + ui(i) * cos_t;
      }

      double new_imag_norm = u_rotated.Imag().Norml2();
      double new_real_norm = u_rotated.Real().Norml2();
      Mpi::Print(" Eigenvector phase rotation: theta = {:.4f} rad, "
                 "|Im|/|Re| before = {:.2e}, after = {:.2e}\n",
                 theta, ui.Norml2() / ur.Norml2(), new_imag_norm / new_real_norm);

      if (new_real_norm > 0.0 && new_imag_norm / new_real_norm > 0.01)
      {
        Mpi::Warning("Eigenvector imaginary component remains significant after phase "
                     "rotation (|Im|/|Re| = {:.2e}). This indicates a truly degenerate "
                     "mode. Mesh sensitivity for degenerate modes requires handling both "
                     "real and imaginary components.\n",
                     new_imag_norm / new_real_norm);
      }
    }

    Mpi::Print("\nComputing mesh sensitivities via analytic differentiation (eigenvalue)\n"
               " lambda = omega^2 = {:.6e}, method = "
#ifdef PALACE_WITH_ENZYME
               "Enzyme reverse-mode AD"
#else
               "forward-mode dual numbers"
#endif
               "\n",
               lambda_val);

    // Extract real part of the phase-rotated eigenvector for the integrator.
    const mfem::Vector &u_real = u_rotated.Real();

    // Build material coefficients matching what SpaceOperator uses for K and M.
    // Stiffness (curl-curl): mu_inv (inverse permeability).
    // Mass (H(curl) Piola): eps (real permittivity).
    // London superconductor term is excluded here (gated by HasLondonDepth check above)
    // because it enters the stiffness matrix as a vector-mass coefficient (adjJt/detJ
    // Piola), not a curl-curl coefficient. Supporting it requires a three-coefficient
    // QFunction context.
    MaterialPropertyCoefficient Q_stiff(mat.MaxCeedAttribute());
    Q_stiff.AddCoefficient(mat.GetAttributeToMaterial(), mat.GetInvPermeability());
    MaterialPropertyCoefficient Q_mass(mat.MaxCeedAttribute());
    Q_mass.AddCoefficient(mat.GetAttributeToMaterial(), mat.GetPermittivityReal());

    // Apply the sensitivity CeedOperator: mesh_nodes → raw sensitivity (L-vector).
    // The CeedOperator element restriction transpose sums contributions from elements
    // sharing a DOF. The output is an L-vector indexed by DofToVDof.
    const mfem::FiniteElementSpace &mesh_fespace = *mesh_nodes.FESpace();
    mfem::Vector raw_sens(mesh_fespace.GetVSize());
    MeshSensitivityIntegrator::Apply(mesh, nd_fespace, &Q_mass, &Q_stiff, lambda_val,
                                     u_real, raw_sens);

    // Compute u^H M u for the eigenvalue-to-frequency conversion factor.
    ComplexVector Au(n);
    Au.UseDevice(true);
    auto M = space_op.GetMassMatrix<ComplexOperator>(Operator::DIAG_ZERO);
    std::complex<double> uMu = 1.0;
    if (M)
    {
      M->Mult(u, Au);
      uMu = linalg::Dot(comm, Au, u);
    }
    // Conversion chain:
    //   raw = d(lambda_nondim)/d(x_nondim)  [from CeedOperator]
    //   d(lambda_nondim)/d(x_mesh) = raw / Lc0   [since x_nondim = x_mesh / Lc0]
    //   d(f_GHz)/d(x_mesh) = d(lambda_nondim)/d(x_mesh) * freq_scale / (uMu * 2*omega * 2*pi)
    // where freq_scale = 1/tc converts nondim freq to GHz, and Lc0 = Lc/L0 converts
    // nondimensional mesh coordinates to mesh file units.
    const double freq_scale =
        iodata.units.GetScaleFactor<Units::ValueType::FREQUENCY>();
    const double Lc0 = iodata.units.GetMeshLengthRelativeScale();
    std::complex<double> factor =
        freq_scale / (Lc0 * uMu * 2.0 * omega * 2.0 * M_PI);
    Mpi::Print(" Eigenvalue conversion: u^H M u = ({:.6e}, {:.6e}), omega = ({:.6e}, "
               "{:.6e})\n"
               " Unit scaling: freq_scale = {:.6e} GHz, Lc0 = {:.6e}\n",
               uMu.real(), uMu.imag(), omega.real(), omega.imag(), freq_scale, Lc0);

    // Build result with ALL vertex sensitivities (not just boundary). The full-mesh
    // sensitivity is needed for the chain contraction dJ/dp = (dJ/dx) · (dx/dp), where
    // dx/dp includes interior mesh deformation that follows the boundary motion.
    SensitivityResult result;
    const int nv = pmesh.GetNV();
    result.num_bdr_nodes = nv;  // All vertices, despite the field name
    result.bdr_node_indices.resize(nv);
    result.dJ_dx.SetSize(sdim * nv);

    for (int vi = 0; vi < nv; vi++)
    {
      result.bdr_node_indices[vi] = vi;
      for (int d = 0; d < sdim; d++)
      {
        const int vdof = mesh_fespace.DofToVDof(vi, d);
        result.dJ_dx(vi * sdim + d) = raw_sens(vdof) * factor.real();
      }
    }
    Mpi::Print(" Full-mesh sensitivity: {:d} vertices, L2 norm = {:.6e} GHz/cm\n",
               nv, result.dJ_dx.Norml2());

    return result;
  }
  else
  {
    // Quadratic EVP (lossy) or London superconductors: analytic sensitivity not yet
    // implemented for the damping term C or the London depth vector-mass contribution.
    // Fall back to finite differences.
    Mpi::Warning("Analytic mesh sensitivity not yet supported for this configuration "
                 "(quadratic EVP or London depth), falling back to finite differences\n");
    auto result = ComputeFDSensitivities(iodata, space_op, u, u, omega, has_damping,
                                         perturbation, "eigenvalue (FD fallback)");

    // Apply the eigenvalue-to-frequency conversion for quadratic EVP.
    ComplexVector Au(n);
    Au.UseDevice(true);
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
    std::complex<double> factor =
        std::complex<double>(0.0, 1.0) / (denom * 2.0 * M_PI);
    const int sdim = space_op.GetNDSpace().GetMesh().Get().SpaceDimension();
    for (int i = 0; i < result.num_bdr_nodes * sdim; i++)
    {
      result.dJ_dx(i) *= factor.real();
    }

    return result;
  }
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
