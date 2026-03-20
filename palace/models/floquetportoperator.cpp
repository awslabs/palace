// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "floquetportoperator.hpp"

#include <cmath>
#include <fmt/core.h>
#include <mfem.hpp>
#include "fem/integrator.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/configfile.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/prettyprint.hpp"

namespace palace
{

using namespace std::complex_literals;

// =============================================================================
// LowRankComplexOperator
// =============================================================================

void LowRankComplexOperator::Mult(const ComplexVector &x, ComplexVector &y) const
{
  y = 0.0;
  AddMult(x, y, 1.0);
}

void LowRankComplexOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                     const std::complex<double> a) const
{
  // Palace uses BILINEAR forms (complex-symmetric, NOT Hermitian). The Floquet DtN
  // operator has the form F = Σ g_k conj(v_k) v_k^T. Applied to x:
  //   F x = Σ g_k conj(v_k) (v_k^T x)
  // where v^T x is the bilinear (NOT Hermitian) inner product.
  for (const auto &t : terms)
  {
    // Bilinear dot product: dot = v^T x = Σ v_j x_j (complex multiplication, no conj)
    // = (v_r · x_r - v_i · x_i) + i(v_r · x_i + v_i · x_r)
    double dot_r =
        linalg::Dot(comm, t.v->Real(), x.Real()) - linalg::Dot(comm, t.v->Imag(), x.Imag());
    double dot_i =
        linalg::Dot(comm, t.v->Real(), x.Imag()) + linalg::Dot(comm, t.v->Imag(), x.Real());

    // Compute s = a * g * dot (complex scalar)
    std::complex<double> s = a * t.g * std::complex<double>(dot_r, dot_i);

    // y += s * conj(v) = (s_r + i s_i)(v_r - i v_i)
    //    Real part: s_r v_r + s_i v_i
    //    Imag part: s_i v_r - s_r v_i
    y.Real().Add(s.real(), t.v->Real());
    y.Real().Add(s.imag(), t.v->Imag());
    y.Imag().Add(s.imag(), t.v->Real());
    y.Imag().Add(-s.real(), t.v->Imag());
  }
}

// =============================================================================
// MassConsistentDtNOperator
// =============================================================================

void MassConsistentDtNOperator::Mult(const ComplexVector &x, ComplexVector &y) const
{
  y = 0.0;
  AddMult(x, y, 1.0);
}

void MassConsistentDtNOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                        const std::complex<double> a) const
{
  // Apply F*x = Σ_i g_i * A_i * x where A_i is a real HypreParMatrix.
  // Since A_i is real: (a * g) * A * x splits into real/imag components.
  for (const auto &t : terms)
  {
    std::complex<double> s = a * t.g;
    if (std::abs(s) < 1e-30)
    {
      continue;
    }

    // tmp = A * x_real
    t.A->Mult(x.Real(), tmp);
    y.Real().Add(s.real(), tmp);
    y.Imag().Add(s.imag(), tmp);

    // tmp = A * x_imag
    t.A->Mult(x.Imag(), tmp);
    y.Real().Add(-s.imag(), tmp);
    y.Imag().Add(s.real(), tmp);
  }
}

namespace
{

void ComputePolarization(const mfem::Vector &kt, const mfem::Vector &port_normal,
                         mfem::Vector &e_te, mfem::Vector &e_tm);

}  // namespace

// =============================================================================
// DenseBoundaryOperator
// =============================================================================

DenseBoundaryOperator::DenseBoundaryOperator(MPI_Comm port_comm, int n_full,
                                             const mfem::Array<int> &bdr_tdof_list_in,
                                             int n_bdr_global)
  : ComplexOperator(n_full), port_comm(port_comm), n_full(n_full),
    n_bdr_global(n_bdr_global), F(Eigen::MatrixXcd::Zero(n_bdr_global, n_bdr_global)),
    x_bdr_local(0), x_bdr_global(n_bdr_global), y_bdr_global(n_bdr_global)
{
  bdr_tdof_list.MakeRef(bdr_tdof_list_in);
  n_bdr_local = bdr_tdof_list.Size();

  // Set up Allgatherv counts and displacements on port_comm.
  if (port_comm != MPI_COMM_NULL)
  {
    int nproc = Mpi::Size(port_comm);
    recv_counts.resize(nproc);
    recv_displs.resize(nproc);
    MPI_Allgather(&n_bdr_local, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, port_comm);
    recv_displs[0] = 0;
    for (int i = 1; i < nproc; i++)
    {
      recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
    }
  }
  x_bdr_local.resize(n_bdr_local);
}

void DenseBoundaryOperator::Mult(const ComplexVector &x, ComplexVector &y) const
{
  y = 0.0;
  AddMult(x, y, 1.0);
}

void DenseBoundaryOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                    const std::complex<double> a) const
{
  if (port_comm == MPI_COMM_NULL)
  {
    return;  // Non-port rank: nothing to do.
  }

  // Restrict x to local boundary DOFs.
  x.Real().HostRead();
  x.Imag().HostRead();
  for (int i = 0; i < n_bdr_local; i++)
  {
    int dof = bdr_tdof_list[i];
    x_bdr_local(i) = std::complex<double>(x.Real()[dof], x.Imag()[dof]);
  }

  // Allgatherv on port_comm to get global boundary DOF vector.
  MPI_Allgatherv(x_bdr_local.data(), n_bdr_local, MPI_DOUBLE_COMPLEX, x_bdr_global.data(),
                 recv_counts.data(), recv_displs.data(), MPI_DOUBLE_COMPLEX, port_comm);

  // Dense matvec in boundary DOF space.
  y_bdr_global.noalias() = F * x_bdr_global;

  // Scatter local portion back to y with scaling a.
  y.Real().HostReadWrite();
  y.Imag().HostReadWrite();
  int offset = recv_displs[Mpi::Rank(port_comm)];
  for (int i = 0; i < n_bdr_local; i++)
  {
    int dof = bdr_tdof_list[i];
    auto val = a * y_bdr_global(offset + i);
    y.Real()[dof] += val.real();
    y.Imag()[dof] += val.imag();
  }
}

// =============================================================================
// FloquetPortData
// =============================================================================

FloquetPortData::FloquetPortData(const config::FloquetPortData &data, const IoData &iodata,
                                 const MaterialOperator &mat_op_ref,
                                 mfem::ParFiniteElementSpace &nd_fespace)
  : excitation(data.excitation), active(data.active), mat_op(mat_op_ref), a1(3), a2(3),
    b1(3), b2(3), k_F(3), port_normal(3), comm(nd_fespace.GetComm())
{
  // Set incident polarization coefficients: E_inc = α_TE ê_TE + α_TM ê_TM.
  const double inv_sqrt2 = 1.0 / std::sqrt(2.0);
  if (data.inc_polarization == "TE")
  {
    inc_alpha_te = 1.0;
    inc_alpha_tm = 0.0;
  }
  else if (data.inc_polarization == "TM")
  {
    inc_alpha_te = 0.0;
    inc_alpha_tm = 1.0;
  }
  else if (data.inc_polarization == "RHC")
  {
    inc_alpha_te = inv_sqrt2;
    inc_alpha_tm = std::complex<double>(0.0, inv_sqrt2);
  }
  else if (data.inc_polarization == "LHC")
  {
    inc_alpha_te = inv_sqrt2;
    inc_alpha_tm = std::complex<double>(0.0, -inv_sqrt2);
  }
  else
  {
    MFEM_ABORT("Unknown IncidentPolarization: " << data.inc_polarization
                                                << ". Expected TE, TM, RHC, or LHC.");
  }

  // Store boundary attributes.
  attr_list.SetSize(static_cast<int>(data.attributes.size()));
  for (int i = 0; i < attr_list.Size(); i++)
  {
    attr_list[i] = data.attributes[i];
  }

  // Floquet ports require periodic boundary conditions in the transverse directions.
  // The periodic mesh provides DOF identification on opposite faces, and the Floquet
  // wave vector (if nonzero) provides the Bloch phase shift.
  const auto &periodic = iodata.boundaries.periodic;
  MFEM_VERIFY(!periodic.boundary_pairs.empty(),
              "FloquetPort requires periodic boundary conditions to be configured under "
              "\"Boundaries\"/\"Periodic\"/\"BoundaryPairs\". At least two periodic "
              "boundary pairs are needed for 3D periodicity.");

  // Extract lattice vectors. If periodic boundary pairs with translations are specified,
  // use them. Otherwise, infer from the port face bounding box (works for axis-aligned
  // rectangular cells with built-in mesh periodicity).
  double mesh_scale = iodata.units.GetMeshLengthRelativeScale();
  if (periodic.boundary_pairs.size() >= 2)
  {
    // The translation component of each affine transform is the lattice vector (in mesh
    // units). Nondimensionalize by dividing by Lc/L0.
    for (int i = 0; i < 3; i++)
    {
      a1(i) = periodic.boundary_pairs[0].affine_transform[i * 4 + 3] / mesh_scale;
      a2(i) = periodic.boundary_pairs[1].affine_transform[i * 4 + 3] / mesh_scale;
    }
  }
  // Fall back to port geometry if translations are missing (e.g., auto-detected periodicity
  // from Gmsh setPeriodic, or no BoundaryPairs at all).
  if (a1.Norml2() < 1e-12 || a2.Norml2() < 1e-12)
  {
    // Infer lattice vectors from the port face bounding box. This works for axis-aligned
    // rectangular periodic cells where the port face spans the full unit cell.
    auto &mesh = *nd_fespace.GetParMesh();
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    Mpi::GlobalMax(1, &bdr_attr_max, comm);
    mfem::Array<int> bdr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list, true);
    mfem::Vector bbmin, bbmax;
    mesh::GetAxisAlignedBoundingBox(mesh, bdr_marker, true, bbmin, bbmax);
    mfem::Vector extent(3);
    for (int i = 0; i < 3; i++)
    {
      extent(i) = bbmax(i) - bbmin(i);
    }

    // The two tangential directions with the largest extents are the lattice directions.
    // Sort dimensions by extent to identify the two periodic directions.
    int dirs[3] = {0, 1, 2};
    std::sort(dirs, dirs + 3, [&](int a, int b) { return extent(a) > extent(b); });

    a1 = 0.0;
    a2 = 0.0;
    a1(dirs[0]) = extent(dirs[0]);
    a2(dirs[1]) = extent(dirs[1]);

    MFEM_VERIFY(a1.Norml2() > 1e-12 && a2.Norml2() > 1e-12,
                "Could not infer lattice vectors from port face bounding box. "
                "Please specify BoundaryPairs with Translation vectors.");

    Mpi::Print(" Floquet port: inferred lattice vectors from port geometry:\n"
               "   a1 = ({:.4e}, {:.4e}, {:.4e})\n"
               "   a2 = ({:.4e}, {:.4e}, {:.4e})\n",
               a1(0), a1(1), a1(2), a2(0), a2(1), a2(2));
  }

  // Bloch wave vector: use the BZ-wrapped value from MaterialOperator so the port is
  // consistent with the bulk ∇_F bilinear form.
  const auto &kF_wrapped = mat_op.GetWaveVector();
  for (int i = 0; i < 3; i++)
  {
    k_F(i) = (i < kF_wrapped.Size()) ? kF_wrapped(i) : 0.0;
  }

  // Compute reciprocal lattice.
  ComputeReciprocalLattice(a1, a2, b1, b2);
  Mpi::Print(" Floquet port reciprocal lattice:\n"
             "   b1 = ({:.4e}, {:.4e}, {:.4e})\n"
             "   b2 = ({:.4e}, {:.4e}, {:.4e})\n",
             b1(0), b1(1), b1(2), b2(0), b2(1), b2(2));

  // Compute the BZ wrapping offset using the BZ-wrapped k_F (not k₀).
  {
    mfem::Vector kF_unwrapped(3);
    for (int i = 0; i < 3; i++)
    {
      kF_unwrapped(i) = periodic.wave_vector[i];
    }
    const auto &kF_bz = mat_op.GetWaveVectorBZ();
    mfem::Vector kF_wrapped(3);
    for (int i = 0; i < 3; i++)
    {
      kF_wrapped(i) = (i < kF_bz.Size()) ? kF_bz(i) : 0.0;
    }
    bz_m = ComputeBZOffset(kF_unwrapped, kF_wrapped, b1, b1 * b1);
    bz_n = ComputeBZOffset(kF_unwrapped, kF_wrapped, b2, b2 * b2);
    if (bz_m != 0 || bz_n != 0)
    {
      Mpi::Print(" Floquet port: BZ wrapping active (shift = {:d}*b1 + {:d}*b2)\n", bz_m,
                 bz_n);
    }
  }

  // Compute port normal and area using geodata utilities.
  {
    auto &mesh = *nd_fespace.GetParMesh();
    int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
    Mpi::GlobalMax(1, &bdr_attr_max, comm);
    mfem::Array<int> bdr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list, true);
    port_normal = mesh::GetSurfaceNormal(mesh, bdr_marker);
    port_area = mesh::GetSurfaceArea(mesh, bdr_marker);
  }

  // Material properties at the port: find the volume element adjacent to any boundary
  // element on the port face and look up its mu_r, eps_r. In parallel, not all ranks may
  // have port boundary elements — use MPI reduction to broadcast the found attribute.
  {
    auto &mesh = *nd_fespace.GetParMesh();
    int port_vol_attr = -1;
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      int battr = mesh.GetBdrAttribute(be);
      bool on_port = false;
      for (int i = 0; i < attr_list.Size(); i++)
      {
        if (battr == attr_list[i])
        {
          on_port = true;
          break;
        }
      }
      if (!on_port)
      {
        continue;
      }
      int f, o;
      mesh.GetBdrElementFace(be, &f, &o);
      int iel1, iel2;
      mesh.GetFaceElements(f, &iel1, &iel2);
      int iel = (iel1 >= 0) ? iel1 : iel2;
      if (iel >= 0)
      {
        port_vol_attr = mesh.GetAttribute(iel);
        break;
      }
    }
    // Look up material properties on the rank that found the boundary element, then
    // broadcast. The CEED attribute map is rank-local (only includes attributes of elements
    // owned by this rank), so mat_op.GetLightSpeedMin(attr) would fail on ranks that don't
    // own elements with that volume attribute.
    double c_local = -1.0;
    double muinv_local = -1.0;
    if (port_vol_attr > 0)
    {
      c_local = mat_op.GetLightSpeedMin(port_vol_attr);
      muinv_local = mat_op.GetInvPermeabilityZZ(port_vol_attr);
    }

    // Broadcast from the rank that found a valid attribute.
    int found_rank = (port_vol_attr > 0) ? Mpi::Rank(comm) : Mpi::Size(comm);
    Mpi::GlobalMin(1, &found_rank, comm);
    MFEM_VERIFY(found_rank < Mpi::Size(comm),
                "Could not find volume element adjacent to Floquet port!");
    MPI_Bcast(&c_local, 1, MPI_DOUBLE, found_rank, comm);
    MPI_Bcast(&muinv_local, 1, MPI_DOUBLE, found_rank, comm);

    MFEM_VERIFY(c_local > 0.0, "Invalid material speed of light at Floquet port!");
    mu_eps_port = 1.0 / (c_local * c_local);
    mu_r_port = 1.0 / muinv_local;

    Mpi::Print("   Port material: mu_r = {:.4e}, mu_r*eps_r = {:.4e}, "
               "c = {:.4e}\n",
               mu_r_port, mu_eps_port, c_local);
  }

  // Determine max diffraction order.
  if (data.max_order >= 0)
  {
    max_order_m = data.max_order;
    max_order_n = data.max_order;
  }
  else
  {
    // Auto-compute: start with a safe default, will be refined in Initialize().
    max_order_m = 3;
    max_order_n = 3;
  }

  // Parse config options.
  use_full_dtn = data.full_dtn;
  use_auxiliary = (data.port_mode == "Auxiliary");
  use_mass_consistent = (data.port_mode == "MassConsistent");

  // Compute boundary DOF infrastructure for DenseBoundaryOperator and preconditioner.
  // Create a sub-communicator for ranks that own port boundary elements (wave port
  // pattern).
  {
    auto &mesh = *nd_fespace.GetParMesh();
    int bdr_attr_max = 0;
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      bdr_attr_max = std::max(bdr_attr_max, mesh.GetBdrAttribute(be));
    }
    Mpi::GlobalMax(1, &bdr_attr_max, nd_fespace.GetComm());
    auto bdr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list, true);
    nd_fespace.GetEssentialTrueDofs(bdr_marker, bdr_tdof_list);

    // Compute global boundary DOF count on the global communicator BEFORE splitting
    // (all ranks must participate in the same collective to avoid deadlock).
    n_bdr_global = bdr_tdof_list.Size();
    Mpi::GlobalSum(1, &n_bdr_global, nd_fespace.GetComm());

    // Create sub-communicator for ranks with port boundary DOFs.
    int color = (bdr_tdof_list.Size() > 0) ? 0 : MPI_UNDEFINED;
    MPI_Comm_split(nd_fespace.GetComm(), color, Mpi::Rank(nd_fespace.GetComm()),
                   &port_comm);

    if (port_comm != MPI_COMM_NULL)
    {
      int nproc = Mpi::Size(port_comm);
      recv_counts.resize(nproc);
      recv_displs.resize(nproc);
      int n_local = bdr_tdof_list.Size();
      MPI_Allgather(&n_local, 1, MPI_INT, recv_counts.data(), 1, MPI_INT, port_comm);
      recv_displs[0] = 0;
      for (int i = 1; i < nproc; i++)
      {
        recv_displs[i] = recv_displs[i - 1] + recv_counts[i - 1];
      }
    }
  }

  // Cap MaxOrder at the mesh Nyquist limit: modes with |B| > π/h can't be resolved
  // by the FEM basis and their Fourier projections decay exponentially (with accurate
  // quadrature) or produce aliasing noise (with default quadrature). Including them is
  // at best wasteful and at worst actively harmful. Per literature (Li et al. 2019):
  // include modes up to ~period/h, over-integrate boundary terms.
  {
    // Compute h_max on port boundary elements. b1, b2 are in reciprocal mesh-units
    // (computed from a1, a2 which use GetVertex coordinates). h_max is in mesh-units
    // (also from GetVertex). So |b| × h_max is dimensionless — no unit conversion needed.
    auto &mesh_tmp = *nd_fespace.GetParMesh();
    int bdr_attr_max_tmp = 0;
    for (int be = 0; be < mesh_tmp.GetNBE(); be++)
    {
      bdr_attr_max_tmp = std::max(bdr_attr_max_tmp, mesh_tmp.GetBdrAttribute(be));
    }
    Mpi::GlobalMax(1, &bdr_attr_max_tmp, nd_fespace.GetComm());
    auto bdr_marker_tmp = mesh::AttrToMarker(bdr_attr_max_tmp, attr_list, true);

    double h_max_local = 0.0;
    for (int be = 0; be < mesh_tmp.GetNBE(); be++)
    {
      if (bdr_marker_tmp[mesh_tmp.GetBdrAttribute(be) - 1])
      {
        // Use the boundary element Jacobian to compute element size.
        // J is spaceDim × (Dim-1) = 3×2 for a surface element in 3D.
        // The max singular value of J gives the element's largest extent (h_max).
        auto *T = mesh_tmp.GetBdrElementTransformation(be);
        T->SetIntPoint(
            &mfem::Geometries.GetCenter(T->GetGeometryType()));
        const mfem::DenseMatrix &J = T->Jacobian();
        double h = J.CalcSingularvalue(0);  // largest singular value
        h_max_local = std::max(h_max_local, h);
      }
    }
    double h_max = h_max_local;
    Mpi::GlobalMax(1, &h_max, nd_fespace.GetComm());

    if (h_max > 0.0)
    {
      // h_max is already in nondimensional (Lc) units since the mesh is
      // nondimensionalized at load time. b1, b2 are also in Lc⁻¹ units.
      // No further unit conversion needed.
      double h_max_nondim = h_max;

      // FE Nyquist: p-th order elements can resolve modes with |B|×h < p×π.
      // b1, b2 are in nondimensional units; h_max_nondim is now consistent.
      int p = nd_fespace.GetMaxElementOrder();
      double b1_norm = b1.Norml2();
      double b2_norm = b2.Norml2();
      int nyquist_m =
          std::max(1, static_cast<int>(std::floor(p * M_PI / (b1_norm * h_max_nondim))));
      int nyquist_n =
          std::max(1, static_cast<int>(std::floor(p * M_PI / (b2_norm * h_max_nondim))));
      Mpi::Print(" Floquet port: h_max = {:.4e} (nondim, Lc units), p = {:d}, "
                 "Nyquist limit = ({:d}, {:d})\n",
                 h_max_nondim, p, nyquist_m, nyquist_n);
      if (max_order_m > nyquist_m || max_order_n > nyquist_n)
      {
        Mpi::Print(" Floquet port: capping MaxOrder from ({:d}, {:d}) to ({:d}, {:d})\n",
                   max_order_m, max_order_n, std::min(max_order_m, nyquist_m),
                   std::min(max_order_n, nyquist_n));
        max_order_m = std::min(max_order_m, nyquist_m);
        max_order_n = std::min(max_order_n, nyquist_n);
      }
    }
  }

  // When FullDtN is active, DISABLE the Nyquist cap. The aliased above-Nyquist
  // modes (with default quadrature) provide additional boundary DOF coverage that
  // prevents the PMC-like reflection from rank-deficient DtN.
  // if (use_full_dtn)
  //{
  //  max_order_m = std::max(max_order_m, 30);
  //  max_order_n = std::max(max_order_n, 30);
  //  Mpi::Print(" Floquet port FullDtN: N_bdr = {:d}, max orders = ({:d}, {:d})\n",
  //             n_bdr_global, max_order_m, max_order_n);
  //}

  // Enumerate diffraction orders and assemble Fourier projections.
  EnumerateOrders();
  AssembleFourierProjections(nd_fespace);
  nd_fespace_ptr = &nd_fespace;  // Store for re-assembly in Initialize.

  Mpi::Print(" Floquet port: {:d} modes ({:d} orders x 2 polarizations), "
             "N_bdr = {:d}, port area = {:.4e}, normal = ({:.3f}, {:.3f}, {:.3f})\n",
             static_cast<int>(modes.size()), static_cast<int>(modes.size()) / 2,
             n_bdr_global, port_area, port_normal(0), port_normal(1), port_normal(2));

  // Parseval diagnostic: measure completeness of the mode basis by testing
  // Σ conj(v)(v^T v_test)/|Γ| against M_bdr v_test. For a complete mode set,
  // the Parseval identity gives Σ conj(v)v^T/|Γ| = M_bdr (boundary mass matrix).
  // The ratio ||Σ ...|| / ||M_bdr v_test|| measures true Parseval completeness.
  if (!modes.empty())
  {
    // Use the (0,0) TE mode's v as test vector.
    const ComplexVector *v_test = nullptr;
    for (const auto &m : modes)
    {
      if (m.m == 0 && m.n == 0 && m.is_te)
      {
        v_test = &m.v;
        break;
      }
    }
    if (v_test)
    {
      int n = v_test->Size();

      // Assemble boundary mass matrix M_bdr for the correct Parseval reference.
      auto &mesh_diag = *nd_fespace.GetParMesh();
      int bdr_attr_max_diag = 0;
      for (int be = 0; be < mesh_diag.GetNBE(); be++)
      {
        bdr_attr_max_diag = std::max(bdr_attr_max_diag, mesh_diag.GetBdrAttribute(be));
      }
      Mpi::GlobalMax(1, &bdr_attr_max_diag, nd_fespace.GetComm());
      auto bdr_marker_diag = mesh::AttrToMarker(bdr_attr_max_diag, attr_list, true);
      mfem::ParBilinearForm m_bdr_form(&nd_fespace);
      m_bdr_form.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(), bdr_marker_diag);
      m_bdr_form.Assemble();
      m_bdr_form.Finalize();
      auto M_bdr_mat = std::unique_ptr<mfem::HypreParMatrix>(m_bdr_form.ParallelAssemble());

      // Compute M_bdr v_test as the Parseval reference.
      Vector Mv_r(n), Mv_i(n);
      M_bdr_mat->Mult(v_test->Real(), Mv_r);
      M_bdr_mat->Mult(v_test->Imag(), Mv_i);
      double Mv_norm =
          std::sqrt(linalg::Dot(comm, Mv_r, Mv_r) + linalg::Dot(comm, Mv_i, Mv_i));

      // Compute rank-1 sum applied to v_test: y = Σ conj(v_i) (v_i^T v_test) / |Γ|
      ComplexVector y_rank1(n);
      y_rank1 = 0.0;
      for (const auto &mode : modes)
      {
        if (!HasFlag(mode.use, FloquetModeUse::Dtn))
        {
          continue;
        }
        // Bilinear dot: v^T v_test
        double dot_r = linalg::Dot(comm, mode.v.Real(), v_test->Real()) -
                       linalg::Dot(comm, mode.v.Imag(), v_test->Imag());
        double dot_i = linalg::Dot(comm, mode.v.Real(), v_test->Imag()) +
                       linalg::Dot(comm, mode.v.Imag(), v_test->Real());
        std::complex<double> dot(dot_r, dot_i);
        std::complex<double> s = dot / port_area;
        // y += s * conj(v) = (s_r+is_i)(v_r-iv_i)
        y_rank1.Real().Add(s.real(), mode.v.Real());
        y_rank1.Real().Add(s.imag(), mode.v.Imag());
        y_rank1.Imag().Add(s.imag(), mode.v.Real());
        y_rank1.Imag().Add(-s.real(), mode.v.Imag());
      }
      // Full complex norm of y (not just real part).
      double y_norm = std::sqrt(linalg::Dot(comm, y_rank1.Real(), y_rank1.Real()) +
                                linalg::Dot(comm, y_rank1.Imag(), y_rank1.Imag()));
      double completeness = (Mv_norm > 0.0) ? y_norm / Mv_norm : 0.0;
      Mpi::Print(" Parseval diagnostic: ||Σ conj(v)(v^Tv_test)/|Γ|||={:.4e} "
                 "||M_bdr v_test||={:.4e} completeness={:.4f}\n",
                 y_norm, Mv_norm, completeness);
    }

    // Cross-polarization diagnostic: check v_TE · v_TM for (0,0) modes.
    // On a tet mesh at k=0 (B=0), these may not be orthogonal, causing
    // energy non-conservation in S-parameter extraction.
    const ComplexVector *v_te = nullptr, *v_tm = nullptr;
    for (const auto &m : modes)
    {
      if (m.m == 0 && m.n == 0)
      {
        if (m.is_te)
        {
          v_te = &m.v;
        }
        else
        {
          v_tm = &m.v;
        }
      }
    }
    if (v_te && v_tm)
    {
      double vte_norm = std::sqrt(linalg::Dot(comm, v_te->Real(), v_te->Real()) +
                                  linalg::Dot(comm, v_te->Imag(), v_te->Imag()));
      double vtm_norm = std::sqrt(linalg::Dot(comm, v_tm->Real(), v_tm->Real()) +
                                  linalg::Dot(comm, v_tm->Imag(), v_tm->Imag()));
      // Bilinear dot: v_TE^T v_TM
      double cross_r = linalg::Dot(comm, v_te->Real(), v_tm->Real()) -
                       linalg::Dot(comm, v_te->Imag(), v_tm->Imag());
      double cross_i = linalg::Dot(comm, v_te->Real(), v_tm->Imag()) +
                       linalg::Dot(comm, v_te->Imag(), v_tm->Real());
      double cross_mag = std::sqrt(cross_r * cross_r + cross_i * cross_i);
      Mpi::Print(" Cross-pol: ||v_TE||={:.4e} ||v_TM||={:.4e} "
                 "|v_TE^T v_TM|={:.4e} cos_angle={:.4f}\n",
                 vte_norm, vtm_norm, cross_mag, cross_mag / (vte_norm * vtm_norm));
    }

    // Normalization diagnostic: check if v^T e(ê) / |Γ| = 1 for the (0,0) TE mode.
    // Here e(ê) is the FE DOF vector obtained by projecting ê_TE onto the ND space.
    // The S-parameter extraction computes c = v^T E / |Γ|. For an exact plane wave
    // E = ê_TE on the port, c should equal 1. If v^T e(ê)/|Γ| ≠ 1, there's a systematic
    // normalization error in the S-parameters.
    if (v_te)
    {
      // Find the (0,0) TE mode to get e_pol.
      const mfem::Vector *e_pol_ptr = nullptr;
      for (const auto &m : modes)
      {
        if (m.m == 0 && m.n == 0 && m.is_te)
        {
          e_pol_ptr = &m.e_pol;
          break;
        }
      }
      if (e_pol_ptr)
      {
        // Project ê_TE onto the ND FE space (volume interpolation).
        mfem::VectorConstantCoefficient e_te_coeff(*e_pol_ptr);
        mfem::ParGridFunction e_te_gf(&nd_fespace);
        e_te_gf.ProjectCoefficient(e_te_coeff);

        // Get true DOF vector.
        int tdof = nd_fespace.GetTrueVSize();
        Vector e_te_tdof(tdof);
        nd_fespace.GetProlongationMatrix()->MultTranspose(e_te_gf, e_te_tdof);

        // v^T e(ê) — bilinear product (v is real for B=0).
        double vTe = linalg::Dot(comm, v_te->Real(), e_te_tdof);

        // Also compute ∫_Γ |ê_TE_h|² dS = e^T M_bdr e (the FE L² norm of ê on boundary).
        auto &mesh_norm = *nd_fespace.GetParMesh();
        int bdr_attr_max_norm = 0;
        for (int be = 0; be < mesh_norm.GetNBE(); be++)
        {
          bdr_attr_max_norm = std::max(bdr_attr_max_norm, mesh_norm.GetBdrAttribute(be));
        }
        Mpi::GlobalMax(1, &bdr_attr_max_norm, nd_fespace.GetComm());
        auto bdr_marker_norm = mesh::AttrToMarker(bdr_attr_max_norm, attr_list, true);
        mfem::ParBilinearForm m_bdr_norm(&nd_fespace);
        m_bdr_norm.AddBoundaryIntegrator(new mfem::VectorFEMassIntegrator(),
                                         bdr_marker_norm);
        m_bdr_norm.Assemble();
        m_bdr_norm.Finalize();
        auto M_bdr_norm =
            std::unique_ptr<mfem::HypreParMatrix>(m_bdr_norm.ParallelAssemble());
        Vector Me(tdof);
        M_bdr_norm->Mult(e_te_tdof, Me);
        double eTMe = linalg::Dot(comm, e_te_tdof, Me);

        Mpi::Print(" Normalization: v^T e(ê)/|Γ| = {:.6f}, "
                   "e^T M_bdr e / |Γ| = {:.6f} (both should be 1.0)\n",
                   vTe / port_area, eTMe / port_area);
      }
    }
  }
}

void FloquetPortData::ComputeReciprocalLattice(const mfem::Vector &a1,
                                               const mfem::Vector &a2, mfem::Vector &b1,
                                               mfem::Vector &b2)
{
  // For 2D periodicity in 3D space, the reciprocal lattice vectors satisfy:
  //   a_i . b_j = 2*pi * delta_ij
  //   b_i is perpendicular to a_j (j != i) and lies in the plane of a1, a2.

  // n = a1 x a2 (normal to the lattice plane)
  mfem::Vector n(3);
  n(0) = a1(1) * a2(2) - a1(2) * a2(1);
  n(1) = a1(2) * a2(0) - a1(0) * a2(2);
  n(2) = a1(0) * a2(1) - a1(1) * a2(0);
  double vol = n.Norml2();  // |a1 x a2| = area of unit cell

  MFEM_VERIFY(vol > 0.0, "Lattice vectors a1, a2 are degenerate (zero cross product)!");

  // b1 = 2*pi * (a2 x n) / |a1 x a2|^2
  // b2 = 2*pi * (n x a1) / |a1 x a2|^2
  mfem::Vector a2xn(3), nxa1(3);
  a2xn(0) = a2(1) * n(2) - a2(2) * n(1);
  a2xn(1) = a2(2) * n(0) - a2(0) * n(2);
  a2xn(2) = a2(0) * n(1) - a2(1) * n(0);

  nxa1(0) = n(1) * a1(2) - n(2) * a1(1);
  nxa1(1) = n(2) * a1(0) - n(0) * a1(2);
  nxa1(2) = n(0) * a1(1) - n(1) * a1(0);

  double vol_sq = vol * vol;
  b1.SetSize(3);
  b2.SetSize(3);
  for (int i = 0; i < 3; i++)
  {
    b1(i) = 2.0 * M_PI * a2xn(i) / vol_sq;
    b2(i) = 2.0 * M_PI * nxa1(i) / vol_sq;
  }

  // Verify: a_i . b_j = 2*pi * delta_ij.
  MFEM_VERIFY(std::abs(a1 * b1 - 2.0 * M_PI) < 1e-10 &&
                  std::abs(a2 * b2 - 2.0 * M_PI) < 1e-10,
              "Reciprocal lattice computation failed: diagonal check!");
  MFEM_VERIFY(std::abs(a1 * b2) < 1e-10 && std::abs(a2 * b1) < 1e-10,
              "Reciprocal lattice computation failed: off-diagonal check!");
}

int FloquetPortData::ComputeBZOffset(const mfem::Vector &kF_unwrapped,
                                     const mfem::Vector &kF_wrapped, const mfem::Vector &b,
                                     double b_sq)
{
  mfem::Vector dkF(kF_unwrapped.Size());
  for (int i = 0; i < kF_unwrapped.Size(); i++)
  {
    dkF(i) = kF_unwrapped(i) - kF_wrapped(i);
  }
  return static_cast<int>(std::round(dkF * b / b_sq));
}

void FloquetPortData::EnumerateOrders()
{
  modes.clear();

  // Enumerate physical mode indices. The range is the union of ±MaxOrder around the
  // specular (0,0) and the BZ-shifted range. This ensures the user always sees their
  // requested diffraction orders, while also including modes near the BZ boundary that
  // are critical for the DtN correction when BZ wrapping is active.
  int m_lo = std::min(-max_order_m, -max_order_m + bz_m);
  int m_hi = std::max(max_order_m, max_order_m + bz_m);
  int n_lo = std::min(-max_order_n, -max_order_n + bz_n);
  int n_hi = std::max(max_order_n, max_order_n + bz_n);
  for (int m = m_lo; m <= m_hi; m++)
  {
    for (int n = n_lo; n <= n_hi; n++)
    {
      // B_mn for the Fourier projection kernel: shifted by BZ offset so that
      // exp(-j B_mn · r) extracts the correct coefficient from E_p (which was solved
      // with the wrapped k_F). Physical mode labels (m, n) are preserved.
      mfem::Vector B_mn(3);
      for (int i = 0; i < 3; i++)
      {
        B_mn(i) = (m - bz_m) * b1(i) + (n - bz_n) * b2(i);
      }

      // Propagation direction k_hat for this order.
      // k_t = B_mn - k_F (transverse wavevector relative to the periodic field)
      // With BZ shift: k_t = (m-bz_m)*b1 + (n-bz_n)*b2 - k_F_wrapped
      //              = m*b1 + n*b2 - k_F_unwrapped (same physical k_t)
      // gamma^2 = omega^2*mu*eps - |k_t|^2, computed in Initialize(omega).

      // Compute TE and TM polarization vectors from the transverse wavevector.
      // Use the physical k_F at the reference frequency for initial assembly.
      // Initialize() will update the polarization when k_F changes with frequency.
      mfem::Vector kt(3);
      {
        const auto &kF_phys = mat_op.HasFloquetFrequencyScaling() ? mat_op.GetWaveVectorBZ()
                                                                  : mat_op.GetWaveVector();
        for (int i = 0; i < 3; i++)
        {
          kt(i) = B_mn(i) - ((i < kF_phys.Size()) ? kF_phys(i) : 0.0);
        }
      }
      mfem::Vector e_te(3), e_tm(3);
      ComputePolarization(kt, port_normal, e_te, e_tm);

      // User-requested modes (±MaxOrder around specular) are for S-parameter output.
      // BZ-shifted modes (±MaxOrder around BZ origin) are for the DtN correction.
      bool in_user_range = (std::abs(m) <= max_order_m && std::abs(n) <= max_order_n);
      bool in_dtn_range =
          (std::abs(m - bz_m) <= max_order_m && std::abs(n - bz_n) <= max_order_n);
      auto use = static_cast<FloquetModeUse>(
          (in_user_range ? static_cast<uint8_t>(FloquetModeUse::Output) : 0) |
          (in_dtn_range ? static_cast<uint8_t>(FloquetModeUse::Dtn) : 0));

      // Add TE and TM modes for this order.
      for (bool is_te : {true, false})
      {
        FloquetMode mode;
        mode.m = m;
        mode.n = n;
        mode.is_te = is_te;
        mode.use = use;
        mode.B_mn = B_mn;
        mode.e_pol = is_te ? e_te : e_tm;
        mode.gamma_sq = 0.0;
        modes.push_back(std::move(mode));
      }
    }
  }
}

void FloquetPortData::AssembleFourierProjections(mfem::ParFiniteElementSpace &nd_fespace)
{
  // For each mode, assemble: v_j = int_Gamma (n x [n x N_j]) . e_p exp(-i B_mn . r) dS
  // The tangential projection n x (n x N_j) extracts the component of N_j tangential to
  // the port face. For ND basis functions on a boundary face, VectorFEBoundaryLFIntegrator
  // already computes int_Gamma (f . N_j) dS where f is a VectorCoefficient. Since ND basis
  // functions on a face are tangential, (n x [n x N_j]) . e_p = N_j . e_p_tangential.
  //
  // So we need: v_j = int_Gamma N_j(r) . [e_p * exp(-i B_mn . r)] dS
  //
  // Split into real and imaginary parts:
  //   Re: int_Gamma N_j . [e_p * cos(B_mn . r)] dS
  //   Im: int_Gamma N_j . [e_p * (-sin(B_mn . r))] dS

  // Create boundary attribute marker. Scan all local boundary elements for the actual max
  // attribute (may include internally-added interface boundary elements).
  auto &mesh = *nd_fespace.GetParMesh();
  int bdr_attr_max = 0;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    bdr_attr_max = std::max(bdr_attr_max, mesh.GetBdrAttribute(be));
  }
  Mpi::GlobalMax(1, &bdr_attr_max, nd_fespace.GetComm());
  mfem::Array<int> bdr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list, true);

  // Compute maximum boundary element diameter for quadrature order estimation.
  // The Fourier integrand exp(-i B . r) oscillates with period 2π/|B|. For modes where
  // |B|×h >> 1, the default quadrature can't resolve the oscillations and produces aliased
  // v vectors. Over-integration eliminates aliasing — above-Nyquist v vectors correctly
  // decay to ~0, and the DtN sum converges exponentially (per literature: Li et al. 2019).
  double h_max_local = 0.0;
  mfem::Geometry::Type bdr_geom = mfem::Geometry::INVALID;
  for (int be = 0; be < mesh.GetNBE(); be++)
  {
    if (bdr_marker[mesh.GetBdrAttribute(be) - 1])
    {
      auto *T = mesh.GetBdrElementTransformation(be);
      T->SetIntPoint(
          &mfem::Geometries.GetCenter(T->GetGeometryType()));
      double h = T->Jacobian().CalcSingularvalue(0);
      h_max_local = std::max(h_max_local, h);
      bdr_geom = mesh.GetBdrElementGeometry(be);
    }
  }
  double h_max = h_max_local;
  Mpi::GlobalMax(1, &h_max, nd_fespace.GetComm());

  int tdof_size = nd_fespace.GetTrueVSize();

  for (auto &mode : modes)
  {
    mode.v.SetSize(tdof_size);
    mode.v = 0.0;
    mode.v.UseDevice(true);

    // Coefficient for real part: e_p * cos(B_mn . r)
    // Coefficient for imag part: e_p * (-sin(B_mn . r))
    // We use a lambda-based VectorCoefficient.

    class FloquetModeCoeff : public mfem::VectorCoefficient
    {
    public:
      const mfem::Vector &e_pol, &B_mn;
      bool is_real;  // true = cos part, false = -sin part

      FloquetModeCoeff(const mfem::Vector &e, const mfem::Vector &B, bool real_part)
        : mfem::VectorCoefficient(3), e_pol(e), B_mn(B), is_real(real_part)
      {
      }

      void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
                const mfem::IntegrationPoint &ip) override
      {
        mfem::Vector x(3);
        T.Transform(ip, x);
        double phase = 0.0;
        for (int i = 0; i < 3; i++)
        {
          phase += B_mn(i) * x(i);
        }
        double scale = is_real ? std::cos(phase) : -std::sin(phase);
        V.SetSize(3);
        for (int i = 0; i < 3; i++)
        {
          V(i) = e_pol(i) * scale;
        }
      }
    };

    FloquetModeCoeff coeff_r(mode.e_pol, mode.B_mn, true);
    FloquetModeCoeff coeff_i(mode.e_pol, mode.B_mn, false);

    // Compute extra quadrature order to resolve oscillatory exp(-i B . r).
    // Need ~2 points per oscillation cycle: extra_order ≈ |B| × h / π.
    // Over-integration ensures accurate v vectors, which is critical for modes
    // with large DtN coefficients (near-cutoff evanescent TM: g ∝ ω²με/|γ|).
    double B_norm = mode.B_mn.Norml2();
    int extra_order = 0; //static_cast<int>(std::ceil(B_norm * h_max / M_PI));
    const mfem::IntegrationRule *ir = nullptr;
    if (extra_order > 0 && bdr_geom != mfem::Geometry::INVALID)
    {
      int fe_order = nd_fespace.GetMaxElementOrder();
      ir = &mfem::IntRules.Get(bdr_geom, 2 * fe_order + extra_order);
      Mpi::Print("    Mode ({:+d},{:+d},{}) |B|={:.2f} extra_order={:d} "
                 "quad_order={:d}\n",
                 mode.m, mode.n, mode.is_te ? "TE" : "TM", B_norm, extra_order,
                 2 * fe_order + extra_order);
    }

    // Assemble real part.
    auto *integ_r = new VectorFEBoundaryLFIntegrator(coeff_r);
    if (ir)
    {
      integ_r->SetIntRule(ir);
    }
    mfem::LinearForm lf_r(&nd_fespace);
    lf_r.AddBoundaryIntegrator(integ_r, bdr_marker);
    lf_r.UseFastAssembly(false);
    lf_r.UseDevice(false);
    lf_r.Assemble();
    lf_r.UseDevice(true);

    // Assemble imaginary part.
    auto *integ_i = new VectorFEBoundaryLFIntegrator(coeff_i);
    if (ir)
    {
      integ_i->SetIntRule(ir);
    }
    mfem::LinearForm lf_i(&nd_fespace);
    lf_i.AddBoundaryIntegrator(integ_i, bdr_marker);
    lf_i.UseFastAssembly(false);
    lf_i.UseDevice(false);
    lf_i.Assemble();
    lf_i.UseDevice(true);

    // Project to true DOFs: v = P^T * lf.
    nd_fespace.GetProlongationMatrix()->MultTranspose(lf_r, mode.v.Real());
    nd_fespace.GetProlongationMatrix()->MultTranspose(lf_i, mode.v.Imag());

    // Mass-consistent mode: assemble A_i[j,k] = ∫_Γ (N_j·ê_i)(N_k·ê_i) dS.
    // The complex exponentials exp(±jB·r) cancel in the bilinear form, so A_i is real
    // and frequency-independent. Uses ê ⊗ ê as a MatrixCoefficient in
    // VectorFEMassIntegrator.
    if (use_mass_consistent && HasFlag(mode.use, FloquetModeUse::Dtn))
    {
      mfem::DenseMatrix ee_mat(3);
      for (int a = 0; a < 3; a++)
      {
        for (int b = 0; b < 3; b++)
        {
          ee_mat(a, b) = mode.e_pol(a) * mode.e_pol(b);
        }
      }
      mfem::MatrixConstantCoefficient ee_coeff(ee_mat);

      auto *integ = new mfem::VectorFEMassIntegrator(ee_coeff);
      if (ir)
      {
        integ->SetIntRule(ir);
      }
      mfem::ParBilinearForm bf(&nd_fespace);
      bf.AddBoundaryIntegrator(integ, bdr_marker);
      bf.Assemble();
      bf.Finalize();
      mode.A_mass.reset(bf.ParallelAssemble());

      // Diagnostic: check if A_mass has non-zero entries by applying to v.
      {
        Vector Av(tdof_size);
        Av = 0.0;
        mode.A_mass->Mult(mode.v.Real(), Av);
        double Av_norm = std::sqrt(linalg::Dot(nd_fespace.GetComm(), Av, Av));
        double v_norm =
            std::sqrt(linalg::Dot(nd_fespace.GetComm(), mode.v.Real(), mode.v.Real()));
        Mpi::Print("    A_mass({:+d},{:+d},{}) ||Av||={:.4e} ||v||={:.4e}\n", mode.m,
                   mode.n, mode.is_te ? "TE" : "TM", Av_norm, v_norm);
      }
    }
  }
}

namespace
{

void ComputePolarization(const mfem::Vector &kt, const mfem::Vector &port_normal,
                         mfem::Vector &e_te, mfem::Vector &e_tm)
{
  double kt_norm = kt.Norml2();
  e_te.SetSize(3);
  e_tm.SetSize(3);

  if (kt_norm > 1e-12)
  {
    mfem::Vector kt_hat(3);
    kt_hat = kt;
    kt_hat *= 1.0 / kt_norm;

    e_te(0) = kt_hat(1) * port_normal(2) - kt_hat(2) * port_normal(1);
    e_te(1) = kt_hat(2) * port_normal(0) - kt_hat(0) * port_normal(2);
    e_te(2) = kt_hat(0) * port_normal(1) - kt_hat(1) * port_normal(0);
    double e_te_norm = e_te.Norml2();
    if (e_te_norm > 1e-12)
    {
      e_te *= 1.0 / e_te_norm;
    }

    e_tm(0) = port_normal(1) * e_te(2) - port_normal(2) * e_te(1);
    e_tm(1) = port_normal(2) * e_te(0) - port_normal(0) * e_te(2);
    e_tm(2) = port_normal(0) * e_te(1) - port_normal(1) * e_te(0);
    double e_tm_norm = e_tm.Norml2();
    if (e_tm_norm > 1e-12)
    {
      e_tm *= 1.0 / e_tm_norm;
    }
  }
  else
  {
    // Normal incidence: degenerate TE/TM. Use consistent reference direction.
    mfem::Vector abs_n(3);
    for (int d = 0; d < 3; d++)
    {
      abs_n(d) = std::abs(port_normal(d));
    }
    mfem::Vector ref(3);
    ref = 0.0;
    int min_idx = 0;
    double min_val = abs_n(0);
    for (int d = 1; d < 3; d++)
    {
      if (abs_n(d) < min_val)
      {
        min_val = abs_n(d);
        min_idx = d;
      }
    }
    ref(min_idx) = 1.0;

    e_te(0) = ref(1) * abs_n(2) - ref(2) * abs_n(1);
    e_te(1) = ref(2) * abs_n(0) - ref(0) * abs_n(2);
    e_te(2) = ref(0) * abs_n(1) - ref(1) * abs_n(0);
    double e_te_norm = e_te.Norml2();
    e_te *= 1.0 / e_te_norm;

    e_tm(0) = abs_n(1) * e_te(2) - abs_n(2) * e_te(1);
    e_tm(1) = abs_n(2) * e_te(0) - abs_n(0) * e_te(2);
    e_tm(2) = abs_n(0) * e_te(1) - abs_n(1) * e_te(0);
    double e_tm_norm = e_tm.Norml2();
    e_tm *= 1.0 / e_tm_norm;
  }
}

}  // namespace

void FloquetPortData::Initialize(double omega)
{
  if (omega == omega0 && omega0 != 0.0)
  {
    return;
  }
  omega0 = omega;

  // Compute k_F at this frequency. When frequency scaling is active, the stored k_F is
  // k₀ = k_F/ω (frequency-independent), so k_F(ω) = ω × k₀. Without scaling, k_F is the
  // fixed Bloch vector and kF_scale = 1.
  double kF_scale = mat_op.HasFloquetFrequencyScaling() ? omega : 1.0;

  // Check BZ wrapping consistency: if k_F(ω) would require a different BZ offset than
  // what was used at construction, the mode enumeration is wrong. Compute the effective
  // k_F(ω) and check its BZ offset.
  if (mat_op.HasFloquetFrequencyScaling())
  {
    mfem::Vector kF_eff(3);
    for (int i = 0; i < 3; i++)
    {
      kF_eff(i) = kF_scale * k_F(i);
    }
    // Wrap to BZ and compute offset.
    double b1_sq = b1 * b1, b2_sq = b2 * b2;
    int new_bz_m = static_cast<int>(std::round((kF_eff * b1) / b1_sq));
    int new_bz_n = static_cast<int>(std::round((kF_eff * b2) / b2_sq));
    // The constructor computed bz_m/bz_n from the UNWRAPPED config k_F. Here we check
    // the effective k_F(ω) — the offset should match after accounting for the reference.
    // Compare with the BZ offset that was used at construction.
    MFEM_VERIFY(new_bz_m == bz_m && new_bz_n == bz_n,
                "Floquet port: k_F("
                    << omega << ") requires BZ offset (" << new_bz_m << "," << new_bz_n
                    << ") different from construction (" << bz_m << "," << bz_n
                    << "). The frequency range may be too wide for the reference "
                       "frequency. Consider using a reference frequency closer to the "
                       "center of the sweep range.");
  }

  bool need_reassemble = false;
  for (auto &mode : modes)
  {
    // gamma_mn^2 = omega^2 * mu_r * eps_r - |B_mn - k_F(ω)|^2
    mfem::Vector kt(3);
    for (int i = 0; i < 3; i++)
    {
      kt(i) = mode.B_mn(i) - kF_scale * k_F(i);
    }
    double kt_sq = kt * kt;
    mode.gamma_sq = omega * omega * mu_eps_port - kt_sq;

    // When k_F scales with frequency, the transverse wavevector direction changes,
    // which rotates the TE/TM polarization vectors. Recompute them.
    if (mat_op.HasFloquetFrequencyScaling())
    {
      mfem::Vector e_te(3), e_tm(3);
      ComputePolarization(kt, port_normal, e_te, e_tm);
      mfem::Vector &new_pol = mode.is_te ? e_te : e_tm;
      // Check if polarization changed (avoid unnecessary reassembly).
      double pol_diff = 0.0;
      for (int d = 0; d < 3; d++)
      {
        pol_diff += std::abs(mode.e_pol(d) - new_pol(d));
      }
      if (pol_diff > 1e-14)
      {
        mode.e_pol = new_pol;
        need_reassemble = true;
      }
    }
  }

  // Re-assemble mode vectors if polarization changed.
  if (need_reassemble && nd_fespace_ptr)
  {
    AssembleFourierProjections(*nd_fespace_ptr);
  }
}

int FloquetPortData::NumPropagatingOrders() const
{
  int count = 0;
  for (const auto &mode : modes)
  {
    if (mode.gamma_sq > 0.0)
    {
      count++;
    }
  }
  return count;
}

std::complex<double>
FloquetPortData::ComputeDtNCorrectionCoeff(const FloquetMode &mode) const
{
  // Compute g_correction = g_full - g_uniform for a single mode.
  // g_uniform uses the TE (0,0) propagation constant as the Robin reference.
  double gamma0 = 0.0;
  for (const auto &m0 : modes)
  {
    if (m0.m == 0 && m0.n == 0 && m0.is_te)
    {
      gamma0 = std::sqrt(std::max(m0.gamma_sq, 0.0));
      break;
    }
  }

  std::complex<double> g_uniform = 1i * gamma0 / (mu_r_port * port_area);
  std::complex<double> g_full = ComputeDtNFullCoeff(mode);
  if (g_full == 0.0)
  {
    return 0.0;
  }

  std::complex<double> g_correction = g_full - g_uniform;

  // Skip negligible corrections.
  if (std::abs(g_correction) < 1e-14 * std::abs(g_full))
  {
    return 0.0;
  }

  return g_correction;
}

std::complex<double> FloquetPortData::ComputeDtNFullCoeff(const FloquetMode &mode) const
{
  // Return the full DtN coefficient g_full (NOT the correction g_full - g_uniform).
  if (mode.gamma_sq > 0.0)
  {
    double gamma = std::sqrt(mode.gamma_sq);
    return mode.is_te
               ? 1i * gamma / (mu_r_port * port_area)
               : 1i * omega0 * omega0 * mu_eps_port / (gamma * mu_r_port * port_area);
  }
  else if (mode.gamma_sq < 0.0)
  {
    double gamma_abs = std::sqrt(-mode.gamma_sq);
    return mode.is_te
               ? gamma_abs / (mu_r_port * port_area)
               : -omega0 * omega0 * mu_eps_port / (gamma_abs * mu_r_port * port_area);
  }
  return 0.0;
}

std::unique_ptr<ComplexOperator> FloquetPortData::GetBoundaryOperator() const
{
  int n = modes.empty() ? 0 : static_cast<int>(modes[0].v.Size());

  // DtN coefficient diagnostic: print propagation constants for all modes.
  // beta = gamma / Lc [rad/m] for propagating, beta = j*gamma_abs / Lc for evanescent.
  // Lc = characteristic length. To convert nondim gamma to physical beta [rad/m]:
  // beta = gamma_nondim / Lc. But we don't have Lc here, so print gamma_nondim
  // (= beta * Lc) and |kt| (= |B - kF|) for COMSOL comparison.
  {
    Mpi::Print(" Floquet port mode propagation constants at omega={:.6e}:\n", omega0);
    for (const auto &mode : modes)
    {
      if (!HasFlag(mode.use, FloquetModeUse::Dtn))
      {
        continue;
      }
      // kt = B - kF(omega) for this mode.
      double kF_scale = mat_op.HasFloquetFrequencyScaling() ? omega0 : 1.0;
      mfem::Vector kt(3);
      for (int d = 0; d < 3; d++)
      {
        kt(d) = mode.B_mn(d) - kF_scale * k_F(d);
      }
      double kt_norm = kt.Norml2();

      // gamma = sqrt(gamma_sq) for propagating, j*|gamma| for evanescent.
      // COMSOL's beta [rad/m] = gamma_nondim / Lc.
      double gamma_re = 0.0, gamma_im = 0.0;
      if (mode.gamma_sq > 0.0)
      {
        gamma_re = std::sqrt(mode.gamma_sq);  // propagating: real beta
      }
      else if (mode.gamma_sq < 0.0)
      {
        gamma_im = -std::sqrt(-mode.gamma_sq);  // evanescent: negative imag (Im(k)<=0)
      }
      auto g_full = ComputeDtNFullCoeff(mode);
      Mpi::Print("  ({:+d},{:+d},{}) gamma_sq={:+.6e} gamma=({:+.6e},{:+.6e}) "
                 "|kt|={:.6e} g_full=({:+.6e},{:+.6e})\n",
                 mode.m, mode.n, mode.is_te ? "TE" : "TM", mode.gamma_sq, gamma_re,
                 gamma_im, kt_norm, g_full.real(), g_full.imag());
    }
  }

  if (use_mass_consistent)
  {
    // Mass-consistent correction: Robin stays in the system for baseline absorption.
    // The correction uses per-mode boundary mass matrices A_i instead of rank-1 conj(v)v^T.
    // Full complex correction (no real-only truncation for evanescent modes) — the
    // mass-consistent A_i should cancel Robin precisely for each mode's polarization
    // direction, avoiding the Cauchy-Schwarz overestimate of rank-1 outer products.
    //
    // Coefficient: g_mc = (g_full - g_uniform) × |Γ| = j(λ_i - γ₀)/μ
    // This removes the 1/|Γ| from the rank-1 coefficient to match A_i's area scaling.
    double gamma0 = 0.0;
    for (const auto &m0 : modes)
    {
      if (m0.m == 0 && m0.n == 0 && m0.is_te)
      {
        gamma0 = std::sqrt(std::max(m0.gamma_sq, 0.0));
        break;
      }
    }
    auto op = std::make_unique<MassConsistentDtNOperator>(n);
    int n_added = 0;
    for (const auto &mode : modes)
    {
      if (!HasFlag(mode.use, FloquetModeUse::Dtn) || !mode.A_mass)
      {
        continue;
      }
      // Full complex correction: g_mc = g_full_mc - Robin_share_per_mode.
      // Robin adds iγ₀/μ × M_bdr to the system. Since Σ_all_modes A_i = N_orders × M_bdr
      // (because A_TE + A_TM = M_bdr for each order), each mode's share of Robin is
      // iγ₀/(μ × N_orders). The correction subtracts this per-mode share.
      int n_orders = static_cast<int>(modes.size()) / 2;  // TE+TM per order
      auto g_mc =
          ComputeDtNFullCoeff(mode) * port_area - 1i * gamma0 / (mu_r_port * n_orders);
      if (std::abs(g_mc) > 1e-14 * std::abs(ComputeDtNFullCoeff(mode) * port_area))
      {
        op->AddTerm(mode.A_mass.get(), g_mc);
        n_added++;
      }
    }
    Mpi::Print(" >> MassConsistent DtN (Robin+full complex): {:d} terms\n", n_added);
    return op;
  }

  // Default: LowRankComplexOperator (rank-1 outer products).
  auto op = std::make_unique<LowRankComplexOperator>(comm, n);
  for (const auto &mode : modes)
  {
    if (!HasFlag(mode.use, FloquetModeUse::Dtn))
    {
      continue;
    }
    auto g = use_full_dtn ? ComputeDtNFullCoeff(mode) : ComputeDtNCorrectionCoeff(mode);
    if (g != 0.0)
    {
      op->AddTerm(&mode.v, g);
    }
  }
  return op;
}

FloquetPortData::IncidentNormalization
FloquetPortData::ComputeIncidentNormalization(double omega) const
{
  IncidentNormalization norm{};
  // TE and TM of the same (0,0) order share the same gamma_sq, so it doesn't matter
  // which polarization we find first.
  for (const auto &mode : modes)
  {
    if (mode.m == 0 && mode.n == 0 && mode.gamma_sq > 0.0)
    {
      norm.gamma_00 = std::sqrt(mode.gamma_sq);
      break;
    }
  }
  MFEM_VERIFY(norm.gamma_00 > 0.0, "Incident Floquet mode is evanescent or not found!");
  norm.lambda_te_00 = norm.gamma_00;
  norm.lambda_tm_00 = omega * omega * mu_eps_port / norm.gamma_00;
  norm.lambda_eff = std::norm(inc_alpha_te) * norm.lambda_te_00 +
                    std::norm(inc_alpha_tm) * norm.lambda_tm_00;
  double p_unit = norm.lambda_eff * port_area / (2.0 * omega * mu_r_port);
  norm.c_inc = 1.0 / std::sqrt(p_unit);
  return norm;
}

std::map<std::tuple<int, int, bool>, std::complex<double>>
FloquetPortData::GetAllSParameters(const GridFunction &E, bool subtract_incident) const
{
  std::map<std::tuple<int, int, bool>, std::complex<double>> result;

  // Power-normalized S-parameter extraction.
  //
  // Step 1: Field amplitude c_mn = v_mn^T E / |Γ| (Fourier coefficient of total field).
  // Step 2: Power normalization: S_mn = √(γ_mn / γ_inc) × c_mn.
  //
  // This ensures Σ |S_mn|² = 1 for energy conservation (lossless). The factor √(γ_mn/γ_inc)
  // accounts for different power-per-amplitude for different diffraction orders
  // (higher-angle orders carry less power per unit amplitude due to oblique propagation).
  //
  // For the incident (0,0) mode: γ_inc = γ_00 = ω √(μ_r ε_r), so √(γ/γ_inc) = 1.
  // For the driving port: S = √(γ/γ_inc) × c - δ_{incident mode} (subtract 1).

  auto norm = ComputeIncidentNormalization(omega0);
  double lambda_te_00 = norm.lambda_te_00;
  double lambda_tm_00 = norm.lambda_tm_00;
  double lambda_eff = norm.lambda_eff;
  double c_inc = norm.c_inc;

  // Restrict E to true DOFs once (shared across all modes).
  const auto *P = E.Real().ParFESpace()->GetProlongationMatrix();
  int tdof_size = modes.empty() ? 0 : static_cast<int>(modes[0].v.Size());
  Vector E_r_tdof(tdof_size), E_i_tdof(tdof_size);
  P->MultTranspose(E.Real(), E_r_tdof);
  if (E.HasImag())
  {
    P->MultTranspose(E.Imag(), E_i_tdof);
  }
  else
  {
    E_i_tdof = 0.0;
  }

  for (const auto &mode : modes)
  {
    if (!HasFlag(mode.use, FloquetModeUse::Output) || mode.gamma_sq <= 0.0)
    {
      continue;  // Only user-requested propagating orders carry S-parameters.
    }

    double gamma_mn = std::sqrt(mode.gamma_sq);

    // Bilinear v^T E = (v_r·E_r - v_i·E_i) + i(v_r·E_i + v_i·E_r).
    double sr = linalg::Dot(comm, mode.v.Real(), E_r_tdof) -
                linalg::Dot(comm, mode.v.Imag(), E_i_tdof);
    double si = linalg::Dot(comm, mode.v.Real(), E_i_tdof) +
                linalg::Dot(comm, mode.v.Imag(), E_r_tdof);

    // Field amplitude: c = v^T E / (c_inc × |Γ|), where c_inc accounts for the
    // unit-power normalization applied to the excitation.
    // Power normalization: S = √(λ_mn / λ_eff) × c, where λ is the DtN eigenvalue.
    // TE: λ = γ. TM: λ = ω²με/γ. λ_eff is the weighted incident eigenvalue.
    double lambda_mn = mode.is_te ? gamma_mn : omega0 * omega0 * mu_eps_port / gamma_mn;
    double power_factor = std::sqrt(lambda_mn / lambda_eff);

    auto key = std::make_tuple(mode.m, mode.n, mode.is_te);
    result[key] = power_factor * std::complex<double>(sr, si) / (c_inc * port_area);
  }

  // Subtract incident field contribution for the driving port.
  // S_{0,0,p} -= √(λ_p / λ_eff) × α_p for each incident polarization component.
  if (subtract_incident)
  {
    for (auto &[key, S] : result)
    {
      auto [m, n, is_te] = key;
      if (IsIncidentMode(m, n, is_te))
      {
        double lambda_p = is_te ? lambda_te_00 : lambda_tm_00;
        S -= std::sqrt(lambda_p / lambda_eff) * GetIncidentAlpha(is_te);
      }
    }
  }

  return result;
}

bool FloquetPortData::AddExcitationVector(double omega, ComplexVector &RHS) const
{
  // Collect the (0,0) TE and TM modes for the incident excitation.
  const FloquetMode *mode_te = nullptr, *mode_tm = nullptr;
  for (const auto &mode : modes)
  {
    if (mode.m == 0 && mode.n == 0)
    {
      if (mode.is_te)
      {
        mode_te = &mode;
      }
      else
      {
        mode_tm = &mode;
      }
    }
  }
  if (!mode_te && !mode_tm)
  {
    return false;
  }

  auto norm = ComputeIncidentNormalization(omega);
  double lambda_te = norm.lambda_te_00;
  double lambda_tm = norm.lambda_tm_00;
  double c_inc = norm.c_inc;

  // f = Σ_p c_inc × 2i × α_p × λ_p / μ_r × conj(v_p)
  // conj(v) appears because the bilinear form uses v^T (not v^H), and the incident
  // field ê exp(+iB·r) projects to conj(v).
  auto add_pol = [&](const FloquetMode *mode, std::complex<double> alpha, double lambda)
  {
    if (!mode || std::abs(alpha) < 1e-14)
    {
      return;
    }
    // s = c_inc × 2i × alpha × lambda / mu_r (complex scalar)
    std::complex<double> s = c_inc * 2.0 * 1i * alpha * lambda / mu_r_port;
    // RHS += s × conj(v) = (s_r + i s_i)(v_r - i v_i)
    RHS.Real().Add(s.real(), mode->v.Real());
    RHS.Real().Add(s.imag(), mode->v.Imag());
    RHS.Imag().Add(s.imag(), mode->v.Real());
    RHS.Imag().Add(-s.real(), mode->v.Imag());
  };
  add_pol(mode_te, inc_alpha_te, lambda_te);
  add_pol(mode_tm, inc_alpha_tm, lambda_tm);
  return true;
}

// =============================================================================
// FloquetPortOperator
// =============================================================================

FloquetPortOperator::FloquetPortOperator(const IoData &iodata,
                                         const MaterialOperator &mat_op_ref,
                                         mfem::ParFiniteElementSpace &nd_fespace)
  : mat_op(mat_op_ref)
{
  const auto &floquet_port_data = iodata.boundaries.floquetport;
  if (floquet_port_data.empty())
  {
    return;
  }
  MFEM_VERIFY(iodata.problem.type == ProblemType::DRIVEN,
              "Floquet port boundaries are only available for frequency domain driven "
              "simulations!");

  Mpi::Print("\nConfiguring {:d} Floquet port boundar{}\n", floquet_port_data.size(),
             (floquet_port_data.size() > 1) ? "ies" : "y");

  for (const auto &[idx, data] : floquet_port_data)
  {
    ports.try_emplace(idx, data, iodata, mat_op_ref, nd_fespace);
  }
}

void FloquetPortOperator::Initialize(double omega)
{
  for (auto &[idx, port] : ports)
  {
    port.Initialize(omega);
  }
}

std::unique_ptr<ComplexOperator> FloquetPortOperator::GetExtraSystemOperator(double omega)
{
  if (ports.empty())
  {
    return nullptr;
  }

  Initialize(omega);

  // Combine boundary operators from all ports.
  // For a single port, return it directly. For multiple, chain with SumComplexOperator.
  // Auxiliary ports skip the low-rank correction (DtN is handled by the coupled system).
  std::unique_ptr<ComplexOperator> combined;
  for (auto &[idx, port] : ports)
  {
    if (!port.active || port.UseAuxiliary())
    {
      continue;
    }
    auto F = port.GetBoundaryOperator();
    if (!combined)
    {
      combined = std::move(F);
    }
    else
    {
      combined = std::make_unique<SumComplexOperator>(std::move(combined), std::move(F));
    }
  }

  return combined;
}

void FloquetPortOperator::AddExtraSystemBdrCoefficients(double omega,
                                                        MaterialPropertyCoefficient &fbr,
                                                        MaterialPropertyCoefficient &fbi,
                                                        bool for_preconditioner)
{
  // Add a full-rank Robin BC (iγ₀/μ) on the Floquet port boundary faces.
  // For FullDtN ports: skip Robin in the system matrix (the dense DtN operator handles
  // absorption). Always include Robin in the preconditioner for the inner solve.
  Initialize(omega);
  for (const auto &[idx, port] : ports)
  {
    if (!port.active)
    {
      continue;
    }
    if ((port.UseFullDtN() || port.UseAuxiliary()) && !for_preconditioner)
    {
      continue;  // FullDtN/Auxiliary: no Robin in system, DtN operator handles everything.
    }
    // γ₀ = TE (0,0) mode propagation constant.
    double gamma0 = 0.0;
    for (const auto &mode : port.GetModes())
    {
      if (mode.m == 0 && mode.n == 0 && mode.is_te)
      {
        gamma0 = std::sqrt(std::max(mode.gamma_sq, 0.0));
        break;
      }
    }
    MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    muinv_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(port.GetAttrList()));
    fbi.AddCoefficient(muinv_func.GetAttributeToMaterial(),
                       muinv_func.GetMaterialProperties(), gamma0);
  }
}

bool FloquetPortOperator::AddExcitationVector(int excitation_idx, double omega,
                                              ComplexVector &RHS)
{
  Initialize(omega);
  bool nnz = false;
  for (auto &[idx, port] : ports)
  {
    if (port.excitation != excitation_idx)
    {
      continue;
    }
    nnz |= port.AddExcitationVector(omega, RHS);
  }
  return nnz;
}

mfem::Array<int> FloquetPortOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, port] : ports)
  {
    attr_list.Append(port.GetAttrList());
  }
  return attr_list;
}

bool FloquetPortOperator::HasAuxiliaryPorts() const
{
  for (const auto &[idx, port] : ports)
  {
    if (port.active && port.UseAuxiliary())
    {
      return true;
    }
  }
  return false;
}

int FloquetPortOperator::GetAuxModeCount(double omega)
{
  Initialize(omega);
  int count = 0;
  for (const auto &[idx, port] : ports)
  {
    if (!port.active || !port.UseAuxiliary())
    {
      continue;
    }
    for (const auto &mode : port.GetModes())
    {
      if (HasFlag(mode.use, FloquetModeUse::Dtn))
      {
        count++;
      }
    }
  }
  return count;
}

void FloquetPortOperator::PopulateAuxModes(double omega, FloquetAuxSystemOperator &aux_op)
{
  Initialize(omega);
  aux_op.ClearModes();
  for (auto &[idx, port] : ports)
  {
    if (!port.active || !port.UseAuxiliary())
    {
      continue;
    }
    for (const auto &mode : port.GetModes())
    {
      if (!HasFlag(mode.use, FloquetModeUse::Dtn))
      {
        continue;
      }
      aux_op.AddMode(&mode.v, port.ComputeDtNFullCoeff(mode));
    }
  }
}

// =============================================================================
// FloquetAuxSystemOperator
// =============================================================================

FloquetAuxSystemOperator::FloquetAuxSystemOperator(MPI_Comm comm,
                                                   const ComplexOperator *A_EE, int n_E,
                                                   int n_s)
  : ComplexOperator(n_E + n_s, n_E + n_s), A_EE(A_EE), comm(comm), n_E(n_E), n_s(n_s),
    E_tmp(n_E), y_E_tmp(n_E)
{
  E_tmp.UseDevice(true);
  y_E_tmp.UseDevice(true);
}

void FloquetAuxSystemOperator::AddMode(const ComplexVector *v, std::complex<double> g_full)
{
  aux_modes.push_back({v, g_full});
}

void FloquetAuxSystemOperator::Mult(const ComplexVector &x, ComplexVector &y) const
{
  y = 0.0;
  AddMult(x, y, 1.0);
}

void FloquetAuxSystemOperator::AddMult(const ComplexVector &x, ComplexVector &y,
                                       const std::complex<double> a) const
{
  const int n_modes = static_cast<int>(aux_modes.size());
  int rank;
  MPI_Comm_rank(comm, &rank);

  // Extract x_E (first n_E entries) into E_tmp.
  {
    const auto *xr = x.Real().Read();
    const auto *xi = x.Imag().Read();
    auto *er = E_tmp.Real().Write();
    auto *ei = E_tmp.Imag().Write();
    std::copy_n(xr, n_E, er);
    std::copy_n(xi, n_E, ei);
  }

  // Read s values from rank 0 and broadcast to all ranks.
  std::vector<double> s_r(n_modes), s_i(n_modes);
  if (rank == 0)
  {
    const auto *xr = x.Real().Read();
    const auto *xi = x.Imag().Read();
    for (int i = 0; i < n_modes; i++)
    {
      s_r[i] = xr[n_E + i];
      s_i[i] = xi[n_E + i];
    }
  }
  MPI_Bcast(s_r.data(), n_modes, MPI_DOUBLE, 0, comm);
  MPI_Bcast(s_i.data(), n_modes, MPI_DOUBLE, 0, comm);

  // E block: y_E += a * (A_EE * x_E + Σ_i g_full_i * s_i * conj(v_i))

  // Apply A_EE * E_tmp → y_E_tmp, then add to y[0:n_E].
  A_EE->Mult(E_tmp, y_E_tmp);
  {
    const auto *yr = y_E_tmp.Real().Read();
    const auto *yi = y_E_tmp.Imag().Read();
    auto *out_r = y.Real().ReadWrite();
    auto *out_i = y.Imag().ReadWrite();
    double ar = a.real(), ai = a.imag();
    for (int j = 0; j < n_E; j++)
    {
      out_r[j] += ar * yr[j] - ai * yi[j];
      out_i[j] += ai * yr[j] + ar * yi[j];
    }
  }

  // A_ES coupling: y_E += a * Σ_i g_full_i * s_i * conj(v_i)
  for (int i = 0; i < n_modes; i++)
  {
    // Complex product: coeff = a * g_full * s  (all complex)
    std::complex<double> s_val(s_r[i], s_i[i]);
    std::complex<double> coeff = a * aux_modes[i].g_full * s_val;

    // Add coeff * conj(v) to y[0:n_E].
    // conj(v) = (v_r, -v_i), so coeff * conj(v) = (coeff_r * v_r + coeff_i * v_i,
    //                                               coeff_i * v_r - coeff_r * v_i)
    double cr = coeff.real(), ci = coeff.imag();
    const auto *vr = aux_modes[i].v->Real().Read();
    const auto *vi = aux_modes[i].v->Imag().Read();
    auto *out_r = y.Real().ReadWrite();
    auto *out_i = y.Imag().ReadWrite();
    for (int j = 0; j < n_E; j++)
    {
      out_r[j] += cr * vr[j] + ci * vi[j];
      out_i[j] += ci * vr[j] - cr * vi[j];
    }
  }

  // s block: y_s[i] += a * (v_i^T * x_E - s_i)  (rank 0 only)
  if (rank == 0)
  {
    auto *out_r = y.Real().ReadWrite();
    auto *out_i = y.Imag().ReadWrite();
    for (int i = 0; i < n_modes; i++)
    {
      // Bilinear dot product: dot = v^T x = Σ v_j x_j (no conjugate)
      double dot_r = linalg::Dot(comm, aux_modes[i].v->Real(), E_tmp.Real()) -
                     linalg::Dot(comm, aux_modes[i].v->Imag(), E_tmp.Imag());
      double dot_i = linalg::Dot(comm, aux_modes[i].v->Real(), E_tmp.Imag()) +
                     linalg::Dot(comm, aux_modes[i].v->Imag(), E_tmp.Real());

      // val = dot - s
      double val_r = dot_r - s_r[i];
      double val_i = dot_i - s_i[i];

      // y_s[i] += a * val
      out_r[n_E + i] += a.real() * val_r - a.imag() * val_i;
      out_i[n_E + i] += a.imag() * val_r + a.real() * val_i;
    }
  }
  else
  {
    // Non-rank-0 ranks still need to participate in the MPI_Allreduce inside linalg::Dot,
    // but discard the result (s entries stay zero on these ranks).
    for (int i = 0; i < n_modes; i++)
    {
      linalg::Dot(comm, aux_modes[i].v->Real(), E_tmp.Real());
      linalg::Dot(comm, aux_modes[i].v->Imag(), E_tmp.Imag());
      linalg::Dot(comm, aux_modes[i].v->Real(), E_tmp.Imag());
      linalg::Dot(comm, aux_modes[i].v->Imag(), E_tmp.Real());
    }
  }
}

}  // namespace palace
