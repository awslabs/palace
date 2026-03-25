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
  // The Floquet DtN operator: F = Σ g_k v_k v_k^H / |Γ| (with 1/|Γ| absorbed into g_k).
  for (const auto &t : terms)
  {
    // Hermitian dot product: dot = v^H x = Σ conj(v_j) x_j
    // = (v_r · x_r + v_i · x_i) + i(v_r · x_i - v_i · x_r)
    double dot_r =
        linalg::Dot(comm, t.v->Real(), x.Real()) + linalg::Dot(comm, t.v->Imag(), x.Imag());
    double dot_i =
        linalg::Dot(comm, t.v->Real(), x.Imag()) - linalg::Dot(comm, t.v->Imag(), x.Real());

    // Compute s = a * g * dot (complex scalar)
    std::complex<double> s = a * t.g * std::complex<double>(dot_r, dot_i);

    // y += s * v = (s_r + i s_i)(v_r + i v_i)
    //    Real part: s_r v_r - s_i v_i
    //    Imag part: s_i v_r + s_r v_i
    y.Real().Add(s.real(), t.v->Real());
    y.Real().Add(-s.imag(), t.v->Imag());
    y.Imag().Add(s.imag(), t.v->Real());
    y.Imag().Add(s.real(), t.v->Imag());
  }
}

namespace
{

void ComputePolarization(const mfem::Vector &kt, const mfem::Vector &port_normal,
                         mfem::Vector &e_te, mfem::Vector &e_tm);

}  // namespace

// =============================================================================
// FloquetPortData
// =============================================================================

FloquetPortData::FloquetPortData(const config::FloquetPortData &data,
                                 const config::PeriodicBoundaryData &periodic_data,
                                 const Units &units, const MaterialOperator &mat_op_ref,
                                 mfem::ParFiniteElementSpace &nd_fespace)
  : excitation(data.excitation), mat_op(mat_op_ref), a1(3), a2(3), b1(3), b2(3), k_F(3),
    port_normal(3), comm(nd_fespace.GetComm())
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
  const auto &periodic = periodic_data;
  MFEM_VERIFY(!periodic.boundary_pairs.empty(),
              "FloquetPort requires periodic boundary conditions to be configured under "
              "\"Boundaries\"/\"Periodic\"/\"BoundaryPairs\". At least two periodic "
              "boundary pairs are needed for 3D periodicity.");

  // Extract lattice vectors. If periodic boundary pairs with translations are specified,
  // use them. Otherwise, infer from the port face bounding box (works for axis-aligned
  // rectangular cells with built-in mesh periodicity).
  double mesh_scale = units.GetMeshLengthRelativeScale();
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
  }

  // Bloch wave vector: use the BZ-wrapped value from MaterialOperator so the port is
  // consistent with the bulk ∇_F bilinear form. The Floquet port must use the same BZ
  // convention as the volume solve for correct Fourier extraction (v^T E).
  const auto &kF_mat = mat_op.GetWaveVector();
  for (int i = 0; i < 3; i++)
  {
    k_F(i) = (i < kF_mat.Size()) ? kF_mat(i) : 0.0;
  }

  // Compute reciprocal lattice (MUST be before BZ offset which uses b1, b2).
  ComputeReciprocalLattice(a1, a2, b1, b2);

  // Store the unwrapped config kF (Lc units) for frequency-dependent BZ offset.
  kF_config_Lc.SetSize(3);
  for (int i = 0; i < 3; i++)
  {
    kF_config_Lc(i) =
        (i < static_cast<int>(periodic.wave_vector.size())) ? periodic.wave_vector[i] : 0.0;
  }

  // Compute the BZ wrapping offset at the reference frequency.
  {
    const auto &kF_bz = mat_op.GetWaveVectorBZ();
    mfem::Vector kF_wrapped(3);
    for (int i = 0; i < 3; i++)
    {
      kF_wrapped(i) = (i < kF_bz.Size()) ? kF_bz(i) : 0.0;
    }
    bz_m = ComputeBZOffset(kF_config_Lc, kF_wrapped, b1, b1 * b1);
    bz_n = ComputeBZOffset(kF_config_Lc, kF_wrapped, b2, b2 * b2);
    if (bz_m != 0 || bz_n != 0)
    {
      Mpi::Print(" Floquet port: BZ wrapping at ref freq (shift = {:d}*b1 + {:d}*b2)\n",
                 bz_m, bz_n);
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
      muinv_local = mat_op.GetInvPermeability(port_vol_attr)(0, 0);
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

  // Cap MaxOrder at the mesh Nyquist limit: p-th order elements can resolve Fourier
  // modes with |B|×h < p×π. Beyond this, projections alias to ~0 with default quadrature.
  {
    auto &mesh = *nd_fespace.GetParMesh();
    int bdr_attr_max = 0;
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      bdr_attr_max = std::max(bdr_attr_max, mesh.GetBdrAttribute(be));
    }
    Mpi::GlobalMax(1, &bdr_attr_max, nd_fespace.GetComm());
    auto bdr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list, true);

    double h_max = 0.0;
    for (int be = 0; be < mesh.GetNBE(); be++)
    {
      if (bdr_marker[mesh.GetBdrAttribute(be) - 1])
      {
        auto *T = mesh.GetBdrElementTransformation(be);
        T->SetIntPoint(&mfem::Geometries.GetCenter(T->GetGeometryType()));
        h_max = std::max(h_max, T->Jacobian().CalcSingularvalue(0));
      }
    }
    Mpi::GlobalMax(1, &h_max, nd_fespace.GetComm());

    if (h_max > 0.0)
    {
      int p = nd_fespace.GetMaxElementOrder();
      int nyquist_m =
          std::max(1, static_cast<int>(std::floor(p * M_PI / (b1.Norml2() * h_max))));
      int nyquist_n =
          std::max(1, static_cast<int>(std::floor(p * M_PI / (b2.Norml2() * h_max))));
      if (max_order_m > nyquist_m || max_order_n > nyquist_n)
      {
        Mpi::Print(" Floquet port: capping MaxOrder from ({:d}, {:d}) to ({:d}, {:d}) "
                   "(Nyquist limit for h_max={:.4e}, p={:d})\n",
                   max_order_m, max_order_n, std::min(max_order_m, nyquist_m),
                   std::min(max_order_n, nyquist_n), h_max, p);
        max_order_m = std::min(max_order_m, nyquist_m);
        max_order_n = std::min(max_order_n, nyquist_n);
      }
    }
  }

  // Enumerate diffraction orders and compute mode vectors.
  EnumerateOrders();
  AssembleFourierProjections(nd_fespace);
  nd_fespace_ptr = &nd_fespace;  // Store for re-assembly in Initialize.

  Mpi::Print(" Floquet port: {:d} modes ({:d} orders x 2 polarizations)\n",
             static_cast<int>(modes.size()), static_cast<int>(modes.size()) / 2);
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
          kt(i) = B_mn(i) + ((i < kF_phys.Size()) ? kF_phys(i) : 0.0);
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

    // Assemble real part: v_r = ∫_Γ N_j · [e_pol cos(B·r)] dS.
    mfem::LinearForm lf_r(&nd_fespace);
    lf_r.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(coeff_r), bdr_marker);
    lf_r.UseFastAssembly(false);
    lf_r.UseDevice(false);
    lf_r.Assemble();
    lf_r.UseDevice(true);

    // Assemble imaginary part: v_i = ∫_Γ N_j · [e_pol (-sin(B·r))] dS.
    mfem::LinearForm lf_i(&nd_fespace);
    lf_i.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(coeff_i), bdr_marker);
    lf_i.UseFastAssembly(false);
    lf_i.UseDevice(false);
    lf_i.Assemble();
    lf_i.UseDevice(true);

    // Project to true DOFs: v = P^T * lf.
    nd_fespace.GetProlongationMatrix()->MultTranspose(lf_r, mode.v.Real());
    nd_fespace.GetProlongationMatrix()->MultTranspose(lf_i, mode.v.Imag());
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

  // BZ consistency: when frequency scaling is active, BZ wrapping is disabled in
  // MaterialOperator (the volume and port both use the unwrapped kF). No zone
  // transitions are needed. For fixed kF (no frequency scaling), BZ wrapping is
  // handled at construction and stays fixed.
  // (No runtime BZ check needed.)

  bool need_reassemble = false;
  for (auto &mode : modes)
  {
    // gamma_mn^2 = omega^2 * mu_r * eps_r - |B_mn + k_F(ω)|^2
    mfem::Vector kt(3);
    for (int i = 0; i < 3; i++)
    {
      kt(i) = mode.B_mn(i) + kF_scale * k_F(i);
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

  // Robin+correction: the system matrix includes a Robin BC (iγ₀/μ × M_bdr) for the
  // (0,0) mode. This low-rank operator adds per-mode corrections g_k - g_uniform so
  // that each mode sees its correct DtN eigenvalue.
  auto op = std::make_unique<LowRankComplexOperator>(comm, n);
  for (const auto &mode : modes)
  {
    if (!HasFlag(mode.use, FloquetModeUse::Dtn))
    {
      continue;
    }
    auto g = ComputeDtNCorrectionCoeff(mode);
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

  // Extract E in the true DOF space. Use GetTrueDofs (extracts owned DOF values)
  // instead of P^T × E_local (which sums shared DOF copies and over-counts in parallel).
  int tdof_size = modes.empty() ? 0 : static_cast<int>(modes[0].v.Size());
  Vector E_r_tdof(tdof_size), E_i_tdof(tdof_size);
  E.Real().GetTrueDofs(E_r_tdof);
  if (E.HasImag())
  {
    E.Imag().GetTrueDofs(E_i_tdof);
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

    // Hermitian v^H E = (v_r·E_r + v_i·E_i) + i(v_r·E_i - v_i·E_r).
    double sr = linalg::Dot(comm, mode.v.Real(), E_r_tdof) +
                linalg::Dot(comm, mode.v.Imag(), E_i_tdof);
    double si = linalg::Dot(comm, mode.v.Real(), E_i_tdof) -
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

FloquetPortOperator::FloquetPortOperator(
    const std::map<int, config::FloquetPortData> &floquetport,
    const config::PeriodicBoundaryData &periodic, ProblemType problem_type,
    const Units &units, const MaterialOperator &mat_op_ref,
    mfem::ParFiniteElementSpace &nd_fespace)
  : mat_op(mat_op_ref)
{
  if (floquetport.empty())
  {
    return;
  }
  MFEM_VERIFY(problem_type == ProblemType::DRIVEN,
              "Floquet port boundaries are only available for frequency domain driven "
              "simulations!");

  Mpi::Print("\nConfiguring {:d} Floquet port boundar{}\n", floquetport.size(),
             (floquetport.size() > 1) ? "ies" : "y");

  for (const auto &[idx, data] : floquetport)
  {
    ports.try_emplace(idx, data, periodic, units, mat_op_ref, nd_fespace);
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
  std::unique_ptr<ComplexOperator> combined;
  for (auto &[idx, port] : ports)
  {
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
                                                        MaterialPropertyCoefficient &fbi)
{
  // Add a full-rank Robin BC (iγ₀/μ) on the Floquet port boundary faces.
  Initialize(omega);
  for (const auto &[idx, port] : ports)
  {
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

}  // namespace palace
