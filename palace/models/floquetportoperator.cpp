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
    double dot_r = linalg::Dot(comm, t.v->Real(), x.Real()) -
                   linalg::Dot(comm, t.v->Imag(), x.Imag());
    double dot_i = linalg::Dot(comm, t.v->Real(), x.Imag()) +
                   linalg::Dot(comm, t.v->Imag(), x.Real());

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
// FloquetPortData
// =============================================================================

FloquetPortData::FloquetPortData(const config::FloquetPortData &data,
                                 const IoData &iodata, const MaterialOperator &mat_op_ref,
                                 mfem::ParFiniteElementSpace &nd_fespace)
  : excitation(data.excitation), active(data.active), mat_op(mat_op_ref),
    a1(3), a2(3), b1(3), b2(3), k_F(3), port_normal(3),
    inc_te(data.inc_polarization == "TE"),
    comm(nd_fespace.GetComm())
{
  // Store boundary attributes.
  attr_list.SetSize(static_cast<int>(data.attributes.size()));
  for (int i = 0; i < attr_list.Size(); i++)
  {
    attr_list[i] = data.attributes[i];
  }

  // Extract lattice vectors from periodic BC translations.
  const auto &periodic = iodata.boundaries.periodic;
  MFEM_VERIFY(!periodic.boundary_pairs.empty(),
              "FloquetPort requires periodic boundary conditions!");
  MFEM_VERIFY(periodic.boundary_pairs.size() >= 2,
              "FloquetPort requires at least 2 periodic boundary pairs "
              "(periodicity in 2 directions)!");

  // The translation component of each affine transform is the lattice vector (in mesh-file
  // units). We must nondimensionalize by dividing by the mesh scaling factor Lc/L0, since
  // the mesh coordinates have been nondimensionalized by that factor.
  double mesh_scale = iodata.units.GetMeshLengthRelativeScale();
  for (int i = 0; i < 3; i++)
  {
    a1(i) = periodic.boundary_pairs[0].affine_transform[i * 4 + 3] / mesh_scale;
    a2(i) = periodic.boundary_pairs[1].affine_transform[i * 4 + 3] / mesh_scale;
  }

  // Bloch wave vector (already nondimensionalized in IoData::NondimensionalizeInputs).
  for (int i = 0; i < 3; i++)
  {
    k_F(i) = periodic.wave_vector[i];
  }

  // Compute reciprocal lattice.
  ComputeReciprocalLattice();

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

  // Enumerate diffraction orders and assemble Fourier projections.
  EnumerateOrders();
  AssembleFourierProjections(nd_fespace);

  Mpi::Print(" Floquet port: {:d} modes ({:d} orders x 2 polarizations), "
             "port area = {:.4e}, normal = ({:.3f}, {:.3f}, {:.3f})\n",
             static_cast<int>(modes.size()), static_cast<int>(modes.size()) / 2,
             port_area, port_normal(0), port_normal(1), port_normal(2));
}

void FloquetPortData::ComputeReciprocalLattice()
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

  // b1 = 2*pi * (a2 x n) / (a1 . (a2 x n))
  // b2 = 2*pi * (n x a1) / (a2 . (n x a1))
  // Simplified: b1 = 2*pi * (a2 x n) / |a1 x a2|^2
  mfem::Vector a2xn(3), nxa1(3);
  a2xn(0) = a2(1) * n(2) - a2(2) * n(1);
  a2xn(1) = a2(2) * n(0) - a2(0) * n(2);
  a2xn(2) = a2(0) * n(1) - a2(1) * n(0);

  nxa1(0) = n(1) * a1(2) - n(2) * a1(1);
  nxa1(1) = n(2) * a1(0) - n(0) * a1(2);
  nxa1(2) = n(0) * a1(1) - n(1) * a1(0);

  double vol_sq = vol * vol;
  for (int i = 0; i < 3; i++)
  {
    b1(i) = 2.0 * M_PI * a2xn(i) / vol_sq;
    b2(i) = 2.0 * M_PI * nxa1(i) / vol_sq;
  }

  // Verify: a_i . b_j = 2*pi * delta_ij.
  double a1b1 = a1 * b1;
  double a1b2 = a1 * b2;
  double a2b1 = a2 * b1;
  double a2b2 = a2 * b2;
  MFEM_VERIFY(std::abs(a1b1 - 2.0 * M_PI) < 1e-10 && std::abs(a2b2 - 2.0 * M_PI) < 1e-10,
              "Reciprocal lattice computation failed: diagonal check!");
  MFEM_VERIFY(std::abs(a1b2) < 1e-10 && std::abs(a2b1) < 1e-10,
              "Reciprocal lattice computation failed: off-diagonal check!");

  Mpi::Print(" Floquet port reciprocal lattice:\n"
             "   b1 = ({:.4e}, {:.4e}, {:.4e})\n"
             "   b2 = ({:.4e}, {:.4e}, {:.4e})\n",
             b1(0), b1(1), b1(2), b2(0), b2(1), b2(2));
}

void FloquetPortData::EnumerateOrders()
{
  modes.clear();

  for (int m = -max_order_m; m <= max_order_m; m++)
  {
    for (int n = -max_order_n; n <= max_order_n; n++)
    {
      // B_mn = m*b1 + n*b2
      mfem::Vector B_mn(3);
      for (int i = 0; i < 3; i++)
      {
        B_mn(i) = m * b1(i) + n * b2(i);
      }

      // Propagation direction k_hat for this order.
      // k_t = B_mn (transverse wavevector in the plane, relative to periodic part)
      // The full transverse wavevector including Bloch is B_mn (since we work with
      // the periodic field E_p, the kernel uses B_mn not B_mn - k_F; see floquet_port.md
      // Eq. 38).
      // But gamma uses: gamma^2 = omega^2*mu*eps - |B_mn - k_F|^2 (Eq. in Section 3).
      // We'll compute gamma in Initialize(omega).

      // Compute TE and TM polarization vectors.
      // The transverse wavevector direction for this order:
      mfem::Vector kt(3);
      for (int i = 0; i < 3; i++)
      {
        kt(i) = B_mn(i) - k_F(i);
      }
      double kt_norm = kt.Norml2();

      mfem::Vector e_te(3), e_tm(3);
      if (kt_norm > 1e-12)
      {
        // TE: e_TE = k_hat_t x n_hat / |k_hat_t x n_hat|
        // where k_hat_t is the transverse wavevector direction (normalized kt)
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

        // TM: e_TM = e_TE x k_hat (full propagation direction)
        // At this point we don't know the full k_hat (need gamma), so use the projection:
        // e_TM is in the plane spanned by kt_hat and n_hat, perpendicular to e_TE.
        // e_TM = n_hat x e_TE (to keep it tangential to port face for now)
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
        // Normal incidence for this order (kt = 0): TE/TM are degenerate.
        // Choose two orthogonal directions in the port plane.
        // Use |n̂| (absolute components) to ensure CONSISTENT polarization at all ports
        // regardless of the outward normal sign (opposite faces have opposite normals).
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

        // e_te = ref x |n_hat| (consistent direction, tangential to port)
        e_te(0) = ref(1) * abs_n(2) - ref(2) * abs_n(1);
        e_te(1) = ref(2) * abs_n(0) - ref(0) * abs_n(2);
        e_te(2) = ref(0) * abs_n(1) - ref(1) * abs_n(0);
        double e_te_norm = e_te.Norml2();
        e_te *= 1.0 / e_te_norm;

        // e_tm = |n_hat| x e_te (consistent with e_te)
        e_tm(0) = abs_n(1) * e_te(2) - abs_n(2) * e_te(1);
        e_tm(1) = abs_n(2) * e_te(0) - abs_n(0) * e_te(2);
        e_tm(2) = abs_n(0) * e_te(1) - abs_n(1) * e_te(0);
        double e_tm_norm = e_tm.Norml2();
        e_tm *= 1.0 / e_tm_norm;
      }

      // Add TE mode.
      {
        FloquetMode mode;
        mode.m = m;
        mode.n = n;
        mode.is_te = true;
        mode.B_mn = B_mn;
        mode.e_pol = e_te;
        mode.gamma_sq = 0.0;
        modes.push_back(std::move(mode));
      }

      // Add TM mode.
      {
        FloquetMode mode;
        mode.m = m;
        mode.n = n;
        mode.is_te = false;
        mode.B_mn = B_mn;
        mode.e_pol = e_tm;
        mode.gamma_sq = 0.0;
        modes.push_back(std::move(mode));
      }
    }
  }
}

void FloquetPortData::AssembleFourierProjections(
    mfem::ParFiniteElementSpace &nd_fespace)
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
  mfem::Array<int> bdr_marker =
      mesh::AttrToMarker(bdr_attr_max, attr_list, true);

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

    // Assemble real part.
    mfem::LinearForm lf_r(&nd_fespace);
    lf_r.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(coeff_r),
                               bdr_marker);
    lf_r.UseFastAssembly(false);
    lf_r.UseDevice(false);
    lf_r.Assemble();

    // Assemble imaginary part.
    mfem::LinearForm lf_i(&nd_fespace);
    lf_i.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(coeff_i),
                               bdr_marker);
    lf_i.UseFastAssembly(false);
    lf_i.UseDevice(false);
    lf_i.Assemble();

    // Project to true DOFs: v = P^T * lf.
    nd_fespace.GetProlongationMatrix()->MultTranspose(lf_r, mode.v.Real());
    nd_fespace.GetProlongationMatrix()->MultTranspose(lf_i, mode.v.Imag());

    // Raw projection vectors (not normalized). The F matrix, excitation, and S-parameter
    // extraction all use consistent formulas with these raw vectors.
  }
}

void FloquetPortData::Initialize(double omega)
{
  if (omega == omega0 && omega0 != 0.0)
  {
    return;
  }
  omega0 = omega;

  for (auto &mode : modes)
  {
    // gamma_mn^2 = omega^2 * mu_r * eps_r - |B_mn - k_F|^2
    mfem::Vector kt(3);
    for (int i = 0; i < 3; i++)
    {
      kt(i) = mode.B_mn(i) - k_F(i);
    }
    double kt_sq = kt * kt;
    mode.gamma_sq = omega * omega * mu_eps_port - kt_sq;
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

std::unique_ptr<LowRankComplexOperator> FloquetPortData::GetBoundaryOperator() const
{
  int n = modes.empty() ? 0 : static_cast<int>(modes[0].v.Size());
  auto op = std::make_unique<LowRankComplexOperator>(comm, n);

  // The full-rank Robin BC (iγ₀/μ × boundary_mass) is already added to the system matrix
  // via AddExtraSystemBdrCoefficients. Here we add only the CORRECTION for modes whose
  // propagation constant γ_mn differs from γ₀:
  //   g_correction = i(γ_mn - γ₀)/(μ|Γ|)  for propagating
  //   g_correction = (-|γ_mn| - iγ₀)/(μ|Γ|) for evanescent

  // γ₀₀ for the (0,0) mode: must include |k_F|² for oblique incidence.
    double kF_sq = k_F * k_F;
    double gamma0 = std::sqrt(omega0 * omega0 * mu_eps_port - kF_sq);

  for (const auto &mode : modes)
  {
    std::complex<double> g_full, g_uniform;
    g_uniform = 1i * gamma0 / (mu_r_port * port_area);

    if (mode.gamma_sq > 0.0)
    {
      double gamma = std::sqrt(mode.gamma_sq);
      if (mode.is_te)
      {
        // TE DtN eigenvalue: gamma (from G_nm tensor, Eq 25 of Floquet_BC.pdf)
        g_full = 1i * gamma / (mu_r_port * port_area);
      }
      else
      {
        // TM DtN eigenvalue: omega^2 * mu_r * eps_r / gamma = k^2 / gamma
        // (the G_nm tensor has different eigenvalues for TE and TM polarizations)
        g_full = 1i * omega0 * omega0 * mu_eps_port / (gamma * mu_r_port * port_area);
      }
    }
    else if (mode.gamma_sq < 0.0)
    {
      double gamma_abs = std::sqrt(-mode.gamma_sq);
      if (mode.is_te)
      {
        // TE evanescent: -|gamma| / (mu * |Gamma|)
        g_full = -gamma_abs / (mu_r_port * port_area);
      }
      else
      {
        // TM evanescent: g = j * omega^2 * eps_r / (gamma * |Gamma|) with gamma = j|gamma|
        //              = j * omega^2 * mu_r * eps_r / (j|gamma| * mu_r * |Gamma|)
        //              = +omega^2 * eps_r / (|gamma| * |Gamma|)  (positive real, inductive)
        g_full = omega0 * omega0 * mu_eps_port /
                 (gamma_abs * mu_r_port * port_area);
      }
    }
    else
    {
      continue;
    }

    std::complex<double> g_correction = g_full - g_uniform;
    if (std::abs(g_correction) < 1e-14 * std::abs(g_full))
    {
      continue;  // Skip if correction is negligible (e.g., (0,0) mode).
    }
    op->AddTerm(&mode.v, g_correction);
  }

  return op;
}

std::map<std::tuple<int, int, bool>, std::complex<double>>
FloquetPortData::GetAllSParameters(const GridFunction &E) const
{
  std::map<std::tuple<int, int, bool>, std::complex<double>> result;

  // Power-normalized S-parameter extraction.
  //
  // Step 1: Field amplitude c_mn = v_mn^T E / |Γ| (Fourier coefficient of total field).
  // Step 2: Power normalization: S_mn = √(γ_mn / γ_inc) × c_mn.
  //
  // This ensures Σ |S_mn|² = 1 for energy conservation (lossless). The factor √(γ_mn/γ_inc)
  // accounts for different power-per-amplitude for different diffraction orders (higher-angle
  // orders carry less power per unit amplitude due to oblique propagation).
  //
  // For the incident (0,0) mode: γ_inc = γ_00 = ω √(μ_r ε_r), so √(γ/γ_inc) = 1.
  // For the driving port: S = √(γ/γ_inc) × c - δ_{incident mode} (subtract 1).

  // Compute γ_inc = γ for the (0,0) mode at this port.
  double kF_sq = k_F * k_F;
  double gamma_inc = std::sqrt(omega0 * omega0 * mu_eps_port - kF_sq);

  // The excitation is normalized to inject 1 W. The incident Fourier amplitude is
  // c_inc = 1/√P_unit where P_unit = γ_inc |Γ| / (2ω μ_r). The S-parameter extraction
  // must divide by c_inc so that S-parameters are ratios of power-normalized amplitudes.
  double p_unit = gamma_inc * port_area / (2.0 * omega0 * mu_r_port);
  double c_inc = 1.0 / std::sqrt(p_unit);

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
    if (mode.gamma_sq <= 0.0)
    {
      continue;  // Only propagating orders carry S-parameters.
    }

    double gamma_mn = std::sqrt(mode.gamma_sq);

    // Bilinear v^T E = (v_r·E_r - v_i·E_i) + i(v_r·E_i + v_i·E_r).
    double sr = linalg::Dot(comm, mode.v.Real(), E_r_tdof) -
                linalg::Dot(comm, mode.v.Imag(), E_i_tdof);
    double si = linalg::Dot(comm, mode.v.Real(), E_i_tdof) +
                linalg::Dot(comm, mode.v.Imag(), E_r_tdof);

    // Field amplitude: c = v^T E / (c_inc × |Γ|), where c_inc accounts for the
    // unit-power normalization applied to the excitation.
    // Power normalization: S = √(λ_mn / λ_inc) × c, where λ is the DtN eigenvalue.
    // TE: λ = γ (propagation constant). TM: λ = ω²με/γ (from the DtN tensor).
    // The incident mode is always TE, so λ_inc = γ_inc.
    double lambda_mn =
        mode.is_te ? gamma_mn : omega0 * omega0 * mu_eps_port / gamma_mn;
    double power_factor = std::sqrt(lambda_mn / gamma_inc);

    auto key = std::make_tuple(mode.m, mode.n, mode.is_te);
    result[key] = power_factor * std::complex<double>(sr, si) / (c_inc * port_area);
  }

  return result;
}

bool FloquetPortData::AddExcitationVector(double omega, ComplexVector &RHS) const
{
  // Find the incident mode (0, 0) with matching polarization.
  const FloquetMode *inc_mode = nullptr;
  for (const auto &mode : modes)
  {
    if (mode.m == 0 && mode.n == 0 && mode.is_te == inc_te)
    {
      inc_mode = &mode;
      break;
    }
  }
  if (!inc_mode)
  {
    return false;
  }

  // From the total-field DtN formulation, the RHS for the incident field is:
  //   f = 2 i gamma_inc / mu_r * conj(v_inc)
  // where v_inc,j = ∫_Γ N_j · ê exp(-iB·r) dS, and conj(v) appears because
  // ∫ E_inc · N_k dS = ∫ ê exp(+iB·r) · N_k dS = conj(v_k) for real N_k, ê.
  // The factor 2 comes from the incident + scattered field decomposition.
  MFEM_VERIFY(inc_mode->gamma_sq > 0.0,
              "Incident Floquet mode is evanescent at this frequency!");
  double gamma_inc = std::sqrt(inc_mode->gamma_sq);

  // Unit-power normalization: scale the incident field so the injected power is 1 W
  // (nondimensional), consistent with lumped and wave port conventions.
  // Power with unit Fourier amplitude: P = γ/(2ωμ_r) × |Γ|.
  // Scale factor: c_inc = 1/√P = √(2ωμ_r / (γ|Γ|)).
  double p_unit = gamma_inc * port_area / (2.0 * omega * mu_r_port);
  double c_inc = 1.0 / std::sqrt(p_unit);

  // f = c_inc × 2i gamma / mu_r × conj(v)
  double c_i = c_inc * 2.0 * gamma_inc / mu_r_port;

  RHS.Real().Add(c_i, inc_mode->v.Imag());
  RHS.Imag().Add(c_i, inc_mode->v.Real());
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

  Mpi::Print("\nConfiguring {:d} Floquet port boundar{}\n",
             floquet_port_data.size(),
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
  std::unique_ptr<ComplexOperator> combined;
  for (auto &[idx, port] : ports)
  {
    if (!port.active)
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
      combined =
          std::make_unique<SumComplexOperator>(std::move(combined), std::move(F));
    }
  }

  return combined;
}

void FloquetPortOperator::AddExtraSystemBdrCoefficients(
    double omega, MaterialPropertyCoefficient &fbr, MaterialPropertyCoefficient &fbi)
{
  // Add a full-rank Robin BC (iγ₀/μ) on the Floquet port boundary faces, matching the
  // wave port pattern. This provides proper absorption for ALL tangential field components,
  // not just the included Floquet modes. The low-rank F operator provides corrections for
  // modes with γ_mn ≠ γ₀.
  Initialize(omega);
  for (const auto &[idx, port] : ports)
  {
    if (!port.active)
    {
      continue;
    }
    // γ₀ for the (0,0) mode at this port.
    // γ₀₀ for the (0,0) mode including |k_F|² for oblique incidence.
    double kF_sq = 0.0;
    for (int i = 0; i < 3; i++) { kF_sq += port.k_F(i) * port.k_F(i); }
    double gamma0 = std::sqrt(omega * omega * port.mu_eps_port - kF_sq);
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
