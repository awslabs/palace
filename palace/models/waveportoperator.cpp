// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "waveportoperator.hpp"
#include <limits>
#include "fem/bilinearform.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/arpack.hpp"
#include "linalg/blockprecond.hpp"
#include "linalg/gmg.hpp"
#include "linalg/hypre.hpp"
#include "linalg/iterative.hpp"
#include "linalg/mumps.hpp"
#include "linalg/rap.hpp"
#include "linalg/solver.hpp"

#include <fmt/ranges.h>
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "fem/interpolator.hpp"
#include "linalg/amg.hpp"
#include "linalg/ams.hpp"
#include "linalg/arpack.hpp"
#include "linalg/blockprecond.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/mumps.hpp"
#include "linalg/slepc.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "models/farfieldboundaryoperator.hpp"
#include "models/materialoperator.hpp"
#include "models/surfaceconductivityoperator.hpp"
#include "models/surfaceimpedanceoperator.hpp"
#include "models/surfacerationalimpedanceoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

namespace
{

// Rank-k, complex-symmetric wave-port operator K = Σ_k g_k s_k s_kᵀ. Unlike the Hermitian
// Floquet DtN (v vᴴ), the outer product is unconjugated (s sᵀ), which preserves the Maxwell
// system operator's complex-symmetry and hence S-parameter reciprocity. Used to add the
// modal correction W_full − W_scalar to the sparse local boundary mass: per port a
// +g·s_full s_fullᵀ term (full modal n×H, incl. ∇ₜEₙ) and a −g·s_scalar s_scalarᵀ term
// (scalar-admittance n×H, matching the local mass). On the port's own mode the sum
// reproduces the full modal n×H; for a TEM mode (E_n≈0) s_full=s_scalar so the correction
// vanishes and the baseline is recovered exactly. s_k lives on the ND true-dof space with
// essential dofs eliminated; g_k = −iω / R_k with R_k the unconjugated modal reaction ∫_Γ
// (n×H_k)·E_k.
class WavePortModalCorrection : public ComplexOperator
{
private:
  struct Term
  {
    std::unique_ptr<ComplexVector> s;
    std::complex<double> g;
  };
  std::vector<Term> terms;
  MPI_Comm comm;

public:
  WavePortModalCorrection(MPI_Comm comm, int n) : ComplexOperator(n), comm(comm) {}

  bool Empty() const { return terms.empty(); }

  void AddTerm(std::unique_ptr<ComplexVector> s, std::complex<double> g)
  {
    terms.push_back({std::move(s), g});
  }

  void Mult(const ComplexVector &x, ComplexVector &y) const override
  {
    y = 0.0;
    AddMult(x, y, 1.0);
  }

  void AddMult(const ComplexVector &x, ComplexVector &y,
               const std::complex<double> a) const override
  {
    for (const auto &t : terms)
    {
      // Unconjugated (complex-symmetric) contraction dot = sᵀx = Σ s_j x_j
      // = (s_r·x_r − s_i·x_i) + i(s_r·x_i + s_i·x_r).
      double dot_r = linalg::Dot(comm, t.s->Real(), x.Real()) -
                     linalg::Dot(comm, t.s->Imag(), x.Imag());
      double dot_i = linalg::Dot(comm, t.s->Real(), x.Imag()) +
                     linalg::Dot(comm, t.s->Imag(), x.Real());
      std::complex<double> c = a * t.g * std::complex<double>(dot_r, dot_i);
      // y += c·s = (c_r + i c_i)(s_r + i s_i).
      y.Real().Add(c.real(), t.s->Real());
      y.Real().Add(-c.imag(), t.s->Imag());
      y.Imag().Add(c.imag(), t.s->Real());
      y.Imag().Add(c.real(), t.s->Imag());
    }
  }
};

}  // namespace

namespace
{

void GetEssentialTrueDofs(mfem::ParGridFunction &E0t, mfem::ParGridFunction &E0n,
                          mfem::ParGridFunction &port_E0t, mfem::ParGridFunction &port_E0n,
                          mfem::ParTransferMap &port_nd_transfer,
                          mfem::ParTransferMap &port_h1_transfer,
                          const mfem::Array<int> &dbc_attr,
                          mfem::Array<int> &port_nd_dbc_tdof_list,
                          mfem::Array<int> &port_h1_dbc_tdof_list)
{
  auto &nd_fespace = *E0t.ParFESpace();
  auto &h1_fespace = *E0n.ParFESpace();
  auto &port_nd_fespace = *port_E0t.ParFESpace();
  auto &port_h1_fespace = *port_E0n.ParFESpace();
  const auto &mesh = *nd_fespace.GetParMesh();

  mfem::Array<int> dbc_marker, nd_dbc_tdof_list, h1_dbc_tdof_list;
  mesh::AttrToMarker(mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0, dbc_attr,
                     dbc_marker);
  nd_fespace.GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);

  Vector tE0t(nd_fespace.GetTrueVSize()), tE0n(h1_fespace.GetTrueVSize());
  tE0t.UseDevice(true);
  tE0n.UseDevice(true);
  tE0t = 0.0;
  tE0n = 0.0;
  linalg::SetSubVector(tE0t, nd_dbc_tdof_list, 1.0);
  linalg::SetSubVector(tE0n, h1_dbc_tdof_list, 1.0);
  E0t.SetFromTrueDofs(tE0t);
  E0n.SetFromTrueDofs(tE0n);
  port_nd_transfer.Transfer(E0t, port_E0t);
  port_h1_transfer.Transfer(E0n, port_E0n);

  Vector port_tE0t(port_nd_fespace.GetTrueVSize()),
      port_tE0n(port_h1_fespace.GetTrueVSize());
  port_tE0t.UseDevice(true);
  port_tE0n.UseDevice(true);
  port_E0t.ParallelProject(port_tE0t);
  port_E0n.ParallelProject(port_tE0n);
  {
    const auto *h_port_tE0t = port_tE0t.HostRead();
    const auto *h_port_tE0n = port_tE0n.HostRead();
    for (int i = 0; i < port_tE0t.Size(); i++)
    {
      if (h_port_tE0t[i] != 0.0)
      {
        port_nd_dbc_tdof_list.Append(i);
      }
    }
    for (int i = 0; i < port_tE0n.Size(); i++)
    {
      if (h_port_tE0n[i] != 0.0)
      {
        port_h1_dbc_tdof_list.Append(i);
      }
    }
  }
}

void GetInitialSpace(const mfem::ParFiniteElementSpace &nd_fespace,
                     const mfem::ParFiniteElementSpace &h1_fespace,
                     const mfem::Array<int> &dbc_tdof_list, ComplexVector &v)
{
  // Initial space which satisfies Dirichlet BCs.
  const int nd_size = nd_fespace.GetTrueVSize(), h1_size = h1_fespace.GetTrueVSize();
  v.SetSize(nd_size + h1_size);
  v.UseDevice(true);
  v = std::complex<double>(1.0, 0.0);
  // linalg::SetRandomReal(nd_fespace.GetComm(), v);
  linalg::SetSubVector(v, nd_size, nd_size + h1_size, 0.0);
  linalg::SetSubVector(v, dbc_tdof_list, 0.0);
}

void Normalize(const GridFunction &S0t, GridFunction &E0t, GridFunction &E0n,
               mfem::LinearForm &sr, mfem::LinearForm &si)
{
  // Normalize grid functions to a chosen polarization direction and unit power, |E x H⋆| ⋅
  // n, integrated over the port surface (+n is the outward mesh normal). The n x H
  // coefficients are updated implicitly as they only store references to the Et, En grid
  // functions. We choose a (rather arbitrary) phase constraint to at least make results for
  // the same port consistent between frequencies/meshes.

  // |E x H⋆| ⋅ n = |E ⋅ (-n x H⋆)|. This also updates the n x H coefficients depending on
  // Et, En. Update linear forms for postprocessing too.
  std::complex<double> dot[2] = {
      {sr * S0t.Real(), si * S0t.Real()},
      {-(sr * E0t.Real()) - (si * E0t.Imag()), -(sr * E0t.Imag()) + (si * E0t.Real())}};
  Mpi::GlobalSum(2, dot, S0t.ParFESpace()->GetComm());
  auto scale = std::abs(dot[0]) / (dot[0] * std::sqrt(std::abs(dot[1])));
  ComplexVector::AXPBY(scale, E0t.Real(), E0t.Imag(), 0.0, E0t.Real(), E0t.Imag());
  ComplexVector::AXPBY(scale, E0n.Real(), E0n.Imag(), 0.0, E0n.Real(), E0n.Imag());
  ComplexVector::AXPBY(scale, sr, si, 0.0, sr, si);

  // This parallel communication is not required since wave port boundaries are true one-
  // sided boundaries.
  // E0t.Real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces for n x H
  // E0t.Imag().ExchangeFaceNbrData();  // coefficients evaluation
  // E0n.Real().ExchangeFaceNbrData();
  // E0n.Imag().ExchangeFaceNbrData();
}

// Helper for BdrSubmeshEVectorCoefficient and BdrSubmeshHVectorCoefficient.
enum class ValueType
{
  REAL,
  IMAG
};

// Return as a vector coefficient the boundary mode electric field.
template <ValueType Type>
class BdrSubmeshEVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const GridFunction &Et, &En;
  const mfem::ParSubMesh &submesh;
  const std::unordered_map<int, int> &submesh_parent_elems;
  mfem::IsoparametricTransformation T_loc;
  const double scaling;

public:
  BdrSubmeshEVectorCoefficient(const GridFunction &Et, const GridFunction &En,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems,
                               double scaling = 1.0)
    : mfem::VectorCoefficient(Et.Real().VectorDim()), Et(Et), En(En), submesh(submesh),
      submesh_parent_elems(submesh_parent_elems), scaling(scaling)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Always do the GridFunction evaluation in the submesh.
    mfem::ElementTransformation *T_submesh = nullptr;
    if (T.mesh == submesh.GetParent())
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                  "BdrSubmeshEVectorCoefficient requires ElementType::BDR_ELEMENT when not "
                  "used on a SubMesh!");
      auto it = submesh_parent_elems.find(T.ElementNo);
      if (it == submesh_parent_elems.end())
      {
        // Just return zero for a parent boundary element not in the submesh.
        V.SetSize(vdim);
        V = 0.0;
        return;
      }
      else
      {
        submesh.GetElementTransformation(it->second, &T_loc);
        T_loc.SetIntPoint(&ip);
        T_submesh = &T_loc;
      }
    }
    else if (T.mesh == &submesh)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                  "BdrSubmeshEVectorCoefficient requires ElementType::ELEMENT when used on "
                  "a SubMesh!");
      T_submesh = &T;
    }
    else
    {
      MFEM_ABORT("Invalid mesh for BdrSubmeshEVectorCoefficient!");
    }

    // Compute Eₜ + n ⋅ Eₙ . The normal returned by GetNormal points out of the
    // computational domain, so we reverse it (direction of propagation is into the domain).
    double normal_data[3];
    mfem::Vector normal(normal_data, vdim);
    BdrGridFunctionCoefficient::GetNormal(*T_submesh, normal);
    if constexpr (Type == ValueType::REAL)
    {
      Et.Real().GetVectorValue(*T_submesh, ip, V);
      auto Vn = En.Real().GetValue(*T_submesh, ip);
      V.Add(-Vn, normal);
    }
    else
    {
      Et.Imag().GetVectorValue(*T_submesh, ip, V);
      auto Vn = En.Imag().GetValue(*T_submesh, ip);
      V.Add(-Vn, normal);
    }
    V *= scaling;
  }
};

// Computes boundary mode n x H, where +n is the outward mesh normal: n x H =
// -1/(iωμ) (ik_n Eₜ + ∇ₜ Eₙ), using the tangential and normal electric field component grid
// functions evaluated on the (single-sided) boundary element. With IncludeGradient=false
// the ∇ₜEₙ longitudinal-gradient term is dropped, leaving the scalar-admittance part
// -1/(iωμ) ik_n Eₜ only, which matches the sparse local boundary mass i·k_n·M_{μ⁻¹} stamped
// on the LHS, and is used to form the modal correction W_full − W_scalar.
template <ValueType Type, bool IncludeGradient = true>
class BdrSubmeshHVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const GridFunction &Et, &En;
  const MaterialOperator &mat_op;
  const mfem::ParSubMesh &submesh;
  const std::unordered_map<int, int> &submesh_parent_elems;
  mfem::IsoparametricTransformation T_loc;
  std::complex<double> kn;
  std::complex<double> omega;

public:
  BdrSubmeshHVectorCoefficient(const GridFunction &Et, const GridFunction &En,
                               const MaterialOperator &mat_op,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems,
                               std::complex<double> kn, std::complex<double> omega)
    : mfem::VectorCoefficient(Et.Real().VectorDim()), Et(Et), En(En), mat_op(mat_op),
      submesh(submesh), submesh_parent_elems(submesh_parent_elems), kn(kn), omega(omega)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    // Always do the GridFunction evaluation in the submesh.
    mfem::ElementTransformation *T_submesh = nullptr;
    if (T.mesh == submesh.GetParent())
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                  "BdrSubmeshHVectorCoefficient requires ElementType::BDR_ELEMENT when not "
                  "used on a SubMesh!");
      auto it = submesh_parent_elems.find(T.ElementNo);
      if (it == submesh_parent_elems.end())
      {
        // Just return zero for a parent boundary element not in the submesh.
        V.SetSize(vdim);
        V = 0.0;
        return;
      }
      else
      {
        submesh.GetElementTransformation(it->second, &T_loc);
        T_loc.SetIntPoint(&ip);
        T_submesh = &T_loc;
      }
    }
    else if (T.mesh == &submesh)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                  "BdrSubmeshHVectorCoefficient requires ElementType::ELEMENT when used on "
                  "a SubMesh!");
      T_submesh = &T;
    }
    else
    {
      MFEM_ABORT("Invalid mesh for BdrSubmeshHVectorCoefficient!");
    }

    // Get the attribute in the neighboring domain element of the parent mesh.
    int attr = [&T, this]()
    {
      int i = -1, o, iel1, iel2;
      if (T.mesh == submesh.GetParent())
      {
        MFEM_ASSERT(
            T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
            "BdrSubmeshHVectorCoefficient requires ElementType::BDR_ELEMENT when not "
            "used on a SubMesh!");
        T.mesh->GetBdrElementFace(T.ElementNo, &i, &o);
      }
      else if (T.mesh == &submesh)
      {
        MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                    "BdrSubmeshHVectorCoefficient requires ElementType::ELEMENT when used "
                    "on a SubMesh!");
        submesh.GetParent()->GetBdrElementFace(submesh.GetParentElementIDMap()[T.ElementNo],
                                               &i, &o);
      }
      else
      {
        MFEM_ABORT("Invalid mesh for BdrSubmeshHVectorCoefficient!");
      }
      submesh.GetParent()->GetFaceElements(i, &iel1, &iel2);
      return submesh.GetParent()->GetAttribute(iel1);
    }();

    // Compute U = -1/i (ik_n Eₜ + ∇ₜ Eₙ) = -kₙEₜ + i∇ₜEₙ (t-gradient in boundary element).
    // The ∇ₜEₙ term is omitted when IncludeGradient=false (scalar-admittance-only n×H). kₙ
    // is the full (complex) propagation constant on the eigenmode path (holomorphic in ω,
    // matches the i·kₙ·M mass); the driven path passes Re{kₙ} to stay consistent with its
    // i·Re{kₙ}·M mass (Im{kₙ} line attenuation carried by the sparse mass, not the modal
    // shape). ω may be complex, so carry both parts of U: with Eₜ = Etr + i·Eti, -kₙEₜ =
    // -(kr·Etr - ki·Eti) - i·(kr·Eti + ki·Etr), and i∇ₜEₙ = -∇ₜEni + i·∇ₜEnr.
    const double kr = kn.real(), ki = kn.imag();
    double Etr_data[3], Eti_data[3];
    mfem::Vector Etr(Etr_data, vdim), Eti(Eti_data, vdim);
    Et.Real().GetVectorValue(*T_submesh, ip, Etr);
    Et.Imag().GetVectorValue(*T_submesh, ip, Eti);
    double Ure_data[3], Uim_data[3];
    mfem::Vector Ure(Ure_data, vdim), Uim(Uim_data, vdim);
    for (int d = 0; d < vdim; d++)
    {
      Ure[d] = -(kr * Etr[d] - ki * Eti[d]);
      Uim[d] = -(kr * Eti[d] + ki * Etr[d]);
    }
    if constexpr (IncludeGradient)
    {
      double dU_data[3];
      mfem::Vector dU(dU_data, vdim);
      En.Imag().GetGradient(*T_submesh, dU);
      Ure -= dU;
      En.Real().GetGradient(*T_submesh, dU);
      Uim += dU;
    }

    // Scale by 1/(ωμ) with μ (real) evaluated in the neighboring element. With complex 1/ω
    // = a + i·b the real/imag parts of n×H = (1/ω)·μ⁻¹·(Ure + i·Uim) couple:
    //   Re = a·μ⁻¹Ure - b·μ⁻¹Uim,  Im = a·μ⁻¹Uim + b·μ⁻¹Ure.
    // Note U uses the scalar kₙ from the cross-section EVP; for a rotated μ tensor that
    // couples axes this reconstruction is approximate.
    double MUre_data[3], MUim_data[3];
    mfem::Vector MUre(MUre_data, vdim), MUim(MUim_data, vdim);
    mat_op.GetInvPermeability(attr).Mult(Ure, MUre);
    mat_op.GetInvPermeability(attr).Mult(Uim, MUim);
    const std::complex<double> inv_omega = 1.0 / omega;
    const double a = inv_omega.real(), b = inv_omega.imag();
    V.SetSize(vdim);
    if constexpr (Type == ValueType::REAL)
    {
      for (int d = 0; d < vdim; d++)
      {
        V[d] = a * MUre[d] - b * MUim[d];
      }
    }
    else
    {
      for (int d = 0; d < vdim; d++)
      {
        V[d] = a * MUim[d] + b * MUre[d];
      }
    }
  }
};

// Distribute a back-transformed mode true-dof vector e = [Eₜ (nd); Eₙ (h1)] into the
// submesh grid functions Et, En (parallel scatter). Mirrors the field distribution in
// WavePortData::Initialize.
void DistributeModeField(ComplexVector &e, int nd_size, int h1_size, GridFunction &Et,
                         GridFunction &En)
{
  e.Real().Read();  // Ensure memory is allocated on device before aliasing
  e.Imag().Read();
  Vector etr(e.Real(), 0, nd_size), eti(e.Imag(), 0, nd_size);
  Vector enr(e.Real(), nd_size, h1_size), eni(e.Imag(), nd_size, h1_size);
  etr.UseDevice(true);
  eti.UseDevice(true);
  enr.UseDevice(true);
  eni.UseDevice(true);
  Et.Real().SetFromTrueDofs(etr);
  Et.Imag().SetFromTrueDofs(eti);
  En.Real().SetFromTrueDofs(enr);
  En.Imag().SetFromTrueDofs(eni);
}

// Unconjugated modal reaction R = sᵀe = ∫_Γ (n×H)·Eₜ for the mode field in (Et, En) at
// (kn, omega), where s is the modal n×H (full n×H with IncludeGradient=true, scalar-
// admittance n×H with false). Global sum over comm. Mirrors the reaction assembly in
// WavePortData::Initialize but for arbitrary complex kn/omega and no normalization.
template <bool IncludeGradient>
std::complex<double>
AssembleReaction(const GridFunction &Et, const GridFunction &En,
                 const MaterialOperator &mat_op, const mfem::ParSubMesh &submesh,
                 const std::unordered_map<int, int> &submesh_parent_elems,
                 std::complex<double> kn, std::complex<double> omega,
                 mfem::ParFiniteElementSpace &fes, MPI_Comm comm)
{
  BdrSubmeshHVectorCoefficient<ValueType::REAL, IncludeGradient> nxHr(
      Et, En, mat_op, submesh, submesh_parent_elems, kn, omega);
  BdrSubmeshHVectorCoefficient<ValueType::IMAG, IncludeGradient> nxHi(
      Et, En, mat_op, submesh, submesh_parent_elems, kn, omega);
  mfem::LinearForm sr(&fes), si(&fes);
  sr.AddDomainIntegrator(new VectorFEDomainLFIntegrator(nxHr));
  si.AddDomainIntegrator(new VectorFEDomainLFIntegrator(nxHi));
  for (auto *lf : {&sr, &si})
  {
    lf->UseFastAssembly(false);
    lf->UseDevice(false);
    lf->Assemble();
    lf->UseDevice(true);
  }
  std::complex<double> R = {sr * Et.Real() - si * Et.Imag(),
                            sr * Et.Imag() + si * Et.Real()};
  Mpi::GlobalSum(1, &R, comm);
  return R;
}

}  // namespace

WavePortData::WavePortData(const config::WavePortData &data,
                           const config::BoundaryData &boundaries,
                           const config::DomainData &domains, ProblemType problem_type,
                           const config::LinearSolverData &linear, const Units &units,
                           const MaterialOperator &mat_op,
                           mfem::ParFiniteElementSpace &nd_fespace,
                           mfem::ParFiniteElementSpace &h1_fespace,
                           const mfem::Array<int> &dbc_attr)
  : mat_op(mat_op), excitation(data.excitation), active(data.active),
    include_in_synthesis(data.include_in_synthesis)
{
  mode_idx = data.mode_idx;
  d_offset = data.d_offset;
  kn0 = 0.0;
  omega0 = 0.0;

  // Construct the SubMesh.
  MFEM_VERIFY(!data.attributes.empty(), "Wave port boundary found with no attributes!");
  const auto &mesh = *nd_fespace.GetParMesh();
  attr_list.Append(data.attributes.data(), data.attributes.size());
  auto port_submesh_ptr = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromBoundary(mesh, attr_list));

  // Add internal boundary elements for edges where the port face intersects other boundary
  // faces (PEC, impedance, conductivity, rational impedance, absorbing).
  // ParSubMesh::CreateFromBoundary only creates boundary elements at the geometric boundary
  // of the selected face region, but internal intersections need boundary elements for the
  // 2D eigenvalue problem BCs.
  {
    std::vector<int> internal_bdr_attrs;
    for (auto a : boundaries.pec.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }
    for (auto a : boundaries.auxpec.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }
    for (const auto &d : boundaries.impedance)
    {
      for (auto a : d.attributes)
      {
        internal_bdr_attrs.push_back(a);
      }
    }
    for (const auto &d : boundaries.conductivity)
    {
      for (auto a : d.attributes)
      {
        internal_bdr_attrs.push_back(a);
      }
    }
    for (const auto &d : boundaries.rational_impedance)
    {
      for (auto a : d.attributes)
      {
        internal_bdr_attrs.push_back(a);
      }
    }
    for (auto a : boundaries.farfield.attributes)
    {
      internal_bdr_attrs.push_back(a);
    }
    mesh::AddSubMeshInternalBoundaryElements(*port_submesh_ptr, attr_list,
                                             internal_bdr_attrs);
  }

  port_mesh = std::make_unique<Mesh>(std::move(port_submesh_ptr));
  port_normal = mesh::GetSurfaceNormal(*port_mesh);

  port_nd_fec = std::make_unique<mfem::ND_FECollection>(nd_fespace.GetMaxElementOrder(),
                                                        port_mesh->Dimension());
  port_h1_fec = std::make_unique<mfem::H1_FECollection>(h1_fespace.GetMaxElementOrder(),
                                                        port_mesh->Dimension());
  port_nd_fespace = std::make_unique<FiniteElementSpace>(*port_mesh, port_nd_fec.get());
  port_h1_fespace = std::make_unique<FiniteElementSpace>(*port_mesh, port_h1_fec.get());

  GridFunction E0t(nd_fespace), E0n(h1_fespace);
  port_E0t = std::make_unique<GridFunction>(*port_nd_fespace, true);
  port_E0n = std::make_unique<GridFunction>(*port_h1_fespace, true);
  port_E = std::make_unique<GridFunction>(*port_nd_fespace, true);
  // ω-fields (ComputeComplexReactions): alloc on all ranks; W assembly runs on the full
  // comm.
  port_Et_omega = std::make_unique<GridFunction>(*port_nd_fespace, true);
  port_En_omega = std::make_unique<GridFunction>(*port_h1_fespace, true);

  port_nd_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0t.Real(), port_E0t->Real()));
  port_h1_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0n.Real(), port_E0n->Real()));

  // Remap submesh attributes so that domain elements get the adjacent volume element
  // attribute (matching material definitions) and boundary edges get the adjacent boundary
  // face attribute (matching BC definitions). This must happen AFTER CreateTransferMap
  // (which depends on the original parent-child attribute relationship). Then rebuild the
  // CEED attribute maps from the remapped submesh.
  {
    auto &port_submesh = static_cast<mfem::ParSubMesh &>(port_mesh->Get());
    mesh::RemapSubMeshAttributes(port_submesh);
    mesh::RemapSubMeshBdrAttributes(port_submesh, attr_list);
  }
  port_mesh->RebuildCeedAttributes();

  // Construct submesh-specific MaterialOperator (uses remapped submesh for CEED attribute
  // mapping) and boundary condition operators. The surface operators are constructed with
  // the parent 3D mesh for bdr_attributes validation (modifying a ParSubMesh's
  // bdr_attributes corrupts MFEM internal state), but use port_mat_op for CEED boundary
  // attribute lookups which correctly reference the remapped submesh.
  port_mat_op = std::make_unique<MaterialOperator>(domains.materials, boundaries.periodic,
                                                   problem_type, *port_mesh);
  port_surf_z_op = std::make_unique<SurfaceImpedanceOperator>(
      boundaries.impedance, boundaries.cracked_attributes, units, *port_mat_op, mesh);
  port_farfield_op = std::make_unique<FarfieldBoundaryOperator>(
      boundaries.farfield, problem_type, *port_mat_op, mesh);
  port_surf_sigma_op = std::make_unique<SurfaceConductivityOperator>(
      boundaries.conductivity, problem_type, units, *port_mat_op, mesh);
  port_surf_rz_op = std::make_unique<SurfaceRationalImpedanceOperator>(
      boundaries.rational_impedance, boundaries.cracked_attributes, problem_type, units,
      *port_mat_op, mesh);

  // Construct mapping from parent (boundary) element indices to submesh (domain)
  // elements.
  {
    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    const mfem::Array<int> &parent_elems = port_submesh.GetParentElementIDMap();
    for (int i = 0; i < parent_elems.Size(); i++)
    {
      submesh_parent_elems[parent_elems[i]] = i;
    }
  }

  // Extract Dirichlet BC true dofs for the port FE spaces.
  {
    mfem::Array<int> port_nd_dbc_tdof_list, port_h1_dbc_tdof_list;
    GetEssentialTrueDofs(E0t.Real(), E0n.Real(), port_E0t->Real(), port_E0n->Real(),
                         *port_nd_transfer, *port_h1_transfer, dbc_attr,
                         port_nd_dbc_tdof_list, port_h1_dbc_tdof_list);
    int nd_tdof_offset = port_nd_fespace->GetTrueVSize();
    port_dbc_tdof_list.Reserve(port_nd_dbc_tdof_list.Size() + port_h1_dbc_tdof_list.Size());
    for (auto tdof : port_nd_dbc_tdof_list)
    {
      port_dbc_tdof_list.Append(tdof);
    }
    for (auto tdof : port_h1_dbc_tdof_list)
    {
      port_dbc_tdof_list.Append(tdof + nd_tdof_offset);
    }
  }

  // Create vector for initial space for eigenvalue solves and eigenmode solution.
  GetInitialSpace(*port_nd_fespace, *port_h1_fespace, port_dbc_tdof_list, v0);
  e0.SetSize(port_nd_fespace->GetTrueVSize() + port_h1_fespace->GetTrueVSize());
  e0.UseDevice(true);

  // Set the shift-and-invert target above the bulk propagation constants of materials
  // present on this port. For anisotropic media the electric and magnetic polarizations
  // can sample different material axes, so max(mu_r * epsilon_r) is not a safe bound;
  // use max(mu_r) * max(epsilon_r) instead.
  mu_eps_max = port_mat_op->GetMaxMuEpsilon() * 1.1;

  // Configure a communicator for the processes which have elements for this port.
  MPI_Comm comm = nd_fespace.GetComm();
  int color = (port_nd_fespace->GetVSize() > 0 || port_h1_fespace->GetVSize() > 0)
                  ? 0
                  : MPI_UNDEFINED;
  MPI_Comm_split(comm, color, Mpi::Rank(comm), &port_comm);
  MFEM_VERIFY((color == 0 && port_comm != MPI_COMM_NULL) ||
                  (color == MPI_UNDEFINED && port_comm == MPI_COMM_NULL),
              "Unexpected error splitting communicator for wave port boundaries!");
  port_root = (color == MPI_UNDEFINED) ? Mpi::Size(comm) : Mpi::Rank(comm);
  Mpi::GlobalMin(1, &port_root, comm);
  MFEM_VERIFY(port_root < Mpi::Size(comm), "No root process found for port!");

  // Configure the boundary mode solver. Matrix assembly is MPI-collective on the FE space
  // communicator (all processes), so the config + construction must happen on all
  // processes. The solver_comm (port_comm) restricts solver setup to port processes only.
  {
    mode_solver = std::make_unique<ModeEigenSolver>(
        *port_mat_op, &port_normal, *port_surf_z_op, *port_farfield_op, *port_surf_sigma_op,
        *port_surf_rz_op, *port_nd_fespace, *port_h1_fespace, port_dbc_tdof_list, mode_idx,
        data.max_size, data.eig_tol, EigenvalueSolver::WhichType::LARGEST_REAL, linear,
        data.eigen_solver, data.verbose, port_comm);
  }

  // Configure port mode sign convention: 1ᵀ Re{-n x H} >= 0 on the "upper-right quadrant"
  // of the wave port boundary, in order to deal with symmetry effectively. The user can
  // override this convention by providing a VoltagePath; see the post-Normalize step in
  // Initialize() which flips the mode sign such that ∫ E · dl > 0 along the path.
  {
    Vector bbmin, bbmax;
    mesh::GetAxisAlignedBoundingBox(*port_mesh, bbmin, bbmax);
    const int dim = port_mesh->SpaceDimension();

    double la = 0.0, lb = 0.0;
    int da = -1, db = -1;
    for (int d = 0; d < dim; d++)
    {
      double diff = bbmax(d) - bbmin(d);
      if (diff > la)
      {
        lb = la;
        la = diff;
        db = da;
        da = d;
      }
      else if (diff > lb)
      {
        lb = diff;
        db = d;
      }
    }
    MFEM_VERIFY(da >= 0 && db >= 0 && da != db,
                "Unexpected wave port geometry for normalization!");
    double ca = 0.5 * (bbmax[da] + bbmin[da]), cb = 0.5 * (bbmax[db] + bbmin[db]);

    auto TDirection = [da, db, ca, cb, dim](const Vector &x, Vector &f)
    {
      MFEM_ASSERT(x.Size() == dim,
                  "Invalid dimension mismatch for wave port mode normalization!");
      f.SetSize(dim);
      if (x[da] >= ca && x[db] >= cb)
      {
        f = 1.0;
      }
      else
      {
        f = 0.0;
      }
    };
    mfem::VectorFunctionCoefficient tfunc(dim, TDirection);
    port_S0t = std::make_unique<GridFunction>(*port_nd_fespace);
    port_S0t->Real().ProjectCoefficient(tfunc);
  }

  // Store voltage path coordinates if provided.
  // Coordinates are already nondimensionalized in IoData::NondimensionalizeInputs.
  if (data.voltage_path.size() >= 2)
  {
    has_voltage_coords = true;
    for (const auto &pt : data.voltage_path)
    {
      mfem::Vector p(pt.size());
      for (int d = 0; d < static_cast<int>(pt.size()); d++)
      {
        p(d) = pt[d];
      }
      voltage_path.push_back(std::move(p));
    }
    voltage_n_samples = data.n_samples;

    // Set up reverse transfer map (port submesh → parent mesh) and a parent-mesh
    // GridFunction to receive the transferred mode field. This enables computing line
    // integrals of the port mode E-field via GSLIB on the 3D parent mesh, since GSLIB
    // requires SpaceDim == Dim and cannot work directly on the 2D-embedded port submesh.
    parent_E0t = std::make_unique<GridFunction>(nd_fespace, true);
    port_nd_transfer_reverse = std::make_unique<mfem::ParTransferMap>(
        mfem::ParSubMesh::CreateTransferMap(port_E0t->Real(), parent_E0t->Real()));

#if defined(MFEM_USE_GSLIB)
    // Build the GSLIB point locator on the (fixed) parent mesh once here.
    // GetExcitationVoltage reuses it for every line integral instead of rebuilding the
    // spatial hash per call — the dominant cost when the mode is evaluated at many
    // frequencies (e.g. synthesis fit).
    auto &parent_mesh = *parent_E0t->Real().FESpace()->GetMesh();
    voltage_gslib_op =
        std::make_unique<mfem::FindPointsGSLIB>(parent_E0t->Real().ParFESpace()->GetComm());
    fem::SetupInterpolator(*voltage_gslib_op, parent_mesh);
#endif
  }

  // Store polarity attributes [high, low]. The config parser enforces that this is
  // mutually exclusive with VoltagePath, and zero means "not set".
  polarity_attributes = data.polarity_attributes;
}

WavePortData::~WavePortData()
{
  // Free the solver before the communicator on which it is based.
  mode_solver.reset();
  if (port_comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&port_comm);
  }
}

void WavePortData::SetSynthesisEigTol(double eig_tol, double ksp_tol)
{
  mode_solver->SetSynthesisTol(eig_tol, ksp_tol);
  // Invalidate the cached real-ω solve (omega0 == 0 initially, and physical ω > 0) so the
  // next Initialize re-solves the cross-section EVP at the tightened tolerance.
  omega0 = -1.0;
}

void WavePortData::Initialize(double omega)
{
  if (omega == omega0)
  {
    return;
  }

  // Solve the generalized eigenvalue problem for the desired wave port mode using the
  // ModeEigenSolver. Frequency-independent matrices were assembled in the constructor.
  const double sigma = -omega * omega * mu_eps_max;
  std::complex<double> lambda;
  {
    const bool has_solver = (port_comm != MPI_COMM_NULL);
    auto result = mode_solver->Solve(omega, sigma, has_solver ? &v0 : nullptr);
    if (has_solver)
    {
      MFEM_VERIFY(result.num_converged >= mode_idx,
                  "Wave port eigensolver did not converge!");
      lambda = mode_solver->GetEigenvalue(mode_idx - 1);
    }
  }
  Mpi::Broadcast(1, &lambda, port_root, port_mesh->GetComm());

  // Extract the eigenmode solution and postprocess. The extracted eigenvalue is λ =
  // 1 / (-k_n² - σ).
  kn0 = std::sqrt(-sigma - 1.0 / lambda);
  omega0 = omega;

  // Separate the computed field out into eₜ and eₙ and transform back to the physical
  // electric field variables Eₜ = eₜ and Eₙ = eₙ / (i·k_n). Order: load raw eigenvector,
  // phase-normalize, then apply the shared VD back-transform.
  {
    if (port_comm != MPI_COMM_NULL)
    {
      mode_solver->GetEigenvector(mode_idx - 1, e0);
      linalg::NormalizePhase(port_comm, e0);
      ComplexVector et_view, en_view;
      mode_assembly::ApplyVDBackTransform(e0, kn0, port_nd_fespace->GetTrueVSize(),
                                          port_h1_fespace->GetTrueVSize(), et_view,
                                          en_view);
    }
    else
    {
      MFEM_ASSERT(e0.Size() == 0,
                  "Unexpected non-empty port FE space in wave port boundary mode solve!");
    }
    e0.Real().Read();  // Ensure memory is allocated on device before aliasing
    e0.Imag().Read();
    Vector e0tr(e0.Real(), 0, port_nd_fespace->GetTrueVSize());
    Vector e0nr(e0.Real(), port_nd_fespace->GetTrueVSize(),
                port_h1_fespace->GetTrueVSize());
    Vector e0ti(e0.Imag(), 0, port_nd_fespace->GetTrueVSize());
    Vector e0ni(e0.Imag(), port_nd_fespace->GetTrueVSize(),
                port_h1_fespace->GetTrueVSize());
    e0tr.UseDevice(true);
    e0nr.UseDevice(true);
    e0ti.UseDevice(true);
    e0ni.UseDevice(true);
    port_E0t->Real().SetFromTrueDofs(e0tr);  // Parallel distribute
    port_E0t->Imag().SetFromTrueDofs(e0ti);
    port_E0n->Real().SetFromTrueDofs(e0nr);
    port_E0n->Imag().SetFromTrueDofs(e0ni);
  }

  // Configure the linear forms for computing S-parameters (projection of the field onto the
  // port mode). Normalize the mode for a chosen polarization direction and unit power,
  // |E x H⋆| ⋅ n, integrated over the port surface (+n is the outward mesh normal).
  {
    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    BdrSubmeshHVectorCoefficient<ValueType::REAL> port_nxH0r_func(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0.real(),
        omega0);
    BdrSubmeshHVectorCoefficient<ValueType::IMAG> port_nxH0i_func(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0.real(),
        omega0);
    {
      port_sr = std::make_unique<mfem::LinearForm>(&port_nd_fespace->Get());
      port_sr->AddDomainIntegrator(new VectorFEDomainLFIntegrator(port_nxH0r_func));
      port_sr->UseFastAssembly(false);
      port_sr->UseDevice(false);
      port_sr->Assemble();
      port_sr->UseDevice(true);
    }
    {
      port_si = std::make_unique<mfem::LinearForm>(&port_nd_fespace->Get());
      port_si->AddDomainIntegrator(new VectorFEDomainLFIntegrator(port_nxH0i_func));
      port_si->UseFastAssembly(false);
      port_si->UseDevice(false);
      port_si->Assemble();
      port_si->UseDevice(true);
    }
    Normalize(*port_S0t, *port_E0t, *port_E0n, *port_sr, *port_si);

    // If the user provided a VoltagePath, use it to fix the mode polarity such that
    // V_exc = ∫ E_mode · dl is real-positive along the path. This ties the wave port
    // mode polarity to a physically meaningful direction (the path direction, like a
    // lumped port's Direction), enabling consistent S-parameter signs in mixed lumped
    // + wave port simulations. PolarityAttributes provides a lightweight (no-GSLIB)
    // alternative when only polarity is needed.
    int sign = 0;
    if (has_voltage_coords)
    {
      auto V_exc = GetExcitationVoltage();
      sign = (V_exc.real() < 0.0) ? -1 : (V_exc.real() > 0.0 ? +1 : 0);
    }
    else if (polarity_attributes[0] != 0 && polarity_attributes[1] != 0)
    {
      sign = GetModePolaritySign(polarity_attributes[0], polarity_attributes[1]);
    }
    if (sign < 0)
    {
      ComplexVector::AXPBY(std::complex<double>(-1.0, 0.0), port_E0t->Real(),
                           port_E0t->Imag(), 0.0, port_E0t->Real(), port_E0t->Imag());
      ComplexVector::AXPBY(std::complex<double>(-1.0, 0.0), port_E0n->Real(),
                           port_E0n->Imag(), 0.0, port_E0n->Real(), port_E0n->Imag());
      ComplexVector::AXPBY(std::complex<double>(-1.0, 0.0), *port_sr, *port_si, 0.0,
                           *port_sr, *port_si);
    }

    // Unconjugated modal reaction R = sᵀe = ∫_Γ (n×H_mode)·E_mode, which scales the rank-1
    // terms of the modal correction W = (−iω/R) s sᵀ. Computed from the same normalized,
    // sign-fixed submesh forms/field that define the excitation n×H, so R matches the 3D-
    // assembled s. Distinct from Normalize's conjugated power sᴴe. Sign-flip invariant (s
    // and E both flip). s = (sr + i·si), e = (E0t.Re + i·E0t.Im).
    modal_reaction = {(*port_sr) * port_E0t->Real() - (*port_si) * port_E0t->Imag(),
                      (*port_sr) * port_E0t->Imag() + (*port_si) * port_E0t->Real()};
    Mpi::GlobalSum(1, &modal_reaction, port_nd_fespace->GetComm());

    // Scalar-admittance reaction R_scalar = s_scalarᵀe from the scalar-only n×H (∇ₜEₙ
    // dropped), scaling the W_scalar term that cancels the local boundary mass's modal
    // action so the correction W_full − W_scalar adds only the ∇ₜEₙ contribution. Assembled
    // from the final (normalized, sign-fixed) mode field so it is consistent with
    // modal_reaction.
    BdrSubmeshHVectorCoefficient<ValueType::REAL, false> port_nxH0r_scalar(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0.real(),
        omega0);
    BdrSubmeshHVectorCoefficient<ValueType::IMAG, false> port_nxH0i_scalar(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0.real(),
        omega0);
    mfem::LinearForm sr_scalar(&port_nd_fespace->Get()), si_scalar(&port_nd_fespace->Get());
    sr_scalar.AddDomainIntegrator(new VectorFEDomainLFIntegrator(port_nxH0r_scalar));
    si_scalar.AddDomainIntegrator(new VectorFEDomainLFIntegrator(port_nxH0i_scalar));
    for (auto *lf : {&sr_scalar, &si_scalar})
    {
      lf->UseFastAssembly(false);
      lf->UseDevice(false);
      lf->Assemble();
      lf->UseDevice(true);
    }
    modal_reaction_scalar = {sr_scalar * port_E0t->Real() - si_scalar * port_E0t->Imag(),
                             sr_scalar * port_E0t->Imag() + si_scalar * port_E0t->Real()};
    Mpi::GlobalSum(1, &modal_reaction_scalar, port_nd_fespace->GetComm());
  }
}

const WavePortData::ModeSolveCache &WavePortData::EvalModeSolve(std::complex<double> omega)
{
  // Solve the cross-section EVP once at complex ω and cache (kₙ, eigenvector), reusing the
  // cached result on a repeat call at the same ω. Mirrors the EVP solve + recovery in
  // Initialize() but skips field reconstruction / normalization and does NOT touch the
  // cached real-ω state (omega0, kn0, port_E0t/E0n). The spectral shift sigma stays REAL —
  // it is a pure algebraic centering of the linearization (exact for any real sigma), so we
  // derive it from the real part of the requested frequency; the cross-section EVP carries
  // the full complex ω.
  if (mode_cache.valid && mode_cache.omega == omega)
  {
    return mode_cache;
  }
  const double omega_ref = omega.real();
  const double sigma = -omega_ref * omega_ref * mu_eps_max;
  const bool has_solver = (port_comm != MPI_COMM_NULL);
  std::complex<double> lambda;
  ComplexVector e_tmp;
  {
    auto result = mode_solver->Solve(omega, sigma, has_solver ? &v0 : nullptr);
    if (has_solver)
    {
      MFEM_VERIFY(result.num_converged >= mode_idx,
                  "Wave port eigensolver did not converge in EvalModeSolve!");
      lambda = mode_solver->GetEigenvalue(mode_idx - 1);
      e_tmp.SetSize(port_nd_fespace->GetTrueVSize() + port_h1_fespace->GetTrueVSize());
      e_tmp.UseDevice(true);
      mode_solver->GetEigenvector(mode_idx - 1, e_tmp);
    }
  }
  Mpi::Broadcast(1, &lambda, port_root, port_mesh->GetComm());
  // k_n = √(−σ − 1/λ_evp). Principal branch gives Re(k_n) ≥ 0 (forward-propagating /
  // forward-decaying sheet), which is the physical convention for the +z mode.
  mode_cache.kn = std::sqrt(-sigma - 1.0 / lambda);
  mode_cache.e = std::move(e_tmp);
  mode_cache.omega = omega;
  mode_cache.valid = true;
  return mode_cache;
}

std::complex<double> WavePortData::SolveKnComplex(std::complex<double> omega)
{
  // Complex-frequency propagation constant, from the shared cached EVP solve. Used by the
  // eigenmode nonlinear solver to evaluate the wave-port BC at the genuinely complex
  // eigenvalue. For real ω this returns the same k_n as Initialize(ω).
  return EvalModeSolve(omega).kn;
}

WavePortData::ComplexReactions
WavePortData::ComputeComplexReactions(std::complex<double> omega)
{
  // Recover kₙ and the eigenvector from the shared cached EVP solve at ω (so kₙ and the
  // reactions come from the same eigenbranch as SolveKnComplex), then on the port ranks
  // recompute the complete mode, scale-lock it to the frozen reference e0, and transport
  // the frozen unit-power reactions by the reaction ratio. Broadcast (kn, R_full, R_scalar)
  // from the port root; falls back to the frozen reactions on ranks/paths without a solver.
  const bool has_solver = (port_comm != MPI_COMM_NULL);
  const auto &solve = EvalModeSolve(omega);
  const std::complex<double> kn = solve.kn;
  // Copy the cached eigenvector: the back-transform below mutates it in place, and the
  // cache must stay pristine for a later SolveKnComplex / ComputeComplexReactions at the
  // same ω.
  ComplexVector e_tmp(solve.e);

  std::complex<double> results[5] = {kn, modal_reaction, modal_reaction_scalar, 0.0, 0.0};
  if (has_solver)
  {
    const int nd = port_nd_fespace->GetTrueVSize();
    const int h1 = port_h1_fespace->GetTrueVSize();
    const auto &submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());

    // Raw-gauge reference reactions of the frozen field e0 at (kn0, ω0). modal_reaction =
    // scale²·R_full_ref with the frozen unit-power scale, so R_full(ω)/R_full_ref
    // transported by modal_reaction stays unit-power. Recomputed only when the reference
    // (ω0) changes.
    if (ref_reaction_omega0 != omega0)
    {
      GridFunction Et0(*port_nd_fespace, true), En0(*port_h1_fespace, true);
      DistributeModeField(e0, nd, h1, Et0, En0);
      R_full_ref =
          AssembleReaction<true>(Et0, En0, mat_op, submesh, submesh_parent_elems,
                                 kn0.real(), omega0, port_nd_fespace->Get(), port_comm);
      R_scalar_ref =
          AssembleReaction<false>(Et0, En0, mat_op, submesh, submesh_parent_elems,
                                  kn0.real(), omega0, port_nd_fespace->Get(), port_comm);
      ref_reaction_omega0 = omega0;
    }

    // Back-transform the complex mode (Eₙ = ẽₙ/(i·kₙ)) and scale-lock it to e0 via the
    // least- squares factor α (α·e_tmp ≈ e0): α = ⟨e_tmp, e0⟩ / ⟨e_tmp, e_tmp⟩. Reactions
    // are quadratic in the field, so the locked reaction α²·R_raw is invariant to the EVP's
    // arbitrary per-solve scale/phase and continuously connected to the ω0 gauge.
    ComplexVector et_view, en_view;
    mode_assembly::ApplyVDBackTransform(e_tmp, kn, nd, h1, et_view, en_view);
    const std::complex<double> num = linalg::Dot(port_comm, e0, e_tmp);  // e_tmpᴴ e0
    const double den = linalg::Dot(port_comm, e_tmp, e_tmp).real();
    const std::complex<double> a2 =
        (den > 0.0) ? (num / den) * (num / den) : std::complex<double>(0.0);

    // Stash the recomputed field (at kn, omega) for the eigenmode full-mode correction. The
    // raw reactions are paired with shape vectors from this same field, so the arbitrary
    // EVP scale/phase cancels in W = (−iω/R) s sᵀ (no gauge-lock needed here, unlike the
    // transported reactions below which feed the frozen-shape synthesis pencil).
    DistributeModeField(e_tmp, nd, h1, *port_Et_omega, *port_En_omega);
    kn_recompute = kn;
    omega_recompute = omega;
    const std::complex<double> R_full_raw = AssembleReaction<true>(
        *port_Et_omega, *port_En_omega, mat_op, submesh, submesh_parent_elems, kn, omega,
        port_nd_fespace->Get(), port_comm);
    const std::complex<double> R_scalar_raw = AssembleReaction<false>(
        *port_Et_omega, *port_En_omega, mat_op, submesh, submesh_parent_elems, kn, omega,
        port_nd_fespace->Get(), port_comm);
    results[3] = R_full_raw;
    results[4] = R_scalar_raw;

    if (std::abs(R_full_ref) > 0.0)
    {
      results[1] = modal_reaction * (a2 * R_full_raw) / R_full_ref;
    }
    if (std::abs(R_scalar_ref) > 0.0)
    {
      results[2] = modal_reaction_scalar * (a2 * R_scalar_raw) / R_scalar_ref;
    }
  }
  Mpi::Broadcast(5, results, port_root, port_mesh->GetComm());
  return {results[0], results[1], results[2], results[3], results[4]};
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeExcitationCoefficientReal(bool include_gradient) const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  if (include_gradient)
  {
    return std::make_unique<
        RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::REAL, true>>>(
        attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems,
        kn0.real(), omega0);
  }
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::REAL, false>>>(
      attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems,
      kn0.real(), omega0);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeExcitationCoefficientImag(bool include_gradient) const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  if (include_gradient)
  {
    return std::make_unique<
        RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::IMAG, true>>>(
        attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems,
        kn0.real(), omega0);
  }
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::IMAG, false>>>(
      attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems,
      kn0.real(), omega0);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetOmegaModeExcitationCoefficientReal(bool include_gradient) const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  if (include_gradient)
  {
    return std::make_unique<
        RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::REAL, true>>>(
        attr_list, *port_Et_omega, *port_En_omega, mat_op, port_submesh,
        submesh_parent_elems, kn_recompute, omega_recompute);
  }
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::REAL, false>>>(
      attr_list, *port_Et_omega, *port_En_omega, mat_op, port_submesh, submesh_parent_elems,
      kn_recompute, omega_recompute);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetOmegaModeExcitationCoefficientImag(bool include_gradient) const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  if (include_gradient)
  {
    return std::make_unique<
        RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::IMAG, true>>>(
        attr_list, *port_Et_omega, *port_En_omega, mat_op, port_submesh,
        submesh_parent_elems, kn_recompute, omega_recompute);
  }
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::IMAG, false>>>(
      attr_list, *port_Et_omega, *port_En_omega, mat_op, port_submesh, submesh_parent_elems,
      kn_recompute, omega_recompute);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeFieldCoefficientReal(double scaling) const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshEVectorCoefficient<ValueType::REAL>>>(
      attr_list, *port_E0t, *port_E0n, port_submesh, submesh_parent_elems, scaling);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeFieldCoefficientImag(double scaling) const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshEVectorCoefficient<ValueType::IMAG>>>(
      attr_list, *port_E0t, *port_E0n, port_submesh, submesh_parent_elems, scaling);
}

double WavePortData::GetExcitationPower() const
{
  // The computed port modes are normalized such that the power integrated over the port is
  // 1: ∫ (E_inc x H_inc⋆) ⋅ n dS = 1.
  return HasExcitation() ? 1.0 : 0.0;
}

std::complex<double> WavePortData::GetPower(GridFunction &E, GridFunction &B) const
{
  // Compute port power, (E x H) ⋅ n = E ⋅ (-n x H), integrated over the port surface using
  // the computed E and H = μ⁻¹ B fields, where +n is the outward mesh normal. The
  // BdrSurfaceCurrentVectorCoefficient computes -n x H for the outward normal. The linear
  // form is reconstructed from scratch each time due to changing H.
  MFEM_VERIFY(E.HasImag() && B.HasImag(),
              "Wave ports expect complex-valued E and B fields in port power "
              "calculation!");
  auto &nd_fespace = *E.ParFESpace();
  const auto &mesh = *nd_fespace.GetParMesh();
  BdrSurfaceCurrentVectorCoefficient nxHr_func(B.Real(), mat_op);
  BdrSurfaceCurrentVectorCoefficient nxHi_func(B.Imag(), mat_op);
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  mfem::Array<int> attr_marker = mesh::AttrToMarker(bdr_attr_max, attr_list);
  std::complex<double> dot;
  {
    mfem::LinearForm pr(&nd_fespace);
    pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHr_func), attr_marker);
    pr.UseFastAssembly(false);
    pr.UseDevice(false);
    pr.Assemble();
    pr.UseDevice(true);
    dot = -(pr * E.Real()) - 1i * (pr * E.Imag());
  }
  {
    mfem::LinearForm pi(&nd_fespace);
    pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHi_func), attr_marker);
    pi.UseFastAssembly(false);
    pi.UseDevice(false);
    pi.Assemble();
    pi.UseDevice(true);
    dot += -(pi * E.Imag()) + 1i * (pi * E.Real());
  }
  Mpi::GlobalSum(1, &dot, nd_fespace.GetComm());
  return dot;
}

std::complex<double> WavePortData::GetSParameter(GridFunction &E) const
{
  // Compute port S-parameter, or the projection of the field onto the port mode:
  // (E x H_inc⋆) ⋅ n = E ⋅ (-n x H_inc⋆), integrated over the port surface.
  MFEM_VERIFY(E.HasImag(),
              "Wave ports expect complex-valued E and B fields in port S-parameter "
              "calculation!");
  port_nd_transfer->Transfer(E.Real(), port_E->Real());
  port_nd_transfer->Transfer(E.Imag(), port_E->Imag());
  std::complex<double> dot(-((*port_sr) * port_E->Real()) - ((*port_si) * port_E->Imag()),
                           -((*port_sr) * port_E->Imag()) + ((*port_si) * port_E->Real()));
  Mpi::GlobalSum(1, &dot, port_nd_fespace->GetComm());
  return dot;
}

int WavePortData::GetModePolaritySign(int high_attr, int low_attr) const
{
  // Without GSLIB, evaluate the modal E-field at the centroids of two port-submesh
  // boundary elements — one with parent attribute high_attr (signal terminal) and one
  // with low_attr (ground terminal) — and return the sign of
  //   Re( E_avg · (c_low - c_high) ),
  // i.e. positive when E points from the signal (high) terminal toward the ground (low)
  // terminal, matching the lumped-port `+R Direction` convention and VoltagePath
  // ordering. After RemapSubMeshBdrAttributes the submesh boundary element's attribute
  // is the parent boundary attribute, so we can match directly against the user ints.
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  const int sdim = port_submesh.SpaceDimension();

  // Local search: for each target attribute, record the first matching boundary
  // element (lowest local index) and compute its centroid and the modal E vector at
  // that centroid (evaluated through the adjacent 2D submesh element).
  auto LocalFind = [&](int target_attr, mfem::Vector &face_centroid, mfem::Vector &E_re,
                       mfem::Vector &E_im, int &found_rank)
  {
    found_rank = std::numeric_limits<int>::max();
    face_centroid.SetSize(sdim);
    E_re.SetSize(sdim);
    E_im.SetSize(sdim);
    face_centroid = 0.0;
    E_re = 0.0;
    E_im = 0.0;
    for (int b = 0; b < port_submesh.GetNBE(); b++)
    {
      if (port_submesh.GetBdrAttribute(b) != target_attr)
      {
        continue;
      }
      mfem::FaceElementTransformations *FT =
          const_cast<mfem::ParSubMesh &>(port_submesh).GetBdrFaceTransformations(b);
      if (FT == nullptr)
      {
        continue;
      }
      // Use the boundary-face centroid (a physical point ON the terminal) as the
      // direction reference, but evaluate E in the interior of the adjacent element
      // since the modal tangential E vanishes on PEC boundaries.
      mfem::IntegrationPoint ip_face;
      ip_face.x = ip_face.y = ip_face.z = 0.5;
      FT->SetAllIntPoints(&ip_face);
      FT->Transform(ip_face, face_centroid);
      mfem::ElementTransformation *T_elem = FT->Elem1;
      mfem::Geometry::Type geom = T_elem->GetGeometryType();
      const mfem::IntegrationPoint &ip_center = mfem::Geometries.GetCenter(geom);
      T_elem->SetIntPoint(&ip_center);
      port_E0t->Real().GetVectorValue(*T_elem, ip_center, E_re);
      port_E0t->Imag().GetVectorValue(*T_elem, ip_center, E_im);
      found_rank = Mpi::Rank(port_submesh.GetComm());
      return;
    }
  };

  mfem::Vector c_high(sdim), c_low(sdim), E_re_high(sdim), E_im_high(sdim);
  mfem::Vector E_re_low(sdim), E_im_low(sdim);
  int rank_high = 0, rank_low = 0;
  LocalFind(high_attr, c_high, E_re_high, E_im_high, rank_high);
  LocalFind(low_attr, c_low, E_re_low, E_im_low, rank_low);

  // Pick the lowest-rank winner for each attribute deterministically, then broadcast
  // its centroid + E values to all ranks.
  MPI_Comm comm = port_submesh.GetComm();
  Mpi::GlobalMin(1, &rank_high, comm);
  Mpi::GlobalMin(1, &rank_low, comm);
  if (rank_high >= Mpi::Size(comm) || rank_low >= Mpi::Size(comm))
  {
    Mpi::Warning(comm,
                 "WavePort PolarityAttributes [{}, {}] not found on the port boundary; "
                 "skipping polarity fix!\n",
                 high_attr, low_attr);
    return 0;
  }
  Mpi::Broadcast(sdim, c_high.GetData(), rank_high, comm);
  Mpi::Broadcast(sdim, E_re_high.GetData(), rank_high, comm);
  Mpi::Broadcast(sdim, E_im_high.GetData(), rank_high, comm);
  Mpi::Broadcast(sdim, c_low.GetData(), rank_low, comm);
  Mpi::Broadcast(sdim, E_re_low.GetData(), rank_low, comm);
  Mpi::Broadcast(sdim, E_im_low.GetData(), rank_low, comm);

  mfem::Vector dir(sdim);
  for (int d = 0; d < sdim; d++)
  {
    dir(d) = c_low(d) - c_high(d);
  }
  // Sign of Re( ((E_high + E_low)/2) · dir ).
  double dot_re = 0.0;
  for (int d = 0; d < sdim; d++)
  {
    dot_re += 0.5 * (E_re_high(d) + E_re_low(d)) * dir(d);
  }
  return (dot_re < 0.0) ? -1 : (dot_re > 0.0 ? +1 : 0);
}

std::complex<double> WavePortData::GetVoltage(GridFunction &E) const
{
  // Compute voltage V = ∫ E · dl along the configured path segments.
  // Uses GSLIB interpolation on the 3D parent mesh E field.
  if (!has_voltage_coords)
  {
    return 0.0;
  }
  MFEM_VERIFY(E.HasImag(),
              "Wave ports expect complex-valued E field in port voltage calculation!");
  std::complex<double> V(0.0, 0.0);
#if defined(MFEM_USE_GSLIB)
  // Reuse the cached point locator: E lives on the same parent mesh the op was Setup on.
  for (std::size_t k = 0; k + 1 < voltage_path.size(); k++)
  {
    V.real(V.real() + fem::ComputeLineIntegral(*voltage_gslib_op, voltage_path[k],
                                               voltage_path[k + 1], E.Real(),
                                               voltage_n_samples));
    V.imag(V.imag() + fem::ComputeLineIntegral(*voltage_gslib_op, voltage_path[k],
                                               voltage_path[k + 1], E.Imag(),
                                               voltage_n_samples));
  }
#else
  MFEM_ABORT("Wave port VoltagePath computation requires MFEM_USE_GSLIB!");
#endif
  return V;
}

std::complex<double> WavePortData::GetExcitationVoltage() const
{
  if (!has_voltage_coords)
  {
    return 0.0;
  }
  // Transfer the port mode tangential E-field from the 2D port submesh back to the 3D
  // parent mesh, then compute the line integral along the voltage path using GSLIB
  // interpolation. Zero the parent field first since SubMeshToParent only writes the
  // mapped DOFs (boundary face DOFs corresponding to the port submesh).
  *parent_E0t = 0.0;
  port_nd_transfer_reverse->Transfer(port_E0t->Real(), parent_E0t->Real());
  port_nd_transfer_reverse->Transfer(port_E0t->Imag(), parent_E0t->Imag());
  // On a nonconforming (AMR-refined) parent mesh the reverse transfer populates only the
  // conforming face DOFs of parent_E0t; hanging-node (slave) DOFs on refined faces crossing
  // the voltage path are left unset, so GSLIB interpolates an under-represented field on the
  // refined child elements and the path integral V collapses. Apply the conforming
  // prolongation (true-dof round trip P·R) to fill slave DOFs consistently before the line
  // integral. On a conforming mesh this is a no-op.
  {
    mfem::ParGridFunction &pr = parent_E0t->Real();
    mfem::ParGridFunction &pi = parent_E0t->Imag();
    mfem::Vector tr, ti;
    pr.GetTrueDofs(tr);
    pr.SetFromTrueDofs(tr);
    pi.GetTrueDofs(ti);
    pi.SetFromTrueDofs(ti);
  }
  std::complex<double> V(0.0, 0.0);
#if defined(MFEM_USE_GSLIB)
  // Reuse the cached point locator (Setup once at construction) — the GSLIB spatial hash
  // depends only on the parent mesh, not the transferred field values. (Line integrals
  // require GSLIB regardless; the cached locator just avoids rebuilding the hash per call.)
  for (std::size_t k = 0; k + 1 < voltage_path.size(); k++)
  {
    V.real(V.real() + fem::ComputeLineIntegral(*voltage_gslib_op, voltage_path[k],
                                               voltage_path[k + 1], parent_E0t->Real(),
                                               voltage_n_samples));
    V.imag(V.imag() + fem::ComputeLineIntegral(*voltage_gslib_op, voltage_path[k],
                                               voltage_path[k + 1], parent_E0t->Imag(),
                                               voltage_n_samples));
  }
#else
  MFEM_ABORT("Wave port VoltagePath computation requires MFEM_USE_GSLIB!");
#endif
  return V;
}

std::complex<double> WavePortData::GetCharacteristicImpedance() const
{
  // Z_PV = (V·V*) / P_mode, where P_mode = ∫(E_mode × H_mode*)·n dS is the full
  // complex Poynting integral of the mode field over the port boundary. The driven
  // solver's Normalize() function normalizes the mode so that |P_mode| = 1, with the
  // mode polarity fixed (via VoltagePath when configured) so that P_mode is real-
  // positive for a propagating mode. Therefore P_mode = 1 and Z_PV reduces to V·V*.
  auto V = GetExcitationVoltage();
  if (std::abs(V) == 0.0)
  {
    return 0.0;
  }
  return V * std::conj(V);
}

WavePortOperator::WavePortOperator(const config::BoundaryData &boundaries,
                                   const config::DomainData &domains,
                                   const config::SolverData &solver,
                                   ProblemType problem_type, const Units &units,
                                   const MaterialOperator &mat_op,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::ParFiniteElementSpace &h1_fespace)
  : suppress_output(false), fc(units.Dimensionalize<Units::ValueType::FREQUENCY>(1.0)),
    kc(1.0 / units.Dimensionalize<Units::ValueType::LENGTH>(1.0))
{
  MFEM_VERIFY(nd_fespace.GetParMesh() == h1_fespace.GetParMesh(),
              "Mesh mismatch in WavePortOperator FE spaces!");
  SetUpBoundaryProperties(boundaries, domains, solver, problem_type, units, mat_op,
                          nd_fespace, h1_fespace);
  PrintBoundaryInfo(units, *nd_fespace.GetParMesh());
}

WavePortOperator::WavePortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::ParFiniteElementSpace &h1_fespace)
  : WavePortOperator(iodata.boundaries, iodata.domains, iodata.solver, iodata.problem.type,
                     iodata.units, mat_op, nd_fespace, h1_fespace)
{
}

void WavePortOperator::SetUpBoundaryProperties(const config::BoundaryData &boundaries,
                                               const config::DomainData &domains,
                                               const config::SolverData &solver,
                                               ProblemType problem_type, const Units &units,
                                               const MaterialOperator &mat_op,
                                               mfem::ParFiniteElementSpace &nd_fespace,
                                               mfem::ParFiniteElementSpace &h1_fespace)
{

  // Check that wave port boundary attributes have been specified correctly.
  const auto &mesh = *nd_fespace.GetParMesh();
  MFEM_VERIFY(boundaries.waveport.empty() || mesh.Dimension() == 3,
              "Wave port boundaries are only available for 3D simulations!");
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  if (!boundaries.waveport.empty())
  {
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), port_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    port_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : boundaries.waveport)
    {
      for (auto attr : data.attributes)
      {
        MFEM_VERIFY(attr > 0 && attr <= bdr_attr_max,
                    "Port boundary attribute tags must be non-negative and correspond to "
                    "boundaries in the mesh!");
        MFEM_VERIFY(bdr_attr_marker[attr - 1],
                    "Unknown port boundary attribute " << attr << "!");
        MFEM_VERIFY(!data.active || !port_marker[attr - 1],
                    "Boundary attribute is assigned to more than one wave port!");
        port_marker[attr - 1] = 1;
      }
    }
  }

  // List of all boundaries which will be marked as essential (Dirichlet) for the purposes
  // of computing wave port modes: PEC and AuxPEC surfaces. Other wave ports are also marked
  // as Dirichlet in case two wave ports are touching and share one or more edges.
  // Impedance, conductivity, and absorbing BCs are handled through their respective
  // boundary operators in the eigenvalue problem.
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(boundaries.pec.attributes.size() +
                                   boundaries.auxpec.attributes.size()));
  for (auto attr : boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;
    }
    dbc_bcs.Append(attr);
  }
  for (auto attr : boundaries.auxpec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;
    }
    dbc_bcs.Append(attr);
  }
  // If user accidentally specifies a surface as both "PEC" and "WavePortPEC", this is fine
  // so allow for duplicates in the attribute list.
  dbc_bcs.Sort();
  dbc_bcs.Unique();

  // Set up wave port data structures.
  for (const auto &[idx, data] : boundaries.waveport)
  {
    mfem::Array<int> port_dbc_bcs(dbc_bcs);
    for (const auto &[other_idx, other_data] : boundaries.waveport)
    {
      if (other_idx == idx || !other_data.active)
      {
        continue;
      }
      for (auto attr : other_data.attributes)
      {
        if (std::binary_search(data.attributes.begin(), data.attributes.end(), attr))
        {
          continue;
        }
        port_dbc_bcs.Append(attr);
      }
    }
    port_dbc_bcs.Sort();
    port_dbc_bcs.Unique();
    ports.try_emplace(idx, data, boundaries, domains, problem_type, solver.linear, units,
                      mat_op, nd_fespace, h1_fespace, port_dbc_bcs);
  }
  MFEM_VERIFY(
      ports.empty() || problem_type == ProblemType::DRIVEN ||
          problem_type == ProblemType::EIGENMODE,
      "Wave port boundaries are only available for frequency domain driven simulations!");
}

void WavePortOperator::PrintBoundaryInfo(const Units &units, const mfem::ParMesh &mesh)
{
  if (ports.empty())
  {
    return;
  }
  fmt::memory_buffer buffer{};
  auto out = fmt::appender{buffer};

  // Print out BC info for all active port attributes.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    for (auto attr : data.GetAttrList())
    {
      fmt::format_to(
          out, " {:d}: Index = {:d}, mode = {:d}, d = {:.3e} m,  n = ({:+.1f})\n", attr,
          idx, data.mode_idx, units.Dimensionalize<Units::ValueType::LENGTH>(data.d_offset),
          fmt::join(data.port_normal, ","));
    }
  }
  if (buffer.size() > 0)
  {
    Mpi::Print("\nConfiguring Robin impedance BC for wave ports at attributes:\n");
    Mpi::Print("{}", fmt::to_string(buffer));
    buffer.clear();
  }

  // Print some information for excited wave ports.
  for (const auto &[idx, data] : ports)
  {
    if (!data.HasExcitation())
    {
      continue;
    }
    for (auto attr : data.GetAttrList())
    {
      fmt::format_to(out, " {:d}: Index = {:d}\n", attr, idx);
    }
  }
  if (buffer.size() > 0)
  {
    Mpi::Print("\nConfiguring wave port excitation source term at attributes:\n");
    Mpi::Print("{}", fmt::to_string(buffer));
  }
}

const WavePortData &WavePortOperator::GetPort(int idx) const
{
  auto it = ports.find(idx);
  MFEM_VERIFY(it != ports.end(), "Unknown wave port index requested!");
  return it->second;
}

mfem::Array<int> WavePortOperator::GetAttrList() const
{
  mfem::Array<int> attr_list;
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    attr_list.Append(data.GetAttrList());
  }
  return attr_list;
}

void WavePortOperator::Initialize(double omega)
{
  bool init = false, first = true;
  for (const auto &[idx, data] : ports)
  {
    init = init || (data.omega0 != omega);
    first = first && (data.omega0 == 0.0);
  }
  if (!init)
  {
    return;
  }
  BlockTimer bt(Timer::WAVE_PORT);
  if (!suppress_output)
  {
    Mpi::Print(
        "\nCalculating boundary modes at wave ports for ω/2π = {:.3e} GHz ({:.3e})\n",
        omega * fc / (2.0 * M_PI), omega);
  }
  for (auto &[idx, data] : ports)
  {
    data.Initialize(omega);
    if (!suppress_output)
    {
      if (first)
      {
        Mpi::Print(" Number of global unknowns for port {:d}:\n"
                   "  H1: {:d}, ND: {:d}\n",
                   idx, data.GlobalTrueH1Size(), data.GlobalTrueNDSize());
      }
      if (data.HasVoltageCoords())
      {
        auto Z_pv = data.GetCharacteristicImpedance().real() * electromagnetics::Z0_;
        Mpi::Print(" Port {:d}, mode {:d}: kₙ = {:.3e}{:+.3e}i m⁻¹, Z_PV = {:.3e} Ω\n", idx,
                   data.mode_idx, data.kn0.real() * kc, data.kn0.imag() * kc, Z_pv);
      }
      else
      {
        Mpi::Print(" Port {:d}, mode {:d}: kₙ = {:.3e}{:+.3e}i m⁻¹\n", idx, data.mode_idx,
                   data.kn0.real() * kc, data.kn0.imag() * kc);
      }
    }
  }
}

void WavePortOperator::AddExtraSystemBdrCoefficients(double omega,
                                                     MaterialPropertyCoefficient &fbr,
                                                     MaterialPropertyCoefficient &fbi)
{
  // Real-ω wave-port BC (driven / boundary-mode path). The contribution looks a lot like
  // the lumped port boundary, except the iω / Z_s coefficient goes to i·k_n / μ where k_n
  // is specific to the port mode at the operating frequency. Only the real (propagating)
  // part of k_n is used here: a real-ω driven/boundary-mode solve models a lossless,
  // energy- carrying port line, so the line-attenuation Im(k_n) is intentionally dropped.
  // (The eigenmode nonlinear solve uses the full complex k_n via the complex overload
  // below.)
  Initialize(omega);
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    AddBoundaryMassBdrCoefficients(idx, fbi, data.kn0.real());
  }
}

void WavePortOperator::AddExtraSystemBdrCoefficients(std::complex<double> omega,
                                                     MaterialPropertyCoefficient &fbr,
                                                     MaterialPropertyCoefficient &fbi)
{
  // Complex-ω wave-port BC for the eigenmode nonlinear solve. The system contribution is
  // A2_wp = i·k_n(ω)·M_{μ⁻¹,p} with the EXACT complex k_n(ω) from the cross-section EVP
  // solved at the genuinely complex eigenvalue (ω = -i·λ). For k_n = β + iα, expanding
  // i·(β + iα) = −α + iβ puts the line-attenuation term −α·M on the real slot (fbr) and
  // the propagating term β·M on the imaginary slot (fbi). SolveKnComplex is
  // side-effect-free (does NOT cache, leaving omega0/kn0/port_E0t untouched), so the real-ω
  // driven/postpro state is preserved. For real ω (imag = 0, α from cross-section loss)
  // this matches the physics of the double overload but additionally carries the
  // attenuation on fbr.
  for (auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    const std::complex<double> kn = data.SolveKnComplex(omega);
    AddBoundaryMassBdrCoefficients(idx, fbr, -kn.imag());
    AddBoundaryMassBdrCoefficients(idx, fbi, kn.real());
  }
}

void WavePortOperator::AddBoundaryMassBdrCoefficients(int port_idx,
                                                      MaterialPropertyCoefficient &fb) const
{
  AddBoundaryMassBdrCoefficients(port_idx, fb, 1.0);
}

void WavePortOperator::AddBoundaryMassBdrCoefficients(int port_idx,
                                                      MaterialPropertyCoefficient &fb,
                                                      double scale) const
{
  // Per-port μ⁻¹ boundary mass coefficient with optional scalar scaling. Pulling this out
  // of AddExtraSystemBdrCoefficients gives the reduced-order model access to the
  // ω-independent operator separately from its per-ω scaling k_n(ω).
  auto it = ports.find(port_idx);
  if (it == ports.end())
  {
    return;
  }
  // This helper exposes the unit boundary mass independently of whether the Robin
  // termination is active. AddExtraSystemBdrCoefficients performs the Active check before
  // calling it, while circuit synthesis also needs the mass of an inactive-but-included
  // port to normalize its node and compute reference data.
  const auto &data = it->second;
  const MaterialOperator &mat_op = data.mat_op;
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  muinv_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(data.GetAttrList()));
  fb.AddCoefficient(muinv_func.GetAttributeToMaterial(), muinv_func.GetMaterialProperties(),
                    scale);
}

double WavePortOperator::GetWavePortKn(int port_idx, double omega)
{
  auto it = ports.find(port_idx);
  MFEM_VERIFY(it != ports.end(),
              "GetWavePortKn called with unknown port index " << port_idx << "!");
  it->second.Initialize(omega);
  return it->second.kn0.real();
}

std::complex<double> WavePortOperator::GetWavePortKnComplex(int port_idx,
                                                            std::complex<double> omega)
{
  auto it = ports.find(port_idx);
  MFEM_VERIFY(it != ports.end(),
              "GetWavePortKnComplex called with unknown port index " << port_idx << "!");
  return it->second.SolveKnComplex(omega);
}

void WavePortOperator::AddExcitationBdrCoefficients(int excitation_idx, double omega,
                                                    SumVectorCoefficient &fbr,
                                                    SumVectorCoefficient &fbi)
{
  // Re/Im{-U_inc} = Re/Im{+2 (-iω) n x H_inc}, which is a function of E_inc as computed by
  // the modal solution (stored as a grid function and coefficient during initialization).
  Initialize(omega);
  for (const auto &[idx, data] : ports)
  {
    if (data.excitation != excitation_idx)
    {
      continue;
    }
    fbr.AddCoefficient(data.GetModeExcitationCoefficientImag(), 2.0 * omega);
    fbi.AddCoefficient(data.GetModeExcitationCoefficientReal(), -2.0 * omega);
  }
}

namespace
{

// Assemble s = ∫_Γ φ·(n×H_mode) on the parent ND true-dof space from an n×H coefficient
// pair, the same n×H the excitation injects (GetModeExcitationCoefficient), so the enforced
// BC, excitation, normalization, and S-projection share one modal calibration. Essential
// (PEC) dofs are eliminated to match the system operator.
std::unique_ptr<ComplexVector> AssembleNxHVector(FiniteElementSpace &nd_fespace,
                                                 const mfem::Array<int> &nd_dbc_tdof_list,
                                                 mfem::VectorCoefficient &cr,
                                                 mfem::VectorCoefficient &ci)
{
  auto s = std::make_unique<ComplexVector>(nd_fespace.GetTrueVSize());
  s->UseDevice(true);
  *s = 0.0;
  for (auto [c, tv] : {std::pair{&cr, &s->Real()}, std::pair{&ci, &s->Imag()}})
  {
    mfem::LinearForm lf(&nd_fespace.Get());
    lf.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(*c));
    lf.UseFastAssembly(false);
    lf.UseDevice(false);
    lf.Assemble();
    lf.UseDevice(true);
    nd_fespace.GetProlongationMatrix()->AddMultTranspose(lf, *tv);
  }
  linalg::SetSubVector(s->Real(), nd_dbc_tdof_list, 0.0);
  linalg::SetSubVector(s->Imag(), nd_dbc_tdof_list, 0.0);
  return s;
}

}  // namespace

std::vector<WavePortOperator::ModalCorrectionTerm>
WavePortOperator::GetModalCorrectionTerms(double omega, FiniteElementSpace &nd_fespace,
                                          const mfem::Array<int> &nd_dbc_tdof_list)
{
  // Assemble the modal correction W = Σ_ports (W_full − W_scalar) over the active wave
  // ports. Needs the current per-port modes and reactions, so refresh them.
  Initialize(omega);
  std::vector<ModalCorrectionTerm> terms;

  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    // The reaction R = sᵀe is generically nonzero for a propagating mode, but can vanish
    // for an evanescent/degenerate mode. Skip the correction for such a port (falling back
    // to the scalar-admittance baseline) rather than dividing by zero.
    if (!(std::abs(data.modal_reaction) > 0.0) ||
        !(std::abs(data.modal_reaction_scalar) > 0.0))
    {
      Mpi::Warning(
          nd_fespace.GetComm(),
          "Wave port {:d} has zero modal reaction; skipping its modal correction!\n", idx);
      continue;
    }
    // W_port = W_full − W_scalar = (−iω/R_full) s_full s_fullᵀ − (−iω/R_scalar) s_scalar
    // s_scalarᵀ. Added to the sparse local mass i·k_n·M (whose modal action equals W_scalar
    // e_mode), the port sees the full modal n×H on its own mode: a matched port satisfies
    // (i·k_n·M + W_port) e_mode = −iω·s_full, consistent with the −2iω·s_full excitation.
    // For a TEM mode s_full=s_scalar and R_full=R_scalar, so W_port ≡ 0 (baseline
    // preserved).
    auto cr = data.GetModeExcitationCoefficientReal(/*include_gradient=*/true);
    auto ci = data.GetModeExcitationCoefficientImag(/*include_gradient=*/true);
    terms.push_back({AssembleNxHVector(nd_fespace, nd_dbc_tdof_list, *cr, *ci),
                     std::complex<double>(0.0, -omega) / data.modal_reaction});
    auto cr_s = data.GetModeExcitationCoefficientReal(/*include_gradient=*/false);
    auto ci_s = data.GetModeExcitationCoefficientImag(/*include_gradient=*/false);
    terms.push_back({AssembleNxHVector(nd_fespace, nd_dbc_tdof_list, *cr_s, *ci_s),
                     std::complex<double>(0.0, omega) / data.modal_reaction_scalar});
  }
  return terms;
}

std::vector<WavePortOperator::ModalCorrectionTerm>
WavePortOperator::GetModalCorrectionTerms(std::complex<double> omega,
                                          FiniteElementSpace &nd_fespace,
                                          const mfem::Array<int> &nd_dbc_tdof_list)
{
  // Complex-ω modal correction for the eigenmode path. Does not call Initialize (ω is
  // complex): recompute the mode at ω (ComputeComplexReactions), then assemble the full n×H
  // and scalar- admittance n×H shape vectors from that recomputed field. The two rank-1
  // terms of W_full − W_scalar are (−iω/R_full) s_full s_fullᵀ and −(−iω/R_scalar) s_scalar
  // s_scalarᵀ, exactly as the real-ω baseline but at complex ω. Pairing each s with the raw
  // reaction of the same field makes W = (−iω/R) s sᵀ invariant to the mode's arbitrary EVP
  // scale/phase, so W tracks the true mode shape at ω (frequency-independent) instead of
  // freezing it at ω0.
  std::vector<ModalCorrectionTerm> terms;
  for (auto &[idx, data] : ports)  // non-const: ComputeComplexReactions re-solves the EVP
  {
    if (!data.active)
    {
      continue;
    }
    if (!(std::abs(data.modal_reaction) > 0.0) ||
        !(std::abs(data.modal_reaction_scalar) > 0.0) || !(std::abs(data.kn0) > 0.0))
    {
      Mpi::Warning(
          nd_fespace.GetComm(),
          "Wave port {:d} has zero modal reaction; skipping its modal correction!\n", idx);
      continue;
    }
    auto smp = SamplePortModalCorrection(data, omega, nd_fespace, nd_dbc_tdof_list);
    if (!smp.active)
    {
      continue;
    }
    terms.push_back({std::move(smp.s_full), smp.g_full});
    terms.push_back({std::move(smp.s_scalar), smp.g_scalar});
  }
  return terms;
}

WavePortOperator::ModalCorrectionSample
WavePortOperator::SamplePortModalCorrection(WavePortData &data, std::complex<double> omega,
                                            FiniteElementSpace &nd_fespace,
                                            const mfem::Array<int> &nd_dbc_tdof_list)
{
  // Recompute the mode at ω (ComputeComplexReactions stashes the ω-field), then assemble
  // the full n×H and scalar-admittance n×H shape vectors from that recomputed field. The
  // two rank-1 terms of W_full − W_scalar are (−iω/R_full) s_full s_fullᵀ and
  // +(iω/R_scalar) s_scalar s_scalarᵀ. Pairing each s with the raw reaction of the same
  // field makes W = (−iω/R) s sᵀ invariant to the mode's arbitrary EVP scale/phase, so W
  // tracks the true mode shape at ω instead of freezing it at ω0.
  ModalCorrectionSample smp;
  if (!(std::abs(data.modal_reaction) > 0.0) ||
      !(std::abs(data.modal_reaction_scalar) > 0.0) || !(std::abs(data.kn0) > 0.0))
  {
    return smp;
  }
  const auto react = data.ComputeComplexReactions(omega);
  if (!(std::abs(react.R_full_raw) > 0.0) || !(std::abs(react.R_scalar_raw) > 0.0))
  {
    return smp;
  }
  auto cr = data.GetOmegaModeExcitationCoefficientReal(/*include_gradient=*/true);
  auto ci = data.GetOmegaModeExcitationCoefficientImag(/*include_gradient=*/true);
  smp.s_full = AssembleNxHVector(nd_fespace, nd_dbc_tdof_list, *cr, *ci);
  auto cr_s = data.GetOmegaModeExcitationCoefficientReal(/*include_gradient=*/false);
  auto ci_s = data.GetOmegaModeExcitationCoefficientImag(/*include_gradient=*/false);
  smp.s_scalar = AssembleNxHVector(nd_fespace, nd_dbc_tdof_list, *cr_s, *ci_s);
  smp.g_full = std::complex<double>(0.0, -1.0) * omega / react.R_full_raw;
  smp.g_scalar = std::complex<double>(0.0, 1.0) * omega / react.R_scalar_raw;
  smp.active = true;
  return smp;
}

std::unique_ptr<ComplexOperator>
WavePortOperator::GetModalCorrectionOperator(double omega, FiniteElementSpace &nd_fespace,
                                             const mfem::Array<int> &nd_dbc_tdof_list)
{
  auto terms = GetModalCorrectionTerms(omega, nd_fespace, nd_dbc_tdof_list);
  if (terms.empty())
  {
    return nullptr;
  }
  auto op = std::make_unique<WavePortModalCorrection>(nd_fespace.GetComm(),
                                                      nd_fespace.GetTrueVSize());
  for (auto &term : terms)
  {
    op->AddTerm(std::move(term.s), term.g);
  }
  return op;
}

std::unique_ptr<ComplexOperator>
WavePortOperator::GetModalCorrectionOperator(std::complex<double> omega,
                                             FiniteElementSpace &nd_fespace,
                                             const mfem::Array<int> &nd_dbc_tdof_list)
{
  auto terms = GetModalCorrectionTerms(omega, nd_fespace, nd_dbc_tdof_list);
  if (terms.empty())
  {
    return nullptr;
  }
  auto op = std::make_unique<WavePortModalCorrection>(nd_fespace.GetComm(),
                                                      nd_fespace.GetTrueVSize());
  for (auto &term : terms)
  {
    op->AddTerm(std::move(term.s), term.g);
  }
  return op;
}

std::vector<int>
WavePortOperator::GetModalCorrectionSynthesisPorts(double omega_ref,
                                                   FiniteElementSpace &nd_fespace,
                                                   const mfem::Array<int> &nd_dbc_tdof_list)
{
  // Enumerate the ports whose modal correction is active at the band. Initialize once at
  // ω_ref and, per port, skip TE/TEM modes where s_full ≈ s_scalar (Eₙ ≈ 0) so synthesis
  // stays exactly at baseline there. The per-ω vectors are recomputed later by
  // SampleModalCorrectionVectors; this only decides which ports participate.
  Initialize(omega_ref);
  std::vector<int> port_idxs;
  for (auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    if (!(std::abs(data.modal_reaction) > 0.0) ||
        !(std::abs(data.modal_reaction_scalar) > 0.0) || !(std::abs(data.kn0) > 0.0))
    {
      Mpi::Warning(nd_fespace.GetComm(),
                   "Wave port {:d} has zero modal reaction; skipping its modal correction "
                   "in synthesis!\n",
                   idx);
      continue;
    }
    auto cr = data.GetModeExcitationCoefficientReal(/*include_gradient=*/true);
    auto ci = data.GetModeExcitationCoefficientImag(/*include_gradient=*/true);
    auto s_full = AssembleNxHVector(nd_fespace, nd_dbc_tdof_list, *cr, *ci);
    auto cr_s = data.GetModeExcitationCoefficientReal(/*include_gradient=*/false);
    auto ci_s = data.GetModeExcitationCoefficientImag(/*include_gradient=*/false);
    auto s_scalar = AssembleNxHVector(nd_fespace, nd_dbc_tdof_list, *cr_s, *ci_s);
    const double norm_full = linalg::Norml2(nd_fespace.GetComm(), *s_full);
    s_full->Add(std::complex<double>(-1.0, 0.0), *s_scalar);
    if (norm_full == 0.0 ||
        linalg::Norml2(nd_fespace.GetComm(), *s_full) <= 1.0e-9 * norm_full)
    {
      continue;
    }
    port_idxs.push_back(idx);
  }
  return port_idxs;
}

WavePortOperator::ModalCorrectionSample
WavePortOperator::SampleModalCorrectionVectors(int port_idx, std::complex<double> omega,
                                               FiniteElementSpace &nd_fespace,
                                               const mfem::Array<int> &nd_dbc_tdof_list)
{
  auto it = ports.find(port_idx);
  MFEM_VERIFY(it != ports.end(),
              "SampleModalCorrectionVectors called with unknown port index " << port_idx
                                                                             << "!");
  return SamplePortModalCorrection(it->second, omega, nd_fespace, nd_dbc_tdof_list);
}

}  // namespace palace
