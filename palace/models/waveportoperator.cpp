// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "waveportoperator.hpp"

#include <tuple>
#include "fem/bilinearform.hpp"
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/arpack.hpp"
#include "linalg/iterative.hpp"
#include "linalg/mumps.hpp"
#include "linalg/rap.hpp"
#include "linalg/slepc.hpp"
#include "linalg/solver.hpp"
#include "linalg/strumpack.hpp"
#include "linalg/superlu.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/timer.hpp"

namespace palace
{

using namespace std::complex_literals;

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

using ComplexHypreParMatrix = std::tuple<std::unique_ptr<mfem::HypreParMatrix>,
                                         std::unique_ptr<mfem::HypreParMatrix>>;
constexpr bool skip_zeros = false;

ComplexHypreParMatrix GetAtt(const MaterialOperator &mat_op,
                             const FiniteElementSpace &nd_fespace,
                             const mfem::Vector &normal, double omega, double sigma)
{
  // Stiffness matrix (shifted): Aₜₜ = (μ⁻¹ ∇ₜ x u, ∇ₜ x v) - ω² (ε u, v) - σ (μ⁻¹ u, v).
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  muinv_func.NormalProjectedCoefficient(normal);
  MaterialPropertyCoefficient epsilon_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetPermittivityReal(), -omega * omega);
  epsilon_func.AddCoefficient(mat_op.GetBdrAttributeToMaterial(),
                              mat_op.GetInvPermeability(), -sigma);
  BilinearForm attr(nd_fespace);
  attr.AddDomainIntegrator<CurlCurlMassIntegrator>(muinv_func, epsilon_func);

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {ParOperator(attr.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
            nullptr};
  }
  MaterialPropertyCoefficient negepstandelta_func(
      mat_op.GetBdrAttributeToMaterial(), mat_op.GetPermittivityImag(), -omega * omega);
  BilinearForm atti(nd_fespace);
  atti.AddDomainIntegrator<VectorFEMassIntegrator>(negepstandelta_func);
  return {ParOperator(attr.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          ParOperator(atti.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble()};
}

ComplexHypreParMatrix GetAtn(const MaterialOperator &mat_op,
                             const FiniteElementSpace &nd_fespace,
                             const FiniteElementSpace &h1_fespace)
{
  // Coupling matrix: Aₜₙ = -(μ⁻¹ ∇ₜ u, v).
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability(), -1.0);
  BilinearForm atn(h1_fespace, nd_fespace);
  atn.AddDomainIntegrator<MixedVectorGradientIntegrator>(muinv_func);
  return {ParOperator(atn.FullAssemble(skip_zeros), h1_fespace, nd_fespace, false)
              .StealParallelAssemble(),
          nullptr};
}

ComplexHypreParMatrix GetAnt(const MaterialOperator &mat_op,
                             const FiniteElementSpace &h1_fespace,
                             const FiniteElementSpace &nd_fespace)
{
  // Coupling matrix: Aₙₜ = -(ε u, ∇ₜ v).
  MaterialPropertyCoefficient epsilon_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetPermittivityReal(), 1.0);

  BilinearForm antr(nd_fespace, h1_fespace);
  antr.AddDomainIntegrator<MixedVectorWeakDivergenceIntegrator>(epsilon_func);

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {ParOperator(antr.FullAssemble(skip_zeros), nd_fespace, h1_fespace, false)
                .StealParallelAssemble(),
            nullptr};
  }
  MaterialPropertyCoefficient negepstandelta_func(mat_op.GetBdrAttributeToMaterial(),
                                                  mat_op.GetPermittivityImag(), 1.0);
  BilinearForm anti(nd_fespace, h1_fespace);
  anti.AddDomainIntegrator<MixedVectorWeakDivergenceIntegrator>(negepstandelta_func);
  return {ParOperator(antr.FullAssemble(skip_zeros), nd_fespace, h1_fespace, false)
              .StealParallelAssemble(),
          ParOperator(anti.FullAssemble(skip_zeros), nd_fespace, h1_fespace, false)
              .StealParallelAssemble()};
}

ComplexHypreParMatrix GetAnn(const MaterialOperator &mat_op,
                             const FiniteElementSpace &h1_fespace,
                             const mfem::Vector &normal)
{
  // Mass matrix: Aₙₙ = -(ε u, v).
  MaterialPropertyCoefficient epsilon_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetPermittivityReal(), -1.0);
  epsilon_func.NormalProjectedCoefficient(normal);
  BilinearForm annr(h1_fespace);
  annr.AddDomainIntegrator<MassIntegrator>(epsilon_func);

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {ParOperator(annr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble(),
            nullptr};
  }
  MaterialPropertyCoefficient negepstandelta_func(mat_op.GetBdrAttributeToMaterial(),
                                                  mat_op.GetPermittivityImag(), -1.0);
  negepstandelta_func.NormalProjectedCoefficient(normal);
  BilinearForm anni(h1_fespace);
  anni.AddDomainIntegrator<MassIntegrator>(negepstandelta_func);
  return {ParOperator(annr.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble(),
          ParOperator(anni.FullAssemble(skip_zeros), h1_fespace).StealParallelAssemble()};
}

ComplexHypreParMatrix GetBtt(const MaterialOperator &mat_op,
                             const FiniteElementSpace &nd_fespace)
{
  // Mass matrix: Bₜₜ = (μ⁻¹ u, v).
  MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                         mat_op.GetInvPermeability());
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator<VectorFEMassIntegrator>(muinv_func);
  return {ParOperator(btt.FullAssemble(skip_zeros), nd_fespace).StealParallelAssemble(),
          nullptr};
}

ComplexHypreParMatrix
GetSystemMatrixA(const mfem::HypreParMatrix *Attr, const mfem::HypreParMatrix *Atti,
                 const mfem::HypreParMatrix *Atnr, const mfem::HypreParMatrix *Atni,
                 const mfem::HypreParMatrix *Antr, const mfem::HypreParMatrix *Anti,
                 const mfem::HypreParMatrix *Annr, const mfem::HypreParMatrix *Anni,
                 const mfem::Array<int> &dbc_tdof_list)
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = λ B e.
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Attr;
  blocks(0, 1) = Atnr;
  blocks(1, 0) = Antr;
  blocks(1, 1) = Annr;
  std::unique_ptr<mfem::HypreParMatrix> Ar(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Ai;
  if (Atti)
  {
    blocks(0, 0) = Atti;
    blocks(0, 1) = Atni;
    blocks(1, 0) = Anti;
    blocks(1, 1) = Anni;
    Ai.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs not associated with this wave port or constrained by
  // Dirichlet BCs.
  Ar->EliminateBC(dbc_tdof_list, Operator::DIAG_ONE);
  if (Ai)
  {
    Ai->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(Ar), std::move(Ai)};
}

ComplexHypreParMatrix GetSystemMatrixB(const mfem::HypreParMatrix *Bttr,
                                       const mfem::HypreParMatrix *Btti,
                                       const mfem::HypreParMatrix *Dnn,
                                       const mfem::Array<int> &dbc_tdof_list)
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = λ B e.
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = Bttr;
  blocks(0, 1) = nullptr;
  blocks(1, 0) = nullptr;
  blocks(1, 1) = Dnn;
  std::unique_ptr<mfem::HypreParMatrix> Br(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> Bi;
  if (Btti)
  {
    blocks(0, 0) = Btti;
    Bi.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs not associated with this wave port or constrained by
  // Dirichlet BCs.
  Br->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  if (Bi)
  {
    Bi->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(Br), std::move(Bi)};
}

void Normalize(const GridFunction &S0t, GridFunction &E0t, GridFunction &E0n,
               mfem::LinearForm &sr, mfem::LinearForm &si)
{
  // Normalize grid functions to a chosen polarization direction and unit power, |E x H⋆| ⋅
  // n, integrated over the port surface (+n is the direction of propagation). The n x H
  // coefficients are updated implicitly as the only store references to the Et, En grid
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

public:
  BdrSubmeshEVectorCoefficient(const GridFunction &Et, const GridFunction &En,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems)
    : mfem::VectorCoefficient(Et.Real().VectorDim()), Et(Et), En(En), submesh(submesh),
      submesh_parent_elems(submesh_parent_elems)
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
  }
};

// Computes boundary mode n x H, where +n is the direction of wave propagation: n x H =
// -1/(iωμ) (ikₙ Eₜ + ∇ₜ Eₙ), using the tangential and normal electric field component grid
// functions evaluated on the (single-sided) boundary element.
template <ValueType Type>
class BdrSubmeshHVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const GridFunction &Et, &En;
  const MaterialOperator &mat_op;
  const mfem::ParSubMesh &submesh;
  const std::unordered_map<int, int> &submesh_parent_elems;
  mfem::IsoparametricTransformation T_loc;
  std::complex<double> kn;
  double omega;

public:
  BdrSubmeshHVectorCoefficient(const GridFunction &Et, const GridFunction &En,
                               const MaterialOperator &mat_op,
                               const mfem::ParSubMesh &submesh,
                               const std::unordered_map<int, int> &submesh_parent_elems,
                               std::complex<double> kn, double omega)
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

    // Compute Re/Im{-1/i (ikₙ Eₜ + ∇ₜ Eₙ)} (t-gradient evaluated in boundary element).
    double U_data[3];
    mfem::Vector U(U_data, vdim);
    if constexpr (Type == ValueType::REAL)
    {
      Et.Real().GetVectorValue(*T_submesh, ip, U);
      U *= -kn.real();

      double dU_data[3];
      mfem::Vector dU(dU_data, vdim);
      En.Imag().GetGradient(*T_submesh, dU);
      U -= dU;
    }
    else
    {
      Et.Imag().GetVectorValue(*T_submesh, ip, U);
      U *= -kn.real();

      double dU_data[3];
      mfem::Vector dU(dU_data, vdim);
      En.Real().GetGradient(*T_submesh, dU);
      U += dU;
    }

    // Scale by 1/(ωμ) with μ evaluated in the neighboring element.
    V.SetSize(U.Size());
    mat_op.GetInvPermeability(attr).Mult(U, V);
    V *= (1.0 / omega);
  }
};

}  // namespace

WavePortData::WavePortData(const config::WavePortData &data,
                           const config::SolverData &solver, const MaterialOperator &mat_op,
                           mfem::ParFiniteElementSpace &nd_fespace,
                           mfem::ParFiniteElementSpace &h1_fespace,
                           const mfem::Array<int> &dbc_attr)
  : mat_op(mat_op)
{
  mode_idx = data.mode_idx;
  d_offset = data.d_offset;
  excitation = data.excitation;
  active = data.active;
  kn0 = 0.0;
  omega0 = 0.0;

  // Construct the SubMesh.
  MFEM_VERIFY(!data.attributes.empty(), "Wave port boundary found with no attributes!");
  const auto &mesh = *nd_fespace.GetParMesh();
  attr_list.Append(data.attributes.data(), data.attributes.size());
  port_mesh = std::make_unique<Mesh>(std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromBoundary(mesh, attr_list)));
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

  port_nd_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0t.Real(), port_E0t->Real()));
  port_h1_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0n.Real(), port_E0n->Real()));

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

  // The operators for the generalized eigenvalue problem are:
  //                [Aₜₜ  Aₜₙ] [eₜ] = -kₙ² [Bₜₜ  0ₜₙ] [eₜ]
  //                [Aₙₜ  Aₙₙ] [eₙ]        [0ₙₜ  0ₙₙ] [eₙ]
  // for the wave port of the given index. The transformed variables are related to the true
  // field by Eₜ = eₜ and Eₙ = eₙ / ikₙ. We will actually solve the shift-and-inverse
  // problem (A - σ B)⁻¹ B e = λ e, with λ = 1 / (-kₙ² - σ).
  // Reference: Vardapetyan and Demkowicz, Full-wave analysis of dielectric waveguides at a
  //            given frequency, Math. Comput. (2003).
  // See also: Halla and Monk, On the analysis of waveguide modes in an electromagnetic
  //           transmission line, arXiv:2302.11994 (2023).
  const double c_max = mat_op.GetLightSpeedMax().Max();
  MFEM_VERIFY(c_max > 0.0 && c_max < mfem::infinity(),
              "Invalid material speed of light detected in WavePortOperator!");
  mu_eps_min = 1.0 / (c_max * c_max) * 0.5;  // Add a safety factor for minimum propagation
                                             // constant possible
  // mu_eps_min = 0.0;  // Use standard inverse transformation to avoid conditioning issues
  //                    // associated with shift
  std::tie(Atnr, Atni) = GetAtn(mat_op, *port_nd_fespace, *port_h1_fespace);
  std::tie(Antr, Anti) = GetAnt(mat_op, *port_h1_fespace, *port_nd_fespace);
  std::tie(Annr, Anni) = GetAnn(mat_op, *port_h1_fespace, port_normal);
  {
    // The HypreParMatrix constructor from a SparseMatrix on each process does not copy
    // the SparseMatrix data, but that's OK since this Dnn is copied in the block system
    // matrix construction.
    Vector d(port_h1_fespace->GetTrueVSize());
    d.UseDevice(false);  // SparseMatrix constructor uses Vector on host
    d = 0.0;
    mfem::SparseMatrix diag(d);
    auto Dnn = std::make_unique<mfem::HypreParMatrix>(
        port_h1_fespace->GetComm(), port_h1_fespace->Get().GlobalTrueVSize(),
        port_h1_fespace->Get().GetTrueDofOffsets(), &diag);
    auto [Bttr, Btti] = GetBtt(mat_op, *port_nd_fespace);
    auto [Br, Bi] = GetSystemMatrixB(Bttr.get(), Btti.get(), Dnn.get(), port_dbc_tdof_list);
    opB = std::make_unique<ComplexWrapperOperator>(std::move(Br), std::move(Bi));
  }

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

  // Configure the eigenvalue problem solver. As for the full 3D case, the system matrices
  // are in general complex and symmetric. We supply the operators to the solver in
  // shift-inverted form and handle the back-transformation externally.
  if (port_comm != MPI_COMM_NULL)
  {
    // Define the linear solver to be used for solving systems associated with the
    // generalized eigenvalue problem.
    auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(port_comm, data.verbose);
    gmres->SetInitialGuess(false);
    gmres->SetRelTol(data.ksp_tol);
    gmres->SetMaxIter(data.ksp_max_its);
    gmres->SetRestartDim(data.ksp_max_its);
    // gmres->SetPrecSide(GmresSolverBase::PrecSide::RIGHT);

    config::LinearSolverData::Type pc_type = solver.linear.type;
    if (pc_type == config::LinearSolverData::Type::SUPERLU)
    {
#if !defined(MFEM_USE_SUPERLU)
      MFEM_ABORT("Solver was not built with SuperLU_DIST support, please choose a "
                 "different solver!");
#endif
    }
    else if (pc_type == config::LinearSolverData::Type::STRUMPACK ||
             pc_type == config::LinearSolverData::Type::STRUMPACK_MP)
    {
#if !defined(MFEM_USE_STRUMPACK)
      MFEM_ABORT("Solver was not built with STRUMPACK support, please choose a "
                 "different solver!");
#endif
    }
    else if (pc_type == config::LinearSolverData::Type::MUMPS)
    {
#if !defined(MFEM_USE_MUMPS)
      MFEM_ABORT("Solver was not built with MUMPS support, please choose a "
                 "different solver!");
#endif
    }
    else  // Default choice
    {
#if defined(MFEM_USE_SUPERLU)
      pc_type = config::LinearSolverData::Type::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
      pc_type = config::LinearSolverData::Type::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
      pc_type = config::LinearSolverData::Type::MUMPS;
#else
#error "Wave port solver requires building with SuperLU_DIST, STRUMPACK, or MUMPS!"
#endif
    }
    auto pc = std::make_unique<MfemWrapperSolver<ComplexOperator>>(
        [&]() -> std::unique_ptr<mfem::Solver>
        {
          if (pc_type == config::LinearSolverData::Type::SUPERLU)
          {
#if defined(MFEM_USE_SUPERLU)
            auto slu = std::make_unique<SuperLUSolver>(
                port_comm, config::LinearSolverData::SymFactType::DEFAULT, false,
                data.verbose - 1);
            // slu->GetSolver().SetColumnPermutation(mfem::superlu::MMD_AT_PLUS_A);
            return slu;
#endif
          }
          else if (pc_type == config::LinearSolverData::Type::STRUMPACK)
          {
#if defined(MFEM_USE_STRUMPACK)
            auto strumpack = std::make_unique<StrumpackSolver>(
                port_comm, config::LinearSolverData::SymFactType::DEFAULT,
                config::LinearSolverData::CompressionType::NONE, 0.0, 0, 0,
                data.verbose - 1);
            // strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::AMD);
            return strumpack;
#endif
          }
          else if (pc_type == config::LinearSolverData::Type::MUMPS)
          {
#if defined(MFEM_USE_MUMPS)
            auto mumps = std::make_unique<MumpsSolver>(
                port_comm, mfem::MUMPSSolver::UNSYMMETRIC,
                config::LinearSolverData::SymFactType::DEFAULT, 0.0, data.verbose - 1);
            // mumps->SetReorderingStrategy(mfem::MUMPSSolver::AMD);
            return mumps;
#endif
          }
          return {};
        }());
    pc->SetSaveAssembled(false);
    ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));

    // Define the eigenvalue solver.
    constexpr int print = 0;
    config::WavePortData::EigenSolverType type = data.eigen_type;
    if (type == config::WavePortData::EigenSolverType::SLEPC)
    {
#if !defined(PALACE_WITH_SLEPC)
      MFEM_ABORT("Solver was not built with SLEPc support, please choose a "
                 "different solver!");
#endif
    }
    else if (type == config::WavePortData::EigenSolverType::ARPACK)
    {
#if !defined(PALACE_WITH_ARPACK)
      MFEM_ABORT("Solver was not built with ARPACK support, please choose a "
                 "different solver!");
#endif
    }
    else  // Default choice
    {
#if defined(PALACE_WITH_SLEPC)
      type = config::WavePortData::EigenSolverType::SLEPC;
#elif defined(PALACE_WITH_ARPACK)
      type = config::WavePortData::EigenSolverType::ARPACK;
#else
#error "Wave port solver requires building with ARPACK or SLEPc!"
#endif
    }
    if (type == config::WavePortData::EigenSolverType::ARPACK)
    {
#if defined(PALACE_WITH_ARPACK)
      eigen = std::make_unique<arpack::ArpackEPSSolver>(port_comm, print);
#endif
    }
    else  // config::WavePortData::EigenSolverType::SLEPC
    {
#if defined(PALACE_WITH_SLEPC)
      auto slepc = std::make_unique<slepc::SlepcEPSSolver>(port_comm, print);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
      eigen = std::move(slepc);
#endif
    }
    eigen->SetNumModes(mode_idx, std::max(2 * mode_idx + 1, 5));
    eigen->SetTol(data.eig_tol);
    eigen->SetLinearSolver(*ksp);

    // We want to ignore evanescent modes (kₙ with large imaginary component). The
    // eigenvalue 1 / (-kₙ² - σ) of the shifted problem will be a large-magnitude negative
    // real number for an eigenvalue kₙ² with real part close to but not below the cutoff σ.
    eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::SMALLEST_REAL);
  }

  // Configure port mode sign convention: 1ᵀ Re{-n x H} >= 0 on the "upper-right quadrant"
  // of the wave port boundary, in order to deal with symmetry effectively.
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
}

WavePortData::~WavePortData()
{
  // Free the solvers before the communicator on which they are based.
  ksp.reset();
  eigen.reset();
  if (port_comm != MPI_COMM_NULL)
  {
    MPI_Comm_free(&port_comm);
  }
}

void WavePortData::Initialize(double omega)
{
  if (omega == omega0)
  {
    return;
  }

  // Construct matrices and solve the generalized eigenvalue problem for the desired wave
  // port mode. The B matrix is operating frequency-independent and has already been
  // constructed.
  std::unique_ptr<ComplexOperator> opA;
  const double sigma = -omega * omega * mu_eps_min;
  {
    auto [Attr, Atti] = GetAtt(mat_op, *port_nd_fespace, port_normal, omega, sigma);
    auto [Ar, Ai] =
        GetSystemMatrixA(Attr.get(), Atti.get(), Atnr.get(), Atni.get(), Antr.get(),
                         Anti.get(), Annr.get(), Anni.get(), port_dbc_tdof_list);
    opA = std::make_unique<ComplexWrapperOperator>(std::move(Ar), std::move(Ai));
  }

  // Configure and solve the (inverse) eigenvalue problem for the desired boundary mode.
  // Linear solves are preconditioned with the real part of the system matrix (ignore loss
  // tangent).
  std::complex<double> lambda;
  if (port_comm != MPI_COMM_NULL)
  {
    ComplexWrapperOperator opP(opA->Real(), nullptr);  // Non-owning constructor
    ksp->SetOperators(*opA, opP);
    eigen->SetOperators(*opB, *opA, EigenvalueSolver::ScaleType::NONE);
    eigen->SetInitialSpace(v0);
    int num_conv = eigen->Solve();
    MFEM_VERIFY(num_conv >= mode_idx, "Wave port eigensolver did not converge!");
    lambda = eigen->GetEigenvalue(mode_idx - 1);
    // Mpi::Print(port_comm, " ... Wave port eigensolver error = {} (bkwd), {} (abs)\n",
    //            eigen->GetError(mode_idx - 1, EigenvalueSolver::ErrorType::BACKWARD),
    //            eigen->GetError(mode_idx - 1, EigenvalueSolver::ErrorType::ABSOLUTE));
  }
  Mpi::Broadcast(1, &lambda, port_root, port_mesh->GetComm());

  // Extract the eigenmode solution and postprocess. The extracted eigenvalue is λ =
  // 1 / (-kₙ² - σ).
  kn0 = std::sqrt(-sigma - 1.0 / lambda);
  omega0 = omega;

  // Separate the computed field out into eₜ and eₙ and and transform back to true
  // electric field variables: Eₜ = eₜ and Eₙ = eₙ / ikₙ.
  {
    if (port_comm != MPI_COMM_NULL)
    {
      eigen->GetEigenvector(mode_idx - 1, e0);
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
    ComplexVector::AXPBY(1.0 / (1i * kn0), e0nr, e0ni, 0.0, e0nr, e0ni);
    port_E0t->Real().SetFromTrueDofs(e0tr);  // Parallel distribute
    port_E0t->Imag().SetFromTrueDofs(e0ti);
    port_E0n->Real().SetFromTrueDofs(e0nr);
    port_E0n->Imag().SetFromTrueDofs(e0ni);
  }

  // Configure the linear forms for computing S-parameters (projection of the field onto the
  // port mode). Normalize the mode for a chosen polarization direction and unit power,
  // |E x H⋆| ⋅ n, integrated over the port surface (+n is the direction of propagation).
  {
    const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
    BdrSubmeshHVectorCoefficient<ValueType::REAL> port_nxH0r_func(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0, omega0);
    BdrSubmeshHVectorCoefficient<ValueType::IMAG> port_nxH0i_func(
        *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0, omega0);
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
  }
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeExcitationCoefficientReal() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::REAL>>>(
      attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0,
      omega0);
}

std::unique_ptr<mfem::VectorCoefficient>
WavePortData::GetModeExcitationCoefficientImag() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshHVectorCoefficient<ValueType::IMAG>>>(
      attr_list, *port_E0t, *port_E0n, mat_op, port_submesh, submesh_parent_elems, kn0,
      omega0);
}

std::unique_ptr<mfem::VectorCoefficient> WavePortData::GetModeFieldCoefficientReal() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshEVectorCoefficient<ValueType::REAL>>>(
      attr_list, *port_E0t, *port_E0n, port_submesh, submesh_parent_elems);
}

std::unique_ptr<mfem::VectorCoefficient> WavePortData::GetModeFieldCoefficientImag() const
{
  const auto &port_submesh = static_cast<const mfem::ParSubMesh &>(port_mesh->Get());
  return std::make_unique<
      RestrictedVectorCoefficient<BdrSubmeshEVectorCoefficient<ValueType::IMAG>>>(
      attr_list, *port_E0t, *port_E0n, port_submesh, submesh_parent_elems);
}

double WavePortData::GetExcitationPower() const
{
  // The computed port modes are normalized such that the power integrated over the port is
  // 1: ∫ (E_inc x H_inc⋆) ⋅ n dS = 1.
  return excitation ? 1.0 : 0.0;
}

std::complex<double> WavePortData::GetPower(GridFunction &E, GridFunction &B) const
{
  // Compute port power, (E x H) ⋅ n = E ⋅ (-n x H), integrated over the port surface using
  // the computed E and H = μ⁻¹ B fields, where +n is the direction of propagation (into the
  // domain). The BdrSurfaceCurrentVectorCoefficient computes -n x H for an outward normal,
  // so we multiply by -1. The linear form is reconstructed from scratch each time due to
  // changing H.
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

WavePortOperator::WavePortOperator(const IoData &iodata, const MaterialOperator &mat_op,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::ParFiniteElementSpace &h1_fespace)
  : suppress_output(false),
    fc(iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, 1.0)),
    kc(1.0 / iodata.DimensionalizeValue(IoData::ValueType::LENGTH, 1.0))
{
  // Set up wave port boundary conditions.
  MFEM_VERIFY(nd_fespace.GetParMesh() == h1_fespace.GetParMesh(),
              "Mesh mismatch in WavePortOperator FE spaces!");
  SetUpBoundaryProperties(iodata, mat_op, nd_fespace, h1_fespace);
  PrintBoundaryInfo(iodata, *nd_fespace.GetParMesh());
}

void WavePortOperator::SetUpBoundaryProperties(const IoData &iodata,
                                               const MaterialOperator &mat_op,
                                               mfem::ParFiniteElementSpace &nd_fespace,
                                               mfem::ParFiniteElementSpace &h1_fespace)
{
  // Check that wave port boundary attributes have been specified correctly.
  const auto &mesh = *nd_fespace.GetParMesh();
  int bdr_attr_max = mesh.bdr_attributes.Size() ? mesh.bdr_attributes.Max() : 0;
  if (!iodata.boundaries.waveport.empty())
  {
    mfem::Array<int> bdr_attr_marker(bdr_attr_max), port_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    port_marker = 0;
    for (auto attr : mesh.bdr_attributes)
    {
      bdr_attr_marker[attr - 1] = 1;
    }
    for (const auto &[idx, data] : iodata.boundaries.waveport)
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

  // List of all boundaries which will be marked as essential for the purposes of computing
  // wave port modes. This includes all PEC surfaces, but may also include others like when
  // a kinetic inductance or other BC is applied for the 3D simulation but should be
  // considered as PEC for the 2D problem. In addition, we mark as Dirichlet boundaries all
  // wave ports other than the wave port being currently considered, in case two wave ports
  // are touching and share one or more edges.
  mfem::Array<int> dbc_bcs;
  dbc_bcs.Reserve(static_cast<int>(iodata.boundaries.pec.attributes.size() +
                                   iodata.boundaries.auxpec.attributes.size() +
                                   iodata.boundaries.conductivity.size()));
  for (auto attr : iodata.boundaries.pec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  for (auto attr : iodata.boundaries.auxpec.attributes)
  {
    if (attr <= 0 || attr > bdr_attr_max)
    {
      continue;  // Can just ignore if wrong
    }
    dbc_bcs.Append(attr);
  }
  for (const auto &data : iodata.boundaries.conductivity)
  {
    for (auto attr : data.attributes)
    {
      if (attr <= 0 || attr > bdr_attr_max)
      {
        continue;  // Can just ignore if wrong
      }
      dbc_bcs.Append(attr);
    }
  }
  // If user accidentally specifies a surface as both "PEC" and "WavePortPEC", this is fine
  // so allow for duplicates in the attribute list.
  dbc_bcs.Sort();
  dbc_bcs.Unique();

  // Set up wave port data structures.
  for (const auto &[idx, data] : iodata.boundaries.waveport)
  {
    mfem::Array<int> port_dbc_bcs(dbc_bcs);
    for (const auto &[other_idx, other_data] : iodata.boundaries.waveport)
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
    ports.try_emplace(idx, data, iodata.solver, mat_op, nd_fespace, h1_fespace,
                      port_dbc_bcs);
  }
  MFEM_VERIFY(
      ports.empty() || iodata.problem.type == config::ProblemData::Type::DRIVEN,
      "Wave port boundaries are only available for frequency domain driven simulations!");
}

void WavePortOperator::PrintBoundaryInfo(const IoData &iodata, const mfem::ParMesh &mesh)
{
  if (ports.empty())
  {
    return;
  }
  fmt::memory_buffer buf{};  // Output buffer & buffer append lambda for cleaner code
  auto to = [&buf](auto fmt, auto &&...args)
  { fmt::format_to(std::back_inserter(buf), fmt, std::forward<decltype(args)>(args)...); };
  using VT = IoData::ValueType;

  // Print out BC info for all active port attributes.
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    for (auto attr : data.GetAttrList())
    {
      to(" {:d}: Index = {:d}, mode = {:d}, d = {:.3e} m,  n = ({:+.1f})\n", attr, idx,
         data.mode_idx, iodata.DimensionalizeValue(VT::LENGTH, data.d_offset),
         fmt::join(data.port_normal, ","));
    }
  }
  if (buf.size() > 0)
  {
    Mpi::Print("\nConfiguring Robin impedance BC for wave ports at attributes:\n");
    Mpi::Print("{}", fmt::to_string(buf));
    buf.clear();
  }

  // Print some information for excited wave ports.
  for (const auto &[idx, data] : ports)
  {
    if (!data.excitation)
    {
      continue;
    }
    for (auto attr : data.GetAttrList())
    {
      to(" {:d}: Index = {:d}\n", attr, idx);
    }
  }
  if (buf.size() > 0)
  {
    Mpi::Print("\nConfiguring wave port excitation source term at attributes:\n");
    Mpi::Print("{}", fmt::to_string(buf));
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
        omega * fc, omega);
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
      Mpi::Print(" Port {:d}, mode {:d}: kₙ = {:.3e}{:+.3e}i m⁻¹\n", idx, data.mode_idx,
                 data.kn0.real() * kc, data.kn0.imag() * kc);
    }
  }
}

void WavePortOperator::AddExtraSystemBdrCoefficients(double omega,
                                                     MaterialPropertyCoefficient &fbr,
                                                     MaterialPropertyCoefficient &fbi)
{
  // Add wave port boundaries to the bilinear form. This looks a lot like the lumped port
  // boundary, except the iω / Z_s coefficient goes to ikₙ / μ where kₙ is specific to the
  // port mode at the given operating frequency (note only the real part of the propagation
  // constant contributes).
  Initialize(omega);
  for (const auto &[idx, data] : ports)
  {
    if (!data.active)
    {
      continue;
    }
    const MaterialOperator &mat_op = data.mat_op;
    MaterialPropertyCoefficient muinv_func(mat_op.GetBdrAttributeToMaterial(),
                                           mat_op.GetInvPermeability());
    muinv_func.RestrictCoefficient(mat_op.GetCeedBdrAttributes(data.GetAttrList()));
    // fbr.AddCoefficient(muinv_func.GetAttributeToMaterial(),
    //                    muinv_func.GetMaterialProperties(),
    //                    -data.kn0.imag());
    fbi.AddCoefficient(muinv_func.GetAttributeToMaterial(),
                       muinv_func.GetMaterialProperties(), data.kn0.real());
  }
}

void WavePortOperator::AddExcitationBdrCoefficients(double omega, SumVectorCoefficient &fbr,
                                                    SumVectorCoefficient &fbi)
{
  // Re/Im{-U_inc} = Re/Im{+2 (-iω) n x H_inc}, which is a function of E_inc as computed by
  // the modal solution (stored as a grid function and coefficient during initialization).
  Initialize(omega);
  for (const auto &[idx, data] : ports)
  {
    if (!data.excitation)
    {
      continue;
    }
    fbr.AddCoefficient(data.GetModeExcitationCoefficientImag(), 2.0 * omega);
    fbi.AddCoefficient(data.GetModeExcitationCoefficientReal(), -2.0 * omega);
  }
}

}  // namespace palace
