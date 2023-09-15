// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "waveportoperator.hpp"

#include <array>
#include <unordered_map>
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
                          const mfem::Array<int> &dbc_marker,
                          mfem::Array<int> &port_nd_dbc_tdof_list,
                          mfem::Array<int> &port_h1_dbc_tdof_list)
{
  mfem::ParFiniteElementSpace &nd_fespace = *E0t.ParFESpace();
  mfem::ParFiniteElementSpace &h1_fespace = *E0n.ParFESpace();
  mfem::ParFiniteElementSpace &port_nd_fespace = *port_E0t.ParFESpace();
  mfem::ParFiniteElementSpace &port_h1_fespace = *port_E0n.ParFESpace();

  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  nd_fespace.GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);

  Vector tE0t(nd_fespace.GetTrueVSize()), tE0n(h1_fespace.GetTrueVSize());
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
  port_E0t.ParallelProject(port_tE0t);
  port_E0n.ParallelProject(port_tE0n);
  for (int i = 0; i < port_tE0t.Size(); i++)
  {
    if (port_tE0t[i] != 0.0)
    {
      port_nd_dbc_tdof_list.Append(i);
    }
  }
  for (int i = 0; i < port_tE0n.Size(); i++)
  {
    if (port_tE0n[i] != 0.0)
    {
      port_h1_dbc_tdof_list.Append(i);
    }
  }
}

constexpr bool skip_zeros = false;

std::unique_ptr<ParOperator> GetBtt(const MaterialOperator &mat_op,
                                    const mfem::ParFiniteElementSpace &nd_fespace)
{
  // Mass matrix: Bₜₜ = (μ⁻¹ u, v).
  constexpr auto MatType = MaterialPropertyType::INV_PERMEABILITY;
  constexpr auto ElemType = MeshElementType::BDR_SUBMESH;
  MaterialPropertyCoefficient<MatType, ElemType> muinv_func(mat_op);
  BilinearForm btt(nd_fespace);
  btt.AddDomainIntegrator(std::make_unique<VectorFEMassIntegrator>(muinv_func));
  return std::make_unique<ParOperator>(btt.FullAssemble(skip_zeros), nd_fespace);
}

std::unique_ptr<ParOperator> GetBtn(const MaterialOperator &mat_op,
                                    const mfem::ParFiniteElementSpace &nd_fespace,
                                    const mfem::ParFiniteElementSpace &h1_fespace)
{
  // Mass matrix: Bₜₙ = (μ⁻¹ ∇ₜ u, v).
  constexpr auto MatType = MaterialPropertyType::INV_PERMEABILITY;
  constexpr auto ElemType = MeshElementType::BDR_SUBMESH;
  MaterialPropertyCoefficient<MatType, ElemType> muinv_func(mat_op);
  BilinearForm btn(h1_fespace, nd_fespace);
  btn.AddDomainIntegrator(std::make_unique<MixedVectorGradientIntegrator>(muinv_func));
  return std::make_unique<ParOperator>(btn.FullAssemble(skip_zeros), h1_fespace, nd_fespace,
                                       false);
}

std::array<std::unique_ptr<ParOperator>, 3>
GetBnn(const MaterialOperator &mat_op, const mfem::ParFiniteElementSpace &h1_fespace)
{
  // Mass matrix: Bₙₙ = (μ⁻¹ ∇ₜ u, ∇ₜ v) - ω² (ε u, v) = Bₙₙ₁ - ω² Bₙₙ₂.
  constexpr auto MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  constexpr auto ElemType = MeshElementType::BDR_SUBMESH;
  MaterialPropertyCoefficient<MatTypeMuInv, ElemType> muinv_func(mat_op);
  BilinearForm bnn1(h1_fespace);
  bnn1.AddDomainIntegrator(std::make_unique<DiffusionIntegrator>(muinv_func));

  constexpr auto MatTypeEpsReal = MaterialPropertyType::PERMITTIVITY_REAL;
  NormalProjectedCoefficient epsilon_func(
      std::make_unique<MaterialPropertyCoefficient<MatTypeEpsReal, ElemType>>(mat_op));
  BilinearForm bnn2r(h1_fespace);
  bnn2r.AddDomainIntegrator(std::make_unique<MassIntegrator>(epsilon_func));

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {std::make_unique<ParOperator>(bnn1.FullAssemble(skip_zeros), h1_fespace),
            std::make_unique<ParOperator>(bnn2r.FullAssemble(skip_zeros), h1_fespace),
            nullptr};
  }
  constexpr auto MatTypeEpsImag = MaterialPropertyType::PERMITTIVITY_IMAG;
  NormalProjectedCoefficient negepstandelta_func(
      std::make_unique<MaterialPropertyCoefficient<MatTypeEpsImag, ElemType>>(mat_op));
  BilinearForm bnn2i(h1_fespace);
  bnn2i.AddDomainIntegrator(std::make_unique<MassIntegrator>(negepstandelta_func));
  return {std::make_unique<ParOperator>(bnn1.FullAssemble(skip_zeros), h1_fespace),
          std::make_unique<ParOperator>(bnn2r.FullAssemble(skip_zeros), h1_fespace),
          std::make_unique<ParOperator>(bnn2i.FullAssemble(skip_zeros), h1_fespace)};
}

std::array<std::unique_ptr<ParOperator>, 3>
GetAtt(const MaterialOperator &mat_op, const mfem::ParFiniteElementSpace &nd_fespace)
{
  // Stiffness matrix: Aₜₜ = (μ⁻¹ ∇ₜ x u, ∇ₜ x v) - ω² (ε u, v) = Aₜₜ₁ - ω² Aₜₜ₂.
  constexpr auto MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  constexpr auto ElemType = MeshElementType::BDR_SUBMESH;
  NormalProjectedCoefficient muinv_func(
      std::make_unique<MaterialPropertyCoefficient<MatTypeMuInv, ElemType>>(mat_op));
  BilinearForm att1(nd_fespace);
  att1.AddDomainIntegrator(std::make_unique<CurlCurlIntegrator>(muinv_func));

  constexpr auto MatTypeEpsReal = MaterialPropertyType::PERMITTIVITY_REAL;
  MaterialPropertyCoefficient<MatTypeEpsReal, ElemType> epsilon_func(mat_op);
  BilinearForm att2r(nd_fespace);
  att2r.AddDomainIntegrator(std::make_unique<VectorFEMassIntegrator>(epsilon_func));

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {std::make_unique<ParOperator>(att1.FullAssemble(skip_zeros), nd_fespace),
            std::make_unique<ParOperator>(att2r.FullAssemble(skip_zeros), nd_fespace),
            nullptr};
  }
  constexpr auto MatTypeEpsImag = MaterialPropertyType::PERMITTIVITY_IMAG;
  MaterialPropertyCoefficient<MatTypeEpsImag, ElemType> negepstandelta_func(mat_op);
  BilinearForm att2i(nd_fespace);
  att2i.AddDomainIntegrator(std::make_unique<VectorFEMassIntegrator>(negepstandelta_func));
  return {std::make_unique<ParOperator>(att1.FullAssemble(skip_zeros), nd_fespace),
          std::make_unique<ParOperator>(att2r.FullAssemble(skip_zeros), nd_fespace),
          std::make_unique<ParOperator>(att2i.FullAssemble(skip_zeros), nd_fespace)};
}

std::array<std::unique_ptr<mfem::HypreParMatrix>, 6>
GetSystemMatrices(std::unique_ptr<ParOperator> Btt, std::unique_ptr<ParOperator> Btn,
                  std::unique_ptr<ParOperator> Bnn1, std::unique_ptr<ParOperator> Bnn2r,
                  std::unique_ptr<ParOperator> Bnn2i, std::unique_ptr<ParOperator> Att1,
                  std::unique_ptr<ParOperator> Att2r, std::unique_ptr<ParOperator> Att2i,
                  const mfem::Array<int> &nd_dbc_tdof_list,
                  const mfem::Array<int> &h1_dbc_tdof_list)
{
  // Construct the 2x2 block matrices for the eigenvalue problem A e = λ B e. We pre-compute
  // the matrices such that:
  //              A = A₁ - ω² A₂, B = A₁ - ω² A₂ + 1/Θ² B₃ - ω²/Θ² B₄.
  std::unique_ptr<mfem::HypreParMatrix> BtnT(Btn->ParallelAssemble().Transpose());

  mfem::Array2D<mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = &Btt->ParallelAssemble();
  blocks(0, 1) = &Btn->ParallelAssemble();
  blocks(1, 0) = BtnT.get();
  blocks(1, 1) = &Bnn1->ParallelAssemble();
  std::unique_ptr<mfem::HypreParMatrix> A1(mfem::HypreParMatrixFromBlocks(blocks));

  auto &Ztt = Btt->ParallelAssemble();
  Ztt *= 0.0;

  blocks = nullptr;
  blocks(0, 0) = &Ztt;
  blocks(1, 1) = &Bnn2r->ParallelAssemble();
  std::unique_ptr<mfem::HypreParMatrix> A2r(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> A2i;
  if (Bnn2i)
  {
    blocks(1, 1) = &Bnn2i->ParallelAssemble();
    A2i.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  auto &Znn = Bnn1->ParallelAssemble();
  Znn *= 0.0;

  blocks = nullptr;
  blocks(0, 0) = &Att1->ParallelAssemble();
  blocks(1, 1) = &Znn;
  std::unique_ptr<mfem::HypreParMatrix> B3(mfem::HypreParMatrixFromBlocks(blocks));

  blocks(0, 0) = &Att2r->ParallelAssemble();
  blocks(1, 1) = &Znn;
  std::unique_ptr<mfem::HypreParMatrix> B4r(mfem::HypreParMatrixFromBlocks(blocks));

  std::unique_ptr<mfem::HypreParMatrix> B4i;
  if (Att2i)
  {
    blocks(0, 0) = &Att2i->ParallelAssemble();
    B4i.reset(mfem::HypreParMatrixFromBlocks(blocks));
  }

  // Eliminate boundary true dofs not associated with this wave port or constrained by
  // Dirichlet BCs. It is not guaranteed that any HypreParMatrix has a full diagonal in its
  // sparsity pattern, so we add a zero diagonal before elimination to guarantee this for A1
  // and B3.
  mfem::Array<int> dbc_tdof_list;
  int nd_tdof_offset = Btt->Height();
  dbc_tdof_list.Reserve(nd_dbc_tdof_list.Size() + h1_dbc_tdof_list.Size());
  for (auto tdof : nd_dbc_tdof_list)
  {
    dbc_tdof_list.Append(tdof);
  }
  for (auto tdof : h1_dbc_tdof_list)
  {
    dbc_tdof_list.Append(tdof + nd_tdof_offset);
  }

  mfem::Vector d(B3->Height());
  d = 0.0;
  mfem::SparseMatrix diag(d);
  mfem::HypreParMatrix Diag(B3->GetComm(), B3->GetGlobalNumRows(), B3->GetRowStarts(),
                            &diag);
  A1.reset(mfem::Add(1.0, *A1, 1.0, Diag));
  B3.reset(mfem::Add(1.0, *B3, 1.0, Diag));

  A1->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  A2r->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  if (A2i)
  {
    A2i->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }
  B3->EliminateBC(dbc_tdof_list, Operator::DIAG_ONE);
  B4r->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  if (B4i)
  {
    B4i->EliminateBC(dbc_tdof_list, Operator::DIAG_ZERO);
  }

  return {std::move(A1), std::move(A2r), std::move(A2i),
          std::move(B3), std::move(B4r), std::move(B4i)};
}

void GetInitialSpace(const mfem::ParFiniteElementSpace &nd_fespace,
                     const mfem::ParFiniteElementSpace &h1_fespace,
                     const mfem::Array<int> &nd_dbc_tdof_list,
                     const mfem::Array<int> &h1_dbc_tdof_list, ComplexVector &v)
{
  // Initial space chosen as such that B v₀ = y₀, with y₀ = [y₀ₜ, 0, ... 0]ᵀ ⟂ null(A)
  // (with Aₜₜ nonsingular). See Lee, Sun, and Cendes, 1991 for reference.
  // Note: When the eigenvalue solver uses a standard ℓ²-inner product instead of B-inner
  // product (since we use a general non-Hermitian solver due to complex symmetric B), then
  // we just use v0 = y0 directly.
  v.SetSize(nd_fespace.GetTrueVSize() + h1_fespace.GetTrueVSize());
  // linalg::SetRandomReal(nd_fespace.GetComm(), v);
  v = std::complex<double>(1.0, 0.0);
  linalg::SetSubVector(v, nd_dbc_tdof_list, 0.0);
  for (int i = nd_fespace.GetTrueVSize();
       i < nd_fespace.GetTrueVSize() + h1_fespace.GetTrueVSize(); i++)
  {
    v.Real()[i] = v.Imag()[i] = 0.0;
  }
}

void NormalizeWithSign(const mfem::ParGridFunction &S0t, mfem::ParComplexGridFunction &E0t,
                       mfem::ParComplexGridFunction &E0n, mfem::LinearForm &sr,
                       mfem::LinearForm &si)
{
  // Normalize grid functions to a chosen polarization direction and unit power, |E x H⋆| ⋅
  // n, integrated over the port surface (+n is the direction of propagation). The n x H
  // coefficients are updated implicitly as the only store references to the Et, En grid
  // functions as well as kₙ, ω. We choose a (rather arbitrary) sign constraint to at least
  // make results for the same port consistent between frequencies/meshes.
  sr = 0.0;
  si = 0.0;
  sr.Assemble();
  si.Assemble();

  // |E x H⋆| ⋅ n = |E ⋅ (-n x H⋆)|
  double sign = sr * S0t;
  std::complex<double> dot(-(sr * E0t.real()) - (si * E0t.imag()),
                           -(sr * E0t.imag()) + (si * E0t.real()));
  std::array<double, 3> data = {sign, dot.real(), dot.imag()};
  Mpi::GlobalSum(3, data.data(), S0t.ParFESpace()->GetComm());
  sign = (data[0] < 0.0) ? -1.0 : 1.0;
  dot = {data[1], data[2]};

  double scale = sign / std::sqrt(std::abs(dot));
  E0t.real() *= scale;  // Updates the n x H coefficients depending on Et, En too
  E0t.imag() *= scale;
  E0n.real() *= scale;
  E0n.imag() *= scale;
  sr *= scale;  // Update linear forms for postprocessing
  si *= scale;

  // This parallel communication is not required since wave port boundaries are true
  // one-sided boundaries.
  // port_E0t->real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
  // port_E0t->imag().ExchangeFaceNbrData();  // for n x H coefficients evaluation
  // port_E0n->real().ExchangeFaceNbrData();
  // port_E0n->imag().ExchangeFaceNbrData();
}

// Computes boundary modal n x H, where +n is the direction of wave propagation: n x H =
// -1/(iωμ) (ikₙ Eₜ + ∇ₜ Eₙ), using the tangential and normal electric field component grid
// functions evaluated on the (single-sided) boundary element. The intent of this vector
// grid function is to be dotted with a function E which is only in the tangential
// component, so the fact that we use the full ∇ Eₙ in the element is fine. We use only the
// real part of kn.
template <bool RealPart>
class BdrSubmeshHVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const mfem::ParComplexGridFunction &Et, &En;
  const MaterialOperator &mat_op;

  mfem::ParSubMesh &submesh;
  const mfem::ParMesh &parent;
  std::unordered_map<int, int> submesh_elem_ids;

  std::complex<double> kn;
  double omega;

  mfem::ParSubMesh &GetSubMesh(mfem::ParMesh &mesh)
  {
    MFEM_ASSERT(
        mfem::ParSubMesh::IsParSubMesh(&mesh),
        "BdrSubmeshHVectorCoefficient requires the input grid function coefficients "
        "to be defined on a SubMesh!");
    mfem::ParSubMesh &submesh = *static_cast<mfem::ParSubMesh *>(&mesh);
    MFEM_ASSERT(submesh.GetFrom() == mfem::SubMesh::From::Boundary,
                "BdrSubmeshHVectorCoefficient requires a SubMesh created using "
                "CreateFromBoundary!");
    return submesh;
  }

public:
  BdrSubmeshHVectorCoefficient(const mfem::ParComplexGridFunction &Et,
                               const mfem::ParComplexGridFunction &En,
                               const MaterialOperator &mat_op)
    : mfem::VectorCoefficient(Et.ParFESpace()->GetParMesh()->SpaceDimension()), Et(Et),
      En(En), mat_op(mat_op), submesh(GetSubMesh(*Et.ParFESpace()->GetParMesh())),
      parent(*submesh.GetParent()), kn(0.0), omega(0.0)
  {
    // Construct mapping from parent (boundary) element indices to submesh (domain)
    // elements.
    const mfem::Array<int> &parent_element_ids = submesh.GetParentElementIDMap();
    for (int i = 0; i < parent_element_ids.Size(); i++)
    {
      submesh_elem_ids[parent_element_ids[i]] = i;
    }
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    mfem::ElementTransformation *submesh_T = nullptr;
    int attr = 0;
    if (T.mesh == &parent)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                  "BdrSubmeshHVectorCoefficient requires ElementType::BDR_ELEMENT when not "
                  "used on a SubMesh!");
      auto it = submesh_elem_ids.find(T.ElementNo);
      if (it == submesh_elem_ids.end())
      {
        // Just return zero for a boundary face not in the submesh.
        V.SetSize(vdim);
        V = 0.0;
        return;
      }
      else
      {
        submesh_T = submesh.GetElementTransformation(it->second);
      }

      int i, o, iel1, iel2;
      parent.GetBdrElementFace(T.ElementNo, &i, &o);
      parent.GetFaceElements(i, &iel1, &iel2);
      attr = parent.GetAttribute(iel1);
    }
    else if (T.mesh == &submesh)
    {
      MFEM_ASSERT(T.ElementType == mfem::ElementTransformation::ELEMENT,
                  "BdrSubmeshHVectorCoefficient requires ElementType::ELEMENT when used on "
                  "a SubMesh!");
      submesh_T = &T;

      int i, o, iel1, iel2;
      parent.GetBdrElementFace(submesh.GetParentElementIDMap()[T.ElementNo], &i, &o);
      parent.GetFaceElements(i, &iel1, &iel2);
      attr = parent.GetAttribute(iel1);
    }
    else
    {
      MFEM_ABORT("Invalid use of BdrSubmeshHVectorCoefficient on an unrecognized mesh!");
    }

    // Compute Re/Im{-1/i (ikₙ Eₜ + ∇ₜ Eₙ)}.
    mfem::Vector U;
    submesh_T->SetIntPoint(&ip);
    if constexpr (RealPart)
    {
      Et.real().GetVectorValue(*submesh_T, ip, U);
      U *= -kn.real();

      mfem::Vector dU;
      En.imag().GetGradient(*submesh_T, dU);
      U -= dU;
    }
    else
    {
      Et.imag().GetVectorValue(*submesh_T, ip, U);
      U *= -kn.real();

      mfem::Vector dU;
      En.real().GetGradient(*submesh_T, dU);
      U += dU;
    }

    // Scale by 1/(ωμ) with μ evaluated in the neighboring element.
    V.SetSize(U.Size());
    mat_op.GetInvPermeability(attr).Mult(U, V);
    V *= (1.0 / omega);
  }

  void SetFrequency(double w, std::complex<double> k)
  {
    omega = w;
    kn = k;
  }
};

}  // namespace

WavePortData::WavePortData(const config::WavePortData &data, const MaterialOperator &mat_op,
                           const mfem::ParFiniteElementSpace &nd_fespace,
                           const mfem::ParFiniteElementSpace &h1_fespace,
                           const mfem::Array<int> &dbc_marker)
{
  excitation = data.excitation;
  mode_idx = data.mode_idx;
  d_offset = data.d_offset;

  // Construct the SubMesh.
  MFEM_VERIFY(!data.attributes.empty(), "Wave port boundary found with no attributes!");
  mfem::ParMesh &mesh = *nd_fespace.GetParMesh();
  attr_list.Reserve(data.attributes.size());
  for (auto attr : data.attributes)
  {
    attr_list.Append(attr);
  }
  mesh::AttrToMarker(nd_fespace.GetParMesh()->bdr_attributes.Size()
                         ? nd_fespace.GetParMesh()->bdr_attributes.Max()
                         : 0,
                     attr_list, attr_marker);
  port_mesh = std::make_unique<mfem::ParSubMesh>(
      mfem::ParSubMesh::CreateFromBoundary(mesh, attr_list));

  int p_nd = nd_fespace.GetMaxElementOrder();
  int p_h1 = h1_fespace.GetMaxElementOrder();
  port_nd_fec = std::make_unique<mfem::ND_FECollection>(p_nd, mesh.Dimension() - 1);
  port_h1_fec = std::make_unique<mfem::H1_FECollection>(p_h1, mesh.Dimension() - 1);
  port_nd_fespace =
      std::make_unique<mfem::ParFiniteElementSpace>(port_mesh.get(), port_nd_fec.get());
  port_h1_fespace =
      std::make_unique<mfem::ParFiniteElementSpace>(port_mesh.get(), port_h1_fec.get());

  mfem::ParGridFunction E0t(const_cast<mfem::ParFiniteElementSpace *>(&nd_fespace)),
      E0n(const_cast<mfem::ParFiniteElementSpace *>(&h1_fespace));
  port_E0t = std::make_unique<mfem::ParComplexGridFunction>(port_nd_fespace.get());
  port_E0n = std::make_unique<mfem::ParComplexGridFunction>(port_h1_fespace.get());

  port_nd_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0t, port_E0t->real()));
  port_h1_transfer = std::make_unique<mfem::ParTransferMap>(
      mfem::ParSubMesh::CreateTransferMap(E0n, port_E0n->real()));

  // Extract Dirichlet BC true dofs for the port FE spaces.
  mfem::Array<int> port_nd_dbc_tdof_list, port_h1_dbc_tdof_list;
  GetEssentialTrueDofs(E0t, E0n, port_E0t->real(), port_E0n->real(), *port_nd_transfer,
                       *port_h1_transfer, dbc_marker, port_nd_dbc_tdof_list,
                       port_h1_dbc_tdof_list);

  // Construct operators for the generalized eigenvalue problem:
  //                [Aₜₜ  0] [eₜ]  = -kₙ² [Bₜₜ   Bₜₙ] [eₜ]
  //                [0   0] [eₙ]        [Bₜₙᵀ  Bₙₙ] [eₙ]
  // for the wave port of the given index. The transformed variables are related to the true
  // field by Eₜ = eₜ/kₙ and Eₙ = ieₙ. This is solved on the global mesh so the result is a
  // grid function over the entire space, not just the port boundary (so that it can be
  // queried from functions which use the global mesh).
  //
  // We will actually solve the shifted problem A e = λ B e, where:
  //                [Bₜₜ   Bₜₙ] [eₜ]   =  λ [Bₜₜ + 1/Θ² Aₜₜ  Bₜₙ] [eₜ]
  //                [Bₜₙᵀ  Bₙₙ] [eₙ]       [Bₜₙᵀ          Bₙₙ] [eₙ] .
  // Here we have λ = Θ²/(Θ²-kₙ²), where Θ² bounds the maximum kₙ² and is taken as Θ² =
  // ω² μₘₐₓ εₘₐₓ over the entire simulation domain.
  // Reference: Lee, Sun, and Cendes, Full-wave analysis of dielectric waveguides using
  //            tangential vector finite elements, IEEE Trans. Microwave Theory Tech.
  //            (1991).
  double c_min = mfem::infinity();
  for (auto attr : mesh.attributes)
  {
    c_min = std::min(c_min, mat_op.GetLightSpeedMin(attr));
  }
  MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
              "Invalid material speed of light detected in WavePortOperator!");
  mu_eps_max = 1.0 / (c_min * c_min);

  // Pre-compute problem matrices such that:
  //            A = A₁ - ω² A₂, B = A₁ - 1 / (μₘ εₘ) B₄ - ω² A₂ + 1/Θ² B₃ .
  {
    std::unique_ptr<mfem::HypreParMatrix> A1, B4r, B4i;
    {
      auto Btt = GetBtt(mat_op, *port_nd_fespace);
      auto Btn = GetBtn(mat_op, *port_nd_fespace, *port_h1_fespace);
      auto [Bnn1, Bnn2r, Bnn2i] = GetBnn(mat_op, *port_h1_fespace);
      auto [Att1, Att2r, Att2i] = GetAtt(mat_op, *port_nd_fespace);

      auto system_mats = GetSystemMatrices(
          std::move(Btt), std::move(Btn), std::move(Bnn1), std::move(Bnn2r),
          std::move(Bnn2i), std::move(Att1), std::move(Att2r), std::move(Att2i),
          port_nd_dbc_tdof_list, port_h1_dbc_tdof_list);
      A1 = std::move(system_mats[0]);
      A2r = std::move(system_mats[1]);
      A2i = std::move(system_mats[2]);
      B3 = std::move(system_mats[3]);
      B4r = std::move(system_mats[4]);
      B4i = std::move(system_mats[5]);
    }

    // Allocate storage for the eigenvalue problem operators. We have sparsity(A2) =
    // sparsity(B3) = sparsity(B4) ⊆ sparsity(A1). Precompute the frequency independent
    // contributions to A and B.
    P = std::make_unique<ComplexWrapperOperator>(
        std::make_unique<mfem::HypreParMatrix>(*A1), nullptr);
    if (A2i)
    {
      A = std::make_unique<ComplexWrapperOperator>(
          std::make_unique<mfem::HypreParMatrix>(*A1),
          std::make_unique<mfem::HypreParMatrix>(*A2i));
      B = std::make_unique<ComplexWrapperOperator>(
          std::make_unique<mfem::HypreParMatrix>(*A1),
          std::make_unique<mfem::HypreParMatrix>(*A2i));

      auto &Br = *static_cast<mfem::HypreParMatrix *>(B->Real());
      Br.Add(-1.0 / mu_eps_max, *B4r);

      auto &Ai = *static_cast<mfem::HypreParMatrix *>(A->Imag());
      auto &Bi = *static_cast<mfem::HypreParMatrix *>(B->Imag());
      Ai *= 0.0;
      Bi *= 0.0;
      Bi.Add(-1.0 / mu_eps_max, *B4i);
    }
    else
    {
      A = std::make_unique<ComplexWrapperOperator>(
          std::make_unique<mfem::HypreParMatrix>(*A1), nullptr);
      B = std::make_unique<ComplexWrapperOperator>(
          std::make_unique<mfem::HypreParMatrix>(*A1), nullptr);

      auto &Br = *static_cast<mfem::HypreParMatrix *>(B->Real());
      Br.Add(-1.0 / mu_eps_max, *B4r);
    }
  }

  // Create vector for initial space for eigenvalue solves (for nullspace of [Aₜₜ  0]
  //                                                                         [0   0] ).
  GetInitialSpace(*port_nd_fespace, *port_h1_fespace, port_nd_dbc_tdof_list,
                  port_h1_dbc_tdof_list, v0);
  e0.SetSize(v0.Size());
  e0t.SetSize(port_nd_fespace->GetTrueVSize());
  e0n.SetSize(port_h1_fespace->GetTrueVSize());

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
    constexpr int ksp_print = 0;
    constexpr double ksp_tol = 1.0e-8;
    constexpr double ksp_max_it = 30;
    auto gmres = std::make_unique<GmresSolver<ComplexOperator>>(port_comm, ksp_print);
    gmres->SetInitialGuess(false);
    gmres->SetRelTol(ksp_tol);
    gmres->SetMaxIter(ksp_max_it);
    gmres->SetRestartDim(ksp_max_it);
    // gmres->SetPrecSide(GmresSolver<ComplexOperator>::PrecSide::RIGHT);

    config::LinearSolverData::Type pc_type;
#if defined(MFEM_USE_SUPERLU)
    pc_type = config::LinearSolverData::Type::SUPERLU;
#elif defined(MFEM_USE_STRUMPACK)
    pc_type = config::LinearSolverData::Type::STRUMPACK;
#elif defined(MFEM_USE_MUMPS)
    pc_type = config::LinearSolverData::Type::MUMPS;
#else
#error "Wave port solver requires building with SuperLU_DIST, STRUMPACK, or MUMPS!"
#endif
    std::unique_ptr<Solver<ComplexOperator>> pc;
    if (pc_type == config::LinearSolverData::Type::SUPERLU)
    {
#if defined(MFEM_USE_SUPERLU)
      auto slu = std::make_unique<SuperLUSolver>(
          port_comm, config::LinearSolverData::SymFactType::DEFAULT, false, ksp_print - 1);
      // slu->GetSolver().SetColumnPermutation(mfem::superlu::NATURAL);
      pc = std::make_unique<WrapperSolver<ComplexOperator>>(std::move(slu));
#endif
    }
    else if (pc_type == config::LinearSolverData::Type::STRUMPACK)
    {
#if defined(MFEM_USE_STRUMPACK)
      auto strumpack = std::make_unique<StrumpackSolver>(
          port_comm, config::LinearSolverData::SymFactType::DEFAULT,
          config::LinearSolverData::CompressionType::NONE, 0.0, 0, 0, ksp_print - 1);
      // strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::NATURAL);
      pc = std::make_unique<WrapperSolver<ComplexOperator>>(std::move(strumpack));
#endif
    }
    else  // config::LinearSolverData::Type::MUMPS
    {
#if defined(MFEM_USE_MUMPS)
      auto mumps = std::make_unique<MumpsSolver>(
          port_comm, mfem::MUMPSSolver::SYMMETRIC_INDEFINITE,
          config::LinearSolverData::SymFactType::DEFAULT, 0.0, ksp_print - 1);
      // mumps->SetReorderingStrategy(mfem::MUMPSSolver::AMD);
      pc = std::make_unique<WrapperSolver<ComplexOperator>>(std::move(mumps));
#endif
    }
    ksp = std::make_unique<ComplexKspSolver>(std::move(gmres), std::move(pc));

    // Define the eigenvalue solver.
    constexpr int print = 0;
    config::EigenSolverData::Type type;
#if defined(PALACE_WITH_SLEPC)
    type = config::EigenSolverData::Type::SLEPC;
#elif defined(PALACE_WITH_ARPACK)
    type = config::EigenSolverData::Type::ARPACK;
#else
#error "Wave port solver requires building with ARPACK or SLEPc!"
#endif
    if (type == config::EigenSolverData::Type::ARPACK)
    {
#if defined(PALACE_WITH_ARPACK)
      eigen = std::make_unique<arpack::ArpackEPSSolver>(port_comm, print);
#endif
    }
    else  // config::EigenSolverData::Type::SLEPC
    {
#if defined(PALACE_WITH_SLEPC)
      auto slepc = std::make_unique<slepc::SlepcEPSSolver>(port_comm, print);
      slepc->SetType(slepc::SlepcEigenvalueSolver::Type::KRYLOVSCHUR);
      slepc->SetProblemType(slepc::SlepcEigenvalueSolver::ProblemType::GEN_NON_HERMITIAN);
      eigen = std::move(slepc);
#endif
    }
    constexpr double tol = 1.0e-6;
    eigen->SetNumModes(mode_idx, std::max(2 * mode_idx + 1, 5));
    eigen->SetTol(tol);
    eigen->SetWhichEigenpairs(EigenvalueSolver::WhichType::LARGEST_MAGNITUDE);
    eigen->SetLinearSolver(*ksp);
  }

  // Coefficients store references to kₙ, ω so they are updated implicitly at each new
  // solve. Also, μ⁻¹ is persistent, so no copy is OK.
  kn0 = 0.0;
  omega0 = 0.0;
  port_nxH0r_func =
      std::make_unique<BdrSubmeshHVectorCoefficient<true>>(*port_E0t, *port_E0n, mat_op);
  port_nxH0i_func =
      std::make_unique<BdrSubmeshHVectorCoefficient<false>>(*port_E0t, *port_E0n, mat_op);
  port_sr = std::make_unique<mfem::LinearForm>(port_nd_fespace.get());
  port_si = std::make_unique<mfem::LinearForm>(port_nd_fespace.get());
  port_sr->AddDomainIntegrator(new VectorFEDomainLFIntegrator(*port_nxH0r_func));
  port_si->AddDomainIntegrator(new VectorFEDomainLFIntegrator(*port_nxH0i_func));
  port_sr->UseFastAssembly(false);
  port_si->UseFastAssembly(false);

  // Configure port mode sign convention: 1ᵀ Re{-n x H} >= 0 on the "upper-right quadrant"
  // of the wave port boundary, in order to deal with symmetry effectively.
  {
    Vector bbmin, bbmax;
    port_mesh->GetBoundingBox(bbmin, bbmax);
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
    port_S0t = std::make_unique<mfem::ParGridFunction>(port_nd_fespace.get());
    port_S0t->ProjectCoefficient(tfunc);
  }
}

WavePortData::~WavePortData()
{
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

  // Use pre-computed matrices to construct and solve the generalized eigenvalue problem for
  // the desired wave port mode.
  double theta2 = mu_eps_max * omega * omega;
  {
    auto &Pr = *static_cast<mfem::HypreParMatrix *>(P->Real());
    Pr *= 0.0;

    auto &Ar = *static_cast<mfem::HypreParMatrix *>(A->Real());
    auto &Br = *static_cast<mfem::HypreParMatrix *>(B->Real());
    Ar.Add(-omega * omega + omega0 * omega0, *A2r);
    Br.Add(-omega * omega + omega0 * omega0, *A2r);
    Br.Add(1.0 / theta2 - ((omega0 == 0.0) ? 0.0 : 1.0 / (mu_eps_max * omega0 * omega0)),
           *B3);
    Pr.Add(1.0, Br);

    if (A2i)
    {
      auto &Ai = *static_cast<mfem::HypreParMatrix *>(A->Imag());
      auto &Bi = *static_cast<mfem::HypreParMatrix *>(B->Imag());
      Ai.Add(-omega * omega + omega0 * omega0, *A2i);
      Bi.Add(-omega * omega + omega0 * omega0, *A2i);
      Pr.Add(1.0, Bi);
    }
  }

  // Configure and solve the eigenvalue problem for the desired boundary mode.
  std::complex<double> lambda;
  if (port_comm != MPI_COMM_NULL)
  {
    ksp->SetOperators(*B, *P);
    eigen->SetOperators(*A, *B, EigenvalueSolver::ScaleType::NONE);
    eigen->SetInitialSpace(v0);
    int num_conv = eigen->Solve();
    MFEM_VERIFY(num_conv >= mode_idx, "Wave port eigensolver did not converge!");
    lambda = eigen->GetEigenvalue(mode_idx - 1);
    // Mpi::Print(port_comm, " ... Wave port eigensolver error = {} (bkwd), {} (abs)\n",
    //            eigen->GetError(mode_idx - 1, EigenvalueSolver::ErrorType::BACKWARD),
    //            eigen->GetError(mode_idx - 1, EigenvalueSolver::ErrorType::ABSOLUTE));
  }
  Mpi::Broadcast(1, &lambda, port_root, B3->GetComm());

  // Extract the eigenmode solution and postprocess. The extracted eigenvalue is λ =
  // Θ² / (Θ² - kₙ²).
  MFEM_VERIFY(lambda.real() > 1.0 / (1.0 - 1.0e-2),
              "Computed wave port mode is or is very close to being evanescent "
                  << "(λ = " << lambda << ")!");
  kn0 = std::sqrt(theta2 - theta2 / lambda);
  omega0 = omega;
  static_cast<BdrSubmeshHVectorCoefficient<true> *>(port_nxH0r_func.get())
      ->SetFrequency(omega0, kn0);
  static_cast<BdrSubmeshHVectorCoefficient<false> *>(port_nxH0i_func.get())
      ->SetFrequency(omega0, kn0);

  // Separate the computed field out into eₜ and eₙ and and transform back to true
  // electric field variables: Eₜ = eₜ/kₙ and Eₙ = ieₙ.
  if (port_comm != MPI_COMM_NULL)
  {
    Vector e0tr, e0ti, e0nr, e0ni;
    eigen->GetEigenvector(mode_idx - 1, e0);
    e0tr.MakeRef(e0.Real(), 0, e0t.Size());
    e0nr.MakeRef(e0.Real(), e0t.Size(), e0n.Size());
    e0ti.MakeRef(e0.Imag(), 0, e0t.Size());
    e0ni.MakeRef(e0.Imag(), e0t.Size(), e0n.Size());
    e0t.Real() = e0tr;
    e0t.Imag() = e0ti;
    e0n.Real() = e0nr;
    e0n.Imag() = e0ni;
    e0t *= 1.0 / kn0;
    e0n *= 1i;
  }
  else
  {
    MFEM_ASSERT(e0.Size() == 0 && e0t.Size() == 0 && e0n.Size() == 0,
                "Unexpected non-empty port FE space in wave port boundary mode solve!");
  }
  port_E0t->real().SetFromTrueDofs(e0t.Real());  // Parallel distribute
  port_E0t->imag().SetFromTrueDofs(e0t.Imag());
  port_E0n->real().SetFromTrueDofs(e0n.Real());
  port_E0n->imag().SetFromTrueDofs(e0n.Imag());

  // Normalize the mode for a chosen polarization direction and unit power, |E x H⋆| ⋅ n,
  // integrated over the port surface (+n is the direction of propagation).
  NormalizeWithSign(*port_S0t, *port_E0t, *port_E0n, *port_sr, *port_si);
}

double WavePortData::GetExcitationPower() const
{
  // The computed port modes are normalized such that the power integrated over the port is
  // 1: ∫ (E_inc x H_inc⋆) ⋅ n dS = 1.
  return excitation ? 1.0 : 0.0;
}

std::complex<double> WavePortData::GetSParameter(mfem::ParComplexGridFunction &E) const
{
  // Compute port S-parameter, or the projection of the field onto the port mode:
  // (E x H_inc⋆) ⋅ n = E ⋅ (-n x H_inc⋆), integrated over the port surface.
  mfem::ParComplexGridFunction port_E(port_nd_fespace.get());
  port_nd_transfer->Transfer(E.real(), port_E.real());
  port_nd_transfer->Transfer(E.imag(), port_E.imag());
  std::complex<double> dot(-((*port_sr) * port_E.real()) - ((*port_si) * port_E.imag()),
                           -((*port_sr) * port_E.imag()) + ((*port_si) * port_E.real()));
  Mpi::GlobalSum(1, &dot, port_nd_fespace->GetComm());
  return dot;
}

std::complex<double> WavePortData::GetPower(mfem::ParComplexGridFunction &E,
                                            mfem::ParComplexGridFunction &B,
                                            const MaterialOperator &mat_op) const
{
  // Compute port power, (E x H) ⋅ n = E ⋅ (-n x H), integrated over the port surface
  // using the computed E and H = μ⁻¹ B fields. The linear form is reconstructed from
  // scratch each time due to changing H. The BdrCurrentVectorCoefficient computes -n x H,
  // where n is an outward normal.
  auto &nd_fespace = *E.ParFESpace();
  BdrCurrentVectorCoefficient nxHr_func(B.real(), mat_op);
  BdrCurrentVectorCoefficient nxHi_func(B.imag(), mat_op);
  mfem::LinearForm pr(&nd_fespace), pi(&nd_fespace);
  pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHr_func), attr_marker);
  pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHi_func), attr_marker);
  pr.UseFastAssembly(false);
  pi.UseFastAssembly(false);
  pr.Assemble();
  pi.Assemble();
  std::complex<double> dot(-(pr * E.real()) - (pi * E.imag()),
                           -(pr * E.imag()) + (pi * E.real()));
  Mpi::GlobalSum(1, &dot, nd_fespace.GetComm());
  return dot;
}

WavePortOperator::WavePortOperator(const IoData &iod, const MaterialOperator &mat,
                                   const mfem::ParFiniteElementSpace &nd_fespace,
                                   const mfem::ParFiniteElementSpace &h1_fespace)
  : iodata(iod), mat_op(mat), suppress_output(false)
{
  // Set up wave port boundary conditions.
  MFEM_VERIFY(nd_fespace.GetParMesh() == h1_fespace.GetParMesh(),
              "Mesh mismatch in WavePortOperator FE spaces!");
  SetUpBoundaryProperties(iodata, mat_op, nd_fespace, h1_fespace);
  PrintBoundaryInfo(iodata, *nd_fespace.GetParMesh());
}

void WavePortOperator::SetUpBoundaryProperties(
    const IoData &iodata, const MaterialOperator &mat_op,
    const mfem::ParFiniteElementSpace &nd_fespace,
    const mfem::ParFiniteElementSpace &h1_fespace)
{
  // Check that wave port boundary attributes have been specified correctly.
  int bdr_attr_max = nd_fespace.GetParMesh()->bdr_attributes.Size()
                         ? nd_fespace.GetParMesh()->bdr_attributes.Max()
                         : 0;
  if (!iodata.boundaries.waveport.empty())
  {
    mfem::Array<int> bdr_attr_marker(bdr_attr_max);
    bdr_attr_marker = 0;
    for (auto attr : nd_fespace.GetParMesh()->bdr_attributes)
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
      }
    }
  }

  // List of all boundaries which will be marked as essential for the purposes of computing
  // wave port modes. This includes all PEC surfaces, but may also include others like when
  // a kinetic inductance or other BC is applied for the 3D simulation but should be
  // considered as PEC for the 2D problem.
  mfem::Array<int> dbc_bcs, dbc_marker;
  dbc_bcs.Reserve(static_cast<int>(iodata.boundaries.pec.attributes.size() +
                                   iodata.boundaries.auxpec.attributes.size()));
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
  // If user accidentally specifies a surface as both "PEC" and "WavePortPEC", this is fine
  // so allow for duplicates in the attribute list.
  dbc_bcs.Sort();
  dbc_bcs.Unique();
  mesh::AttrToMarker(bdr_attr_max, dbc_bcs, dbc_marker);

  // Set up wave port data structures.
  for (const auto &[idx, data] : iodata.boundaries.waveport)
  {
    ports.try_emplace(idx, data, mat_op, nd_fespace, h1_fespace, dbc_marker);
  }
  MFEM_VERIFY(
      ports.empty() || iodata.problem.type == config::ProblemData::Type::DRIVEN,
      "Wave port boundaries are only available for frequency domain driven simulations!");

  // Mark selected boundary attributes from the mesh for wave ports.
  port_marker.SetSize(bdr_attr_max);
  port_marker = 0;
  for (const auto &[idx, data] : ports)
  {
    for (int i = 0; i < data.GetMarker().Size(); i++)
    {
      MFEM_VERIFY(!(port_marker[i] && data.GetMarker()[i]),
                  "Boundary attribute is assigned to more than one wave port!");
      port_marker[i] = port_marker[i] || data.GetMarker()[i];
    }
  }
}

void WavePortOperator::PrintBoundaryInfo(const IoData &iodata, mfem::ParMesh &mesh)
{
  // Print out BC info for all port attributes.
  if (ports.empty())
  {
    return;
  }
  Mpi::Print("\nConfiguring Robin impedance BC for wave ports at attributes:\n");
  for (const auto &[idx, data] : ports)
  {
    for (int i = 0; i < data.GetMarker().Size(); i++)
    {
      if (!data.GetMarker()[i])
      {
        continue;
      }
      const int attr = i + 1;
      mfem::Vector nor;
      mesh::GetSurfaceNormal(mesh, attr, nor);
      Mpi::Print(
          " {:d}: Index = {:d}, mode = {:d}, d = {:.3e} m", attr, idx, data.GetModeIndex(),
          iodata.DimensionalizeValue(IoData::ValueType::LENGTH, data.GetOffsetDistance()));
      if (mesh.SpaceDimension() == 3)
      {
        Mpi::Print(", n = ({:+.1f}, {:+.1f}, {:+.1f})", nor(0), nor(1), nor(2));
      }
      else
      {
        Mpi::Print(", n = ({:+.1f}, {:+.1f})", nor(0), nor(1));
      }
      Mpi::Print("\n");
    }
  }

  // Print some information for excited wave ports.
  bool first = true;
  for (const auto &[idx, data] : ports)
  {
    if (!data.IsExcited())
    {
      continue;
    }
    if (first)
    {
      Mpi::Print("\nConfiguring wave port excitation source term at attributes:\n");
      first = false;
    }
    for (int i = 0; i < data.GetMarker().Size(); i++)
    {
      if (data.GetMarker()[i])
      {
        Mpi::Print(" {:d}: Index = {:d}\n", i + 1, idx);
      }
    }
  }
}

const WavePortData &WavePortOperator::GetPort(int idx) const
{
  auto it = ports.find(idx);
  MFEM_VERIFY(it != ports.end(), "Unknown wave port index requested!");
  return it->second;
}

void WavePortOperator::Initialize(double omega)
{
  bool init = false, first = true;
  for (const auto &[idx, data] : ports)
  {
    init = init || (data.GetOperatingFrequency() != omega);
    first = first && (data.GetOperatingFrequency() == 0.0);
  }
  if (!init)
  {
    return;
  }
  BlockTimer bt(Timer::WAVEPORT);
  if (!suppress_output)
  {
    const double freq = iodata.DimensionalizeValue(IoData::ValueType::FREQUENCY, omega);
    Mpi::Print(
        "\nCalculating boundary modes at wave ports for ω/2π = {:.3e} GHz ({:.3e})\n", freq,
        omega);
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
      double k0 = 1.0 / iodata.DimensionalizeValue(IoData::ValueType::LENGTH, 1.0);
      Mpi::Print(" Port {:d}, mode {:d}: kₙ = {:.3e}{:+.3e}i m⁻¹\n", idx,
                 data.GetModeIndex(), k0 * data.GetPropagationConstant().real(),
                 k0 * data.GetPropagationConstant().imag());
    }
  }
}

void WavePortOperator::AddExtraSystemBdrCoefficients(double omega,
                                                     SumMatrixCoefficient &fbr,
                                                     SumMatrixCoefficient &fbi)
{
  // Add wave port boundaries to the bilinear form. This looks a lot like the lumped port
  // boundary, except the iω / Z_s coefficient goes to ikₙ / μ where kₙ is specific to the
  // port mode at the given operating frequency (note only the real part of the propagation
  // constant contributes).
  Initialize(omega);
  for (auto &[idx, data] : ports)
  {
    constexpr auto MatType = MaterialPropertyType::INV_PERMEABILITY;
    constexpr auto ElemType = MeshElementType::BDR_ELEMENT;
    fbi.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType, ElemType>>(
                           mat_op, data.GetPropagationConstant().real()),
                       data.GetMarker());
  }
}

void WavePortOperator::AddExcitationBdrCoefficients(double omega, SumVectorCoefficient &fbr,
                                                    SumVectorCoefficient &fbi)
{
  // Re{-U_inc} = Re{+2 (-iω) n x H_inc}, which is a function of E_inc as computed by the
  // modal solution (stored as a grid function and coefficient during initialization).
  // Likewise for the imaginary part.
  Initialize(omega);
  for (auto &[idx, data] : ports)
  {
    if (!data.IsExcited())
    {
      continue;
    }
    fbr.AddCoefficient(std::make_unique<mfem::ScalarVectorProductCoefficient>(
                           2.0 * omega, data.GetModeCoefficientImag()),
                       data.GetMarker());
    fbi.AddCoefficient(std::make_unique<mfem::ScalarVectorProductCoefficient>(
                           -2.0 * omega, data.GetModeCoefficientReal()),
                       data.GetMarker());
  }
}

}  // namespace palace
