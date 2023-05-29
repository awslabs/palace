// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "waveportoperator.hpp"

#include <array>
#include <tuple>
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

namespace palace
{

using namespace std::complex_literals;

namespace
{

void GetEssentialTrueDofs(mfem::ParFiniteElementSpace &nd_fespace,
                          mfem::ParFiniteElementSpace &h1_fespace,
                          const mfem::Array<int> &attr_marker,
                          const mfem::Array<int> &dbc_marker,
                          mfem::Array<int> &nd_dbc_tdof_list,
                          mfem::Array<int> &h1_dbc_tdof_list,
                          HYPRE_BigInt attr_tdof_sizes[2])
{
  // Mark all ND and H1 dofs which are not on the port, and then mark PEC boundaries on
  // the port as well.
  mfem::Array<int> nd_attr_tdof_list, h1_attr_tdof_list;
  nd_fespace.GetEssentialTrueDofs(attr_marker, nd_attr_tdof_list);
  h1_fespace.GetEssentialTrueDofs(attr_marker, h1_attr_tdof_list);
  nd_fespace.GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  attr_tdof_sizes[0] = nd_attr_tdof_list.Size();
  attr_tdof_sizes[1] = h1_attr_tdof_list.Size();
  Mpi::GlobalSum(2, attr_tdof_sizes, nd_fespace.GetComm());

  mfem::Array<int> nd_dbc_tdof_marker(nd_fespace.GetTrueVSize()),
      h1_dbc_tdof_marker(h1_fespace.GetTrueVSize());
  nd_dbc_tdof_marker = 1;
  h1_dbc_tdof_marker = 1;
  for (auto tdof : nd_attr_tdof_list)
  {
    nd_dbc_tdof_marker[tdof] = 0;
  }
  for (auto tdof : nd_dbc_tdof_list)
  {
    nd_dbc_tdof_marker[tdof] = 1;
  }
  for (auto tdof : h1_attr_tdof_list)
  {
    h1_dbc_tdof_marker[tdof] = 0;
  }
  for (auto tdof : h1_dbc_tdof_list)
  {
    h1_dbc_tdof_marker[tdof] = 1;
  }

  // Convert back to a list.
  nd_dbc_tdof_list.DeleteAll();
  nd_dbc_tdof_list.Reserve(nd_fespace.GetTrueVSize());
  for (int i = 0; i < nd_dbc_tdof_marker.Size(); i++)
  {
    if (nd_dbc_tdof_marker[i])
    {
      nd_dbc_tdof_list.Append(i);
    }
  }
  h1_dbc_tdof_list.DeleteAll();
  h1_dbc_tdof_list.Reserve(h1_fespace.GetTrueVSize());
  for (int i = 0; i < h1_dbc_tdof_marker.Size(); i++)
  {
    if (h1_dbc_tdof_marker[i])
    {
      h1_dbc_tdof_list.Append(i);
    }
  }
}

constexpr int skip_zeros = 0;

std::unique_ptr<ParOperator> GetBtt(const MaterialOperator &mat_op,
                                    mfem::ParFiniteElementSpace &nd_fespace,
                                    mfem::Array<int> &attr_marker)
{
  // Mass matrix: Bₜₜ = (μ⁻¹ u, v).
  constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_PERMEABILITY;
  BdrMaterialPropertyCoefficient<MatType> muinv_func(mat_op);
  auto btt = std::make_unique<mfem::SymmetricBilinearForm>(&nd_fespace);
  btt->AddBoundaryIntegrator(new mfem::MixedVectorMassIntegrator(muinv_func), attr_marker);
  btt->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
  btt->Assemble(skip_zeros);
  btt->Finalize(skip_zeros);
  return std::make_unique<ParOperator>(std::move(btt), nd_fespace);
}

std::unique_ptr<ParOperator> GetBtn(const MaterialOperator &mat_op,
                                    mfem::ParFiniteElementSpace &nd_fespace,
                                    mfem::ParFiniteElementSpace &h1_fespace,
                                    mfem::Array<int> &attr_marker)
{
  // Mass matrix: Bₜₙ = (μ⁻¹ ∇ₜ u, v).
  constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_PERMEABILITY;
  BdrMaterialPropertyCoefficient<MatType> muinv_func(mat_op);
  auto btn = std::make_unique<mfem::MixedBilinearForm>(&h1_fespace, &nd_fespace);
  btn->AddBoundaryIntegrator(new mfem::MixedVectorGradientIntegrator(muinv_func),
                             attr_marker);
  btn->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
  btn->Assemble(skip_zeros);
  btn->Finalize(skip_zeros);
  return std::make_unique<ParOperator>(std::move(btn), h1_fespace, nd_fespace, false);
}

std::array<std::unique_ptr<ParOperator>, 3> GetBnn(const MaterialOperator &mat_op,
                                                   mfem::ParFiniteElementSpace &h1_fespace,
                                                   mfem::Array<int> &attr_marker)
{
  // Mass matrix: Bₙₙ = (μ⁻¹ ∇ₜ u, ∇ₜ v) - ω² (ε u, v) = Bₙₙ₁ - ω² Bₙₙ₂.
  constexpr MaterialPropertyType MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  BdrMaterialPropertyCoefficient<MatTypeMuInv> muinv_func(mat_op);
  auto bnn1 = std::make_unique<mfem::SymmetricBilinearForm>(&h1_fespace);
  bnn1->AddBoundaryIntegrator(new mfem::MixedGradGradIntegrator(muinv_func), attr_marker);
  bnn1->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
  bnn1->Assemble(skip_zeros);
  bnn1->Finalize(skip_zeros);

  constexpr MaterialPropertyType MatTypeEpsReal = MaterialPropertyType::PERMITTIVITY_REAL;
  NormalProjectedCoefficient epsilon_func(
      std::make_unique<BdrMaterialPropertyCoefficient<MatTypeEpsReal>>(mat_op));
  auto bnn2r = std::make_unique<mfem::SymmetricBilinearForm>(&h1_fespace);
  bnn2r->AddBoundaryIntegrator(new mfem::MixedScalarMassIntegrator(epsilon_func),
                               attr_marker);
  bnn2r->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
  bnn2r->Assemble(skip_zeros);
  bnn2r->Finalize(skip_zeros);

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {std::make_unique<ParOperator>(std::move(bnn1), h1_fespace),
            std::make_unique<ParOperator>(std::move(bnn2r), h1_fespace), nullptr};
  }
  constexpr MaterialPropertyType MatTypeEpsImag = MaterialPropertyType::PERMITTIVITY_IMAG;
  NormalProjectedCoefficient negepstandelta_func(
      std::make_unique<BdrMaterialPropertyCoefficient<MatTypeEpsImag>>(mat_op));
  auto bnn2i = std::make_unique<mfem::SymmetricBilinearForm>(&h1_fespace);
  bnn2i->AddBoundaryIntegrator(new mfem::MixedScalarMassIntegrator(negepstandelta_func),
                               attr_marker);
  bnn2i->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
  bnn2i->Assemble(skip_zeros);
  bnn2i->Finalize(skip_zeros);
  return {std::make_unique<ParOperator>(std::move(bnn1), h1_fespace),
          std::make_unique<ParOperator>(std::move(bnn2r), h1_fespace),
          std::make_unique<ParOperator>(std::move(bnn2i), h1_fespace)};
}

std::array<std::unique_ptr<ParOperator>, 3> GetAtt(const MaterialOperator &mat_op,
                                                   mfem::ParFiniteElementSpace &nd_fespace,
                                                   mfem::Array<int> &attr_marker)
{
  // Stiffness matrix: Aₜₜ = (μ⁻¹ ∇ₜ x u, ∇ₜ x v) - ω² (ε u, v) = Aₜₜ₁ - ω² Aₜₜ₂.
  constexpr MaterialPropertyType MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  NormalProjectedCoefficient muinv_func(
      std::make_unique<BdrMaterialPropertyCoefficient<MatTypeMuInv>>(mat_op));
  auto att1 = std::make_unique<mfem::SymmetricBilinearForm>(&nd_fespace);
  att1->AddBoundaryIntegrator(new mfem::CurlCurlIntegrator(muinv_func), attr_marker);
  att1->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
  att1->Assemble(skip_zeros);
  att1->Finalize(skip_zeros);

  constexpr MaterialPropertyType MatTypeEpsReal = MaterialPropertyType::PERMITTIVITY_REAL;
  BdrMaterialPropertyCoefficient<MatTypeEpsReal> epsilon_func(mat_op);
  auto att2r = std::make_unique<mfem::SymmetricBilinearForm>(&nd_fespace);
  att2r->AddBoundaryIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func),
                               attr_marker);
  att2r->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
  att2r->Assemble(skip_zeros);
  att2r->Finalize(skip_zeros);

  // Contribution for loss tangent: ε -> ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {std::make_unique<ParOperator>(std::move(att1), nd_fespace),
            std::make_unique<ParOperator>(std::move(att2r), nd_fespace), nullptr};
  }
  constexpr MaterialPropertyType MatTypeEpsImag = MaterialPropertyType::PERMITTIVITY_IMAG;
  BdrMaterialPropertyCoefficient<MatTypeEpsImag> negepstandelta_func(mat_op);
  auto att2i = std::make_unique<mfem::SymmetricBilinearForm>(&nd_fespace);
  att2i->AddBoundaryIntegrator(new mfem::MixedVectorMassIntegrator(negepstandelta_func),
                               attr_marker);
  att2i->SetAssemblyLevel(mfem::AssemblyLevel::LEGACY);
  att2i->Assemble(skip_zeros);
  att2i->Finalize(skip_zeros);
  return {std::make_unique<ParOperator>(std::move(att1), nd_fespace),
          std::make_unique<ParOperator>(std::move(att2r), nd_fespace),
          std::make_unique<ParOperator>(std::move(att2i), nd_fespace)};
}

std::array<std::unique_ptr<mfem::HypreParMatrix>, 6>
GetSystemMatrices(std::unique_ptr<ParOperator> Btt, std::unique_ptr<ParOperator> Btn,
                  std::unique_ptr<ParOperator> Bnn1, std::unique_ptr<ParOperator> Bnn2r,
                  std::unique_ptr<ParOperator> Bnn2i, std::unique_ptr<ParOperator> Att1,
                  std::unique_ptr<ParOperator> Att2r, std::unique_ptr<ParOperator> Att2i,
                  mfem::Array<int> &nd_dbc_tdof_list, mfem::Array<int> &h1_dbc_tdof_list)
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

  // Eliminate boundary tdofs not associated with this wave port or constrained by Dirichlet
  // BCs. It is not guaranteed that any HypreParMatrix has a full diagonal in its sparsity
  // pattern, so we add a zero diagonal before elimination to guarantee this for A1 and B3.
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

void GetInitialSpace(mfem::ParFiniteElementSpace &nd_fespace,
                     mfem::ParFiniteElementSpace &h1_fespace,
                     const mfem::Array<int> &nd_dbc_tdof_list,
                     const mfem::Array<int> &h1_dbc_tdof_list, ComplexVector &v)
{
  // Initial space chosen as such that B v₀ = y₀, with y₀ = [y₀ₜ, 0, ... 0]ᵀ ⟂ null(A)
  // (with Aₜₜ nonsingular). See Lee, Sun, and Cendes, 1991 for reference.
  // Note: When the eigenvalue solver uses a standard ℓ²-inner product instead of B-inner
  // product (since we use a general non-Hermitian solver due to complex symmetric B), then
  // we just use v0 = y0 directly.
  v.SetSize(nd_fespace.GetTrueVSize() + h1_fespace.GetTrueVSize());
  linalg::SetRandom(nd_fespace.GetComm(), v);
  // v = std::complex<double>(1.0, 0.0);
  v.SetSubVector(nd_dbc_tdof_list, 0.0);
  for (int i = nd_fespace.GetTrueVSize();
       i < nd_fespace.GetTrueVSize() + h1_fespace.GetTrueVSize(); i++)
  {
    v.Real()[i] = v.Imag()[i] = 0.0;
  }
}

}  // namespace

// Computes boundary modal n x H, where +n is the direction of wave propagation: n x H =
// -1/(iωμ) (ikₙ Eₜ + ∇ₜ Eₙ), using the tangential and normal electric field component grid
// functions evaluated on the (single-sided) boundary element. The intent of this vector
// grid function is to be dotted with a function E which is only in the tangential
// component, so the fact that we use the full ∇ Eₙ in the element is fine. We use only the
// real part of kn.
class BdrHVectorCoefficient : public mfem::VectorCoefficient
{
private:
  const mfem::ParComplexGridFunction &gridfunc_t, &gridfunc_n;
  const MaterialOperator &mat_op;
  const bool imaginary;
  std::complex<double> kn;
  double omega;

public:
  BdrHVectorCoefficient(const mfem::ParComplexGridFunction &Et,
                        const mfem::ParComplexGridFunction &En, const MaterialOperator &op,
                        bool imag)
    : mfem::VectorCoefficient(Et.ParFESpace()->GetParMesh()->SpaceDimension()),
      gridfunc_t(Et), gridfunc_n(En), mat_op(op), imaginary(imag), kn(0.0), omega(0.0)
  {
  }

  void Eval(mfem::Vector &V, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override
  {
    MFEM_VERIFY(T.ElementType == mfem::ElementTransformation::BDR_ELEMENT,
                "Unexpected element type in BdrHVectorCoefficient!");
    MFEM_VERIFY(gridfunc_t.ParFESpace()->GetParMesh() == T.mesh &&
                    gridfunc_n.ParFESpace()->GetParMesh() == T.mesh,
                "Invalid mesh for BdrHVectorCoefficient!");

    // This coefficient is only to be used on true exterior boundaries.
    int i, o;
    int iel1, iel2, info1, info2;
    const mfem::Mesh &mesh = *T.mesh;
    mesh.GetBdrElementFace(T.ElementNo, &i, &o);
    mesh.GetFaceElements(i, &iel1, &iel2);
    mesh.GetFaceInfos(i, &info1, &info2);
    if (info2 >= 0)
    {
      // Just return for an non-true boundary face.
      V.SetSize(vdim);
      V = 0.0;
      return;
    }

    // Compute Re/Im{-1/i (ikₙ Eₜ + ∇ₜ Eₙ)}.
    T.SetIntPoint(&ip);
    mfem::Vector U;
    if (!imaginary)
    {
      gridfunc_t.real().GetVectorValue(T, ip, U);
      U *= -kn.real();

      mfem::Vector dU;
      gridfunc_n.imag().GetGradient(T, dU);
      U -= dU;
    }
    else
    {
      gridfunc_t.imag().GetVectorValue(T, ip, U);
      U *= -kn.real();

      mfem::Vector dU;
      gridfunc_n.real().GetGradient(T, dU);
      U += dU;
    }

    // Scale by 1/(ωμ) with μ evaluated in the neighboring element.
    V.SetSize(U.Size());
    mat_op.GetInvPermeability(mesh.GetAttribute(iel1)).Mult(U, V);
    V *= (1.0 / omega);
  }

  void SetFrequency(double w, std::complex<double> k)
  {
    omega = w;
    kn = k;
  }
};

WavePortData::WavePortData(const config::WavePortData &data, const MaterialOperator &mat_op,
                           mfem::ParFiniteElementSpace &nd_fespace,
                           mfem::ParFiniteElementSpace &h1_fespace,
                           const mfem::Array<int> &dbc_marker)
{
  excitation = data.excitation;
  mode_idx = data.mode_idx;
  d_offset = data.d_offset;
  MFEM_VERIFY(!data.attributes.empty(), "Wave port boundary found with no attributes!");
  mesh::AttrToMarker(nd_fespace.GetParMesh()->bdr_attributes.Max(), data.attributes,
                     attr_marker);

  // Construct operators for the generalized eigenvalue problem:
  //                [Aₜₜ  0] [eₜ]  = -kₙ² [Bₜₜ   Bₜₙ]  [eₜ]
  //                [0   0] [eₙ]        [Bₜₙᵀ  Bₙₙ] [eₙ]
  // for the wave port of the given index. The transformed variables are related to the true
  // field by Eₜ = eₜ/kₙ and Eₙ = ieₙ. This is solved on the global mesh so the result is a
  // grid function over the entire space, not just the port boundary (so that it can be
  // queried from functions which use the global mesh).
  //
  // We will actually solve the shifted problem A e = λ B e, where (see Lee, Sun, and
  // Cendes, 1991):
  //                [Bₜₜ   Bₜₙ]  [eₜ]  =  λ [Bₜₜ + 1/Θ² Aₜₜ  Bₜₙ] [eₜ]
  //                [Bₜₙᵀ  Bₙₙ] [eₙ]       [Bₜₙᵀ          Bₙₙ] [eₙ] .
  // Here we have λ = Θ²/(Θ²-kₙ²), where Θ² bounds the maximum kₙ² and is taken as Θ² =
  // ω² μₘₐₓ εₘₐₓ over the entire simulation domain.
  double c_min = mfem::infinity();
  for (auto attr : nd_fespace.GetParMesh()->attributes)
  {
    c_min = std::min(c_min, mat_op.GetLightSpeedMin(attr));
  }
  MFEM_VERIFY(c_min > 0.0 && c_min < mfem::infinity(),
              "Invalid material speed of light detected in WavePortOperator!");
  mu_eps_max = 1.0 / (c_min * c_min);

  // Pre-compute problem matrices such that:
  //            A = A₁ - ω² A₂, B = A₁ - 1 / (μₘ εₘ) B₄ - ω² A₂ + 1/Θ² B₃ .
  mfem::Array<int> nd_dbc_tdof_list, h1_dbc_tdof_list;
  GetEssentialTrueDofs(nd_fespace, h1_fespace, attr_marker, dbc_marker, nd_dbc_tdof_list,
                       h1_dbc_tdof_list, attr_tdof_sizes);
  {
    std::unique_ptr<mfem::HypreParMatrix> A1, B4r, B4i;
    {
      auto Btt = GetBtt(mat_op, nd_fespace, attr_marker);
      auto Btn = GetBtn(mat_op, nd_fespace, h1_fespace, attr_marker);
      auto [Bnn1, Bnn2r, Bnn2i] = GetBnn(mat_op, h1_fespace, attr_marker);
      auto [Att1, Att2r, Att2i] = GetAtt(mat_op, nd_fespace, attr_marker);

      auto system_mats = GetSystemMatrices(
          std::move(Btt), std::move(Btn), std::move(Bnn1), std::move(Bnn2r),
          std::move(Bnn2i), std::move(Att1), std::move(Att2r), std::move(Att2i),
          nd_dbc_tdof_list, h1_dbc_tdof_list);
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
  GetInitialSpace(nd_fespace, h1_fespace, nd_dbc_tdof_list, h1_dbc_tdof_list, v0);
  e0.SetSize(v0.Size());
  e0t.SetSize(nd_fespace.GetTrueVSize());
  e0n.SetSize(h1_fespace.GetTrueVSize());

  // Configure the eigenvalue problem solver. As for the full 3D case, the system matrices
  // are in general complex and symmetric. We supply the operators to the solver in
  // shift-inverted form and handle the back-transformation externally.
  {
    // Define the linear solver to be used for solving systems associated with the
    // generalized eigenvalue problem.
    constexpr int ksp_print = 0;
    constexpr double ksp_tol = 1.0e-8;
    constexpr double ksp_max_it = 100;
    auto gmres =
        std::make_unique<GmresSolver<ComplexOperator>>(nd_fespace.GetComm(), ksp_print);
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
      auto slu =
          std::make_unique<SuperLUSolver>(nd_fespace.GetComm(), 0, false, ksp_print - 1);
      slu->GetSolver().SetColumnPermutation(mfem::superlu::NATURAL);
      pc = std::make_unique<WrapperSolver<ComplexOperator>>(std::move(slu));
#endif
    }
    else if (pc_type == config::LinearSolverData::Type::STRUMPACK)
    {
#if defined(MFEM_USE_STRUMPACK)
      auto strumpack = std::make_unique<StrumpackSolver>(nd_fespace.GetComm(), 0,
                                                         strumpack::CompressionType::NONE,
                                                         0.0, 0, 0, ksp_print - 1);
      strumpack->SetReorderingStrategy(strumpack::ReorderingStrategy::NATURAL);
      pc = std::make_unique<WrapperSolver<ComplexOperator>>(std::move(strumpack));
#endif
    }
    else  // config::LinearSolverData::Type::MUMPS
    {
#if defined(MFEM_USE_MUMPS)
      auto mumps = std::make_unique<MumpsSolver>(nd_fespace.GetComm(),
                                                 mfem::MUMPSSolver::SYMMETRIC_INDEFINITE, 0,
                                                 0.0, ksp_print - 1);
      mumps->SetReorderingStrategy(mfem::MUMPSSolver::AMD);
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
      eigen = std::make_unique<arpack::ArpackEPSSolver>(nd_fespace.GetComm(), print);
#endif
    }
    else  // config::EigenSolverData::Type::SLEPC
    {
#if defined(PALACE_WITH_SLEPC)
      auto slepc = std::make_unique<slepc::SlepcEPSSolver>(nd_fespace.GetComm(), print);
      slepc->SetType(slepc::SlepcEigenSolver::Type::KRYLOVSCHUR);
      slepc->SetProblemType(slepc::SlepcEigenSolver::ProblemType::GEN_NON_HERMITIAN);
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
  E0t = std::make_unique<mfem::ParComplexGridFunction>(&nd_fespace);
  E0n = std::make_unique<mfem::ParComplexGridFunction>(&h1_fespace);
  nxH0r_func = std::make_unique<BdrHVectorCoefficient>(*E0t, *E0n, mat_op, false);
  nxH0i_func = std::make_unique<BdrHVectorCoefficient>(*E0t, *E0n, mat_op, true);
  sr = std::make_unique<mfem::ParLinearForm>(&nd_fespace);
  si = std::make_unique<mfem::ParLinearForm>(&nd_fespace);
  sr->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(*nxH0r_func), attr_marker);
  si->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(*nxH0i_func), attr_marker);
  sr->UseFastAssembly(false);
  si->UseFastAssembly(false);
  ones = std::make_unique<mfem::ParGridFunction>(&nd_fespace);
  *ones = 1.0;
}  // namespace palace

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
    Br.Add(1.0 / theta2 - (omega0 == 0.0 ? 0.0 : 1.0 / (mu_eps_max * omega0 * omega0)),
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
  ksp->SetOperators(*B, *P);
  eigen->SetOperators(*A, *B, EigenvalueSolver::ScaleType::NONE);
  eigen->SetInitialSpace(v0);
  int num_conv = eigen->Solve();
  MFEM_VERIFY(num_conv >= mode_idx, "Wave port eigensolver did not converge!");
  std::complex<double> lambda = eigen->GetEigenvalue(mode_idx - 1);
  // Mpi::Print(" ... Wave port eigensolver error = {} (bkwd), {} (abs)\n",
  //            eigen->GetError(mode_idx - 1, EigenvalueSolver::ErrorType::BACKWARD),
  //            eigen->GetError(mode_idx - 1, EigenvalueSolver::ErrorType::ABSOLUTE));

  // Extract the eigenmode solution and postprocess. The extracted eigenvalue is λ =
  // Θ² / (Θ² - kₙ²).
  MFEM_VERIFY(lambda.real() > 1.0 / (1.0 - 1.0e-2),
              "Computed wave port mode is or is very close to being evanescent "
                  << "(λ = " << lambda << ")!");
  kn0 = std::sqrt(theta2 - theta2 / lambda);
  omega0 = omega;
  static_cast<BdrHVectorCoefficient *>(nxH0r_func.get())->SetFrequency(omega0, kn0);
  static_cast<BdrHVectorCoefficient *>(nxH0i_func.get())->SetFrequency(omega0, kn0);

  // Separate the computed field out into eₜ and eₙ and and transform back to true electric
  // field variables: Eₜ = eₜ/kₙ and Eₙ = ieₙ.
  eigen->GetEigenvector(mode_idx - 1, e0);
  {
    Vector e0tr, e0ti, e0nr, e0ni;
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
    E0t->real().SetFromTrueDofs(e0t.Real());  // Parallel distribute
    E0t->imag().SetFromTrueDofs(e0t.Imag());
    E0n->real().SetFromTrueDofs(e0n.Real());
    E0n->imag().SetFromTrueDofs(e0n.Imag());
  }

  // Normalize grid functions to a chosen polarization direction and unit power, |E x H⋆| ⋅
  // n, integrated over the port surface (+n is the direction of propagation). The n x H
  // coefficients are updated implicitly as the only store references to the Et, En grid
  // functions as well as kₙ, ω. We choose a (rather arbitrary) sign constraint to at least
  // make results for the same port consistent between frequencies/meshes.
  {
    // |E x H⋆| ⋅ n = |E ⋅ (-n x H⋆)|
    *sr = 0.0;
    *si = 0.0;
    sr->Assemble();
    si->Assemble();
    double sign = ((*sr)(*ones) > 0.0) ? 1.0 : -1.0;
    std::complex<double> s0(-(*sr)(E0t->real()) - (*si)(E0t->imag()),
                            -(*sr)(E0t->imag()) + (*si)(E0t->real()));
    double scale = sign / std::sqrt(std::abs(s0));
    E0t->real() *= scale;  // This updates the n x H coefficients depending on Et, En too
    E0t->imag() *= scale;
    E0n->real() *= scale;
    E0n->imag() *= scale;
    *sr *= scale;  // Update linear forms for postprocessing
    *si *= scale;
  }

  // This parallel communication is not required since wave port boundaries are true
  // one-sided boundaries.
  // E0t->real().ExchangeFaceNbrData();  // Ready for parallel comm on shared faces
  // E0t->imag().ExchangeFaceNbrData();  // for n x H coefficients evaluation
  // E0n->real().ExchangeFaceNbrData();
  // E0n->imag().ExchangeFaceNbrData();
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
  return {-(*sr)(E.real()) - (*si)(E.imag()), -(*sr)(E.imag()) + (*si)(E.real())};
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
  mfem::ParLinearForm pr(&nd_fespace), pi(&nd_fespace);
  pr.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHr_func), attr_marker);
  pi.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(nxHi_func), attr_marker);
  pr.UseFastAssembly(false);
  pi.UseFastAssembly(false);
  pr.Assemble();
  pi.Assemble();
  return {pr(E.real()) + pi(E.imag()), pr(E.imag()) - pi(E.real())};
}

WavePortOperator::WavePortOperator(const IoData &iod, const MaterialOperator &mat,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::ParFiniteElementSpace &h1_fespace)
  : iodata(iod), mat_op(mat), suppress_output(false)
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
  int bdr_attr_max = nd_fespace.GetParMesh()->bdr_attributes.Max();
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
    constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_PERMEABILITY;
    fbi.AddCoefficient(std::make_unique<BdrMaterialPropertyCoefficient<MatType>>(
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
                           2.0 * omega, *data.GetModeCoefficientImag()),
                       data.GetMarker());
    fbi.AddCoefficient(std::make_unique<mfem::ScalarVectorProductCoefficient>(
                           -2.0 * omega, *data.GetModeCoefficientReal()),
                       data.GetMarker());
  }
}

}  // namespace palace
