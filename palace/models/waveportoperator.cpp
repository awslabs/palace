// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "waveportoperator.hpp"

#include <optional>
#include "fem/coefficient.hpp"
#include "fem/integrator.hpp"
#include "linalg/arpack.hpp"
#include "linalg/hypre.hpp"
#include "linalg/slepc.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

constexpr int skip_zeros = 0;

inline mfem::HypreParMatrix GetBtt(const MaterialOperator &mat_op,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::Array<int> &attr_marker)
{
  // Mass matrix: Bₜₜ = (μ⁻¹ u, v).
  constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_PERMEABILITY;
  MaterialPropertyCoefficient<MatType> muinv_func(mat_op);
  mfem::ParBilinearForm btt(&nd_fespace);
  btt.AddBoundaryIntegrator(new mfem::MixedVectorMassIntegrator(muinv_func), attr_marker);
  // btt.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  btt.Assemble(skip_zeros);
  btt.Finalize(skip_zeros);
  return *btt.ParallelAssemble();
}

inline mfem::HypreParMatrix GetBtn(const MaterialOperator &mat_op,
                                   mfem::ParFiniteElementSpace &nd_fespace,
                                   mfem::ParFiniteElementSpace &h1_fespace,
                                   mfem::Array<int> &attr_marker)
{
  // Mass matrix: Bₜₙ = (μ⁻¹ ∇ₜ u, v).
  constexpr MaterialPropertyType MatType = MaterialPropertyType::INV_PERMEABILITY;
  MaterialPropertyCoefficient<MatType> muinv_func(mat_op);
  mfem::ParMixedBilinearForm btn(&h1_fespace, &nd_fespace);
  btn.AddBoundaryIntegrator(new mfem::MixedVectorGradientIntegrator(muinv_func),
                            attr_marker);
  // btn.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  btn.Assemble(skip_zeros);
  btn.Finalize(skip_zeros);
  return *btn.ParallelAssemble();
}

struct Bnn
{
  mfem::HypreParMatrix Bnn1;
  mfem::HypreParMatrix Bnn2r;
  std::optional<mfem::HypreParMatrix> Bnn2i;
};

inline Bnn GetBnn(const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &h1_fespace,
                  mfem::Array<int> &attr_marker)
{
  // Mass matrix: Bₙₙ = (μ⁻¹ ∇ₜ u, ∇ₜ v) - ω² (ε u, v) = Bₙₙ₁ - ω² Bₙₙ₂.
  constexpr MaterialPropertyType MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  MaterialPropertyCoefficient<MatTypeMuInv> muinv_func(mat_op);
  mfem::ParBilinearForm bnn1(&h1_fespace);
  bnn1.AddBoundaryIntegrator(new mfem::MixedGradGradIntegrator(muinv_func), attr_marker);
  // bnn1.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  bnn1.Assemble(skip_zeros);
  bnn1.Finalize(skip_zeros);

  constexpr MaterialPropertyType MatTypeEpsReal = MaterialPropertyType::PERMITTIVITY_REAL;
  NormalProjectedCoefficient epsilon_func(
      std::make_unique<MaterialPropertyCoefficient<MatTypeEpsReal>>(mat_op));
  mfem::ParBilinearForm bnn2r(&h1_fespace);
  bnn2r.AddBoundaryIntegrator(new mfem::MixedScalarMassIntegrator(epsilon_func),
                              attr_marker);
  // bnn2r.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  bnn2r.Assemble(skip_zeros);
  bnn2r.Finalize(skip_zeros);

  // Contribution for loss tangent: ε => ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {*bnn1.ParallelAssemble(), *bnn2r.ParallelAssemble()};
  }
  constexpr MaterialPropertyType MatTypeEpsImag = MaterialPropertyType::PERMITTIVITY_IMAG;
  NormalProjectedCoefficient negepstandelta_func(
      std::make_unique<MaterialPropertyCoefficient<MatTypeEpsImag>>(mat_op));
  mfem::ParBilinearForm bnn2i(&h1_fespace);
  bnn2i.AddBoundaryIntegrator(new mfem::MixedScalarMassIntegrator(negepstandelta_func),
                              attr_marker);
  // bnn2i.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  bnn2i.Assemble(skip_zeros);
  bnn2i.Finalize(skip_zeros);
  return {*bnn1.ParallelAssemble(), *bnn2r.ParallelAssemble(), *bnn2i.ParallelAssemble()};
}

struct Att
{
  mfem::HypreParMatrix Att1;
  mfem::HypreParMatrix Att2r;
  std::optional<mfem::HypreParMatrix> Att2i;
};

inline Att GetAtt(const MaterialOperator &mat_op, mfem::ParFiniteElementSpace &nd_fespace,
                  mfem::Array<int> &attr_marker)
{
  // Stiffness matrix: Aₜₜ = (μ⁻¹ ∇ₜ x u, ∇ₜ x v) - ω² (ε u, v) = Aₜₜ₁ - ω² Aₜₜ₂.
  constexpr MaterialPropertyType MatTypeMuInv = MaterialPropertyType::INV_PERMEABILITY;
  NormalProjectedCoefficient muinv_func(
      std::make_unique<MaterialPropertyCoefficient<MatTypeMuInv>>(mat_op));
  mfem::ParBilinearForm att1(&nd_fespace);
  att1.AddBoundaryIntegrator(new mfem::CurlCurlIntegrator(muinv_func), attr_marker);
  // att1.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  att1.Assemble(skip_zeros);
  att1.Finalize(skip_zeros);

  constexpr MaterialPropertyType MatTypeEpsReal = MaterialPropertyType::PERMITTIVITY_REAL;
  MaterialPropertyCoefficient<MatTypeEpsReal> epsilon_func(mat_op);
  mfem::ParBilinearForm att2r(&nd_fespace);
  att2r.AddBoundaryIntegrator(new mfem::MixedVectorMassIntegrator(epsilon_func),
                              attr_marker);
  // att2r.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  att2r.Assemble(skip_zeros);
  att2r.Finalize(skip_zeros);

  // Contribution for loss tangent: ε => ε * (1 - i tan(δ)).
  if (!mat_op.HasLossTangent())
  {
    return {*att1.ParallelAssemble(), *att2r.ParallelAssemble()};
  }
  constexpr MaterialPropertyType MatTypeEpsImag = MaterialPropertyType::PERMITTIVITY_IMAG;
  MaterialPropertyCoefficient<MatTypeEpsImag> negepstandelta_func(mat_op);
  mfem::ParBilinearForm att2i(&nd_fespace);
  att2i.AddBoundaryIntegrator(new mfem::MixedVectorMassIntegrator(negepstandelta_func),
                              attr_marker);
  // att2i.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  att2i.Assemble(skip_zeros);
  att2i.Finalize(skip_zeros);
  return {*att1.ParallelAssemble(), *att2r.ParallelAssemble(), *att2i.ParallelAssemble()};
}

inline mfem::HypreParMatrix GetZ(mfem::ParFiniteElementSpace &fespace)
{
  // Zero matrix on ND or H1 space dofs.
  mfem::ParBilinearForm z(&fespace);
  // z.SetAssemblyLevel(mfem::AssemblyLevel::FULL);
  z.Assemble(skip_zeros);
  z.Finalize(skip_zeros);
  return *z.ParallelAssemble();
}

struct SystemMatrices
{
  petsc::PetscParMatrix A1;
  petsc::PetscParMatrix A2;
  petsc::PetscParMatrix B3;
  petsc::PetscParMatrix B4;
};

SystemMatrices
GetSystemMatrices(const mfem::HypreParMatrix &Att1, const mfem::HypreParMatrix &Att2r,
                  const std::optional<mfem::HypreParMatrix> &Att2i,
                  const mfem::HypreParMatrix &Btt, const mfem::HypreParMatrix &Btn,
                  const mfem::HypreParMatrix &Bnn1, const mfem::HypreParMatrix &Bnn2r,
                  const std::optional<mfem::HypreParMatrix> &Bnn2i,
                  const mfem::HypreParMatrix &Ztt, const mfem::HypreParMatrix &Znn,
                  const mfem::Array<int> &nd_tdof_list,
                  const mfem::Array<int> &h1_tdof_list, int nd_tdof_offset)
{
  // Construct the 2x2 block matrices for the eigenvalue problem. We pre-compute the
  // eigenvalue problem matrices such that:
  //              A = A₁ - ω² A₂, B = A + 1/Θ² B₃ - ω²/Θ² B₄.
  mfem::Array2D<const mfem::HypreParMatrix *> blocks(2, 2);
  blocks(0, 0) = &Btt;
  blocks(0, 1) = &Btn;
  blocks(1, 0) = Btn.Transpose();
  blocks(1, 1) = &Bnn1;
  std::unique_ptr<mfem::HypreParMatrix> hA1s(mfem::HypreParMatrixFromBlocks(blocks));
  auto A1s = petsc::PetscAijMatrix(*hA1s);

  blocks = nullptr;
  blocks(0, 0) = &Ztt;
  blocks(1, 1) = &Bnn2r;
  std::unique_ptr<mfem::HypreParMatrix> hA2r(mfem::HypreParMatrixFromBlocks(blocks));
  auto A2s = [&]()
  {
    if (!Bnn2i)
    {
      return petsc::PetscAijMatrix(*hA2r);
    }
    blocks(1, 1) = &*Bnn2i;
    std::unique_ptr<mfem::HypreParMatrix> hA2i(mfem::HypreParMatrixFromBlocks(blocks));
    return petsc::PetscAijMatrix(*hA2r, *hA2i);
  }();

  blocks = nullptr;
  blocks(0, 0) = &Att1;
  blocks(1, 1) = &Znn;
  std::unique_ptr<mfem::HypreParMatrix> hB3s(mfem::HypreParMatrixFromBlocks(blocks));
  auto B3s = petsc::PetscAijMatrix(*hB3s);

  blocks = nullptr;
  blocks(0, 0) = &Att2r;
  blocks(1, 1) = &Znn;
  std::unique_ptr<mfem::HypreParMatrix> hB4r(mfem::HypreParMatrixFromBlocks(blocks));
  auto B4s = [&]()
  {
    if (!Att2i)
    {
      return petsc::PetscAijMatrix(*hB4r);
    }
    blocks(0, 0) = &*Att2i;
    std::unique_ptr<mfem::HypreParMatrix> hB4i(mfem::HypreParMatrixFromBlocks(blocks));
    return petsc::PetscAijMatrix(*hB4r, *hB4i);
  }();

  // Consolidate list of local ND and H1 true dofs before extracting the respective
  // submatrices. The matrix is still distributed over the same number of processors,
  // though some are empty (PETSc handles this).
  mfem::Array<int> tdof_list;
  tdof_list.Reserve(nd_tdof_list.Size() + h1_tdof_list.Size());
  for (auto tdof : nd_tdof_list)
  {
    tdof_list.Append(tdof);
  }
  for (auto tdof : h1_tdof_list)
  {
    tdof_list.Append(tdof + nd_tdof_offset);
  }
  return {*A1s.GetSubMatrix(tdof_list, tdof_list), *A2s.GetSubMatrix(tdof_list, tdof_list),
          *B3s.GetSubMatrix(tdof_list, tdof_list), *B4s.GetSubMatrix(tdof_list, tdof_list)};
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
    if (imaginary)
    {
      gridfunc_t.imag().GetVectorValue(T, ip, V);
      V *= -kn.real();

      mfem::Vector Vn;
      gridfunc_n.real().GetGradient(T, Vn);
      V += Vn;
    }
    else
    {
      gridfunc_t.real().GetVectorValue(T, ip, V);
      V *= -kn.real();

      mfem::Vector Vn;
      gridfunc_n.imag().GetGradient(T, Vn);
      V -= Vn;
    }

    // Scale by 1/(ωμ) with μ evaluated in the neighboring element.
    mfem::Vector t(V.Size());
    V *= (1.0 / omega);
    mat_op.GetInvPermeability(mesh.GetAttribute(iel1)).Mult(V, t);
    V = std::move(t);
  }

  void SetFrequency(double w, std::complex<double> k)
  {
    omega = w;
    kn = k;
  }
};

WavePortData::WavePortData(const config::WavePortData &data, const MaterialOperator &mat_op,
                           const mfem::Array<int> &dbc_marker,
                           mfem::ParFiniteElementSpace &nd_fespace,
                           mfem::ParFiniteElementSpace &h1_fespace)
{
  excitation = data.excitation;
  mode_idx = data.mode_idx;
  d_offset = data.d_offset;
  MFEM_VERIFY(!data.attributes.empty(), "Wave port boundary found with no attributes!");
  mesh::AttrToMarker(nd_fespace.GetParMesh()->bdr_attributes.Max(), data.attributes,
                     attr_marker);

  // Construct operators for the generalized eigenvalue problem:
  //                [Aₜₜ  0] [eₜ]  = -kₙ² [Bₜₜ   Bₜₙ] [eₜ]
  //                [0   0] [eₙ]        [Bₜₙᵀ  Bₙₙ] [eₙ]
  // for the wave port of the given index. The transformed variables are related to the true
  // field by Eₜ = eₜ/kₙ and Eₙ = ieₙ. This is solved on the global mesh so the result is a
  // grid function over the entire space, not just the port boundary (so that it can be
  // queried from functions which use the global mesh).
  GetTrueDofs(dbc_marker, nd_fespace, h1_fespace, nd_attr_tdof_list, h1_attr_tdof_list);

  // Construct the system matrices. We will actually solve the shifted problem:
  //                [Bₜₜ   Bₜₙ]  [eₜ]  =  λ [Bₜₜ + 1/Θ² Aₜₜ  Bₜₙ] [eₜ]
  //                [Bₜₙᵀ  Bₙₙ] [eₙ]       [Bₜₙᵀ          Bₙₙ] [eₙ]
  // (see Lee, Sun, and Cendes, 1991). Here we have λ = Θ²/(Θ²-kₙ²), where Θ² bounds the
  // maximum kₙ² and is taken as ω² μₘₐₓ εₘₐₓ over the entire simulation domain.
  double cmin = mfem::infinity();
  for (auto attr : nd_fespace.GetParMesh()->attributes)
  {
    double s = mat_op.GetLightSpeedMin(attr);
    if (s < cmin)
    {
      cmin = s;
    }
  }
  MFEM_VERIFY(cmin > 0.0, "Invalid material speed of light detected in WavePortOperator!");
  muepsmax = 1.0 / (cmin * cmin);

  // Pre-compute problem matrices such that:
  //                A = A₁ - ω² A₂, B = A + 1/Θ² B₃ - ω²/Θ² B₄.
  // First, create parallel objects and then gather to matrices and vectors to root.
  {
    const auto &Btt = GetBtt(mat_op, nd_fespace, attr_marker);
    const auto &Btn = GetBtn(mat_op, nd_fespace, h1_fespace, attr_marker);
    const auto &[Bnn1, Bnn2r, Bnn2i] = GetBnn(mat_op, h1_fespace, attr_marker);
    const auto &[Att1, Att2r, Att2i] = GetAtt(mat_op, nd_fespace, attr_marker);
    const auto &Ztt = GetZ(nd_fespace);
    const auto &Znn = GetZ(h1_fespace);
    auto system_mat =
        GetSystemMatrices(Att1, Att2r, Att2i, Btt, Btn, Bnn1, Bnn2r, Bnn2i, Ztt, Znn,
                          nd_attr_tdof_list, h1_attr_tdof_list, nd_fespace.GetTrueVSize());
    A1 = std::make_unique<petsc::PetscParMatrix>(std::move(system_mat.A1));
    A2 = std::make_unique<petsc::PetscParMatrix>(std::move(system_mat.A2));
    B3 = std::make_unique<petsc::PetscParMatrix>(std::move(system_mat.B3));
    B4 = std::make_unique<petsc::PetscParMatrix>(std::move(system_mat.B4));
  }

  // Configure sequential vector and scatter from parallel. The original vector is created
  // to be compatible with the parallel matrix, and the scatter creates a sequential vector
  // compatible with the sequential matrix. Then, gather matrices so eigenvalue problem can
  // be solved sequentially without communication. A1/A2/B3/B4 = nullptr if !root.
  {
    bool root = Mpi::Root(A1->GetComm());
    e = std::make_unique<petsc::PetscParVector>(*A1);
    scatter =
        std::make_unique<petsc::PetscScatter>(petsc::PetscScatter::Type::TO_ZERO, *e, e0);
    A1 = A1->GetSequentialMatrix(root);
    A2 = A2->GetSequentialMatrix(root);
    B3 = B3->GetSequentialMatrix(root);
    B4 = B4->GetSequentialMatrix(root);
  }
  if (A1)
  {
    // sparsity(A2) ⊆ sparsity(A1), sparsity(B4) ⊆ sparsity(B3) ⊆ sparsity(A)
    A = std::make_unique<petsc::PetscParMatrix>(*A1);
    B = std::make_unique<petsc::PetscParMatrix>(*A1);
    A->SetSymmetric();
    B->SetSymmetric();
    A1->SetSymmetric();
    A2->SetSymmetric();
    B3->SetSymmetric();
    B4->SetSymmetric();
  }

  // Create vector for initial space (initially parallel, then scattered to root).
  {
    petsc::PetscParVector y(*e);
    GetInitialSpace(nd_attr_tdof_list.Size(), h1_attr_tdof_list.Size(), y);
    y0 = std::make_unique<petsc::PetscParVector>(*e0);
    scatter->Forward(y, *y0);
  }

  // Coefficients store references to kₙ, ω so they are updated implicitly at each new
  // solve. Also, μ⁻¹ is persistent, so no copy is OK.
  kn0 = 0.0;
  omega0 = 0.0;
  E0t = std::make_unique<mfem::ParComplexGridFunction>(&nd_fespace);
  E0n = std::make_unique<mfem::ParComplexGridFunction>(&h1_fespace);
  nxH0r_func = std::make_unique<BdrHVectorCoefficient>(*E0t, *E0n, mat_op, false);
  nxH0i_func = std::make_unique<BdrHVectorCoefficient>(*E0t, *E0n, mat_op, true);

  // Configure the eigenvalue problem solver. As for the full 3D case, the system matrices
  // are in general complex and symmetric. We supply the operators to the solver in
  // shift-inverted form and handle the back- transformation externally.
  if (A)
  {
    // Define the linear solver to be used for solving systems associated with the
    // generalized eigenvalue problem. We use PETSc's sequential sparse solvers.
    int print = 0;
    ksp = std::make_unique<KspSolver>(A->GetComm(), print, "port_");
    ksp->SetType(KspSolver::Type::CHOLESKY);  // Symmetric indefinite factorization
    ksp->SetOperator(*B);

    // Define the eigenvalue solver.
    config::EigenSolverData::Type type = config::EigenSolverData::Type::DEFAULT;
#if defined(PALACE_WITH_ARPACK) && defined(PALACE_WITH_SLEPC)
    if (type == config::EigenSolverData::Type::DEFAULT)
    {
      type = config::EigenSolverData::Type::SLEPC;
    }
#elif defined(PALACE_WITH_ARPACK)
    if (type == config::EigenSolverData::Type::SLEPC)
    {
      Mpi::Warning("SLEPc eigensolver not available, using ARPACK!\n");
    }
    type = config::EigenSolverData::Type::ARPACK;
#elif defined(PALACE_WITH_SLEPC)
    if (type == config::EigenSolverData::Type::ARPACK)
    {
      Mpi::Warning("ARPACK eigensolver not available, using SLEPc!\n");
    }
    type = config::EigenSolverData::Type::SLEPC;
#else
#error "Wave port solver requires building with ARPACK or SLEPc!"
#endif
    if (type == config::EigenSolverData::Type::ARPACK)
    {
#if defined(PALACE_WITH_ARPACK)
      eigen = std::unique_ptr<EigenSolverBase>(new arpack::ArpackEPSSolver(print));
#endif
    }
    else  // config::EigenSolverData::Type::SLEPC
    {
#if defined(PALACE_WITH_SLEPC)
      eigen =
          std::unique_ptr<EigenSolverBase>(new slepc::SlepcEPSSolver(A->GetComm(), print));
      auto *slepc = dynamic_cast<slepc::SlepcEPSSolver *>(eigen.get());
      slepc->SetProblemType(slepc::SlepcEigenSolver::ProblemType::GEN_NON_HERMITIAN);
      slepc->SetType(slepc::SlepcEigenSolver::Type::KRYLOVSCHUR);
#endif
    }
    constexpr double tol = 1.0e-6;
    eigen->SetLinearSolver(*ksp);
    eigen->SetWhichEigenpairs(EigenSolverBase::WhichType::LARGEST_MAGNITUDE);
    eigen->SetNumModes(mode_idx, std::max(2 * mode_idx + 1, 5));
    eigen->SetTol(tol);
  }
}

void WavePortData::GetTrueDofs(const mfem::Array<int> &dbc_marker,
                               mfem::ParFiniteElementSpace &nd_fespace,
                               mfem::ParFiniteElementSpace &h1_fespace,
                               mfem::Array<int> &nd_tdof_list,
                               mfem::Array<int> &h1_tdof_list)
{
  // Ensures no duplicates in the attribute list for this port index (this would imply a
  // mistake in the configuration file). We can, however, have multiple unique ports with
  // shared boundary attributes.
  nd_fespace.GetEssentialTrueDofs(attr_marker, nd_tdof_list);
  h1_fespace.GetEssentialTrueDofs(attr_marker, h1_tdof_list);
  int nd_tdofs = nd_tdof_list.Size();
  int h1_tdofs = h1_tdof_list.Size();

  // Mark all ND and H1 dofs on the port, then unmark PEC boundaries.
  mfem::Array<int> nd_tdof_marker(nd_fespace.GetTrueVSize()),
      h1_tdof_marker(h1_fespace.GetTrueVSize()), nd_dbc_tdof_list, h1_dbc_tdof_list;
  nd_tdof_marker = 0;
  h1_tdof_marker = 0;
  nd_fespace.GetEssentialTrueDofs(dbc_marker, nd_dbc_tdof_list);
  h1_fespace.GetEssentialTrueDofs(dbc_marker, h1_dbc_tdof_list);
  for (auto tdof : nd_tdof_list)
  {
    nd_tdof_marker[tdof] = 1;
  }
  for (auto tdof : nd_dbc_tdof_list)
  {
    nd_tdof_marker[tdof] = 0;
  }
  for (auto tdof : h1_tdof_list)
  {
    h1_tdof_marker[tdof] = 1;
  }
  for (auto tdof : h1_dbc_tdof_list)
  {
    h1_tdof_marker[tdof] = 0;
  }

  // Convert back to a list.
  nd_tdof_list.DeleteAll();
  nd_tdof_list.Reserve(nd_tdofs);
  for (int i = 0; i < nd_tdof_marker.Size(); i++)
  {
    if (nd_tdof_marker[i])
    {
      nd_tdof_list.Append(i);
    }
  }
  h1_tdof_list.DeleteAll();
  h1_tdof_list.Reserve(h1_tdofs);
  for (int i = 0; i < h1_tdof_marker.Size(); i++)
  {
    if (h1_tdof_marker[i])
    {
      h1_tdof_list.Append(i);
    }
  }
}

void WavePortData::GetInitialSpace(int nt, int nn, petsc::PetscParVector &y0)
{
  // Initial space chosen as such that B v₀ = y₀, with y₀ = [y₀ₜ, 0, ... 0]ᵀ ⟂ null(A)
  // (with Aₜₜ nonsingular). See Lee, Sun, and Cendes, 1991 for reference.
  // Note: When the eigenvalue solver uses a standard ℓ²-inner product instead of B-inner
  // product(since we use a general non-Hermitian solver due to complex symmetric B), then
  // we just use v0 = y0 directly.
  MFEM_VERIFY(y0.GetSize() == nt + nn, "Invalid vector size!");
  y0.SetRandomReal();
  PetscScalar *py0 = y0.GetArray();
  // for (int i = 0; i < nt; i++)     { py0[i] = 1.0; }
  for (int i = nt; i < nt + nn; i++)
  {
    py0[i] = 0.0;
  }
  y0.RestoreArray(py0);
}

std::complex<double> WavePortData::Solve(petsc::PetscParVector &y0,
                                         petsc::PetscParVector &e0,
                                         petsc::PetscParVector &e,
                                         petsc::PetscScatter &scatter)
{
  double eig[2];
  if (A)  // Only on root
  {
    // The y0 and e0 vectors are still parallel vectors, but with all data on root. We want
    // true sequential vectors.
    PetscScalar *pe0 = e0.GetArray();
    petsc::PetscParVector e0s(e0.GetSize(), pe0);

    // Set starting vector.
    {
      PetscScalar *py0 = y0.GetArray();
      petsc::PetscParVector y0s(y0.GetSize(), py0);
      eigen->SetInitialSpace(y0s);
      y0.RestoreArray(py0);
    }

#if 0
    // Alternatively, use B-orthogonal initial space. Probably want to call SetBMat for
    // the eigensolver in this case.
    {
      PetscScalar *py0 = y0.GetArray();
      petsc::PetscParVector y0s(y0.GetSize(), py0);
      petsc::PetscParVector v0s(y0s);
      ksp->Mult(y0s, v0s);
      eigen->SetInitialSpace(v0s);
      y0.RestoreArray(py0);
    }
#endif

    // Solve (operators have been set in constructor).
    int num_conv = 0;
    eigen->SetOperators(*A, *B, EigenSolverBase::ScaleType::NONE);
    num_conv = eigen->Solve();
    MFEM_VERIFY(num_conv >= mode_idx, "Wave port eigensolver did not converge!");
    eigen->GetEigenvalue(mode_idx - 1, eig[0], eig[1]);
    eigen->GetEigenvector(mode_idx - 1, e0s);
    e0.RestoreArray(pe0);
  }

  // Scatter the result to all processors.
  scatter.Reverse(e0, e);
  Mpi::Broadcast(2, eig, 0, e.GetComm());
  return {eig[0], eig[1]};
}

void WavePortData::Initialize(double omega)
{
  if (omega == omega0)
  {
    return;
  }

  // Use pre-computed matrices to construct and solve the generalized eigenvalue problem for
  // the desired wave port mode.
  double theta2 = muepsmax * omega * omega;
  if (A)
  {
    MFEM_VERIFY(A1 && A2 && B3 && B4 && A && B,
                "Boundary mode eigenvalue problem operators uninitialized for solve!");
    A->Scale(0.0);
    A->AXPY(1.0, *A1, petsc::PetscParMatrix::NNZStructure::SAME);
    A->AXPY(-omega * omega, *A2, petsc::PetscParMatrix::NNZStructure::SUBSET);
    B->Scale(0.0);
    B->AXPY(1.0, *A, petsc::PetscParMatrix::NNZStructure::SAME);
    B->AXPY(1.0 / theta2, *B3, petsc::PetscParMatrix::NNZStructure::SUBSET);
    B->AXPY(-omega * omega / theta2, *B4, petsc::PetscParMatrix::NNZStructure::SUBSET);
  }

  // Configure and solve the eigenvalue problem for the desired boundary mode.
  std::complex<double> lambda = Solve(*y0, *e0, *e, *scatter);

  // Extract the eigenmode solution and postprocess. The extracted eigenvalue is λ =
  // Θ²/(Θ²-kₙ²).
  MFEM_VERIFY(lambda.real() > 1.0 / (1.0 - 1.0e-2),
              "Computed wave port mode is or is very close to being evanescent "
                  << "(λ = " << lambda << ")!");
  kn0 = std::sqrt(theta2 - theta2 / lambda);
  omega0 = omega;
  dynamic_cast<BdrHVectorCoefficient &>(*nxH0r_func).SetFrequency(omega0, kn0);
  dynamic_cast<BdrHVectorCoefficient &>(*nxH0i_func).SetFrequency(omega0, kn0);

  mfem::Vector etr(nd_attr_tdof_list.Size()), eti(nd_attr_tdof_list.Size()),
      enr(h1_attr_tdof_list.Size()), eni(h1_attr_tdof_list.Size());
  MFEM_VERIFY(e->GetSize() == etr.Size() + enr.Size(),
              "Unexpected vector size in wave port eigenmode solver!");
  e->GetToVectors(etr, eti, 0, nd_attr_tdof_list.Size());
  e->GetToVectors(enr, eni, nd_attr_tdof_list.Size(),
                  nd_attr_tdof_list.Size() + h1_attr_tdof_list.Size());

  // Re-expand from restricted boundary dofs to true dofs and transform back to true
  // electric field variables: Eₜ = eₜ/kₙ and Eₙ = ieₙ.
  auto &nd_fespace = *E0t->ParFESpace();
  auto &h1_fespace = *E0n->ParFESpace();
  mfem::Vector E0tr(nd_fespace.GetTrueVSize()), E0ti(nd_fespace.GetTrueVSize()),
      E0nr(h1_fespace.GetTrueVSize()), E0ni(h1_fespace.GetTrueVSize());
  E0tr = 0.0;
  E0ti = 0.0;
  E0nr = 0.0;
  E0ni = 0.0;
  std::complex<double> ookn = 1.0 / kn0;
  for (int i = 0; i < nd_attr_tdof_list.Size(); i++)
  {
    E0tr(nd_attr_tdof_list[i]) = ookn.real() * etr(i) - ookn.imag() * eti(i);
    E0ti(nd_attr_tdof_list[i]) = ookn.imag() * etr(i) + ookn.real() * eti(i);
  }
  for (int i = 0; i < h1_attr_tdof_list.Size(); i++)
  {
    E0nr(h1_attr_tdof_list[i]) = -eni(i);
    E0ni(h1_attr_tdof_list[i]) = enr(i);
  }
  E0t->real().SetFromTrueDofs(E0tr);  // Parallel distribute
  E0t->imag().SetFromTrueDofs(E0ti);
  E0n->real().SetFromTrueDofs(E0nr);
  E0n->imag().SetFromTrueDofs(E0ni);

  // Normalize grid functions to a chosen polarization direction and unit power, |E x H⋆| ⋅
  // n, integrated over the port surface (+n is the direction of propagation). The n x H
  // coefficients are updated implicitly as the only store references to the Et, En grid
  // functions as well as kₙ, ω.
  {
    // Choose a (rather arbitrary) sign constraint: @ t = 0, 1ᵀ E > 0 when integrated over
    // the port surface. This at least makes results for the same port consistent between
    // frequencies/meshes.
    mfem::Vector ones(nd_fespace.GetParMesh()->SpaceDimension());
    ones = 1.0;
    mfem::VectorConstantCoefficient tdir(ones);
    mfem::ConstantCoefficient ndir(1.0);
    mfem::ParLinearForm sut(&nd_fespace), sun(&h1_fespace);
    sut.AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(tdir), attr_marker);
    sun.AddBoundaryIntegrator(new BoundaryLFIntegrator(ndir), attr_marker);
    sut.UseFastAssembly(false);
    sun.UseFastAssembly(false);
    sut.Assemble();
    sun.Assemble();
    if (sut(E0t->real()) + sun(E0n->real()) < 0.0)
    {
      E0t->real().Neg();  // This updates the n x H coefficients depending on Et, En
      E0t->imag().Neg();
      E0n->real().Neg();
      E0n->imag().Neg();
    }
  }
  {
    // |E x H⋆| ⋅ n = |E ⋅ (-n x H⋆)|
    sr = std::make_unique<mfem::ParLinearForm>(&nd_fespace);
    si = std::make_unique<mfem::ParLinearForm>(&nd_fespace);
    sr->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(*nxH0r_func), attr_marker);
    si->AddBoundaryIntegrator(new VectorFEBoundaryLFIntegrator(*nxH0i_func), attr_marker);
    sr->UseFastAssembly(false);
    si->UseFastAssembly(false);
    sr->Assemble();
    si->Assemble();
    std::complex<double> s0(-(*sr)(E0t->real()) - (*si)(E0t->imag()),
                            -(*sr)(E0t->imag()) + (*si)(E0t->real()));
    double scale = 1.0 / std::sqrt(std::abs(s0));
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
                                            const MaterialOperator &mat_op,
                                            const std::map<int, int> &local_to_shared) const
{
  // Compute port power, (E x H) ⋅ n = E ⋅ (-n x H), integrated over the port surface
  // using the computed E and H = μ⁻¹ B fields. The linear form is reconstructed from
  // scratch each time due to changing H. The BdrCurrentVectorCoefficient computes -n x H,
  // where n is an outward normal.
  auto &nd_fespace = *E.ParFESpace();
  BdrCurrentVectorCoefficient nxHr_func(B.real(), mat_op, local_to_shared);
  BdrCurrentVectorCoefficient nxHi_func(B.imag(), mat_op, local_to_shared);
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
    ports.try_emplace(idx, data, mat_op, dbc_marker, nd_fespace, h1_fespace);
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
        // Print header at first solve.
        if (data.GetA() && data.GetB())
        {
          Mpi::Print(" Number of global unknowns for port {:d}: {}\n", idx,
                     data.GetA()->GetGlobalNumRows());
          Mpi::Print("  A: NNZ = {:d}, norm = {:e}\n", data.GetA()->NNZ(),
                     data.GetA()->NormF());
          Mpi::Print("  B: NNZ = {:d}, norm = {:e}\n", data.GetB()->NNZ(),
                     data.GetB()->NormF());
        }
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
    fbi.AddCoefficient(std::make_unique<MaterialPropertyCoefficient<MatType>>(
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
