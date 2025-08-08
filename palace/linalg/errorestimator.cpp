// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"

#include <limits>
#include "fem/bilinearform.hpp"
#include "fem/integrator.hpp"
#include "fem/libceed/ceed.hpp"
#include "fem/libceed/coefficient.hpp"
#include "fem/libceed/integrator.hpp"
#include "linalg/amg.hpp"
#include "linalg/densematrix.hpp"
#include "linalg/gmg.hpp"
#include "linalg/iterative.hpp"
#include "linalg/jacobi.hpp"
#include "linalg/rap.hpp"
#include "models/materialoperator.hpp"
#include "utils/communication.hpp"
#include "utils/diagnostic.hpp"
#include "utils/omp.hpp"
#include "utils/timer.hpp"

PalacePragmaDiagnosticPush
PalacePragmaDiagnosticDisableUnused

#include "fem/qfunctions/hcurlhdiv_error_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

template <typename OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&a,
                           const FiniteElementSpace &trial_fespace,
                           const FiniteElementSpace &test_fespace);

template <>
auto BuildLevelParOperator<Operator>(std::unique_ptr<Operator> &&a,
                                     const FiniteElementSpace &trial_fespace,
                                     const FiniteElementSpace &test_fespace)
{
  return std::make_unique<ParOperator>(std::move(a), trial_fespace, test_fespace, false);
}

template <>
auto BuildLevelParOperator<ComplexOperator>(std::unique_ptr<Operator> &&a,
                                            const FiniteElementSpace &trial_fespace,
                                            const FiniteElementSpace &test_fespace)
{
  return std::make_unique<ComplexParOperator>(std::move(a), nullptr, trial_fespace,
                                              test_fespace, false);
}

template <typename OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&a, const FiniteElementSpace &fespace)
{
  return BuildLevelParOperator<OperType>(std::move(a), fespace, fespace);
}

template <typename OperType>
auto ConfigureLinearSolver(const FiniteElementSpaceHierarchy &fespaces, double tol,
                           int max_it, int print, bool use_mg)
{
  // The system matrix for the projection is real, SPD and diagonally dominant.
  std::unique_ptr<Solver<OperType>> pc;
  if (!use_mg)
  {
    // Use eigenvalue estimate to compute optimal Jacobi damping parameter.
    pc = std::make_unique<JacobiSmoother<OperType>>(fespaces.GetFinestFESpace().GetComm(),
                                                    0.0);
  }
  else
  {
    auto amg = std::make_unique<BoomerAmgSolver>(1, 1, true, 0);
    amg->SetStrengthThresh(0.8);  // More coarsening to save memory
    if (fespaces.GetNumLevels() > 1)
    {
      const int mg_smooth_order = 2;  // Smooth order independent of FE space order
      pc = std::make_unique<GeometricMultigridSolver<OperType>>(
          fespaces.GetFinestFESpace().GetComm(),
          std::make_unique<MfemWrapperSolver<OperType>>(std::move(amg)),
          fespaces.GetProlongationOperators(), nullptr, 1, 1, mg_smooth_order, 1.0, 0.0,
          true);
    }
    else
    {
      pc = std::make_unique<MfemWrapperSolver<OperType>>(std::move(amg));
    }
  }
  auto pcg =
      std::make_unique<CgSolver<OperType>>(fespaces.GetFinestFESpace().GetComm(), print);
  pcg->SetInitialGuess(false);
  pcg->SetRelTol(tol);
  pcg->SetAbsTol(std::numeric_limits<double>::epsilon());
  pcg->SetMaxIter(max_it);
  return std::make_unique<BaseKspSolver<OperType>>(std::move(pcg), std::move(pc));
}

}  // namespace

template <typename VecType>
FluxProjector<VecType>::FluxProjector(const MaterialPropertyCoefficient &coeff,
                                      const FiniteElementSpaceHierarchy &smooth_fespaces,
                                      const FiniteElementSpace &rhs_fespace, double tol,
                                      int max_it, int print, bool use_mg)
{
  BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
  const auto &smooth_fespace = smooth_fespaces.GetFinestFESpace();
  {
    constexpr bool skip_zeros = false;
    BilinearForm m(smooth_fespace);
    m.AddDomainIntegrator<VectorFEMassIntegrator>();
    if (!use_mg)
    {
      M = BuildLevelParOperator<OperType>(m.Assemble(skip_zeros), smooth_fespace);
    }
    else
    {
      auto m_vec = m.Assemble(smooth_fespaces, skip_zeros);
      auto M_mg =
          std::make_unique<BaseMultigridOperator<OperType>>(smooth_fespaces.GetNumLevels());
      for (std::size_t l = 0; l < smooth_fespaces.GetNumLevels(); l++)
      {
        const auto &fespace_l = smooth_fespaces.GetFESpaceAtLevel(l);
        M_mg->AddOperator(BuildLevelParOperator<OperType>(std::move(m_vec[l]), fespace_l));
      }
      M = std::move(M_mg);
    }
  }
  {
    // Flux operator is always partially assembled.
    BilinearForm flux(rhs_fespace, smooth_fespace);
    flux.AddDomainIntegrator<VectorFEMassIntegrator>(coeff);
    Flux = BuildLevelParOperator<OperType>(flux.PartialAssemble(), rhs_fespace,
                                           smooth_fespace);
  }
  ksp = ConfigureLinearSolver<OperType>(smooth_fespaces, tol, max_it, print, use_mg);
  ksp->SetOperators(*M, *M);
  rhs.SetSize(smooth_fespace.GetTrueVSize());
  rhs.UseDevice(true);
}

template <typename VecType>
void FluxProjector<VecType>::Mult(const VecType &x, VecType &y) const
{
  BlockTimer bt(Timer::SOLVE_ESTIMATOR);
  MFEM_ASSERT(x.Size() == Flux->Width() && y.Size() == rhs.Size(),
              "Invalid vector dimensions for FluxProjector::Mult!");
  // Mpi::Print(" Computing smooth flux recovery (projection) for error estimation\n");
  Flux->Mult(x, rhs);
  ksp->Mult(rhs, y);
}

namespace
{

template <typename VecType>
Vector ComputeErrorEstimates(const VecType &F, VecType &F_gf, VecType &G, VecType &G_gf,
                             const FiniteElementSpace &fespace,
                             const FiniteElementSpace &smooth_fespace,
                             const FluxProjector<VecType> &projector,
                             const ceed::Operator &integ_op)
{
  // Compute the projection of the discontinuous flux onto the smooth finite element space
  // (recovery) and populate the corresponding grid functions.
  BlockTimer bt(Timer::ESTIMATION);
  projector.Mult(F, G);
  if constexpr (std::is_same<VecType, ComplexVector>::value)
  {
    fespace.GetProlongationMatrix()->Mult(F.Real(), F_gf.Real());
    fespace.GetProlongationMatrix()->Mult(F.Imag(), F_gf.Imag());
    smooth_fespace.GetProlongationMatrix()->Mult(G.Real(), G_gf.Real());
    smooth_fespace.GetProlongationMatrix()->Mult(G.Imag(), G_gf.Imag());
  }
  else
  {
    fespace.GetProlongationMatrix()->Mult(F, F_gf);
    smooth_fespace.GetProlongationMatrix()->Mult(G, G_gf);
  }

  // Use libCEED operators to perform the error estimate integration over each element.
  const auto &mesh = fespace.GetMesh();
  Vector estimates(mesh.GetNE());
  estimates.UseDevice(true);
  estimates = 0.0;
  PalacePragmaOmp(parallel if (ceed::internal::NumCeeds() > 1))
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];

    // We need to update the state of the underlying libCEED vectors to indicate that the
    // data has changed. Each thread has it's own vector, referencing the same underlying
    // data.
    CeedVector F_gf_vec, G_gf_vec;
    {
      CeedInt nsub_ops;
      CeedOperator *sub_ops;
      PalaceCeedCall(
          ceed, CeedCompositeOperatorGetNumSub(integ_op[utils::GetThreadNum()], &nsub_ops));
      PalaceCeedCall(
          ceed, CeedCompositeOperatorGetSubList(integ_op[utils::GetThreadNum()], &sub_ops));
      MFEM_ASSERT(nsub_ops > 0, "Unexpected empty libCEED composite operator!");
      CeedOperatorField field;
      PalaceCeedCall(ceed, CeedOperatorGetFieldByName(sub_ops[0], "u_1", &field));
      PalaceCeedCall(ceed, CeedOperatorFieldGetVector(field, &F_gf_vec));
      PalaceCeedCall(ceed, CeedOperatorGetFieldByName(sub_ops[0], "u_2", &field));
      PalaceCeedCall(ceed, CeedOperatorFieldGetVector(field, &G_gf_vec));
      if constexpr (std::is_same<VecType, ComplexVector>::value)
      {
        ceed::InitCeedVector(F_gf.Real(), ceed, &F_gf_vec, false);
        ceed::InitCeedVector(G_gf.Real(), ceed, &G_gf_vec, false);
      }
      else
      {
        ceed::InitCeedVector(F_gf, ceed, &F_gf_vec, false);
        ceed::InitCeedVector(G_gf, ceed, &G_gf_vec, false);
      }
    }

    // Each thread writes to non-overlapping entries of the estimates vector.
    CeedVector estimates_vec;
    ceed::InitCeedVector(estimates, ceed, &estimates_vec);

    // Do the integration (both input vectors are passive). For the complex case, add sum of
    // squares of real and imaginary parts to the estimates before square root.
    PalaceCeedCall(ceed,
                   CeedOperatorApplyAdd(integ_op[utils::GetThreadNum()], CEED_VECTOR_NONE,
                                        estimates_vec, CEED_REQUEST_IMMEDIATE));
    if constexpr (std::is_same<VecType, ComplexVector>::value)
    {
      ceed::InitCeedVector(F_gf.Imag(), ceed, &F_gf_vec, false);
      ceed::InitCeedVector(G_gf.Imag(), ceed, &G_gf_vec, false);
      PalaceCeedCall(ceed,
                     CeedOperatorApplyAdd(integ_op[utils::GetThreadNum()], CEED_VECTOR_NONE,
                                          estimates_vec, CEED_REQUEST_IMMEDIATE));
    }

    // Cleanup.
    PalaceCeedCall(ceed, CeedVectorDestroy(&estimates_vec));
  }

  return estimates;
}

}  // namespace

template <typename VecType>
GradFluxErrorEstimator<VecType>::GradFluxErrorEstimator(
    const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
    FiniteElementSpaceHierarchy &rt_fespaces, double tol, int max_it, int print,
    bool use_mg)
  : nd_fespace(nd_fespace), rt_fespace(rt_fespaces.GetFinestFESpace()),
    projector(MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                          mat_op.GetPermittivityReal()),
              rt_fespaces, nd_fespace, tol, max_it, print, use_mg),
    integ_op(nd_fespace.GetMesh().GetNE(), nd_fespace.GetVSize()),
    E_gf(nd_fespace.GetVSize()), D(rt_fespace.GetTrueVSize()), D_gf(rt_fespace.GetVSize())
{
  E_gf.UseDevice(true);
  D.UseDevice(true);
  D_gf.UseDevice(true);

  // Construct the libCEED operator used for integrating the element-wise error. The
  // discontinuous flux is ε E = ε ∇V.
  const auto &mesh = nd_fespace.GetMesh();
  PalacePragmaOmp(parallel if (ceed::internal::NumCeeds() > 1))
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];
    for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
    {
      // Only integrate over domain elements (not on the boundary).
      if (mfem::Geometry::Dimension[geom] < mesh.Dimension())
      {
        continue;
      }

      // Create libCEED vector wrappers for use with libCEED operators.
      CeedVector E_gf_vec, D_gf_vec;
      if constexpr (std::is_same<VecType, ComplexVector>::value)
      {
        ceed::InitCeedVector(E_gf.Real(), ceed, &E_gf_vec);
        ceed::InitCeedVector(D_gf.Real(), ceed, &D_gf_vec);
      }
      else
      {
        ceed::InitCeedVector(E_gf, ceed, &E_gf_vec);
        ceed::InitCeedVector(D_gf, ceed, &D_gf_vec);
      }

      // Construct mesh element restriction for elements of this element geometry type.
      CeedElemRestriction mesh_elem_restr;
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                               ceed, static_cast<CeedInt>(data.indices.size()), 1, 1,
                               mesh.GetNE(), mesh.GetNE(), CEED_MEM_HOST, CEED_USE_POINTER,
                               data.indices.data(), &mesh_elem_restr));

      // Element restriction and basis objects for inputs.
      CeedElemRestriction nd_restr =
          nd_fespace.GetCeedElemRestriction(ceed, geom, data.indices);
      CeedElemRestriction rt_restr =
          rt_fespace.GetCeedElemRestriction(ceed, geom, data.indices);
      CeedBasis nd_basis = nd_fespace.GetCeedBasis(ceed, geom);
      CeedBasis rt_basis = rt_fespace.GetCeedBasis(ceed, geom);

      // Construct coefficient for discontinuous flux, then smooth flux.
      auto mat_sqrtepsilon = linalg::MatrixSqrt(mat_op.GetPermittivityReal());
      auto mat_invsqrtepsilon = linalg::MatrixPow(mat_op.GetPermittivityReal(), -0.5);
      MaterialPropertyCoefficient sqrtepsilon_func(mat_op.GetAttributeToMaterial(),
                                                   mat_sqrtepsilon);
      MaterialPropertyCoefficient invsqrtepsilon_func(mat_op.GetAttributeToMaterial(),
                                                      mat_invsqrtepsilon);
      auto ctx =
          ceed::PopulateCoefficientContext(mesh.SpaceDimension(), &sqrtepsilon_func,
                                           mesh.SpaceDimension(), &invsqrtepsilon_func);

      // Assemble the libCEED operator. Inputs: E (for discontinuous flux), then smooth
      // flux.
      ceed::CeedQFunctionInfo info;
      info.assemble_q_data = false;
      switch (10 * mesh.SpaceDimension() + mesh.Dimension())
      {
        case 22:
          info.apply_qf = f_apply_hcurlhdiv_error_22;
          info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlhdiv_error_22_loc);
          break;
        case 33:
          info.apply_qf = f_apply_hcurlhdiv_error_33;
          info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hcurlhdiv_error_33_loc);
          break;
        default:
          MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                     << mesh.Dimension() << ", " << mesh.SpaceDimension()
                     << ") for GradFluxErrorEstimator!");
      }
      info.trial_ops = info.test_ops = ceed::EvalMode::Interp;

      CeedOperator sub_op;
      ceed::AssembleCeedElementErrorIntegrator(
          info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed, E_gf_vec,
          D_gf_vec, nd_restr, rt_restr, nd_basis, rt_basis, mesh_elem_restr, data.geom_data,
          data.geom_data_restr, &sub_op);
      integ_op.AddSubOperator(sub_op);  // Sub-operator owned by ceed::Operator

      // Element restriction and passive input vectors are owned by the operator.
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_elem_restr));
      PalaceCeedCall(ceed, CeedVectorDestroy(&E_gf_vec));
      PalaceCeedCall(ceed, CeedVectorDestroy(&D_gf_vec));
    }
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  integ_op.Finalize();
}

template <typename VecType>
void GradFluxErrorEstimator<VecType>::AddErrorIndicator(const VecType &E, double Et,
                                                        ErrorIndicator &indicator) const
{
  auto estimates =
      ComputeErrorEstimates(E, E_gf, D, D_gf, nd_fespace, rt_fespace, projector, integ_op);
  linalg::Sqrt(estimates, (Et > 0.0) ? 0.5 / Et : 1.0);  // Correct factor of 1/2 in energy
  indicator.AddIndicator(estimates);
}

template <typename VecType>
CurlFluxErrorEstimator<VecType>::CurlFluxErrorEstimator(
    const MaterialOperator &mat_op, FiniteElementSpace &rt_fespace,
    FiniteElementSpaceHierarchy &nd_fespaces, double tol, int max_it, int print,
    bool use_mg)
  : rt_fespace(rt_fespace), nd_fespace(nd_fespaces.GetFinestFESpace()),
    projector(MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                          mat_op.GetInvPermeability()),
              nd_fespaces, rt_fespace, tol, max_it, print, use_mg),
    integ_op(nd_fespace.GetMesh().GetNE(), rt_fespace.GetVSize()),
    B_gf(rt_fespace.GetVSize()), H(nd_fespace.GetTrueVSize()), H_gf(nd_fespace.GetVSize())
{
  B_gf.UseDevice(true);
  H.UseDevice(true);
  H_gf.UseDevice(true);

  // Construct the libCEED operator used for integrating the element-wise error. The
  // discontinuous flux is μ⁻¹ B ≃ μ⁻¹ ∇ × E.
  const auto &mesh = rt_fespace.GetMesh();
  PalacePragmaOmp(parallel if (ceed::internal::NumCeeds() > 1))
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];
    for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
    {
      // Only integrate over domain elements (not on the boundary).
      if (mfem::Geometry::Dimension[geom] < mesh.Dimension())
      {
        continue;
      }

      // Create libCEED vector wrappers for use with libCEED operators.
      CeedVector B_gf_vec, H_gf_vec;
      if constexpr (std::is_same<VecType, ComplexVector>::value)
      {
        ceed::InitCeedVector(B_gf.Real(), ceed, &B_gf_vec);
        ceed::InitCeedVector(H_gf.Real(), ceed, &H_gf_vec);
      }
      else
      {
        ceed::InitCeedVector(B_gf, ceed, &B_gf_vec);
        ceed::InitCeedVector(H_gf, ceed, &H_gf_vec);
      }

      // Construct mesh element restriction for elements of this element geometry type.
      CeedElemRestriction mesh_elem_restr;
      PalaceCeedCall(ceed, CeedElemRestrictionCreate(
                               ceed, static_cast<CeedInt>(data.indices.size()), 1, 1,
                               mesh.GetNE(), mesh.GetNE(), CEED_MEM_HOST, CEED_USE_POINTER,
                               data.indices.data(), &mesh_elem_restr));

      // Element restriction and basis objects for inputs.
      CeedElemRestriction rt_restr =
          rt_fespace.GetCeedElemRestriction(ceed, geom, data.indices);
      CeedElemRestriction nd_restr =
          nd_fespace.GetCeedElemRestriction(ceed, geom, data.indices);
      CeedBasis rt_basis = rt_fespace.GetCeedBasis(ceed, geom);
      CeedBasis nd_basis = nd_fespace.GetCeedBasis(ceed, geom);

      // Construct coefficient for discontinuous flux, then smooth flux.
      auto mat_invsqrtmu = linalg::MatrixSqrt(mat_op.GetInvPermeability());
      auto mat_sqrtmu = linalg::MatrixPow(mat_op.GetInvPermeability(), -0.5);
      MaterialPropertyCoefficient invsqrtmu_func(mat_op.GetAttributeToMaterial(),
                                                 mat_invsqrtmu);
      MaterialPropertyCoefficient sqrtmu_func(mat_op.GetAttributeToMaterial(), mat_sqrtmu);
      auto ctx = ceed::PopulateCoefficientContext(mesh.SpaceDimension(), &invsqrtmu_func,
                                                  mesh.SpaceDimension(), &sqrtmu_func);

      // Assemble the libCEED operator. Inputs: B (for discontinuous flux), then smooth
      // flux. Currently only supports 3D, since curl in 2D requires special treatment.
      ceed::CeedQFunctionInfo info;
      info.assemble_q_data = false;
      switch (10 * mesh.SpaceDimension() + mesh.Dimension())
      {
        case 33:
          info.apply_qf = f_apply_hdivhcurl_error_33;
          info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_hdivhcurl_error_33_loc);
          break;
        default:
          MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                     << mesh.Dimension() << ", " << mesh.SpaceDimension()
                     << ") for CurlFluxErrorEstimator!");
      }
      info.trial_ops = info.test_ops = ceed::EvalMode::Interp;

      CeedOperator sub_op;
      ceed::AssembleCeedElementErrorIntegrator(
          info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed, B_gf_vec,
          H_gf_vec, rt_restr, nd_restr, rt_basis, nd_basis, mesh_elem_restr, data.geom_data,
          data.geom_data_restr, &sub_op);
      integ_op.AddSubOperator(sub_op);  // Sub-operator owned by ceed::Operator

      // Element restriction and passive input vectors are owned by the operator.
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_elem_restr));
      PalaceCeedCall(ceed, CeedVectorDestroy(&B_gf_vec));
      PalaceCeedCall(ceed, CeedVectorDestroy(&H_gf_vec));
    }
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  integ_op.Finalize();
}

template <typename VecType>
void CurlFluxErrorEstimator<VecType>::AddErrorIndicator(const VecType &B, double Et,
                                                        ErrorIndicator &indicator) const
{
  auto estimates =
      ComputeErrorEstimates(B, B_gf, H, H_gf, rt_fespace, nd_fespace, projector, integ_op);
  linalg::Sqrt(estimates, (Et > 0.0) ? 0.5 / Et : 1.0);  // Correct factor of 1/2 in energy
  indicator.AddIndicator(estimates);
}

template <typename VecType>
TimeDependentFluxErrorEstimator<VecType>::TimeDependentFluxErrorEstimator(
    const MaterialOperator &mat_op, FiniteElementSpaceHierarchy &nd_fespaces,
    FiniteElementSpaceHierarchy &rt_fespaces, double tol, int max_it, int print,
    bool use_mg)
  : grad_estimator(mat_op, nd_fespaces.GetFinestFESpace(), rt_fespaces, tol, max_it, print,
                   use_mg),
    curl_estimator(mat_op, rt_fespaces.GetFinestFESpace(), nd_fespaces, tol, max_it, print,
                   use_mg)
{
}

template <typename VecType>
void TimeDependentFluxErrorEstimator<VecType>::AddErrorIndicator(
    const VecType &E, const VecType &B, double Et, ErrorIndicator &indicator) const
{
  auto grad_estimates =
      ComputeErrorEstimates(E, grad_estimator.E_gf, grad_estimator.D, grad_estimator.D_gf,
                            grad_estimator.nd_fespace, grad_estimator.rt_fespace,
                            grad_estimator.projector, grad_estimator.integ_op);
  auto curl_estimates =
      ComputeErrorEstimates(B, curl_estimator.B_gf, curl_estimator.H, curl_estimator.H_gf,
                            curl_estimator.rt_fespace, curl_estimator.nd_fespace,
                            curl_estimator.projector, curl_estimator.integ_op);
  grad_estimates += curl_estimates;  // Sum of squares
  linalg::Sqrt(grad_estimates,
               (Et > 0.0) ? 0.5 / Et : 1.0);  // Correct factor of 1/2 in energy
  indicator.AddIndicator(grad_estimates);
}

template class FluxProjector<Vector>;
template class FluxProjector<ComplexVector>;
template class GradFluxErrorEstimator<Vector>;
template class GradFluxErrorEstimator<ComplexVector>;
template class CurlFluxErrorEstimator<Vector>;
template class CurlFluxErrorEstimator<ComplexVector>;
template class TimeDependentFluxErrorEstimator<Vector>;
template class TimeDependentFluxErrorEstimator<ComplexVector>;

}  // namespace palace
