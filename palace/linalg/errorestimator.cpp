// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "errorestimator.hpp"

#include <cstring>
#include <functional>
#include <limits>
#include <string>
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

#include "fem/qfunctions/22/hcurlhdiv_error_22_qf.h"
#include "fem/qfunctions/hcurlhdiv_complex_error_qf.h"
#include "fem/qfunctions/hcurlhdiv_error_qf.h"
#include "fem/qfunctions/l2h1_error_qf.h"

PalacePragmaDiagnosticPop

namespace palace
{

namespace
{

template <OperatorType OperType>
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

template <OperatorType OperType>
auto BuildLevelParOperator(std::unique_ptr<Operator> &&a, const FiniteElementSpace &fespace)
{
  return BuildLevelParOperator<OperType>(std::move(a), fespace, fespace);
}

template <OperatorType OperType>
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
FluxProjector<VecType>::FluxProjector(
    const MaterialPropertyCoefficient &coeff,
    const FiniteElementSpaceHierarchy &smooth_fespaces,
    const FiniteElementSpace &rhs_fespace, double tol, int max_it, int print, bool use_mg,
    const std::vector<MaterialPropertyCoefficient> &coeff_k)
{
  BlockTimer bt(Timer::CONSTRUCT_ESTIMATOR);
  const auto &smooth_fespace = smooth_fespaces.GetFinestFESpace();
  // Detect scalar FE spaces (H1, L2) vs vector FE spaces (ND, RT) based on map type.
  const int dim = smooth_fespace.Dimension();
  const auto map_type = smooth_fespace.GetFEColl().GetMapType(dim);
  const bool scalar_flux =
      (map_type == mfem::FiniteElement::VALUE || map_type == mfem::FiniteElement::INTEGRAL);
  {
    constexpr bool skip_zeros = false;
    BilinearForm m(smooth_fespace);
    if (scalar_flux)
    {
      m.AddDomainIntegrator<MassIntegrator>();
    }
    else
    {
      m.AddDomainIntegrator<VectorFEMassIntegrator>();
    }
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
  auto AssembleFlux = [&](const MaterialPropertyCoefficient &f)
  {
    // Flux operators are always partially assembled. Every process assembles each term,
    // empty where it has none of its elements, so that the parallel operators match.
    BilinearForm flux(rhs_fespace, smooth_fespace);
    if (!f.empty())
    {
      if (scalar_flux)
      {
        flux.AddDomainIntegrator<MassIntegrator>(f);
      }
      else
      {
        flux.AddDomainIntegrator<VectorFEMassIntegrator>(f);
      }
    }
    return BuildLevelParOperator<OperType>(flux.PartialAssemble(), rhs_fespace,
                                           smooth_fespace);
  };
  Flux = AssembleFlux(coeff);
  for (const auto &f : coeff_k)
  {
    Flux_k.push_back(AssembleFlux(f));
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

template <typename VecType>
void FluxProjector<VecType>::Mult(const VecType &x, VecType &y,
                                  const std::vector<std::complex<double>> &w) const
  requires std::is_same_v<VecType, ComplexVector>
{
  BlockTimer bt(Timer::SOLVE_ESTIMATOR);
  MFEM_VERIFY(w.size() == Flux_k.size(), "Invalid number of flux projection weights!");
  MFEM_ASSERT(x.Size() == Flux->Width() && y.Size() == rhs.Size(),
              "Invalid vector dimensions for FluxProjector::Mult!");
  Flux->Mult(x, rhs);
  for (std::size_t k = 0; k < Flux_k.size(); k++)
  {
    if (w[k] != 0.0)
    {
      Flux_k[k]->AddMult(x, rhs, w[k]);
    }
  }
  ksp->Mult(rhs, y);
}

CeedContextList::~CeedContextList()
{
  for (auto &[ceed, c] : ctx)
  {
    PalaceCeedCall(ceed, CeedQFunctionContextDestroy(&c));
  }
}

void CeedContextList::Append(Ceed ceed, CeedQFunctionContext c)
{
  PalacePragmaOmp(critical(CeedContextListAppend))
  {
    ctx.emplace_back(ceed, c);
  }
}

void CeedContextList::Update(const void *data, std::size_t size) const
{
  for (const auto &[ceed, c] : ctx)
  {
    std::size_t ctx_size;
    PalaceCeedCall(ceed, CeedQFunctionContextGetContextSize(c, &ctx_size));
    MFEM_VERIFY(ctx_size == size, "Invalid size for a QFunction context update!");
    void *ctx_data;
    PalaceCeedCall(ceed, CeedQFunctionContextGetData(c, CEED_MEM_HOST, &ctx_data));
    std::memcpy(ctx_data, data, size);
    PalaceCeedCall(ceed, CeedQFunctionContextRestoreData(c, &ctx_data));
  }
}

namespace
{

// Integrate the elemental errors of the libCEED operator, whose passive inputs u_1, u_2,
// ... are re-pointed at the given grid functions (shared by every element geometry
// sub-operator of a thread), accumulating into the estimates.
void IntegrateErrors(const ceed::Operator &integ_op, const std::vector<const Vector *> &u,
                     Vector &estimates)
{
  PalacePragmaOmp(parallel if (ceed::internal::NumCeeds() > 1))
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];

    // We need to update the state of the underlying libCEED vectors to indicate that the
    // data has changed. Each thread has its own set of vectors, referencing the same
    // underlying data, and every sub-operator of the thread's composite (one per element
    // geometry type) shares that set, so re-pointing the vectors of the first sub-operator
    // updates all of them.
    std::vector<CeedVector> u_vec(u.size());
    {
      CeedInt nsub_ops;
      CeedOperator *sub_ops;
      PalaceCeedCall(
          ceed, CeedOperatorCompositeGetNumSub(integ_op[utils::GetThreadNum()], &nsub_ops));
      PalaceCeedCall(
          ceed, CeedOperatorCompositeGetSubList(integ_op[utils::GetThreadNum()], &sub_ops));
      MFEM_ASSERT(nsub_ops > 0, "Unexpected empty libCEED composite operator!");
      for (std::size_t k = 0; k < u.size(); k++)
      {
        CeedOperatorField field;
        PalaceCeedCall(ceed,
                       CeedOperatorGetFieldByName(
                           sub_ops[0], ("u_" + std::to_string(k + 1)).c_str(), &field));
        PalaceCeedCall(ceed, CeedOperatorFieldGetVector(field, &u_vec[k]));
        ceed::InitCeedVector(*u[k], ceed, &u_vec[k], false);
      }
    }

    // Each thread writes to non-overlapping entries of the estimates vector.
    CeedVector estimates_vec;
    ceed::InitCeedVector(estimates, ceed, &estimates_vec);
    PalaceCeedCall(ceed,
                   CeedOperatorApplyAdd(integ_op[utils::GetThreadNum()], CEED_VECTOR_NONE,
                                        estimates_vec, CEED_REQUEST_IMMEDIATE));

    // Cleanup (the field vectors are reference counted, the operator keeps its own
    // references).
    PalaceCeedCall(ceed, CeedVectorDestroy(&estimates_vec));
    for (auto &v : u_vec)
    {
      PalaceCeedCall(ceed, CeedVectorDestroy(&v));
    }
  }
}

// Context of the complex constitutive flux error QFunctions: the ε', ε'', and W⁻¹ tables.
std::vector<CeedIntScalar> FrequencyDependentContext(const MaterialOperator &mat_op,
                                                     const mfem::DenseTensor &eps_r,
                                                     const mfem::DenseTensor &eps_i,
                                                     const mfem::DenseTensor &winv)
{
  const int dim = mat_op.SpaceDimension();
  MaterialPropertyCoefficient eps_r_func(mat_op.GetAttributeToMaterial(), eps_r),
      eps_i_func(mat_op.GetAttributeToMaterial(), eps_i),
      winv_func(mat_op.GetAttributeToMaterial(), winv);
  auto ctx = ceed::PopulateCoefficientContext(dim, &eps_r_func, dim, &eps_i_func);
  const auto ctx_w = ceed::PopulateCoefficientContext(dim, &winv_func);
  ctx.insert(ctx.end(), ctx_w.begin(), ctx_w.end());
  return ctx;
}

Vector ZeroEstimates(const Mesh &mesh)
{
  Vector estimates(mesh.GetNE());
  estimates.UseDevice(true);
  estimates = 0.0;
  return estimates;
}

// Elemental flux errors with the passive inputs (discontinuous, smooth flux), for each of
// the real and imaginary parts of complex fields (the squares add).
template <typename VecType>
Vector IntegrateRealErrors(const Mesh &mesh, const ceed::Operator &integ_op,
                           const VecType &F_gf, const VecType &G_gf)
{
  auto estimates = ZeroEstimates(mesh);
  if (mesh.GetNE() == 0)  // Possible for a submesh, with an empty libCEED operator
  {
    return estimates;
  }
  if constexpr (std::is_same_v<VecType, ComplexVector>)
  {
    IntegrateErrors(integ_op, {&F_gf.Real(), &G_gf.Real()}, estimates);
    IntegrateErrors(integ_op, {&F_gf.Imag(), &G_gf.Imag()}, estimates);
  }
  else
  {
    IntegrateErrors(integ_op, {&F_gf, &G_gf}, estimates);
  }
  return estimates;
}

template <typename VecType>
void Prolong(const FiniteElementSpace &fespace, const VecType &x, VecType &x_gf)
{
  if constexpr (std::is_same_v<VecType, ComplexVector>)
  {
    fespace.GetProlongationMatrix()->Mult(x.Real(), x_gf.Real());
    fespace.GetProlongationMatrix()->Mult(x.Imag(), x_gf.Imag());
  }
  else
  {
    fespace.GetProlongationMatrix()->Mult(x, x_gf);
  }
}

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
  Prolong(fespace, F, F_gf);
  Prolong(smooth_fespace, G, G_gf);

  // Use libCEED operators to perform the error estimate integration over each element.
  return IntegrateRealErrors(fespace.GetMesh(), integ_op, F_gf, G_gf);
}

// The energy norms of the recovered flux of the last estimate, with a zero discontinuous
// flux input. This overwrites the field workspace F_gf.
template <typename VecType>
Vector RecoveredNorms(const Mesh &mesh, const ceed::Operator &integ_op, VecType &F_gf,
                      const VecType &G_gf)
{
  F_gf = 0.0;
  return IntegrateRealErrors(mesh, integ_op, F_gf, G_gf);
}

// The scale of the squared indicators normalized by the total energy, with the factor of
// 1/2 in the energy. A nonpositive or nonfinite energy (possible for a general
// PoleResidue model with a negative stored energy) is replaced by the energy norm of the
// recovered fluxes, so that the indicators stay positive and comparable across samples.
double IndicatorScale(double Et, const std::function<double()> &recovered_energy,
                      MPI_Comm comm)
{
  if (Et > 0.0 && std::isfinite(Et))
  {
    return 0.5 / Et;
  }
  const double Er = recovered_energy();
  Mpi::Warning(comm,
               "Error indicator normalization: total field energy {:.3e} is not positive "
               "(a frequency-dependent material with a negative stored energy?); using the "
               "energy norm of the recovered fluxes, {:.3e}, instead!\n",
               Et, Er);
  return (Er > 0.0) ? 0.5 / Er : 1.0;
}

double GlobalSum(const Vector &v, MPI_Comm comm)
{
  double sum = linalg::LocalSum(v);
  Mpi::GlobalSum(1, &sum, comm);
  return sum;
}

}  // namespace

template <typename VecType>
GradFluxErrorEstimator<VecType>::GradFluxErrorEstimator(
    const MaterialOperator &mat_op, FiniteElementSpace &nd_fespace,
    FiniteElementSpaceHierarchy &rt_fespaces, double tol, int max_it, int print,
    bool use_mg)
  : nd_fespace(nd_fespace), rt_fespace(rt_fespaces.GetFinestFESpace()), mat_op(mat_op),
    fd_materials(
        [&mat_op]()
        {
          std::vector<std::size_t> m;
          for (std::size_t i = 0; i < mat_op.NumFrequencyDependentPermittivityMaterials();
               i++)
          {
            if (mat_op.HasFrequencyDependentPermittivityModel(i))
            {
              m.push_back(i);
            }
          }
          return m;
        }()),
    fd_loss_tangent(!fd_materials.empty() && mat_op.HasLossTangent()),
    projector(MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                          mat_op.GetPermittivityReal()),
              rt_fespaces, nd_fespace, tol, max_it, print, use_mg,
              [&]()
              {
                // The weighted terms of the complex constitutive flux: Im{ε} of the
                // constant materials (weight i), then the unit flux of each
                // frequency-dependent material (weight ε(iω) − ε∞).
                std::vector<MaterialPropertyCoefficient> coeff_k;
                if (fd_loss_tangent)
                {
                  coeff_k.emplace_back(mat_op.GetAttributeToMaterial(),
                                       mat_op.GetPermittivityImag());
                }
                for (auto m : fd_materials)
                {
                  auto &f = coeff_k.emplace_back(mat_op.MaxCeedAttribute());
                  const auto attrs = mat_op.GetCeedAttributes(
                      mat_op.GetFrequencyDependentPermittivityAttributes(m));
                  if (attrs.Size() > 0)
                  {
                    f.AddMaterialProperty(attrs, 1.0);
                  }
                }
                return coeff_k;
              }()),
    integ_op(nd_fespace.GetMesh().GetNE(), nd_fespace.GetVSize()),
    fd_warned(mat_op.NumFrequencyDependentPermittivityMaterials(), false),
    E_gf(nd_fespace.GetVSize()), D(rt_fespace.GetTrueVSize()), D_gf(rt_fespace.GetVSize())
{
  E_gf.UseDevice(true);
  D.UseDevice(true);
  D_gf.UseDevice(true);
  const bool complex_route = !fd_materials.empty();
  if constexpr (!std::is_same_v<VecType, ComplexVector>)
  {
    MFEM_VERIFY(!complex_route, "Real-valued flux error estimation does not support "
                                "frequency-dependent permittivity!");
  }

  // Construct the libCEED operator used for integrating the element-wise error. The
  // discontinuous flux is ε E = ε ∇V.
  const auto &mesh = nd_fespace.GetMesh();
  PalacePragmaOmp(parallel if (ceed::internal::NumCeeds() > 1))
  {
    Ceed ceed = ceed::internal::GetCeedObjects()[utils::GetThreadNum()];

    // Create libCEED vector wrappers for use with libCEED operators. One set of vectors
    // per Ceed context is shared by the sub-operators of every element geometry type, so
    // that the estimate can re-point the passive inputs once through any sub-operator (in
    // particular between the real and imaginary parts of complex fields).
    CeedVector E_gf_vec, D_gf_vec, E_gf_vec_i = nullptr, D_gf_vec_i = nullptr;
    if constexpr (std::is_same_v<VecType, ComplexVector>)
    {
      ceed::InitCeedVector(E_gf.Real(), ceed, &E_gf_vec);
      ceed::InitCeedVector(D_gf.Real(), ceed, &D_gf_vec);
      if (complex_route)
      {
        ceed::InitCeedVector(E_gf.Imag(), ceed, &E_gf_vec_i);
        ceed::InitCeedVector(D_gf.Imag(), ceed, &D_gf_vec_i);
      }
    }
    else
    {
      ceed::InitCeedVector(E_gf, ceed, &E_gf_vec);
      ceed::InitCeedVector(D_gf, ceed, &D_gf_vec);
    }
    for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
    {
      // Only integrate over domain elements (not on the boundary).
      if (mfem::Geometry::Dimension[geom] < mesh.Dimension())
      {
        continue;
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

      ceed::CeedQFunctionInfo info;
      info.assemble_q_data = false;
      info.trial_ops = info.test_ops = ceed::EvalMode::Interp;
      CeedOperator sub_op;
      if (complex_route)
      {
        // Complex constitutive flux: context tables ε', ε'', W⁻¹, set at each frequency
        // (here with ε∞ and W = Re{ε}, only for the layout).
        const auto ctx = FrequencyDependentContext(
            mat_op, mat_op.GetPermittivityReal(), mat_op.GetPermittivityImag(),
            linalg::MatrixPow(mat_op.GetPermittivityReal(), -1.0));
        switch (10 * mesh.SpaceDimension() + mesh.Dimension())
        {
          case 22:
            info.apply_qf = f_apply_hcurlhdiv_complex_error_22;
            info.apply_qf_path =
                PalaceQFunctionRelativePath(f_apply_hcurlhdiv_complex_error_22_loc);
            break;
          case 33:
            info.apply_qf = f_apply_hcurlhdiv_complex_error_33;
            info.apply_qf_path =
                PalaceQFunctionRelativePath(f_apply_hcurlhdiv_complex_error_33_loc);
            break;
          default:
            MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                       << mesh.Dimension() << ", " << mesh.SpaceDimension()
                       << ") for GradFluxErrorEstimator!");
        }
        CeedQFunctionContext sub_ctx;
        ceed::AssembleCeedElementErrorIntegrator(
            info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed,
            {{E_gf_vec, nd_restr, nd_basis},
             {E_gf_vec_i, nd_restr, nd_basis},
             {D_gf_vec, rt_restr, rt_basis},
             {D_gf_vec_i, rt_restr, rt_basis}},
            mesh_elem_restr, data.geom_data, data.geom_data_restr, &sub_op, &sub_ctx);
        integ_ctx.Append(ceed, sub_ctx);
      }
      else
      {
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
        switch (10 * mesh.SpaceDimension() + mesh.Dimension())
        {
          case 22:
            info.apply_qf = f_apply_hcurlhdiv_error_22;
            info.apply_qf_path =
                PalaceQFunctionRelativePath(f_apply_hcurlhdiv_error_22_loc);
            break;
          case 33:
            info.apply_qf = f_apply_hcurlhdiv_error_33;
            info.apply_qf_path =
                PalaceQFunctionRelativePath(f_apply_hcurlhdiv_error_33_loc);
            break;
          default:
            MFEM_ABORT("Invalid value of (dim, space_dim) = ("
                       << mesh.Dimension() << ", " << mesh.SpaceDimension()
                       << ") for GradFluxErrorEstimator!");
        }
        ceed::AssembleCeedElementErrorIntegrator(
            info, (void *)ctx.data(), ctx.size() * sizeof(CeedIntScalar), ceed, E_gf_vec,
            D_gf_vec, nd_restr, rt_restr, nd_basis, rt_basis, mesh_elem_restr,
            data.geom_data, data.geom_data_restr, &sub_op);
      }
      integ_op.AddSubOperator(sub_op);  // Sub-operator owned by ceed::Operator

      // Element restriction is owned by the operator.
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_elem_restr));
    }

    // Passive input vectors are owned (shared) by the sub-operators.
    for (auto *v : {&E_gf_vec, &D_gf_vec, &E_gf_vec_i, &D_gf_vec_i})
    {
      if (*v)
      {
        PalaceCeedCall(ceed, CeedVectorDestroy(v));
      }
    }
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  integ_op.Finalize();
}

template <typename VecType>
void GradFluxErrorEstimator<VecType>::SetFrequency(double omega) const
  requires std::is_same_v<VecType, ComplexVector>
{
  if (omega == fd_omega)
  {
    return;
  }
  MFEM_VERIFY(omega > 0.0, "Error estimation for a frequency-dependent material requires "
                           "a positive frequency!");
  // The positive stored-energy permittivity of each model weights the indicator. A
  // general model with W(ω) ≤ 0 keeps the solve and its reported (signed) energy, and the
  // indicator alone uses ε∞ > 0.
  std::vector<double> W(mat_op.NumFrequencyDependentPermittivityMaterials(), 0.0);
  for (auto m : fd_materials)
  {
    W[m] = mat_op.GetFrequencyDependentStoredEnergyPermittivity(m, omega);
    if (W[m] <= 0.0)
    {
      const double eps_inf = mat_op.GetFrequencyDependentPermittivityHighFrequency(m);
      if (!fd_warned[m])
      {
        Mpi::Warning(mat_op.GetMesh().GetComm(),
                     "Material {:d} has a nonpositive stored-energy permittivity {:.3e} at "
                     "ω = {:.3e} (a general PoleResidue model with a negative stored "
                     "energy); the error indicator uses ε∞ = {:.3e} for it at this and any "
                     "later such frequency!\n",
                     m, W[m], omega, eps_inf);
        fd_warned[m] = true;
      }
      W[m] = eps_inf;
    }
  }
  auto winv = mat_op.MapFrequencyDependentMaterials(
      linalg::MatrixPow(mat_op.GetPermittivityReal(), -1.0),
      [&W](std::size_t m, mfem::DenseMatrix &M)
      {
        M = 0.0;
        for (int d = 0; d < M.Height(); d++)
        {
          M(d, d) = 1.0 / W[m];
        }
      });
  // A process without elements (of a submesh) has no contexts, and no material table to
  // build them from.
  if (!integ_ctx.Empty())
  {
    const auto ctx =
        FrequencyDependentContext(mat_op, mat_op.EvaluatePermittivityReal(omega),
                                  mat_op.EvaluatePermittivityImag(omega), winv);
    integ_ctx.Update(ctx.data(), ctx.size() * sizeof(CeedIntScalar));
  }
  fd_omega = omega;
}

template <typename VecType>
Vector GradFluxErrorEstimator<VecType>::ComputeEstimates(const VecType &E) const
{
  MFEM_VERIFY(fd_materials.empty(),
              "Error estimation for a frequency-dependent material requires the frequency "
              "of the field!");
  return ComputeErrorEstimates(E, E_gf, D, D_gf, nd_fespace, rt_fespace, projector,
                               integ_op);
}

template <typename VecType>
Vector GradFluxErrorEstimator<VecType>::ComputeEstimates(const VecType &E,
                                                         double omega) const
  requires std::is_same_v<VecType, ComplexVector>
{
  if (fd_materials.empty())
  {
    return ComputeEstimates(E);
  }
  BlockTimer bt(Timer::ESTIMATION);
  SetFrequency(omega);

  // Recover the complete constitutive flux D = ε(iω) E: the constant Re{ε} flux, plus the
  // weighted terms.
  std::vector<std::complex<double>> w;
  if (fd_loss_tangent)
  {
    w.emplace_back(0.0, 1.0);
  }
  for (auto m : fd_materials)
  {
    w.push_back(mat_op.EvaluateFrequencyDependentPermittivity(m, {0.0, omega}) -
                mat_op.GetFrequencyDependentPermittivityHighFrequency(m));
  }
  projector.Mult(E, D, w);
  Prolong(nd_fespace, E, E_gf);
  Prolong(rt_fespace, D, D_gf);

  const auto &mesh = nd_fespace.GetMesh();
  auto estimates = ZeroEstimates(mesh);
  if (mesh.GetNE() > 0)
  {
    IntegrateErrors(integ_op, {&E_gf.Real(), &E_gf.Imag(), &D_gf.Real(), &D_gf.Imag()},
                    estimates);
  }
  return estimates;
}

template <typename VecType>
Vector GradFluxErrorEstimator<VecType>::RecoveredFluxNorms() const
{
  const auto &mesh = nd_fespace.GetMesh();
  if constexpr (std::is_same_v<VecType, ComplexVector>)
  {
    if (!fd_materials.empty())
    {
      E_gf = 0.0;
      auto norms = ZeroEstimates(mesh);
      if (mesh.GetNE() > 0)
      {
        IntegrateErrors(integ_op, {&E_gf.Real(), &E_gf.Imag(), &D_gf.Real(), &D_gf.Imag()},
                        norms);
      }
      return norms;
    }
  }
  return RecoveredNorms(mesh, integ_op, E_gf, D_gf);
}

template <typename VecType>
void GradFluxErrorEstimator<VecType>::AddErrorIndicator(const VecType &E, double Et,
                                                        ErrorIndicator &indicator) const
{
  auto estimates = ComputeEstimates(E);
  linalg::Sqrt(estimates, (Et > 0.0) ? 0.5 / Et : 1.0);  // Correct factor of 1/2 in energy
  indicator.AddIndicator(estimates);
}

template <typename VecType>
void GradFluxErrorEstimator<VecType>::AddErrorIndicator(const VecType &E, double Et,
                                                        double omega,
                                                        ErrorIndicator &indicator) const
  requires std::is_same_v<VecType, ComplexVector>
{
  auto estimates = ComputeEstimates(E, omega);
  const auto comm = nd_fespace.GetComm();
  linalg::Sqrt(
      estimates,
      IndicatorScale(
          Et, [&]() { return 0.5 * GlobalSum(RecoveredFluxNorms(), comm); }, comm));
  indicator.AddIndicator(estimates);
}

template <typename VecType>
CurlFluxErrorEstimator<VecType>::CurlFluxErrorEstimator(
    const MaterialOperator &mat_op, FiniteElementSpace &rt_fespace,
    FiniteElementSpaceHierarchy &nd_fespaces, double tol, int max_it, int print,
    bool use_mg)
  : rt_fespace(rt_fespace), nd_fespace(nd_fespaces.GetFinestFESpace()),
    projector(MaterialPropertyCoefficient(mat_op.GetAttributeToMaterial(),
                                          mat_op.GetCurlCurlInvPermeability()),
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

    // Create libCEED vector wrappers for use with libCEED operators. One pair of vectors
    // per Ceed context is shared by the sub-operators of every element geometry type, so
    // that ComputeErrorEstimates can re-point the passive inputs once through any
    // sub-operator (in particular between the real and imaginary parts of complex fields).
    CeedVector B_gf_vec, H_gf_vec;
    if constexpr (std::is_same_v<VecType, ComplexVector>)
    {
      ceed::InitCeedVector(B_gf.Real(), ceed, &B_gf_vec);
      ceed::InitCeedVector(H_gf.Real(), ceed, &H_gf_vec);
    }
    else
    {
      ceed::InitCeedVector(B_gf, ceed, &B_gf_vec);
      ceed::InitCeedVector(H_gf, ceed, &H_gf_vec);
    }
    for (const auto &[geom, data] : mesh.GetCeedGeomFactorData(ceed))
    {
      // Only integrate over domain elements (not on the boundary).
      if (mfem::Geometry::Dimension[geom] < mesh.Dimension())
      {
        continue;
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
      // In 2D, the curl is scalar so we use scalar (1x1) √μ⁻¹ and √μ coefficients.
      const bool scalar_curl = (mesh.Dimension() == 2);
      const auto &muinv_tensor =
          scalar_curl ? mat_op.GetCurlCurlInvPermeability() : mat_op.GetInvPermeability();
      auto mat_invsqrtmu = linalg::MatrixSqrt(muinv_tensor);
      auto mat_sqrtmu = linalg::MatrixPow(muinv_tensor, -0.5);
      MaterialPropertyCoefficient invsqrtmu_func(mat_op.GetAttributeToMaterial(),
                                                 mat_invsqrtmu);
      MaterialPropertyCoefficient sqrtmu_func(mat_op.GetAttributeToMaterial(), mat_sqrtmu);
      const int coeff_dim = scalar_curl ? 1 : mesh.SpaceDimension();
      auto ctx = ceed::PopulateCoefficientContext(coeff_dim, &invsqrtmu_func, coeff_dim,
                                                  &sqrtmu_func);

      // Assemble the libCEED operator. Inputs: B (for discontinuous flux), then smooth
      // flux.
      ceed::CeedQFunctionInfo info;
      info.assemble_q_data = false;
      switch (10 * mesh.SpaceDimension() + mesh.Dimension())
      {
        case 22:
          // 2D scalar curl: L2 (B) × H1 (smooth H) error
          info.apply_qf = f_apply_l2h1_error;
          info.apply_qf_path = PalaceQFunctionRelativePath(f_apply_l2h1_error_loc);
          break;
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

      // Element restriction is owned by the operator.
      PalaceCeedCall(ceed, CeedElemRestrictionDestroy(&mesh_elem_restr));
    }

    // Passive input vectors are owned (shared) by the sub-operators.
    PalaceCeedCall(ceed, CeedVectorDestroy(&B_gf_vec));
    PalaceCeedCall(ceed, CeedVectorDestroy(&H_gf_vec));
  }

  // Finalize the operator (call CeedOperatorCheckReady).
  integ_op.Finalize();
}

template <typename VecType>
Vector CurlFluxErrorEstimator<VecType>::ComputeEstimates(const VecType &B) const
{
  return ComputeErrorEstimates(B, B_gf, H, H_gf, rt_fespace, nd_fespace, projector,
                               integ_op);
}

template <typename VecType>
Vector CurlFluxErrorEstimator<VecType>::RecoveredFluxNorms() const
{
  return RecoveredNorms(rt_fespace.GetMesh(), integ_op, B_gf, H_gf);
}

template <typename VecType>
void CurlFluxErrorEstimator<VecType>::AddErrorIndicator(const VecType &B, double Et,
                                                        ErrorIndicator &indicator) const
{
  auto estimates = ComputeEstimates(B);
  linalg::Sqrt(estimates, (Et > 0.0) ? 0.5 / Et : 1.0);  // Correct factor of 1/2 in energy
  indicator.AddIndicator(estimates);
}

namespace
{

// Sum of the squared gradient and curl flux indicators, normalized by the total energy.
void AddFluxErrorIndicator(Vector grad_estimates, const Vector &curl_estimates,
                           double scale, ErrorIndicator &indicator)
{
  grad_estimates += curl_estimates;
  linalg::Sqrt(grad_estimates, scale);  // Correct factor of 1/2 in energy
  indicator.AddIndicator(grad_estimates);
}

}  // namespace

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
  AddFluxErrorIndicator(grad_estimator.ComputeEstimates(E),
                        curl_estimator.ComputeEstimates(B), (Et > 0.0) ? 0.5 / Et : 1.0,
                        indicator);
}

template <typename VecType>
void TimeDependentFluxErrorEstimator<VecType>::AddErrorIndicator(
    const VecType &E, const VecType &B, double Et, double omega,
    ErrorIndicator &indicator) const
  requires std::is_same_v<VecType, ComplexVector>
{
  auto grad_estimates = grad_estimator.ComputeEstimates(E, omega);
  auto curl_estimates = curl_estimator.ComputeEstimates(B);
  const auto comm = grad_estimator.nd_fespace.GetComm();
  const double scale = IndicatorScale(
      Et,
      [&]()
      {
        return 0.5 * (GlobalSum(grad_estimator.RecoveredFluxNorms(), comm) +
                      GlobalSum(curl_estimator.RecoveredFluxNorms(), comm));
      },
      comm);
  AddFluxErrorIndicator(std::move(grad_estimates), curl_estimates, scale, indicator);
}

template <typename VecType>
BoundaryModeFluxErrorEstimator<VecType>::BoundaryModeFluxErrorEstimator(
    const MaterialOperator &mat_op, FiniteElementSpaceHierarchy &nd_fespaces,
    FiniteElementSpaceHierarchy &rt_fespaces, FiniteElementSpace &curl_fespace,
    FiniteElementSpaceHierarchy &h1_fespaces, double tol, int max_it, int print,
    bool use_mg)
  : grad_estimator(mat_op, nd_fespaces.GetFinestFESpace(), rt_fespaces, tol, max_it, print,
                   use_mg),
    curl_estimator(mat_op, curl_fespace, h1_fespaces, tol, max_it, print, use_mg)
{
}

template <typename VecType>
void BoundaryModeFluxErrorEstimator<VecType>::AddErrorIndicator(
    const VecType &E, const VecType &B, double Et, ErrorIndicator &indicator) const
{
  AddFluxErrorIndicator(grad_estimator.ComputeEstimates(E),
                        curl_estimator.ComputeEstimates(B), (Et > 0.0) ? 0.5 / Et : 1.0,
                        indicator);
}

template <typename VecType>
void BoundaryModeFluxErrorEstimator<VecType>::AddErrorIndicator(
    const VecType &E, const VecType &B, double Et, double omega,
    ErrorIndicator &indicator) const
  requires std::is_same_v<VecType, ComplexVector>
{
  auto grad_estimates = grad_estimator.ComputeEstimates(E, omega);
  auto curl_estimates = curl_estimator.ComputeEstimates(B);
  const auto comm = grad_estimator.nd_fespace.GetComm();
  const double scale = IndicatorScale(
      Et,
      [&]()
      {
        return 0.5 * (GlobalSum(grad_estimator.RecoveredFluxNorms(), comm) +
                      GlobalSum(curl_estimator.RecoveredFluxNorms(), comm));
      },
      comm);
  AddFluxErrorIndicator(std::move(grad_estimates), curl_estimates, scale, indicator);
}

template class FluxProjector<Vector>;
template class FluxProjector<ComplexVector>;
template class GradFluxErrorEstimator<Vector>;
template class GradFluxErrorEstimator<ComplexVector>;
template class CurlFluxErrorEstimator<Vector>;
template class CurlFluxErrorEstimator<ComplexVector>;
template class TimeDependentFluxErrorEstimator<Vector>;
template class TimeDependentFluxErrorEstimator<ComplexVector>;
template class BoundaryModeFluxErrorEstimator<Vector>;
template class BoundaryModeFluxErrorEstimator<ComplexVector>;

}  // namespace palace
