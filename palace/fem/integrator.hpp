// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_INTEGRATOR_HPP
#define PALACE_FEM_INTEGRATOR_HPP

#include <vector>
#include <mfem.hpp>

// Forward declarations of libCEED objects.
typedef struct Ceed_private *Ceed;
typedef struct CeedOperator_private *CeedOperator;

namespace palace
{

//
// Classes which implement or extend bilinear and linear form integrators.
//

namespace fem
{

// Helper functions for creating an integration rule to exactly integrate polynomials of
// order 2p + w + q_extra.
inline int GetDefaultIntegrationOrder(const mfem::FiniteElement &trial_fe,
                                      const mfem::FiniteElement &test_fe,
                                      const mfem::ElementTransformation &T, int q_extra = 0)
{

  // //XX TODO DEBUG
  // std::cout << "Mesh geometry order: " << T.OrderW() << "\n";

  return trial_fe.GetOrder() + test_fe.GetOrder() + T.OrderW() + q_extra;
}

inline int GetDefaultIntegrationOrder(const mfem::FiniteElementSpace &trial_fespace,
                                      const mfem::FiniteElementSpace &test_fespace,
                                      int q_extra = 0)
{
  // Every process is guaranteed to have at least one element, and assumes no variable order
  // spaces are used.
  mfem::Mesh &mesh = *trial_fespace.GetMesh();
  const mfem::FiniteElement &trial_fe = *trial_fespace.GetFE(0);
  const mfem::FiniteElement &test_fe = *test_fespace.GetFE(0);
  const mfem::ElementTransformation &T = *mesh.GetElementTransformation(0);
  return GetDefaultIntegrationOrder(trial_fe, test_fe, T, q_extra);
}

}  // namespace fem

// Base class for libCEED-based bilinear form integrators.
class BilinearFormIntegrator
{
public:
  virtual ~BilinearFormIntegrator() = default;

  virtual void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) = 0;

  virtual void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                                const mfem::FiniteElementSpace &test_fespace,
                                const mfem::IntegrationRule &ir,
                                const std::vector<int> &indices, Ceed ceed,
                                CeedOperator *op, CeedOperator *op_t) = 0;
};

// Integrator for a(u, v) = (Q u, v) for H1 elements (also for vector (H1)ᵈ spaces).
class MassIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  MassIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  MassIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  MassIntegrator(mfem::VectorCoefficient &VQ) : Q(nullptr), VQ(&VQ), MQ(nullptr) {}
  MassIntegrator(mfem::MatrixCoefficient &MQ) : Q(nullptr), VQ(nullptr), MQ(&MQ) {}

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Q u, v) for vector finite elements.
class VectorFEMassIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  VectorFEMassIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  VectorFEMassIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  VectorFEMassIntegrator(mfem::VectorCoefficient &VQ) : Q(nullptr), VQ(&VQ), MQ(nullptr) {}
  VectorFEMassIntegrator(mfem::MatrixCoefficient &MQ) : Q(nullptr), VQ(nullptr), MQ(&MQ) {}

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Q curl u, curl v) for Nedelec elements.
class CurlCurlIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  CurlCurlIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  CurlCurlIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  CurlCurlIntegrator(mfem::VectorCoefficient &VQ) : Q(nullptr), VQ(&VQ), MQ(nullptr) {}
  CurlCurlIntegrator(mfem::MatrixCoefficient &MQ) : Q(nullptr), VQ(nullptr), MQ(&MQ) {}

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Qc curl u, curl v) + (Qm u, v) for Nedelec elements.
class CurlCurlMassIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Qc, *Qm;
  mfem::VectorCoefficient *VQc, *VQm;
  mfem::MatrixCoefficient *MQc, *MQm;

public:
  CurlCurlMassIntegrator()
    : Qc(nullptr), Qm(nullptr), VQc(nullptr), VQm(nullptr), MQc(nullptr), MQm(nullptr)
  {
  }
  CurlCurlMassIntegrator(mfem::Coefficient &Qc, mfem::Coefficient &Qm)
    : Qc(&Qc), Qm(&Qm), VQc(nullptr), VQm(nullptr), MQc(nullptr), MQm(nullptr)
  {
  }
  CurlCurlMassIntegrator(mfem::Coefficient &Qc, mfem::VectorCoefficient &VQm)
    : Qc(&Qc), Qm(nullptr), VQc(nullptr), VQm(&VQm), MQc(nullptr), MQm(nullptr)
  {
  }
  CurlCurlMassIntegrator(mfem::Coefficient &Qc, mfem::MatrixCoefficient &MQm)
    : Qc(&Qc), Qm(nullptr), VQc(nullptr), VQm(nullptr), MQc(nullptr), MQm(&MQm)
  {
  }
  CurlCurlMassIntegrator(mfem::VectorCoefficient &VQc, mfem::Coefficient &Qm)
    : Qc(nullptr), Qm(&Qm), VQc(&VQc), VQm(nullptr), MQc(nullptr), MQm(nullptr)
  {
  }
  CurlCurlMassIntegrator(mfem::VectorCoefficient &VQc, mfem::VectorCoefficient &VQm)
    : Qc(nullptr), Qm(nullptr), VQc(&VQc), VQm(&VQm), MQc(nullptr), MQm(nullptr)
  {
  }
  CurlCurlMassIntegrator(mfem::VectorCoefficient &VQc, mfem::MatrixCoefficient &MQm)
    : Qc(nullptr), Qm(nullptr), VQc(&VQc), VQm(nullptr), MQc(nullptr), MQm(&MQm)
  {
  }
  CurlCurlMassIntegrator(mfem::MatrixCoefficient &MQc, mfem::Coefficient &Qm)
    : Qc(nullptr), Qm(&Qm), VQc(nullptr), VQm(nullptr), MQc(&MQc), MQm(nullptr)
  {
  }
  CurlCurlMassIntegrator(mfem::MatrixCoefficient &MQc, mfem::VectorCoefficient &VQm)
    : Qc(nullptr), Qm(nullptr), VQc(nullptr), VQm(&VQm), MQc(&MQc), MQm(nullptr)
  {
  }
  CurlCurlMassIntegrator(mfem::MatrixCoefficient &MQc, mfem::MatrixCoefficient &MQm)
    : Qc(nullptr), Qm(nullptr), VQc(nullptr), VQm(nullptr), MQc(&MQc), MQm(&MQm)
  {
  }

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Q grad u, grad v) for H1 elements.
class DiffusionIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  DiffusionIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  DiffusionIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  DiffusionIntegrator(mfem::VectorCoefficient &VQ) : Q(nullptr), VQ(&VQ), MQ(nullptr) {}
  DiffusionIntegrator(mfem::MatrixCoefficient &MQ) : Q(nullptr), VQ(nullptr), MQ(&MQ) {}

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Qd grad u, grad v) + (Qm u, v) for H1 elements.
class DiffusionMassIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Qd, *Qm;
  mfem::VectorCoefficient *VQd;
  mfem::MatrixCoefficient *MQd;

public:
  DiffusionMassIntegrator() : Qd(nullptr), Qm(nullptr), VQd(nullptr), MQd(nullptr) {}
  DiffusionMassIntegrator(mfem::Coefficient &Qd, mfem::Coefficient &Qm)
    : Qd(&Qd), Qm(&Qm), VQd(nullptr), MQd(nullptr)
  {
  }
  DiffusionMassIntegrator(mfem::VectorCoefficient &VQd, mfem::Coefficient &Qm)
    : Qd(nullptr), Qm(&Qm), VQd(&VQd), MQd(nullptr)
  {
  }
  DiffusionMassIntegrator(mfem::MatrixCoefficient &MQd, mfem::Coefficient &Qm)
    : Qd(nullptr), Qm(&Qm), VQd(nullptr), MQd(&MQd)
  {
  }

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Q div u, div v) for Raviart-Thomas elements.
class DivDivIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;

public:
  DivDivIntegrator() : Q(nullptr) {}
  DivDivIntegrator(mfem::Coefficient &Q) : Q(&Q) {}

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Qd div u, div v) + (Qm u, v) for Raviart-Thomas elements.
class DivDivMassIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Qd, *Qm;
  mfem::VectorCoefficient *VQm;
  mfem::MatrixCoefficient *MQm;

public:
  DivDivMassIntegrator() : Qd(nullptr), Qm(nullptr), VQm(nullptr), MQm(nullptr) {}
  DivDivMassIntegrator(mfem::Coefficient &Qd, mfem::Coefficient &Qm)
    : Qd(&Qd), Qm(&Qm), VQm(nullptr), MQm(nullptr)
  {
  }
  DivDivMassIntegrator(mfem::Coefficient &Qd, mfem::VectorCoefficient &VQm)
    : Qd(&Qd), Qm(nullptr), VQm(&VQm), MQm(nullptr)
  {
  }
  DivDivMassIntegrator(mfem::Coefficient &Qd, mfem::MatrixCoefficient &MQm)
    : Qd(&Qd), Qm(nullptr), VQm(nullptr), MQm(&MQm)
  {
  }

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Q grad u, v) for u in H1 and v in H(curl).
class MixedVectorGradientIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  MixedVectorGradientIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  MixedVectorGradientIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  MixedVectorGradientIntegrator(mfem::VectorCoefficient &VQ)
    : Q(nullptr), VQ(&VQ), MQ(nullptr)
  {
  }
  MixedVectorGradientIntegrator(mfem::MatrixCoefficient &MQ)
    : Q(nullptr), VQ(nullptr), MQ(&MQ)
  {
  }

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = -(Q u, grad v) for u in H(curl) and v in H1.
class MixedVectorWeakDivergenceIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  MixedVectorWeakDivergenceIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  MixedVectorWeakDivergenceIntegrator(mfem::Coefficient &Q)
    : Q(&Q), VQ(nullptr), MQ(nullptr)
  {
  }
  MixedVectorWeakDivergenceIntegrator(mfem::VectorCoefficient &VQ)
    : Q(nullptr), VQ(&VQ), MQ(nullptr)
  {
  }
  MixedVectorWeakDivergenceIntegrator(mfem::MatrixCoefficient &MQ)
    : Q(nullptr), VQ(nullptr), MQ(&MQ)
  {
  }

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Q curl u, v) for u in H(curl) and v in H(div).
class MixedVectorCurlIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  MixedVectorCurlIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  MixedVectorCurlIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  MixedVectorCurlIntegrator(mfem::VectorCoefficient &VQ) : Q(nullptr), VQ(&VQ), MQ(nullptr)
  {
  }
  MixedVectorCurlIntegrator(mfem::MatrixCoefficient &MQ) : Q(nullptr), VQ(nullptr), MQ(&MQ)
  {
  }

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Integrator for a(u, v) = (Q u, curl v) for u in H(div) and v in H(curl).
class MixedVectorWeakCurlIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  MixedVectorWeakCurlIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  MixedVectorWeakCurlIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  MixedVectorWeakCurlIntegrator(mfem::VectorCoefficient &VQ)
    : Q(nullptr), VQ(&VQ), MQ(nullptr)
  {
  }
  MixedVectorWeakCurlIntegrator(mfem::MatrixCoefficient &MQ)
    : Q(nullptr), VQ(nullptr), MQ(&MQ)
  {
  }

  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override;
};

// Base class for all discrete interpolators.
class DiscreteInterpolator : public BilinearFormIntegrator
{
public:
  void Assemble(const mfem::FiniteElementSpace &trial_fespace,
                const mfem::FiniteElementSpace &test_fespace,
                const mfem::IntegrationRule &ir, const std::vector<int> &indices, Ceed ceed,
                CeedOperator *op, CeedOperator *op_t) override;

  void AssembleBoundary(const mfem::FiniteElementSpace &trial_fespace,
                        const mfem::FiniteElementSpace &test_fespace,
                        const mfem::IntegrationRule &ir, const std::vector<int> &indices,
                        Ceed ceed, CeedOperator *op, CeedOperator *op_t) override
  {
    MFEM_ABORT("Boundary assembly is not implemented for DiscreteInterpolator objects!");
  }
};

// Interpolator for the identity map, where the domain space is a subspace of the range
// space (discrete embedding matrix).
using IdentityInterpolator = DiscreteInterpolator;

// Interpolator for the discrete gradient map from an H1 space to an H(curl) space.
using GradientInterpolator = DiscreteInterpolator;

// Interpolator for the discrete curl map from an H(curl) space to an H(div) space.
using CurlInterpolator = DiscreteInterpolator;

// Similar to MFEM's VectorFEBoundaryTangentLFIntegrator for ND spaces, but instead of
// computing (n x f, v), this just computes (f, v). Also eliminates the a and b quadrature
// parameters and uses fem::GetDefaultIntegrationOrder instead.
class VectorFEBoundaryLFIntegrator : public mfem::LinearFormIntegrator
{
private:
  mfem::VectorCoefficient &Q;
  mfem::DenseMatrix vshape;
  mfem::Vector f_loc, f_hat;
  int q_order;

public:
  VectorFEBoundaryLFIntegrator(mfem::VectorCoefficient &QG, int q_order = -1)
    : Q(QG), f_loc(QG.GetVDim()), q_order(q_order)
  {
  }

  void AssembleRHSElementVect(const mfem::FiniteElement &fe, mfem::ElementTransformation &T,
                              mfem::Vector &elvect) override;
};

// Similar to MFEM's BoundaryLFIntegrator for H1 spaces, but eliminates the a and b
// quadrature parameters and uses fem::GetDefaultIntegrationOrder instead.
class BoundaryLFIntegrator : public mfem::LinearFormIntegrator
{
private:
  mfem::Coefficient &Q;
  mfem::Vector shape;
  int q_order;

public:
  BoundaryLFIntegrator(mfem::Coefficient &QG, int q_order = -1) : Q(QG), q_order(q_order) {}

  void AssembleRHSElementVect(const mfem::FiniteElement &fe, mfem::ElementTransformation &T,
                              mfem::Vector &elvect) override;
};

using VectorFEDomainLFIntegrator = VectorFEBoundaryLFIntegrator;
using DomainLFIntegrator = BoundaryLFIntegrator;

}  // namespace palace

#endif  // PALACE_FEM_INTEGRATOR_HPP
