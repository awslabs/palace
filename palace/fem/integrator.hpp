// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_INTEGRATOR_HPP
#define PALACE_FEM_INTEGRATOR_HPP

#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"

namespace palace
{

//
// Classes which implement or extend bilinear and linear form integrators.
//

namespace fem
{

// Helper functions for creating an integration rule to exactly integrate polynomials of
// order 2 * p_trial + order(|J|) + q_extra.
struct DefaultIntegrationOrder
{
  inline static int p_trial = 1;
  inline static bool q_order_jac = true;
  inline static int q_order_extra_pk = 0;
  inline static int q_order_extra_qk = 0;
  static int Get(const mfem::IsoparametricTransformation &T);
  static int Get(const mfem::ElementTransformation &T);
  static int Get(const mfem::Mesh &mesh, mfem::Geometry::Type geom);
};

}  // namespace fem

// XX TODO WIP FOR NOW, NO COEFFICIENTS IN ASSEMBLY
// XX TODO WIP INITIALIZE ALL MEMBERS IN CONSTRUCTORS (MAP TYPE VARIABLES)

// XX TODO CONST FOR ASSEMBLE METHOD?

// XX TODO SHARED COEFFICIENT ASSEMBLY METHOD?

// Base class for libCEED-based bilinear form integrators.
class BilinearFormIntegrator
{
public:
  virtual ~BilinearFormIntegrator() = default;

  virtual void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                        CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                        CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) = 0;

  virtual void SetMapTypes(int trial_type, int test_type) {}
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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
};

// Integrator for a(u, v) = (Q u, v) for vector finite elements.
class VectorFEMassIntegrator : public BilinearFormIntegrator
{
protected:
  int trial_map_type, test_map_type;

  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  VectorFEMassIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  VectorFEMassIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  VectorFEMassIntegrator(mfem::VectorCoefficient &VQ) : Q(nullptr), VQ(&VQ), MQ(nullptr) {}
  VectorFEMassIntegrator(mfem::MatrixCoefficient &MQ) : Q(nullptr), VQ(nullptr), MQ(&MQ) {}

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;

  void SetMapTypes(int trial_type, int test_type) override
  {
    trial_map_type = trial_type;
    test_map_type = test_type;
  }
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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
};

// Integrator for a(u, v) = (Qc curl u, curl v) + (Qm u, v) for Nedelec elements.
class CurlCurlMassIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Qc, *Qm;
  mfem::VectorCoefficient *VQc, *VQm;
  mfem::MatrixCoefficient *MQc, *MQm;

public:
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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
};

// Integrator for a(u, v) = (Qd grad u, grad v) + (Qm u, v) for H1 elements.
class DiffusionMassIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Qd, *Qm;
  mfem::VectorCoefficient *VQd;
  mfem::MatrixCoefficient *MQd;

public:
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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
};

// Integrator for a(u, v) = (Q div u, div v) for Raviart-Thomas elements.
class DivDivIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;

public:
  DivDivIntegrator() : Q(nullptr) {}
  DivDivIntegrator(mfem::Coefficient &Q) : Q(&Q) {}

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
};

// Integrator for a(u, v) = (Qd div u, div v) + (Qm u, v) for Raviart-Thomas elements.
class DivDivMassIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Qd, *Qm;
  mfem::VectorCoefficient *VQm;
  mfem::MatrixCoefficient *MQm;

public:
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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
};

// Integrator for a(u, v) = (Q curl u, v) for u in H(curl) and v in H(div).
class MixedVectorCurlIntegrator : public BilinearFormIntegrator
{
protected:
  int trial_map_type, test_map_type;

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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;

  void SetMapTypes(int trial_type, int test_type) override
  {
    trial_map_type = trial_type;
    test_map_type = test_type;
  }
};

// Integrator for a(u, v) = (Q u, curl v) for u in H(div) and v in H(curl).
class MixedVectorWeakCurlIntegrator : public BilinearFormIntegrator
{
protected:
  int trial_map_type, test_map_type;

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

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;

  void SetMapTypes(int trial_type, int test_type) override
  {
    trial_map_type = trial_type;
    test_map_type = test_type;
  }
};

// Integrator for a(u, v) = (Q grad u, v) for u in H1 and v in (H1)ᵈ.
class GradientIntegrator : public BilinearFormIntegrator
{
protected:
  mfem::Coefficient *Q;
  mfem::VectorCoefficient *VQ;
  mfem::MatrixCoefficient *MQ;

public:
  GradientIntegrator() : Q(nullptr), VQ(nullptr), MQ(nullptr) {}
  GradientIntegrator(mfem::Coefficient &Q) : Q(&Q), VQ(nullptr), MQ(nullptr) {}
  GradientIntegrator(mfem::VectorCoefficient &VQ) : Q(nullptr), VQ(&VQ), MQ(nullptr) {}
  GradientIntegrator(mfem::MatrixCoefficient &MQ) : Q(nullptr), VQ(nullptr), MQ(&MQ) {}

  void Assemble(const ceed::CeedGeomFactorData &geom_data, Ceed ceed,
                CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedOperator *op) override;
};

// Base class for all discrete interpolators.
class DiscreteInterpolator
{
public:
  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis interp_basis, CeedOperator *op, CeedOperator *op_t);
};

// Interpolator for the identity map, where the domain space is a subspace of the range
// space (discrete embedding matrix).
using IdentityInterpolator = DiscreteInterpolator;

// Interpolator for the discrete gradient map from an H1 space to an H(curl) space.
using GradientInterpolator = DiscreteInterpolator;

// Interpolator for the discrete curl map from an H(curl) space to an H(div) space.
using CurlInterpolator = DiscreteInterpolator;

// Interpolator for the discrete divergence map from an H(div) space to an L2 space.
using DivergenceInterpolator = DiscreteInterpolator;

// Similar to MFEM's VectorFEBoundaryTangentLFIntegrator for ND spaces, but instead of
// computing (n x f, v), this just computes (f, v). Also eliminates the a and b quadrature
// parameters and uses fem::DefaultIntegrationOrder instead.
class VectorFEBoundaryLFIntegrator : public mfem::LinearFormIntegrator
{
private:
  mfem::VectorCoefficient &Q;
  mfem::DenseMatrix vshape;
  mfem::Vector f_loc, f_hat;

public:
  VectorFEBoundaryLFIntegrator(mfem::VectorCoefficient &QG) : Q(QG) {}

  void AssembleRHSElementVect(const mfem::FiniteElement &fe, mfem::ElementTransformation &T,
                              mfem::Vector &elvect) override;
};

// Similar to MFEM's BoundaryLFIntegrator for H1 spaces, but eliminates the a and b
// quadrature parameters and uses fem::DefaultIntegrationOrder instead.
class BoundaryLFIntegrator : public mfem::LinearFormIntegrator
{
private:
  mfem::Coefficient &Q;
  mfem::Vector shape;

public:
  BoundaryLFIntegrator(mfem::Coefficient &QG) : Q(QG) {}

  void AssembleRHSElementVect(const mfem::FiniteElement &fe, mfem::ElementTransformation &T,
                              mfem::Vector &elvect) override;
};

using VectorFEDomainLFIntegrator = VectorFEBoundaryLFIntegrator;
using DomainLFIntegrator = BoundaryLFIntegrator;

}  // namespace palace

#endif  // PALACE_FEM_INTEGRATOR_HPP
