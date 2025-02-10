// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_FEM_INTEGRATOR_HPP
#define PALACE_FEM_INTEGRATOR_HPP

#include <mfem.hpp>
#include "fem/libceed/ceed.hpp"

namespace palace
{

class MaterialPropertyCoefficient;

//
// Classes which implement or extend bilinear and linear form integrators. In doc strings u
// refers to the trial function, and v the test function.
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

// Base class for libCEED-based bilinear form integrators.
class BilinearFormIntegrator
{
protected:
  const MaterialPropertyCoefficient *Q;
  bool assemble_q_data;
  bool transpose;

public:
  BilinearFormIntegrator(const MaterialPropertyCoefficient *Q = nullptr,
                         const bool transpose = false)
    : Q(Q), assemble_q_data(false), transpose(transpose)
  {
  }
  BilinearFormIntegrator(const MaterialPropertyCoefficient &Q, const bool transpose = false)
    : Q(&Q), assemble_q_data(false), transpose(transpose)
  {
  }
  virtual ~BilinearFormIntegrator() = default;

  virtual void Assemble(Ceed ceed, CeedElemRestriction trial_restr,
                        CeedElemRestriction test_restr, CeedBasis trial_basis,
                        CeedBasis test_basis, CeedVector geom_data,
                        CeedElemRestriction geom_data_restr, CeedOperator *op) const = 0;

  virtual void SetMapTypes(int trial_type, int test_type) {}

  void AssembleQuadratureData() { assemble_q_data = true; }
};

// Integrator for a(u, v) = (Q u, v) for H1 elements (also for vector (H1)ᵈ spaces).
class MassIntegrator : public BilinearFormIntegrator
{
public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
};

// Integrator for a(u, v) = (Q u, v) for vector finite elements.
class VectorFEMassIntegrator : public BilinearFormIntegrator
{
protected:
  int trial_map_type = mfem::FiniteElement::UNKNOWN_MAP_TYPE;
  int test_map_type = mfem::FiniteElement::UNKNOWN_MAP_TYPE;

public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;

  void SetMapTypes(int trial_type, int test_type) override
  {
    trial_map_type = trial_type;
    test_map_type = test_type;
  }
};

// Integrator for a(u, v) = (Q grad u, grad v) for H1 elements.
class DiffusionIntegrator : public BilinearFormIntegrator
{
public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
};

// Integrator for a(u, v) = (Q curl u, curl v) for Nedelec elements.
class CurlCurlIntegrator : public BilinearFormIntegrator
{
public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
};

// Integrator for a(u, v) = (Q div u, div v) for Raviart-Thomas elements.
class DivDivIntegrator : public BilinearFormIntegrator
{
public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
};

// Integrator for a(u, v) = (Qd grad u, grad v) + (Qm u, v) for H1 elements.
class DiffusionMassIntegrator : public BilinearFormIntegrator
{
protected:
  const MaterialPropertyCoefficient *Q_mass;
  bool transpose_mass;

public:
  using BilinearFormIntegrator::BilinearFormIntegrator;
  DiffusionMassIntegrator(const MaterialPropertyCoefficient &Q,
                          const MaterialPropertyCoefficient &Q_mass,
                          const bool transpose = false, const bool transpose_mass = false)
    : BilinearFormIntegrator(Q, transpose), Q_mass(&Q_mass), transpose_mass(transpose_mass)
  {
  }

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
};

// Integrator for a(u, v) = (Qc curl u, curl v) + (Qm u, v) for Nedelec elements.
class CurlCurlMassIntegrator : public BilinearFormIntegrator
{
protected:
  const MaterialPropertyCoefficient *Q_mass;
  bool transpose_mass;

public:
  using BilinearFormIntegrator::BilinearFormIntegrator;
  CurlCurlMassIntegrator(const MaterialPropertyCoefficient &Q,
                         const MaterialPropertyCoefficient &Q_mass,
                         const bool transpose = false, const bool transpose_mass = false)
    : BilinearFormIntegrator(Q, transpose), Q_mass(&Q_mass), transpose_mass(transpose_mass)
  {
  }

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
};

// Integrator for a(u, v) = (Qd div u, div v) + (Qm u, v) for Raviart-Thomas elements.
class DivDivMassIntegrator : public BilinearFormIntegrator
{
protected:
  const MaterialPropertyCoefficient *Q_mass;
  bool transpose_mass;

public:
  using BilinearFormIntegrator::BilinearFormIntegrator;
  DivDivMassIntegrator(const MaterialPropertyCoefficient &Q,
                       const MaterialPropertyCoefficient &Q_mass,
                       const bool transpose = false, const bool transpose_mass = false)
    : BilinearFormIntegrator(Q, transpose), Q_mass(&Q_mass), transpose_mass(transpose_mass)
  {
  }

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
};

// Integrator for a(u, v) = (Q grad u, v) for u in H1 and v in H(curl) or H(div).
class MixedVectorGradientIntegrator : public BilinearFormIntegrator
{
protected:
  int trial_map_type = mfem::FiniteElement::UNKNOWN_MAP_TYPE;
  int test_map_type = mfem::FiniteElement::UNKNOWN_MAP_TYPE;

public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;

  void SetMapTypes(int trial_type, int test_type) override
  {
    trial_map_type = trial_type;
    test_map_type = test_type;
  }
};

// Integrator for a(u, v) = -(Q u, grad v) for u in H(curl) and v in H1.
class MixedVectorWeakDivergenceIntegrator : public BilinearFormIntegrator
{
public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
};

// Integrator for a(u, v) = (Q curl u, v) for u in H(curl) and v in H(div) or H(curl).
class MixedVectorCurlIntegrator : public BilinearFormIntegrator
{
protected:
  int trial_map_type = mfem::FiniteElement::UNKNOWN_MAP_TYPE;
  int test_map_type = mfem::FiniteElement::UNKNOWN_MAP_TYPE;

public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;

  void SetMapTypes(int trial_type, int test_type) override
  {
    trial_map_type = trial_type;
    test_map_type = test_type;
  }
};

// Integrator for a(u, v) = -(Q u, curl v) for u in H(div) or H(curl) and v in H(curl).
class MixedVectorWeakCurlIntegrator : public BilinearFormIntegrator
{
protected:
  int trial_map_type = mfem::FiniteElement::UNKNOWN_MAP_TYPE;
  int test_map_type = mfem::FiniteElement::UNKNOWN_MAP_TYPE;

public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;

  void SetMapTypes(int trial_type, int test_type) override
  {
    trial_map_type = trial_type;
    test_map_type = test_type;
  }
};

// Integrator for a(u, v) = (Q grad u, v) for u in H1 and v in (H1)ᵈ.
class GradientIntegrator : public BilinearFormIntegrator
{
public:
  using BilinearFormIntegrator::BilinearFormIntegrator;

  void Assemble(Ceed ceed, CeedElemRestriction trial_restr, CeedElemRestriction test_restr,
                CeedBasis trial_basis, CeedBasis test_basis, CeedVector geom_data,
                CeedElemRestriction geom_data_restr, CeedOperator *op) const override;
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
