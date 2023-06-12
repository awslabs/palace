// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "materialoperator.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

// Compute matrix functions for symmetric real-valued 2x2 or 3x3 matrices. Returns the
// matrix U * f(Λ) * U' for input U * Λ * U'
// Reference: Deledalle et al., Closed-form expressions of the eigen decomposition of 2x2
//            and 3x3 Hermitian matrices, HAL hal-01501221 (2017).
mfem::DenseMatrix MatrixFunction(const mfem::DenseMatrix &M,
                                 const std::function<double(const double &)> &functor)
{
  MFEM_ASSERT(M.Height() == M.Width(),
              "MatrixFunction only available for square matrices!");
  const auto N = M.Height();
  constexpr auto tol = 10.0 * std::numeric_limits<double>::epsilon();
  for (int i = 0; i < N; i++)
  {
    for (int j = i + 1; j < N; j++)
    {
      MFEM_VERIFY(std::abs(M(i, j) - M(j, i)) < tol,
                  "MatrixFunction only available for symmetric matrices!");
    }
  }
  mfem::DenseMatrix Mout(N, N);
  Mout = 0.0;
  if (N == 2)
  {
    MFEM_ABORT("2x2 MatrixFunction is not implemented yet!");
  }
  else if (N == 3)
  {
    // Need to specialize based on the number of zeros and their locations.
    const auto &a = M(0, 0), &b = M(1, 1), &c = M(2, 2);
    const auto &d = M(0, 1), &e = M(1, 2), &f = M(0, 2);
    const bool d_non_zero = std::abs(d) > tol;
    const bool e_non_zero = std::abs(e) > tol;
    const bool f_non_zero = std::abs(f) > tol;
    if (!d_non_zero && !e_non_zero && !f_non_zero)
    {
      // a 0 0
      // 0 b 0
      // 0 0 c
      for (int i = 0; i < 3; i++)
      {
        Mout(i, i) = functor(M(i, i));
      }
      return Mout;
    }
    if (d_non_zero && !e_non_zero && !f_non_zero)
    {
      // a d 0
      // d b 0
      // 0 0 c
      const double disc = std::sqrt(a * a - 2.0 * a * b + b * b + 4.0 * d * d);
      const double lambda1 = c;
      const double lambda2 = (a + b - disc) / 2.0;
      const double lambda3 = (a + b + disc) / 2.0;
      const mfem::Vector v1{{0.0, 0.0, 1.0}};
      const mfem::Vector v2{{-(-a + b + disc) / (2.0 * d), 1.0, 0.0}};
      const mfem::Vector v3{{-(-a + b - disc) / (2.0 * d), 1.0, 0.0}};
      AddMult_a_VVt(functor(lambda1), v1, Mout);
      AddMult_a_VVt(functor(lambda2), v2, Mout);
      AddMult_a_VVt(functor(lambda3), v3, Mout);
      return Mout;
    }
    if (!d_non_zero && e_non_zero && !f_non_zero)
    {
      // a 0 0
      // 0 b e
      // 0 e c
      const double disc = std::sqrt(b * b - 2.0 * b * c + c * c + 4.0 * e * e);
      const double lambda1 = a;
      const double lambda2 = 0.5 * (b + c - disc);
      const double lambda3 = 0.5 * (b + c + disc);
      const mfem::Vector v1{{1.0, 0.0, 0.0}};
      const mfem::Vector v2{{0.0, -(-b + c + disc) / (2.0 * e), 1.0}};
      const mfem::Vector v3{{0.0, -(-b + c - disc) / (2.0 * e), 1.0}};
      AddMult_a_VVt(functor(lambda1), v1, Mout);
      AddMult_a_VVt(functor(lambda2), v2, Mout);
      AddMult_a_VVt(functor(lambda3), v3, Mout);
      return Mout;
    }
    if (!d_non_zero && !e_non_zero && f_non_zero)
    {
      // a 0 f
      // 0 b 0
      // f 0 c
      const double disc = std::sqrt(a * a - 2.0 * a * c + c * c + 4.0 * f * f);
      const double lambda1 = b;
      const double lambda2 = 0.5 * (a + c - disc);
      const double lambda3 = 0.5 * (a + c + disc);
      const mfem::Vector v1{{0.0, 1.0, 0.0}};
      const mfem::Vector v2{{-(-a + c + disc) / (2.0 * f), 0.0, 1.0}};
      const mfem::Vector v3{{-(-a + c - disc) / (2.0 * f), 0.0, 1.0}};
      AddMult_a_VVt(functor(lambda1), v1, Mout);
      AddMult_a_VVt(functor(lambda2), v2, Mout);
      AddMult_a_VVt(functor(lambda3), v3, Mout);
      return Mout;
    }
    if ((!d_non_zero && e_non_zero && f_non_zero) ||
        (d_non_zero && !e_non_zero && f_non_zero) ||
        (d_non_zero && e_non_zero && !f_non_zero))
    {
      MFEM_ABORT("This nonzero pattern is not currently supported for MatrixFunction!");
    }
    // General case for all nonzero:
    // a d f
    // d b e
    // f e c
    const double a2 = a * a, b2 = b * b, c2 = c * c, d2 = d * d, e2 = e * e, f2 = f * f;
    const double a2mbmc = 2.0 * a - b - c;
    const double b2mamc = 2.0 * b - a - c;
    const double c2mamb = 2.0 * c - a - b;
    const double x1 = a2 + b2 + c2 - a * b - b * c + 3.0 * (d2 + e2 + f2);
    const double x2 = -(a2mbmc * b2mamc * c2mamb) +
                      9.0 * (c2mamb * d2 + b2mamc * f2 + a2mbmc * e2) - 54.0 * d * e * f;
    const double phi = std::atan2(std::sqrt(4.0 * x1 * x1 * x1 - x2 * x2), x2);
    const double lambda1 = (a + b + c - 2.0 * std::sqrt(x1) * std::cos(phi / 3.0)) / 3.0;
    const double lambda2 =
        (a + b + c + 2.0 * std::sqrt(x1) * std::cos((phi - M_PI) / 3.0)) / 3.0;
    const double lambda3 =
        (a + b + c + 2.0 * std::sqrt(x1) * std::cos((phi + M_PI) / 3.0)) / 3.0;

    auto SafeDivide = [&](double x, double y)
    {
      if (std::abs(x) <= tol)
      {
        return 0.0;
      }
      if (std::abs(x) >= tol && std::abs(y) <= tol)
      {
        MFEM_ABORT("Logic error: Zero denominator with nonzero numerator!");
        return 0.0;
      }
      return x / y;
    };
    const double m1 = SafeDivide(d * (c - lambda1) - e * f, f * (b - lambda1) - d * e);
    const double m2 = SafeDivide(d * (c - lambda2) - e * f, f * (b - lambda2) - d * e);
    const double m3 = SafeDivide(d * (c - lambda3) - e * f, f * (b - lambda3) - d * e);
    const double l1mcmem1 = lambda1 - c - e * m1;
    const double l2mcmem2 = lambda2 - c - e * m2;
    const double l3mcmem3 = lambda3 - c - e * m3;
    const double n1 = 1.0 + m1 * m1 + SafeDivide(std::pow(l1mcmem1, 2), f2);
    const double n2 = 1.0 + m2 * m2 + SafeDivide(std::pow(l2mcmem2, 2), f2);
    const double n3 = 1.0 + m3 * m3 + SafeDivide(std::pow(l3mcmem3, 2), f2);
    const double tlambda1 = functor(lambda1) / n1;
    const double tlambda2 = functor(lambda2) / n2;
    const double tlambda3 = functor(lambda3) / n3;

    const double at = (tlambda1 * l1mcmem1 * l1mcmem1 + tlambda2 * l2mcmem2 * l2mcmem2 +
                       tlambda3 * l3mcmem3 * l3mcmem3) /
                      f2;
    const double bt = tlambda1 * m1 * m1 + tlambda2 * m2 * m2 + tlambda3 * m3 * m3;
    const double ct = tlambda1 + tlambda2 + tlambda3;
    const double dt =
        (tlambda1 * m1 * l1mcmem1 + tlambda2 * m2 * l2mcmem2 + tlambda3 * m3 * l3mcmem3) /
        f;
    const double et = tlambda1 * m1 + tlambda2 * m2 + tlambda3 * m3;
    const double ft = (tlambda1 * l1mcmem1 + tlambda2 * l2mcmem2 + tlambda3 * l3mcmem3) / f;
    Mout(0, 0) = at;
    Mout(0, 1) = dt;
    Mout(0, 2) = ft;
    Mout(1, 0) = dt;
    Mout(1, 1) = bt;
    Mout(1, 2) = et;
    Mout(2, 0) = ft;
    Mout(2, 1) = et;
    Mout(2, 2) = ct;
    return Mout;
  }
  else
  {
    MFEM_ABORT("MatrixFunction only supports 2x2 or 3x3 matrices!");
  }
  return Mout;
}

mfem::DenseMatrix MatrixSqrt(const mfem::DenseMatrix &M)
{
  return MatrixFunction(M, [](auto s) { return std::sqrt(s); });
}

template <std::size_t N>
bool IsValid(const config::SymmetricMatrixData<N> &data)
{
  // All the coefficients are nonzero.
  bool valid =
      std::all_of(data.s.begin(), data.s.end(), [](auto d) { return std::abs(d) > 0.0; });

  // All the vectors are normalized.
  constexpr auto tol = 1.0e-6;
  auto UnitNorm = [&](const std::array<double, N> &x)
  {
    double s = -1.0;
    for (const auto &i : x)
    {
      s += std::pow(i, 2);
    }
    return std::abs(s) < tol;
  };
  valid &= std::all_of(data.v.begin(), data.v.end(), UnitNorm);

  // All the vectors are orthogonal.
  for (std::size_t i1 = 0; i1 < N; i1++)
  {
    const auto &v1 = data.v.at(i1);
    for (std::size_t i2 = i1 + 1; i2 < N; i2++)
    {
      const auto &v2 = data.v.at(i2);
      double s = 0.0;
      for (std::size_t j = 0; j < N; j++)
      {
        s += v1[j] * v2[j];
      }
      valid &= std::abs(s) < tol;
    }
  }
  return valid;
}

template <std::size_t N>
bool IsIsotropic(const config::SymmetricMatrixData<N> &data)
{
  bool valid = true;
  for (std::size_t i = 0; i < N; i++)
  {
    for (std::size_t j = 0; j < N; j++)
    {
      if (i == j)
      {
        valid &= data.v[i][j] == 1.0;
      }
      else
      {
        valid &= data.v[i][j] == 0.0;
      }
    }
  }
  return valid;
}

template <std::size_t N>
bool IsIdentity(const config::SymmetricMatrixData<N> &data)
{
  auto valid = std::all_of(data.s.begin(), data.s.end(), [](auto d) { return d == 1.0; });
  valid &= IsIsotropic(data);
  return valid;
}

template <std::size_t N>
mfem::DenseMatrix ToDenseMatrix(const config::SymmetricMatrixData<N> &data)
{
  mfem::DenseMatrix M(N, N);
  mfem::Vector V(N);
  for (std::size_t i = 0; i < N; i++)
  {
    for (std::size_t j = 0; j < N; j++)
    {
      V(j) = data.v[i][j];
    }
    AddMult_a_VVt(data.s[i], V, M);
  }
  return M;
}

}  // namespace

MaterialOperator::MaterialOperator(const IoData &iodata, const mfem::ParMesh &mesh)
  : sdim(mesh.SpaceDimension())
{
  SetUpMaterialProperties(iodata, mesh);
}

void MaterialOperator::SetUpMaterialProperties(const IoData &iodata,
                                               const mfem::ParMesh &mesh)
{
  // Check that material attributes have been specified correctly. The mesh attributes may
  // be non-contiguous and when no material attribute is specified the elements are deleted
  // from the mesh so as to not cause problems.
  MFEM_VERIFY(!iodata.domains.materials.empty(), "Materials must be non-empty!");
  int attr_max = mesh.attributes.Max();
  mfem::Array<int> attr_marker(attr_max);
  attr_marker = 0;
  for (auto attr : mesh.attributes)
  {
    attr_marker[attr - 1] = 1;
  }
  for (const auto &data : iodata.domains.materials)
  {
    for (auto attr : data.attributes)
    {
      MFEM_VERIFY(
          attr > 0 && attr <= attr_max,
          "Material attribute tags must be non-negative and correspond to attributes "
          "in the mesh!");
      MFEM_VERIFY(attr_marker[attr - 1], "Unknown material attribute " << attr << "!");
    }
  }

  // Set up material properties of the different domain regions, represented with piece-wise
  // constant matrix-valued coefficients for the relative permeability and permittivity,
  // and other material properties.
  mat_muinv.resize(attr_max, mfem::DenseMatrix(sdim));
  mat_epsilon.resize(attr_max, mfem::DenseMatrix(sdim));
  mat_epsilon_imag.resize(attr_max, mfem::DenseMatrix(sdim));
  mat_epsilon_abs.resize(attr_max, mfem::DenseMatrix(sdim));
  mat_invz0.resize(attr_max, mfem::DenseMatrix(sdim));
  mat_c0.resize(attr_max, mfem::DenseMatrix(sdim));
  mat_sigma.resize(attr_max, mfem::DenseMatrix(sdim));
  mat_invLondon.resize(attr_max, mfem::DenseMatrix(sdim));
  mat_c0_min.resize(attr_max, 0.0);
  mat_c0_max.resize(attr_max, 0.0);
  for (const auto &data : iodata.domains.materials)
  {
    if (iodata.problem.type == config::ProblemData::Type::ELECTROSTATIC)
    {
      MFEM_VERIFY(IsValid(data.epsilon_r), "Material has no valid permittivity defined!");
      if (!IsIdentity(data.mu_r) || IsValid(data.sigma) || std::abs(data.lambda_L) > 0.0)
      {
        Mpi::Warning("Electrostatic problem type does not account for material "
                     "permeability, electrical conductivity, or London depth!\n");
      }
    }
    else if (iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC)
    {
      MFEM_VERIFY(IsValid(data.mu_r), "Material has no valid permeability defined!");
      if (!IsIdentity(data.epsilon_r) || IsValid(data.tandelta) || IsValid(data.sigma) ||
          std::abs(data.lambda_L) > 0.0)
      {
        Mpi::Warning(
            "Magnetostatic problem type does not account for material permittivity, loss "
            "tangent, electrical conductivity, or London depth!\n");
      }
    }
    else
    {
      MFEM_VERIFY(IsValid(data.mu_r) && IsValid(data.epsilon_r),
                  "Material has no valid permeability or no valid permittivity defined!");
      if (iodata.problem.type == config::ProblemData::Type::TRANSIENT)
      {
        MFEM_VERIFY(!IsValid(data.tandelta),
                    "Transient problem type does not support material loss tangent, use "
                    "electrical conductivity instead!");
      }
      else
      {
        MFEM_VERIFY(!(IsValid(data.tandelta) && IsValid(data.sigma)),
                    "Material loss model should probably use only one of loss tangent or "
                    "electrical conductivity!");
      }
    }
    for (auto attr : data.attributes)
    {
      MFEM_VERIFY(
          mat_c0_min.at(attr - 1) == 0.0 && mat_c0_max.at(attr - 1) == 0.0,
          "Detected multiple definitions of material properties for domain attribute "
              << attr << "!");

      // Compute the inverse of the input permeability matrix.
      mfem::DenseMatrix mu_r = ToDenseMatrix(data.mu_r);
      mfem::DenseMatrixInverse(mu_r, true).GetInverseMatrix(mat_muinv.at(attr - 1));

      // Material permittivity: Im{ϵ} = - ϵ * tan(δ)
      mfem::DenseMatrix T(sdim, sdim);
      mat_epsilon.at(attr - 1) = ToDenseMatrix(data.epsilon_r);
      Mult(mat_epsilon.at(attr - 1), ToDenseMatrix(data.tandelta), T);
      T *= -1.0;
      mat_epsilon_imag.at(attr - 1) = T;

      // ϵ * √(I + tan(δ) * tan(δ)ᵀ)
      MultAAt(ToDenseMatrix(data.tandelta), T);
      for (int i = 0; i < T.Height(); i++)
      {
        T(i, i) += 1.0;
      }
      Mult(mat_epsilon.at(attr - 1), MatrixSqrt(T), mat_epsilon_abs.at(attr - 1));

      // √μ⁻¹ ϵ
      Mult(mat_muinv.at(attr - 1), mat_epsilon.at(attr - 1), mat_invz0.at(attr - 1));
      mat_invz0.at(attr - 1) = MatrixSqrt(mat_invz0.at(attr - 1));

      // (√μ ϵ)⁻¹
      mfem::DenseMatrixInverse(mat_epsilon.at(attr - 1), true).GetInverseMatrix(T);
      Mult(mat_muinv.at(attr - 1), T, mat_c0.at(attr - 1));
      mat_c0.at(attr - 1) = MatrixSqrt(mat_c0.at(attr - 1));
      mat_c0_min.at(attr - 1) = mat_c0.at(attr - 1).CalcSingularvalue(sdim - 1);
      mat_c0_max.at(attr - 1) = mat_c0.at(attr - 1).CalcSingularvalue(0);

      // Electrical conductivity, σ
      mat_sigma.at(attr - 1) = ToDenseMatrix(data.sigma);

      // λ⁻² * μ⁻¹
      mat_invLondon.at(attr - 1) = mat_muinv.at(attr - 1);
      mat_invLondon.at(attr - 1) *=
          std::abs(data.lambda_L) > 0.0 ? std::pow(data.lambda_L, -2.0) : 0.0;
    }
  }

  // Mark selected material attributes from the mesh as having certain local properties.
  mfem::Array<int> losstan_mats, conductivity_mats, london_mats;
  losstan_mats.Reserve(attr_max);
  conductivity_mats.Reserve(attr_max);
  london_mats.Reserve(attr_max);
  for (int i = 0; i < attr_max; i++)
  {
    if (mat_epsilon_imag.at(i).MaxMaxNorm() > 0.0)
    {
      losstan_mats.Append(i + 1);  // Markers are 1-based
    }
    if (mat_sigma.at(i).MaxMaxNorm() > 0.0)
    {
      conductivity_mats.Append(i + 1);
    }
    if (mat_invLondon.at(i).MaxMaxNorm() > 0.0)
    {
      london_mats.Append(i + 1);
    }
  }
  mesh::AttrToMarker(attr_max, losstan_mats, losstan_marker);
  mesh::AttrToMarker(attr_max, conductivity_mats, conductivity_marker);
  mesh::AttrToMarker(attr_max, london_mats, london_marker);
}

}  // namespace palace
