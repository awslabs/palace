// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "materialoperator.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include "fem/mesh.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"
#include "utils/omp.hpp"

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

MaterialOperator::MaterialOperator(const IoData &iodata, Mesh &mesh)
  : loc_attr(mesh.GetAttributeGlobalToLocal()),
    loc_bdr_attr(mesh.GetBdrAttributeGlobalToLocal()),
    local_to_shared(mesh.GetLocalToSharedFaceMap())
{
  SetUpMaterialProperties(iodata, mesh);
}

void MaterialOperator::SetUpMaterialProperties(const IoData &iodata, mfem::ParMesh &mesh)
{
  // Check that material attributes have been specified correctly. The mesh attributes may
  // be non-contiguous and when no material attribute is specified the elements are deleted
  // from the mesh so as to not cause problems.
  MFEM_VERIFY(!iodata.domains.materials.empty(), "Materials must be non-empty!");
  {
    int attr_max = mesh.attributes.Size() ? mesh.attributes.Max() : 0;
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
  }

  // Set up material properties of the different domain regions, represented with element-
  // wise constant matrix-valued coefficients for the relative permeability, permittivity,
  // and other material properties.
  mfem::Array<int> mat_marker(iodata.domains.materials.size());
  mat_marker = 0;
  int nmats = 0;
  for (std::size_t i = 0; i < iodata.domains.materials.size(); i++)
  {
    const auto &data = iodata.domains.materials[i];
    for (auto attr : data.attributes)
    {
      auto it = loc_attr.find(attr);
      if (it != loc_attr.end())
      {
        mat_marker[i] = 1;
        nmats++;
        break;
      }
    }
  }

  attr_mat.SetSize(loc_attr.size());  // XX TODO WIP...
  // int attr_max = 0;
  // for (auto &[attr, l_attr] : loc_attr)
  // {
  //   attr_max = std::max(attr_max, attr);
  // }
  // attr_mat.SetSize(attr_max);
  attr_mat = -1;

  const int sdim = mesh.SpaceDimension();
  mat_muinv.SetSize(sdim, sdim, nmats);
  mat_epsilon.SetSize(sdim, sdim, nmats);
  mat_epsilon_imag.SetSize(sdim, sdim, nmats);
  mat_epsilon_abs.SetSize(sdim, sdim, nmats);
  mat_invz0.SetSize(sdim, sdim, nmats);
  mat_c0.SetSize(sdim, sdim, nmats);
  mat_sigma.SetSize(sdim, sdim, nmats);
  mat_invLondon.SetSize(sdim, sdim, nmats);
  mat_c0_min.SetSize(nmats);
  mat_c0_max.SetSize(nmats);

  int count = 0;
  for (std::size_t i = 0; i < iodata.domains.materials.size(); i++)
  {
    if (!mat_marker[i])
    {
      continue;
    }
    const auto &data = iodata.domains.materials[i];
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

    // Map all attributes to this material property index.
    for (auto attr : data.attributes)
    {
      auto it = loc_attr.find(attr);
      if (it != loc_attr.end())
      {
        MFEM_VERIFY(
            attr_mat[it->second - 1] < 0,
            "Detected multiple definitions of material properties for domain attribute "
                << attr << "!");
        attr_mat[it->second - 1] = count;
      }
    }

    // Compute the inverse of the input permeability matrix.
    mfem::DenseMatrix mu_r = ToDenseMatrix(data.mu_r);
    mfem::DenseMatrixInverse(mu_r, true).GetInverseMatrix(mat_muinv(count));

    // Material permittivity: Re{ε} = ε, Im{ε} = -ε * tan(δ)
    mfem::DenseMatrix T(sdim, sdim);
    mat_epsilon(count) = ToDenseMatrix(data.epsilon_r);
    Mult(mat_epsilon(count), ToDenseMatrix(data.tandelta), T);
    T *= -1.0;
    mat_epsilon_imag(count) = T;
    if (mat_epsilon_imag(count).MaxMaxNorm() > 0.0)
    {
      for (auto attr : data.attributes)
      {
        losstan_attr.Append(attr);
      }
    }

    // ε * √(I + tan(δ) * tan(δ)ᵀ)
    MultAAt(ToDenseMatrix(data.tandelta), T);
    for (int d = 0; d < T.Height(); d++)
    {
      T(d, d) += 1.0;
    }
    Mult(mat_epsilon(count), MatrixSqrt(T), mat_epsilon_abs(count));

    // √μ⁻¹ ε
    Mult(mat_muinv(count), mat_epsilon(count), mat_invz0(count));
    mat_invz0(count) = MatrixSqrt(mat_invz0(count));

    // (√μ ε)⁻¹
    mfem::DenseMatrixInverse(mat_epsilon(count), true).GetInverseMatrix(T);
    Mult(mat_muinv(count), T, mat_c0(count));
    mat_c0(count) = MatrixSqrt(mat_c0(count));
    mat_c0_min[count] = mat_c0(count).CalcSingularvalue(sdim - 1);
    mat_c0_max[count] = mat_c0(count).CalcSingularvalue(0);

    // Electrical conductivity, σ
    mat_sigma(count) = ToDenseMatrix(data.sigma);
    if (mat_sigma(count).MaxMaxNorm() > 0.0)
    {
      for (auto attr : data.attributes)
      {
        conductivity_attr.Append(attr);
      }
    }

    // λ⁻² * μ⁻¹
    mat_invLondon(count) = mat_muinv(count);
    mat_invLondon(count) *=
        std::abs(data.lambda_L) > 0.0 ? std::pow(data.lambda_L, -2.0) : 0.0;
    if (mat_invLondon(count).MaxMaxNorm() > 0.0)
    {
      for (auto attr : data.attributes)
      {
        london_attr.Append(attr);
      }
    }

    count++;
  }

  // Set up map from boundary attributes (local) to material indices based on the
  // neighboring elements.

  bdr_attr_mat.SetSize(loc_bdr_attr.size());  // XX TODO WIP...
  // int bdr_attr_max = 0;
  // for (auto &[attr, l_attr] : loc_bdr_attr)
  // {
  //   bdr_attr_max = std::max(bdr_attr_max, attr);
  // }
  // bdr_attr_mat.SetSize(bdr_attr_max);
  bdr_attr_mat = -1;

  for (int i = 0; i < mesh.GetNBE(); i++)
  {
    const int attr = mesh.GetBdrAttribute(i);
    if (bdr_attr_mat[loc_bdr_attr.at(attr) - 1] >= 0)
    {
      continue;
    }

    // See also BdrGridFunctionCoefficient::GetElementTransformations.
    int f, o;
    int iel1, iel2, info1, info2;
    mesh.GetBdrElementFace(i, &f, &o);
    mesh.GetFaceElements(f, &iel1, &iel2);
    mesh.GetFaceInfos(f, &info1, &info2);

    mfem::FaceElementTransformations *FET;
    if (info2 >= 0 && iel2 < 0)
    {
      // Face is shared with another subdomain.
      const int &ishared = local_to_shared.at(f);
      FET = mesh.GetSharedFaceTransformations(ishared);
    }
    else
    {
      // Face is either internal to the subdomain, or a true one-sided boundary.
      FET = mesh.GetFaceElementTransformations(f);
    }

    auto *T1 = &FET->GetElement1Transformation();
    auto *T2 = (info2 >= 0) ? &FET->GetElement2Transformation() : nullptr;

    // For internal boundaries, use the element which corresponds to the vacuum domain, or
    // at least the one with the higher speed of light.
    const int nbr_attr =
        (T2 && GetLightSpeedMin(T2->Attribute) > GetLightSpeedMax(T1->Attribute))
            ? T2->Attribute
            : T1->Attribute;
    bdr_attr_mat[loc_bdr_attr.at(attr) - 1] = attr_mat[loc_attr.at(nbr_attr) - 1];
  }
}

MaterialPropertyCoefficient::MaterialPropertyCoefficient(
    const mfem::Array<int> &attr_mat_, const mfem::DenseTensor &mat_coeff_, double a)
  : attr_mat(attr_mat_), mat_coeff(mat_coeff_)
{
  for (int k = 0; k < mat_coeff.SizeK(); k++)
  {
    mat_coeff(k) *= a;
  }
}

namespace
{

void UpdateProperty(mfem::DenseTensor &mat_coeff, int k, double coeff, double a)
{
  // Constant diagonal coefficient.
  if (mat_coeff.SizeI() == 0 && mat_coeff.SizeJ() == 0)
  {
    // Initialize the coefficient material properties.
    MFEM_VERIFY(k == 0 && mat_coeff.SizeK() == 1,
                "Unexpected initial size for MaterialPropertyCoefficient!");
    mat_coeff.SetSize(1, 1, mat_coeff.SizeK());
    mat_coeff(0, 0, k) = a * coeff;
  }
  else
  {
    MFEM_VERIFY(mat_coeff.SizeI() == mat_coeff.SizeJ(),
                "Invalid dimensions for MaterialPropertyCoefficient update!");
    for (int i = 0; i < mat_coeff.SizeI(); i++)
    {
      mat_coeff(i, i, k) += a * coeff;
    }
  }
}

void UpdateProperty(mfem::DenseTensor &mat_coeff, int k, const mfem::DenseMatrix &coeff,
                    double a)
{
  if (mat_coeff.SizeI() == 0 && mat_coeff.SizeJ() == 0)
  {
    // Initialize the coefficient material properties.
    MFEM_VERIFY(k == 0 && mat_coeff.SizeK() == 1,
                "Unexpected initial size for MaterialPropertyCoefficient!");
    mat_coeff.SetSize(coeff.Height(), coeff.Width(), mat_coeff.SizeK());
    mat_coeff(k).Set(a, coeff);
  }
  else if (coeff.Height() == mat_coeff.SizeI() && coeff.Width() == mat_coeff.SizeJ())
  {
    // Add as full matrix.
    mat_coeff(k).Add(a, coeff);
  }
  else if (coeff.Height() == 1 && coeff.Width() == 1)
  {
    // Add as diagonal.
    UpdateProperty(mat_coeff, k, coeff(0, 0), a);
  }
  else if (mat_coeff.SizeI() == 1 && mat_coeff.SizeJ() == 1)
  {
    // Convert to matrix coefficient and previous data add as diagonal.
    mfem::DenseTensor mat_coeff_scalar(mat_coeff);
    mat_coeff.SetSize(coeff.Height(), coeff.Width(), mat_coeff_scalar.SizeK());
    mat_coeff = 0.0;
    for (int l = 0; l < mat_coeff.SizeK(); l++)
    {
      UpdateProperty(mat_coeff, l, mat_coeff_scalar(0, 0, l), 1.0);
    }
    mat_coeff(k).Add(a, coeff);
  }
  else
  {
    MFEM_ABORT("Invalid dimensions when updating material property at index " << k << "!");
  }
}

bool Equals(const mfem::DenseMatrix &mat_coeff, double coeff, double a)
{
  MFEM_VERIFY(mat_coeff.Height() == mat_coeff.Width(),
              "Invalid dimensions for MaterialPropertyCoefficient update!");
  constexpr double tol = 1.0e-9;
  for (int i = 0; i < mat_coeff.Height(); i++)
  {
    if (std::abs(mat_coeff(i, i) - a * coeff) >= tol * std::abs(mat_coeff(i, i)))
    {
      return false;
    }
    for (int j = 0; j < mat_coeff.Width(); j++)
    {
      if (j != i && std::abs(mat_coeff(i, j)) > 0.0)
      {
        return false;
      }
    }
  }
  return true;
}

bool Equals(const mfem::DenseMatrix &mat_coeff, const mfem::DenseMatrix &coeff, double a)
{
  if (coeff.Height() == 1 && coeff.Width() == 1)
  {
    return Equals(mat_coeff, coeff(0, 0), a);
  }
  else
  {
    constexpr double tol = 1.0e-9;
    mfem::DenseMatrix T(mat_coeff);
    T.Add(-a, coeff);
    return (T.MaxMaxNorm() < tol * mat_coeff.MaxMaxNorm());
  }
}

}  // namespace

void MaterialPropertyCoefficient::AddCoefficient(const mfem::Array<int> &attr_mat_,
                                                 const mfem::DenseTensor &mat_coeff_,
                                                 double a)
{
  if (empty())
  {
    attr_mat = attr_mat_;
    mat_coeff = mat_coeff_;
    for (int k = 0; k < mat_coeff.SizeK(); k++)
    {
      mat_coeff(k) *= a;
    }
  }
  else if (attr_mat_ == attr_mat)
  {
    MFEM_VERIFY(mat_coeff_.SizeK() == mat_coeff.SizeK(),
                "Invalid dimensions for MaterialPropertyCoefficient::AddCoefficient!");
    for (int k = 0; k < mat_coeff.SizeK(); k++)
    {
      UpdateProperty(mat_coeff, k, mat_coeff_(k), a);
    }
  }
  else
  {
    for (int k = 0; k < mat_coeff_.SizeK(); k++)
    {
      // Get list of all attributes which use this material property.
      mfem::Array<int> attr_list;
      attr_list.Reserve(attr_mat_.Size());
      for (int i = 0; i < attr_mat_.Size(); i++)
      {
        if (attr_mat_[i] == k)
        {
          attr_list.Append(i + 1);
        }
      }

      // Add or update the material property.
      AddMaterialProperty(attr_list, mat_coeff_(k), a);
    }
  }
}

template <typename T>
void MaterialPropertyCoefficient::AddMaterialProperty(const mfem::Array<int> &attr_list,
                                                      const T &coeff, double a)
{
  // Preprocess the attribute list. If any of the given attributes already have material
  // properties assigned, then they all need to point to the same material and it is
  // updated in place. Otherwise a new material is added for these attributes.
  int mat_idx = -1, attr_max = attr_mat.Size();
  for (auto attr : attr_list)
  {
    if (mat_idx < 0)
    {
      mat_idx = (attr > attr_mat.Size()) ? -1 : attr_mat[attr - 1];
    }
    else
    {
      MFEM_VERIFY(attr <= attr_mat.Size() && mat_idx == attr_mat[attr - 1],
                  "All attributes for AddMaterialProperty must correspond to the same "
                  "existing material if it exists!");
    }
    attr_max = std::max(attr, attr_max);
  }

  if (mat_idx < 0)
  {
    // Check if we can reuse an existing material.
    for (int k = 0; k < mat_coeff.SizeK(); k++)
    {
      if (Equals(mat_coeff(k), coeff, a))
      {
        mat_idx = k;
        break;
      }
    }
    if (mat_idx < 0)
    {
      // Append a new material and assign the attributes to it.
      const mfem::DenseTensor mat_coeff_backup(mat_coeff);
      mat_coeff.SetSize(mat_coeff_backup.SizeI(), mat_coeff_backup.SizeJ(),
                        mat_coeff_backup.SizeK() + 1);
      for (int k = 0; k < mat_coeff_backup.SizeK(); k++)
      {
        mat_coeff(k) = mat_coeff_backup(k);
      }
      mat_idx = mat_coeff.SizeK() - 1;
    }
    mat_coeff(mat_idx) = 0.0;  // Zero out so we can add

    // Copy the previous attribute materials, initialize no material to all new ones, then
    // populate.
    attr_mat.SetSize(attr_max, -1);
    for (auto attr : attr_list)
    {
      attr_mat[attr - 1] = mat_idx;
    }
  }
  UpdateProperty(mat_coeff, mat_idx, coeff, a);
}

void MaterialPropertyCoefficient::RestrictCoefficient(const mfem::Array<int> &attr_list)
{
  // Create a new material property coefficient with materials corresponding to only the
  // unique ones in the given attribute list.
  const mfem::Array<int> attr_mat_orig(attr_mat);
  const mfem::DenseTensor mat_coeff_orig(mat_coeff);
  attr_mat = -1;
  mat_coeff.SetSize(mat_coeff_orig.SizeI(), mat_coeff_orig.SizeJ(), 0);
  for (auto attr : attr_list)
  {
    if (attr_mat[attr - 1] >= 0)
    {
      // Attribute has already been processed.
      continue;
    }

    // Find all attributes in restricted list of attributes which map to this material index
    // and process them together.
    const int orig_mat_idx = attr_mat_orig[attr - 1];
    const int new_mat_idx = mat_coeff.SizeK();
    for (auto attr2 : attr_list)
    {
      if (attr_mat_orig[attr2 - 1] == orig_mat_idx)
      {
        attr_mat[attr2 - 1] = new_mat_idx;
      }
    }

    // Append the new material property.
    const mfem::DenseTensor mat_coeff_backup(mat_coeff);
    mat_coeff.SetSize(mat_coeff_backup.SizeI(), mat_coeff_backup.SizeJ(),
                      mat_coeff_backup.SizeK() + 1);
    for (int k = 0; k < mat_coeff_backup.SizeK(); k++)
    {
      mat_coeff(k) = mat_coeff_backup(k);
    }
    mat_coeff(new_mat_idx) = mat_coeff_orig(orig_mat_idx);
  }
}

void MaterialPropertyCoefficient::NormalProjectedCoefficient(const mfem::Vector &normal)
{
  mfem::DenseTensor mat_coeff_backup(mat_coeff);
  mat_coeff.SetSize(1, 1, mat_coeff_backup.SizeK());
  for (int k = 0; k < mat_coeff_backup.SizeK(); k++)
  {
    mat_coeff(k) = mat_coeff_backup(k).InnerProduct(normal, normal);
  }
}

template void MaterialPropertyCoefficient::AddMaterialProperty(const mfem::Array<int> &,
                                                               const mfem::DenseMatrix &,
                                                               double);
template void MaterialPropertyCoefficient::AddMaterialProperty(const mfem::Array<int> &,
                                                               const double &, double);

}  // namespace palace
