// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "materialoperator.hpp"

#include <cmath>
#include <functional>
#include <limits>
#include "linalg/densematrix.hpp"
#include "utils/communication.hpp"
#include "utils/geodata.hpp"
#include "utils/iodata.hpp"

namespace palace
{

namespace
{

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

MaterialOperator::MaterialOperator(const IoData &iodata, const Mesh &mesh) : mesh(mesh)
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
  const auto &loc_attr = this->mesh.GetCeedAttributes();
  mfem::Array<int> mat_marker(iodata.domains.materials.size());
  mat_marker = 0;
  int nmats = 0;
  for (std::size_t i = 0; i < iodata.domains.materials.size(); i++)
  {
    const auto &data = iodata.domains.materials[i];
    for (auto attr : data.attributes)
    {
      if (loc_attr.find(attr) != loc_attr.end())
      {
        mat_marker[i] = 1;
        nmats++;
        break;
      }
    }
  }
  attr_mat.SetSize(loc_attr.size());
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
  mat_muinvkx.SetSize(sdim, sdim, nmats);
  mat_kxTmuinvkx.SetSize(sdim, sdim, nmats);
  mat_kx.SetSize(sdim, sdim, nmats);
  has_losstan_attr = has_conductivity_attr = has_london_attr = has_wave_attr = false;

  // Set up Floquet wave vector for periodic meshes with phase-delay constraints.
  SetUpFloquetWaveVector(iodata, mesh);

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
        Mpi::Warning(
            "Electrostatic problem type does not account for material permeability,\n"
            "electrical conductivity, or London depth!\n");
      }
    }
    else if (iodata.problem.type == config::ProblemData::Type::MAGNETOSTATIC)
    {
      MFEM_VERIFY(IsValid(data.mu_r), "Material has no valid permeability defined!");
      if (!IsIdentity(data.epsilon_r) || IsValid(data.tandelta) || IsValid(data.sigma) ||
          std::abs(data.lambda_L) > 0.0)
      {
        Mpi::Warning(
            "Magnetostatic problem type does not account for material permittivity,\n"
            "loss tangent, electrical conductivity, or London depth!\n");
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
    mfem::DenseMatrix mat_mu = ToDenseMatrix(data.mu_r);
    mfem::DenseMatrixInverse(mat_mu, true).GetInverseMatrix(mat_muinv(count));

    // Material permittivity: Re{ε} = ε, Im{ε} = -ε * tan(δ)
    mfem::DenseMatrix T(sdim, sdim);
    mat_epsilon(count) = ToDenseMatrix(data.epsilon_r);
    Mult(mat_epsilon(count), ToDenseMatrix(data.tandelta), T);
    T *= -1.0;
    mat_epsilon_imag(count) = T;
    if (mat_epsilon_imag(count).MaxMaxNorm() > 0.0)
    {
      has_losstan_attr = true;
    }

    // ε * √(I + tan(δ) * tan(δ)ᵀ)
    MultAAt(ToDenseMatrix(data.tandelta), T);
    for (int d = 0; d < T.Height(); d++)
    {
      T(d, d) += 1.0;
    }
    Mult(mat_epsilon(count), linalg::MatrixSqrt(T), mat_epsilon_abs(count));

    // √(μ⁻¹ ε)
    Mult(mat_muinv(count), mat_epsilon(count), mat_invz0(count));
    mat_invz0(count) = linalg::MatrixSqrt(mat_invz0(count));

    // √((μ ε)⁻¹)
    Mult(mat_mu, mat_epsilon(count), T);
    mat_c0(count) = linalg::MatrixPow(T, -0.5);
    mat_c0_min[count] = linalg::SingularValueMin(mat_c0(count));
    mat_c0_max[count] = linalg::SingularValueMax(mat_c0(count));

    // Electrical conductivity, σ
    mat_sigma(count) = ToDenseMatrix(data.sigma);
    if (mat_sigma(count).MaxMaxNorm() > 0.0)
    {
      has_conductivity_attr = true;
    }

    // λ⁻² * μ⁻¹
    mat_invLondon(count) = mat_muinv(count);
    mat_invLondon(count) *=
        std::abs(data.lambda_L) > 0.0 ? std::pow(data.lambda_L, -2.0) : 0.0;
    if (mat_invLondon(count).MaxMaxNorm() > 0.0)
    {
      has_london_attr = true;
    }

    // μ⁻¹ [k x]
    Mult(mat_muinv(count), wave_vector_cross, mat_muinvkx(count));

    // [k x]^T μ⁻¹ [k x]
    T.Transpose(wave_vector_cross);
    Mult(T, mat_muinvkx(count), mat_kxTmuinvkx(count));

    // [k x]
    mat_kx(count) = wave_vector_cross;

    count++;
  }
  bool has_attr[4] = {has_losstan_attr, has_conductivity_attr, has_london_attr,
                      has_wave_attr};
  Mpi::GlobalOr(4, has_attr, mesh.GetComm());
  has_losstan_attr = has_attr[0];
  has_conductivity_attr = has_attr[1];
  has_london_attr = has_attr[2];
  has_wave_attr = has_attr[3];
}

void MaterialOperator::SetUpFloquetWaveVector(const IoData &iodata,
                                              const mfem::ParMesh &mesh)
{
  const int sdim = mesh.SpaceDimension();
  const double tol = std::numeric_limits<double>::epsilon();

  // Sum Floquet wave vector over periodic boundary pairs.
  mfem::Vector wave_vector(sdim), local_wave_vector(sdim);
  wave_vector = 0.0;
  for (const auto &data : iodata.boundaries.periodic)
  {
    MFEM_VERIFY(data.wave_vector.size() == sdim,
                "Floquet wave vector size must equal the spatial dimension.");
    std::copy(data.wave_vector.begin(), data.wave_vector.end(),
              local_wave_vector.GetData());
    wave_vector += local_wave_vector;
  }

  // Get Floquet wave vector specified outside of periodic boundary definitions.
  const auto &data = iodata.boundaries.floquet;
  MFEM_VERIFY(data.wave_vector.size() == sdim,
              "Floquet wave vector size must equal the spatial dimension.");
  std::copy(data.wave_vector.begin(), data.wave_vector.end(), local_wave_vector.GetData());
  wave_vector += local_wave_vector;
  has_wave_attr = (wave_vector.Norml2() > tol);

  MFEM_VERIFY(!has_wave_attr || iodata.problem.type == config::ProblemData::Type::DRIVEN ||
                  iodata.problem.type == config::ProblemData::Type::EIGENMODE,
              "Quasi-periodic Floquet boundary conditions are only available for "
              " frequency domain driven or eigenmode simulations!");
  MFEM_VERIFY(!has_wave_attr || sdim == 3,
              "Quasi-periodic Floquet periodic boundary conditions are only available "
              " in 3D!");

  // Get mesh dimensions in x/y/z coordinates
  mfem::Vector bbmin, bbmax;
  mesh::GetAxisAlignedBoundingBox(mesh, bbmin, bbmax);
  bbmax -= bbmin;

  // Ensure Floquet wave vector components are in range [-π/L, π/L].
  for (int i = 0; i < sdim; i++)
  {
    if (wave_vector[i] > M_PI / bbmax[i])
    {
      wave_vector[i] =
          -M_PI / bbmax[i] + fmod(wave_vector[i] + M_PI / bbmax[i], 2 * M_PI / bbmax[i]);
    }
    else if (wave_vector[i] < M_PI / bbmax[i])
    {
      wave_vector[i] =
          M_PI / bbmax[i] + fmod(wave_vector[i] - M_PI / bbmax[i], 2 * M_PI / bbmax[i]);
    }
  }

  // Matrix representation of cross product with wave vector
  // [k x] = | 0  -k3  k2|
  //         | k3  0  -k1|
  //         |-k2  k1  0 |
  wave_vector_cross.SetSize(3);
  wave_vector_cross = 0.0;
  wave_vector_cross(0, 1) = -wave_vector[2];
  wave_vector_cross(0, 2) = wave_vector[1];
  wave_vector_cross(1, 0) = wave_vector[2];
  wave_vector_cross(1, 2) = -wave_vector[0];
  wave_vector_cross(2, 0) = -wave_vector[1];
  wave_vector_cross(2, 1) = wave_vector[0];
}

mfem::Array<int> MaterialOperator::GetBdrAttributeToMaterial() const
{
  // Construct map from all (contiguous) local libCEED boundary attributes to the material
  // index in the neighboring element.
  mfem::Array<int> bdr_attr_mat(mesh.MaxCeedBdrAttribute());
  bdr_attr_mat = -1;
  for (const auto &[attr, bdr_attr_map] : mesh.GetCeedBdrAttributes())
  {
    for (auto it = bdr_attr_map.begin(); it != bdr_attr_map.end(); ++it)
    {
      MFEM_ASSERT(it->second > 0 && it->second <= bdr_attr_mat.Size(),
                  "Invalid libCEED boundary attribute " << it->second << "!");
      bdr_attr_mat[it->second - 1] = AttrToMat(it->first);
    }
  }
  return bdr_attr_mat;
}

MaterialPropertyCoefficient::MaterialPropertyCoefficient(int attr_max)
{
  attr_mat.SetSize(attr_max);
  attr_mat = -1;
}

MaterialPropertyCoefficient::MaterialPropertyCoefficient(
    const mfem::Array<int> &attr_mat_, const mfem::DenseTensor &mat_coeff_, double a)
  : attr_mat(attr_mat_), mat_coeff(mat_coeff_)
{
  *this *= a;
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
    MFEM_VERIFY(attr_mat_.Size() == attr_mat.Size(),
                "Invalid resize of attribute to material property map in "
                "MaterialPropertyCoefficient::AddCoefficient!");
    attr_mat = attr_mat_;
    mat_coeff = mat_coeff_;
    *this *= a;
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
  if (attr_list.Size() == 0)
  {
    // No attributes, nothing to add.
    return;
  }

  int mat_idx = -1;
  for (auto attr : attr_list)
  {
    MFEM_VERIFY(attr <= attr_mat.Size(),
                "Out of bounds access for attribute "
                    << attr << " in MaterialPropertyCoefficient::AddMaterialProperty!");
    if (mat_idx < 0)
    {
      mat_idx = attr_mat[attr - 1];
    }
    else
    {
      MFEM_VERIFY(mat_idx == attr_mat[attr - 1],
                  "All attributes for MaterialPropertyCoefficient::AddMaterialProperty "
                  "must correspond to the same "
                  "existing material if it exists!");
    }
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

    // Assign all attributes to this new material.
    for (auto attr : attr_list)
    {
      attr_mat[attr - 1] = mat_idx;
    }
  }
  UpdateProperty(mat_coeff, mat_idx, coeff, a);
}

MaterialPropertyCoefficient &MaterialPropertyCoefficient::operator*=(double a)
{
  for (int k = 0; k < mat_coeff.SizeK(); k++)
  {
    mat_coeff(k) *= a;
  }
  return *this;
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
  for (int k = 0; k < mat_coeff.SizeK(); k++)
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
