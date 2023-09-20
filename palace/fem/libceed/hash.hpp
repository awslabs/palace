// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HASH_HPP
#define PALACE_LIBCEED_HASH_HPP

#include <functional>
#include <mfem.hpp>

namespace palace::ceed
{

// Base case for combining hashes.
inline void CeedHashCombine(std::size_t &seed) {}

// See for example https://onlinelibrary.wiley.com/doi/abs/10.1002/asi.10170, the source
// of https://www.boost.org/doc/libs/1_35_0/doc/html/boost/hash_combine_id241013.html.
template <typename T, typename... U>
inline void CeedHashCombine(std::size_t &seed, const T &v, const U &...args)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  (CeedHashCombine(seed, args), ...);
}

namespace internal
{

struct ElementKey
{
  mfem::Geometry::Type type;
  mfem::Array<int> vertices;
  ElementKey(const mfem::Element &el)
    : type(el.GetGeometryType()), vertices(el.GetNVertices())
  {
    for (int i = 0; i < el.GetNVertices(); i++)
    {
      vertices[i] = el.GetVertices()[i];
    }
  }
  bool operator==(const ElementKey &k) const
  {
    if (type != k.type || vertices.Size() != k.vertices.Size())
    {
      return false;
    }
    for (int i = 0; i < vertices.Size(); i++)
    {
      if (vertices[i] != k.vertices[i])
      {
        return false;
      }
    }
    return true;
  }
  bool operator!=(const ElementKey &k) const { return !(*this == k); }
};

struct FiniteElementKey
{
  mfem::Geometry::Type type;
  int order, P;
  int space, range_type, map_type, deriv_type, deriv_range_type, deriv_map_type;
  FiniteElementKey(const mfem::FiniteElement &fe)
    : type(fe.GetGeomType()), order(fe.GetOrder()), P(fe.GetDof()), space(fe.Space()),
      range_type(fe.GetRangeType()), map_type(fe.GetMapType()),
      deriv_type(fe.GetDerivType()), deriv_range_type(fe.GetDerivRangeType()),
      deriv_map_type(fe.GetDerivMapType())
  {
  }
  bool operator==(const FiniteElementKey &k) const
  {
    return (type == k.type && order == k.order && P == k.P && space == k.space &&
            range_type == k.range_type && map_type == k.map_type &&
            deriv_type == k.deriv_type && deriv_range_type == k.deriv_range_type &&
            deriv_map_type == k.deriv_map_type);
  }
  bool operator!=(const FiniteElementKey &k) const { return !(*this == k); }
};

struct IntegrationRuleKey
{
  int order, Q;
  IntegrationRuleKey(const mfem::IntegrationRule &ir)
    : order(ir.GetOrder()), Q(ir.GetNPoints())
  {
  }
  bool operator==(const IntegrationRuleKey &k) const
  {
    return (order == k.order && Q == k.Q);
  }
  bool operator!=(const IntegrationRuleKey &k) const { return !(*this == k); }
};

struct BasisKey
{
  Ceed ceed;
  FiniteElementKey fe;
  IntegrationRuleKey ir;
  int ncomp;
  BasisKey(Ceed ceed, const mfem::FiniteElementSpace &fespace,
           const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir)
    : ceed(ceed), fe(fe), ir(ir), ncomp(fespace.GetVDim())
  {
  }
  bool operator==(const BasisKey &k) const
  {
    return (ceed == k.ceed && fe == k.fe && ir == k.ir && ncomp == k.ncomp);
  }
};

struct BasisHash
{
  std::size_t operator()(const BasisKey &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.ceed, k.fe, k.ir, k.ncomp);
    return hash;
  }
};

struct InterpBasisKey
{
  Ceed ceed;
  FiniteElementKey trial_fe, test_fe;
  int ncomp;
  InterpBasisKey(Ceed ceed, const mfem::FiniteElementSpace &trial_fespace,
                 const mfem::FiniteElementSpace &test_fespace,
                 const mfem::FiniteElement &trial_fe, const mfem::FiniteElement &test_fe)
    : ceed(ceed), trial_fe(trial_fe), test_fe(test_fe), ncomp(trial_fespace.GetVDim())
  {
  }
  bool operator==(const InterpBasisKey &k) const
  {
    return (ceed == k.ceed && trial_fe == k.trial_fe && test_fe == k.test_fe &&
            ncomp == k.ncomp);
  }
};

struct InterpBasisHash
{
  std::size_t operator()(const InterpBasisKey &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.ceed, k.trial_fe, k.test_fe, k.ncomp);
    return hash;
  }
};

struct RestrKey
{
  Ceed ceed;
  const mfem::Mesh *mesh;
  ElementKey el;
  FiniteElementKey fe;
  mfem::Array<int> dofs;
  int ne, ncomp, lsize;
  mfem::Ordering::Type ordering;
  bool unique_interp_restr, unique_interp_range_restr;
  RestrKey(Ceed ceed, const mfem::FiniteElementSpace &fespace, const mfem::Element &el,
           const mfem::FiniteElement &fe, const mfem::Array<int> &dofs, std::size_t ne,
           bool unique_interp_restr, bool unique_interp_range_restr)
    : ceed(ceed), mesh(fespace.GetMesh()), el(el), fe(fe), dofs(dofs), ne(ne),
      ncomp(fespace.GetVDim()), lsize(fespace.GetVDim() * fespace.GetNDofs()),
      ordering(fespace.GetOrdering()), unique_interp_restr(unique_interp_restr),
      unique_interp_range_restr(unique_interp_range_restr)
  {
  }
  bool operator==(const RestrKey &k) const
  {
    if (ceed != k.ceed || mesh != k.mesh || el != k.el || fe != k.fe ||
        dofs.Size() != k.dofs.Size() || ne != k.ne || ncomp != k.ncomp ||
        lsize != k.lsize || ordering != k.ordering ||
        unique_interp_restr != k.unique_interp_restr ||
        unique_interp_range_restr != k.unique_interp_range_restr)
    {
      return false;
    }
    for (int i = 0; i < dofs.Size(); i++)
    {
      if (dofs[i] != k.dofs[i])
      {
        return false;
      }
    }
    return true;
  }
};

struct RestrHash
{
  std::size_t operator()(const RestrKey &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.ceed, k.mesh, k.el, k.fe, k.dofs.Size(), k.ne, k.ncomp, k.lsize,
                    k.ordering, k.unique_interp_restr, k.unique_interp_range_restr);
    return hash;
  }
};

}  // namespace internal

}  // namespace palace::ceed

namespace std
{

template <>
struct hash<palace::ceed::internal::ElementKey>
{
  std::size_t operator()(const palace::ceed::internal::ElementKey &k) const noexcept
  {
    std::size_t hash = 0;
    palace::ceed::CeedHashCombine(hash, k.type, k.vertices.Size());
    return hash;
  }
};

template <>
struct hash<palace::ceed::internal::FiniteElementKey>
{
  std::size_t operator()(const palace::ceed::internal::FiniteElementKey &k) const noexcept
  {
    std::size_t hash = 0;
    palace::ceed::CeedHashCombine(hash, k.type, k.order, k.P, k.space, k.range_type,
                                  k.map_type, k.deriv_type, k.deriv_range_type,
                                  k.deriv_map_type);
    return hash;
  }
};

template <>
struct hash<palace::ceed::internal::IntegrationRuleKey>
{
  std::size_t operator()(const palace::ceed::internal::IntegrationRuleKey &k) const noexcept
  {
    std::size_t hash = 0;
    palace::ceed::CeedHashCombine(hash, k.order, k.Q);
    return hash;
  }
};

}  // namespace std

#endif  // PALACE_LIBCEED_HASH_HPP
