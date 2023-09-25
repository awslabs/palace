// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HASH_HPP
#define PALACE_LIBCEED_HASH_HPP

#include <functional>
#include <string>
#include <utility>
#include <mfem.hpp>
#include "fem/fespace.hpp"

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
};

using FiniteElementPairKey = std::pair<FiniteElementKey, FiniteElementKey>;

struct FiniteElementPairHash
{
  std::size_t operator()(const FiniteElementPairKey &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.first, k.second);
    return hash;
  }
};

struct BasisKey
{
  Ceed ceed;
  FiniteElementKey fe;
  int qorder, nqpts, ncomp;
  BasisKey(Ceed ceed, const mfem::ParFiniteElementSpace &fespace,
           const mfem::FiniteElement &fe, const mfem::IntegrationRule &ir)
    : ceed(ceed), fe(fe), qorder(ir.GetOrder()), nqpts(ir.GetNPoints()),
      ncomp(fespace.GetVDim())
  {
  }
  bool operator==(const BasisKey &k) const
  {
    return (ceed == k.ceed && fe == k.fe && qorder == k.qorder && nqpts == k.nqpts &&
            ncomp == k.ncomp);
  }
};

struct BasisHash
{
  std::size_t operator()(const BasisKey &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.ceed, k.fe, k.qorder, k.nqpts, k.ncomp);
    return hash;
  }
};

struct InterpBasisKey
{
  Ceed ceed;
  FiniteElementKey trial_fe, test_fe;
  int ncomp;
  InterpBasisKey(Ceed ceed, const mfem::ParFiniteElementSpace &trial_fespace,
                 const mfem::ParFiniteElementSpace &test_fespace,
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
  std::size_t fespace, first_elem;
  bool use_bdr, unique_interp_restr, unique_interp_range_restr;
  RestrKey(Ceed ceed, const FiniteElementSpace &fespace, std::size_t first_elem,
           bool use_bdr, bool unique_interp_restr, bool unique_interp_range_restr)
    : ceed(ceed), fespace(fespace.GetId()), first_elem(first_elem), use_bdr(use_bdr),
      unique_interp_restr(unique_interp_restr),
      unique_interp_range_restr(unique_interp_range_restr)
  {
  }
  bool operator==(const RestrKey &k) const
  {
    return (ceed == k.ceed && fespace == k.fespace && first_elem == k.first_elem &&
            use_bdr == k.use_bdr && unique_interp_restr == k.unique_interp_restr &&
            unique_interp_range_restr == k.unique_interp_range_restr);
  }
};

struct RestrHash
{
  std::size_t operator()(const RestrKey &k) const
  {
    std::size_t hash = 0;
    CeedHashCombine(hash, k.ceed, k.fespace, k.first_elem, k.use_bdr, k.unique_interp_restr,
                    k.unique_interp_range_restr);
    return hash;
  }
};

}  // namespace internal

}  // namespace palace::ceed

namespace std
{

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

}  // namespace std

#endif  // PALACE_LIBCEED_HASH_HPP
