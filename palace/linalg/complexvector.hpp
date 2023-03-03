// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_COMPLEX_VECTOR_HPP
#define PALACE_COMPLEX_VECTOR_HPP

#include <type_traits>
#include <utility>
#include <mfem.hpp>

namespace palace
{

// A wrapper for a pair of mfem vectors, convenient when working with complex values
struct ComplexVector
{
  ComplexVector(int n) : real(n), imag(n) {}

  ComplexVector() = default;
  ComplexVector(const ComplexVector &) = default;
  ComplexVector(ComplexVector &&) = default;
  ComplexVector &operator=(const ComplexVector &) = default;
  ComplexVector &operator=(ComplexVector &&) = default;
  ~ComplexVector() = default;

  mfem::Vector real, imag;

  template <size_t I>
  auto &get() &
  {
    if constexpr (I == 0)
      return real;
    else if constexpr (I == 1)
      return imag;
  }

  template <size_t I>
  const auto &get() const &
  {
    if constexpr (I == 0)
      return real;
    else if constexpr (I == 1)
      return imag;
  }

  template <size_t I>
  auto &&get() &&
  {
    if constexpr (I == 0)
      return std::move(real);
    else if constexpr (I == 1)
      return std::move(imag);
  }
};

}  // namespace palace

// helpers to allow for structured binding
namespace std
{
template <>
struct tuple_size<palace::ComplexVector> : integral_constant<size_t, 2>
{
};
template <>
struct tuple_element<0, palace::ComplexVector>
{
  using type = mfem::Vector;
};
template <>
struct tuple_element<1, palace::ComplexVector>
{
  using type = mfem::Vector;
};
}  // namespace std

#endif  // PALACE_COMPLEX_VECTOR_HPP