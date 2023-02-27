// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_COMPLEX_VECTOR_HPP
#define PALACE_COMPLEX_VECTOR_HPP

#include <utility>
#include <type_traits>
#include <mfem.hpp>

namespace palace {

// A wrapper for a pair of mfem vectors, convenient when working with complex values
struct ComplexVector
{
    mfem::Vector real, imag;

    template <size_t I>
    auto& get() & {
        if constexpr (I == 0) return real;
        else if constexpr (I == 1) return imag;
    }

    template <size_t I>
    auto const& get() const& {
        if constexpr (I == 0) return real;
        else if constexpr (I == 1) return imag;
    }

    template <size_t I>
    auto&& get() && {
        if constexpr (I == 0) return std::move(real);
        else if constexpr (I == 1) return std::move(imag);
    }
};

}

// helpers to allow for structured binding
namespace std
{
  template<> struct tuple_size<palace::ComplexVector> : integral_constant<size_t, 2> {};
  template<> struct tuple_element<0, palace::ComplexVector> { using type = mfem::Vector; };
  template<> struct tuple_element<1, palace::ComplexVector> { using type = mfem::Vector; };
}

#endif // PALACE_COMPLEX_VECTOR_HPP