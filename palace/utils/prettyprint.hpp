// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PRETTY_PRINT_HPP
#define PALACE_UTILS_PRETTY_PRINT_HPP

#include <string>
#include <string_view>
#include <type_traits>
#include <vector>
#include <fmt/core.h>
#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace::utils
{

//
// Utility functions for formatted printing.
//

namespace internal
{

constexpr std::size_t max_width = 80;

template <typename T>
inline std::size_t GetSize(const T &v)
{
  return v.size();
}

template <typename T>
inline std::size_t GetSize(const mfem::Array<T> &v)
{
  return v.Size();
}

inline std::size_t PrePrint(MPI_Comm comm, std::size_t w, std::size_t wv, std::size_t lead)
{
  auto end = w + 2 + wv + 1;  // Consider space for trailing comma
  if (w > 0 && end > max_width - lead)
  {
    Mpi::Print(comm, ",\n{}", std::string(lead, ' '));  // Line break
    w = 0;
  }
  if (w)
  {
    Mpi::Print(comm, ", ");
    return w + 2;
  }
  else
  {
    Mpi::Print(comm, " ");
    return w + 1;
  }
}

}  // namespace internal

// Fixed column width wrapped printing for the contents of an array, with range notation for
// integral types.
template <template <typename...> class Container, typename T, typename... U>
inline void PrettyPrint(const Container<T, U...> &data, T scale,
                        const std::string_view prefix = "", MPI_Comm comm = MPI_COMM_WORLD)
{
  std::size_t w = 0, lead = prefix.length();
  Mpi::Print(comm, fmt::runtime(prefix));
  auto i = data.begin();
  while (i != data.end())
  {
    if constexpr (std::is_integral<T>::value)
    {
      auto j = i;
      if (scale == 1)
      {
        while ((std::next(j) != data.end()) && *std::next(j) == (*j) + 1)
        {
          j++;
        }
      }
      if (i == j)
      {
        auto wi = 1 + static_cast<T>(std::log10((*i) + 1));
        w = internal::PrePrint(comm, w, wi, lead) + wi;
        Mpi::Print(comm, "{:d}", (*i) * scale);
      }
      else
      {
        auto wi =
            3 + static_cast<T>(std::log10((*i) + 1)) + static_cast<T>(std::log10((*j) + 1));
        w = internal::PrePrint(comm, w, wi, lead) + wi;
        Mpi::Print(comm, "{:d}-{:d}", (*i) * scale, (*j) * scale);
      }
      i = std::next(j);
    }
    else
    {
      constexpr auto pv = 3;       // Value precision
      constexpr auto wv = pv + 6;  // Total printed width of a value
      w = internal::PrePrint(comm, w, wv, lead) + wv;
      Mpi::Print(comm, "{:.{}e}", (*i) * scale, pv);
      i++;
    }
  }
  Mpi::Print(comm, "\n");
}

template <template <typename...> class Container, typename T, typename... U>
inline void PrettyPrint(const Container<T, U...> &data, const std::string &prefix = "",
                        MPI_Comm comm = MPI_COMM_WORLD)
{
  PrettyPrint(data, T(1), prefix, comm);
}

}  // namespace palace::utils

#endif  // PALACE_UTILS_PRETTY_PRINT_HPP
