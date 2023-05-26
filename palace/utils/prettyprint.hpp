// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_PRETTY_PRINT_HPP
#define PALACE_UTILS_PRETTY_PRINT_HPP

#include <string>
#include <type_traits>
#include <vector>
#include <mfem.hpp>
#include "utils/communication.hpp"

namespace palace::utils
{

//
// Utility functions for formatted printing.
//

namespace internal
{

constexpr std::size_t max_width = 60;

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

// Fixed column width wrapped printing with range notation for the contents of a marker
// array.
template <template <typename...> class Container, typename T, typename... U>
inline void PrettyPrintMarker(const Container<T, U...> &data,
                              const std::string &prefix = "",
                              MPI_Comm comm = MPI_COMM_WORLD)
{
  static_assert(std::is_integral<T>::value,
                "PrettyPrintMarker requires containers with an integral type marker!");
  std::size_t i = 0, w = 0, lead = prefix.length();
  Mpi::Print(comm, prefix);
  while (i < internal::GetSize(data))
  {
    if (data[i])
    {
      auto j = i;
      while ((j + 1 < internal::GetSize(data)) && data[j + 1])
      {
        j++;
      }
      if (i == j)
      {
        auto wi = 1 + static_cast<int>(std::log10(i + 1));
        w = internal::PrePrint(comm, w, wi, lead) + wi;
        Mpi::Print(comm, "{:d}", i + 1);
        i++;
      }
      else
      {
        auto wi =
            3 + static_cast<int>(std::log10(i + 1)) + static_cast<int>(std::log10(j + 1));
        w = internal::PrePrint(comm, w, wi, lead) + wi;
        Mpi::Print(comm, "{:d}-{:d}", i + 1, j + 1);
        i = j + 1;
      }
    }
    else
    {
      i++;
    }
  }
  Mpi::Print(comm, "\n");
}

// Fixed column width wrapped printing for the contents of an array.
template <template <typename...> class Container, typename T, typename... U>
inline void PrettyPrint(const Container<T, U...> &data, T scale = T(1.0),
                        const std::string &prefix = "", MPI_Comm comm = MPI_COMM_WORLD)
{
  constexpr int pv = 3;       // Value precision
  constexpr int wv = pv + 6;  // Total printed width of a value
  std::size_t w = 0, lead = prefix.length();
  Mpi::Print(comm, prefix);
  for (const auto &v : data)
  {
    w = internal::PrePrint(comm, w, wv, lead) + wv;
    Mpi::Print(comm, "{:.{}e}", v * scale, pv);
  }
  Mpi::Print(comm, "\n");
}

}  // namespace palace::utils

#endif  // PALACE_UTILS_PRETTY_PRINT_HPP
