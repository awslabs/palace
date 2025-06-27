// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_UTILS_VIEW_HELPER_HPP
#define PALACE_UTILS_VIEW_HELPER_HPP

#include <cstddef>
#include <functional>
#include <iterator>
#include <utility>

namespace palace
{

// Naive C++17 helper class to implement "enum", "drop" view (in that order). It only does
// forward iteration (no random access / bidirectional iteration). Replace in C++20/C++23 or
// with a range library.
//
// Only used in drivensolver.cpp but included as separate header for unit tests.
//
// Assume Container is sized iterator.
template <typename Container>
class EnumerateView
{
public:
  using container_t = std::remove_reference_t<Container>;
  container_t *range_;
  std::size_t drop_;

  using it_t = decltype(std::declval<container_t>().begin());  //::iterator;

public:
  template <typename QualifiedContainer,
            typename = std::enable_if_t<
                std::is_same_v<std::remove_reference_t<QualifiedContainer>, Container>>>
  EnumerateView(QualifiedContainer &&range, std::size_t drop = 0)
    : range_(std::addressof(range)), drop_(drop)
  {
  }

  // Custom Iterator with step increment
  class EnumerateViewIt
  {
  private:
    it_t it_;
    it_t it_end_;
    std::size_t counter_;

  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type =
        std::pair<std::size_t, typename std::iterator_traits<it_t>::value_type>;
    using difference_type = typename std::iterator_traits<it_t>::difference_type;
    using reference =
        std::pair<std::size_t, typename std::iterator_traits<it_t>::reference>;
    using pointer = reference *;

    // Precondition: drop and it location must match correctly.
    EnumerateViewIt(it_t it, it_t it_end, std::size_t drop)
      : it_(it), it_end_(it_end), counter_(drop)
    {
    }

    reference operator*() { return std::make_pair(counter_, std::ref(*it_)); }

    // Increment counter and range iterator.
    auto &operator++()
    {
      ++it_;
      ++counter_;
      return *this;
    }

    auto operator++(int)
    {
      EnumerateViewIt temp = *this;
      ++(*this);
      return temp;
    }

    auto operator==(const EnumerateViewIt &other) const { return it_ == other.it_; }
    auto operator!=(const EnumerateViewIt &other) const { return it_ != other.it_; }
  };

  auto begin()
  {
    // Note: std::advance beyond end is undefined behavior. Clamp to end.
    if (drop_ > range_->size())
    {
      return end();
    }
    auto it_ = range_->begin();
    std::advance(it_, drop_);
    return EnumerateViewIt(it_, range_->end(), drop_);
  }

  auto end() { return EnumerateViewIt(range_->end(), range_->end(), range_->size()); }
};

// Deduction guides so it knows how tot get Container from QualifiedContainer.
template <typename QualifiedContainer>
EnumerateView(QualifiedContainer &&, std::size_t)
    -> EnumerateView<std::remove_reference_t<QualifiedContainer>>;

template <typename QualifiedContainer>
EnumerateView(QualifiedContainer &&)
    -> EnumerateView<std::remove_reference_t<QualifiedContainer>>;

}  // namespace palace

#endif  // PALACE_UTILS_VIEW_HELPER_HPP
