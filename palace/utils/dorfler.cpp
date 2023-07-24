// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#include "dorfler.hpp"

#include <algorithm>
#include <numeric>
#include "utils/communication.hpp"

namespace palace::utils
{

double ComputeDorflerThreshold(double fraction, const mfem::Vector &e)
{
  // Pre compute the sort and partial sum to make evaluating a candidate
  // partition very fast.
  std::vector<double> estimates(e.begin(), e.end());
  std::sort(estimates.begin(), estimates.end());

  MFEM_ASSERT(estimates.front() >= 0, "Indicators must be non-negative");

  std::vector<double> sum;
  sum.reserve(estimates.size());
  std::partial_sum(estimates.begin(), estimates.end(), std::back_inserter(sum));

  auto pivot = std::distance(
      sum.begin(), std::lower_bound(sum.begin(), sum.end(), (1 - fraction) * sum.back()));

  double error_threshold = estimates[pivot];

  auto marked = [&estimates, &sum](double e) -> std::pair<std::size_t, double>
  {
    const auto lb = std::lower_bound(estimates.begin(), estimates.end(), e);
    const auto elems_marked = std::distance(lb, estimates.end());
    const double error_unmarked =
        lb != estimates.begin() ? sum[sum.size() - elems_marked - 1] : 0;
    const double error_marked = sum.back() - error_unmarked;
    return {elems_marked, error_marked};
  };

  // Each rank has computed a different threshold, if a given rank has lots of
  // low error elements, their value will be lower and if a rank has high error,
  // their value will be higher. Using the value from the low error rank will
  // give too many elements, and using the value from the high error rank will
  // give too few. The correct threshold value will be an intermediate
  // between the min and max over ranks.

  auto comm = Mpi::World();

  // Send initial information to all
  const double total_error = [&]()
  {
    double total_error = sum.back();
    Mpi::GlobalSum(1, &total_error, comm);
    return total_error;
  }();
  const double max_indicator = [&]()
  {
    double max_indicator = estimates.back();
    Mpi::GlobalSum(1, &max_indicator, comm);
    return max_indicator;
  }();

  double min_threshold = error_threshold;
  double max_threshold = error_threshold;
  Mpi::GlobalMin(1, &min_threshold, comm);
  Mpi::GlobalMax(1, &max_threshold, comm);
  const std::size_t total_elem = [&]()
  {
    std::size_t tmp = estimates.size();
    Mpi::GlobalSum(1, &tmp, comm);
    return tmp;
  }();

  MFEM_ASSERT(min_threshold <= max_threshold,
              "min: " << min_threshold << " max " << max_threshold);

  auto [elem_marked, error_marked] = marked(error_threshold);

  // Keep track of the number of elements marked by the threshold bounds. If the
  // top and bottom values are equal, there's no point further bisecting.
  auto [max_elem_marked, max_error_marked] = marked(min_threshold);
  auto [min_elem_marked, min_error_marked] = marked(max_threshold);
  Mpi::GlobalSum(1, &min_elem_marked, comm);
  Mpi::GlobalSum(1, &max_elem_marked, comm);
  Mpi::GlobalSum(1, &min_error_marked, comm);
  Mpi::GlobalSum(1, &max_error_marked, comm);

  constexpr int maxiter = 100;  // Maximum limit to prevent runaway.
  for (int i = 0; i < maxiter; ++i)
  {
    error_threshold = (min_threshold + max_threshold) / 2;

    std::tie(elem_marked, error_marked) = marked(error_threshold);

    // All processors need the values used for the stopping criteria.
    Mpi::GlobalSum(1, &elem_marked, comm);
    Mpi::GlobalSum(1, &error_marked, comm);

    MFEM_ASSERT(elem_marked > 0, "Some elements must have been marked");
    MFEM_ASSERT(error_marked > 0, "Some error must have been marked");

    const auto candidate_fraction = error_marked / total_error;

    Mpi::Debug("Threshold: {:e} < {:e} < {:e}, Marked Elems: {} <= {} <= {}\n",
               min_threshold, error_threshold, max_threshold, min_elem_marked, elem_marked,
               max_elem_marked);

    // Set the tolerance based off of the largest local indicator value. These
    // tolerance values are chosen based on testing, opt not to expose them.
    constexpr double frac_tol = 2 * std::numeric_limits<double>::epsilon();
    const double error_tol = 2 * std::numeric_limits<double>::epsilon() * max_indicator;
    if (std::abs(max_threshold - min_threshold) < error_tol ||
        std::abs(candidate_fraction - fraction) < frac_tol ||
        max_elem_marked == min_elem_marked)
    {
      // Candidate fraction matches to tolerance, or the number of marked
      // elements is no longer changing.
      Mpi::Debug(
          "ΔFraction: {:.3e}, Tol {:.3e}, ΔThreshold: {:.3e}, Tol {:.3e},  ΔElements: {}\n",
          candidate_fraction - fraction, frac_tol, max_threshold - min_threshold, error_tol,
          max_elem_marked - min_elem_marked);
      break;
    }

    // Update in preparation for next iteration. The logic here looks inverted
    // compared to a usual binary search, because a smaller value marks a larger
    // number of elements and thus a greater fraction of error.
    if (candidate_fraction > fraction)
    {
      // This candidate marked too much, raise the lower bound.
      min_threshold = error_threshold;
      max_elem_marked = elem_marked;
      max_error_marked = error_marked;
    }
    else if (candidate_fraction < fraction)
    {
      // This candidate marked too little, lower the upper bound.
      max_threshold = error_threshold;
      min_elem_marked = elem_marked;
      min_error_marked = error_marked;
    }
  }

  // Always choose the lower threshold value, thereby marking the larger number
  // of elements and fraction of the total error. Would rather over mark than
  // under mark, as Dorfler marking is the smallest set that covers AT LEAST the
  // specified fraction of the error.

  error_threshold = min_threshold;
  elem_marked = max_elem_marked;
  error_marked = max_error_marked;

  Mpi::Print("Indicator threshold {:.3e} marked {} of {} elements and {:.3f}% error\n",
             error_threshold, elem_marked, total_elem, 100 * error_marked / total_error);

  MFEM_ASSERT(error_threshold > 0, "error_threshold must be positive");
  return error_threshold;
}

}  // namespace palace::utils
