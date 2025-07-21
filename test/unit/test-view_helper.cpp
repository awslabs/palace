#include <iterator>
#include <map>
#include <vector>
#include <fmt/format.h>
#include <catch2/catch_test_macros.hpp>
#include <catch2/benchmark/catch_benchmark_all.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include "utils/view_helper.hpp"

using namespace palace;

TEST_CASE("EnumerateViewVectorBasic", "[enumerateview]")
{
  // Vector no off-set
  {
    std::vector<double> test_vector = {1.0, 2.0, 4.0, 8.0};
    std::vector<std::pair<int, double>> out_vector{};

    for (const auto &v_it : EnumerateView(test_vector))
    {
      out_vector.emplace_back(v_it);
    }
    CHECK(out_vector ==
          std::vector<std::pair<int, double>>{{0, 1.0}, {1, 2.0}, {2, 4.0}, {3, 8.0}});
  }

  // Vector with off-set
  {
    std::vector<double> test_vector = {1.0, 2.0, 4.0, 8.0};
    std::vector<std::pair<int, double>> out_vector{};

    auto view = EnumerateView(test_vector, 2);
    CHECK(view.begin() != view.end());

    for (const auto &v_it : EnumerateView(test_vector, 2))
    {
      out_vector.emplace_back(v_it);
    }
    CHECK(out_vector == std::vector<std::pair<int, double>>{
                            {2, 4.0},
                            {3, 8.0},
                        });
  }

  // Vector with off-set at end
  {
    std::vector<double> test_vector = {1.0, 2.0, 4.0, 8.0};
    std::vector<std::pair<int, double>> out_vector{};

    auto view = EnumerateView(test_vector, 4);
    CHECK(view.begin() == view.end());

    for (const auto &v_it : EnumerateView(test_vector, 4))
    {
      out_vector.emplace_back(v_it);
    }
    CHECK(out_vector == std::vector<std::pair<int, double>>{});
  }

  // Vector with off-set beyond end
  {
    std::vector<double> test_vector = {1.0, 2.0, 4.0, 8.0};
    std::vector<std::pair<int, double>> out_vector{};

    auto view = EnumerateView(test_vector, 8);
    REQUIRE(view.begin() == view.end());

    for (const auto &v_it : EnumerateView(test_vector, 8))
    {
      out_vector.emplace_back(v_it);
    }
    CHECK(out_vector == std::vector<std::pair<int, double>>{});
  }
}

TEST_CASE("EnumerateViewVectorTypeConstChecks", "[enumerateview]")
{
  {
    std::vector<double> test_vector = {1.0, 2.0, 4.0, 8.0};
    auto view = EnumerateView(test_vector);

    STATIC_CHECK(std::is_same_v<decltype(view), EnumerateView<std::vector<double>>>);
    STATIC_CHECK(std::is_same_v<decltype(view)::container_t, std::vector<double>>);

    STATIC_CHECK(std::is_same_v<decltype(std::begin(view))::value_type,
                                std::pair<std::size_t, double>>);
    STATIC_CHECK(std::is_same_v<decltype(std::begin(view))::reference,
                                std::pair<std::size_t, double &>>);

    for (auto [i_it, v_it] : EnumerateView(test_vector))
    {
      STATIC_CHECK(std::is_same_v<decltype(i_it), std::size_t>);
      STATIC_CHECK(std::is_same_v<decltype(v_it), double &>);
      v_it *= 2;
    }
    CHECK(test_vector == std::vector{2.0, 4.0, 8.0, 16.0});

    // Having const in loop will make the tuple can by value iterator const but not the
    // reference to the double, which we can still modify.
    for (const auto &[i_it, v_it] : EnumerateView(test_vector))
    {
      STATIC_CHECK(std::is_same_v<decltype(i_it), const std::size_t>);
      STATIC_CHECK(std::is_same_v<decltype(v_it), double &>);
      v_it *= 2;
    }
    CHECK(test_vector == std::vector{4.0, 8.0, 16.0, 32.0});
  }

  {
    const std::vector<double> test_vector{1.0, 2.0, 4.0, 8.0};
    auto view = EnumerateView(test_vector);

    STATIC_CHECK(std::is_same_v<decltype(view), EnumerateView<const std::vector<double>>>);
    STATIC_CHECK(std::is_same_v<decltype(view)::container_t, const std::vector<double>>);

    STATIC_CHECK(std::is_same_v<decltype(std::begin(view))::value_type,
                                std::pair<std::size_t, double>>);
    STATIC_CHECK(std::is_same_v<decltype(std::begin(view))::reference,
                                std::pair<std::size_t, const double &>>);

    for (auto [i_it, v_it] : EnumerateView(test_vector))
    {
      STATIC_CHECK(std::is_same_v<decltype(i_it), std::size_t>);
      STATIC_CHECK(std::is_same_v<decltype(v_it), const double &>);
    }
  }

  {
    struct my_struct
    {
      double test;

      bool operator==(const my_struct &other) const { return test == other.test; }
    };
    std::vector<my_struct> test_vector{{1.0}, {2.0}, {4.0}, {8.0}};
    auto view = EnumerateView(test_vector);

    STATIC_CHECK(std::is_same_v<decltype(view), EnumerateView<std::vector<my_struct>>>);
    STATIC_CHECK(std::is_same_v<decltype(view)::container_t, std::vector<my_struct>>);

    STATIC_CHECK(std::is_same_v<decltype(std::begin(view))::value_type,
                                std::pair<std::size_t, my_struct>>);
    STATIC_CHECK(std::is_same_v<decltype(std::begin(view))::reference,
                                std::pair<std::size_t, my_struct &>>);

    for (auto [i_it, v_it] : EnumerateView(test_vector))
    {
      STATIC_CHECK(std::is_same_v<decltype(i_it), std::size_t>);
      STATIC_CHECK(std::is_same_v<decltype(v_it), my_struct &>);
      v_it.test *= 2;
    }
    CHECK(test_vector == std::vector<my_struct>{{2.0}, {4.0}, {8.0}, {16.0}});
  }
}

TEST_CASE("EnumerateViewMapBasic", "[enumerateview]")
{
  // Map no off-set
  {
    std::map<int, double> test_map{{1, 1.0}, {2, 2.0}, {3, 4.0}, {4, 8.0}};
    std::vector<std::pair<int, std::pair<int, double>>> out_vector;
    for (const auto &v_it : EnumerateView(test_map))
    {
      out_vector.emplace_back(v_it);
    }
    CHECK(out_vector ==
          decltype(out_vector){{0, {1, 1.0}}, {1, {2, 2.0}}, {2, {3, 4.0}}, {3, {4, 8.0}}});
  }

  // Map with off-set
  {
    std::map<int, double> test_map{{1, 1.0}, {2, 2.0}, {3, 4.0}, {4, 8.0}};
    std::vector<std::pair<int, std::pair<int, double>>> out_vector;

    auto view = EnumerateView(test_map, 2);
    CHECK(view.begin() != view.end());

    for (const auto &v_it : EnumerateView(test_map, 2))
    {
      out_vector.emplace_back(v_it);
    }
    CHECK(out_vector == decltype(out_vector){{2, {3, 4.0}}, {3, {4, 8.0}}});
  }

  // Map with off-set at end
  {
    std::map<int, double> test_map{{1, 1.0}, {2, 2.0}, {3, 4.0}, {4, 8.0}};
    std::vector<std::pair<int, std::pair<int, double>>> out_vector;

    auto view = EnumerateView(test_map, 4);
    CHECK(view.begin() == view.end());

    for (const auto &v_it : EnumerateView(test_map, 4))
    {
      out_vector.emplace_back(v_it);
    }
    CHECK(out_vector == decltype(out_vector){});
  }

  // Map with off-set beyond end
  {
    std::map<int, double> test_map{{1, 1.0}, {2, 2.0}, {3, 4.0}, {4, 8.0}};
    std::vector<std::pair<int, std::pair<int, double>>> out_vector;

    auto view = EnumerateView(test_map, 8);
    REQUIRE(view.begin() == view.end());

    for (const auto &v_it : EnumerateView(test_map, 8))
    {
      out_vector.emplace_back(v_it);
    }
    CHECK(out_vector == decltype(out_vector){});
  }
}
