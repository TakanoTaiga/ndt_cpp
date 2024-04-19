//   Copyright 2024 Shota Minami
//
//   Licensed under the Apache License, Version 2.0 (the "License");
//   you may not use this file except in compliance with the License.
//   You may obtain a copy of the License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the License is distributed on an "AS IS" BASIS,
//   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//   See the License for the specific language governing permissions and
//   limitations under the License.

#ifndef FLATKDTREE_FLATKDTREE_H_
#define FLATKDTREE_FLATKDTREE_H_

#include <algorithm>
#include <array>
#include <cstddef>
#include <iterator>
#include <type_traits>

namespace kdtree {

template <typename T, typename Enabler = void>
struct default_point_policy;

template <typename Float, std::size_t Size>
struct default_point_policy<std::array<Float, Size>, typename std::enable_if<std::is_floating_point<Float>::value>::type>
{
  using point_type = std::array<Float, Size>;
  using distance_type = Float;

  static constexpr std::size_t dimension = Size;

  template <std::size_t Index>
  auto element_compare(const point_type &p, const point_type &q) const -> bool
  {
    return p[Index] < q[Index];
  }

  template <std::size_t Index>
  auto element_distance(const point_type &p, const point_type &q) const -> distance_type
  {
    return (p[Index] - q[Index]) * (p[Index] - q[Index]);
  }

  auto distance(const point_type &p, const point_type &q) const -> distance_type
  {
    distance_type dist = 0;
    for (std::size_t i = 0; i < Size; ++i) {
      dist += (p[i] - q[i]) * (p[i] - q[i]);
    }
    return dist;
  }
};

namespace internal {
  template <typename Iterator>
  using iter_value_type = typename std::iterator_traits<Iterator>::value_type;

  template <int Dimension = 0, typename RandomAccessIterator, typename Policy>
  auto construct_recursive(
    RandomAccessIterator first, RandomAccessIterator last, const Policy &policy) -> void
  {
    using point_type = typename Policy::point_type;

    if (first == last) {
      return;
    }

    const RandomAccessIterator middle = first + std::distance(first, last) / 2;

    std::nth_element(first, middle, last, [&](const point_type &p, const point_type &q) {
      return policy.template element_compare<Dimension>(p, q);
    });

    static constexpr std::size_t NextDimension = (Dimension + 1) % Policy::dimension;
    construct_recursive<NextDimension>(first, middle, policy);
    construct_recursive<NextDimension>(middle + 1, last, policy);
  }

  template <
    std::size_t Dimension = 0,
    typename RandomAccessIterator1, typename RandomAccessIterator2, typename RandomAccessIterator3,
    typename PointPolicy = default_point_policy<iter_value_type<RandomAccessIterator1>>>
  auto search_knn_recursive(
    RandomAccessIterator1 first, RandomAccessIterator1 last,
    RandomAccessIterator2 out_point, RandomAccessIterator3 out_distance,
    std::size_t k, std::size_t n, const internal::iter_value_type<RandomAccessIterator1> &query,
    const PointPolicy &policy = {})
    -> std::size_t
  {
    using distance_type = typename PointPolicy::distance_type;

    if (first == last) {
      return n;
    }

    const RandomAccessIterator1 middle = first + std::distance(first, last) / 2;
    distance_type distance = policy.distance(query, *middle);

    if (n < k) {
      // Insert the point to the heap.
      std::size_t i = n;
      while (i > 0) {
        std::size_t p = (i - 1) >> 1;
        if (distance < out_distance[p]) {
          break;
        }

        out_distance[i] = std::move(out_distance[p]);
        out_point[i] = std::move(out_point[p]);
        i = p;
      }
      out_distance[i] = std::move(distance);
      out_point[i] = *middle;
      ++n;
    } else if (distance < *out_distance) {
      // Replace the root of the heap.
      std::size_t p = 0;
      for (std::size_t i = 1; i < n; i = (p << 1) | 1) {
        if (i + 1 < n and out_distance[i] < out_distance[i + 1]) {
          ++i;
        }

        if (out_distance[i] < distance) {
          break;
        }

        out_distance[p] = std::move(out_distance[i]);
        out_point[p] = std::move(out_point[i]);
        p = i;
      }
      out_distance[p] = std::move(distance);
      out_point[p] = *middle;
    }

    const auto search = [&](RandomAccessIterator1 first, RandomAccessIterator1 last) {
      constexpr auto next_dimension = (Dimension + 1) % PointPolicy::dimension;
      return search_knn_recursive<next_dimension>(first, last, out_point, out_distance, k, n, query, policy);
    };

    if (policy.template element_compare<Dimension>(query, *middle)) {
      n = search(first, middle);

      if (policy.template element_distance<Dimension>(*middle, query) < *out_distance) {
        n = search(middle + 1, last);
      }
    } else {
      n = search(middle + 1, last);

      if (policy.template element_distance<Dimension>(*middle, query) < *out_distance) {
        n = search(first, middle);
      }
    }

    return n;
  }
}

/**
 * Construct a k-d tree on a range of points.
 *
 * @tparam RandomAccessIterator The type of the iterator.
 * @tparam PointPolicy The policy for the point type.
 * @param first The iterator to the first element.
 * @param last The iterator to the last element.
 * @param policy The policy for the point type.
 */
template <
  typename RandomAccessIterator,
  typename PointPolicy = default_point_policy<internal::iter_value_type<RandomAccessIterator>>>
auto construct(RandomAccessIterator first, RandomAccessIterator last, const PointPolicy &policy = {})
  -> void
{
  using point_type = typename PointPolicy::point_type;

  static_assert(std::is_same<internal::iter_value_type<RandomAccessIterator>, point_type>::value,
    "kd-tree point type is not compatible with the policy");

  internal::construct_recursive(first, last, policy);
}

/**
 * Search k-nearest neighbors.
 *
 * @tparam RandomAccessIterator1 The type of the iterator for the points.
 * @tparam RandomAccessIterator2 The type of the iterator for the output points.
 * @tparam RandomAccessIterator3 The type of the iterator for the output distances.
 * @tparam PointPolicy The policy for the point type.
 * @param first The iterator to the first element.
 * @param last The iterator to the last element.
 * @param out_point The iterator to the first element of the output points.
 * @param out_distance The iterator to the first element of the output distances.
 * @param k The number of neighbors to search.
 * @param query The query point.
 * @param policy The policy for the point type.
 * @return The number of found neighbors.
 */
template <
  typename RandomAccessIterator1, typename RandomAccessIterator2, typename RandomAccessIterator3,
  typename PointPolicy = default_point_policy<internal::iter_value_type<RandomAccessIterator1>>>
auto search_knn(
  RandomAccessIterator1 first, RandomAccessIterator1 last,
  RandomAccessIterator2 out_point, RandomAccessIterator3 out_distance,
  std::size_t k, const internal::iter_value_type<RandomAccessIterator1> &query,
  const PointPolicy &policy = {})
  -> std::size_t
{
  using point_type = typename PointPolicy::point_type;
  using distance_type = typename PointPolicy::distance_type;

  static_assert(std::is_same<internal::iter_value_type<RandomAccessIterator1>, point_type>::value,
    "kd-tree point type is not compatible with the policy");
  static_assert(std::is_same<internal::iter_value_type<RandomAccessIterator2>, point_type>::value,
    "output point type is not compatible with the policy");
  static_assert(std::is_same<internal::iter_value_type<RandomAccessIterator3>, distance_type>::value,
    "output distance type is not compatible with the policy");

  return internal::search_knn_recursive(first, last, out_point, out_distance, k, 0, query, policy);
}

} // namespace kdtree

#endif // FLATKDTREE_FLATKDTREE_H_
