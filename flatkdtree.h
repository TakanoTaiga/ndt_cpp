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
#include <utility>

namespace kdtree {

namespace trait {
  template <typename T, typename Enabler = void>
  struct access
  { };

  template <typename T, std::size_t N>
  struct access<std::array<T, N>, typename std::enable_if<std::is_arithmetic<T>::value>::type>
  {
    static constexpr std::size_t dimension = N;

    template <std::size_t I>
    static auto get(const std::array<T, N> &t) -> T
    {
      return t[I];
    }
  };
}

namespace internal {
  template <std::size_t I, typename T>
  auto get(const T &t) -> decltype(trait::access<typename std::decay<T>::type>::template get<I>(t))
  {
    return trait::access<typename std::decay<T>::type>::template get<I>(t);
  }

  namespace detail {
    template <typename S, typename T, std::size_t I>
    struct distance
    {
      static auto apply(const S &s, const T &t) -> decltype((get<I>(s) - get<I>(t)) * (get<I>(s) - get<I>(t)))
      {
        return (get<I>(s) - get<I>(t)) * (get<I>(s) - get<I>(t)) + distance<S, T, I - 1>::apply(s, t);
      }
    };

    template <typename S, typename T>
    struct distance<S, T, 0>
    {
      static auto apply(const S &s, const T &t) -> decltype((get<0>(s) - get<0>(t)) * (get<0>(s) - get<0>(t)))
      {
        return (get<0>(s) - get<0>(t)) * (get<0>(s) - get<0>(t));
      }
    };
  }

  template <typename S, typename T>
  auto distance(const S &s, const T &t) -> decltype(detail::distance<S, T, trait::access<S>::dimension - 1>::apply(s, t))
  {
    static_assert(trait::access<S>::dimension == trait::access<T>::dimension, "dimension mismatch");
    return detail::distance<S, T, trait::access<S>::dimension - 1>::apply(s, t);
  }

  template <std::size_t D, typename I>
  auto do_construct(I first, I last) -> void
  {
    using P = typename std::iterator_traits<I>::value_type;

    const auto middle = first + std::distance(first, last) / 2;

    std::nth_element(first, middle, last, [&](const P &l, const P &r) {
      return get<D>(l) < get<D>(r);
    });

    constexpr auto ND = (D + 1) % trait::access<P>::dimension;
    if (first != middle) {
      do_construct<ND>(first, middle);
    }
    if (middle + 1 != last) {
      do_construct<ND>(middle + 1, last);
    }
  }

  template <std::size_t D, typename I1, typename I2, typename I3, typename Q>
  auto do_search_knn(I1 first, I1 last, I2 out_point, I3 out_distance, std::size_t k, std::size_t n, const Q &query)
    -> std::size_t
  {
    using P = typename std::iterator_traits<I1>::value_type;

    const auto middle = first + std::distance(first, last) / 2;
    auto distance = internal::distance(query, *middle);

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

    constexpr auto ND = (D + 1) % trait::access<Q>::dimension;
    if (get<D>(query) < get<D>(*middle)) {
      if (first != middle) {
        n = do_search_knn<ND>(first, middle, out_point, out_distance, k, n, query);
      }
      if (middle + 1 != last and (get<D>(query) - get<D>(*middle)) * (get<D>(query) - get<D>(*middle)) < *out_distance) {
        n = do_search_knn<ND>(middle + 1, last, out_point, out_distance, k, n, query);
      }
    } else {
      if (middle + 1 != last) {
        n = do_search_knn<ND>(middle + 1, last, out_point, out_distance, k, n, query);
      }
      if (first != middle and (get<D>(query) - get<D>(*middle)) * (get<D>(query) - get<D>(*middle)) < *out_distance) {
        n = do_search_knn<ND>(first, middle, out_point, out_distance, k, n, query);
      }
    }

    return n;
  }
}

template <typename I>
auto construct(I first, I last) -> void
{
  if (first != last) {
    internal::do_construct<0>(first, last);
  }
}

template <typename I1, typename I2, typename I3, typename Q>
auto search_knn(I1 first, I1 last, I2 out_point, I3 out_distance, std::size_t k, const Q &query) -> std::size_t
{
  if (first != last) {
    return internal::do_search_knn<0>(first, last, out_point, out_distance, k, 0, query);
  }
  return 0;
}

} // namespace kdtree

#endif // FLATKDTREE_FLATKDTREE_H_
