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
#include <tuple>
#include <type_traits>
#include <utility>

namespace kdtree {

namespace trait {
  template <typename T, std::size_t I, typename Enabler = void>
  struct access
  { };

  template <typename T>
  struct access<T, 0, typename std::enable_if<std::is_arithmetic<T>::value>::type>
  {
    static T get(const T &t)
    {
      return t;
    }
  };

  namespace detail {
    template <typename T, std::size_t = std::tuple_size<T>::value>
    using require_tuple = void;
  }

  template <typename T, std::size_t I>
  struct access<T, I, detail::require_tuple<T>>
  {
    static auto get(const T &t) -> typename std::tuple_element<I, T>::type
    {
      using std::get;
      return get<I>(t);
    }
  };

  template <typename T, typename Enabler = void>
  struct dimension
  { };

  template <typename T>
  struct dimension<T, typename std::enable_if<std::is_arithmetic<T>::value>::type>
  {
    static constexpr std::size_t value = 1;
  };

  template <typename T>
  struct dimension<T, detail::require_tuple<T>>
  {
    static constexpr std::size_t value = std::tuple_size<T>::value;
  };

  template <typename S, typename T, typename Enabler = void>
  struct compare
  {
    static bool apply(const S &lhs, const T &rhs)
    {
      return lhs < rhs;
    }
  };

  namespace detail {
    template <typename R, typename S, typename T, std::size_t L>
    struct default_squared_distance
    {
      static R apply(const S &s, const T &t)
      {
        R d = static_cast<R>(access<S, L>::get(s)) - static_cast<R>(access<T, L>::get(t));
        return d * d + default_squared_distance<R, S, T, L - 1>::apply(s, t);
      }
    };

    template <typename R, typename S, typename T>
    struct default_squared_distance<R, S, T, 0>
    {
      static R apply(const S &s, const T &t)
      {
        R d = static_cast<R>(access<S, 0>::get(s)) - static_cast<R>(access<T, 0>::get(t));
        return d * d;
      }
    };
  }

  template <typename R, typename S, typename T, typename Enabler = void>
  struct squared_distance
  {
    static R apply(const S &s, const T &t)
    {
      return detail::default_squared_distance<R, S, T, dimension<S>::value - 1>::apply(s, t);
    }
  };
}

namespace internal {
  template <std::size_t I, typename T>
  auto get(const T &t) -> decltype(trait::access<T, I>::get(t)) { return trait::access<T, I>::get(t); }

  template <typename T>
  constexpr std::size_t dimension() { return trait::dimension<T>::value; }

  template <typename S, typename T>
  bool compare(const S &lhs, const T &rhs) { return trait::compare<S, T>::apply(lhs, rhs); }

  template <typename R, typename S, typename T>
  R squared_distance(const S &s, const T &t) { return trait::squared_distance<R, S, T>::apply(s, t); }

  template <std::size_t L, typename I>
  void do_construct(I first, I last)
  {
    using P = typename std::iterator_traits<I>::value_type;

    const I middle = first + std::distance(first, last) / 2;
    std::nth_element(first, middle, last, [&](const P &l, const P &r) { return compare(get<L>(l), get<L>(r)); });

    constexpr std::size_t NL = (L + 1) % dimension<P>();
    if (first != middle) {
      do_construct<NL>(first, middle);
    }
    if (middle + 1 != last) {
      do_construct<NL>(middle + 1, last);
    }
  }

  template <std::size_t L, typename I, typename OP, typename OD, typename Q>
  std::size_t do_search_knn(I first, I last, OP out_point, OD out_distance, std::size_t k, std::size_t n, const Q &query)
  {
    using P = typename std::iterator_traits<I>::value_type;
    using D = typename std::iterator_traits<OD>::value_type;

    const I middle = first + std::distance(first, last) / 2;
    const D distance = squared_distance<D>(query, *middle);

    if (n < k) {
      // Insert the point to the heap.
      std::size_t i = n;
      while (i > 0) {
        const std::size_t p = (i - 1) >> 1;
        if (compare(distance, out_distance[p])) {
          break;
        }

        out_distance[i] = std::move(out_distance[p]);
        out_point[i] = std::move(out_point[p]);
        i = p;
      }
      out_distance[i] = std::move(distance);
      out_point[i] = *middle;
      ++n;
    } else if (compare(distance, *out_distance)) {
      // Replace the root of the heap.
      std::size_t p = 0;
      for (std::size_t i = 1; i < n; i = (p << 1) | 1) {
        if (i + 1 < n and compare(out_distance[i], out_distance[i + 1])) {
          ++i;
        }

        if (compare(out_distance[i], distance)) {
          break;
        }

        out_distance[p] = std::move(out_distance[i]);
        out_point[p] = std::move(out_point[i]);
        p = i;
      }
      out_distance[p] = std::move(distance);
      out_point[p] = *middle;
    }

    constexpr std::size_t NL = (L + 1) % dimension<P>();
    if (compare(get<L>(query), get<L>(*middle))) {
      if (first != middle) {
        n = do_search_knn<NL>(first, middle, out_point, out_distance, k, n, query);
      }
      if (middle + 1 != last and compare(squared_distance<D>(get<L>(query), get<L>(*middle)), *out_distance)) {
        n = do_search_knn<NL>(middle + 1, last, out_point, out_distance, k, n, query);
      }
    } else {
      if (middle + 1 != last) {
        n = do_search_knn<NL>(middle + 1, last, out_point, out_distance, k, n, query);
      }
      if (first != middle and compare(squared_distance<D>(get<L>(query), get<L>(*middle)), *out_distance)) {
        n = do_search_knn<NL>(first, middle, out_point, out_distance, k, n, query);
      }
    }
    return n;
  }
}

template <typename I>
void construct(I first, I last)
{
  if (first != last) {
    internal::do_construct<0>(first, last);
  }
}

template <typename I, typename OP, typename OD, typename Q>
std::size_t search_knn(I first, I last, OP out_point, OD out_distance, std::size_t k, const Q &query)
{
  if (first != last) {
    return internal::do_search_knn<0>(first, last, out_point, out_distance, k, 0, query);
  }
  return 0;
}

} // namespace kdtree

#endif // FLATKDTREE_FLATKDTREE_H_
