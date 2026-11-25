/**
 * Shadow header for icpx + SYCL device compilation vs system libfmt 9.
 * On the NVPTX device pass, long double keeps sizeof==16 but digits==53,
 * so stock dragonbox::float_info has no match. Provide a uint128 specialization.
 */
#pragma once

#include_next <fmt/format.h>

#if defined(__SYCL_DEVICE_ONLY__)
namespace fmt {
namespace detail {
namespace dragonbox {
template <>
struct float_info<long double> {
  using carrier_uint = detail::uint128_t;
  static const int exponent_bits = 15;
  static const int kappa = 2;
  static const int big_divisor = 1000;
  static const int small_divisor = 100;
  static const int min_k = -292;
  static const int max_k = 326;
  static const int shorter_interval_tie_lower_threshold = -77;
  static const int shorter_interval_tie_upper_threshold = -77;
};
} // namespace dragonbox
} // namespace detail
} // namespace fmt
#endif
