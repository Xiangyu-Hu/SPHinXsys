#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace ophelie_client_nc {

// The PCHIP curves use the normalized coordinate of the base profile.
struct Curve {
    const double *knots;
    const std::array<double, 4> *coefficients;
    int segment_count;

    double evaluate(double depth_fraction) const {
        const auto *first = knots;
        const auto *last = knots + segment_count;
        const int index = std::clamp(
            static_cast<int>(std::upper_bound(first, last, depth_fraction) - first) - 1,
            0, segment_count - 1);
        const double local = depth_fraction - knots[index];
        const auto &c = coefficients[index];
        return ((c[0] * local + c[1]) * local + c[2]) * local + c[3];
    }
};

inline constexpr double axis_knots[] = {0, 0.1111111111111111, 0.3888888888888889, 0.66666666666666663, 1};
inline constexpr std::array<double, 4> axis_coefficients[] = {
    {-6130.7513955497134, 2369.9751709340953, 60.217714285713939, 790},
    {-18784.24996043081, 7179.1056224615622, 359.81399688033832, 817.53999999999996},
    {3183.9082456731326, -2026.9725571314257, 0, 1068.8199999999999},
    {150.77066555583755, -418.95668612813375, -389.08006746349275, 980.65999999999997},
};
inline constexpr Curve axis_curve{axis_knots, axis_coefficients, 4};

inline constexpr double amp_knots[] = {0, 0.1111111111111111, 0.22222222222222221, 0.30555555555555558, 0.3888888888888889, 0.5, 0.66666666666666663, 0.83333333333333337, 1};
inline constexpr std::array<double, 4> amp_coefficients[] = {
    {-32805, -1215, 1485, 455},
    {975.095541401275, -1323.3439490445862, 0, 560},
    {214016.2938823874, -32739.149755591723, -257.96178343949043, 545},
    {-31794.755686271008, 4759.3304157163775, -1255.8139534883719, 420},
    {9553.7622523471528, -1547.6824901530892, -1124.9829609996216, 330},
    {28053.379441244939, -4285.1470989716754, -1115.0693568726358, 199},
    {-906.12371134020543, 701.04123711340139, -205.67010309278339, 24},
    {450.00000000000097, 29.999999999999801, -47.499999999999993, 5},
};
inline constexpr Curve amp_curve{amp_knots, amp_coefficients, 8};

inline constexpr double hot_r_knots[] = {0, 0.1111111111111111, 0.22222222222222221, 0.30555555555555558, 0.3888888888888889, 0.5, 0.66666666666666663, 0.83333333333333337, 1};
inline constexpr std::array<double, 4> hot_r_coefficients[] = {
    {5.8319999999999785, -2.2679999999999945, -0.09000000000000058, 0.76000000000000001},
    {16.836649681528584, -4.1387388535031704, -0.37800000000000017, 0.72999999999999998},
    {-1.3715629593335614, -0.43602155805978643, -0.67414012738853502, 0.66000000000000003},
    {18.198733031674141, -2.2919457013574576, -0.77538461538461523, 0.59999999999999998},
    {44.094245524296703, -4.3752429667519213, -0.77823529411764714, 0.53000000000000003},
    {0.093913043478260669, 0.32869565217391356, -0.11739130434782623, 0.45000000000000001},
    {0, 0, 0, 0.44},
    {0, 0, 0, 0.44},
};
inline constexpr Curve hot_r_curve{hot_r_knots, hot_r_coefficients, 8};

inline constexpr double wall_knots[] = {0, 0.1111111111111111, 0.22222222222222221, 0.30555555555555558, 0.3888888888888889, 0.5, 0.66666666666666663, 0.83333333333333337, 1};
inline constexpr std::array<double, 4> wall_coefficients[] = {
    {-260.35714285714351, -173.5714285714285, 202.5, 805},
    {890.15961691939071, -272.47805267358302, 154.28571428571428, 825},
    {-4794.6368715083627, 319.10614525139511, 126.70391061452511, 840},
    {-5760.0000000000064, 240.00000000000057, 80, 850},
    {2535.6521739130426, -686.73913043478251, 0, 855},
    {664.0993788819892, -298.50931677018662, -58.695652173913054, 850},
    {3209.142857142851, -637.71428571428476, -102.85714285714283, 835},
    {432.00000000000097, 35.999999999999801, -47.999999999999993, 815},
};
inline constexpr Curve wall_curve{wall_knots, wall_coefficients, 8};

// Arguments are in metres. Use only on initial liquid particles.
// The measured depths remain fixed at 20, 70 and 120 mm below the surface.
// This profile is intended for a melt deeper than 0.150 m.
inline double temperature_K(double radius_m, double depth_m,
                            double melt_radius_m = 0.1415, double melt_height_m = 0.230) {
    const double r = std::clamp(radius_m / melt_radius_m, 0.0, 1.0);
    constexpr double join_depth_m = 0.150;
    constexpr double base_tail_m = 0.030;
    const double depth = std::clamp(depth_m, 0.0, melt_height_m);
    double base_depth = depth;
    if (depth > join_depth_m) {
        const double tail_m = melt_height_m - join_depth_m;
        const double tail_fraction = (depth - join_depth_m) / tail_m;
        // Match the slope at the join when stretching the lower melt.
        base_depth = join_depth_m + base_tail_m *
            (1.0 - std::pow(1.0 - tail_fraction, tail_m / base_tail_m));
    }
    const double d = base_depth / (join_depth_m + base_tail_m);
    const double s = r * r;
    const double rc = hot_r_curve.evaluate(d);
    const double sc = rc * rc;
    const double a = 1.5;
    const double b = a * (1.0 - sc) / sc;
    const double ring = std::pow(s / sc, a) * std::pow((1.0 - s) / (1.0 - sc), b);
    const double centre = axis_curve.evaluate(d);
    return 273.15 + centre + amp_curve.evaluate(d) * ring +
           (wall_curve.evaluate(d) - centre) * std::pow(s, 6.0);
}

} // namespace ophelie_client_nc
