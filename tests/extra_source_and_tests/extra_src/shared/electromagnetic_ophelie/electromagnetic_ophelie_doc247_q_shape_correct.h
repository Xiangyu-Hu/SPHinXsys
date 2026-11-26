/**
 * @file electromagnetic_ophelie_doc247_q_shape_correct.h
 * @brief Doc247 ASSUMED Q-shape correct: load r–z seed, blend live EM Q, rescale ∫Q dV → P_target.
 *
 * Used by French NC (--fig6-q-correct) and stirring (--fig10-q-correct). Engineering mix only —
 * not literature PASS / FRENCH_VALIDATED.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC247_Q_SHAPE_CORRECT_H
#define ELECTROMAGNETIC_OPHELIE_DOC247_Q_SHAPE_CORRECT_H

#include "base_data_type.h"
#include "base_particles.h"
#include "sphinxsys_variable.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct Doc247RzQGrid
{
    std::vector<double> r_axis;
    std::vector<double> z_axis;
    std::vector<double> Q;
    double Q_min = 0;
    double Q_max = 0;
    size_t n_nodes = 0;
    bool ok = false;
    bool looks_normalized_unit_box = false;

    size_t nr() const { return r_axis.size(); }
    size_t nz() const { return z_axis.size(); }
    double at(size_t ir, size_t iz) const { return Q[ir + nr() * iz]; }

    double sample(double r, double z) const
    {
        if (!ok || nr() < 2 || nz() < 2)
            return Q.empty() ? 0.0 : Q.front();
        r = std::max(r_axis.front(), std::min(r_axis.back(), r));
        z = std::max(z_axis.front(), std::min(z_axis.back(), z));
        auto upper = [](const std::vector<double> &ax, double v) {
            return static_cast<size_t>(
                std::max<std::ptrdiff_t>(1, std::lower_bound(ax.begin(), ax.end(), v) - ax.begin()));
        };
        const size_t ir1 = std::min(upper(r_axis, r), nr() - 1);
        const size_t iz1 = std::min(upper(z_axis, z), nz() - 1);
        const size_t ir0 = ir1 - 1;
        const size_t iz0 = iz1 - 1;
        const double r0 = r_axis[ir0], r1 = r_axis[ir1];
        const double z0 = z_axis[iz0], z1 = z_axis[iz1];
        const double tr = (r1 > r0) ? (r - r0) / (r1 - r0) : 0.0;
        const double tz = (z1 > z0) ? (z - z0) / (z1 - z0) : 0.0;
        const double q00 = at(ir0, iz0), q10 = at(ir1, iz0), q01 = at(ir0, iz1), q11 = at(ir1, iz1);
        return (1.0 - tr) * (1.0 - tz) * q00 + tr * (1.0 - tz) * q10 + (1.0 - tr) * tz * q01 + tr * tz * q11;
    }
};

inline bool hostLoadDoc247RzQGrid(const std::string &csv_path, Doc247RzQGrid &grid, std::string &err)
{
    grid = Doc247RzQGrid{};
    std::ifstream in(csv_path);
    if (!in)
    {
        err = "cannot open init-Q-rz-csv";
        return false;
    }
    std::string line;
    if (!std::getline(in, line))
    {
        err = "empty init-Q-rz-csv";
        return false;
    }
    auto strip = [](std::string s) {
        while (!s.empty() && (s.back() == '\r' || s.back() == ' ' || s.back() == '\t'))
            s.pop_back();
        size_t i = 0;
        while (i < s.size() && (s[i] == ' ' || s[i] == '\t'))
            ++i;
        return s.substr(i);
    };
    auto lower = [&](std::string s) {
        s = strip(s);
        for (char &c : s)
            c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    };
    std::vector<std::string> hdr;
    {
        std::stringstream ss(line);
        std::string cell;
        while (std::getline(ss, cell, ','))
            hdr.push_back(lower(cell));
    }
    int i_r = -1, i_z = -1, i_Q = -1;
    auto pick_r = [](const std::string &h) -> int {
        if (h == "r_m" || h == "r")
            return 2;
        if (h == "r_over_r" || h == "x_unit")
            return 1;
        return 0;
    };
    auto pick_z = [](const std::string &h) -> int {
        if (h == "z_m" || h == "z")
            return 2;
        if (h == "z_over_h" || h == "z_or_y_unit")
            return 1;
        return 0;
    };
    int pri_r = 0, pri_z = 0;
    for (int i = 0; i < static_cast<int>(hdr.size()); ++i)
    {
        const std::string &h = hdr[static_cast<size_t>(i)];
        const int pr = pick_r(h);
        if (pr > pri_r)
        {
            pri_r = pr;
            i_r = i;
        }
        const int pz = pick_z(h);
        if (pz > pri_z)
        {
            pri_z = pz;
            i_z = i;
        }
        if (i_Q < 0 && (h == "q_norm" || h == "q_wpm3" || h == "q_shape" || h == "q" ||
                        h == "value_normalized" || h == "value"))
            i_Q = i;
    }
    if (i_r < 0 || i_z < 0 || i_Q < 0)
    {
        err = "init-Q-rz-csv needs r_m,z_m,Q_norm (or Q_Wpm3 / value_normalized)";
        return false;
    }
    std::vector<double> rs, zs, Qs;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::vector<std::string> cells;
        std::string cell;
        while (std::getline(ss, cell, ','))
            cells.push_back(strip(cell));
        if (static_cast<int>(cells.size()) <= std::max(i_r, std::max(i_z, i_Q)))
            continue;
        const double r = std::atof(cells[static_cast<size_t>(i_r)].c_str());
        const double z = std::atof(cells[static_cast<size_t>(i_z)].c_str());
        const double Qv = std::atof(cells[static_cast<size_t>(i_Q)].c_str());
        if (!(std::isfinite(r) && std::isfinite(z) && std::isfinite(Qv)) || r < 0.0 || Qv < 0.0)
        {
            err = "non-finite/negative Q or bad r in init-Q-rz-csv";
            return false;
        }
        rs.push_back(r);
        zs.push_back(z);
        Qs.push_back(Qv);
    }
    if (rs.empty())
    {
        err = "no usable rows in init-Q-rz-csv";
        return false;
    }
    auto uniq_sorted = [](std::vector<double> v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end(),
                            [](double a, double b) { return std::fabs(a - b) < 1e-12; }),
                v.end());
        return v;
    };
    grid.r_axis = uniq_sorted(rs);
    grid.z_axis = uniq_sorted(zs);
    if (grid.nr() < 2 || grid.nz() < 2 || grid.nr() * grid.nz() != rs.size())
    {
        err = "init-Q-rz-csv is not a complete regular r-z grid";
        return false;
    }
    grid.Q.assign(grid.nr() * grid.nz(), std::numeric_limits<double>::quiet_NaN());
    grid.Q_min = Qs.front();
    grid.Q_max = Qs.front();
    for (size_t k = 0; k < rs.size(); ++k)
    {
        const auto ir = static_cast<size_t>(
            std::lower_bound(grid.r_axis.begin(), grid.r_axis.end(), rs[k]) - grid.r_axis.begin());
        const auto iz = static_cast<size_t>(
            std::lower_bound(grid.z_axis.begin(), grid.z_axis.end(), zs[k]) - grid.z_axis.begin());
        if (ir >= grid.nr() || iz >= grid.nz() || std::fabs(grid.r_axis[ir] - rs[k]) > 1e-12 ||
            std::fabs(grid.z_axis[iz] - zs[k]) > 1e-12)
        {
            err = "init-Q-rz-csv grid axis mismatch";
            return false;
        }
        grid.Q[ir + grid.nr() * iz] = Qs[k];
        grid.Q_min = std::min(grid.Q_min, Qs[k]);
        grid.Q_max = std::max(grid.Q_max, Qs[k]);
    }
    for (double v : grid.Q)
    {
        if (!std::isfinite(v))
        {
            err = "init-Q-rz-csv missing grid node";
            return false;
        }
    }
    grid.n_nodes = rs.size();
    grid.looks_normalized_unit_box =
        grid.r_axis.front() >= -1e-12 && grid.z_axis.front() >= -1e-12 && grid.r_axis.back() > 0.95 &&
        grid.r_axis.back() <= 1.0 + 1e-6 && grid.z_axis.back() > 0.95 && grid.z_axis.back() <= 1.0 + 1e-6;
    grid.ok = true;
    return true;
}

/**
 * Blend live EM particle Q with an r–z shape seed; rescale ∫Q dV → P_target.
 * Returns realized power after correction. Log tag is caller-chosen (Fig6 / Fig10).
 */
inline Real hostBlendEmQWithRzShape(BaseParticles &particles, const std::string &q_field,
                                    const Doc247RzQGrid &shape, const Vecd &center, Real glass_radius,
                                    Real glass_half_height, Real P_target, Real blend_shape, bool verbose,
                                    const char *log_tag = "Q-shape-correct")
{
    const Real b = std::max(Real(0), std::min(Real(1), blend_shape));
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, q_field);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *q = particles.getVariableDataByName<Real>(q_field);
    const size_t n = particles.TotalRealParticles();
    const Real R = std::max(glass_radius, TinyReal);
    const Real H = std::max(Real(2) * glass_half_height, TinyReal);
    const Real z_bottom = center[2] - glass_half_height;
    const bool unit = shape.looks_normalized_unit_box;

    StdVec<Real> q_fig(n, Real(0));
    Real P_em = Real(0);
    Real P_fig_raw = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        P_em += q[i] * vol[i];
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        const double z = static_cast<double>(pos[i][2]);
        double q_raw = 0.0;
        if (unit)
        {
            const double rn = std::max(0.0, std::min(1.0, r / static_cast<double>(R)));
            const double zn =
                std::max(0.0, std::min(1.0, (z - static_cast<double>(z_bottom)) / static_cast<double>(H)));
            q_raw = shape.sample(rn, zn);
        }
        else
        {
            q_raw = shape.sample(r, z);
        }
        q_fig[i] = static_cast<Real>(q_raw);
        P_fig_raw += q_fig[i] * vol[i];
    }
    const Real scale_em = (P_em > TinyReal && P_target > TinyReal) ? (P_target / P_em) : Real(0);
    const Real scale_fig = (P_fig_raw > TinyReal && P_target > TinyReal) ? (P_target / P_fig_raw) : Real(0);
    Real P_mix = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        q[i] = (Real(1) - b) * (q[i] * scale_em) + b * (q_fig[i] * scale_fig);
        P_mix += q[i] * vol[i];
    }
    const Real scale_out = (P_mix > TinyReal && P_target > TinyReal) ? (P_target / P_mix) : Real(0);
    Real P_out = Real(0);
    Real q_max = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        q[i] *= scale_out;
        P_out += q[i] * vol[i];
        q_max = std::max(q_max, q[i]);
    }
    syncVariableToDevice<Real>(particles, q_field);
    if (verbose)
    {
        std::cout << "[ophelie][doc247][" << log_tag << "] blend=" << b << " unit_box=" << (unit ? 1 : 0)
                  << " P_em=" << P_em << " P_fig_raw=" << P_fig_raw << " P_out=" << P_out << " Q_max=" << q_max
                  << "\n";
    }
    return P_out;
}

/** Backward-compatible alias used by NC frozen_q. */
inline Real hostBlendEmQWithFig6Shape(BaseParticles &particles, const std::string &q_field,
                                      const Doc247RzQGrid &shape, const Vecd &center, Real glass_radius,
                                      Real glass_half_height, Real P_target, Real blend_fig6, bool verbose)
{
    return hostBlendEmQWithRzShape(particles, q_field, shape, center, glass_radius, glass_half_height, P_target,
                                   blend_fig6, verbose, "Fig6Q-correct");
}

inline bool hostDoc247QShapeCorrectDue(int update_count_after, int every)
{
    if (every <= 1)
        return true;
    return update_count_after > 0 && (update_count_after % every) == 0;
}

inline bool hostFig6QCorrectDue(int update_count_after, int every)
{
    return hostDoc247QShapeCorrectDue(update_count_after, every);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_DOC247_Q_SHAPE_CORRECT_H
