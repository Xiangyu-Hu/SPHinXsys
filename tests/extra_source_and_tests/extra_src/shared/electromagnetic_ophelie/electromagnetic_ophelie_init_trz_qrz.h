/**
 * @file electromagnetic_ophelie_init_trz_qrz.h
 * @brief Load axisymmetric T(r,z) / Q(r,z) CSV seeds and apply to glass particles.
 *
 * Used by Doc248 client stirring (C2 / C2r) for ASSUMED French Fig.8a / Fig.10 remapped ICs.
 * Not literature PASS.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_INIT_TRZ_QRZ_H
#define ELECTROMAGNETIC_OPHELIE_INIT_TRZ_QRZ_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_doc247_q_shape_correct.h"
#include "electromagnetic_ophelie_joule_to_heat_one_way.h"

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

struct Doc248RzTemperatureGrid
{
    std::vector<double> r_axis;
    std::vector<double> z_axis;
    std::vector<double> T;
    double T_min = 0;
    double T_max = 0;
    size_t n_nodes = 0;
    bool ok = false;

    size_t nr() const { return r_axis.size(); }
    size_t nz() const { return z_axis.size(); }
    double at(size_t ir, size_t iz) const { return T[ir + nr() * iz]; }

    double sample(double r, double z) const
    {
        if (!ok || nr() < 2 || nz() < 2)
            return T.empty() ? 0.0 : T.front();
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
        const double t00 = at(ir0, iz0), t10 = at(ir1, iz0), t01 = at(ir0, iz1), t11 = at(ir1, iz1);
        return (1.0 - tr) * (1.0 - tz) * t00 + tr * (1.0 - tz) * t10 + (1.0 - tr) * tz * t01 + tr * tz * t11;
    }
};

inline bool hostLoadDoc248RzTemperatureGrid(const std::string &csv_path, Doc248RzTemperatureGrid &grid,
                                            std::string &err)
{
    grid = Doc248RzTemperatureGrid{};
    std::ifstream in(csv_path);
    if (!in)
    {
        err = "cannot open init-T-rz-csv";
        return false;
    }
    std::string line;
    if (!std::getline(in, line))
    {
        err = "empty init-T-rz-csv";
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
    int i_r = -1, i_z = -1, i_T = -1;
    for (int i = 0; i < static_cast<int>(hdr.size()); ++i)
    {
        const std::string &h = hdr[static_cast<size_t>(i)];
        if (h == "r_m" || h == "r")
            i_r = i;
        else if (h == "z_m" || h == "z")
            i_z = i;
        else if (h == "t_k" || h == "t" || h == "temperature_k")
            i_T = i;
    }
    if (i_r < 0 || i_z < 0 || i_T < 0)
    {
        err = "init-T-rz-csv needs r_m,z_m,T_K";
        return false;
    }
    std::vector<double> rs, zs, Ts;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        std::stringstream ss(line);
        std::vector<std::string> cells;
        std::string cell;
        while (std::getline(ss, cell, ','))
            cells.push_back(strip(cell));
        if (static_cast<int>(cells.size()) <= std::max(i_r, std::max(i_z, i_T)))
            continue;
        const double r = std::atof(cells[static_cast<size_t>(i_r)].c_str());
        const double z = std::atof(cells[static_cast<size_t>(i_z)].c_str());
        const double Tv = std::atof(cells[static_cast<size_t>(i_T)].c_str());
        if (!(std::isfinite(r) && std::isfinite(z) && std::isfinite(Tv)) || r < 0.0 || !(Tv > 0.0))
        {
            err = "non-finite/non-positive T or bad r in init-T-rz-csv";
            return false;
        }
        rs.push_back(r);
        zs.push_back(z);
        Ts.push_back(Tv);
    }
    if (rs.empty())
    {
        err = "no usable rows in init-T-rz-csv";
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
        err = "init-T-rz-csv is not a complete regular r-z grid";
        return false;
    }
    grid.T.assign(grid.nr() * grid.nz(), std::numeric_limits<double>::quiet_NaN());
    grid.T_min = Ts.front();
    grid.T_max = Ts.front();
    for (size_t k = 0; k < rs.size(); ++k)
    {
        const auto ir = static_cast<size_t>(
            std::lower_bound(grid.r_axis.begin(), grid.r_axis.end(), rs[k]) - grid.r_axis.begin());
        const auto iz = static_cast<size_t>(
            std::lower_bound(grid.z_axis.begin(), grid.z_axis.end(), zs[k]) - grid.z_axis.begin());
        if (ir >= grid.nr() || iz >= grid.nz() || std::fabs(grid.r_axis[ir] - rs[k]) > 1e-12 ||
            std::fabs(grid.z_axis[iz] - zs[k]) > 1e-12)
        {
            err = "init-T-rz-csv grid axis mismatch";
            return false;
        }
        grid.T[ir + grid.nr() * iz] = Ts[k];
        grid.T_min = std::min(grid.T_min, Ts[k]);
        grid.T_max = std::max(grid.T_max, Ts[k]);
    }
    for (double v : grid.T)
    {
        if (!std::isfinite(v))
        {
            err = "init-T-rz-csv missing grid node";
            return false;
        }
    }
    grid.n_nodes = rs.size();
    grid.ok = true;
    return true;
}

/**
 * Sample absolute (r,z) seed (z in solver coords). Sets Temperature and thermal_deltaT = T - t_ref.
 * Radial origin = crucible axis at `center` (xy).
 */
inline void hostApplyInitTRzAbsolute(BaseParticles &particles, const Vecd &center, Real t_thermal_ref,
                                     const Doc248RzTemperatureGrid &grid, const char *log_tag = "init-T-rz")
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(particles, kOphelieThermalDeltaTField);
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *temp = particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const size_t n = particles.TotalRealParticles();
    double Tmin = 1.0e300, Tmax = -1.0e300, Tsum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const double dx = static_cast<double>(pos[i][0] - center[0]);
        const double dy = static_cast<double>(pos[i][1] - center[1]);
        const double r = std::sqrt(dx * dx + dy * dy);
        const double z = static_cast<double>(pos[i][2]);
        const double Ti = grid.sample(r, z);
        temp[i] = static_cast<Real>(Ti);
        delta_t[i] = static_cast<Real>(Ti) - t_thermal_ref;
        Tmin = std::min(Tmin, Ti);
        Tmax = std::max(Tmax, Ti);
        Tsum += Ti;
    }
    syncVariableToDevice<Real>(particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(particles, kOphelieThermalDeltaTField);
    std::cout << "[ophelie][" << log_tag << "] n_nodes=" << grid.n_nodes << " nr=" << grid.nr()
              << " nz=" << grid.nz() << " T_csv=[" << grid.T_min << "," << grid.T_max << "] T_particle=["
              << Tmin << "," << Tmax << "] T_mean=" << (n > 0 ? Tsum / static_cast<double>(n) : 0.0)
              << " T_thermal_ref=" << t_thermal_ref << "\n";
}

/**
 * Assign prescribed shaped volumetric Q from rz CSV; rescale ∫Q dV → P_target.
 * Physical (r_m,z_m) seeds preferred; unit-box seeds use live R,H,z_bottom.
 * Returns realized absorbed power after scale.
 */
inline Real hostAssignPrescribedShapeJouleHeat(BaseParticles &particles, const std::string &q_field,
                                               const Doc247RzQGrid &shape, const Vecd &center,
                                               Real glass_radius, Real glass_half_height, Real P_target,
                                               const char *log_tag = "init-Q-rz")
{
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

    Real P_raw = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
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
        q[i] = static_cast<Real>(q_raw);
        P_raw += q[i] * vol[i];
    }
    const Real scale = (P_raw > TinyReal && P_target > TinyReal) ? (P_target / P_raw) : Real(0);
    Real P_after = Real(0);
    Real q_max = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        q[i] *= scale;
        P_after += q[i] * vol[i];
        q_max = std::max(q_max, q[i]);
    }
    syncVariableToDevice<Real>(particles, q_field);
    std::cout << "[ophelie][" << log_tag << "] n=" << n << " unit_box=" << (unit ? 1 : 0)
              << " nodes=" << shape.n_nodes << " Qshape=[" << shape.Q_min << "," << shape.Q_max << "]"
              << " P_raw=" << P_raw << " scale=" << scale << " P_abs=" << P_after << " Q_max=" << q_max
              << "\n";
    return P_after;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_INIT_TRZ_QRZ_H
