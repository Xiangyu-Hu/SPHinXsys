#ifndef ELECTROMAGNETIC_OPHELIE_THERMO_EM_COUPLING_H
#define ELECTROMAGNETIC_OPHELIE_THERMO_EM_COUPLING_H

/**
 * Periodic thermo–EM coupling helpers (Jacoutot CES 2008 Fig. 4 dataflow).
 *
 * EM may still be solved in 3D; σ and Q are reduced to axisymmetric (r,z) so the
 * coupling interface matches the paper even before a 2D OPHELIE solver exists.
 *
 * Q remapping never silently rescales. Conservative remap is opt-in and logged.
 */

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_literature_parameters.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct FrenchCylindricalFrame
{
    Vecd origin = Vecd::Zero();
    int axis = 2; // 0/1/2 = x/y/z
    Real radius = 0.25;
    Real z_min = 0.0;
    Real z_max = 0.185;
};

inline void frenchCylindricalRZ(const Vecd &pos, const FrenchCylindricalFrame &frame, Real &r, Real &z)
{
    z = pos[frame.axis];
    Real r2 = 0.0;
    for (int d = 0; d < 3; ++d)
    {
        if (d == frame.axis)
        {
            continue;
        }
        const Real dr = pos[d] - frame.origin[d];
        r2 += dr * dr;
    }
    r = std::sqrt(r2);
}

inline FrenchCylindricalFrame makeFrenchCylinderFrameFromGlass(const Vecd &center, Real radius, Real half_height,
                                                              FrenchVerticalAxis axis)
{
    FrenchCylindricalFrame frame;
    frame.origin = center;
    frame.axis = frenchVerticalAxisIndex(axis);
    frame.radius = radius;
    frame.z_min = center[frame.axis] - half_height;
    frame.z_max = center[frame.axis] + half_height;
    return frame;
}

struct AxisymmetricRZGrid
{
    FrenchCylindricalFrame frame;
    size_t n_r = 32;
    size_t n_z = 24;
    StdVec<Real> value;
    StdVec<Real> weight;
    StdVec<size_t> count;
    StdVec<uint8_t> valid;

    size_t index(size_t ir, size_t iz) const { return ir + n_r * iz; }
    Real dr() const { return frame.radius / static_cast<Real>(std::max<size_t>(n_r, 1)); }
    Real dz() const
    {
        return (frame.z_max - frame.z_min) / static_cast<Real>(std::max<size_t>(n_z, 1));
    }

    void allocate()
    {
        const size_t n = n_r * n_z;
        value.assign(n, Real(0));
        weight.assign(n, Real(0));
        count.assign(n, 0);
        valid.assign(n, 0);
    }

    bool binIndex(Real r, Real z, size_t &ir, size_t &iz) const
    {
        if (r < 0.0 || r > frame.radius * Real(1.000001) || z < frame.z_min || z > frame.z_max)
        {
            if (r > frame.radius && r <= frame.radius * Real(1.05))
            {
                r = frame.radius;
            }
            else if (r < 0.0 || r > frame.radius)
            {
                return false;
            }
            if (z < frame.z_min || z > frame.z_max)
            {
                return false;
            }
        }
        ir = static_cast<size_t>(std::min(static_cast<Real>(n_r - 1),
                                          std::floor(r / (dr() + TinyReal))));
        iz = static_cast<size_t>(std::min(static_cast<Real>(n_z - 1),
                                          std::floor((z - frame.z_min) / (dz() + TinyReal))));
        return true;
    }

    Real sample(Real r, Real z) const
    {
        size_t ir = 0;
        size_t iz = 0;
        if (!binIndex(r, z, ir, iz) || !valid[index(ir, iz)])
        {
            return Real(0);
        }
        return value[index(ir, iz)];
    }
};

inline void azimuthalDeposit(AxisymmetricRZGrid &grid, const Vecd *pos, const Real *field, const Real *vol, size_t n)
{
    grid.allocate();
    for (size_t i = 0; i < n; ++i)
    {
        Real r = 0.0;
        Real z = 0.0;
        frenchCylindricalRZ(pos[i], grid.frame, r, z);
        size_t ir = 0;
        size_t iz = 0;
        if (!grid.binIndex(r, z, ir, iz))
        {
            continue;
        }
        if (!std::isfinite(field[i]) || vol[i] <= TinyReal)
        {
            continue;
        }
        const size_t id = grid.index(ir, iz);
        grid.value[id] += field[i] * vol[i];
        grid.weight[id] += vol[i];
        grid.count[id] += 1;
    }
    for (size_t id = 0; id < grid.value.size(); ++id)
    {
        if (grid.weight[id] > TinyReal && grid.count[id] > 0)
        {
            grid.value[id] /= grid.weight[id];
            grid.valid[id] = 1;
        }
        else
        {
            grid.value[id] = Real(0);
            grid.valid[id] = 0;
        }
    }
}

/** Empty bins inherit a neighbor / global mean. Never silently write σ=0. */
inline void fillEmptyAxisymmetricBins(AxisymmetricRZGrid &grid, Real fallback)
{
    if (!std::isfinite(fallback) || fallback <= TinyReal)
    {
        fallback = Real(1);
    }
    bool any_valid = false;
    Real mean = 0.0;
    size_t n_valid = 0;
    for (size_t id = 0; id < grid.value.size(); ++id)
    {
        if (grid.valid[id])
        {
            any_valid = true;
            mean += grid.value[id];
            ++n_valid;
        }
    }
    if (any_valid)
    {
        fallback = mean / static_cast<Real>(n_valid);
    }
    for (size_t pass = 0; pass < 8; ++pass)
    {
        bool changed = false;
        for (size_t iz = 0; iz < grid.n_z; ++iz)
        {
            for (size_t ir = 0; ir < grid.n_r; ++ir)
            {
                const size_t id = grid.index(ir, iz);
                if (grid.valid[id])
                {
                    continue;
                }
                Real acc = 0.0;
                size_t n_n = 0;
                const int nbr[4][2] = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
                for (const auto &d : nbr)
                {
                    const int jr = static_cast<int>(ir) + d[0];
                    const int jz = static_cast<int>(iz) + d[1];
                    if (jr < 0 || jz < 0 || jr >= static_cast<int>(grid.n_r) || jz >= static_cast<int>(grid.n_z))
                    {
                        continue;
                    }
                    const size_t nid = grid.index(static_cast<size_t>(jr), static_cast<size_t>(jz));
                    if (grid.valid[nid])
                    {
                        acc += grid.value[nid];
                        ++n_n;
                    }
                }
                if (n_n > 0)
                {
                    grid.value[id] = acc / static_cast<Real>(n_n);
                    grid.valid[id] = 1;
                    changed = true;
                }
            }
        }
        if (!changed)
        {
            break;
        }
    }
    for (size_t id = 0; id < grid.value.size(); ++id)
    {
        if (!grid.valid[id])
        {
            grid.value[id] = fallback;
            grid.valid[id] = 1;
        }
    }
}

inline Real integrateAxisymmetricPower(const AxisymmetricRZGrid &grid)
{
    Real power = 0.0;
    const Real pi = std::acos(Real(-1));
    const Real dr = grid.dr();
    const Real dz = grid.dz();
    for (size_t iz = 0; iz < grid.n_z; ++iz)
    {
        for (size_t ir = 0; ir < grid.n_r; ++ir)
        {
            const Real r0 = static_cast<Real>(ir) * dr;
            const Real r1 = r0 + dr;
            const Real area = pi * (r1 * r1 - r0 * r0);
            power += grid.value[grid.index(ir, iz)] * area * dz;
        }
    }
    return power;
}

inline void mapAxisymmetricToParticles(const AxisymmetricRZGrid &grid, const Vecd *pos, StdVec<Real> &out, size_t n)
{
    out.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        Real r = 0.0;
        Real z = 0.0;
        frenchCylindricalRZ(pos[i], grid.frame, r, z);
        out[i] = grid.sample(r, z);
    }
}

inline Real underRelaxScalar(Real old_v, Real new_v, Real alpha)
{
    alpha = std::min(Real(1), std::max(Real(0), alpha));
    return (Real(1) - alpha) * old_v + alpha * new_v;
}

struct FrenchPowerMapBook
{
    Real p_em_recon = 0.0;
    Real p_axisymmetric_grid = 0.0;
    Real p_euler_sample = 0.0;
    Real p_sph_particles = 0.0;
    Real remap_scale = 1.0;
    Real p_euler_before_remap = 0.0;
    Real p_euler_after_remap = 0.0;
    bool conservative_remap_applied = false;

    Real relErr(Real a, Real b) const { return std::abs(a - b) / (std::abs(b) + TinyReal); }
};

inline void logFrenchPowerMapBook(const FrenchPowerMapBook &book)
{
    std::cout << "[ophelie][q-map] P_em_recon=" << book.p_em_recon
              << " P_axisymmetric=" << book.p_axisymmetric_grid
              << " P_euler_sample=" << book.p_euler_sample << " P_sph=" << book.p_sph_particles
              << " rel_axi=" << book.relErr(book.p_axisymmetric_grid, book.p_em_recon)
              << " rel_euler=" << book.relErr(book.p_euler_sample, book.p_em_recon)
              << " rel_sph=" << book.relErr(book.p_sph_particles, book.p_em_recon)
              << " conservative_remap=" << (book.conservative_remap_applied ? 1 : 0);
    if (book.conservative_remap_applied)
    {
        std::cout << " scale=" << book.remap_scale << " P_before=" << book.p_euler_before_remap
                  << " P_after=" << book.p_euler_after_remap;
    }
    std::cout << std::endl;
}

struct ThermoEMCouplingState
{
    Real last_update_time = -1.0e300;
    size_t update_count = 0;
    Real calibrated_current_per_loop = 0.0;
    Real generator_power_w = 400000.0;
    Real target_glass_absorbed_power_w = 50000.0;
    Real reconstructed_glass_power_w = 0.0;
    Real coil_current_per_loop = 0.0;
    bool current_frozen = false;
    StdVec<Real> sigma_em;
    AxisymmetricRZGrid sigma_rz;
    AxisymmetricRZGrid q_rz;
    FrenchPowerMapBook last_power;
    Real last_phi_eq_res = 0.0;
    Real last_sigma_min = 0.0;
    Real last_sigma_max = 0.0;
    Real last_sigma_mean = 0.0;
};

inline bool shouldUpdateEM(Real physical_time, Real interval_s, Real last_update_time)
{
    if (interval_s <= TinyReal)
    {
        return true;
    }
    if (last_update_time < Real(-1.0e100))
    {
        return true;
    }
    return physical_time + TinyReal >= last_update_time + interval_s;
}

inline bool evaluateSigmaFromTemperature(const OphelieTemperatureLaw &sigma_law, const Real *temperature, size_t n,
                                         StdVec<Real> &sigma_raw, std::string &error)
{
    sigma_raw.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        const Real t = temperature[i];
        if (!std::isfinite(t) || t <= TinyReal)
        {
            error = "invalid temperature at particle " + std::to_string(i) + " T=" + std::to_string(t);
            return false;
        }
        const Real s = sigma_law.evaluate(t);
        if (!std::isfinite(s) || s <= TinyReal)
        {
            error = "invalid sigma(T) at particle " + std::to_string(i) + " T=" + std::to_string(t) +
                    " sigma=" + std::to_string(s);
            return false;
        }
        sigma_raw[i] = s;
    }
    return true;
}

inline void sigmaFieldStats(const StdVec<Real> &sigma, const Real *vol, size_t n, Real &smin, Real &smax, Real &smean)
{
    smin = std::numeric_limits<Real>::max();
    smax = 0.0;
    Real acc = 0.0;
    Real w = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        smin = std::min(smin, sigma[i]);
        smax = std::max(smax, sigma[i]);
        acc += sigma[i] * vol[i];
        w += vol[i];
    }
    smean = acc / (w + TinyReal);
}

inline Real hostParticlePower(const Real *q, const Real *vol, size_t n)
{
    Real p = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        p += q[i] * vol[i];
    }
    return p;
}

inline void hostAssignSigmaField(BaseParticles &particles, const std::string &sigma_field, const StdVec<Real> &sigma)
{
    syncVariableToHost<Real>(particles, sigma_field);
    Real *s = particles.getVariableDataByName<Real>(sigma_field);
    const size_t n = particles.TotalRealParticles();
    if (sigma.size() != n)
    {
        std::cerr << "[ophelie][thermo-em] sigma size mismatch n=" << n << " sigma=" << sigma.size() << std::endl;
        return;
    }
    for (size_t i = 0; i < n; ++i)
    {
        s[i] = sigma[i];
    }
    syncVariableToDevice<Real>(particles, sigma_field);
}

/**
 * Host-side melt diagnostics so SSH runs do not need VTP / ParaView.
 * Outer / center split matches the Q skin diagnostic: r > 0.75 R vs r < 0.25 R.
 */
struct FrenchMeltSpatialStats
{
    Real t_mean = 0.0;
    Real t_std = 0.0;
    Real t_outer = 0.0;
    Real t_center = 0.0;
    Real t_outer_over_center = 0.0;
    Real t_bottom = 0.0;
    Real t_top = 0.0;
    Real u_rms = 0.0;
    Real u_r_rms = 0.0;
    Real u_theta_rms = 0.0;
    Real u_z_rms = 0.0;
    Real u_max = 0.0;
    Real volume = 0.0;
};

inline FrenchMeltSpatialStats computeFrenchMeltSpatialStats(const Vecd *pos, const Vecd *vel, const Real *temperature,
                                                            const Real *vol, size_t n,
                                                            const FrenchCylindricalFrame &frame)
{
    FrenchMeltSpatialStats s;
    const int ax = frame.axis;
    Vecd e_z = Vecd::Zero();
    e_z[ax] = Real(1);
    const Real r_outer = Real(0.75) * frame.radius;
    const Real r_center = Real(0.25) * frame.radius;
    const Real h = std::max(frame.z_max - frame.z_min, TinyReal);
    const Real z_bot = frame.z_min + Real(0.25) * h;
    const Real z_top = frame.z_min + Real(0.75) * h;

    double t_acc = 0.0;
    double t2_acc = 0.0;
    Real u2_acc = 0.0;
    Real ur2_acc = 0.0;
    Real uth2_acc = 0.0;
    Real uz2_acc = 0.0;
    double tout_acc = 0.0;
    double tout_w = 0.0;
    double tctr_acc = 0.0;
    double tctr_w = 0.0;
    double tbot_acc = 0.0;
    double tbot_w = 0.0;
    double ttop_acc = 0.0;
    double ttop_w = 0.0;
    double w = 0.0;

    for (size_t i = 0; i < n; ++i)
    {
        if (vol[i] <= TinyReal || !std::isfinite(temperature[i]))
        {
            continue;
        }
        const double vi = static_cast<double>(vol[i]);
        const double ti = static_cast<double>(temperature[i]);
        w += vi;
        t_acc += ti * vi;
        t2_acc += ti * ti * vi;

        Real r = 0.0;
        Real z = 0.0;
        frenchCylindricalRZ(pos[i], frame, r, z);
        if (r >= r_outer)
        {
            tout_acc += ti * vi;
            tout_w += vi;
        }
        if (r <= r_center)
        {
            tctr_acc += ti * vi;
            tctr_w += vi;
        }
        if (z <= z_bot)
        {
            tbot_acc += ti * vi;
            tbot_w += vi;
        }
        if (z >= z_top)
        {
            ttop_acc += ti * vi;
            ttop_w += vi;
        }

        if (vel == nullptr)
        {
            continue;
        }
        const Real un = vel[i].norm();
        if (!std::isfinite(un))
        {
            continue;
        }
        s.u_max = std::max(s.u_max, un);
        u2_acc += un * un * vi;
        const Real uz = vel[i].dot(e_z);
        uz2_acc += uz * uz * vi;
        Vecd e_r = pos[i] - frame.origin;
        e_r[ax] = Real(0);
        const Real rxy = e_r.norm();
        if (rxy > TinyReal)
        {
            e_r /= rxy;
            const Real ur = vel[i].dot(e_r);
            const Real uth = vel[i].dot(e_z.cross(e_r));
            ur2_acc += ur * ur * vi;
            uth2_acc += uth * uth * vi;
        }
    }

    const double w_safe = w + static_cast<double>(TinyReal);
    const double t_mean_d = t_acc / w_safe;
    s.volume = static_cast<Real>(w);
    s.t_mean = static_cast<Real>(t_mean_d);
    s.t_std = static_cast<Real>(std::sqrt(std::max(t2_acc / w_safe - t_mean_d * t_mean_d, 0.0)));
    s.t_outer = static_cast<Real>(tout_acc / (tout_w + static_cast<double>(TinyReal)));
    s.t_center = static_cast<Real>(tctr_acc / (tctr_w + static_cast<double>(TinyReal)));
    s.t_outer_over_center = s.t_outer / (s.t_center + TinyReal);
    s.t_bottom = static_cast<Real>(tbot_acc / (tbot_w + static_cast<double>(TinyReal)));
    s.t_top = static_cast<Real>(ttop_acc / (ttop_w + static_cast<double>(TinyReal)));
    s.u_rms = std::sqrt(u2_acc / static_cast<Real>(w_safe));
    s.u_r_rms = std::sqrt(ur2_acc / static_cast<Real>(w_safe));
    s.u_theta_rms = std::sqrt(uth2_acc / static_cast<Real>(w_safe));
    s.u_z_rms = std::sqrt(uz2_acc / static_cast<Real>(w_safe));
    return s;
}

inline FrenchMeltSpatialStats hostFrenchMeltSpatialStats(BaseParticles &particles, const FrenchCylindricalFrame &frame,
                                                         const std::string &temperature_field)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, temperature_field);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_field);
    return computeFrenchMeltSpatialStats(pos, vel, temperature, vol, particles.TotalRealParticles(), frame);
}

inline void writeFrenchSpatialStatsCsvHeader(const std::string &path)
{
    std::ofstream out(path);
    out << "time,T_mean,T_std,T_outer,T_center,T_outer_over_center,T_bottom,T_top,"
           "U_rms,U_r_rms,U_theta_rms,U_z_rms,U_max\n";
}

inline void appendFrenchSpatialStatsCsv(const std::string &path, Real time, const FrenchMeltSpatialStats &s)
{
    std::ofstream out(path, std::ios::app);
    out << time << "," << s.t_mean << "," << s.t_std << "," << s.t_outer << "," << s.t_center << ","
        << s.t_outer_over_center << "," << s.t_bottom << "," << s.t_top << "," << s.u_rms << "," << s.u_r_rms << ","
        << s.u_theta_rms << "," << s.u_z_rms << "," << s.u_max << "\n";
}

inline void writeAxisymmetricCsv(const std::string &path, const AxisymmetricRZGrid &grid, const char *value_name)
{
    std::ofstream out(path);
    out << "ir,iz,r_m,z_m,count,weight," << value_name << ",valid\n";
    const Real dr = grid.dr();
    const Real dz = grid.dz();
    for (size_t iz = 0; iz < grid.n_z; ++iz)
    {
        for (size_t ir = 0; ir < grid.n_r; ++ir)
        {
            const size_t id = grid.index(ir, iz);
            const Real r = (static_cast<Real>(ir) + Real(0.5)) * dr;
            const Real z = grid.frame.z_min + (static_cast<Real>(iz) + Real(0.5)) * dz;
            out << ir << "," << iz << "," << r << "," << z << "," << grid.count[id] << "," << grid.weight[id] << ","
                << grid.value[id] << "," << static_cast<int>(grid.valid[id]) << "\n";
        }
    }
}

inline void writeFluidTemperatureRzCsv(const std::string &path, BaseParticles &particles,
                                       const FrenchCylindricalFrame &frame, const std::string &temperature_field)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, temperature_field);
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *temperature = particles.getVariableDataByName<Real>(temperature_field);
    AxisymmetricRZGrid t_rz;
    t_rz.frame = frame;
    azimuthalDeposit(t_rz, pos, temperature, vol, particles.TotalRealParticles());
    writeAxisymmetricCsv(path, t_rz, "T_K");
}

inline void writeThermoEMCouplingCsvHeader(const std::string &path)
{
    std::ofstream out(path);
    out << "update,physical_time_s,em_control,coil_current,generator_power,target_glass_absorbed_power,"
           "reconstructed_glass_power,P_axisymmetric,P_euler_sample,P_sph,rel_axi,rel_euler,rel_sph,"
           "sigma_min,sigma_max,sigma_mean,phi_eq_res,remap_scale,conservative_remap\n";
}

inline void appendThermoEMCouplingCsv(const std::string &path, const ThermoEMCouplingState &state, Real physical_time,
                                      EMControlMode mode)
{
    std::ofstream out(path, std::ios::app);
    const auto &b = state.last_power;
    out << state.update_count << "," << physical_time << "," << emControlModeName(mode) << ","
        << state.coil_current_per_loop << "," << state.generator_power_w << ","
        << state.target_glass_absorbed_power_w << "," << state.reconstructed_glass_power_w << ","
        << b.p_axisymmetric_grid << "," << b.p_euler_sample << "," << b.p_sph_particles << ","
        << b.relErr(b.p_axisymmetric_grid, b.p_em_recon) << "," << b.relErr(b.p_euler_sample, b.p_em_recon) << ","
        << b.relErr(b.p_sph_particles, b.p_em_recon) << "," << state.last_sigma_min << "," << state.last_sigma_max
        << "," << state.last_sigma_mean << "," << state.last_phi_eq_res << "," << b.remap_scale << ","
        << (b.conservative_remap_applied ? 1 : 0) << "\n";
}

inline Real frenchEnergyResidual(Real dE_dt, Real joule_power, Real wall_loss, Real free_loss, Real bottom_loss)
{
    return dE_dt - (joule_power - wall_loss - free_loss - bottom_loss);
}

inline void writeFrenchEnergyBudgetCsvHeader(const std::string &path)
{
    std::ofstream out(path);
    out << "time,thermal_energy,dE_dt,joule_power,wall_heat_loss,free_surface_heat_loss,bottom_heat_loss,"
           "energy_residual,energy_residual_rel,T_min,T_mean,T_max,U_max,U_buoy,thermal_diffusion\n";
}

inline void appendFrenchEnergyBudgetCsv(const std::string &path, Real time, Real thermal_energy, Real dE_dt,
                                        Real joule_power, Real wall_loss, Real free_loss, Real bottom_loss,
                                        Real t_min = 0.0, Real t_mean = 0.0, Real t_max = 0.0, Real u_max = 0.0,
                                        Real u_buoy = 0.0, int thermal_diffusion = 0)
{
    const Real residual = frenchEnergyResidual(dE_dt, joule_power, wall_loss, free_loss, bottom_loss);
    const Real residual_rel = residual / (std::abs(joule_power) + TinyReal);
    std::ofstream out(path, std::ios::app);
    out << time << "," << thermal_energy << "," << dE_dt << "," << joule_power << "," << wall_loss << "," << free_loss
        << "," << bottom_loss << "," << residual << "," << residual_rel << "," << t_min << "," << t_mean << ","
        << t_max << "," << u_max << "," << u_buoy << "," << thermal_diffusion << "\n";
}

/**
 * Build σ_em for the next EM solve from SPH temperature (or a constant).
 * Returns false on invalid T/σ (caller must abort; do not continue with NaN).
 */
inline bool buildUnderRelaxedAxisymmetricSigma(const FrenchLiteratureParameters &lit, const FrenchCylindricalFrame &frame,
                                              const Vecd *pos, const Real *vol, const Real *temperature, size_t n,
                                              ThermoEMCouplingState &state, StdVec<Real> &sigma_em_out,
                                              std::string &error)
{
    StdVec<Real> sigma_raw(n, Real(0));
    if (lit.constant_sigma)
    {
        const Real s0 = lit.material.sigma.evaluate(lit.t_initial_k);
        for (size_t i = 0; i < n; ++i)
        {
            sigma_raw[i] = s0;
        }
    }
    else
    {
        if (!evaluateSigmaFromTemperature(lit.material.sigma, temperature, n, sigma_raw, error))
        {
            return false;
        }
    }

    state.sigma_rz.frame = frame;
    azimuthalDeposit(state.sigma_rz, pos, sigma_raw.data(), vol, n);
    Real smin = 0, smax = 0, smean = 0;
    sigmaFieldStats(sigma_raw, vol, n, smin, smax, smean);
    fillEmptyAxisymmetricBins(state.sigma_rz, smean > TinyReal ? smean : Real(1));

    StdVec<Real> sigma_bar;
    mapAxisymmetricToParticles(state.sigma_rz, pos, sigma_bar, n);
    if (state.sigma_em.size() != n)
    {
        state.sigma_em = sigma_bar;
    }
    sigma_em_out.resize(n);
    for (size_t i = 0; i < n; ++i)
    {
        sigma_em_out[i] = underRelaxScalar(state.sigma_em[i], sigma_bar[i], lit.sigma_under_relaxation);
        if (!std::isfinite(sigma_em_out[i]) || sigma_em_out[i] <= TinyReal)
        {
            error = "sigma under-relaxation produced invalid value";
            return false;
        }
    }
    state.sigma_em = sigma_em_out;
    sigmaFieldStats(sigma_em_out, vol, n, state.last_sigma_min, state.last_sigma_max, state.last_sigma_mean);
    return true;
}

/**
 * Fluid T → azimuthal σ_bar(r,z) → under-relaxed σ on the (possibly different) EM particle set.
 * Empty r–z bins are filled; invalid T/σ aborts via `error`.
 */
inline bool mapFluidTemperatureToUnderRelaxedEmSigma(const FrenchLiteratureParameters &lit,
                                                     const FrenchCylindricalFrame &frame, const Vecd *pos_fluid,
                                                     const Real *vol_fluid, const Real *t_fluid, size_t n_fluid,
                                                     const Vecd *pos_em, const Real *vol_em, size_t n_em,
                                                     ThermoEMCouplingState &state, StdVec<Real> &sigma_on_em,
                                                     std::string &error)
{
    StdVec<Real> sigma_raw;
    if (lit.constant_sigma)
    {
        const Real s0 = lit.material.sigma.evaluate(lit.t_initial_k);
        sigma_raw.assign(n_fluid, s0);
    }
    else if (!evaluateSigmaFromTemperature(lit.material.sigma, t_fluid, n_fluid, sigma_raw, error))
    {
        return false;
    }

    state.sigma_rz.frame = frame;
    azimuthalDeposit(state.sigma_rz, pos_fluid, sigma_raw.data(), vol_fluid, n_fluid);
    Real smin = 0, smax = 0, smean = 0;
    sigmaFieldStats(sigma_raw, vol_fluid, n_fluid, smin, smax, smean);
    fillEmptyAxisymmetricBins(state.sigma_rz, smean > TinyReal ? smean : Real(1));

    StdVec<Real> sigma_bar;
    mapAxisymmetricToParticles(state.sigma_rz, pos_em, sigma_bar, n_em);
    if (state.sigma_em.size() != n_em)
    {
        state.sigma_em = sigma_bar;
    }
    sigma_on_em.resize(n_em);
    for (size_t i = 0; i < n_em; ++i)
    {
        sigma_on_em[i] = underRelaxScalar(state.sigma_em[i], sigma_bar[i], lit.sigma_under_relaxation);
        if (!std::isfinite(sigma_on_em[i]) || sigma_on_em[i] <= TinyReal)
        {
            error = "sigma under-relaxation produced invalid value";
            return false;
        }
    }
    state.sigma_em = sigma_on_em;
    sigmaFieldStats(sigma_on_em, vol_em, n_em, state.last_sigma_min, state.last_sigma_max, state.last_sigma_mean);
    return true;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_THERMO_EM_COUPLING_H
