#ifndef ELECTROMAGNETIC_OPHELIE_THERMO_EM_COUPLING_H
#define ELECTROMAGNETIC_OPHELIE_THERMO_EM_COUPLING_H

/**
 * Periodic thermo–EM coupling helpers (Jacoutot CES 2008 Fig. 4 dataflow).
 *
 * Production Q path (Batch D / P1-02 — do not overstate as paper-style axi→3D):
 *
 *   3D T → sigma(T) → under-relaxed averaged sigma(r,z) → 3D EM solve
 *   → original 3D Q deposited on Euler grid → 3D fluid heating
 *
 * The r–z azimuthal average of Q is recorded for diagnostics / power books only.
 * It is not currently rotated back into 3D for heating. An explicit axi-Q mode
 * would be a separate opt-in.
 *
 * Q remapping never silently rescales. Conservative remap is opt-in and logged.
 * A single under-relaxed sigma update per em_update_interval introduces lag;
 * nondivergence alone does not prove coupling convergence.
 */

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_literature_parameters.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
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

inline Real hostParticlePower(const Real *q, const Real *vol, size_t n)
{
    Real p = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        p += q[i] * vol[i];
    }
    return p;
}

struct FrenchPowerMapBook
{
    /** Schema: 2 = explicit P_em_particles / P_grid_at_em / P_heat_fluid (P_sph kept as alias). */
    int schema_version = 2;
    Real p_em_recon = 0.0;
    Real p_axisymmetric_grid = 0.0;
    /** Grid sampled at EM particle positions with EM volumes (doc 23: P_grid_at_em). */
    Real p_euler_sample = 0.0;
    /** Integral of EM-body Q·V — NOT the moving fluid heating power. Alias of p_em_particles. */
    Real p_sph_particles = 0.0;
    /** Same as p_sph_particles; unambiguous name for EM-body Q integral. */
    Real p_em_particles = 0.0;
    /** Fluid-body JouleHeat·V after grid sampling (doc 23: P_heat_fluid). */
    Real p_heat_fluid = 0.0;
    /** Alias used in energy books: last sampled fluid heating, not EM recon. */
    Real p_heat_fluid_sampled = 0.0;
    Real remap_scale = 1.0;
    Real p_euler_before_remap = 0.0;
    Real p_euler_after_remap = 0.0;
    bool conservative_remap_applied = false;
    bool fluid_power_recorded = false;

    Real relErr(Real a, Real b) const { return std::abs(a - b) / (std::abs(b) + TinyReal); }
};

inline void logFrenchPowerMapBook(const FrenchPowerMapBook &book)
{
    std::cout << "[ophelie][q-map] schema=" << book.schema_version << " P_em_recon=" << book.p_em_recon
              << " P_axisymmetric=" << book.p_axisymmetric_grid
              << " P_grid_at_em=" << book.p_euler_sample << " P_em_particles=" << book.p_em_particles
              << " P_sph=" << book.p_sph_particles << " (alias_em_particles)"
              << " P_heat_fluid=" << book.p_heat_fluid
              << " fluid_recorded=" << (book.fluid_power_recorded ? 1 : 0)
              << " rel_axi=" << book.relErr(book.p_axisymmetric_grid, book.p_em_recon)
              << " rel_grid_em=" << book.relErr(book.p_euler_sample, book.p_em_recon)
              << " rel_em_particles=" << book.relErr(book.p_em_particles, book.p_em_recon)
              << " rel_fluid=" << book.relErr(book.p_heat_fluid, book.p_em_recon)
              << " conservative_remap=" << (book.conservative_remap_applied ? 1 : 0);
    if (book.conservative_remap_applied)
    {
        std::cout << " scale=" << book.remap_scale << " P_before=" << book.p_euler_before_remap
                  << " P_after=" << book.p_euler_after_remap;
    }
    std::cout << std::endl;
}

/**
 * Integrated volumetric heating on a body: sum_i Q_i V_i [W].
 * For fluid after SampleJouleHeatFromGridCK, this is the heat source actually seen by T += Q dt/(rho cp).
 */
inline Real hostBodyJouleHeatPower(BaseParticles &particles, const std::string &joule_field = "JouleHeat")
{
    syncVariableToHost<Real>(particles, joule_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *q = particles.getVariableDataByName<Real>(joule_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    return hostParticlePower(q, vol, particles.TotalRealParticles());
}

/**
 * Cumulative applied-energy book (doc 40 R01 / doc 41).
 *   R_E = H(t) - H(t0) - E_Q_applied + E_boundary_out
 * E_Q_applied integrates the fluid sampled source used by T += Q dt/(ρ cp).
 * E_boundary_out integrates the same discrete face flux the Robin kernel applies.
 * Clamp/imposed-T energy is not included until those operators report it.
 */
struct FrenchAppliedEnergyLedger
{
    double h0 = 0.0;
    double e_q_applied = 0.0;
    double e_boundary_out = 0.0;
    bool h0_set = false;
    /** Sticky: any non-finite accumulate input sets this; gate must fail (doc 45 F03). */
    bool invalid_accumulate = false;
    Real p_heat_fluid_live = 0.0;
    Real p_loss_live = 0.0;
    Real p_side_live = 0.0;
    Real p_bottom_live = 0.0;
    Real p_free_live = 0.0;
    size_t fluid_sample_count = 0;
    Real last_out_of_grid_sampling = 0.0;
    Real max_out_of_grid_sampling = 0.0;

    void setH0(Real h)
    {
        h0 = static_cast<double>(h);
        h0_set = true;
    }

    void accumulateQ(Real p_heat_fluid, Real dt)
    {
        if (!(std::isfinite(p_heat_fluid) && std::isfinite(dt)) || dt <= Real(0))
        {
            invalid_accumulate = true;
            return;
        }
        e_q_applied += static_cast<double>(p_heat_fluid) * static_cast<double>(dt);
    }

    void accumulateBoundary(Real p_loss, Real dt)
    {
        if (!(std::isfinite(p_loss) && std::isfinite(dt)) || dt <= Real(0))
        {
            invalid_accumulate = true;
            return;
        }
        e_boundary_out += static_cast<double>(p_loss) * static_cast<double>(dt);
    }

    void noteLiveLoss(Real side_w, Real bottom_w, Real free_w)
    {
        p_side_live = side_w;
        p_bottom_live = bottom_w;
        p_free_live = free_w;
        p_loss_live = side_w + bottom_w + free_w;
    }

    Real residualH(Real h_now) const
    {
        return static_cast<Real>((static_cast<double>(h_now) - h0) - e_q_applied + e_boundary_out);
    }

    Real residualHRel(Real h_now) const
    {
        const double dh = static_cast<double>(h_now) - h0;
        const double residual = dh - e_q_applied + e_boundary_out;
        const double scale = std::max({std::abs(e_q_applied), std::abs(e_boundary_out), std::abs(dh), 1.0});
        return static_cast<Real>(std::abs(residual) / scale);
    }

    /** Residual relative to |predicted net storage| = |E_Q - E_boundary|; floor when near zero (doc 45 F06). */
    Real residualHRelVsNet(Real h_now) const
    {
        const double dh = static_cast<double>(h_now) - h0;
        const double residual = dh - e_q_applied + e_boundary_out;
        const double net_pred = e_q_applied - e_boundary_out;
        const double scale = std::max(std::abs(net_pred), 1.0);
        return static_cast<Real>(std::abs(residual) / scale);
    }
};

struct ThermoEMCouplingState
{
    Real last_update_time = -1.0e300;
    /** Next nominal refresh time (Doc 47 §5); avoids drift from last+interval. */
    Real next_due_time = -1.0e300;
    /** Max (t_refresh − last_update) seen at periodic refreshes; end-of-run age folded in at gate. */
    Real max_field_age_s = 0.0;
    /** Max (t_refresh − next_due) when a refresh fires late. */
    Real max_nominal_lateness_s = 0.0;
    bool freshness_measured = false;
    size_t refresh_solver_fail_count = 0;
    /** Total coupling CSV rows (includes the t=0 initialization book). */
    size_t update_count = 0;
    /** Periodic refreshes after initialization only (doc 23 P0-02). */
    size_t periodic_refresh_count = 0;
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
    Real max_out_of_grid = 0.0;
    Real last_out_of_grid = 0.0;
    FrenchAppliedEnergyLedger ledger;
};

/** Record fluid heating after grid→particle sampling; updates last_power and logs. */
inline void recordFrenchFluidHeatPower(ThermoEMCouplingState &state, BaseParticles &fluid_particles,
                                       const std::string &joule_field = "JouleHeat")
{
    state.last_power.p_heat_fluid = hostBodyJouleHeatPower(fluid_particles, joule_field);
    state.last_power.p_heat_fluid_sampled = state.last_power.p_heat_fluid;
    state.last_power.fluid_power_recorded = true;
    state.last_power.schema_version = 2;
    state.ledger.p_heat_fluid_live = state.last_power.p_heat_fluid;
    ++state.ledger.fluid_sample_count;
    std::cout << "[ophelie][q-map] P_heat_fluid=" << state.last_power.p_heat_fluid
              << " rel_fluid_vs_em=" << state.last_power.relErr(state.last_power.p_heat_fluid, state.last_power.p_em_recon)
              << " rel_fluid_vs_grid_em="
              << state.last_power.relErr(state.last_power.p_heat_fluid, state.last_power.p_euler_sample)
              << " fluid_samples=" << state.ledger.fluid_sample_count << std::endl;
}

/** Same as recordFrenchFluidHeatPower, plus out-of-grid history during fluid sampling (doc 40 R01). */
inline void recordFrenchFluidHeatSample(ThermoEMCouplingState &state, BaseParticles &fluid_particles, Real out_of_grid,
                                        const std::string &joule_field = "JouleHeat")
{
    recordFrenchFluidHeatPower(state, fluid_particles, joule_field);
    if (out_of_grid >= Real(0) && std::isfinite(out_of_grid))
    {
        state.ledger.last_out_of_grid_sampling = out_of_grid;
        state.ledger.max_out_of_grid_sampling = std::max(state.ledger.max_out_of_grid_sampling, out_of_grid);
        state.last_out_of_grid = out_of_grid;
        state.max_out_of_grid = std::max(state.max_out_of_grid, out_of_grid);
    }
}

/**
 * Doc 47 §5: fire when physical_time reaches the nominal next-due tick (not last+interval).
 * `next_due_time` < -1e100 means "not armed" → allow (caller should arm after t=0 init).
 */
inline bool shouldUpdateEM(Real physical_time, Real interval_s, Real next_due_time)
{
    if (interval_s <= TinyReal)
    {
        return true;
    }
    if (next_due_time < Real(-1.0e100))
    {
        return true;
    }
    return physical_time + TinyReal >= next_due_time;
}

/** After t=0 initialization book: last_update=0, first periodic due at interval. */
inline void armFrenchPeriodicRefreshSchedule(ThermoEMCouplingState &state, Real interval_s)
{
    state.last_update_time = Real(0);
    state.next_due_time = (interval_s > TinyReal) ? interval_s : Real(0);
    state.max_field_age_s = Real(0);
    state.max_nominal_lateness_s = Real(0);
    state.freshness_measured = false;
    state.refresh_solver_fail_count = 0;
}

/**
 * Record a completed periodic refresh: field age, nominal lateness, advance next-due.
 * If late by more than one interval, advance next_due past t_now (do not invent extra refreshes).
 */
inline void noteFrenchPeriodicRefresh(ThermoEMCouplingState &state, Real t_now, Real interval_s, bool solver_ok)
{
    if (state.last_update_time > Real(-1.0e100))
    {
        state.max_field_age_s = std::max(state.max_field_age_s, t_now - state.last_update_time);
    }
    if (state.next_due_time > Real(-1.0e100))
    {
        state.max_nominal_lateness_s =
            std::max(state.max_nominal_lateness_s, t_now - state.next_due_time);
    }
    state.last_update_time = t_now;
    state.freshness_measured = true;
    if (!solver_ok)
    {
        ++state.refresh_solver_fail_count;
    }
    if (interval_s > TinyReal)
    {
        if (state.next_due_time < Real(-1.0e100))
        {
            state.next_due_time = t_now + interval_s;
        }
        else
        {
            state.next_due_time += interval_s;
            while (state.next_due_time <= t_now + TinyReal)
            {
                state.next_due_time += interval_s;
            }
        }
    }
}

/** End-of-run field age since last successful refresh (for gate). */
inline Real frenchCurrentFieldAge(const ThermoEMCouplingState &state, Real physical_time)
{
    if (state.last_update_time < Real(-1.0e100))
    {
        return Real(-1);
    }
    return physical_time - state.last_update_time;
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

/**
 * Time-averaged axisymmetric meridional (U_r, U_z) diagnostic for Round 3C.
 * Accumulates volume-weighted signed components on an (r,z) grid so a slow roll
 * can be separated from instantaneous WCSPH noise / U_max spikes.
 */
struct FrenchMeridionalFlowAccumulator
{
    FrenchCylindricalFrame frame;
    size_t n_r = 16;
    size_t n_z = 12;
    StdVec<double> ur_acc;
    StdVec<double> uz_acc;
    StdVec<double> w_acc;
    size_t n_samples = 0;
    Real t_first = -1.0;
    Real t_last = 0.0;
    double u_rms_sum = 0.0;
    double ur_rms_sum = 0.0;
    double uz_rms_sum = 0.0;
    double circ_proxy_sum = 0.0; // volume mean of U_z * (2 r/R - 1)

    void reset(const FrenchCylindricalFrame &f, size_t nr = 16, size_t nz = 12)
    {
        frame = f;
        n_r = std::max<size_t>(nr, 1);
        n_z = std::max<size_t>(nz, 1);
        const size_t n = n_r * n_z;
        ur_acc.assign(n, 0.0);
        uz_acc.assign(n, 0.0);
        w_acc.assign(n, 0.0);
        n_samples = 0;
        t_first = -1.0;
        t_last = 0.0;
        u_rms_sum = ur_rms_sum = uz_rms_sum = circ_proxy_sum = 0.0;
    }

    size_t index(size_t ir, size_t iz) const { return ir + n_r * iz; }

    bool binIndex(Real r, Real z, size_t &ir, size_t &iz) const
    {
        if (r < 0.0 || z < frame.z_min || z > frame.z_max)
        {
            return false;
        }
        const Real R = std::max(frame.radius, TinyReal);
        const Real H = std::max(frame.z_max - frame.z_min, TinyReal);
        ir = static_cast<size_t>(std::min(std::floor(r / R * static_cast<Real>(n_r)), static_cast<Real>(n_r - 1)));
        iz = static_cast<size_t>(
            std::min(std::floor((z - frame.z_min) / H * static_cast<Real>(n_z)), static_cast<Real>(n_z - 1)));
        return true;
    }

    void accumulate(const Vecd *pos, const Vecd *vel, const Real *vol, size_t n, Real time)
    {
        if (pos == nullptr || vel == nullptr || vol == nullptr || n == 0)
        {
            return;
        }
        const int ax = frame.axis;
        Vecd e_z = Vecd::Zero();
        e_z[ax] = Real(1);
        const Real R = std::max(frame.radius, TinyReal);

        double u2 = 0.0;
        double ur2 = 0.0;
        double uz2 = 0.0;
        double circ = 0.0;
        double wtot = 0.0;

        for (size_t i = 0; i < n; ++i)
        {
            if (vol[i] <= TinyReal)
            {
                continue;
            }
            Real r = 0.0;
            Real z = 0.0;
            frenchCylindricalRZ(pos[i], frame, r, z);
            size_t ir = 0;
            size_t iz = 0;
            if (!binIndex(r, z, ir, iz))
            {
                continue;
            }
            const Real uz = vel[i].dot(e_z);
            Vecd e_r = pos[i] - frame.origin;
            e_r[ax] = Real(0);
            const Real rxy = e_r.norm();
            Real ur = Real(0);
            if (rxy > TinyReal)
            {
                e_r /= rxy;
                ur = vel[i].dot(e_r);
            }
            const double vi = static_cast<double>(vol[i]);
            const size_t id = index(ir, iz);
            ur_acc[id] += static_cast<double>(ur) * vi;
            uz_acc[id] += static_cast<double>(uz) * vi;
            w_acc[id] += vi;

            const Real un = vel[i].norm();
            if (std::isfinite(un))
            {
                u2 += static_cast<double>(un * un) * vi;
                ur2 += static_cast<double>(ur * ur) * vi;
                uz2 += static_cast<double>(uz * uz) * vi;
                circ += static_cast<double>(uz) * (Real(2) * r / R - Real(1)) * vi;
                wtot += vi;
            }
        }

        if (wtot > static_cast<double>(TinyReal))
        {
            u_rms_sum += std::sqrt(u2 / wtot);
            ur_rms_sum += std::sqrt(ur2 / wtot);
            uz_rms_sum += std::sqrt(uz2 / wtot);
            circ_proxy_sum += circ / wtot;
            ++n_samples;
            if (t_first < Real(0))
            {
                t_first = time;
            }
            t_last = time;
        }
    }

    Real meanMeridionalSpeed() const
    {
        double num = 0.0;
        double den = 0.0;
        for (size_t id = 0; id < w_acc.size(); ++id)
        {
            if (w_acc[id] <= TinyReal)
            {
                continue;
            }
            const double ur = ur_acc[id] / w_acc[id];
            const double uz = uz_acc[id] / w_acc[id];
            num += std::sqrt(ur * ur + uz * uz) * w_acc[id];
            den += w_acc[id];
        }
        return static_cast<Real>(num / (den + TinyReal));
    }

    Real meanCircProxy() const
    {
        return n_samples > 0 ? static_cast<Real>(circ_proxy_sum / static_cast<double>(n_samples)) : Real(0);
    }

    Real meanInstUrms() const
    {
        return n_samples > 0 ? static_cast<Real>(u_rms_sum / static_cast<double>(n_samples)) : Real(0);
    }

    Real meanInstUrRms() const
    {
        return n_samples > 0 ? static_cast<Real>(ur_rms_sum / static_cast<double>(n_samples)) : Real(0);
    }

    Real meanInstUzRms() const
    {
        return n_samples > 0 ? static_cast<Real>(uz_rms_sum / static_cast<double>(n_samples)) : Real(0);
    }

    /** Volume-weighted mean U_z on time-averaged bins with r ∈ [r_lo, r_hi). */
    Real meanUzRadiusBand(Real r_lo, Real r_hi) const
    {
        double num = 0.0;
        double den = 0.0;
        const Real dr = frame.radius / static_cast<Real>(n_r);
        for (size_t iz = 0; iz < n_z; ++iz)
        {
            for (size_t ir = 0; ir < n_r; ++ir)
            {
                const Real r = (static_cast<Real>(ir) + Real(0.5)) * dr;
                if (r < r_lo || r >= r_hi)
                    continue;
                const size_t id = index(ir, iz);
                if (w_acc[id] <= TinyReal)
                    continue;
                num += (uz_acc[id] / w_acc[id]) * w_acc[id];
                den += w_acc[id];
            }
        }
        return static_cast<Real>(num / (den + TinyReal));
    }

    void writeCsv(const std::string &path) const
    {
        std::ofstream out(path);
        out << "ir,iz,r_m,z_m,weight,U_r_mean,U_z_mean,U_merid,valid\n";
        const Real dr = frame.radius / static_cast<Real>(n_r);
        const Real dz = (frame.z_max - frame.z_min) / static_cast<Real>(n_z);
        for (size_t iz = 0; iz < n_z; ++iz)
        {
            for (size_t ir = 0; ir < n_r; ++ir)
            {
                const size_t id = index(ir, iz);
                const Real r = (static_cast<Real>(ir) + Real(0.5)) * dr;
                const Real z = frame.z_min + (static_cast<Real>(iz) + Real(0.5)) * dz;
                const int valid = w_acc[id] > TinyReal ? 1 : 0;
                const double ur = valid ? ur_acc[id] / w_acc[id] : 0.0;
                const double uz = valid ? uz_acc[id] / w_acc[id] : 0.0;
                out << ir << "," << iz << "," << r << "," << z << "," << w_acc[id] << "," << ur << "," << uz << ","
                    << std::sqrt(ur * ur + uz * uz) << "," << valid << "\n";
            }
        }
    }

    void logSummary(const char *tag) const
    {
        std::cout << "[ophelie][meridional] tag=" << tag << " samples=" << n_samples << " window=[" << t_first << ","
                  << t_last << "] U_merid_mean=" << meanMeridionalSpeed() << " circ_proxy=" << meanCircProxy()
                  << " U_rms_inst_avg=" << meanInstUrms() << " U_r_rms_inst_avg=" << meanInstUrRms()
                  << " U_z_rms_inst_avg=" << meanInstUzRms() << "\n";
    }
};

inline void hostAccumulateFrenchMeridionalFlow(FrenchMeridionalFlowAccumulator &acc, BaseParticles &particles,
                                               Real time)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    acc.accumulate(pos, vel, vol, particles.TotalRealParticles(), time);
}

/** Doc97 P3: one-shot volume-weighted band / circ scalars (for non-overlapping time windows). */
struct FrenchMeridionalInstantBands
{
    Real uz_axis = 0;
    Real uz_rim = 0;
    Real circ_proxy = 0;
    Real u_rms = 0;
    Real u_merid = 0;
    double w_axis = 0;
    double w_rim = 0;
    double w_tot = 0;
};

inline FrenchMeridionalInstantBands hostFrenchMeridionalInstantBands(BaseParticles &particles,
                                                                    const FrenchCylindricalFrame &frame,
                                                                    Real r_axis_hi, Real r_rim_lo)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    const int ax = frame.axis;
    Vecd e_z = Vecd::Zero();
    e_z[ax] = Real(1);
    const Real R = std::max(frame.radius, TinyReal);
    FrenchMeridionalInstantBands out;
    double num_ax = 0, num_rim = 0, circ = 0, u2 = 0, umer = 0;
    for (size_t i = 0; i < n; ++i)
    {
        if (vol[i] <= TinyReal)
            continue;
        Real r = 0, z = 0;
        frenchCylindricalRZ(pos[i], frame, r, z);
        if (r < 0.0 || z < frame.z_min || z > frame.z_max)
            continue;
        const Real uz = vel[i].dot(e_z);
        Vecd e_r = pos[i] - frame.origin;
        e_r[ax] = Real(0);
        const Real rxy = e_r.norm();
        Real ur = Real(0);
        if (rxy > TinyReal)
        {
            e_r /= rxy;
            ur = vel[i].dot(e_r);
        }
        const double vi = static_cast<double>(vol[i]);
        out.w_tot += vi;
        u2 += static_cast<double>(vel[i].squaredNorm()) * vi;
        umer += std::sqrt(static_cast<double>(ur * ur + uz * uz)) * vi;
        circ += static_cast<double>(uz) * (Real(2) * r / R - Real(1)) * vi;
        if (r < r_axis_hi)
        {
            num_ax += static_cast<double>(uz) * vi;
            out.w_axis += vi;
        }
        if (r >= r_rim_lo)
        {
            num_rim += static_cast<double>(uz) * vi;
            out.w_rim += vi;
        }
    }
    if (out.w_axis > static_cast<double>(TinyReal))
        out.uz_axis = static_cast<Real>(num_ax / out.w_axis);
    if (out.w_rim > static_cast<double>(TinyReal))
        out.uz_rim = static_cast<Real>(num_rim / out.w_rim);
    if (out.w_tot > static_cast<double>(TinyReal))
    {
        out.circ_proxy = static_cast<Real>(circ / out.w_tot);
        out.u_rms = static_cast<Real>(std::sqrt(u2 / out.w_tot));
        out.u_merid = static_cast<Real>(umer / out.w_tot);
    }
    return out;
}

inline void appendFrenchMeridionalTimeseriesCsv(const std::string &path, Real time,
                                                const FrenchMeridionalInstantBands &b, bool write_header)
{
    std::ofstream out(path, write_header ? std::ios::out : std::ios::app);
    if (write_header)
        out << "time_s,Uz_axis,Uz_rim,circ_proxy,U_rms,U_merid,w_axis,w_rim,w_tot\n";
    out << time << "," << b.uz_axis << "," << b.uz_rim << "," << b.circ_proxy << "," << b.u_rms << "," << b.u_merid
        << "," << b.w_axis << "," << b.w_rim << "," << b.w_tot << "\n";
}

/** Doc97 P3: time-weighted mean over a half-open physical window [t0,t1). */
struct FrenchMeridionalWindowMean
{
    Real t0 = 0;
    Real t1 = 0;
    Real uz_axis = 0;
    Real uz_rim = 0;
    Real circ_proxy = 0;
    Real u_rms = 0;
    Real u_merid = 0;
    Real weight_s = 0;
    int n_samples = 0;
    bool ok = false;
};

inline std::vector<FrenchMeridionalWindowMean> frenchMeridionalWindowsFromTimeseriesCsv(
    const std::string &path, const std::vector<std::pair<Real, Real>> &windows)
{
    std::vector<FrenchMeridionalWindowMean> out;
    out.reserve(windows.size());
    for (const auto &w : windows)
    {
        FrenchMeridionalWindowMean m;
        m.t0 = w.first;
        m.t1 = w.second;
        out.push_back(m);
    }
    std::ifstream in(path);
    if (!in)
        return out;
    std::string line;
    if (!std::getline(in, line))
        return out; // header
    struct Sample
    {
        Real t = 0;
        Real uz_axis = 0;
        Real uz_rim = 0;
        Real circ = 0;
        Real u_rms = 0;
        Real u_merid = 0;
    };
    std::vector<Sample> samples;
    while (std::getline(in, line))
    {
        if (line.empty())
            continue;
        Sample s;
        char c = 0;
        std::stringstream ss(line);
        if (!(ss >> s.t >> c >> s.uz_axis >> c >> s.uz_rim >> c >> s.circ >> c >> s.u_rms >> c >> s.u_merid))
            continue;
        samples.push_back(s);
    }
    if (samples.size() < 2)
        return out;
    for (size_t i = 0; i + 1 < samples.size(); ++i)
    {
        const Sample &a = samples[i];
        const Real t_lo = a.t;
        const Real t_hi = samples[i + 1].t;
        if (!(t_hi > t_lo))
            continue;
        for (auto &m : out)
        {
            const Real lo = std::max(t_lo, m.t0);
            const Real hi = std::min(t_hi, m.t1);
            if (!(hi > lo))
                continue;
            const double dt = static_cast<double>(hi - lo);
            m.weight_s = static_cast<Real>(static_cast<double>(m.weight_s) + dt);
            m.uz_axis = static_cast<Real>(static_cast<double>(m.uz_axis) + dt * static_cast<double>(a.uz_axis));
            m.uz_rim = static_cast<Real>(static_cast<double>(m.uz_rim) + dt * static_cast<double>(a.uz_rim));
            m.circ_proxy = static_cast<Real>(static_cast<double>(m.circ_proxy) + dt * static_cast<double>(a.circ));
            m.u_rms = static_cast<Real>(static_cast<double>(m.u_rms) + dt * static_cast<double>(a.u_rms));
            m.u_merid = static_cast<Real>(static_cast<double>(m.u_merid) + dt * static_cast<double>(a.u_merid));
            ++m.n_samples;
        }
    }
    for (auto &m : out)
    {
        if (m.weight_s > TinyReal && m.n_samples > 0)
        {
            m.uz_axis /= m.weight_s;
            m.uz_rim /= m.weight_s;
            m.circ_proxy /= m.weight_s;
            m.u_rms /= m.weight_s;
            m.u_merid /= m.weight_s;
            m.ok = true;
        }
    }
    return out;
}

inline void writeFrenchMeridionalWindowsJsonArray(std::ostream &js,
                                                  const std::vector<FrenchMeridionalWindowMean> &wins)
{
    js << "  \"meridional_windows\": [\n";
    for (size_t i = 0; i < wins.size(); ++i)
    {
        const auto &w = wins[i];
        js << "    {\"t0\": " << w.t0 << ", \"t1\": " << w.t1 << ", \"ok\": " << (w.ok ? 1 : 0)
           << ", \"n_samples\": " << w.n_samples << ", \"weight_s\": " << w.weight_s
           << ", \"Uz_axis\": " << w.uz_axis << ", \"Uz_rim\": " << w.uz_rim
           << ", \"circ_proxy\": " << w.circ_proxy << ", \"U_rms\": " << w.u_rms
           << ", \"U_merid\": " << w.u_merid << "}" << (i + 1 < wins.size() ? "," : "") << "\n";
    }
    js << "  ]";
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
           "reconstructed_glass_power,P_axisymmetric,P_grid_at_em,P_em_particles,P_sph,P_heat_fluid,"
           "rel_axi,rel_grid_em,rel_em_particles,rel_fluid,rel_sph,"
           "sigma_min,sigma_max,sigma_mean,phi_eq_res,remap_scale,conservative_remap,"
           "fluid_recorded,schema_version,periodic_refresh_count\n";
}

inline void appendThermoEMCouplingCsv(const std::string &path, const ThermoEMCouplingState &state, Real physical_time,
                                      EMControlMode mode)
{
    std::ofstream out(path, std::ios::app);
    const auto &b = state.last_power;
    out << state.update_count << "," << physical_time << "," << emControlModeName(mode) << ","
        << state.coil_current_per_loop << "," << state.generator_power_w << ","
        << state.target_glass_absorbed_power_w << "," << state.reconstructed_glass_power_w << ","
        << b.p_axisymmetric_grid << "," << b.p_euler_sample << "," << b.p_em_particles << "," << b.p_sph_particles
        << "," << b.p_heat_fluid << "," << b.relErr(b.p_axisymmetric_grid, b.p_em_recon) << ","
        << b.relErr(b.p_euler_sample, b.p_em_recon) << "," << b.relErr(b.p_em_particles, b.p_em_recon) << ","
        << b.relErr(b.p_heat_fluid, b.p_em_recon) << "," << b.relErr(b.p_sph_particles, b.p_em_recon) << ","
        << state.last_sigma_min << "," << state.last_sigma_max << "," << state.last_sigma_mean << ","
        << state.last_phi_eq_res << "," << b.remap_scale << "," << (b.conservative_remap_applied ? 1 : 0) << ","
        << (b.fluid_power_recorded ? 1 : 0) << "," << b.schema_version << "," << state.periodic_refresh_count
        << "\n";
}

/** Honest smoke vs validation gates (doc 23 P0-02, doc 40 R04). */
struct FrenchCouplingRunGates
{
    bool startup_ok = false;
    bool completed_requested_time = false;
    bool periodic_refresh_ok = false;
    bool em_solver_ok = false;
    bool q_transfer_ok = false;
    bool flow_sanity_ok = false;
    bool thermal_energy_ok = false;
    bool thermal_finite_ok = false;
    bool freshness_ok = false;
    bool smoke_passed = false;
    bool validation_exercised = false;
    bool validation_passed = false;
    size_t expected_periodic_refreshes = 0;
    size_t periodic_refresh_count = 0;
    std::string fail_reasons;
};

inline constexpr Real kFrenchValidationEnergyRelTol = Real(0.02);
inline constexpr size_t kFrenchValidationMinPeriodicRefreshes = 3;
/**
 * Absolute end-time tolerance [s] for "completed requested time" (Doc 47 §5).
 * Not a fraction of the run duration.
 */
inline constexpr Real kFrenchValidationEndTimeTol = Real(1.0e-3);

struct FrenchCouplingGateInputs
{
    bool periodic_coupling = false;
    bool zero_q = false;
    bool constant_thermal_rheology = true;
    bool require_u_max_positive = true;
    bool energy_residual_available = false;
    bool budget_exhausted = false;
    bool diverged = false;
    bool fluid_power_recorded = false;
    bool ledger_invalid_accumulate = false;
    bool out_of_grid_measured = true;
    Real end_time = 0.0;
    Real em_interval = 0.0;
    Real physical_time = 0.0;
    Real p_em = 0.0;
    Real p_heat_fluid = 0.0;
    Real phi_eq = 0.0;
    Real power_rel_err = 0.0;
    /** <0 means not measured (natural host path). */
    Real q_sample_rel_l2 = Real(-1);
    Real n_out_of_grid_init = 0.0;
    Real n_out_of_grid_now = 0.0;
    Real max_out_of_grid = 0.0;
    Real u_max = 0.0;
    Real c0 = 0.0;
    Real t_mean = 0.0;
    Real t_max = 0.0;
    Real escaped_frac = 0.0;
    Real escaped_frac_max = 1.0;
    Real cumulative_energy_residual_rel = 0.0;
    size_t update_count = 0;
    size_t periodic_refresh_count = 0;
    /** Doc 47 §5 freshness; <0 means not measured (do not certify validation). */
    Real max_field_age_s = Real(-1);
    Real max_nominal_lateness_s = Real(-1);
    size_t refresh_solver_fail_count = 0;
    bool freshness_measured = false;
};

/**
 * Doc 47 §5 exclusive-endpoint schedule: count k with k*interval < end_time
 * (refresh exactly at end_time is optional). Init at t=0 is separate and not counted.
 * No blanket +1 slack — that previously accepted a missing *interior* tick.
 *
 * Use kFrenchValidationEndTimeTol (not TinyReal): with float Real, TinyReal≪ulp(end)
 * so `end - TinyReal == end` and floor(end/interval) wrongly includes the endpoint.
 */
inline size_t frenchExpectedPeriodicRefreshes(Real end_time, Real interval_s)
{
    if (interval_s <= TinyReal || end_time <= TinyReal)
    {
        return 0;
    }
    const Real exclusive_end = end_time - kFrenchValidationEndTimeTol;
    if (exclusive_end <= TinyReal)
    {
        return 0;
    }
    const Real n = std::floor(exclusive_end / interval_s);
    return n > Real(0) ? static_cast<size_t>(n) : 0;
}

inline void appendGateFail(std::string &reasons, const char *tag)
{
    if (!reasons.empty())
    {
        reasons += ",";
    }
    reasons += tag;
}

inline void evaluateFrenchCouplingRunGates(const FrenchCouplingGateInputs &in, FrenchCouplingRunGates &g)
{
    g = FrenchCouplingRunGates{};
    g.periodic_refresh_count = in.periodic_refresh_count;
    g.expected_periodic_refreshes =
        in.periodic_coupling ? frenchExpectedPeriodicRefreshes(in.end_time, in.em_interval) : 0;

    const bool power_present = in.zero_q || (std::isfinite(in.p_em) && in.p_em > TinyReal);
    g.startup_ok = power_present && (!in.periodic_coupling || in.update_count >= 1);
    g.completed_requested_time =
        !in.budget_exhausted && in.physical_time + kFrenchValidationEndTimeTol >= in.end_time;

    // Doc 47 F03: minimum exercise count AND exclusive-endpoint schedule (no slack).
    const size_t min_refresh =
        std::min(kFrenchValidationMinPeriodicRefreshes, g.expected_periodic_refreshes);
    const bool enough_refreshes =
        !in.periodic_coupling ||
        (g.expected_periodic_refreshes == 0
             ? (in.periodic_refresh_count >= kFrenchValidationMinPeriodicRefreshes)
             : (in.periodic_refresh_count >= g.expected_periodic_refreshes &&
                in.periodic_refresh_count >= min_refresh));
    g.periodic_refresh_ok = enough_refreshes;

    g.em_solver_ok = in.zero_q || (std::isfinite(in.p_em) &&
                                   (in.periodic_coupling ? (std::isfinite(in.phi_eq) && in.phi_eq < Real(1.0e-2))
                                                         : (in.power_rel_err < Real(1.0e-2))));

    const bool q_l2_ok = in.q_sample_rel_l2 < Real(0) ||
                         (std::isfinite(in.q_sample_rel_l2) && in.q_sample_rel_l2 < Real(0.2));
    // Doc 45 F03: nonzero EM with zero/unrecorded fluid heat is a transfer failure.
    const bool transfer_power_ok =
        in.zero_q ||
        (in.fluid_power_recorded && std::isfinite(in.p_heat_fluid) && in.p_heat_fluid > TinyReal);
    const bool out_of_grid_ok =
        !in.out_of_grid_measured ||
        (in.n_out_of_grid_now < Real(1) && in.max_out_of_grid < Real(1));
    g.q_transfer_ok = q_l2_ok && out_of_grid_ok &&
                      (!in.periodic_coupling || transfer_power_ok);

    // Doc 47 §5: unmeasured transfer / out-of-grid must not certify validation_passed.
    const bool transfer_criteria_measured =
        !in.periodic_coupling ||
        (in.out_of_grid_measured && in.q_sample_rel_l2 >= Real(0) && in.fluid_power_recorded);

    g.flow_sanity_ok = !in.diverged && std::isfinite(in.u_max) && in.u_max <= in.c0 &&
                       in.escaped_frac <= in.escaped_frac_max;
    if (in.require_u_max_positive)
    {
        g.flow_sanity_ok = g.flow_sanity_ok && in.u_max > TinyReal;
    }

    g.thermal_finite_ok = std::isfinite(in.t_mean) && std::isfinite(in.t_max);
    const bool energy_closed =
        !in.energy_residual_available ||
        (std::isfinite(in.cumulative_energy_residual_rel) &&
         in.cumulative_energy_residual_rel < kFrenchValidationEnergyRelTol && !in.ledger_invalid_accumulate);
    // Smoke: finite T. Validation: applied-energy identity when the ledger ran.
    g.thermal_energy_ok = g.thermal_finite_ok && energy_closed;

    // Doc 47 §5: age/lateness bounded by interval, plus a small discrete-step margin
    // (advection dt can push age slightly past interval; Tiny/1e-3 alone is too tight).
    const Real step_margin = std::max(kFrenchValidationEndTimeTol, Real(0.01) * in.em_interval);
    const Real age_budget = in.em_interval + step_margin;
    const bool age_ok = in.max_field_age_s < Real(0) ||
                        (std::isfinite(in.max_field_age_s) && in.max_field_age_s <= age_budget);
    const bool late_ok = in.max_nominal_lateness_s < Real(0) ||
                         (std::isfinite(in.max_nominal_lateness_s) && in.max_nominal_lateness_s <= age_budget);
    const bool sticky_solver_ok = in.refresh_solver_fail_count == 0;
    g.freshness_ok = !in.periodic_coupling ||
                     (in.freshness_measured && age_ok && late_ok && sticky_solver_ok);

    if (!g.startup_ok)
        appendGateFail(g.fail_reasons, "startup");
    if (!g.completed_requested_time)
        appendGateFail(g.fail_reasons, "incomplete_time");
    if (in.periodic_coupling && !g.periodic_refresh_ok)
        appendGateFail(g.fail_reasons, "periodic_refresh");
    if (!g.em_solver_ok)
        appendGateFail(g.fail_reasons, "em_solver");
    if (!g.q_transfer_ok)
        appendGateFail(g.fail_reasons, "q_transfer");
    if (in.periodic_coupling && !transfer_criteria_measured)
        appendGateFail(g.fail_reasons, "transfer_unmeasured");
    if (in.periodic_coupling && !g.freshness_ok)
        appendGateFail(g.fail_reasons, "freshness");
    if (!g.flow_sanity_ok)
        appendGateFail(g.fail_reasons, "flow");
    if (!g.thermal_finite_ok)
        appendGateFail(g.fail_reasons, "thermal_finite");
    if (in.energy_residual_available && !energy_closed)
        appendGateFail(g.fail_reasons, "energy_residual");
    if (in.ledger_invalid_accumulate)
        appendGateFail(g.fail_reasons, "ledger_invalid");
    if (!in.constant_thermal_rheology)
        appendGateFail(g.fail_reasons, "variable_cp_not_validatable");

    g.smoke_passed = g.startup_ok && g.em_solver_ok && g.flow_sanity_ok && g.thermal_finite_ok && q_l2_ok &&
                     (!in.out_of_grid_measured || in.n_out_of_grid_init < Real(1));
    g.validation_exercised = g.completed_requested_time && enough_refreshes && in.energy_residual_available;
    g.validation_passed = g.smoke_passed && g.validation_exercised && g.q_transfer_ok && g.thermal_energy_ok &&
                          transfer_criteria_measured && g.freshness_ok && !in.budget_exhausted &&
                          in.constant_thermal_rheology && !in.ledger_invalid_accumulate;
}

/** Instantaneous residual using the heat actually sampled onto the fluid (doc 40 R01). */
inline Real frenchEnergyResidual(Real dE_dt, Real joule_power_applied_fluid, Real wall_loss, Real free_loss,
                                 Real bottom_loss)
{
    return dE_dt - (joule_power_applied_fluid - wall_loss - free_loss - bottom_loss);
}

inline bool verifyFrenchAppliedEnergyResidualIdentity(std::string &message)
{
    FrenchAppliedEnergyLedger ledger;
    ledger.setH0(Real(0));
    ledger.e_q_applied = 120.0;
    ledger.e_boundary_out = 20.0;
    if (std::abs(ledger.residualH(Real(100))) > Real(1.0e-9))
    {
        message = "applied_energy_residual_identity_failed";
        return false;
    }
    const Real inst = frenchEnergyResidual(Real(80), Real(100), Real(15), Real(5), Real(0));
    if (std::abs(inst) > Real(1.0e-9))
    {
        message = "instantaneous_applied_residual_identity_failed";
        return false;
    }
    message = "applied_energy_residual_identity_ok";
    return true;
}

/** Doc 45/47 F03: compiled-in negative tests for the shared gate (must fail false-pass paths). */
inline bool verifyFrenchCouplingGateNegativeTests(std::string &message)
{
    auto base_ok = []() {
        FrenchCouplingGateInputs in;
        in.periodic_coupling = true;
        in.constant_thermal_rheology = true;
        in.energy_residual_available = true;
        in.fluid_power_recorded = true;
        in.out_of_grid_measured = true;
        in.end_time = Real(30);
        in.em_interval = Real(0.5);
        in.physical_time = Real(30);
        in.p_em = Real(60000);
        in.p_heat_fluid = Real(59000);
        in.phi_eq = Real(1.0e-4);
        in.q_sample_rel_l2 = Real(0.01);
        in.u_max = Real(0.1);
        in.c0 = Real(10);
        in.t_mean = Real(1473);
        in.t_max = Real(1480);
        in.cumulative_energy_residual_rel = Real(0.01);
        in.update_count = 61;
        // Exclusive-endpoint: end=30, interval=0.5 → expected 59 (t=0.5..29.5).
        in.periodic_refresh_count = 59;
        in.freshness_measured = true;
        in.max_field_age_s = Real(0.5);
        in.max_nominal_lateness_s = Real(0);
        in.refresh_solver_fail_count = 0;
        return in;
    };

    {
        FrenchCouplingGateInputs in = base_ok();
        in.periodic_refresh_count = 3;
        in.update_count = 4;
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed)
        {
            message = "gate_false_pass_under_refresh_schedule";
            return false;
        }
    }
    {
        // Doc 47: missing interior tick must fail (no blanket slack).
        FrenchCouplingGateInputs in = base_ok();
        in.end_time = Real(30.25);
        in.periodic_refresh_count = 59; // exclusive expected = 60
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed || g.periodic_refresh_ok)
        {
            message = "gate_false_pass_missing_interior_refresh";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        in.p_heat_fluid = Real(0);
        in.fluid_power_recorded = true;
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed)
        {
            message = "gate_false_pass_zero_fluid_heat";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        in.out_of_grid_measured = false;
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed)
        {
            message = "gate_false_pass_out_of_grid_unmeasured";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        in.q_sample_rel_l2 = Real(-1);
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed)
        {
            message = "gate_false_pass_q_l2_unmeasured";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        in.physical_time = Real(29.7); // was accepted under 99% horizon
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed || g.completed_requested_time)
        {
            message = "gate_false_pass_incomplete_end_time";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        in.ledger_invalid_accumulate = true;
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed)
        {
            message = "gate_false_pass_sticky_ledger_invalid";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        in.freshness_measured = false;
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed)
        {
            message = "gate_false_pass_freshness_unmeasured";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        in.max_field_age_s = Real(2.0); // > interval+tol
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed || g.freshness_ok)
        {
            message = "gate_false_pass_stale_field";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        in.refresh_solver_fail_count = 1;
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.validation_passed)
        {
            message = "gate_false_pass_sticky_refresh_solver_fail";
            return false;
        }
    }
    {
        FrenchCouplingGateInputs in = base_ok();
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (!g.validation_passed)
        {
            message = "gate_true_pass_baseline_failed:" + g.fail_reasons;
            return false;
        }
        if (g.expected_periodic_refreshes != 59)
        {
            message = "gate_expected_refresh_count_not_exclusive_endpoint";
            return false;
        }
    }
    {
        // end=2, interval=0.5 → expected 3 (matches Batch A/B smoke).
        FrenchCouplingGateInputs in = base_ok();
        in.end_time = Real(2);
        in.physical_time = Real(2);
        in.periodic_refresh_count = 3;
        in.update_count = 4;
        in.max_field_age_s = Real(0.5);
        FrenchCouplingRunGates g;
        evaluateFrenchCouplingRunGates(in, g);
        if (g.expected_periodic_refreshes != 3 || !g.validation_passed)
        {
            message = "gate_smoke_end2_exclusive_failed:" + g.fail_reasons;
            return false;
        }
    }
    {
        // Doc 54: nominal next-due schedule (no last+interval drift).
        ThermoEMCouplingState sched;
        armFrenchPeriodicRefreshSchedule(sched, Real(0.5));
        if (sched.next_due_time < Real(0.49) || sched.next_due_time > Real(0.51))
        {
            message = "schedule_arm_next_due_failed";
            return false;
        }
        if (!shouldUpdateEM(Real(0.5), Real(0.5), sched.next_due_time))
        {
            message = "schedule_should_fire_at_due_failed";
            return false;
        }
        if (shouldUpdateEM(Real(0.49), Real(0.5), sched.next_due_time))
        {
            message = "schedule_premature_fire_failed";
            return false;
        }
        noteFrenchPeriodicRefresh(sched, Real(0.52), Real(0.5), true);
        if (sched.max_nominal_lateness_s < Real(0.019) || sched.max_field_age_s < Real(0.51))
        {
            message = "schedule_lateness_or_age_failed";
            return false;
        }
        if (sched.next_due_time < Real(0.99) || sched.next_due_time > Real(1.01))
        {
            message = "schedule_advance_next_due_failed";
            return false;
        }
        if (shouldUpdateEM(Real(0.75), Real(0.5), sched.next_due_time))
        {
            message = "schedule_mid_interval_fire_failed";
            return false;
        }
    }
    message = "gate_negative_tests_ok";
    return true;
}

struct FrenchEnergyBudgetRowExtras
{
    Real p_em_recon = 0.0;
    Real p_heat_fluid = 0.0;
    Real e_q_applied = 0.0;
    Real e_boundary_out = 0.0;
    Real residual_H = 0.0;
    Real residual_H_rel = 0.0;
    Real residual_H_rel_vs_net = 0.0;
    Real sum_mass = 0.0;
    Real sum_vol = 0.0;
    Real sum_rho0_vol = 0.0;
    /** Doc 68 F07: denser mass/ρ/V + parallel mass-based H (diagnostic only). */
    Real sum_rho_vol = 0.0;
    Real rho_min = 0.0;
    Real rho_max = 0.0;
    Real rho_mean = 0.0;
    Real max_abs_mass_minus_rho_vol = 0.0;
    Real rms_mass_minus_rho_vol = 0.0;
    Real rel_l2_mass_vs_rho_vol = 0.0;
    Real thermal_energy_mass = 0.0; // H_m = Σ m cp ΔT; ledger still uses rho0·V H
};

inline void writeFrenchEnergyBudgetCsvHeader(const std::string &path)
{
    std::ofstream out(path);
    out << "time,thermal_energy,dE_dt,joule_power,wall_heat_loss,free_surface_heat_loss,bottom_heat_loss,"
           "energy_residual,energy_residual_rel,T_min,T_mean,T_max,U_max,U_buoy,thermal_diffusion,"
           "mu_min,mu_mean,mu_max,cp_min,cp_mean,cp_max,k_min,k_mean,k_max,constant_thermal_rheology,"
           "p_em_recon,p_heat_fluid,e_q_applied,e_boundary_out,energy_residual_H,energy_residual_H_rel,"
           "energy_residual_H_rel_vs_net,sum_mass,sum_vol,sum_rho0_vol,sum_rho_vol,rho_min,rho_max,rho_mean,"
           "max_abs_mass_minus_rho_vol,rms_mass_minus_rho_vol,rel_l2_mass_vs_rho_vol,thermal_energy_mass,schema\n";
}

inline void appendFrenchEnergyBudgetCsv(const std::string &path, Real time, Real thermal_energy, Real dE_dt,
                                        Real joule_power, Real wall_loss, Real free_loss, Real bottom_loss,
                                        Real t_min = 0.0, Real t_mean = 0.0, Real t_max = 0.0, Real u_max = 0.0,
                                        Real u_buoy = 0.0, int thermal_diffusion = 0,
                                        const OphelieThermalRheologyStats *rheology = nullptr,
                                        int constant_thermal_rheology = 1,
                                        const FrenchEnergyBudgetRowExtras *extras = nullptr)
{
    const Real residual = frenchEnergyResidual(dE_dt, joule_power, wall_loss, free_loss, bottom_loss);
    const Real residual_rel = residual / (std::abs(joule_power) + TinyReal);
    OphelieThermalRheologyStats zero_stats;
    const OphelieThermalRheologyStats &r = rheology != nullptr ? *rheology : zero_stats;
    FrenchEnergyBudgetRowExtras z;
    const FrenchEnergyBudgetRowExtras &x = extras != nullptr ? *extras : z;
    std::ofstream out(path, std::ios::app);
    out << time << "," << thermal_energy << "," << dE_dt << "," << joule_power << "," << wall_loss << "," << free_loss
        << "," << bottom_loss << "," << residual << "," << residual_rel << "," << t_min << "," << t_mean << ","
        << t_max << "," << u_max << "," << u_buoy << "," << thermal_diffusion << "," << r.mu_min << "," << r.mu_mean
        << "," << r.mu_max << "," << r.cp_min << "," << r.cp_mean << "," << r.cp_max << "," << r.k_min << ","
        << r.k_mean << "," << r.k_max << "," << constant_thermal_rheology << "," << x.p_em_recon << ","
        << x.p_heat_fluid << "," << x.e_q_applied << "," << x.e_boundary_out << "," << x.residual_H << ","
        << x.residual_H_rel << "," << x.residual_H_rel_vs_net << "," << x.sum_mass << "," << x.sum_vol << ","
        << x.sum_rho0_vol << "," << x.sum_rho_vol << "," << x.rho_min << "," << x.rho_max << "," << x.rho_mean << ","
        << x.max_abs_mass_minus_rho_vol << "," << x.rms_mass_minus_rho_vol << "," << x.rel_l2_mass_vs_rho_vol << ","
        << x.thermal_energy_mass << ",5\n";
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
