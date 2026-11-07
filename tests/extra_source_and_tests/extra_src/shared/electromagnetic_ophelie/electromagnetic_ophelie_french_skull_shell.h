#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_SKULL_SHELL_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_SKULL_SHELL_H

/**
 * Round 3B Level-2: fixed solid skull shell parameters and CLI.
 *
 * The resolved conduction body is implemented on the existing WallBoundary SolidBody
 * (see docs/ophelie/17_ROUND3B_LEVEL2_FIXED_SKULL_SHELL_PLAN.md). This header owns
 * parameters, CLI parsing, and mutual-exclusion guards vs Level-1 effective resistance.
 *
 * OPEN_PARAMETER: delta_m, k_skull (frozen-glass ~1 W/(m·K)).
 */

#include "electromagnetic_ophelie_french_skull_resistance.h"

#include <cstring>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct FrenchSkullShellParams
{
    bool enable = false;
    /** Shell thickness [m]. Also used to size WallBoundary when shell path is active. */
    Real delta_m = Real(0.005);
    Real k_skull = Real(1.0);
    Real h_cool_side = Real(300);
    Real h_cool_bottom = Real(35);
    Real t_cool = Real(300);
    Real rho_skull = Real(2300); // frozen glass density OPEN_PARAMETER
    Real cp_skull = Real(1100);  // near low-T Fig. 2 OPEN_PARAMETER
};

inline bool isFrenchSkullShellCliOption(const char *arg)
{
    return std::strcmp(arg, "--skull-shell") == 0 || std::strncmp(arg, "--skull-shell-delta=", 20) == 0 ||
           std::strncmp(arg, "--skull-shell-k=", 16) == 0 ||
           std::strncmp(arg, "--skull-shell-h-cool-side=", 26) == 0 ||
           std::strncmp(arg, "--skull-shell-h-cool-bottom=", 28) == 0 ||
           std::strncmp(arg, "--skull-shell-t-cool=", 21) == 0;
}

inline bool tryApplyFrenchSkullShellCli(const char *arg, FrenchSkullShellParams &shell)
{
    if (std::strcmp(arg, "--skull-shell") == 0)
    {
        shell.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-shell-delta=", 20) == 0)
    {
        shell.delta_m = static_cast<Real>(std::atof(arg + 20));
        shell.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-shell-k=", 16) == 0)
    {
        shell.k_skull = static_cast<Real>(std::atof(arg + 16));
        shell.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-shell-h-cool-side=", 26) == 0)
    {
        shell.h_cool_side = static_cast<Real>(std::atof(arg + 26));
        shell.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-shell-h-cool-bottom=", 28) == 0)
    {
        shell.h_cool_bottom = static_cast<Real>(std::atof(arg + 28));
        shell.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-shell-t-cool=", 21) == 0)
    {
        shell.t_cool = static_cast<Real>(std::atof(arg + 21));
        shell.enable = true;
        return true;
    }
    return false;
}

/** Reject Level-1 effective skull / legacy balance when Level-2 shell is requested. */
inline bool assertFrenchSkullShellExclusive(const FrenchSkullShellParams &shell,
                                            const FrenchSkullResistanceParams &level1, bool balance_heat_loss,
                                            std::string &error)
{
    if (!shell.enable)
    {
        return true;
    }
    if (level1.enable)
    {
        error = "refuse --skull-shell together with --skull-resistance (Level-1). Use shell alone.";
        return false;
    }
    if (balance_heat_loss)
    {
        error = "refuse --skull-shell together with --balance-heat-loss.";
        return false;
    }
    if (!(shell.delta_m > TinyReal && shell.k_skull > TinyReal && shell.h_cool_side > TinyReal &&
          shell.h_cool_bottom > TinyReal && shell.t_cool > TinyReal))
    {
        error = "invalid --skull-shell delta/k/h_cool/t_cool";
        return false;
    }
    return true;
}

inline void logFrenchSkullShellParams(const FrenchSkullShellParams &shell, Real wall_thickness_used)
{
    if (!shell.enable)
    {
        return;
    }
    std::cout << "[ophelie][skull-shell] Level-2 fixed solid skull shell (WallBoundary conduction path)\n"
              << "  delta=" << shell.delta_m << " m  wall_thickness_used=" << wall_thickness_used
              << " m  k_skull=" << shell.k_skull << " W/m/K\n"
              << "  h_cool side/bottom=" << shell.h_cool_side << "/" << shell.h_cool_bottom
              << "  T_cool=" << shell.t_cool << " K\n"
              << "  rho/cp OPEN_PARAMETER=" << shell.rho_skull << "/" << shell.cp_skull << "\n"
              << "  note: melt side/bottom Level-1 resistance is off; free-crust may remain\n";
}

/**
 * Host Ra / Pr / velocity-scale estimate for Round 3C detectability (natural case).
 *
 * Scales:
 * - U_inertial ~ sqrt(g β ΔT H)  (low-Pr free-fall)
 * - U_viscous  ~ g β ΔT H² / ν   (high-Pr Stokes / creeping buoyancy)
 *
 * Compare U_viscous to the known WCSPH residual floor (~2 cm/s at dp=15 mm) before
 * deciding on a finer natural mesh.
 *
 * @param mu_hot optional μ(T0+ΔT); if >0, prints μ contrast for μ(T) awareness.
 * @param u_noise_hint documented residual U_rms floor from prior natural runs.
 */
inline void logFrenchNaturalRaPrEstimate(Real radius, Real height, Real delta_t, Real rho, Real mu, Real cp, Real k,
                                         Real beta, Real g = Real(9.81), Real mu_hot = Real(-1),
                                         Real u_noise_hint = Real(0.02))
{
    const Real dT = std::abs(delta_t);
    const Real H = std::max(height, TinyReal);
    const Real nu = mu / (rho + TinyReal);
    const Real alpha = k / (rho * cp + TinyReal);
    const Real ra = g * std::abs(beta) * dT * H * H * H / (nu * alpha + TinyReal);
    const Real pr = nu / (alpha + TinyReal);
    const Real u_inertial = std::sqrt(std::max(g * std::abs(beta) * dT * H, Real(0)));
    const Real u_viscous = g * std::abs(beta) * dT * H * H / (nu + TinyReal);
    const Real ratio = u_viscous / (u_noise_hint + TinyReal);
    const char *verdict = (ratio >= Real(1)) ? "signal_may_exceed_noise"
                          : (ratio >= Real(0.2)) ? "marginal"
                                                 : "buried_defer_finer_mesh";
    std::cout << "[ophelie][natural-ra] R=" << radius << " H=" << H << " dT_proxy=" << dT
              << " mu=" << mu << " nu=" << nu << " alpha=" << alpha
              << " Ra=" << ra << " Pr=" << pr
              << " U_inertial~" << u_inertial << " U_viscous~" << u_viscous << " m/s"
              << " U_noise_hint=" << u_noise_hint
              << " U_visc/U_noise=" << ratio
              << " verdict=" << verdict;
    if (mu_hot > TinyReal)
    {
        std::cout << " mu(T0+dT)=" << mu_hot << " mu_ratio=" << (mu / (mu_hot + TinyReal));
    }
    std::cout << " (OPEN: dT_proxy not paper figure ΔT)\n";
}

/**
 * Effective wall conductivity so that a hydro wall of thickness L has the same
 * conduction resistance as thermal OPEN_PARAMETER delta / k_skull:
 *   R = L / k_eff = delta / k_skull  ⇒  k_eff = k_skull * L / delta
 */
inline Real frenchSkullShellEffectiveWallConductivity(Real k_skull, Real wall_thickness, Real delta_m)
{
    const Real delta = std::max(delta_m, TinyReal);
    return k_skull * std::max(wall_thickness, TinyReal) / delta;
}

/**
 * Initialize wall temperature with a linear profile from melt-side ~T_melt to outer T_cool,
 * so Contact sees an immediate radial gradient on coarse (thick hydro) walls.
 */
inline void hostInitializeFrenchSkullShellWallTemperatureProfile(BaseParticles &wall_particles,
                                                                 const OphelieFrenchReducedCaseParams &french,
                                                                 Real wall_thickness, Real t_melt, Real t_cool,
                                                                 Real t_initial)
{
    const size_t n = wall_particles.TotalRealParticles();
    if (n == 0)
    {
        return;
    }
    syncVariableToHost<Vecd>(wall_particles, "Position");
    Real *temperature = wall_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = wall_particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Vecd *pos = wall_particles.getVariableDataByName<Vecd>("Position");
    const Vecd &c = french.glass_center;
    const Real z_lo = c[2] - french.glass_half_height;
    const Real L = std::max(wall_thickness, TinyReal);
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - c[0];
        const Real dy = pos[i][1] - c[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real dist_side = std::max(r - french.glass_radius, Real(0));
        const Real dist_bottom = std::max(z_lo - pos[i][2], Real(0));
        const Real dist = std::max(dist_side, dist_bottom);
        const Real xi = std::min(dist / L, Real(1));
        const Real t = t_melt + (t_cool - t_melt) * xi;
        temperature[i] = t;
        delta_t[i] = t - t_initial;
    }
    syncVariableToDevice<Real>(wall_particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(wall_particles, kOphelieThermalDeltaTField);
}

/**
 * Dirichlet T_cool on the coolant-facing outer layer of WallBoundary particles
 * (outer ~1.5 dp of the wall thickness). Updates Temperature and OphelieThermalDeltaT.
 */
inline size_t hostApplyFrenchSkullShellOuterCoolDirichlet(BaseParticles &wall_particles,
                                                          const OphelieFrenchReducedCaseParams &french,
                                                          Real wall_thickness, Real t_cool, Real t_initial)
{
    const size_t n = wall_particles.TotalRealParticles();
    if (n == 0)
    {
        return 0;
    }
    syncVariableToHost<Vecd>(wall_particles, "Position");
    syncVariableToHost<Real>(wall_particles, kOphelieTemperatureField);
    syncVariableToHost<Real>(wall_particles, kOphelieThermalDeltaTField);
    const Vecd *pos = wall_particles.getVariableDataByName<Vecd>("Position");
    Real *temperature = wall_particles.getVariableDataByName<Real>(kOphelieTemperatureField);
    Real *delta_t = wall_particles.getVariableDataByName<Real>(kOphelieThermalDeltaTField);
    const Vecd &c = french.glass_center;
    const Real z_lo = c[2] - french.glass_half_height;
    const Real outer_band = std::max(Real(1.5) * french.dp, Real(0.25) * wall_thickness);
    size_t n_cool = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - c[0];
        const Real dy = pos[i][1] - c[1];
        const Real r = std::sqrt(dx * dx + dy * dy);
        const Real dist_from_inner_side = r - french.glass_radius;
        const Real dist_from_inner_bottom = z_lo - pos[i][2];
        const Real dist_from_inner = std::max(dist_from_inner_side, dist_from_inner_bottom);
        if (dist_from_inner < wall_thickness - outer_band)
        {
            continue;
        }
        temperature[i] = t_cool;
        delta_t[i] = t_cool - t_initial;
        ++n_cool;
    }
    syncVariableToDevice<Real>(wall_particles, kOphelieTemperatureField);
    syncVariableToDevice<Real>(wall_particles, kOphelieThermalDeltaTField);
    return n_cool;
}

/**
 * Instantaneous Contact thermal power [W] from DiffusionRelaxationCK transfer field.
 * Field name is TemperatureTransferWith{SourceBodyName} on the source body particles.
 * Convention: P = sum(transfer_i * V_i) is dE/dt of that body from Contact (negative ⇒ heat leaving).
 */
inline Real hostOphelieContactThermalPowerW(BaseParticles &particles, const std::string &transfer_field)
{
    syncVariableToHost<Real>(particles, transfer_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *transfer = particles.getVariableDataByName<Real>(transfer_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();
    Real power = Real(0);
    for (size_t i = 0; i < n; ++i)
    {
        power += transfer[i] * vol[i];
    }
    return power;
}

inline std::string frenchSkullShellMeltContactTransferField(const std::string &melt_body_name)
{
    return std::string(kOphelieTemperatureField) + "TransferWith" + melt_body_name;
}

inline std::string frenchSkullShellWallContactTransferField(const std::string &wall_body_name)
{
    return std::string(kOphelieTemperatureField) + "TransferWith" + wall_body_name;
}

/**
 * Melt→wall contact heat-loss power [W] (>=0 when melt loses heat).
 * Also returns wall-side contact dE/dt for consistency check.
 */
inline void hostFrenchSkullShellContactHeatPowers(BaseParticles &melt_particles, BaseParticles &wall_particles,
                                                  const std::string &melt_body_name,
                                                  const std::string &wall_body_name, Real &melt_loss_w,
                                                  Real &wall_gain_w)
{
    const Real melt_dE_dt =
        hostOphelieContactThermalPowerW(melt_particles, frenchSkullShellMeltContactTransferField(melt_body_name));
    const Real wall_dE_dt =
        hostOphelieContactThermalPowerW(wall_particles, frenchSkullShellWallContactTransferField(wall_body_name));
    melt_loss_w = -melt_dE_dt; // heat leaving melt
    wall_gain_w = wall_dE_dt;  // heat entering wall
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_SKULL_SHELL_H
