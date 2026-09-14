#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_LITERATURE_PARAMETERS_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_LITERATURE_PARAMETERS_H

#include "electromagnetic_ophelie_french_material_laws.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_stirring_geometry.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * Vertical cylinder axis. CAD for stirring was delivered as -y; meshes are rotated to +z.
 * Do not assume z without checking this field.
 */
enum class FrenchVerticalAxis
{
    X = 0,
    Y = 1,
    Z = 2
};

enum class FrenchLiteratureRegime
{
    NaturalConvection,
    MechanicallyStirred
};

enum class EMControlMode
{
    FixedCoilCurrent,
    FixedAbsorbedPower
};

enum class ThermoEMCouplingMode
{
    Off,
    Periodic
};

/**
 * Unified French-literature knobs. Existing case structs remain the storage for geometry;
 * this record documents the two regimes and OPEN_PARAMETER items that must not be silently
 * overwritten.
 *
 * OPEN_PARAMETER (do not guess away):
 * - Natural frequency in the frozen-Q case is 282 kHz, not 300 kHz.
 * - Natural melt H=0.185 m vs stirred H=0.21 m (different papers / EREBUS vs CEP natural).
 * - rho0: Table-1 2750 vs stirred-run docs 2800 kg/m^3.
 * - sigma(1473 K): Table-1 16 S/m vs thesis IV.10 ~24.75 S/m.
 * - Generator 400 kW is supply-side, never a Joule calibration target.
 * - Np characteristic diameter / exact Nu definition: not frozen here.
 */
struct FrenchLiteratureParameters
{
    FrenchLiteratureRegime regime = FrenchLiteratureRegime::NaturalConvection;
    FrenchVerticalAxis vertical_axis = FrenchVerticalAxis::Z;

    Real frequency_hz = 300000.0;
    Real t_initial_k = 1473.0;
    Real generator_power_w = 400000.0;          // supply-side; not a calibration target
    Real target_glass_absorbed_power_w = 50000.0;
    Real stirrer_rpm = 0.0;

    EMControlMode em_control = EMControlMode::FixedCoilCurrent;
    ThermoEMCouplingMode coupling = ThermoEMCouplingMode::Off;
    Real em_update_interval_s = 10.0;
    Real sigma_under_relaxation = 0.3;
    bool enable_aind = false;
    bool enable_thermal_diffusion = false;
    bool q_conservative_remap = false;
    bool constant_sigma = false;
    bool zero_q = false;
    bool enable_boussinesq = true;

    OphelieFrenchGlassMaterialLaws material = makeFrenchNaturalCep2008MaterialLaws();
};

inline const char *frenchLiteratureRegimeName(FrenchLiteratureRegime regime)
{
    return regime == FrenchLiteratureRegime::MechanicallyStirred ? "mechanically_stirred" : "natural_convection";
}

inline const char *emControlModeName(EMControlMode mode)
{
    return mode == EMControlMode::FixedAbsorbedPower ? "fixed-power" : "fixed-current";
}

inline const char *thermoEMCouplingModeName(ThermoEMCouplingMode mode)
{
    return mode == ThermoEMCouplingMode::Periodic ? "periodic" : "off";
}

inline int frenchVerticalAxisIndex(FrenchVerticalAxis axis)
{
    return static_cast<int>(axis);
}

inline FrenchLiteratureParameters makeFrenchNaturalLiteratureParameters()
{
    FrenchLiteratureParameters p;
    p.regime = FrenchLiteratureRegime::NaturalConvection;
    p.vertical_axis = FrenchVerticalAxis::Z;
    p.frequency_hz = 282000.0; // OPEN_PARAMETER vs 300 kHz
    p.target_glass_absorbed_power_w = 50000.0;
    p.stirrer_rpm = 0.0;
    p.material = makeFrenchNaturalCep2008MaterialLaws();
    p.em_control = EMControlMode::FixedCoilCurrent;
    p.enable_aind = false;
    return p;
}

inline FrenchLiteratureParameters makeFrenchStirredLiteratureParameters()
{
    FrenchLiteratureParameters p;
    p.regime = FrenchLiteratureRegime::MechanicallyStirred;
    p.vertical_axis = FrenchVerticalAxis::Z; // CAD rotated; see french_stirring_geometry.h
    p.frequency_hz = 300000.0;
    p.target_glass_absorbed_power_w = 60000.0;
    p.stirrer_rpm = 10.0;
    p.material = makeJacoutotThesisIV10MaterialLaws();
    p.em_control = EMControlMode::FixedCoilCurrent;
    p.enable_aind = false;
    return p;
}

inline Real frenchBuoyancySpeedScale(Real gravity, Real beta, Real delta_t, Real height)
{
    const Real g_beta_dt_h = gravity * std::abs(beta) * std::abs(delta_t) * std::max(height, Real(0));
    return std::sqrt(std::max(g_beta_dt_h, Real(0)));
}

inline bool isFrenchLiteratureCouplingCliOption(const char *arg)
{
    return std::strncmp(arg, "--thermo-em-coupling=", 21) == 0 ||
           std::strncmp(arg, "--em-update-interval=", 21) == 0 ||
           std::strncmp(arg, "--sigma-under-relaxation=", 25) == 0 ||
           std::strncmp(arg, "--em-control=", 13) == 0 ||
           std::strncmp(arg, "--target-absorbed-power=", 24) == 0 || std::strncmp(arg, "--aind=", 7) == 0 ||
           std::strcmp(arg, "--enable-thermal-diffusion") == 0 || std::strcmp(arg, "--no-thermal-diffusion") == 0 ||
           std::strcmp(arg, "--constant-sigma") == 0 || std::strcmp(arg, "--zero-q") == 0 ||
           std::strcmp(arg, "--q-conservative-remap") == 0 || std::strcmp(arg, "--no-boussinesq") == 0;
}

/**
 * Shared coupling CLI. Returns true if the argument was consumed.
 * `thermal_diffusion_explicit` is set when the user passed enable/disable so periodic
 * mode can default-on conduction without overriding an explicit `--no-thermal-diffusion`.
 */
inline bool tryApplyFrenchLiteratureCli(const char *arg, FrenchLiteratureParameters &p,
                                        bool &thermal_diffusion_explicit)
{
    if (std::strncmp(arg, "--thermo-em-coupling=", 21) == 0)
    {
        const std::string mode(arg + 21);
        p.coupling = (mode == "periodic") ? ThermoEMCouplingMode::Periodic : ThermoEMCouplingMode::Off;
        return true;
    }
    if (std::strncmp(arg, "--em-update-interval=", 21) == 0)
    {
        p.em_update_interval_s = static_cast<Real>(std::atof(arg + 21));
        return true;
    }
    if (std::strncmp(arg, "--sigma-under-relaxation=", 25) == 0)
    {
        p.sigma_under_relaxation = static_cast<Real>(std::atof(arg + 25));
        return true;
    }
    if (std::strncmp(arg, "--em-control=", 13) == 0)
    {
        const std::string mode(arg + 13);
        p.em_control = (mode == "fixed-power") ? EMControlMode::FixedAbsorbedPower : EMControlMode::FixedCoilCurrent;
        return true;
    }
    if (std::strncmp(arg, "--target-absorbed-power=", 24) == 0)
    {
        p.target_glass_absorbed_power_w = static_cast<Real>(std::atof(arg + 24));
        return true;
    }
    if (std::strncmp(arg, "--aind=", 7) == 0)
    {
        const std::string mode(arg + 7);
        p.enable_aind = (mode != "off");
        return true;
    }
    if (std::strcmp(arg, "--enable-thermal-diffusion") == 0)
    {
        p.enable_thermal_diffusion = true;
        thermal_diffusion_explicit = true;
        return true;
    }
    if (std::strcmp(arg, "--no-thermal-diffusion") == 0)
    {
        p.enable_thermal_diffusion = false;
        thermal_diffusion_explicit = true;
        return true;
    }
    if (std::strcmp(arg, "--constant-sigma") == 0)
    {
        p.constant_sigma = true;
        return true;
    }
    if (std::strcmp(arg, "--zero-q") == 0)
    {
        p.zero_q = true;
        return true;
    }
    if (std::strcmp(arg, "--q-conservative-remap") == 0)
    {
        p.q_conservative_remap = true;
        return true;
    }
    if (std::strcmp(arg, "--no-boussinesq") == 0)
    {
        p.enable_boussinesq = false;
        return true;
    }
    return false;
}

/** Periodic coupling turns SPH conduction on unless the user passed `--no-thermal-diffusion`. */
inline void finalizeFrenchPeriodicLiteratureDefaults(FrenchLiteratureParameters &p, bool thermal_diffusion_explicit)
{
    if (p.coupling == ThermoEMCouplingMode::Periodic && !thermal_diffusion_explicit)
    {
        p.enable_thermal_diffusion = true;
    }
}

inline void logFrenchLiteratureParameters(const FrenchLiteratureParameters &p)
{
    std::cout << "[ophelie][french-lit] regime=" << frenchLiteratureRegimeName(p.regime)
              << " vertical_axis=" << frenchVerticalAxisIndex(p.vertical_axis)
              << " f_hz=" << p.frequency_hz << " T0=" << p.t_initial_k
              << " generator_power=" << p.generator_power_w
              << " target_glass_absorbed_power=" << p.target_glass_absorbed_power_w
              << " stirrer_rpm=" << p.stirrer_rpm << " em_control=" << emControlModeName(p.em_control)
              << " coupling=" << thermoEMCouplingModeName(p.coupling)
              << " em_update_interval_s=" << p.em_update_interval_s
              << " sigma_under_relaxation=" << p.sigma_under_relaxation
              << " aind=" << (p.enable_aind ? "on" : "off")
              << " thermal_diffusion=" << (p.enable_thermal_diffusion ? 1 : 0)
              << " q_conservative_remap=" << (p.q_conservative_remap ? 1 : 0) << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_LITERATURE_PARAMETERS_H
