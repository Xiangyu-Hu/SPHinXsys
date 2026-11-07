#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_SKULL_RESISTANCE_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_SKULL_RESISTANCE_H

/**
 * Round 3B Level 1: effective fixed skull thermal resistance (not a resolved skull shape).
 *
 * Series model at cooled side/bottom walls:
 *   R_total = delta_skull / k_skull + 1 / h_cool
 *   h_eff   = 1 / R_total
 *   q_wall  = h_eff * (T_melt - T_cool)
 *
 * Sources / OPEN_PARAMETER notes:
 * - Skull thickness: Jacoutot CES/CEP — millimetre to centimetre scale (OPEN_PARAMETER).
 * - k_skull: frozen-glass conductivity ~1 W/(m·K) near low-T Fig. 2 NUM (OPEN_PARAMETER).
 * - h_cool: coolant-side coefficient (metal + water). CEP sensitivity often 100–300 W/(m²·K)
 *   on the side; bottom skull formation cited above ~35 W/(m²·K). Defaults keep the current
 *   production h_side/h_bottom as the coolant-side values before adding skull conduction.
 *
 * Free-surface convection/radiation is left unchanged in Level-1 wall skull unless
 * `--free-crust-resistance` enables the Level-1b companion (conduction + surface conv/rad at T_s).
 * `--balance-heat-loss` remains an engineering legacy knob and must not be mixed into this path.
 */

#include "electromagnetic_ophelie_thermal_diffusion_one_way.h"

#include <cmath>
#include <cstring>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct FrenchSkullResistanceParams
{
    bool enable = false;
    /** Side skull thickness [m]. Default 5 mm (OPEN_PARAMETER within mm–cm). */
    Real delta_side_m = Real(0.005);
    /** Bottom skull thickness [m]. Default 5 mm. */
    Real delta_bottom_m = Real(0.005);
    /** Frozen-glass conductivity [W/(m·K)]. */
    Real k_skull = Real(1.0);
    /**
     * Coolant-side coefficients [W/(m²·K)] before skull conduction.
     * Negative means "use whatever is already on the thermal BC options".
     */
    Real h_cool_side = Real(-1);
    Real h_cool_bottom = Real(-1);
    Real t_cool = Real(-1); // <0 -> keep thermal BC t_cool
};

inline Real frenchSkullSeriesResistance(Real delta_m, Real k_skull, Real h_cool)
{
    const Real r_cond = std::max(delta_m, Real(0)) / (std::max(k_skull, TinyReal));
    const Real r_cool = Real(1) / (std::max(h_cool, TinyReal));
    return r_cond + r_cool;
}

inline Real frenchSkullEffectiveHeatTransfer(Real delta_m, Real k_skull, Real h_cool)
{
    return Real(1) / (frenchSkullSeriesResistance(delta_m, k_skull, h_cool) + TinyReal);
}

struct FrenchZeroDHeatBudget
{
    Real area_side = 0.0;
    Real area_bottom = 0.0;
    Real area_free = 0.0;
    Real h_side_eff = 0.0;
    Real h_bottom_eff = 0.0;
    Real h_free = 0.0;
    Real emissivity = 0.0;
    Real t_melt = 0.0;
    Real t_cool = 0.0;
    Real t_ambient = 0.0;
    Real t_rad_ambient = 0.0;
    Real stefan_boltzmann = Real(5.670374419e-8);
    Real side_loss_w = 0.0;
    Real bottom_loss_w = 0.0;
    Real free_conv_loss_w = 0.0;
    Real free_rad_loss_w = 0.0;
    Real total_loss_w = 0.0;
    Real target_joule_w = 0.0;
    Real loss_over_joule = 0.0;
};

/** Cylinder melt estimate: side = 2π R H, bottom = π R², free = π R². */
inline FrenchZeroDHeatBudget estimateFrenchZeroDHeatBudget(Real radius, Real height, Real t_melt,
                                                          const OphelieThermalDiffusionOneWayOptions &bc,
                                                          Real target_joule_w)
{
    FrenchZeroDHeatBudget b;
    b.area_side = Real(2) * Real(M_PI) * radius * height;
    b.area_bottom = Real(M_PI) * radius * radius;
    b.area_free = b.area_bottom;
    b.h_side_eff = bc.h_side;
    b.h_bottom_eff = bc.h_bottom;
    b.h_free = bc.h_free;
    b.emissivity = bc.emissivity;
    b.t_melt = t_melt;
    b.t_cool = bc.t_cool;
    b.t_ambient = bc.t_ambient;
    b.t_rad_ambient = bc.t_rad_ambient;
    b.stefan_boltzmann = bc.stefan_boltzmann;
    b.target_joule_w = target_joule_w;

    b.side_loss_w = b.h_side_eff * (t_melt - b.t_cool) * b.area_side;
    b.bottom_loss_w = b.h_bottom_eff * (t_melt - b.t_cool) * b.area_bottom;
    {
        Real q_total = Real(0);
        Real q_conv = Real(0);
        Real q_rad = Real(0);
        frenchFreeSurfaceHeatFlux(t_melt, b.h_free, b.emissivity, b.t_ambient, b.t_rad_ambient, b.stefan_boltzmann,
                                  bc.enable_free_crust ? Real(1) : Real(0), bc.free_crust_delta_m, bc.free_crust_k,
                                  q_total, q_conv, q_rad);
        b.free_conv_loss_w = q_conv * b.area_free;
        b.free_rad_loss_w = q_rad * b.area_free;
        b.total_loss_w = b.side_loss_w + b.bottom_loss_w + q_total * b.area_free;
    }
    b.loss_over_joule = b.total_loss_w / (std::abs(target_joule_w) + TinyReal);
    return b;
}

/**
 * Apply skull series resistance onto side/bottom h of `bc`.
 * Returns false if parameters are invalid. Set `log` to print the mapping.
 */
inline bool applyFrenchSkullResistanceToThermalBc(OphelieThermalDiffusionOneWayOptions &bc,
                                                  const FrenchSkullResistanceParams &skull, std::string &error,
                                                  bool log = true)
{
    if (!skull.enable)
    {
        return true;
    }
    if (!(skull.delta_side_m >= Real(0) && skull.delta_bottom_m >= Real(0) && skull.k_skull > TinyReal))
    {
        error = "invalid skull delta/k_skull";
        return false;
    }
    const Real h_cool_side = skull.h_cool_side > TinyReal ? skull.h_cool_side : bc.h_side;
    const Real h_cool_bottom = skull.h_cool_bottom > TinyReal ? skull.h_cool_bottom : bc.h_bottom;
    if (!(h_cool_side > TinyReal && h_cool_bottom > TinyReal))
    {
        error = "invalid h_cool for skull series model";
        return false;
    }
    if (skull.t_cool > TinyReal)
    {
        bc.t_cool = skull.t_cool;
    }

    const Real h_side_eff = frenchSkullEffectiveHeatTransfer(skull.delta_side_m, skull.k_skull, h_cool_side);
    const Real h_bottom_eff =
        frenchSkullEffectiveHeatTransfer(skull.delta_bottom_m, skull.k_skull, h_cool_bottom);
    if (log)
    {
        std::cout << "[ophelie][skull-resistance] Level-1 effective fixed skull (not resolved shape)\n"
                  << "  side:   delta=" << skull.delta_side_m << " m  k_skull=" << skull.k_skull
                  << " W/m/K  h_cool=" << h_cool_side << " -> h_eff=" << h_side_eff << " W/m^2/K\n"
                  << "  bottom: delta=" << skull.delta_bottom_m << " m  k_skull=" << skull.k_skull
                  << " W/m/K  h_cool=" << h_cool_bottom << " -> h_eff=" << h_bottom_eff << " W/m^2/K\n"
                  << "  free surface: use --free-crust-resistance for Level-1b companion\n";
    }
    bc.h_side = h_side_eff;
    bc.h_bottom = h_bottom_eff;
    return true;
}

/**
 * Level-1b: enable free-surface crust conduction in series with surface conv+rad.
 * OPEN_PARAMETER thickness; solves outer T_s by Newton (see frenchFreeSurfaceHeatFlux).
 */
inline bool applyFrenchFreeCrustToThermalBc(OphelieThermalDiffusionOneWayOptions &bc, bool enable, Real delta_m,
                                            Real k_crust, std::string &error, bool log = true)
{
    if (!enable)
    {
        return true;
    }
    if (!(delta_m >= Real(0) && k_crust > TinyReal))
    {
        error = "invalid free-crust delta/k";
        return false;
    }
    bc.enable_free_crust = true;
    bc.free_crust_delta_m = delta_m;
    bc.free_crust_k = k_crust;
    if (log)
    {
        std::cout << "[ophelie][free-crust] Level-1b free-surface crust (not resolved geometry)\n"
                  << "  delta=" << delta_m << " m  k_crust=" << k_crust
                  << " W/m/K  (outer T_s from conv+rad balance)\n";
    }
    return true;
}

inline void logFrenchZeroDHeatBudget(const FrenchZeroDHeatBudget &b, const char *tag)
{
    std::cout << "[ophelie][0d-heat] " << tag << " T_melt=" << b.t_melt << " target_P=" << b.target_joule_w
              << " W\n"
              << "  areas side/bottom/free=" << b.area_side << "/" << b.area_bottom << "/" << b.area_free
              << " m^2\n"
              << "  h_side/bottom/free=" << b.h_side_eff << "/" << b.h_bottom_eff << "/" << b.h_free << "\n"
              << "  loss side/bottom/free_conv/free_rad=" << b.side_loss_w << "/" << b.bottom_loss_w << "/"
              << b.free_conv_loss_w << "/" << b.free_rad_loss_w << " W\n"
              << "  total_loss=" << b.total_loss_w << " W  loss/P_joule=" << b.loss_over_joule << std::endl;
}

/**
 * Host regression for the skull formula and 0D budget monotonicity.
 * Returns true if Table-1 anchors and thicker skull ⇒ lower side loss.
 */
inline bool verifyFrenchSkullResistanceRegression(std::string &message)
{
    const Real h_cool = Real(300);
    const Real k = Real(1);
    const Real h_thin = frenchSkullEffectiveHeatTransfer(Real(0.002), k, h_cool);
    const Real h_thick = frenchSkullEffectiveHeatTransfer(Real(0.020), k, h_cool);
    const Real h_none = frenchSkullEffectiveHeatTransfer(Real(0), k, h_cool);
    if (!(std::abs(h_none - h_cool) / h_cool < Real(1.0e-12)))
    {
        message = "zero-thickness skull must recover h_cool";
        return false;
    }
    if (!(h_thick + TinyReal < h_thin && h_thin + TinyReal < h_cool))
    {
        message = "thicker skull must reduce h_eff monotonically";
        return false;
    }

    OphelieThermalDiffusionOneWayOptions bc;
    bc.h_side = h_cool;
    bc.h_bottom = Real(35);
    bc.h_free = Real(20);
    bc.emissivity = Real(0.8);
    bc.t_cool = Real(300);
    bc.t_ambient = Real(300);
    bc.t_rad_ambient = Real(300);

    FrenchSkullResistanceParams skull;
    skull.enable = true;
    skull.delta_side_m = Real(0.005);
    skull.delta_bottom_m = Real(0.005);
    skull.k_skull = k;
    std::string err;
    if (!applyFrenchSkullResistanceToThermalBc(bc, skull, err, false))
    {
        message = err;
        return false;
    }
    const FrenchZeroDHeatBudget with_skull =
        estimateFrenchZeroDHeatBudget(Real(0.25), Real(0.21), Real(1473), bc, Real(60000));
    bc.h_side = h_cool;
    bc.h_bottom = Real(35);
    const FrenchZeroDHeatBudget bare =
        estimateFrenchZeroDHeatBudget(Real(0.25), Real(0.21), Real(1473), bc, Real(60000));
    if (!(with_skull.side_loss_w + TinyReal < bare.side_loss_w &&
          with_skull.total_loss_w + TinyReal < bare.total_loss_w))
    {
        message = "0D budget with skull must reduce side/total loss vs bare walls";
        return false;
    }

    // Level-1b: free crust must reduce free-surface loss vs bare free surface.
    bc.h_side = with_skull.h_side_eff;
    bc.h_bottom = with_skull.h_bottom_eff;
    bc.enable_free_crust = false;
    const FrenchZeroDHeatBudget skull_no_crust =
        estimateFrenchZeroDHeatBudget(Real(0.25), Real(0.21), Real(1473), bc, Real(60000));
    bc.enable_free_crust = true;
    bc.free_crust_delta_m = Real(0.005);
    bc.free_crust_k = Real(1.0);
    const FrenchZeroDHeatBudget skull_with_crust =
        estimateFrenchZeroDHeatBudget(Real(0.25), Real(0.21), Real(1473), bc, Real(60000));
    if (!(skull_with_crust.total_loss_w + TinyReal < skull_no_crust.total_loss_w &&
          (skull_with_crust.free_conv_loss_w + skull_with_crust.free_rad_loss_w) <
              Real(0.5) * (skull_no_crust.free_conv_loss_w + skull_no_crust.free_rad_loss_w)))
    {
        message = "0D budget with free crust must cut free-surface loss by >50% (δ=5mm)";
        return false;
    }
    message = "skull_resistance_regression_ok";
    return true;
}

inline bool isFrenchSkullCliOption(const char *arg)
{
    return std::strcmp(arg, "--skull-resistance") == 0 || std::strncmp(arg, "--skull-delta=", 14) == 0 ||
           std::strncmp(arg, "--skull-delta-side=", 19) == 0 ||
           std::strncmp(arg, "--skull-delta-bottom=", 21) == 0 || std::strncmp(arg, "--skull-k=", 10) == 0 ||
           std::strncmp(arg, "--skull-h-cool-side=", 20) == 0 ||
           std::strncmp(arg, "--skull-h-cool-bottom=", 22) == 0 || std::strncmp(arg, "--skull-t-cool=", 15) == 0 ||
           std::strcmp(arg, "--free-crust-resistance") == 0 || std::strncmp(arg, "--free-crust-delta=", 19) == 0 ||
           std::strncmp(arg, "--free-crust-k=", 15) == 0;
}

inline bool tryApplyFrenchSkullCli(const char *arg, FrenchSkullResistanceParams &skull)
{
    if (std::strcmp(arg, "--skull-resistance") == 0)
    {
        skull.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-delta=", 14) == 0)
    {
        const Real d = static_cast<Real>(std::atof(arg + 14));
        skull.delta_side_m = d;
        skull.delta_bottom_m = d;
        skull.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-delta-side=", 19) == 0)
    {
        skull.delta_side_m = static_cast<Real>(std::atof(arg + 19));
        skull.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-delta-bottom=", 21) == 0)
    {
        skull.delta_bottom_m = static_cast<Real>(std::atof(arg + 21));
        skull.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-k=", 10) == 0)
    {
        skull.k_skull = static_cast<Real>(std::atof(arg + 10));
        skull.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-h-cool-side=", 20) == 0)
    {
        skull.h_cool_side = static_cast<Real>(std::atof(arg + 20));
        skull.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-h-cool-bottom=", 22) == 0)
    {
        skull.h_cool_bottom = static_cast<Real>(std::atof(arg + 22));
        skull.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--skull-t-cool=", 15) == 0)
    {
        skull.t_cool = static_cast<Real>(std::atof(arg + 15));
        skull.enable = true;
        return true;
    }
    return false;
}

struct FrenchFreeCrustParams
{
    bool enable = false;
    Real delta_m = Real(0.005);
    Real k_crust = Real(1.0);
};

inline bool tryApplyFrenchFreeCrustCli(const char *arg, FrenchFreeCrustParams &crust)
{
    if (std::strcmp(arg, "--free-crust-resistance") == 0)
    {
        crust.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--free-crust-delta=", 19) == 0)
    {
        crust.delta_m = static_cast<Real>(std::atof(arg + 19));
        crust.enable = true;
        return true;
    }
    if (std::strncmp(arg, "--free-crust-k=", 15) == 0)
    {
        crust.k_crust = static_cast<Real>(std::atof(arg + 15));
        crust.enable = true;
        return true;
    }
    return false;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_SKULL_RESISTANCE_H
