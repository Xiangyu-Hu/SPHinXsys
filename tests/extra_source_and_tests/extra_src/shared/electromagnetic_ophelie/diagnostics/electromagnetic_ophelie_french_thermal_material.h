#ifndef ELECTROMAGNETIC_OPHELIE_FRENCH_THERMAL_MATERIAL_H
#define ELECTROMAGNETIC_OPHELIE_FRENCH_THERMAL_MATERIAL_H

#include "electromagnetic_ophelie_joule_to_heat_one_way.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Thermal material preset for one-way / diffusion prototypes. */
enum class OphelieThermalMaterialPreset
{
    /** MMS-friendly defaults (T₀=300 K, moderate k). */
    ReducedPrototype,
    /** Jacoutot et al., Chem. Eng. Process. 47 (2008) Table 1 @ 1473 K. */
    JacoutotTable1_1473K,
};

/**
 * Jacoutot 2008 Table 1 (molten glass @ 1473 K):
 *   ρ = 2750 kg/m³, Cp = 1150 J/(kg·K), λ = 4 W/(m·K), σ = 16 S/m (EM, separate).
 */
inline void applyOphelieThermalMaterialPreset(OphelieJouleHeatOneWayMaterialProps &material,
                                              OphelieThermalMaterialPreset preset)
{
    switch (preset)
    {
    case OphelieThermalMaterialPreset::JacoutotTable1_1473K:
        material.rho = 2750.0;
        material.cp = 1150.0;
        material.k = 4.0;
        material.t_initial = 1473.0;
        break;
    case OphelieThermalMaterialPreset::ReducedPrototype:
    default:
        material.rho = 2500.0;
        material.cp = 1200.0;
        material.k = 1.0;
        material.t_initial = 300.0;
        break;
    }
}

inline const char *ophelieThermalMaterialPresetName(OphelieThermalMaterialPreset preset)
{
    switch (preset)
    {
    case OphelieThermalMaterialPreset::JacoutotTable1_1473K:
        return "jacoutot_table1_1473K";
    case OphelieThermalMaterialPreset::ReducedPrototype:
    default:
        return "reduced_prototype";
    }
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FRENCH_THERMAL_MATERIAL_H
