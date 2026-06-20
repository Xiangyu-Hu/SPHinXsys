/**
 * @file rh200_material_preset.h
 * @brief RH200 flow/thermal material presets: demo vs French-literature-aligned (Jacoutot 2008 Table 1).
 */
#ifndef RH200_MATERIAL_PRESET_H
#define RH200_MATERIAL_PRESET_H

#include <cstring>
#include <iostream>
#include <string>

namespace SPH
{
namespace rh200
{

enum class Rh200MaterialPreset
{
    DemoCurrent,
    /** Jacoutot rho/cp/k/T0 + demo mu (plan Task 6-B, stability-friendly). */
    FrenchLikeStable,
    /** Jacoutot thermal + mu=4 Pa·s @ 1473 K (plan Task 6-C). */
    FrenchLikeFull,
};

struct Rh200FlowMaterialParams
{
    Real rho_glass = Real(2500);
    Real cp_glass = Real(1200);
    Real k_glass = Real(0.8);
    Real mu_glass = Real(0.05);
    Real t_initial_glass = Real(373.15);
    Real rho_rotor = Real(1000);
    Real cp_rotor = Real(1200);
    Real k_rotor_solid = Real(16);
};

inline Rh200MaterialPreset parseRh200MaterialPreset(const std::string &text)
{
    if (text == "rh200-french-like" || text == "french-like" || text == "french_like" ||
        text == "rh200-french-like-full" || text == "french-like-full" || text == "french_like_full")
    {
        return Rh200MaterialPreset::FrenchLikeFull;
    }
    if (text == "rh200-french-like-stable" || text == "french-like-stable" || text == "french_like_stable" ||
        text == "rh200-french-like-thermal" || text == "french-like-thermal")
    {
        return Rh200MaterialPreset::FrenchLikeStable;
    }
    return Rh200MaterialPreset::DemoCurrent;
}

inline const char *rh200MaterialPresetName(Rh200MaterialPreset preset)
{
    switch (preset)
    {
    case Rh200MaterialPreset::FrenchLikeStable:
        return "rh200-french-like-stable";
    case Rh200MaterialPreset::FrenchLikeFull:
        return "rh200-french-like-full";
    case Rh200MaterialPreset::DemoCurrent:
    default:
        return "rh200-demo-current";
    }
}

inline void applyJacoutot1473KThermalToRh200FlowMaterial(Rh200FlowMaterialParams &material)
{
    /** Jacoutot et al., Chem. Eng. Process. 47 (2008) Table 1 @ 1473 K. */
    material.rho_glass = Real(2750);
    material.cp_glass = Real(1150);
    material.k_glass = Real(4.0);
    material.t_initial_glass = Real(1473.0);
    material.cp_rotor = material.cp_glass;
}

inline Rh200FlowMaterialParams makeRh200FlowMaterialParams(Rh200MaterialPreset preset)
{
    Rh200FlowMaterialParams material;
    switch (preset)
    {
    case Rh200MaterialPreset::FrenchLikeStable:
        applyJacoutot1473KThermalToRh200FlowMaterial(material);
        material.mu_glass = Real(0.05);
        break;
    case Rh200MaterialPreset::FrenchLikeFull:
        applyJacoutot1473KThermalToRh200FlowMaterial(material);
        material.mu_glass = Real(4.0);
        break;
    case Rh200MaterialPreset::DemoCurrent:
    default:
        break;
    }
    return material;
}

inline void applyRh200MaterialPresetToEmDefaults(Rh200MaterialPreset preset, Real &sigma0, Real &frequency_hz,
                                                 Real &target_power)
{
    if (preset == Rh200MaterialPreset::DemoCurrent)
    {
        return;
    }
    sigma0 = Real(16);
    frequency_hz = Real(300000);
    target_power = Real(50000);
}

inline void logRh200MaterialPreset(const Rh200MaterialPreset preset, const Rh200FlowMaterialParams &material,
                                   Real sigma0, Real frequency_hz, Real target_power)
{
    std::cout << "[rh200] material preset=" << rh200MaterialPresetName(preset) << " rho_glass=" << material.rho_glass
              << " cp_glass=" << material.cp_glass << " k_glass=" << material.k_glass
              << " mu_glass=" << material.mu_glass << " T0_glass=" << material.t_initial_glass
              << " rho_rotor=" << material.rho_rotor << " cp_rotor=" << material.cp_rotor
              << " k_rotor=" << material.k_rotor_solid << " em: sigma=" << sigma0 << " f=" << frequency_hz
              << " Hz target_power=" << target_power << " W" << std::endl;
}

} // namespace rh200
} // namespace SPH

#endif // RH200_MATERIAL_PRESET_H
