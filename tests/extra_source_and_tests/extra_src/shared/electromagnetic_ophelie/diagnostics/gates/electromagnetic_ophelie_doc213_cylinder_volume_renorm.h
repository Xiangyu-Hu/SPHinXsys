/**
 * @file electromagnetic_ophelie_doc213_cylinder_volume_renorm.h
 * @brief Doc213: match discrete ΣV and ΣV r² to continuum cylinder moments.
 *
 * Lattice packing overstates ΣV and ΣV r² → continuous power floor (G0b/G6).
 * Apply once after particle generation, before operators / power integrals:
 *   V_i' = α V_i (1 + β r_i²)
 * choosing α,β so ΣV' = π R² H and ΣV' r² = π R⁴ H / 2.
 * Mass_i scaled by V_i'/V_i so ρ = Mass/V stays consistent.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_DOC213_CYLINDER_VOLUME_RENORM_H
#define ELECTROMAGNETIC_OPHELIE_DOC213_CYLINDER_VOLUME_RENORM_H

#include "electromagnetic_ophelie_device_sync.h"

#include <cmath>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieDoc213CylinderVolumeRenormReport
{
    size_t n = 0;
    Real sum_v_before = 0.0;
    Real sum_v_r2_before = 0.0;
    Real sum_v_r4_before = 0.0;
    Real vol_exact = 0.0;
    Real mom_r2_exact = 0.0;
    Real alpha = 1.0;
    Real beta = 0.0;
    Real sum_v_after = 0.0;
    Real sum_v_r2_after = 0.0;
    Real rel_vol_before = 0.0;
    Real rel_mom_before = 0.0;
    Real rel_vol_after = 0.0;
    Real rel_mom_after = 0.0;
    bool applied = false;
    std::string mode = "off";
};

/**
 * @param mode "total" → uniform α only (match ΣV); "moment" → α,β match ΣV and ΣVr².
 */
inline OphelieDoc213CylinderVolumeRenormReport hostOphelieDoc213CylinderVolumeRenorm(
    BaseParticles &particles, const Vecd &center, Real radius, Real height, const std::string &mode)
{
    OphelieDoc213CylinderVolumeRenormReport rep;
    rep.mode = mode;
    if (mode != "total" && mode != "moment")
    {
        return rep;
    }

    const size_t n = particles.TotalRealParticles();
    rep.n = n;
    rep.vol_exact = Pi * radius * radius * height;
    rep.mom_r2_exact = Pi * radius * radius * radius * radius * height / Real(2);

    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, "Mass");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *mass = particles.getVariableDataByName<Real>("Mass");

    double s0 = 0, s2 = 0, s4 = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r2 = dx * dx + dy * dy;
        const double v = static_cast<double>(vol[i]);
        s0 += v;
        s2 += v * static_cast<double>(r2);
        s4 += v * static_cast<double>(r2) * static_cast<double>(r2);
    }
    rep.sum_v_before = static_cast<Real>(s0);
    rep.sum_v_r2_before = static_cast<Real>(s2);
    rep.sum_v_r4_before = static_cast<Real>(s4);
    rep.rel_vol_before = static_cast<Real>((s0 - rep.vol_exact) / (rep.vol_exact + TinyReal));
    rep.rel_mom_before = static_cast<Real>((s2 - rep.mom_r2_exact) / (rep.mom_r2_exact + TinyReal));

    Real alpha = Real(1);
    Real beta = Real(0);
    if (mode == "total")
    {
        alpha = rep.vol_exact / (rep.sum_v_before + TinyReal);
    }
    else
    {
        // β = (M S0 − V S2) / (V S4 − M S2); α = V / (S0 + β S2)
        const double Ve = static_cast<double>(rep.vol_exact);
        const double Me = static_cast<double>(rep.mom_r2_exact);
        const double num = Me * s0 - Ve * s2;
        const double den = Ve * s4 - Me * s2;
        if (std::abs(den) < 1.0e-30)
        {
            std::cerr << "[ophelie][doc213] volume_renorm moment: singular den, fallback to total"
                      << std::endl;
            alpha = rep.vol_exact / (rep.sum_v_before + TinyReal);
            beta = Real(0);
            rep.mode = "total_fallback";
        }
        else
        {
            beta = static_cast<Real>(num / den);
            alpha = static_cast<Real>(Ve / (s0 + static_cast<double>(beta) * s2 + TinyReal));
        }
        if (!(alpha > Real(0)))
        {
            std::cerr << "[ophelie][doc213] volume_renorm moment: alpha<=0, fallback to total"
                      << std::endl;
            alpha = rep.vol_exact / (rep.sum_v_before + TinyReal);
            beta = Real(0);
            rep.mode = "total_fallback";
        }
    }
    rep.alpha = alpha;
    rep.beta = beta;

    double s0a = 0, s2a = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const Real dx = pos[i][0] - center[0];
        const Real dy = pos[i][1] - center[1];
        const Real r2 = dx * dx + dy * dy;
        const Real factor = alpha * (Real(1) + beta * r2);
        const Real v_old = vol[i];
        const Real v_new = v_old * factor;
        vol[i] = v_new;
        if (mass != nullptr && std::abs(v_old) > TinyReal)
        {
            mass[i] *= (v_new / v_old);
        }
        s0a += static_cast<double>(v_new);
        s2a += static_cast<double>(v_new) * static_cast<double>(r2);
    }
    rep.sum_v_after = static_cast<Real>(s0a);
    rep.sum_v_r2_after = static_cast<Real>(s2a);
    rep.rel_vol_after = static_cast<Real>((s0a - rep.vol_exact) / (rep.vol_exact + TinyReal));
    rep.rel_mom_after = static_cast<Real>((s2a - rep.mom_r2_exact) / (rep.mom_r2_exact + TinyReal));
    rep.applied = true;

    syncVariableToDevice<Real>(particles, "VolumetricMeasure");
    if (mass != nullptr)
    {
        syncVariableToDevice<Real>(particles, "Mass");
    }

    std::cout << "[ophelie][doc213] volume_renorm mode=" << rep.mode << " n=" << n
              << " alpha=" << alpha << " beta=" << beta
              << " rel_vol " << rep.rel_vol_before << "->" << rep.rel_vol_after
              << " rel_mom_r2 " << rep.rel_mom_before << "->" << rep.rel_mom_after << std::endl;
    return rep;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif
