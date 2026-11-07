/**
 * @file electromagnetic_ophelie_viscosity_operator_diag.h
 * @brief Doc177 P0: fixed-geometry explicit ViscousForce vs implicit DampingPairwise compare.
 */
#ifndef ELECTROMAGNETIC_OPHELIE_VISCOSITY_OPERATOR_DIAG_H
#define ELECTROMAGNETIC_OPHELIE_VISCOSITY_OPERATOR_DIAG_H

#include "base_particles.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_french_material_laws.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct ViscosityOperatorDiagRow
{
    size_t original_id = 0;
    size_t index = 0;
    size_t n_classic_nb = 0;
    Real mu = 0;
    Real mass = 0;
    Real vol = 0;
    Real ax_exp = 0, ay_exp = 0, az_exp = 0;
    Real ax_imp = 0, ay_imp = 0, az_imp = 0;
    Real a_exp_norm = 0, a_imp_norm = 0;
    Real rel_err = 0;
};

/**
 * Compare a_exp = ViscousForce/mass to a_imp ≈ (v_after - v_before)/dt from a user-provided
 * implicit exec callback. Caller must:
 *  - freeze geometry (no position update)
 *  - zero ForcePrior viscous history before explicit force
 *  - provide the same dt used inside the implicit InteractionSplit exec
 */
template <typename ExplicitExec, typename ImplicitExec>
inline bool hostCompareViscosityOperators(
    BaseParticles &particles, ExplicitExec &&exec_explicit, ImplicitExec &&exec_implicit, Real dt,
    const std::string &csv_path, size_t max_rows, Real &mean_rel_err, Real &max_rel_err)
{
    if (!(dt > TinyReal))
        return false;

    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Vecd>(particles, "ViscousForce");
    syncVariableToHost<Real>(particles, "Mass");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, kOphelieDynamicViscosityField);
    syncVariableToHost<UnsignedInt>(particles, "OriginalID");

    Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real *mu = particles.getVariableDataByName<Real>(kOphelieDynamicViscosityField);
    UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
    const size_t n = particles.TotalRealParticles();

    std::vector<Vecd> v0(n);
    for (size_t i = 0; i < n; ++i)
        v0[i] = vel[i];

    // Explicit force → acceleration
    exec_explicit();
    syncVariableToHost<Vecd>(particles, "ViscousForce");
    Vecd *vf = particles.getVariableDataByName<Vecd>("ViscousForce");
    std::vector<Vecd> a_exp(n);
    for (size_t i = 0; i < n; ++i)
        a_exp[i] = vf[i] / std::max(mass[i], TinyReal);

    // Restore velocity, run implicit over dt
    for (size_t i = 0; i < n; ++i)
        vel[i] = v0[i];
    syncVariableToDevice<Vecd>(particles, "Velocity");
    exec_implicit(dt);
    syncVariableToHost<Vecd>(particles, "Velocity");

    mean_rel_err = 0;
    max_rel_err = 0;
    size_t counted = 0;
    std::vector<ViscosityOperatorDiagRow> rows;
    rows.reserve(std::min(max_rows, n));

    const size_t stride = std::max<size_t>(1, n / std::max<size_t>(1, max_rows));
    for (size_t i = 0; i < n; i += stride)
    {
        const Vecd a_imp = (vel[i] - v0[i]) / dt;
        const Real ne = a_exp[i].norm();
        const Real ni = a_imp.norm();
        const Real denom = std::max(ne, TinyReal);
        const Real rel = (a_imp - a_exp[i]).norm() / denom;
        mean_rel_err += rel;
        max_rel_err = std::max(max_rel_err, rel);
        ++counted;

        ViscosityOperatorDiagRow r;
        r.original_id = static_cast<size_t>(oid[i]);
        r.index = i;
        r.mu = mu[i];
        r.mass = mass[i];
        r.vol = vol[i];
        r.ax_exp = a_exp[i][0];
        r.ay_exp = a_exp[i][1];
        r.az_exp = a_exp[i][2];
        r.ax_imp = a_imp[0];
        r.ay_imp = a_imp[1];
        r.az_imp = a_imp[2];
        r.a_exp_norm = ne;
        r.a_imp_norm = ni;
        r.rel_err = rel;
        if (rows.size() < max_rows)
            rows.push_back(r);
    }
    if (counted > 0)
        mean_rel_err /= static_cast<Real>(counted);

    // restore velocity
    for (size_t i = 0; i < n; ++i)
        vel[i] = v0[i];
    syncVariableToDevice<Vecd>(particles, "Velocity");

    std::ofstream out(csv_path);
    if (!out)
        return false;
    out << "original_id,index,mu,mass,vol,ax_exp,ay_exp,az_exp,ax_imp,ay_imp,az_imp,"
           "a_exp_norm,a_imp_norm,rel_err\n";
    for (const auto &r : rows)
    {
        out << r.original_id << "," << r.index << "," << r.mu << "," << r.mass << "," << r.vol << ","
            << r.ax_exp << "," << r.ay_exp << "," << r.az_exp << "," << r.ax_imp << "," << r.ay_imp << ","
            << r.az_imp << "," << r.a_exp_norm << "," << r.a_imp_norm << "," << r.rel_err << "\n";
    }
    std::cout << "[doc177][viscosity-diag] n_sample=" << counted << " mean_rel_err=" << mean_rel_err
              << " max_rel_err=" << max_rel_err << " dt=" << dt << " csv=" << csv_path << std::endl;
    return true;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif
