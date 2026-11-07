/**
 * @file electromagnetic_ophelie_field_export.h
 * @brief Doc183: lightweight host CSV/NPZ-ready particle field dump (no VTP).
 *
 * Syncs required variables to host, writes one CSV + sidecar JSON meta.
 * Does not modify particle state. mu_raw is evaluated read-only from material laws
 * when provided; mu_used comes from caller (particle field pointer) or
 * constant_mu_fallback when the pointer is null (constant-property F1).
 */
#ifndef ELECTROMAGNETIC_OPHELIE_FIELD_EXPORT_H
#define ELECTROMAGNETIC_OPHELIE_FIELD_EXPORT_H

#include "base_particles.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_french_material_laws.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieFieldExportResult
{
    bool ok = false;
    std::string path_csv;
    std::string path_meta;
    size_t n = 0;
    Real time_s = 0;
    Real u_max = 0;
    Real u_rms = 0;
    std::string error;
};

/**
 * Write particle field CSV at current host-synced state.
 * @param mu_used_or_null host pointer to OphelieDynamicViscosity if registered; else nullptr
 * @param constant_mu_fallback used when mu_used_or_null is nullptr (must be >0 then)
 */
inline OphelieFieldExportResult hostExportOphelieParticleFieldCsv(
    BaseParticles &particles, const std::string &out_dir, Real time_s, size_t step,
    const std::string &tag, const std::string &phase_note, bool sorted_flag,
    const OphelieFrenchGlassMaterialLaws *laws_for_mu_raw, Real mu_cap_for_note,
    bool variable_properties, const Real *mu_used_or_null, Real constant_mu_fallback,
    const std::string &coord_note = "x,y,z with +z up; gravity -z",
    const std::string &q_used_field = "OphelieThermalSourceQ")
{
    OphelieFieldExportResult r;
    r.time_s = time_s;
    std::error_code ec;
    std::filesystem::create_directories(out_dir, ec);
    if (ec)
    {
        r.error = "mkdir_failed:" + ec.message();
        return r;
    }

    if (mu_used_or_null == nullptr && !(constant_mu_fallback > Real(0)))
    {
        r.error = "missing_mu_field_and_no_constant_mu_fallback";
        return r;
    }

    // Required production fields for Doc183 coupling export.
    syncVariableToHost<UnsignedInt>(particles, "OriginalID");
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, "Velocity");
    syncVariableToHost<Real>(particles, "Temperature");
    syncVariableToHost<Real>(particles, "Density");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Real>(particles, q_used_field);

    const UnsignedInt *oid = particles.getVariableDataByName<UnsignedInt>("OriginalID");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *vel = particles.getVariableDataByName<Vecd>("Velocity");
    const Real *T = particles.getVariableDataByName<Real>("Temperature");
    const Real *rho = particles.getVariableDataByName<Real>("Density");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *Q = particles.getVariableDataByName<Real>(q_used_field);
    const size_t n = particles.TotalRealParticles();
    r.n = n;

    double u2sum = 0, umax = 0;
    for (size_t i = 0; i < n; ++i)
    {
        const double u2 = static_cast<double>(vel[i].squaredNorm());
        u2sum += u2;
        umax = std::max(umax, std::sqrt(u2));
    }
    r.u_max = static_cast<Real>(umax);
    r.u_rms = static_cast<Real>(std::sqrt(u2sum / static_cast<double>(std::max<size_t>(n, 1))));

    {
        std::vector<char> seen;
        UnsignedInt max_id = 0;
        for (size_t i = 0; i < n; ++i)
            max_id = std::max(max_id, oid[i]);
        if (max_id < n * 4 + 16)
        {
            seen.assign(static_cast<size_t>(max_id) + 1, 0);
            for (size_t i = 0; i < n; ++i)
            {
                if (seen[oid[i]])
                {
                    r.error = "non_unique_OriginalID";
                    return r;
                }
                seen[oid[i]] = 1;
            }
        }
    }

    std::ostringstream stem;
    stem << "field_" << tag << "_t" << std::fixed << std::setprecision(6) << time_s;
    r.path_csv = out_dir + "/" + stem.str() + ".csv";
    r.path_meta = out_dir + "/" + stem.str() + ".meta.json";

    {
        std::ofstream out(r.path_csv);
        if (!out)
        {
            r.error = "open_csv_failed";
            return r;
        }
        out << "OriginalID,time_s,x_m,y_m,z_m,ux_m_s,uy_m_s,uz_m_s,T_K,rho_kg_m3,vol_m3,"
               "mu_raw_Pa_s,mu_used_Pa_s,Q_used_W_m3\n";
        out.setf(std::ios::scientific);
        out.precision(10);
        for (size_t i = 0; i < n; ++i)
        {
            Real mu_u = mu_used_or_null != nullptr ? mu_used_or_null[i] : constant_mu_fallback;
            Real mu_raw = mu_u;
            if (laws_for_mu_raw != nullptr && variable_properties)
                mu_raw = std::max(laws_for_mu_raw->mu.evaluate(T[i]), TinyReal);
            const Real q = Q[i];
            if (!std::isfinite(static_cast<double>(T[i])) || !std::isfinite(static_cast<double>(vel[i][0])))
            {
                r.error = "non_finite_field";
                return r;
            }
            out << oid[i] << "," << time_s << "," << pos[i][0] << "," << pos[i][1] << "," << pos[i][2] << ","
                << vel[i][0] << "," << vel[i][1] << "," << vel[i][2] << "," << T[i] << "," << rho[i] << ","
                << vol[i] << "," << mu_raw << "," << mu_u << "," << q << "\n";
        }
    }

    {
        std::ofstream js(r.path_meta);
        js << "{\n"
           << "  \"schema\": \"doc183_ophelie_field_export_v1\",\n"
           << "  \"time_s\": " << time_s << ",\n"
           << "  \"step\": " << step << ",\n"
           << "  \"tag\": \"" << tag << "\",\n"
           << "  \"phase_note\": \"" << phase_note << "\",\n"
           << "  \"n_particles\": " << n << ",\n"
           << "  \"U_max\": " << r.u_max << ",\n"
           << "  \"U_rms\": " << r.u_rms << ",\n"
           << "  \"particle_sorted_flag\": " << (sorted_flag ? "true" : "false") << ",\n"
           << "  \"coord_note\": \"" << coord_note << "\",\n"
           << "  \"Q_field\": \"" << q_used_field << " (applied / normalized fluid source)\",\n"
           << "  \"mu_raw_note\": \"read-only material_laws.mu.evaluate(T); not written back\",\n"
           << "  \"mu_used_note\": \""
           << (mu_used_or_null != nullptr ? "OphelieDynamicViscosity after clamp"
                                          : "constant_mu_fallback (no particle viscosity field)")
           << "\",\n"
           << "  \"mu_field_present\": " << (mu_used_or_null != nullptr ? "true" : "false") << ",\n"
           << "  \"mu_cap\": " << mu_cap_for_note << ",\n"
           << "  \"csv\": \"" << r.path_csv << "\"\n"
           << "}\n";
    }
    r.ok = true;
    return r;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_FIELD_EXPORT_H
