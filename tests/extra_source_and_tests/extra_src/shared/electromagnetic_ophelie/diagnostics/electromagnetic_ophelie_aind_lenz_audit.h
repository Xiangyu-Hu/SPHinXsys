#ifndef ELECTROMAGNETIC_OPHELIE_AIND_LENZ_AUDIT_H
#define ELECTROMAGNETIC_OPHELIE_AIND_LENZ_AUDIT_H

#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_biot_savart.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Vol-weighted inner product treating (real, imag) as in-phase complex vector components. */
inline Real hostVolWeightedComplexVecdInnerProduct(BaseParticles &particles, const Vecd *u_real, const Vecd *u_imag,
                                                   const Vecd *v_real, const Vecd *v_imag, size_t n)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sum += vol[i] * (u_real[i].dot(v_real[i]) + u_imag[i].dot(v_imag[i]));
    }
    return sum;
}

inline Real hostVolWeightedComplexVecdNormSquared(BaseParticles &particles, const Vecd *u_real, const Vecd *u_imag,
                                                  size_t n)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sum += vol[i] * (u_real[i].squaredNorm() + u_imag[i].squaredNorm());
    }
    return sum;
}

inline Real hostVolWeightedScalarInnerProduct(BaseParticles &particles, const Real *u, const Real *v, size_t n)
{
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        sum += vol[i] * u[i] * v[i];
    }
    return sum;
}

struct OphelieAIndLenzAuditRecord
{
    Real a_coil_vol_norm = 0.0;
    Real a_ind_vol_norm = 0.0;
    Real b_coil_vol_norm = 0.0;
    Real b_ind_vol_norm = 0.0;
    Real b_total_vol_norm = 0.0;
    Real a_ind_over_a_coil = 0.0;
    Real b_ind_over_b_coil = 0.0;
    Real b_total_over_b_coil = 0.0;
    /** Vol-weighted cos(angle) between complex-vector B fields; -1=opposing, +1=aligned. */
    Real corr_b_ind_b_coil = 0.0;
    Real corr_b_total_b_coil = 0.0;
    Real corr_a_ind_a_coil = 0.0;
    Real signed_dot_b_ind_b_coil = 0.0;
    /** Dominant z component correlation (TEAM7 / uniform-Bz drive). */
    Real corr_b_ind_z_b_coil_z = 0.0;
    Real mean_b_ind_z = 0.0;
    Real mean_b_coil_z = 0.0;
    Real mean_b_ind_imag_z = 0.0;
    Real mean_b_ind_real_z = 0.0;
    /** Imag-chain phasor coupling: dot(B_ind_imag, B_coil_real). Negative => opposing quadrature. */
    Real signed_dot_b_ind_imag_b_coil_real = 0.0;
    Real corr_phasor_b_ind_imag_b_coil_real = 0.0;
    /** Analytic uniform B_drive_z reference (rotational-A box); 0 if unused. */
    Real b_drive_ref_z = 0.0;
    Real corr_b_ind_z_b_drive_z = 0.0;
    /** True when induced B is not aligned with coil B (Lenz shielding direction). */
    bool lenz_opposes_coil_b = false;
    /** Rotational benchmark: B_ind_z opposes uniform B0 z drive. */
    bool lenz_opposes_uniform_b_drive = false;
    bool lenz_audit_passed = false;
};

inline OphelieAIndLenzAuditRecord computeOphelieAIndLenzAuditFromFields(
    BaseParticles &particles, const OphelieGlassFieldNames &names, size_t n, Real b_drive_ref_z = 0.0)
{
    OphelieAIndLenzAuditRecord rec;
    rec.b_drive_ref_z = b_drive_ref_z;

    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Vecd>(particles, names.a_coil_imag);
    syncVariableToHost<Vecd>(particles, names.a_ind_real);
    syncVariableToHost<Vecd>(particles, names.a_ind_imag);
    syncVariableToHost<Vecd>(particles, names.b_coil_real);
    syncVariableToHost<Vecd>(particles, names.b_coil_imag);
    syncVariableToHost<Vecd>(particles, names.b_ind_real);
    syncVariableToHost<Vecd>(particles, names.b_ind_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Vecd *a_coil_r = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Vecd *a_coil_i = particles.getVariableDataByName<Vecd>(names.a_coil_imag);
    const Vecd *a_ind_r = particles.getVariableDataByName<Vecd>(names.a_ind_real);
    const Vecd *a_ind_i = particles.getVariableDataByName<Vecd>(names.a_ind_imag);
    const Vecd *b_coil_r = particles.getVariableDataByName<Vecd>(names.b_coil_real);
    const Vecd *b_coil_i = particles.getVariableDataByName<Vecd>(names.b_coil_imag);
    const Vecd *b_ind_r = particles.getVariableDataByName<Vecd>(names.b_ind_real);
    const Vecd *b_ind_i = particles.getVariableDataByName<Vecd>(names.b_ind_imag);

    rec.a_coil_vol_norm = hostComplexVecdPairVolWeightedNorm(particles, names.a_coil_real, names.a_coil_imag, n);
    rec.a_ind_vol_norm = hostComplexVecdPairVolWeightedNorm(particles, names.a_ind_real, names.a_ind_imag, n);
    rec.b_coil_vol_norm = hostComplexVecdPairVolWeightedNorm(particles, names.b_coil_real, names.b_coil_imag, n);
    rec.b_ind_vol_norm = hostComplexVecdPairVolWeightedNorm(particles, names.b_ind_real, names.b_ind_imag, n);

    rec.a_ind_over_a_coil = rec.a_ind_vol_norm / (rec.a_coil_vol_norm + TinyReal);
    rec.b_ind_over_b_coil = rec.b_ind_vol_norm / (rec.b_coil_vol_norm + TinyReal);

    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real b_total_sq = 0.0;
    Real b_total_dot_b_coil = 0.0;
    Real b_ind_z_sq = 0.0;
    Real b_ind_imag_z_sq = 0.0;
    Real b_coil_z_sq = 0.0;
    Real b_ind_z_dot_b_coil_z = 0.0;
    Real b_ind_imag_z_dot_drive = 0.0;
    Real vol_sum = 0.0;
    Real b_ind_z_sum = 0.0;
    Real b_coil_z_sum = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd b_total_r = b_coil_r[i] + b_ind_r[i];
        const Vecd b_total_i = b_coil_i[i] + b_ind_i[i];
        b_total_sq += vol[i] * (b_total_r.squaredNorm() + b_total_i.squaredNorm());
        b_total_dot_b_coil +=
            vol[i] * (b_total_r.dot(b_coil_r[i]) + b_total_i.dot(b_coil_i[i]));

        const Real b_ind_z_r = b_ind_r[i][2];
        const Real b_ind_z_i = b_ind_i[i][2];
        const Real b_coil_z = b_coil_r[i][2];
        b_ind_z_sq += vol[i] * (b_ind_z_r * b_ind_z_r + b_ind_z_i * b_ind_z_i);
        b_ind_imag_z_sq += vol[i] * b_ind_z_i * b_ind_z_i;
        b_coil_z_sq += vol[i] * b_coil_z * b_coil_z;
        b_ind_z_dot_b_coil_z += vol[i] * (b_ind_z_r * b_coil_z);
        if (std::abs(b_drive_ref_z) > TinyReal)
        {
            b_ind_imag_z_dot_drive += vol[i] * b_ind_z_i * b_drive_ref_z;
        }
        vol_sum += vol[i];
        b_ind_z_sum += vol[i] * b_ind_z_r;
        b_coil_z_sum += vol[i] * b_coil_z;
        rec.mean_b_ind_imag_z += vol[i] * b_ind_z_i;
        rec.mean_b_ind_real_z += vol[i] * b_ind_z_r;
    }
    rec.b_total_vol_norm = std::sqrt(b_total_sq);
    rec.b_total_over_b_coil = rec.b_total_vol_norm / (rec.b_coil_vol_norm + TinyReal);
    if (vol_sum > TinyReal)
    {
        rec.mean_b_ind_z = b_ind_z_sum / vol_sum;
        rec.mean_b_coil_z = b_coil_z_sum / vol_sum;
        rec.mean_b_ind_imag_z /= vol_sum;
        rec.mean_b_ind_real_z /= vol_sum;
    }

    rec.signed_dot_b_ind_b_coil =
        hostVolWeightedComplexVecdInnerProduct(particles, b_ind_r, b_ind_i, b_coil_r, b_coil_i, n);
    rec.signed_dot_b_ind_imag_b_coil_real = 0.0;
    Real b_ind_imag_only_sq = 0.0;
    Real b_coil_real_only_sq = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        rec.signed_dot_b_ind_imag_b_coil_real += vol[i] * b_ind_i[i].dot(b_coil_r[i]);
        b_ind_imag_only_sq += vol[i] * b_ind_i[i].squaredNorm();
        b_coil_real_only_sq += vol[i] * b_coil_r[i].squaredNorm();
    }
    rec.corr_phasor_b_ind_imag_b_coil_real = rec.signed_dot_b_ind_imag_b_coil_real /
                                             (std::sqrt(b_ind_imag_only_sq * b_coil_real_only_sq) + TinyReal);
    const Real b_ind_sq = hostVolWeightedComplexVecdNormSquared(particles, b_ind_r, b_ind_i, n);
    const Real b_coil_sq = hostVolWeightedComplexVecdNormSquared(particles, b_coil_r, b_coil_i, n);
    rec.corr_b_ind_b_coil = rec.signed_dot_b_ind_b_coil / (std::sqrt(b_ind_sq * b_coil_sq) + TinyReal);
    rec.corr_b_total_b_coil =
        b_total_dot_b_coil / (std::sqrt(b_total_sq * b_coil_sq) + TinyReal);

    const Real a_dot =
        hostVolWeightedComplexVecdInnerProduct(particles, a_ind_r, a_ind_i, a_coil_r, a_coil_i, n);
    const Real a_ind_sq = hostVolWeightedComplexVecdNormSquared(particles, a_ind_r, a_ind_i, n);
    const Real a_coil_sq = hostVolWeightedComplexVecdNormSquared(particles, a_coil_r, a_coil_i, n);
    rec.corr_a_ind_a_coil = a_dot / (std::sqrt(a_ind_sq * a_coil_sq) + TinyReal);

    rec.corr_b_ind_z_b_coil_z = b_ind_z_dot_b_coil_z / (std::sqrt(b_ind_z_sq * b_coil_z_sq) + TinyReal);
    if (std::abs(b_drive_ref_z) > TinyReal)
    {
        const Real drive_sq = b_drive_ref_z * b_drive_ref_z * vol_sum;
        rec.corr_b_ind_z_b_drive_z = b_ind_imag_z_dot_drive / (std::sqrt(b_ind_imag_z_sq * drive_sq) + TinyReal);
        rec.lenz_opposes_uniform_b_drive = rec.corr_b_ind_z_b_drive_z < Real(0);
    }

    if (rec.b_coil_vol_norm > TinyReal)
    {
        if (std::abs(rec.corr_b_ind_b_coil) > Real(0.15))
        {
            rec.lenz_opposes_coil_b = rec.corr_b_ind_b_coil < Real(0);
        }
        else
        {
            rec.lenz_opposes_coil_b = rec.corr_phasor_b_ind_imag_b_coil_real < Real(0);
        }
        rec.lenz_audit_passed = rec.lenz_opposes_coil_b;
    }
    else if (std::abs(b_drive_ref_z) > TinyReal)
    {
        rec.lenz_opposes_coil_b = rec.corr_b_ind_z_b_drive_z < Real(0);
        rec.lenz_audit_passed = rec.lenz_opposes_uniform_b_drive;
    }
    else
    {
        rec.lenz_audit_passed = rec.lenz_opposes_coil_b;
    }
    return rec;
}

template <class ExecutionPolicy>
inline OphelieAIndLenzAuditRecord execOphelieGlassSelfInducedBiotAndLenzAudit(
    SolidBody &glass_body, const OphelieGlassFieldNames &names, OphelieParameters &params,
    const std::string &j_real_field, const std::string &j_imag_field, Real b_drive_ref_z = 0.0)
{
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    hostZeroVecdField(particles, names.b_ind_real, n);
    hostZeroVecdField(particles, names.b_ind_imag, n);

    StateDynamics<ExecutionPolicy, ComputeOphelieGlassSelfInducedBiotSavartCK> compute_aind(
        glass_body, names, params, j_real_field, j_imag_field);
    compute_aind.exec();

    return computeOphelieAIndLenzAuditFromFields(particles, names, n, b_drive_ref_z);
}

inline void printOphelieAIndLenzAudit(const OphelieAIndLenzAuditRecord &rec, const std::string &context_tag = "")
{
    std::cout << "[ophelie] P4 A_ind/Lenz audit";
    if (!context_tag.empty())
    {
        std::cout << " (" << context_tag << ")";
    }
    std::cout << ": Aind/Acoil=" << rec.a_ind_over_a_coil << " Bind/Bcoil=" << rec.b_ind_over_b_coil
              << " Btotal/Bcoil=" << rec.b_total_over_b_coil << " corr(Bind,Bcoil)=" << rec.corr_b_ind_b_coil
              << " corr(Btotal,Bcoil)=" << rec.corr_b_total_b_coil               << " corr(Aind,Acoil)=" << rec.corr_a_ind_a_coil
              << " corr_phasor(Bind_i,Bcoil_r)=" << rec.corr_phasor_b_ind_imag_b_coil_real
              << " corr(Bind_z,Bcoil_z)=" << rec.corr_b_ind_z_b_coil_z << " mean_Bind_z_r=" << rec.mean_b_ind_real_z
              << " mean_Bind_z_i=" << rec.mean_b_ind_imag_z << " mean_Bcoil_z=" << rec.mean_b_coil_z;
    if (std::abs(rec.b_drive_ref_z) > TinyReal)
    {
        std::cout << " Bdrive_z=" << rec.b_drive_ref_z << " corr(Bind_z,Bdrive_z)=" << rec.corr_b_ind_z_b_drive_z
                  << " lenz_vs_drive=" << (rec.lenz_opposes_uniform_b_drive ? 1 : 0);
    }
    std::cout << " lenz_opposes_coil=" << (rec.lenz_opposes_coil_b ? 1 : 0)
              << " lenz_pass=" << (rec.lenz_audit_passed ? 1 : 0) << std::endl;
}

inline void writeOphelieAIndLenzAuditCsv(const std::string &path, const OphelieAIndLenzAuditRecord &rec,
                                         const std::string &context_tag = "")
{
    const fs::path parent = fs::path(path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(path);
    if (!out)
    {
        std::cerr << "[ophelie] could not write Lenz audit CSV: " << path << std::endl;
        return;
    }
    out << "context,a_ind_over_a_coil,b_ind_over_b_coil,b_total_over_b_coil,corr_b_ind_b_coil,corr_b_total_b_coil,"
           "corr_a_ind_a_coil,corr_phasor_b_ind_imag_b_coil_real,signed_dot_b_ind_b_coil,signed_dot_b_ind_imag_b_coil_real,"
           "corr_b_ind_z_b_coil_z,mean_b_ind_z,mean_b_coil_z,mean_b_ind_imag_z,b_drive_ref_z,"
           "corr_b_ind_z_b_drive_z,lenz_opposes_coil_b,lenz_opposes_uniform_b_drive,lenz_audit_passed\n";
    out << context_tag << "," << rec.a_ind_over_a_coil << "," << rec.b_ind_over_b_coil << "," << rec.b_total_over_b_coil
        << "," << rec.corr_b_ind_b_coil << "," << rec.corr_b_total_b_coil << "," << rec.corr_a_ind_a_coil << ","
        << rec.corr_phasor_b_ind_imag_b_coil_real << "," << rec.signed_dot_b_ind_b_coil << ","
        << rec.signed_dot_b_ind_imag_b_coil_real << "," << rec.corr_b_ind_z_b_coil_z << "," << rec.mean_b_ind_z << ","
        << rec.mean_b_coil_z << "," << rec.mean_b_ind_imag_z << "," << rec.b_drive_ref_z << ","
        << rec.corr_b_ind_z_b_drive_z << "," << (rec.lenz_opposes_coil_b ? 1 : 0) << ","
        << (rec.lenz_opposes_uniform_b_drive ? 1 : 0) << "," << (rec.lenz_audit_passed ? 1 : 0) << "\n";
    std::cout << "[ophelie] P4 A_ind/Lenz audit CSV: " << path << std::endl;
}

inline OphelieAIndLenzAuditRecord runTeam7AIndLenzAuditAfterOneWay(
    SolidBody &plate_body, const OphelieGlassFieldNames &names, const std::string &j_imag_field,
    const std::string &csv_path, const std::string &context_tag)
{
    BaseParticles &particles = plate_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    OphelieAIndLenzAuditRecord rec = computeOphelieAIndLenzAuditFromFields(particles, names, n, 0.0);
    printOphelieAIndLenzAudit(rec, context_tag);
    if (!csv_path.empty())
    {
        writeOphelieAIndLenzAuditCsv(csv_path, rec, context_tag);
    }
    return rec;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif
