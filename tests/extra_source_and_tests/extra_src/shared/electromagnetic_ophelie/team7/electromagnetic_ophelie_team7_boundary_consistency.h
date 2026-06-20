#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_BOUNDARY_CONSISTENCY_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_BOUNDARY_CONSISTENCY_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux_operator_audit.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_team7_boundary_normal.h"
#include "electromagnetic_ophelie_team7_edge_recon_boundary.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct Team7BoundaryConsistencyRecord
{
    std::string partition_name;
    std::string test_name;
    Real metric_value = 0.0;
    Real threshold = 0.0;
    bool passed = false;
    std::string notes;
};

inline bool team7BoundaryConsistencyPassLowerIsBetter(Real metric, Real threshold)
{
    return metric <= threshold;
}

inline bool team7BoundaryConsistencyPassHigherIsBetter(Real metric, Real threshold)
{
    return metric >= threshold;
}

inline StdVec<Team7BoundaryConsistencyRecord> computeTeam7BoundaryConsistencyAudit(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const BoundingBoxd &plate_bbox, Real z_lower_m, Real z_upper_m, Real partition_skin_h_m, Real boundary_width_m,
    const StdVec<Team7EdgeReconBoundaryPartitionCompare> *boundary_compare = nullptr)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.e_imag);
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *e_edge = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    const Vecd *j_edge = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Vecd *e_imag = particles.getVariableDataByName<Vecd>(names.e_imag);
    const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
    const Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const bool have_sphinx_normals = team7ParticlesHaveSphinxNormals(particles);
    const Real omega = params.omega();

    const Team7PlateHoleGrid hole_grid = buildTeam7PlateHoleGridFromParticles(particles, plate_bbox, partition_skin_h_m);
    const StdVec<Team7EdgeFluxPartitionKind> kinds = team7EdgeFluxPartitionKindsForAudit();

    struct Accum
    {
        size_t count = 0;
        size_t boundary_count = 0;
        Real vol_sum = 0.0;
        Real grad_phi_sq = 0.0;
        Real omega_a_sq = 0.0;
        Real edge_mis_sq = 0.0;
        Real edge_norm_sq = 0.0;
        Real e_imag_mis_sq = 0.0;
        Real e_imag_sq = 0.0;
        Real j_sigma_mis_sq = 0.0;
        Real j_sq = 0.0;
        Real j_normal_sq = 0.0;
        Real j_tangent_sq = 0.0;
        Real normal_angle_sum_deg = 0.0;
        Real normal_angle_max_deg = 0.0;
    };
    StdVec<Accum> acc(kinds.size());
    for (size_t k = 0; k < kinds.size(); ++k)
    {
        acc[k] = Accum{};
    }

    const Real *signed_distance =
        have_sphinx_normals ? particles.getVariableDataByName<Real>(kTeam7SignedDistance) : nullptr;

    for (size_t i = 0; i < n; ++i)
    {
        const Team7ParticlePartitionInfo pinfo = classifyTeam7ParticlePartition(
            pos[i], plate_bbox, partition_skin_h_m, z_lower_m, z_upper_m, hole_grid, boundary_width_m,
            signed_distance != nullptr ? signed_distance[i] : Real(0), signed_distance != nullptr);
        StdVec<Team7EdgeFluxPartitionKind> targets;
        team7ParticlePartitionTargets(pinfo, targets);
        const Real v = vol[i];
        const Vecd omega_a = omega * a_coil_real[i];
        const Vecd e_em = -grad_phi[i] - omega_a;
        const Vecd j_sigma_e = sigma[i] * e_edge[i];

        Vecd sphinx_normal = Vecd::Zero();
        const bool on_boundary_shell =
            have_sphinx_normals && getTeam7BoundaryNormal(particles, i, boundary_width_m, sphinx_normal);

        for (Team7EdgeFluxPartitionKind target : targets)
        {
            size_t idx = 0;
            for (; idx < kinds.size(); ++idx)
            {
                if (kinds[idx] == target)
                {
                    break;
                }
            }
            Accum &a = acc[idx];
            ++a.count;
            a.vol_sum += v;
            if (on_boundary_shell)
            {
                ++a.boundary_count;
                const Vecd n_analytic =
                    team7PartitionNormal(pinfo.surface_kind, pos[i], plate_bbox, hole_grid);
                const Real cos_angle = std::clamp(sphinx_normal.dot(n_analytic), Real(-1), Real(1));
                const Real angle_deg = std::acos(std::abs(cos_angle)) * 180.0 / Pi;
                a.normal_angle_sum_deg += angle_deg;
                a.normal_angle_max_deg = std::max(a.normal_angle_max_deg, angle_deg);
            }
            a.grad_phi_sq += v * grad_phi[i].squaredNorm();
            a.omega_a_sq += v * omega_a.squaredNorm();
            const Vecd edge_delta = e_edge[i] - e_em;
            a.edge_mis_sq += v * edge_delta.squaredNorm();
            a.edge_norm_sq += v * e_edge[i].squaredNorm();
            const Vecd e_imag_delta = e_imag[i] - e_em;
            a.e_imag_mis_sq += v * e_imag_delta.squaredNorm();
            a.e_imag_sq += v * e_imag[i].squaredNorm();
            a.j_sigma_mis_sq += v * (j_edge[i] - j_sigma_e).squaredNorm();
            a.j_sq += v * j_edge[i].squaredNorm();
            Vecd n_for_jn = on_boundary_shell ? sphinx_normal
                                                : team7PartitionNormal(pinfo.surface_kind, pos[i], plate_bbox, hole_grid);
            const Real jn = j_edge[i].dot(n_for_jn);
            const Real jt_sq = (j_edge[i] - jn * n_for_jn).squaredNorm();
            a.j_normal_sq += v * jn * jn;
            a.j_tangent_sq += v * jt_sq;
        }
    }

    StdVec<Team7BoundaryConsistencyRecord> records;
    const auto add_record = [&](const std::string &partition, const std::string &test_name, Real metric, Real threshold,
                                bool pass, const std::string &notes)
    {
        records.push_back({partition, test_name, metric, threshold, pass, notes});
    };

    for (size_t k = 0; k < kinds.size(); ++k)
    {
        const std::string partition = team7EdgeFluxPartitionName(kinds[k]);
        const Accum &a = acc[k];
        if (a.count == 0)
        {
            continue;
        }
        const Real grad_phi_vol = std::sqrt(a.grad_phi_sq);
        const Real omega_a_vol = std::sqrt(a.omega_a_sq);
        const Real omega_over_grad = grad_phi_vol > TinyReal ? omega_a_vol / grad_phi_vol : 0.0;
        const Real e_edge_em_mis = std::sqrt(a.edge_mis_sq) / (std::sqrt(a.edge_norm_sq) + TinyReal);
        const Real e_imag_emf_mis = std::sqrt(a.e_imag_mis_sq) / (std::sqrt(a.e_imag_sq) + TinyReal);
        const Real j_vol = std::sqrt(a.j_sq);
        const Real j_sigma_mis = j_vol > TinyReal ? std::sqrt(a.j_sigma_mis_sq) / j_vol : 0.0;
        const Real jn_vol = std::sqrt(a.j_normal_sq);
        const Real jt_vol = std::sqrt(a.j_tangent_sq);
        const Real jn_over_jt = jn_vol / (jt_vol + TinyReal);
        const Real normal_mean_deg =
            a.boundary_count > 0 ? a.normal_angle_sum_deg / static_cast<Real>(a.boundary_count) : 0.0;

        add_record(partition, "emf_omega_a_over_grad_phi", omega_over_grad, 5.0,
                   team7BoundaryConsistencyPassLowerIsBetter(omega_over_grad, 5.0),
                   "EMF balance: |omegaA|/|gradPhi| (TEAM7 one-way ~2)");
        add_record(partition, "e_edge_em_mismatch_raw", e_edge_em_mis, 0.10,
                   team7BoundaryConsistencyPassLowerIsBetter(e_edge_em_mis, 0.10),
                   "||E_edge-(-gradPhi-omegaA)||/||E_edge||");
        add_record(partition, "e_imag_emf_mismatch", e_imag_emf_mis, 0.10,
                   team7BoundaryConsistencyPassLowerIsBetter(e_imag_emf_mis, 0.10),
                   "||E_imag-(-gradPhi-omegaA)||/||E_imag|| particle-gradient path");
        add_record(partition, "j_sigma_e_edge_closure", j_sigma_mis, 0.01,
                   team7BoundaryConsistencyPassLowerIsBetter(j_sigma_mis, 0.01), "||J_edge-sigma*E_edge||/||J_edge||");
        const bool skip_no_flux_partition =
            kinds[k] == Team7EdgeFluxPartitionKind::TrueInterior ||
            kinds[k] == Team7EdgeFluxPartitionKind::All;
        if (!skip_no_flux_partition)
        {
            add_record(partition, "no_flux_jn_over_jt_raw", jn_over_jt, 0.10,
                       team7BoundaryConsistencyPassLowerIsBetter(jn_over_jt, 0.10),
                       "no-flux diagnostic: |Jn|/|Jt| using SPHinXsys n on shell");
        }
        if (a.boundary_count > 0 && kinds[k] != Team7EdgeFluxPartitionKind::TrueInterior)
        {
            add_record(partition, "sphinx_normal_mean_angle_deg", normal_mean_deg, 20.0,
                       team7BoundaryConsistencyPassLowerIsBetter(normal_mean_deg, 20.0),
                       "mean angle(SPHinXsys n, partition analytic n) on boundary shell");
            add_record(partition, "sphinx_normal_max_angle_deg", a.normal_angle_max_deg, 45.0,
                       team7BoundaryConsistencyPassLowerIsBetter(a.normal_angle_max_deg, 45.0),
                       "max angle on boundary shell");
        }
    }

    if (boundary_compare != nullptr)
    {
        for (const Team7EdgeReconBoundaryPartitionCompare &compare : *boundary_compare)
        {
            if (compare.has_project_normal)
            {
                add_record(compare.partition_name, "e_edge_em_mismatch_project_normal",
                           compare.project_normal.e_edge_em_mismatch, 0.10,
                           team7BoundaryConsistencyPassLowerIsBetter(compare.project_normal.e_edge_em_mismatch, 0.10),
                           "P5.1 post-hoc E projection diagnostic");
                add_record(compare.partition_name, "no_flux_jn_over_jt_project_normal",
                           compare.project_normal.jn_over_jt, 0.01,
                           team7BoundaryConsistencyPassLowerIsBetter(compare.project_normal.jn_over_jt, 0.01),
                           "expect ~0 after E<-E-n(n·E)");
            }
            if (compare.has_tangent_ls)
            {
                add_record(compare.partition_name, "e_edge_em_mismatch_tangent_ls",
                           compare.tangent_ls.e_edge_em_mismatch, 0.10,
                           team7BoundaryConsistencyPassLowerIsBetter(compare.tangent_ls.e_edge_em_mismatch, 0.10),
                           "P5.2 tangent-plane LS diagnostic");
                add_record(compare.partition_name, "e_edge_em_mismatch_tangent_pair",
                           compare.tangent_ls.e_edge_em_mismatch_tangent_pair, 0.10,
                           team7BoundaryConsistencyPassLowerIsBetter(compare.tangent_ls.e_edge_em_mismatch_tangent_pair,
                                                                     0.10),
                           "fair compare: E_edge vs projected E_em on tangent plane");
                add_record(compare.partition_name, "em_tangential_mismatch_tangent_ls",
                           compare.tangent_ls.em_tangential_mismatch, 0.10,
                           team7BoundaryConsistencyPassLowerIsBetter(compare.tangent_ls.em_tangential_mismatch, 0.10),
                           "tangential EMF component mismatch");
                if (compare.partition_name != "true_interior" && compare.partition_name != "all")
                {
                    add_record(compare.partition_name, "no_flux_jn_over_jt_tangent_ls", compare.tangent_ls.jn_over_jt,
                               0.01, team7BoundaryConsistencyPassLowerIsBetter(compare.tangent_ls.jn_over_jt, 0.01),
                               "expect ~0 after tangent LS");
                }
            }
        }
    }

    add_record("reference_box", "constant_A_gauge_cancellation_J_ratio", Real(3.6e-5), Real(1.0e-3), true,
               "P1a high_sigma box: J_phi/J_edge~0 expected (HIGH_SIGMA_EDGE_FLUX_P1A_LOG.md)");
    add_record("reference_box", "rotational_A_uniform_B_J_ratio", Real(0.91), Real(0.80), true,
               "P1a high_sigma box: J retained ~91% (not TEAM7 geometry)");
    return records;
}

inline void printTeam7BoundaryConsistencyReport(const StdVec<Team7BoundaryConsistencyRecord> &records)
{
    size_t n_pass = 0;
    size_t n_fail = 0;
    std::cout << "[team7] P5.3 boundary operator consistency audit:" << std::endl;
    for (const Team7BoundaryConsistencyRecord &record : records)
    {
        if (record.partition_name == "reference_box")
        {
            continue;
        }
        if (record.passed)
        {
            ++n_pass;
        }
        else
        {
            ++n_fail;
        }
        std::cout << "  " << (record.passed ? "PASS" : "FAIL") << " " << record.partition_name << " "
                  << record.test_name << " metric=" << record.metric_value << " threshold=" << record.threshold;
        if (!record.notes.empty())
        {
            std::cout << " (" << record.notes << ")";
        }
        std::cout << std::endl;
    }
    std::cout << "[team7] P5.3 TEAM7 plate summary: pass=" << n_pass << " fail=" << n_fail << std::endl;
    for (const Team7BoundaryConsistencyRecord &record : records)
    {
        if (record.partition_name != "reference_box")
        {
            continue;
        }
        std::cout << "  REF " << record.test_name << " metric=" << record.metric_value
                  << " (documented box benchmark, not re-run here)" << std::endl;
    }
}

inline void writeTeam7BoundaryConsistencyCsv(const std::string &output_path,
                                             const StdVec<Team7BoundaryConsistencyRecord> &records)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    out << "partition,test_name,metric,threshold,passed,notes\n";
    for (const Team7BoundaryConsistencyRecord &record : records)
    {
        out << record.partition_name << "," << record.test_name << "," << record.metric_value << ","
            << record.threshold << "," << (record.passed ? 1 : 0) << ",\"" << record.notes << "\"\n";
    }
    std::cout << "[team7] P5.3 boundary consistency CSV: " << output_path << std::endl;
}

template <class ExecutionPolicy>
inline void runTeam7BoundaryConsistencyAudit(SolidBody &plate_body, Inner<> &plate_inner,
                                           const OphelieGlassFieldNames &plate_names, OphelieParameters &params,
                                           const BoundingBoxd &plate_bbox, Real z_lower_m, Real z_upper_m, Real dp,
                                           const std::string &csv_path,
                                           const StdVec<Team7EdgeReconBoundaryPartitionCompare> *boundary_compare =
                                               nullptr)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(plate_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(plate_inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(plate_inner,
                                                                                                        plate_names);
    compute_grad_phi.exec();

    const Real boundary_width_m = params.edge_recon_boundary_width_factor_ * dp;
    const StdVec<Team7BoundaryConsistencyRecord> records = computeTeam7BoundaryConsistencyAudit(
        plate_body.getBaseParticles(), plate_names, params, plate_bbox, z_lower_m, z_upper_m, dp, boundary_width_m,
        boundary_compare);
    printTeam7BoundaryConsistencyReport(records);
    writeTeam7BoundaryConsistencyCsv(csv_path, records);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_BOUNDARY_CONSISTENCY_H
