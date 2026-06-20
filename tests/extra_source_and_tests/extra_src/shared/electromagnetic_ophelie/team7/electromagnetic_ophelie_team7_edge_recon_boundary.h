#ifndef ELECTROMAGNETIC_OPHELIE_TEAM7_EDGE_RECON_BOUNDARY_H
#define ELECTROMAGNETIC_OPHELIE_TEAM7_EDGE_RECON_BOUNDARY_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_edge_flux_operator_audit.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_team7_boundary_normal.h"
#include "interaction_ck.h"
#include "update_body_relation.h"
#include "vector_functions.h"

#include <Eigen/Dense>
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

inline Vecd projectOphelieEfieldNoNormalFlux(const Vecd &e_field, const Vecd &unit_normal)
{
    return e_field - unit_normal * unit_normal.dot(e_field);
}

inline void buildOphelieSurfaceTangentBasis(const Vecd &unit_normal, Vecd &t1_out, Vecd &t2_out)
{
    const Vecd reference_axis =
        std::abs(unit_normal[2]) < Real(0.9) ? Vecd(0.0, 0.0, 1.0) : Vecd(1.0, 0.0, 0.0);
    t1_out = unit_normal.cross(reference_axis);
    const Real t1_norm = t1_out.norm();
    t1_out = t1_norm > TinyReal ? t1_out / t1_norm : Vecd(1.0, 0.0, 0.0);
    t2_out = unit_normal.cross(t1_out);
    const Real t2_norm = t2_out.norm();
    t2_out = t2_norm > TinyReal ? t2_out / t2_norm : Vecd(0.0, 1.0, 0.0);
}

struct Team7EdgeReconBoundaryVariantMetrics
{
    size_t n_particles = 0;
    size_t n_boundary_shell = 0;
    Real vol_sum_m3 = 0.0;
    Real j_vol_norm = 0.0;
    Real e_edge_em_mismatch = 0.0;
    /** E_edge vs E_em with both projected onto tangent plane (fair no-flux compare). */
    Real e_edge_em_mismatch_tangent_pair = 0.0;
    /** E_edge projected; E_em raw (legacy unfair compare). */
    Real e_edge_em_mismatch_projected_edge = 0.0;
    /** Both E_edge and E_em projected with n·E=0. */
    Real e_edge_em_mismatch_projected_pair = 0.0;
    Real em_normal_mismatch = 0.0;
    Real em_tangential_mismatch = 0.0;
    Real j_normal_l2 = 0.0;
    Real j_tangential_l2 = 0.0;
    Real jn_over_jt = 0.0;
    Real e_normal_l2 = 0.0;
    Real p_recon_w = 0.0;
    Real tangent_ls_fallback_fraction = 0.0;
};

inline Real team7VolWeightedRelativeMismatch(Real mis_sq, Real ref_norm_sq)
{
    return std::sqrt(mis_sq) / (std::sqrt(ref_norm_sq) + TinyReal);
}

inline void team7AccumulateEmfComponentMismatches(Real vol, const Vecd &e_edge, const Vecd &e_em, const Vecd &unit_normal,
                                                  Real &normal_mis_sq, Real &tangential_mis_sq)
{
    const Real en_edge = e_edge.dot(unit_normal);
    const Real en_em = e_em.dot(unit_normal);
    const Vecd et_edge = e_edge - en_edge * unit_normal;
    const Vecd et_em = e_em - en_em * unit_normal;
    const Real normal_delta = en_edge - en_em;
    const Vecd tangential_delta = et_edge - et_em;
    normal_mis_sq += vol * normal_delta * normal_delta;
    tangential_mis_sq += vol * tangential_delta.squaredNorm();
}

struct Team7EdgeReconBoundaryPartitionCompare
{
    std::string partition_name;
    Team7EdgeReconBoundaryVariantMetrics raw;
    Team7EdgeReconBoundaryVariantMetrics project_normal;
    Team7EdgeReconBoundaryVariantMetrics tangent_ls;
    bool has_project_normal = false;
    bool has_tangent_ls = false;
};

inline Team7EdgeReconBoundaryVariantMetrics finalizeTeam7EdgeReconBoundaryVariantMetrics(
    size_t count, size_t boundary_count, Real vol_sum, Real j_sq, Real edge_mis_sq, Real edge_norm_sq,
    Real j_normal_sq, Real j_tangent_sq, Real e_normal_sq, Real p_recon_sum, Real fallback_vol = 0.0,
    Real tangent_pair_mis_sq = 0.0, Real projected_edge_mis_sq = 0.0, Real projected_pair_mis_sq = 0.0,
    Real normal_mis_sq = 0.0, Real tangential_mis_sq = 0.0)
{
    Team7EdgeReconBoundaryVariantMetrics metrics;
    metrics.n_particles = count;
    metrics.n_boundary_shell = boundary_count;
    metrics.vol_sum_m3 = vol_sum;
    metrics.j_vol_norm = std::sqrt(j_sq);
    metrics.e_edge_em_mismatch = team7VolWeightedRelativeMismatch(edge_mis_sq, edge_norm_sq);
    metrics.e_edge_em_mismatch_tangent_pair = team7VolWeightedRelativeMismatch(tangent_pair_mis_sq, edge_norm_sq);
    metrics.e_edge_em_mismatch_projected_edge = team7VolWeightedRelativeMismatch(projected_edge_mis_sq, edge_norm_sq);
    metrics.e_edge_em_mismatch_projected_pair = team7VolWeightedRelativeMismatch(projected_pair_mis_sq, edge_norm_sq);
    metrics.em_normal_mismatch = team7VolWeightedRelativeMismatch(normal_mis_sq, edge_norm_sq);
    metrics.em_tangential_mismatch = team7VolWeightedRelativeMismatch(tangential_mis_sq, edge_norm_sq);
    metrics.j_normal_l2 = std::sqrt(j_normal_sq);
    metrics.j_tangential_l2 = std::sqrt(j_tangent_sq);
    metrics.jn_over_jt = metrics.j_normal_l2 / (metrics.j_tangential_l2 + TinyReal);
    metrics.e_normal_l2 = std::sqrt(e_normal_sq);
    metrics.p_recon_w = p_recon_sum;
    metrics.tangent_ls_fallback_fraction = vol_sum > TinyReal ? fallback_vol / vol_sum : 0.0;
    return metrics;
}

inline void ensureTeam7TangentLsDiagnosticFields(BaseParticles &particles)
{
    syncVariableToDevice<Vecd>(particles, kTeam7NormalDirection);
    syncVariableToDevice<Real>(particles, kTeam7SignedDistance);
    syncVariableToDevice<Vecd>(particles, kTeam7EEdgeTangentLsDiag);
    syncVariableToDevice<Vecd>(particles, kTeam7JEdgeTangentLsDiag);
    syncVariableToDevice<Real>(particles, kTeam7EdgeTangentLsFallback);
}

template <typename... RelationTypes>
class ReconstructOphelieEdgeFluxTangentLsBoundaryDiagnosticCK;

template <template <typename...> class RelationType, typename... Parameters>
class ReconstructOphelieEdgeFluxTangentLsBoundaryDiagnosticCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ReconstructOphelieEdgeFluxTangentLsBoundaryDiagnosticCK(RelationType<Parameters...> &relation,
                                                            const OphelieGlassFieldNames &names,
                                                            const OphelieEdgeFluxComponent &component, Real omega,
                                                            Real pair_weight_regularization,
                                                            Real recon_condition_threshold, Real boundary_width_m,
                                                            Real solver_local_rhs_scale = Real(1),
                                                            bool use_distance_tangent_norm = false)
        : BaseInteraction(relation), omega_(omega), pair_weight_regularization_(pair_weight_regularization),
          a_sign_(component.a_sign), recon_condition_threshold_(recon_condition_threshold),
          boundary_width_m_(boundary_width_m), solver_local_rhs_scale_(solver_local_rhs_scale),
          use_distance_tangent_norm_(use_distance_tangent_norm),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_phi_(this->particles_->template getVariableByName<Real>(component.phi_field)),
          dv_active_a_(this->particles_->template getVariableByName<Vecd>(component.active_a_field)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_signed_distance_(this->particles_->template getVariableByName<Real>(kTeam7SignedDistance)),
          dv_normal_(this->particles_->template getVariableByName<Vecd>(kTeam7NormalDirection)),
          dv_e_edge_recon_(this->particles_->template getVariableByName<Vecd>(component.e_recon_field)),
          dv_j_edge_recon_(this->particles_->template getVariableByName<Vecd>(component.j_recon_field)),
          dv_e_tangent_(this->particles_->template getVariableByName<Vecd>(kTeam7EEdgeTangentLsDiag)),
          dv_j_tangent_(this->particles_->template getVariableByName<Vecd>(kTeam7JEdgeTangentLsDiag)),
          dv_fallback_(this->particles_->template getVariableByName<Real>(kTeam7EdgeTangentLsFallback))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              phi_(encloser.dv_phi_->DelegatedData(ex_policy)),
              active_a_(encloser.dv_active_a_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              signed_distance_(encloser.dv_signed_distance_->DelegatedData(ex_policy)),
              normal_direction_(encloser.dv_normal_->DelegatedData(ex_policy)),
              e_edge_recon_(encloser.dv_e_edge_recon_->DelegatedData(ex_policy)),
              j_edge_recon_(encloser.dv_j_edge_recon_->DelegatedData(ex_policy)),
              e_tangent_(encloser.dv_e_tangent_->DelegatedData(ex_policy)),
              j_tangent_(encloser.dv_j_tangent_->DelegatedData(ex_policy)),
              fallback_(encloser.dv_fallback_->DelegatedData(ex_policy)), omega_(encloser.omega_),
              a_sign_(encloser.a_sign_), pair_weight_regularization_(encloser.pair_weight_regularization_),
              recon_condition_threshold_(encloser.recon_condition_threshold_),
              boundary_width_m_(encloser.boundary_width_m_), solver_local_rhs_scale_(encloser.solver_local_rhs_scale_),
              use_distance_tangent_norm_(encloser.use_distance_tangent_norm_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_[index_i];
            const Real sigma_i = sigma_[index_i];
            const Vecd a_i = active_a_[index_i];
            const Real signed_distance = signed_distance_[index_i];
            const bool on_boundary_shell = std::abs(signed_distance) <= boundary_width_m_;
            if (!on_boundary_shell)
            {
                e_tangent_[index_i] = e_edge_recon_[index_i];
                j_tangent_[index_i] = j_edge_recon_[index_i];
                fallback_[index_i] = 0.0;
                return;
            }

            Vecd unit_normal = normal_direction_[index_i];
            const Real normal_norm = unit_normal.norm();
            if (normal_norm <= TinyReal)
            {
                e_tangent_[index_i] = e_edge_recon_[index_i];
                j_tangent_[index_i] = j_edge_recon_[index_i];
                fallback_[index_i] = 1.0;
                return;
            }
            unit_normal /= normal_norm;

            Vecd t1 = Vecd::Zero();
            Vecd t2 = Vecd::Zero();
            buildOphelieSurfaceTangentBasis(unit_normal, t1, t2);

            Eigen::Matrix2d m_acc = Eigen::Matrix2d::Zero();
            Eigen::Vector2d b_acc = Eigen::Vector2d::Zero();
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Vecd xj_minus_xi = -r_ij_vec;
                const Real distance = r_ij_vec.norm();
                if (distance <= TinyReal)
                {
                    continue;
                }
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real c_ij = computeOphelieEdgeFluxPairConductance(
                    sigma_i, sigma_[index_j], dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                    reference_smoothing_length_);
                const Real edge_drop = computeOphelieEdgeFluxEdgeDrop(phi_i, phi_[index_j], a_i, active_a_[index_j],
                                                                      xj_minus_xi, omega_, a_sign_);
                const Vecd r_tangent = r_ij_vec - unit_normal * unit_normal.dot(r_ij_vec);
                const Real distance_t = r_tangent.norm();
                if (distance_t <= TinyReal)
                {
                    continue;
                }
                const Vecd e_hat_t = r_tangent / distance_t;
                const Eigen::Vector2d u(static_cast<double>(e_hat_t.dot(t1)), static_cast<double>(e_hat_t.dot(t2)));
                const double directional_e =
                    static_cast<double>(edge_drop) /
                    static_cast<double>(use_distance_tangent_norm_ ? distance_t : distance) *
                    static_cast<double>(solver_local_rhs_scale_);
                const double w_ij = static_cast<double>(std::abs(c_ij));
                m_acc += w_ij * u * u.transpose();
                b_acc += w_ij * directional_e * u;
            }

            const double condition_proxy = m_acc.norm();
            if (condition_proxy > static_cast<double>(recon_condition_threshold_))
            {
                const double tikhonov_eps = static_cast<double>(SqrtEps) * condition_proxy * condition_proxy;
                const Eigen::Vector2d e2 = inverseTikhonov(m_acc, tikhonov_eps) * b_acc;
                Vecd e_i = static_cast<Real>(e2[0]) * t1 + static_cast<Real>(e2[1]) * t2;
                if (solver_local_rhs_scale_ < Real(1) - TinyReal)
                {
                    e_i /= solver_local_rhs_scale_;
                }
                if (!e_i.allFinite())
                {
                    e_tangent_[index_i] = e_edge_recon_[index_i];
                    j_tangent_[index_i] = j_edge_recon_[index_i];
                    fallback_[index_i] = 1.0;
                }
                else
                {
                    e_tangent_[index_i] = e_i;
                    j_tangent_[index_i] = sigma_i * e_i;
                    fallback_[index_i] = 0.0;
                }
            }
            else
            {
                e_tangent_[index_i] = e_edge_recon_[index_i];
                j_tangent_[index_i] = j_edge_recon_[index_i];
                fallback_[index_i] = 1.0;
            }
        }

      protected:
        Real *Vol_;
        Real *phi_;
        Vecd *active_a_;
        Real *sigma_;
        Real *signed_distance_;
        Vecd *normal_direction_;
        Vecd *e_edge_recon_;
        Vecd *j_edge_recon_;
        Vecd *e_tangent_;
        Vecd *j_tangent_;
        Real *fallback_;
        Real omega_;
        Real a_sign_;
        Real pair_weight_regularization_;
        Real recon_condition_threshold_;
        Real boundary_width_m_;
        Real solver_local_rhs_scale_;
        bool use_distance_tangent_norm_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    Real pair_weight_regularization_;
    Real a_sign_;
    Real recon_condition_threshold_;
    Real boundary_width_m_;
    Real solver_local_rhs_scale_;
    bool use_distance_tangent_norm_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Vecd> *dv_active_a_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_signed_distance_;
    DiscreteVariable<Vecd> *dv_normal_;
    DiscreteVariable<Vecd> *dv_e_edge_recon_;
    DiscreteVariable<Vecd> *dv_j_edge_recon_;
    DiscreteVariable<Vecd> *dv_e_tangent_;
    DiscreteVariable<Vecd> *dv_j_tangent_;
    DiscreteVariable<Real> *dv_fallback_;
};

template <typename... Parameters>
class ReconstructOphelieEdgeFluxTangentLsBoundaryDiagnosticCK<Inner<Parameters...>>
    : public ReconstructOphelieEdgeFluxTangentLsBoundaryDiagnosticCK<Base, Inner<Parameters...>>
{
  public:
    ReconstructOphelieEdgeFluxTangentLsBoundaryDiagnosticCK(Inner<Parameters...> &inner_relation,
                                                            const OphelieGlassFieldNames &names,
                                                            const OphelieEdgeFluxComponent &component, Real omega,
                                                            Real pair_weight_regularization,
                                                            Real recon_condition_threshold, Real boundary_width_m,
                                                            Real solver_local_rhs_scale = Real(1),
                                                            bool use_distance_tangent_norm = false)
        : ReconstructOphelieEdgeFluxTangentLsBoundaryDiagnosticCK<Base, Inner<Parameters...>>(
              inner_relation, names, component, omega, pair_weight_regularization, recon_condition_threshold,
              boundary_width_m, solver_local_rhs_scale, use_distance_tangent_norm)
    {
    }
};

inline StdVec<Team7EdgeReconBoundaryPartitionCompare> computeTeam7EdgeReconBoundaryCompareAudit(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const BoundingBoxd &plate_bbox, Real z_lower_m, Real z_upper_m, Real partition_skin_h_m, Real boundary_width_m,
    bool include_project_normal, bool include_tangent_ls)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    if (include_tangent_ls)
    {
        syncVariableToHost<Vecd>(particles, kTeam7EEdgeTangentLsDiag);
        syncVariableToHost<Vecd>(particles, kTeam7JEdgeTangentLsDiag);
        syncVariableToHost<Real>(particles, kTeam7EdgeTangentLsFallback);
    }

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *e_edge = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    const Vecd *j_edge = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
    const Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *e_tangent = include_tangent_ls ? particles.getVariableDataByName<Vecd>(kTeam7EEdgeTangentLsDiag) : nullptr;
    const Vecd *j_tangent = include_tangent_ls ? particles.getVariableDataByName<Vecd>(kTeam7JEdgeTangentLsDiag) : nullptr;
    const Real *tangent_fallback =
        include_tangent_ls ? particles.getVariableDataByName<Real>(kTeam7EdgeTangentLsFallback) : nullptr;
    const bool have_sphinx_normals = team7ParticlesHaveSphinxNormals(particles);
    const Real omega = params.omega();

    const Team7PlateHoleGrid hole_grid = buildTeam7PlateHoleGridFromParticles(particles, plate_bbox, partition_skin_h_m);
    const StdVec<Team7EdgeFluxPartitionKind> kinds = team7EdgeFluxPartitionKindsForAudit();
    StdVec<Team7EdgeReconBoundaryPartitionCompare> records(kinds.size());
    for (size_t k = 0; k < kinds.size(); ++k)
    {
        records[k].partition_name = team7EdgeFluxPartitionName(kinds[k]);
        records[k].has_project_normal = include_project_normal;
        records[k].has_tangent_ls = include_tangent_ls;
    }

    struct Accum
    {
        size_t count = 0;
        size_t boundary_count = 0;
        Real vol_sum = 0.0;
        Real j_sq = 0.0;
        Real edge_mis_sq = 0.0;
        Real edge_norm_sq = 0.0;
        Real tangent_pair_mis_sq = 0.0;
        Real projected_edge_mis_sq = 0.0;
        Real projected_pair_mis_sq = 0.0;
        Real normal_mis_sq = 0.0;
        Real tangential_mis_sq = 0.0;
        Real j_normal_sq = 0.0;
        Real j_tangent_sq = 0.0;
        Real e_normal_sq = 0.0;
        Real p_recon_sum = 0.0;
        Real fallback_vol = 0.0;
    };
    StdVec<Accum> raw_acc(records.size());
    StdVec<Accum> project_acc(records.size());
    StdVec<Accum> tangent_acc(records.size());

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

        Vecd sphinx_normal = Vecd::Zero();
        const bool on_boundary_shell =
            have_sphinx_normals && getTeam7BoundaryNormal(particles, i, boundary_width_m, sphinx_normal);

        const Vecd e_raw = e_edge[i];
        const Vecd j_raw = j_edge[i];
        Vecd e_project = e_raw;
        Vecd j_project = j_raw;
        if (include_project_normal && on_boundary_shell)
        {
            e_project = projectOphelieEfieldNoNormalFlux(e_raw, sphinx_normal);
            j_project = sigma[i] * e_project;
        }
        const Vecd e_tangent_i = include_tangent_ls && e_tangent != nullptr ? e_tangent[i] : e_raw;
        const Vecd j_tangent_i = include_tangent_ls && j_tangent != nullptr ? j_tangent[i] : j_raw;

        Vecd n_for_jn = Vecd::Zero();
        if (on_boundary_shell)
        {
            n_for_jn = sphinx_normal;
        }
        else
        {
            n_for_jn = team7PartitionNormal(pinfo.surface_kind, pos[i], plate_bbox, hole_grid);
        }

        const Vecd e_em_project =
            on_boundary_shell ? projectOphelieEfieldNoNormalFlux(e_em, n_for_jn) : e_em;
        const Vecd e_raw_tangent =
            on_boundary_shell ? projectOphelieEfieldNoNormalFlux(e_raw, n_for_jn) : e_raw;
        const Vecd e_em_tangent =
            on_boundary_shell ? projectOphelieEfieldNoNormalFlux(e_em, n_for_jn) : e_em;

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

            const auto accumulate = [&](Accum &a, const Vecd &e_var, const Vecd &j_var, bool count_boundary,
                                          bool count_fallback, const Vecd &e_em_ref, const Vecd &e_em_proj_ref,
                                          const Vecd &e_em_tangent_ref)
            {
                ++a.count;
                a.vol_sum += v;
                if (count_boundary)
                {
                    ++a.boundary_count;
                }
                if (count_fallback)
                {
                    a.fallback_vol += v;
                }
                a.j_sq += v * j_var.squaredNorm();
                const Vecd e_var_projected =
                    on_boundary_shell ? projectOphelieEfieldNoNormalFlux(e_var, n_for_jn) : e_var;
                a.edge_mis_sq += v * (e_var - e_em_ref).squaredNorm();
                a.edge_norm_sq += v * e_var.squaredNorm();
                a.tangent_pair_mis_sq += v * (e_var - e_em_tangent_ref).squaredNorm();
                a.projected_edge_mis_sq += v * (e_var_projected - e_em_ref).squaredNorm();
                a.projected_pair_mis_sq += v * (e_var_projected - e_em_proj_ref).squaredNorm();
                team7AccumulateEmfComponentMismatches(v, e_var, e_em_ref, n_for_jn, a.normal_mis_sq, a.tangential_mis_sq);
                const Real jn = j_var.dot(n_for_jn);
                const Real jt_sq = (j_var - jn * n_for_jn).squaredNorm();
                a.j_normal_sq += v * jn * jn;
                a.j_tangent_sq += v * jt_sq;
                const Real en = e_var.dot(n_for_jn);
                a.e_normal_sq += v * en * en;
                a.p_recon_sum += Real(0.5) * j_var.dot(e_var) * v;
            };

            accumulate(raw_acc[idx], e_raw, j_raw, on_boundary_shell, false, e_em, e_em_project, e_em_tangent);
            if (include_project_normal)
            {
                accumulate(project_acc[idx], e_project, j_project, on_boundary_shell, false, e_em, e_em_project,
                           e_em_tangent);
            }
            if (include_tangent_ls)
            {
                const bool used_tangent_fallback =
                    tangent_fallback != nullptr && tangent_fallback[i] > Real(0.5);
                accumulate(tangent_acc[idx], e_tangent_i, j_tangent_i, on_boundary_shell, used_tangent_fallback, e_em,
                           e_em_project, e_em_tangent);
            }
        }
    }

    for (size_t k = 0; k < records.size(); ++k)
    {
        const Accum &ra = raw_acc[k];
        records[k].raw = finalizeTeam7EdgeReconBoundaryVariantMetrics(
            ra.count, ra.boundary_count, ra.vol_sum, ra.j_sq, ra.edge_mis_sq, ra.edge_norm_sq, ra.j_normal_sq,
            ra.j_tangent_sq, ra.e_normal_sq, ra.p_recon_sum, 0.0, ra.tangent_pair_mis_sq, ra.projected_edge_mis_sq,
            ra.projected_pair_mis_sq, ra.normal_mis_sq, ra.tangential_mis_sq);
        if (include_project_normal)
        {
            const Accum &pa = project_acc[k];
            records[k].project_normal = finalizeTeam7EdgeReconBoundaryVariantMetrics(
                pa.count, pa.boundary_count, pa.vol_sum, pa.j_sq, pa.edge_mis_sq, pa.edge_norm_sq, pa.j_normal_sq,
                pa.j_tangent_sq, pa.e_normal_sq, pa.p_recon_sum, 0.0, pa.tangent_pair_mis_sq, pa.projected_edge_mis_sq,
                pa.projected_pair_mis_sq, pa.normal_mis_sq, pa.tangential_mis_sq);
        }
        if (include_tangent_ls)
        {
            const Accum &ta = tangent_acc[k];
            records[k].tangent_ls = finalizeTeam7EdgeReconBoundaryVariantMetrics(
                ta.count, ta.boundary_count, ta.vol_sum, ta.j_sq, ta.edge_mis_sq, ta.edge_norm_sq, ta.j_normal_sq,
                ta.j_tangent_sq, ta.e_normal_sq, ta.p_recon_sum, ta.fallback_vol, ta.tangent_pair_mis_sq,
                ta.projected_edge_mis_sq, ta.projected_pair_mis_sq, ta.normal_mis_sq, ta.tangential_mis_sq);
        }
    }
    return records;
}

inline void printTeam7EdgeReconBoundaryCompareReport(const StdVec<Team7EdgeReconBoundaryPartitionCompare> &records,
                                                     OphelieEdgeReconBoundaryMode mode, Real boundary_width_m)
{
    std::cout << "[team7] P5 edge-recon boundary diagnostic mode=" << ophelieEdgeReconBoundaryModeName(mode)
              << " boundary_width_m=" << boundary_width_m
              << " (diagnostic-only; production E/J unchanged)" << std::endl;
    for (const Team7EdgeReconBoundaryPartitionCompare &record : records)
    {
        std::cout << "  " << record.partition_name << " raw: e_edge_em_mis=" << record.raw.e_edge_em_mismatch
                  << " em_tan_mis=" << record.raw.em_tangential_mismatch
                  << " em_proj_pair=" << record.raw.e_edge_em_mismatch_projected_pair
                  << " Jn/Jt=" << record.raw.jn_over_jt << " |Jn|=" << record.raw.j_normal_l2
                  << " P_recon=" << record.raw.p_recon_w;
        if (record.has_project_normal)
        {
            std::cout << " | project-normal: e_edge_em_mis=" << record.project_normal.e_edge_em_mismatch
                      << " em_tan_mis=" << record.project_normal.em_tangential_mismatch
                      << " Jn/Jt=" << record.project_normal.jn_over_jt
                      << " |Jn|=" << record.project_normal.j_normal_l2
                      << " P_recon=" << record.project_normal.p_recon_w;
        }
        if (record.has_tangent_ls)
        {
            std::cout << " | tangent-ls: e_edge_em_mis=" << record.tangent_ls.e_edge_em_mismatch
                      << " em_tan_mis=" << record.tangent_ls.em_tangential_mismatch
                      << " em_proj_pair=" << record.tangent_ls.e_edge_em_mismatch_projected_pair
                      << " Jn/Jt=" << record.tangent_ls.jn_over_jt << " |Jn|=" << record.tangent_ls.j_normal_l2
                      << " P_recon=" << record.tangent_ls.p_recon_w
                      << " fallback_frac=" << record.tangent_ls.tangent_ls_fallback_fraction;
        }
        std::cout << " (n_boundary=" << record.raw.n_boundary_shell << ")" << std::endl;
    }
}

inline void writeTeam7EdgeReconBoundaryCompareCsv(const std::string &output_path,
                                                  const StdVec<Team7EdgeReconBoundaryPartitionCompare> &records,
                                                  OphelieEdgeReconBoundaryMode mode, Real boundary_width_m,
                                                  Real partition_skin_h_m)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    out << "partition,boundary_mode,boundary_width_m,partition_skin_h_m,variant,n_particles,n_boundary_shell,"
           "vol_sum_m3,j_vol_norm,e_edge_em_mismatch,e_edge_em_mismatch_tangent_pair,e_edge_em_mismatch_projected_pair,"
           "em_normal_mismatch,em_tangential_mismatch,j_normal_l2,j_tangential_l2,Jn_over_Jt,e_normal_l2,p_recon_w,"
           "tangent_ls_fallback_fraction\n";
    for (const Team7EdgeReconBoundaryPartitionCompare &record : records)
    {
        const auto write_variant = [&](const std::string &variant, const Team7EdgeReconBoundaryVariantMetrics &m)
        {
            out << record.partition_name << "," << ophelieEdgeReconBoundaryModeName(mode) << "," << boundary_width_m
                << "," << partition_skin_h_m << "," << variant << "," << m.n_particles << "," << m.n_boundary_shell
                << "," << m.vol_sum_m3 << "," << m.j_vol_norm << "," << m.e_edge_em_mismatch << ","
                << m.e_edge_em_mismatch_tangent_pair << "," << m.e_edge_em_mismatch_projected_pair << ","
                << m.em_normal_mismatch << "," << m.em_tangential_mismatch << "," << m.j_normal_l2 << ","
                << m.j_tangential_l2 << "," << m.jn_over_jt << "," << m.e_normal_l2 << "," << m.p_recon_w << ","
                << m.tangent_ls_fallback_fraction << "\n";
        };
        write_variant("raw", record.raw);
        if (record.has_project_normal)
        {
            write_variant("project-normal", record.project_normal);
        }
        if (record.has_tangent_ls)
        {
            write_variant("tangent-ls", record.tangent_ls);
        }
    }
    std::cout << "[team7] P5 edge-recon boundary compare CSV: " << output_path << std::endl;
}

inline void printTeam7EdgeReconBoundaryHoleSummary(const StdVec<Team7EdgeReconBoundaryPartitionCompare> &records)
{
    const Team7EdgeReconBoundaryPartitionCompare *hole = nullptr;
    const Team7EdgeReconBoundaryPartitionCompare *boundary_shell = nullptr;
    const Team7EdgeReconBoundaryPartitionCompare *true_interior = nullptr;
    for (const Team7EdgeReconBoundaryPartitionCompare &record : records)
    {
        if (record.partition_name == "hole_lateral")
        {
            hole = &record;
        }
        if (record.partition_name == "boundary_shell_all")
        {
            boundary_shell = &record;
        }
        if (record.partition_name == "true_interior")
        {
            true_interior = &record;
        }
    }
    if (hole == nullptr)
    {
        return;
    }
    std::cout << "[team7] P5 hole_lateral boundary delta:"
              << " e_edge_em raw=" << hole->raw.e_edge_em_mismatch
              << " em_tan=" << hole->raw.em_tangential_mismatch;
    if (hole->has_project_normal)
    {
        std::cout << " project-normal=" << hole->project_normal.e_edge_em_mismatch
                  << " em_tan=" << hole->project_normal.em_tangential_mismatch;
    }
    if (hole->has_tangent_ls)
    {
        std::cout << " tangent-ls=" << hole->tangent_ls.e_edge_em_mismatch
                  << " em_tan=" << hole->tangent_ls.em_tangential_mismatch
                  << " fallback_frac=" << hole->tangent_ls.tangent_ls_fallback_fraction;
    }
    std::cout << " |Jn| raw=" << hole->raw.j_normal_l2;
    if (hole->has_project_normal)
    {
        std::cout << " project-normal=" << hole->project_normal.j_normal_l2;
    }
    if (hole->has_tangent_ls)
    {
        std::cout << " tangent-ls=" << hole->tangent_ls.j_normal_l2;
    }
    if (boundary_shell != nullptr)
    {
        std::cout << " (boundary_shell Jn/Jt raw=" << boundary_shell->raw.jn_over_jt;
        if (boundary_shell->has_tangent_ls)
        {
            std::cout << " tangent-ls=" << boundary_shell->tangent_ls.jn_over_jt;
        }
        std::cout << ")";
    }
    if (true_interior != nullptr)
    {
        std::cout << " (true_interior e_edge_em raw=" << true_interior->raw.e_edge_em_mismatch;
        if (true_interior->has_tangent_ls)
        {
            std::cout << " tangent-ls=" << true_interior->tangent_ls.e_edge_em_mismatch;
        }
        std::cout << ")";
    }
    std::cout << std::endl;
}

template <class ExecutionPolicy>
inline StdVec<Team7EdgeReconBoundaryPartitionCompare> runTeam7EdgeReconBoundaryDiagnostic(
    SolidBody &plate_body, Inner<> &plate_inner, const OphelieGlassFieldNames &plate_names, OphelieParameters &params,
    const BoundingBoxd &plate_bbox, Real z_lower_m, Real z_upper_m, Real dp, const std::string &csv_path)
{
    if (params.edge_recon_boundary_mode_ == OphelieEdgeReconBoundaryMode::None)
    {
        return {};
    }
    if (ophelieEdgeReconBoundaryModeUsesProductionClosure(params.edge_recon_boundary_mode_))
    {
        std::cout << "[team7] edge-recon boundary diagnostic skipped: mode="
                  << ophelieEdgeReconBoundaryModeName(params.edge_recon_boundary_mode_)
                  << " is production-only (use --team7-boundary-consistency-audit=1)" << std::endl;
        return {};
    }
    if (!team7ParticlesHaveSphinxNormals(plate_body.getBaseParticles()))
    {
        std::cerr << "[team7] edge-recon boundary diagnostic skipped: NormalDirection missing (regenerate reload with "
                     "--relax=1)"
                  << std::endl;
        return {};
    }

    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(plate_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(plate_inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(plate_inner,
                                                                                                        plate_names);
    compute_grad_phi.exec();

    const Real boundary_width_m = params.edge_recon_boundary_width_factor_ * dp;
    const bool include_project_normal =
        params.edge_recon_boundary_mode_ == OphelieEdgeReconBoundaryMode::ProjectNormal ||
        params.edge_recon_boundary_mode_ == OphelieEdgeReconBoundaryMode::TangentLs;
    const bool include_tangent_ls = params.edge_recon_boundary_mode_ == OphelieEdgeReconBoundaryMode::TangentLs;

    if (include_tangent_ls)
    {
        ensureTeam7TangentLsDiagnosticFields(plate_body.getBaseParticles());
        const OphelieEdgeFluxComponent imag_component = makeOphelieEdgeFluxImagComponent(plate_names, params);
        const Real solver_local_rhs_scale = ophelieEdgeFluxEffectiveSolverLocalRhsScale(params);
        const bool use_distance_tangent_norm =
            params.tangent_ls_distance_norm_ == OphelieTangentLsDistanceNorm::Tangent;
        InteractionDynamicsCK<ExecutionPolicy, ReconstructOphelieEdgeFluxTangentLsBoundaryDiagnosticCK<Inner<>>>
            reconstruct_tangent_ls(plate_inner, plate_names, imag_component, params.omega(),
                                   params.pair_weight_regularization_, params.edge_recon_condition_threshold_,
                                   boundary_width_m, solver_local_rhs_scale, use_distance_tangent_norm);
        reconstruct_tangent_ls.exec();
        std::cout << "[team7] P5.2 tangent-ls boundary reconstruction finished (diagnostic fields only)"
                  << " distance_norm=" << ophelieTangentLsDistanceNormName(params.tangent_ls_distance_norm_)
                  << std::endl;
    }

    const StdVec<Team7EdgeReconBoundaryPartitionCompare> records = computeTeam7EdgeReconBoundaryCompareAudit(
        plate_body.getBaseParticles(), plate_names, params, plate_bbox, z_lower_m, z_upper_m, dp, boundary_width_m,
        include_project_normal, include_tangent_ls);
    printTeam7EdgeReconBoundaryCompareReport(records, params.edge_recon_boundary_mode_, boundary_width_m);
    printTeam7EdgeReconBoundaryHoleSummary(records);
    writeTeam7EdgeReconBoundaryCompareCsv(csv_path, records, params.edge_recon_boundary_mode_, boundary_width_m, dp);
    return records;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_TEAM7_EDGE_RECON_BOUNDARY_H
