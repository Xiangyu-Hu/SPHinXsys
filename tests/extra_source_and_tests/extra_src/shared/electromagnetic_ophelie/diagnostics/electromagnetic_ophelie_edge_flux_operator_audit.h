#ifndef ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_OPERATOR_AUDIT_H
#define ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_OPERATOR_AUDIT_H

#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_parameters.h"
#include "interaction_ck.h"
#include "update_body_relation.h"
#include "vector_functions.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

namespace
{
constexpr const char *kOphelieConductanceSumR2 = "OphelieConductanceSumR2";
constexpr const char *kOphelieConductanceCMin = "OphelieConductanceCMin";
constexpr const char *kOphelieConductanceCMax = "OphelieConductanceCMax";
constexpr const char *kOphelieConductanceCMean = "OphelieConductanceCMean";
constexpr const char *kOphelieConductanceMoment = "OphelieConductanceMoment";
} // namespace

/** Per-particle conductance moment M_i = Σ_j C_ij r_ij ⊗ r_ij (host-side eigen analysis). */
struct OphelieConductanceMomentParticleStats
{
    Real trace = 0.0;
    Real eig1 = 0.0;
    Real eig2 = 0.0;
    Real eig3 = 0.0;
    Real anisotropy = 0.0;
    Real trace_over_sigma_vol = 0.0;
};

struct OphelieEdgeFluxConductanceAudit
{
    Real c_ij_global_min = 0.0;
    Real c_ij_global_max = 0.0;
    Real c_ij_global_mean = 0.0;
    size_t pair_count = 0;
    Real conductance_sum_r2_min = 0.0;
    Real conductance_sum_r2_max = 0.0;
    Real conductance_sum_r2_mean = 0.0;
    Real sigma_vol_min = 0.0;
    Real sigma_vol_max = 0.0;
    Real sigma_vol_mean = 0.0;
    Real conductance_ratio_min = 0.0;
    Real conductance_ratio_max = 0.0;
    Real conductance_ratio_mean = 0.0;
    Real conductance_ratio_median = 0.0;
    Real moment_trace_mean = 0.0;
    Real moment_trace_median = 0.0;
    Real moment_eig1_mean = 0.0;
    Real moment_eig2_mean = 0.0;
    Real moment_eig3_mean = 0.0;
    Real moment_anisotropy_mean = 0.0;
    Real moment_anisotropy_median = 0.0;
    Real typical_spacing_m = 0.0;
    Real reference_inv_h3 = 0.0;
    Real ratio_median_over_inv_h3 = 0.0;
    Real trace_sum_r2_max_rel_err = 0.0;
    size_t n_particles = 0;
};

enum class Team7EdgeFluxPartitionKind
{
    /** Signed-distance far from boundary; use for bulk EMF metrics. */
    TrueInterior,
    /** |SignedDistance| <= boundary shell width. */
    BoundaryShellAll,
    /** Near two or more surface classes (edge/corner). */
    CornerOrEdge,
    /** Legacy: geometric interior bucket (may include mislabeled shell particles). */
    Interior,
    TopSurface,
    BottomSurface,
    ExteriorLateral,
    HoleLateral,
    All
};

struct Team7ParticlePartitionInfo
{
    Team7EdgeFluxPartitionKind surface_kind = Team7EdgeFluxPartitionKind::Interior;
    bool on_boundary_shell = false;
    bool is_corner_or_edge = false;
    bool is_true_interior = false;
};

struct Team7PlateHoleGrid
{
    Real x_lower_m = 0.0;
    Real y_lower_m = 0.0;
    Real cell_h_m = 0.0;
    size_t nx = 0;
    size_t ny = 0;
    StdVec<bool> occupied;
    StdVec<bool> exterior_void;
    StdVec<bool> hole_void;
};

inline Team7PlateHoleGrid buildTeam7PlateHoleGridFromParticles(BaseParticles &particles, const BoundingBoxd &plate_bbox,
                                                               Real cell_h_m)
{
    Team7PlateHoleGrid grid;
    if (cell_h_m <= TinyReal)
    {
        return grid;
    }
    grid.x_lower_m = plate_bbox.lower_[0];
    grid.y_lower_m = plate_bbox.lower_[1];
    grid.cell_h_m = cell_h_m;
    const Real x_span = plate_bbox.upper_[0] - plate_bbox.lower_[0];
    const Real y_span = plate_bbox.upper_[1] - plate_bbox.lower_[1];
    grid.nx = std::max(static_cast<size_t>(1), static_cast<size_t>(std::ceil(x_span / cell_h_m)));
    grid.ny = std::max(static_cast<size_t>(1), static_cast<size_t>(std::ceil(y_span / cell_h_m)));
    const size_t n_cells = grid.nx * grid.ny;
    grid.occupied.assign(n_cells, false);
    grid.exterior_void.assign(n_cells, false);
    grid.hole_void.assign(n_cells, false);

    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const size_t n = particles.TotalRealParticles();
    auto cell_index = [&](Real x, Real y) -> size_t
    {
        size_t ix = static_cast<size_t>(std::floor((x - grid.x_lower_m) / cell_h_m));
        size_t iy = static_cast<size_t>(std::floor((y - grid.y_lower_m) / cell_h_m));
        ix = std::min(ix, grid.nx - 1);
        iy = std::min(iy, grid.ny - 1);
        return iy * grid.nx + ix;
    };
    for (size_t i = 0; i < n; ++i)
    {
        grid.occupied[cell_index(pos[i][0], pos[i][1])] = true;
    }
    StdVec<size_t> queue;
    auto try_push = [&](size_t ix, size_t iy)
    {
        if (ix >= grid.nx || iy >= grid.ny)
        {
            return;
        }
        const size_t idx = iy * grid.nx + ix;
        if (grid.occupied[idx] || grid.exterior_void[idx])
        {
            return;
        }
        grid.exterior_void[idx] = true;
        queue.push_back(idx);
    };
    for (size_t ix = 0; ix < grid.nx; ++ix)
    {
        try_push(ix, 0);
        try_push(ix, grid.ny - 1);
    }
    for (size_t iy = 0; iy < grid.ny; ++iy)
    {
        try_push(0, iy);
        try_push(grid.nx - 1, iy);
    }
    for (size_t head = 0; head < queue.size(); ++head)
    {
        const size_t idx = queue[head];
        const size_t ix = idx % grid.nx;
        const size_t iy = idx / grid.nx;
        if (ix > 0)
        {
            try_push(ix - 1, iy);
        }
        if (ix + 1 < grid.nx)
        {
            try_push(ix + 1, iy);
        }
        if (iy > 0)
        {
            try_push(ix, iy - 1);
        }
        if (iy + 1 < grid.ny)
        {
            try_push(ix, iy + 1);
        }
    }
    for (size_t idx = 0; idx < n_cells; ++idx)
    {
        grid.hole_void[idx] = !grid.occupied[idx] && !grid.exterior_void[idx];
    }
    return grid;
}

inline Real team7DistanceToHoleVoidBoundaryM(const Team7PlateHoleGrid &grid, Real x_m, Real y_m)
{
    if (grid.nx == 0 || grid.ny == 0)
    {
        return std::numeric_limits<Real>::max();
    }
    Real min_dist = std::numeric_limits<Real>::max();
    for (size_t iy = 0; iy < grid.ny; ++iy)
    {
        for (size_t ix = 0; ix < grid.nx; ++ix)
        {
            const size_t idx = iy * grid.nx + ix;
            if (!grid.hole_void[idx])
            {
                continue;
            }
            const Real cx = grid.x_lower_m + (static_cast<Real>(ix) + Real(0.5)) * grid.cell_h_m;
            const Real cy = grid.y_lower_m + (static_cast<Real>(iy) + Real(0.5)) * grid.cell_h_m;
            const Real dx = x_m - cx;
            const Real dy = y_m - cy;
            min_dist = std::min(min_dist, std::sqrt(dx * dx + dy * dy));
        }
    }
    return min_dist;
}

struct Team7EdgeFluxPartitionRecord
{
    std::string partition_name;
    size_t n_particles = 0;
    Real vol_sum_m3 = 0.0;
    Real j_imag_vol_norm = 0.0;
    Real e_imag_vol_norm = 0.0;
    Real b_ind_vol_norm = 0.0;
    Real b_coil_vol_norm = 0.0;
    Real bind_over_bcoil = 0.0;
    Real grad_phi_vol_norm = 0.0;
    Real omega_a_vol_norm = 0.0;
    Real e_edge_em_mismatch = 0.0;
    Real edge_residual_vol_norm = 0.0;
    Real conductance_ratio_mean = 0.0;
    Real j_normal_l2 = 0.0;
    Real j_tangential_l2 = 0.0;
    Real jn_over_jt = 0.0;
    Real e_normal_l2 = 0.0;
    Real moment_trace_mean = 0.0;
    Real moment_eig1_mean = 0.0;
    Real moment_eig2_mean = 0.0;
    Real moment_eig3_mean = 0.0;
    Real moment_anisotropy_mean = 0.0;
    Real moment_trace_over_sigma_vol_mean = 0.0;
    /** |ωA|_vol / |∇φ|_vol — EMF balance partition indicator. */
    Real omega_a_over_grad_phi = 0.0;
    /** Vol-fraction of particles with edge LS fallback (ill-conditioned neighbor set). */
    Real edge_recon_fallback_fraction = 0.0;
    /** Vol-weighted RMS pair edge_drop (imag chain). */
    Real edge_drop_l2 = 0.0;
    /** |J − σE|_vol / |J|_vol — σE=J closure on partition. */
    Real j_sigma_e_mismatch = 0.0;
    /** |J_edge_recon − J|_vol / |J|_vol — edge recon vs stored J. */
    Real j_edge_recon_mismatch = 0.0;
};

template <typename... RelationTypes>
class ComputeOphelieEdgeFluxConductanceAuditCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieEdgeFluxConductanceAuditCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieEdgeFluxConductanceAuditCK(RelationType<Parameters...> &relation, const OphelieGlassFieldNames &names,
                                             Real pair_weight_regularization)
        : BaseInteraction(relation), pair_weight_regularization_(pair_weight_regularization),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_sum_r2_(this->particles_->template getVariableByName<Real>(kOphelieConductanceSumR2)),
          dv_c_min_(this->particles_->template getVariableByName<Real>(kOphelieConductanceCMin)),
          dv_c_max_(this->particles_->template getVariableByName<Real>(kOphelieConductanceCMax)),
          dv_c_mean_(this->particles_->template getVariableByName<Real>(kOphelieConductanceCMean)),
          dv_moment_(this->particles_->template getVariableByName<Matd>(kOphelieConductanceMoment))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              sum_r2_(encloser.dv_sum_r2_->DelegatedData(ex_policy)),
              c_min_(encloser.dv_c_min_->DelegatedData(ex_policy)),
              c_max_(encloser.dv_c_max_->DelegatedData(ex_policy)),
              c_mean_(encloser.dv_c_mean_->DelegatedData(ex_policy)),
              moment_(encloser.dv_moment_->DelegatedData(ex_policy)),
              pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real sigma_i = sigma_[index_i];
            Real sum_r2 = 0.0;
            Matd moment = Matd::Zero();
            Real c_min = std::numeric_limits<Real>::max();
            Real c_max = 0.0;
            Real c_sum = 0.0;
            Real neighbor_count = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Real distance = r_ij_vec.norm();
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real c_ij = computeOphelieEdgeFluxPairConductance(
                    sigma_i, sigma_[index_j], dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                    reference_smoothing_length_);
                sum_r2 += c_ij * distance_sq;
                moment += c_ij * r_ij_vec * r_ij_vec.transpose();
                c_min = std::min(c_min, c_ij);
                c_max = std::max(c_max, c_ij);
                c_sum += c_ij;
                neighbor_count += 1.0;
            }
            sum_r2_[index_i] = sum_r2;
            moment_[index_i] = moment;
            c_min_[index_i] = neighbor_count > 0.0 ? c_min : 0.0;
            c_max_[index_i] = c_max;
            c_mean_[index_i] = neighbor_count > 0.0 ? c_sum / neighbor_count : 0.0;
        }

      protected:
        Real *Vol_;
        Real *sigma_;
        Real *sum_r2_;
        Real *c_min_;
        Real *c_max_;
        Real *c_mean_;
        Matd *moment_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real pair_weight_regularization_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_sum_r2_;
    DiscreteVariable<Real> *dv_c_min_;
    DiscreteVariable<Real> *dv_c_max_;
    DiscreteVariable<Real> *dv_c_mean_;
    DiscreteVariable<Matd> *dv_moment_;
};

template <typename... Parameters>
class ComputeOphelieEdgeFluxConductanceAuditCK<Inner<Parameters...>>
    : public ComputeOphelieEdgeFluxConductanceAuditCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieEdgeFluxConductanceAuditCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names,
                                             Real pair_weight_regularization)
        : ComputeOphelieEdgeFluxConductanceAuditCK<Base, Inner<Parameters...>>(inner_relation, names,
                                                                               pair_weight_regularization)
    {
    }
};

inline Matd ophelieZeroMatd()
{
    return Matd::Zero();
}

inline void registerOphelieConductanceAuditFields(BaseParticles &particles)
{
    particles.registerStateVariable<Real>(kOphelieConductanceSumR2, Real(0));
    particles.registerStateVariable<Real>(kOphelieConductanceCMin, Real(0));
    particles.registerStateVariable<Real>(kOphelieConductanceCMax, Real(0));
    particles.registerStateVariable<Real>(kOphelieConductanceCMean, Real(0));
    particles.registerStateVariable<Matd>(kOphelieConductanceMoment, ophelieZeroMatd());
}

inline OphelieConductanceMomentParticleStats computeOphelieConductanceMomentParticleStats(
    const Matd &moment, Real sum_r2, Real sigma_vol)
{
    OphelieConductanceMomentParticleStats stats;
    stats.trace = moment.trace();
    const Vec3d eigs = getPrincipalValuesFromMatrix(moment);
    stats.eig1 = eigs[0];
    stats.eig2 = eigs[1];
    stats.eig3 = eigs[2];
    stats.anisotropy = (stats.eig1 - stats.eig3) / (stats.trace + TinyReal);
    stats.trace_over_sigma_vol = stats.trace / (sigma_vol + TinyReal);
    (void)sum_r2;
    return stats;
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxConductanceAudit evaluateOphelieEdgeFluxConductanceAudit(
    SolidBody &body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params)
{
    registerOphelieConductanceAuditFields(body.getBaseParticles());
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieEdgeFluxConductanceAuditCK<Inner<>>> compute_audit(
        inner, names, params.pair_weight_regularization_);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_audit.exec();

    BaseParticles &particles = body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Real>(particles, kOphelieConductanceSumR2);
    syncVariableToHost<Real>(particles, kOphelieConductanceCMin);
    syncVariableToHost<Real>(particles, kOphelieConductanceCMax);
    syncVariableToHost<Real>(particles, kOphelieConductanceCMean);
    syncVariableToHost<Matd>(particles, kOphelieConductanceMoment);
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *sum_r2 = particles.getVariableDataByName<Real>(kOphelieConductanceSumR2);
    const Real *c_min = particles.getVariableDataByName<Real>(kOphelieConductanceCMin);
    const Real *c_max = particles.getVariableDataByName<Real>(kOphelieConductanceCMax);
    const Real *c_mean = particles.getVariableDataByName<Real>(kOphelieConductanceCMean);
    const Matd *moment = particles.getVariableDataByName<Matd>(kOphelieConductanceMoment);
    const Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OphelieEdgeFluxConductanceAudit audit;
    audit.n_particles = n;
    audit.c_ij_global_min = std::numeric_limits<Real>::max();
    audit.c_ij_global_max = 0.0;
    Real c_mean_sum = 0.0;
    Real sum_r2_min = std::numeric_limits<Real>::max();
    Real sum_r2_max = 0.0;
    Real sum_r2_total = 0.0;
    Real sigma_vol_min = std::numeric_limits<Real>::max();
    Real sigma_vol_max = 0.0;
    Real sigma_vol_total = 0.0;
    Real ratio_min = std::numeric_limits<Real>::max();
    Real ratio_max = 0.0;
    Real ratio_total = 0.0;
    StdVec<Real> ratios;
    StdVec<Real> traces;
    StdVec<Real> anisotropies;
    Real eig1_total = 0.0;
    Real eig2_total = 0.0;
    Real eig3_total = 0.0;
    Real spacing_total = 0.0;
    Real trace_err_max = 0.0;
    ratios.reserve(n);
    traces.reserve(n);
    anisotropies.reserve(n);

    for (size_t i = 0; i < n; ++i)
    {
        if (c_min[i] > 0.0)
        {
            audit.c_ij_global_min = std::min(audit.c_ij_global_min, c_min[i]);
            audit.c_ij_global_max = std::max(audit.c_ij_global_max, c_max[i]);
            c_mean_sum += c_mean[i];
            ++audit.pair_count;
        }
        sum_r2_min = std::min(sum_r2_min, sum_r2[i]);
        sum_r2_max = std::max(sum_r2_max, sum_r2[i]);
        sum_r2_total += sum_r2[i];
        const Real sigma_vol = sigma[i] * vol[i];
        sigma_vol_min = std::min(sigma_vol_min, sigma_vol);
        sigma_vol_max = std::max(sigma_vol_max, sigma_vol);
        sigma_vol_total += sigma_vol;
        const Real ratio = sum_r2[i] / (sigma_vol + TinyReal);
        ratio_min = std::min(ratio_min, ratio);
        ratio_max = std::max(ratio_max, ratio);
        ratio_total += ratio;
        ratios.push_back(ratio);

        const OphelieConductanceMomentParticleStats mstats =
            computeOphelieConductanceMomentParticleStats(moment[i], sum_r2[i], sigma_vol);
        traces.push_back(mstats.trace);
        anisotropies.push_back(mstats.anisotropy);
        eig1_total += mstats.eig1;
        eig2_total += mstats.eig2;
        eig3_total += mstats.eig3;
        spacing_total += std::cbrt(vol[i]);
        if (sum_r2[i] > TinyReal)
        {
            trace_err_max = std::max(trace_err_max, std::abs(mstats.trace - sum_r2[i]) / sum_r2[i]);
        }
    }

    if (audit.pair_count == 0)
    {
        audit.c_ij_global_min = 0.0;
    }
    audit.c_ij_global_mean = audit.pair_count > 0 ? c_mean_sum / static_cast<Real>(audit.pair_count) : 0.0;
    audit.conductance_sum_r2_min = n > 0 ? sum_r2_min : 0.0;
    audit.conductance_sum_r2_max = sum_r2_max;
    audit.conductance_sum_r2_mean = n > 0 ? sum_r2_total / static_cast<Real>(n) : 0.0;
    audit.sigma_vol_min = n > 0 ? sigma_vol_min : 0.0;
    audit.sigma_vol_max = sigma_vol_max;
    audit.sigma_vol_mean = n > 0 ? sigma_vol_total / static_cast<Real>(n) : 0.0;
    audit.conductance_ratio_min = n > 0 ? ratio_min : 0.0;
    audit.conductance_ratio_max = ratio_max;
    audit.conductance_ratio_mean = n > 0 ? ratio_total / static_cast<Real>(n) : 0.0;
    if (!ratios.empty())
    {
        std::sort(ratios.begin(), ratios.end());
        audit.conductance_ratio_median = ratios[ratios.size() / 2];
    }
    if (!traces.empty())
    {
        std::sort(traces.begin(), traces.end());
        audit.moment_trace_mean = traces.empty() ? 0.0 : std::accumulate(traces.begin(), traces.end(), Real(0)) /
                                                           static_cast<Real>(traces.size());
        audit.moment_trace_median = traces[traces.size() / 2];
    }
    if (!anisotropies.empty())
    {
        std::sort(anisotropies.begin(), anisotropies.end());
        audit.moment_anisotropy_mean =
            std::accumulate(anisotropies.begin(), anisotropies.end(), Real(0)) / static_cast<Real>(anisotropies.size());
        audit.moment_anisotropy_median = anisotropies[anisotropies.size() / 2];
    }
    audit.moment_eig1_mean = n > 0 ? eig1_total / static_cast<Real>(n) : 0.0;
    audit.moment_eig2_mean = n > 0 ? eig2_total / static_cast<Real>(n) : 0.0;
    audit.moment_eig3_mean = n > 0 ? eig3_total / static_cast<Real>(n) : 0.0;
    audit.typical_spacing_m = n > 0 ? spacing_total / static_cast<Real>(n) : 0.0;
    audit.reference_inv_h3 =
        audit.typical_spacing_m > TinyReal ? 1.0 / (audit.typical_spacing_m * audit.typical_spacing_m * audit.typical_spacing_m)
                                           : 0.0;
    audit.ratio_median_over_inv_h3 =
        audit.reference_inv_h3 > TinyReal ? audit.conductance_ratio_median / audit.reference_inv_h3 : 0.0;
    audit.trace_sum_r2_max_rel_err = trace_err_max;
    return audit;
}

inline void printOphelieEdgeFluxConductanceAudit(const OphelieEdgeFluxConductanceAudit &audit)
{
    std::cout << "[ophelie] P4.1 conductance audit: C_ij[min,max,mean]=" << audit.c_ij_global_min << ","
              << audit.c_ij_global_max << "," << audit.c_ij_global_mean
              << " sum(C*r^2)[min,max,mean]=" << audit.conductance_sum_r2_min << "," << audit.conductance_sum_r2_max
              << "," << audit.conductance_sum_r2_mean << " sigma*Vol[min,max,mean]=" << audit.sigma_vol_min << ","
              << audit.sigma_vol_max << "," << audit.sigma_vol_mean
              << " ratio[min,max,mean,median]=" << audit.conductance_ratio_min << "," << audit.conductance_ratio_max
              << "," << audit.conductance_ratio_mean << "," << audit.conductance_ratio_median
              << " moment_trace[mean,median]=" << audit.moment_trace_mean << "," << audit.moment_trace_median
              << " eig[mean]=" << audit.moment_eig1_mean << "," << audit.moment_eig2_mean << ","
              << audit.moment_eig3_mean << " anisotropy[mean,median]=" << audit.moment_anisotropy_mean << ","
              << audit.moment_anisotropy_median << " h_typ=" << audit.typical_spacing_m
              << " ratio_med/inv_h3=" << audit.ratio_median_over_inv_h3
              << " trace_err_max=" << audit.trace_sum_r2_max_rel_err << std::endl;
}

inline void writeOphelieEdgeFluxConductanceAuditCsv(const std::string &output_path,
                                                    const OphelieEdgeFluxConductanceAudit &audit)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    out << "metric,value\n";
    out << "c_ij_global_min," << audit.c_ij_global_min << "\n";
    out << "c_ij_global_max," << audit.c_ij_global_max << "\n";
    out << "c_ij_global_mean," << audit.c_ij_global_mean << "\n";
    out << "conductance_sum_r2_min," << audit.conductance_sum_r2_min << "\n";
    out << "conductance_sum_r2_max," << audit.conductance_sum_r2_max << "\n";
    out << "conductance_sum_r2_mean," << audit.conductance_sum_r2_mean << "\n";
    out << "sigma_vol_min," << audit.sigma_vol_min << "\n";
    out << "sigma_vol_max," << audit.sigma_vol_max << "\n";
    out << "sigma_vol_mean," << audit.sigma_vol_mean << "\n";
    out << "conductance_ratio_min," << audit.conductance_ratio_min << "\n";
    out << "conductance_ratio_max," << audit.conductance_ratio_max << "\n";
    out << "conductance_ratio_mean," << audit.conductance_ratio_mean << "\n";
    out << "conductance_ratio_median," << audit.conductance_ratio_median << "\n";
    out << "moment_trace_mean," << audit.moment_trace_mean << "\n";
    out << "moment_trace_median," << audit.moment_trace_median << "\n";
    out << "moment_eig1_mean," << audit.moment_eig1_mean << "\n";
    out << "moment_eig2_mean," << audit.moment_eig2_mean << "\n";
    out << "moment_eig3_mean," << audit.moment_eig3_mean << "\n";
    out << "moment_anisotropy_mean," << audit.moment_anisotropy_mean << "\n";
    out << "moment_anisotropy_median," << audit.moment_anisotropy_median << "\n";
    out << "typical_spacing_m," << audit.typical_spacing_m << "\n";
    out << "reference_inv_h3," << audit.reference_inv_h3 << "\n";
    out << "ratio_median_over_inv_h3," << audit.ratio_median_over_inv_h3 << "\n";
    out << "trace_sum_r2_max_rel_err," << audit.trace_sum_r2_max_rel_err << "\n";
    out << "n_particles," << audit.n_particles << "\n";
    std::cout << "[ophelie] P4.1 conductance audit CSV: " << output_path << std::endl;
}

inline Team7EdgeFluxPartitionKind classifyTeam7EdgeFluxPartition(
    const Vecd &pos, const BoundingBoxd &plate_bbox, Real skin_h_m, Real z_lower_m, Real z_upper_m,
    const Team7PlateHoleGrid &hole_grid)
{
    const Real x = pos[0];
    const Real y = pos[1];
    const Real z = pos[2];
    const Real dist_x_lo = x - plate_bbox.lower_[0];
    const Real dist_x_hi = plate_bbox.upper_[0] - x;
    const Real dist_y_lo = y - plate_bbox.lower_[1];
    const Real dist_y_hi = plate_bbox.upper_[1] - y;
    const bool near_top = z >= z_upper_m - skin_h_m;
    const bool near_bottom = z <= z_lower_m + skin_h_m;
    const bool near_exterior =
        dist_x_lo < skin_h_m || dist_x_hi < skin_h_m || dist_y_lo < skin_h_m || dist_y_hi < skin_h_m;
    const Real hole_dist = team7DistanceToHoleVoidBoundaryM(hole_grid, x, y);
    const bool near_hole = hole_dist < skin_h_m;
    if (near_top)
    {
        return Team7EdgeFluxPartitionKind::TopSurface;
    }
    if (near_bottom)
    {
        return Team7EdgeFluxPartitionKind::BottomSurface;
    }
    if (near_hole)
    {
        return Team7EdgeFluxPartitionKind::HoleLateral;
    }
    if (near_exterior)
    {
        return Team7EdgeFluxPartitionKind::ExteriorLateral;
    }
    return Team7EdgeFluxPartitionKind::Interior;
}

inline Team7ParticlePartitionInfo classifyTeam7ParticlePartition(
    const Vecd &pos, const BoundingBoxd &plate_bbox, Real skin_h_m, Real z_lower_m, Real z_upper_m,
    const Team7PlateHoleGrid &hole_grid, Real boundary_width_m, Real signed_distance_m, bool have_signed_distance)
{
    Team7ParticlePartitionInfo info;
    info.surface_kind =
        classifyTeam7EdgeFluxPartition(pos, plate_bbox, skin_h_m, z_lower_m, z_upper_m, hole_grid);

    const Real x = pos[0];
    const Real y = pos[1];
    const Real z = pos[2];
    const Real dist_x_lo = x - plate_bbox.lower_[0];
    const Real dist_x_hi = plate_bbox.upper_[0] - x;
    const Real dist_y_lo = y - plate_bbox.lower_[1];
    const Real dist_y_hi = plate_bbox.upper_[1] - y;
    const bool near_top = z >= z_upper_m - skin_h_m;
    const bool near_bottom = z <= z_lower_m + skin_h_m;
    const bool near_exterior =
        dist_x_lo < skin_h_m || dist_x_hi < skin_h_m || dist_y_lo < skin_h_m || dist_y_hi < skin_h_m;
    const Real hole_dist = team7DistanceToHoleVoidBoundaryM(hole_grid, x, y);
    const bool near_hole = hole_dist < skin_h_m;
    const int n_surface_flags =
        static_cast<int>(near_top) + static_cast<int>(near_bottom) + static_cast<int>(near_hole) +
        static_cast<int>(near_exterior);
    info.is_corner_or_edge = n_surface_flags >= 2;

    if (have_signed_distance && boundary_width_m > TinyReal)
    {
        info.on_boundary_shell = std::abs(signed_distance_m) <= boundary_width_m;
    }
    else
    {
        info.on_boundary_shell = n_surface_flags >= 1;
    }
    info.is_true_interior = info.surface_kind == Team7EdgeFluxPartitionKind::Interior && !info.on_boundary_shell &&
                            !info.is_corner_or_edge;
    return info;
}

inline void team7ParticlePartitionTargets(const Team7ParticlePartitionInfo &info,
                                          StdVec<Team7EdgeFluxPartitionKind> &targets_out)
{
    targets_out.clear();
    if (info.is_true_interior)
    {
        targets_out.push_back(Team7EdgeFluxPartitionKind::TrueInterior);
    }
    if (info.on_boundary_shell)
    {
        targets_out.push_back(Team7EdgeFluxPartitionKind::BoundaryShellAll);
    }
    if (info.is_corner_or_edge)
    {
        targets_out.push_back(Team7EdgeFluxPartitionKind::CornerOrEdge);
    }
    targets_out.push_back(info.surface_kind);
    targets_out.push_back(Team7EdgeFluxPartitionKind::All);
}

inline StdVec<Team7EdgeFluxPartitionKind> team7EdgeFluxPartitionKindsForAudit()
{
    return {Team7EdgeFluxPartitionKind::TrueInterior,     Team7EdgeFluxPartitionKind::BoundaryShellAll,
            Team7EdgeFluxPartitionKind::CornerOrEdge,   Team7EdgeFluxPartitionKind::Interior,
            Team7EdgeFluxPartitionKind::TopSurface,   Team7EdgeFluxPartitionKind::BottomSurface,
            Team7EdgeFluxPartitionKind::ExteriorLateral, Team7EdgeFluxPartitionKind::HoleLateral,
            Team7EdgeFluxPartitionKind::All};
}

inline Vecd team7PartitionNormal(Team7EdgeFluxPartitionKind kind, const Vecd &pos, const BoundingBoxd &plate_bbox,
                                 const Team7PlateHoleGrid &hole_grid)
{
    switch (kind)
    {
    case Team7EdgeFluxPartitionKind::TopSurface:
        return Vecd(0.0, 0.0, 1.0);
    case Team7EdgeFluxPartitionKind::BottomSurface:
        return Vecd(0.0, 0.0, -1.0);
    case Team7EdgeFluxPartitionKind::HoleLateral:
    {
        Real best_dist = std::numeric_limits<Real>::max();
        Vecd toward = Vecd::Zero();
        for (size_t iy = 0; iy < hole_grid.ny; ++iy)
        {
            for (size_t ix = 0; ix < hole_grid.nx; ++ix)
            {
                const size_t idx = iy * hole_grid.nx + ix;
                if (!hole_grid.hole_void[idx])
                {
                    continue;
                }
                const Real cx = hole_grid.x_lower_m + (static_cast<Real>(ix) + Real(0.5)) * hole_grid.cell_h_m;
                const Real cy = hole_grid.y_lower_m + (static_cast<Real>(iy) + Real(0.5)) * hole_grid.cell_h_m;
                const Vecd c(cx, cy, pos[2]);
                const Real d = (c - pos).norm();
                if (d < best_dist)
                {
                    best_dist = d;
                    toward = c - pos;
                }
            }
        }
        return toward.norm() > TinyReal ? toward.normalized() : Vecd(1.0, 0.0, 0.0);
    }
    case Team7EdgeFluxPartitionKind::ExteriorLateral:
    {
        const Real dist_x_lo = pos[0] - plate_bbox.lower_[0];
        const Real dist_x_hi = plate_bbox.upper_[0] - pos[0];
        const Real dist_y_lo = pos[1] - plate_bbox.lower_[1];
        const Real dist_y_hi = plate_bbox.upper_[1] - pos[1];
        const Real min_dist = std::min({dist_x_lo, dist_x_hi, dist_y_lo, dist_y_hi});
        if (min_dist == dist_x_lo)
        {
            return Vecd(-1.0, 0.0, 0.0);
        }
        if (min_dist == dist_x_hi)
        {
            return Vecd(1.0, 0.0, 0.0);
        }
        if (min_dist == dist_y_lo)
        {
            return Vecd(0.0, -1.0, 0.0);
        }
        return Vecd(0.0, 1.0, 0.0);
    }
    default:
        return Vecd(1.0, 0.0, 0.0);
    }
}

inline const char *team7EdgeFluxPartitionName(Team7EdgeFluxPartitionKind kind)
{
    switch (kind)
    {
    case Team7EdgeFluxPartitionKind::TrueInterior:
        return "true_interior";
    case Team7EdgeFluxPartitionKind::BoundaryShellAll:
        return "boundary_shell_all";
    case Team7EdgeFluxPartitionKind::CornerOrEdge:
        return "corner_or_edge";
    case Team7EdgeFluxPartitionKind::Interior:
        return "interior";
    case Team7EdgeFluxPartitionKind::TopSurface:
        return "top_surface";
    case Team7EdgeFluxPartitionKind::BottomSurface:
        return "bottom_surface";
    case Team7EdgeFluxPartitionKind::ExteriorLateral:
        return "exterior_lateral";
    case Team7EdgeFluxPartitionKind::HoleLateral:
        return "hole_lateral";
    default:
        return "all";
    }
}

inline StdVec<Team7EdgeFluxPartitionRecord> computeTeam7EdgeFluxPartitionAudit(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const std::string &j_imag_field, const BoundingBoxd &plate_bbox, Real z_lower_m, Real z_upper_m, Real skin_h_m)
{
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, j_imag_field);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.e_imag);
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Vecd>(particles, names.a_coil_real);
    syncVariableToHost<Vecd>(particles, names.b_ind_imag);
    syncVariableToHost<Vecd>(particles, names.b_coil_real);
    syncVariableToHost<Real>(particles, names.edge_flux_residual_imag);
    syncVariableToHost<Real>(particles, names.edge_recon_fallback);
    syncVariableToHost<Real>(particles, names.edge_drop_sq_mean);
    syncVariableToHost<Real>(particles, kOphelieConductanceSumR2);
    syncVariableToHost<Matd>(particles, kOphelieConductanceMoment);
    syncVariableToHost<Real>(particles, names.sigma);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *j_imag = particles.getVariableDataByName<Vecd>(j_imag_field);
    const Vecd *j_edge_recon = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Vecd *e_imag = particles.getVariableDataByName<Vecd>(names.e_imag);
    const Vecd *e_edge = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
    const Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    const Vecd *b_ind = particles.getVariableDataByName<Vecd>(names.b_ind_imag);
    const Vecd *b_coil = particles.getVariableDataByName<Vecd>(names.b_coil_real);
    const Real *edge_res = particles.getVariableDataByName<Real>(names.edge_flux_residual_imag);
    const Real *fallback = particles.getVariableDataByName<Real>(names.edge_recon_fallback);
    const Real *edge_drop_sq = particles.getVariableDataByName<Real>(names.edge_drop_sq_mean);
    const Real *sum_r2 = particles.getVariableDataByName<Real>(kOphelieConductanceSumR2);
    const Matd *moment = particles.getVariableDataByName<Matd>(kOphelieConductanceMoment);
    const Real *sigma = particles.getVariableDataByName<Real>(names.sigma);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *signed_distance = particles.getVariableDataByName<Real>("SignedDistance");
    const bool have_signed_distance = signed_distance != nullptr;
    const Real boundary_width_m = Real(2) * skin_h_m;
    const Real omega = params.omega();

    const Team7PlateHoleGrid hole_grid = buildTeam7PlateHoleGridFromParticles(particles, plate_bbox, skin_h_m);
    size_t hole_void_cells = 0;
    for (bool is_hole : hole_grid.hole_void)
    {
        if (is_hole)
        {
            ++hole_void_cells;
        }
    }
    std::cout << "[ophelie] TEAM7 hole grid: nx=" << hole_grid.nx << " ny=" << hole_grid.ny
              << " hole_void_cells=" << hole_void_cells << std::endl;

    const StdVec<Team7EdgeFluxPartitionKind> kinds = team7EdgeFluxPartitionKindsForAudit();
    StdVec<Team7EdgeFluxPartitionRecord> records(kinds.size());
    for (size_t k = 0; k < kinds.size(); ++k)
    {
        records[k].partition_name = team7EdgeFluxPartitionName(kinds[k]);
    }

    struct Accum
    {
        size_t count = 0;
        Real vol_sum = 0.0;
        Real j_sq = 0.0;
        Real e_sq = 0.0;
        Real b_ind_sq = 0.0;
        Real b_coil_sq = 0.0;
        Real grad_phi_sq = 0.0;
        Real omega_a_sq = 0.0;
        Real edge_mis_sq = 0.0;
        Real edge_norm_sq = 0.0;
        Real edge_res_sq = 0.0;
        Real ratio_sum = 0.0;
        Real j_normal_sq = 0.0;
        Real j_tangent_sq = 0.0;
        Real e_normal_sq = 0.0;
        Real moment_trace_sum = 0.0;
        Real moment_eig1_sum = 0.0;
        Real moment_eig2_sum = 0.0;
        Real moment_eig3_sum = 0.0;
        Real moment_anisotropy_sum = 0.0;
        Real moment_trace_over_sigma_vol_sum = 0.0;
        Real fallback_vol = 0.0;
        Real edge_drop_sq_sum = 0.0;
        Real j_sigma_mis_sq = 0.0;
        Real j_edge_mis_sq = 0.0;
    };
    StdVec<Accum> acc(records.size());

    for (size_t i = 0; i < n; ++i)
    {
        const Team7ParticlePartitionInfo pinfo = classifyTeam7ParticlePartition(
            pos[i], plate_bbox, skin_h_m, z_lower_m, z_upper_m, hole_grid, boundary_width_m,
            have_signed_distance ? signed_distance[i] : Real(0), have_signed_distance);
        const Team7EdgeFluxPartitionKind kind = pinfo.surface_kind;
        const Vecd n_hat = team7PartitionNormal(kind, pos[i], plate_bbox, hole_grid);
        const Vecd j_vec = j_imag[i];
        const Vecd j_sigma_e = sigma[i] * e_imag[i];
        const Vecd j_edge_delta = j_edge_recon[i] - j_vec;
        const Real jn = j_vec.dot(n_hat);
        const Real jt_sq = (j_vec - jn * n_hat).squaredNorm();
        const Real en = e_imag[i].dot(n_hat);
        StdVec<Team7EdgeFluxPartitionKind> targets;
        team7ParticlePartitionTargets(pinfo, targets);
        const Real v = vol[i];
        const Vecd omega_a = omega * a_coil_real[i];
        const Vecd e_em = -grad_phi[i] - omega_a;
        const Vecd edge_delta = e_edge[i] - e_em;
        const Real sigma_vol = sigma[i] * v;
        const Real ratio = sum_r2[i] / (sigma_vol + TinyReal);
        const OphelieConductanceMomentParticleStats mstats =
            computeOphelieConductanceMomentParticleStats(moment[i], sum_r2[i], sigma_vol);

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
            a.j_sq += v * j_imag[i].squaredNorm();
            a.e_sq += v * e_imag[i].squaredNorm();
            a.b_ind_sq += v * b_ind[i].squaredNorm();
            a.b_coil_sq += v * b_coil[i].squaredNorm();
            a.grad_phi_sq += v * grad_phi[i].squaredNorm();
            a.omega_a_sq += v * omega_a.squaredNorm();
            a.edge_mis_sq += v * edge_delta.squaredNorm();
            a.edge_norm_sq += v * e_edge[i].squaredNorm();
            a.edge_res_sq += v * edge_res[i] * edge_res[i];
            a.ratio_sum += ratio;
            a.j_normal_sq += v * jn * jn;
            a.j_tangent_sq += v * jt_sq;
            a.e_normal_sq += v * en * en;
            a.moment_trace_sum += mstats.trace;
            a.moment_eig1_sum += mstats.eig1;
            a.moment_eig2_sum += mstats.eig2;
            a.moment_eig3_sum += mstats.eig3;
            a.moment_anisotropy_sum += mstats.anisotropy;
            a.moment_trace_over_sigma_vol_sum += mstats.trace_over_sigma_vol;
            if (fallback[i] > 0.5)
            {
                a.fallback_vol += v;
            }
            a.edge_drop_sq_sum += v * edge_drop_sq[i];
            a.j_sigma_mis_sq += v * (j_vec - j_sigma_e).squaredNorm();
            a.j_edge_mis_sq += v * j_edge_delta.squaredNorm();
        }
    }

    for (size_t k = 0; k < records.size(); ++k)
    {
        const Accum &a = acc[k];
        records[k].n_particles = a.count;
        records[k].vol_sum_m3 = a.vol_sum;
        records[k].j_imag_vol_norm = std::sqrt(a.j_sq);
        records[k].e_imag_vol_norm = std::sqrt(a.e_sq);
        records[k].b_ind_vol_norm = std::sqrt(a.b_ind_sq);
        records[k].b_coil_vol_norm = std::sqrt(a.b_coil_sq);
        records[k].bind_over_bcoil = records[k].b_ind_vol_norm / (records[k].b_coil_vol_norm + TinyReal);
        records[k].grad_phi_vol_norm = std::sqrt(a.grad_phi_sq);
        records[k].omega_a_vol_norm = std::sqrt(a.omega_a_sq);
        records[k].e_edge_em_mismatch = std::sqrt(a.edge_mis_sq) / (std::sqrt(a.edge_norm_sq) + TinyReal);
        records[k].edge_residual_vol_norm = std::sqrt(a.edge_res_sq);
        records[k].conductance_ratio_mean = a.count > 0 ? a.ratio_sum / static_cast<Real>(a.count) : 0.0;
        records[k].j_normal_l2 = std::sqrt(a.j_normal_sq);
        records[k].j_tangential_l2 = std::sqrt(a.j_tangent_sq);
        records[k].jn_over_jt = records[k].j_normal_l2 / (records[k].j_tangential_l2 + TinyReal);
        records[k].e_normal_l2 = std::sqrt(a.e_normal_sq);
        records[k].moment_trace_mean = a.count > 0 ? a.moment_trace_sum / static_cast<Real>(a.count) : 0.0;
        records[k].moment_eig1_mean = a.count > 0 ? a.moment_eig1_sum / static_cast<Real>(a.count) : 0.0;
        records[k].moment_eig2_mean = a.count > 0 ? a.moment_eig2_sum / static_cast<Real>(a.count) : 0.0;
        records[k].moment_eig3_mean = a.count > 0 ? a.moment_eig3_sum / static_cast<Real>(a.count) : 0.0;
        records[k].moment_anisotropy_mean = a.count > 0 ? a.moment_anisotropy_sum / static_cast<Real>(a.count) : 0.0;
        records[k].moment_trace_over_sigma_vol_mean =
            a.count > 0 ? a.moment_trace_over_sigma_vol_sum / static_cast<Real>(a.count) : 0.0;
        records[k].omega_a_over_grad_phi =
            records[k].grad_phi_vol_norm > TinyReal ? records[k].omega_a_vol_norm / records[k].grad_phi_vol_norm : 0.0;
        records[k].edge_recon_fallback_fraction = a.vol_sum > TinyReal ? a.fallback_vol / a.vol_sum : 0.0;
        records[k].edge_drop_l2 = a.vol_sum > TinyReal ? std::sqrt(a.edge_drop_sq_sum / a.vol_sum) : 0.0;
        records[k].j_sigma_e_mismatch =
            records[k].j_imag_vol_norm > TinyReal ? std::sqrt(a.j_sigma_mis_sq) / records[k].j_imag_vol_norm : 0.0;
        records[k].j_edge_recon_mismatch =
            records[k].j_imag_vol_norm > TinyReal ? std::sqrt(a.j_edge_mis_sq) / records[k].j_imag_vol_norm : 0.0;
    }
    return records;
}

inline void printTeam7EdgeFluxPartitionAuditReport(const StdVec<Team7EdgeFluxPartitionRecord> &records, Real skin_h_m)
{
    std::cout << "[ophelie] P4.4 edge-flux partition audit (skin_h_m=" << skin_h_m << "):" << std::endl;
    for (const Team7EdgeFluxPartitionRecord &r : records)
    {
        std::cout << "  " << r.partition_name << " n=" << r.n_particles << " Bind/B=" << r.bind_over_bcoil
                  << " |J|_vol=" << r.j_imag_vol_norm << " e_edge_em_mis=" << r.e_edge_em_mismatch
                  << " |edge_res|_vol=" << r.edge_residual_vol_norm << " C_ratio_mean=" << r.conductance_ratio_mean
                  << " trace/sigmaVol=" << r.moment_trace_over_sigma_vol_mean
                  << " aniso=" << r.moment_anisotropy_mean << " omegaA/gradPhi=" << r.omega_a_over_grad_phi
                  << " fallback_frac=" << r.edge_recon_fallback_fraction << " edge_drop_l2=" << r.edge_drop_l2
                  << " j_sigmaE_mis=" << r.j_sigma_e_mismatch << " j_edge_mis=" << r.j_edge_recon_mismatch
                  << " Jn/Jt=" << r.jn_over_jt << " |Jn|=" << r.j_normal_l2 << " |Jt|=" << r.j_tangential_l2
                  << std::endl;
    }
}

inline void writeTeam7EdgeFluxPartitionAuditCsv(const std::string &output_path,
                                                const StdVec<Team7EdgeFluxPartitionRecord> &records, Real skin_h_m)
{
    namespace fs = std::filesystem;
    const fs::path parent = fs::path(output_path).parent_path();
    if (!parent.empty())
    {
        fs::create_directories(parent);
    }
    std::ofstream out(output_path);
    out << "partition,skin_h_m,n_particles,vol_sum_m3,j_imag_vol_norm,e_imag_vol_norm,b_ind_vol_norm,"
           "b_coil_vol_norm,Bind_over_Bcoil,grad_phi_vol_norm,omega_a_vol_norm,e_edge_em_mismatch,"
           "edge_residual_vol_norm,conductance_ratio_mean,j_normal_l2,j_tangential_l2,Jn_over_Jt,e_normal_l2,"
           "moment_trace_mean,moment_eig1_mean,moment_eig2_mean,moment_eig3_mean,moment_anisotropy_mean,"
           "moment_trace_over_sigma_vol_mean,omega_a_over_grad_phi,edge_recon_fallback_fraction,edge_drop_l2,"
           "j_sigma_e_mismatch,j_edge_recon_mismatch\n";
    for (const Team7EdgeFluxPartitionRecord &r : records)
    {
        out << r.partition_name << "," << skin_h_m << "," << r.n_particles << "," << r.vol_sum_m3 << ","
            << r.j_imag_vol_norm << "," << r.e_imag_vol_norm << "," << r.b_ind_vol_norm << "," << r.b_coil_vol_norm
            << "," << r.bind_over_bcoil << "," << r.grad_phi_vol_norm << "," << r.omega_a_vol_norm << ","
            << r.e_edge_em_mismatch << "," << r.edge_residual_vol_norm << "," << r.conductance_ratio_mean << ","
            << r.j_normal_l2 << "," << r.j_tangential_l2 << "," << r.jn_over_jt << "," << r.e_normal_l2 << ","
            << r.moment_trace_mean << "," << r.moment_eig1_mean << "," << r.moment_eig2_mean << ","
            << r.moment_eig3_mean << "," << r.moment_anisotropy_mean << ","
            << r.moment_trace_over_sigma_vol_mean << "," << r.omega_a_over_grad_phi << ","
            << r.edge_recon_fallback_fraction << "," << r.edge_drop_l2 << "," << r.j_sigma_e_mismatch << ","
            << r.j_edge_recon_mismatch << "\n";
    }
    std::cout << "[ophelie] P4.4 partition audit CSV: " << output_path << std::endl;
}

inline void printTeam7HoleLateralEdgeReconSummary(const StdVec<Team7EdgeFluxPartitionRecord> &records)
{
    const Team7EdgeFluxPartitionRecord *hole = nullptr;
    const Team7EdgeFluxPartitionRecord *interior = nullptr;
    for (const Team7EdgeFluxPartitionRecord &r : records)
    {
        if (r.partition_name == "hole_lateral")
        {
            hole = &r;
        }
        if (r.partition_name == "interior")
        {
            interior = &r;
        }
    }
    if (hole == nullptr || interior == nullptr)
    {
        return;
    }
    std::cout << "[ophelie] P4.5 hole_lateral edge-recon vs interior:"
              << " e_edge_em_mis=" << hole->e_edge_em_mismatch << " vs " << interior->e_edge_em_mismatch
              << " omegaA/gradPhi=" << hole->omega_a_over_grad_phi << " vs " << interior->omega_a_over_grad_phi
              << " fallback_frac=" << hole->edge_recon_fallback_fraction << " vs "
              << interior->edge_recon_fallback_fraction << " edge_drop_l2=" << hole->edge_drop_l2 << " vs "
              << interior->edge_drop_l2 << " j_sigmaE_mis=" << hole->j_sigma_e_mismatch << " vs "
              << interior->j_sigma_e_mismatch << " Bind/B=" << hole->bind_over_bcoil << " vs "
              << interior->bind_over_bcoil << std::endl;
}

template <class ExecutionPolicy>
inline void runTeam7EdgeFluxOperatorAudit(SolidBody &plate_body, Inner<> &plate_inner,
                                          const OphelieGlassFieldNames &plate_names, OphelieParameters &params,
                                          const std::string &j_imag_field, const BoundingBoxd &plate_bbox,
                                          Real z_lower_m, Real z_upper_m, Real skin_h_m,
                                          const std::string &conductance_csv_path,
                                          const std::string &partition_csv_path)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(plate_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(plate_inner);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(plate_inner,
                                                                                                        plate_names);
    compute_grad_phi.exec();
    (void)evaluateOphelieEdgeFluxResidualForComponent<ExecutionPolicy>(
        plate_body, plate_inner, plate_names, makeOphelieEdgeFluxImagComponent(plate_names, params), params);
    (void)evaluateOphelieEdgeFluxImagEdgeDropMetrics<ExecutionPolicy>(plate_body, plate_inner, plate_names, params);

    const OphelieEdgeFluxConductanceAudit conductance_audit =
        evaluateOphelieEdgeFluxConductanceAudit<ExecutionPolicy>(plate_body, plate_inner, plate_names, params);
    printOphelieEdgeFluxConductanceAudit(conductance_audit);
    writeOphelieEdgeFluxConductanceAuditCsv(conductance_csv_path, conductance_audit);

    const StdVec<Team7EdgeFluxPartitionRecord> partition_records = computeTeam7EdgeFluxPartitionAudit(
        plate_body.getBaseParticles(), plate_names, params, j_imag_field, plate_bbox, z_lower_m, z_upper_m, skin_h_m);
    printTeam7EdgeFluxPartitionAuditReport(partition_records, skin_h_m);
    printTeam7HoleLateralEdgeReconSummary(partition_records);
    writeTeam7EdgeFluxPartitionAuditCsv(partition_csv_path, partition_records, skin_h_m);
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_OPERATOR_AUDIT_H
