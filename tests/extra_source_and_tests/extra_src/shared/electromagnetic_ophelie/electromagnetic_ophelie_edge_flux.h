#ifndef ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_H
#define ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux_boundary_closure.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_parameters.h"
#include "interaction_ck.h"
#include "update_body_relation.h"
#include "vector_functions.h"

#include <cmath>
#include <iostream>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

struct OphelieEdgeFluxResidualMetrics
{
    Real edge_res_l2 = 0.0;
    Real edge_res_linf = 0.0;
    Real edge_res_volume_integral = 0.0;
    Real global_conservation = 0.0;
    int a_flux_sign = 1;
};

struct OphelieEdgeFluxDiagnosticReport
{
    OphelieEdgeFluxResidualMetrics level0;
    OphelieEdgeFluxResidualMetrics post_phi;
    Real edge_res_red_l2 = 0.0;
};

inline bool ophelieUseEdgeFluxElectromotiveRhs(const OphelieParameters &params)
{
    return params.ophelie_current_form_ == OphelieCurrentFormKind::EdgeFlux;
}

/** Active A_real for edge-flux kernels: coil-only (A_coil) or total (A_coil+A_ind). */
inline const std::string &getOphelieActiveARealFieldName(const OphelieGlassFieldNames &names,
                                                         const OphelieParameters &params)
{
    return params.use_a_total_for_edge_flux_ ? names.a_src_real : names.a_coil_real;
}

/** Active A_imag for real-chain edge-flux: coil-only or total A_imag. */
inline const std::string &getOphelieActiveAImagFieldName(const OphelieGlassFieldNames &names,
                                                           const OphelieParameters &params)
{
    return params.use_a_total_for_edge_flux_ ? names.a_src_imag : names.a_coil_imag;
}

/** One scalar edge-flux phasor component (imag or real chain). */
struct OphelieEdgeFluxComponent
{
    std::string phi_field;
    std::string phi_rhs_field;
    std::string phi_lhs_field;
    std::string active_a_field;
    std::string edge_residual_field;
    std::string e_recon_field;
    std::string j_recon_field;
    std::string joule_heat_recon_field;
    Real a_sign = Real(1.0);
    const char *chain_label = "imag";
};

using OphelieEdgeFluxEdgeDropComponent = OphelieEdgeFluxComponent;

inline OphelieEdgeFluxComponent makeOphelieEdgeFluxImagComponent(const OphelieGlassFieldNames &names,
                                                                 const OphelieParameters &params)
{
    OphelieEdgeFluxComponent component;
    component.phi_field = names.phi_imag;
    component.phi_rhs_field = names.phi_rhs_imag;
    component.phi_lhs_field = names.phi_lhs_imag;
    component.active_a_field = getOphelieActiveARealFieldName(names, params);
    component.edge_residual_field = names.edge_flux_residual_imag;
    component.e_recon_field = names.e_edge_recon_imag;
    component.j_recon_field = names.j_edge_recon_imag;
    component.joule_heat_recon_field = names.joule_heat_edge_recon_imag;
    component.a_sign = params.edge_flux_imag_a_sign_;
    component.chain_label = "imag";
    return component;
}

inline OphelieEdgeFluxComponent makeOphelieEdgeFluxRealComponent(const OphelieGlassFieldNames &names,
                                                                 const OphelieParameters &params)
{
    OphelieEdgeFluxComponent component;
    component.phi_field = names.phi_real;
    component.phi_rhs_field = names.phi_rhs_real;
    component.phi_lhs_field = names.phi_lhs_real;
    component.active_a_field = getOphelieActiveAImagFieldName(names, params);
    component.edge_residual_field = names.edge_flux_residual_real;
    component.e_recon_field = names.e_edge_recon_real;
    component.j_recon_field = names.j_edge_recon_real;
    component.joule_heat_recon_field = names.joule_heat_edge_recon_real;
    component.a_sign = Real(-1.0);
    component.chain_label = "real";
    return component;
}

inline OphelieEdgeFluxComponent makeOphelieEdgeFluxImagEdgeDropComponent(const OphelieGlassFieldNames &names,
                                                                         const OphelieParameters &params)
{
    return makeOphelieEdgeFluxImagComponent(names, params);
}

inline OphelieEdgeFluxComponent makeOphelieEdgeFluxRealEdgeDropComponent(const OphelieGlassFieldNames &names,
                                                                         const OphelieParameters &params)
{
    return makeOphelieEdgeFluxRealComponent(names, params);
}

/** J_imag source for A_ind Biot–Savart: edge-flux uses JEdgeReconImag, else primary JImag. */
inline const std::string &getOphelieAIndJImagFieldName(const OphelieGlassFieldNames &names,
                                                       const OphelieParameters &params)
{
    return ophelieUseEdgeFluxElectromotiveRhs(params) ? names.j_edge_recon_imag : names.j_imag;
}

inline const std::string &getOphelieAIndJRealFieldName(const OphelieGlassFieldNames &names,
                                                       const OphelieParameters &params)
{
    return ophelieUseEdgeFluxElectromotiveRhs(params) ? names.j_edge_recon_real : names.j_real;
}

class CopyOpheliePrimaryEJQToParticleDiagnosticCK : public LocalDynamics
{
  public:
    CopyOpheliePrimaryEJQToParticleDiagnosticCK(SPHBody &sph_body, const OphelieGlassFieldNames &names)
        : LocalDynamics(sph_body),
          dv_e_imag_(particles_->template getVariableByName<Vecd>(names.e_imag)),
          dv_j_imag_(particles_->template getVariableByName<Vecd>(names.j_imag)),
          dv_joule_heat_(particles_->template getVariableByName<Real>(names.joule_heat)),
          dv_e_diag_(particles_->template getVariableByName<Vecd>(names.e_imag_particle_diag)),
          dv_j_diag_(particles_->template getVariableByName<Vecd>(names.j_imag_particle_diag)),
          dv_q_diag_(particles_->template getVariableByName<Real>(names.joule_heat_particle_diag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : e_imag_(encloser.dv_e_imag_->DelegatedData(ex_policy)),
              j_imag_(encloser.dv_j_imag_->DelegatedData(ex_policy)),
              joule_heat_(encloser.dv_joule_heat_->DelegatedData(ex_policy)),
              e_diag_(encloser.dv_e_diag_->DelegatedData(ex_policy)),
              j_diag_(encloser.dv_j_diag_->DelegatedData(ex_policy)),
              q_diag_(encloser.dv_q_diag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            e_diag_[index_i] = e_imag_[index_i];
            j_diag_[index_i] = j_imag_[index_i];
            q_diag_[index_i] = joule_heat_[index_i];
        }

      protected:
        Vecd *e_imag_;
        Vecd *j_imag_;
        Real *joule_heat_;
        Vecd *e_diag_;
        Vecd *j_diag_;
        Real *q_diag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_e_imag_;
    DiscreteVariable<Vecd> *dv_j_imag_;
    DiscreteVariable<Real> *dv_joule_heat_;
    DiscreteVariable<Vecd> *dv_e_diag_;
    DiscreteVariable<Vecd> *dv_j_diag_;
    DiscreteVariable<Real> *dv_q_diag_;
};

class CopyOphelieEdgeReconToPrimaryEJQCK : public LocalDynamics
{
  public:
    CopyOphelieEdgeReconToPrimaryEJQCK(SPHBody &sph_body, const OphelieGlassFieldNames &names, bool complex_mode)
        : LocalDynamics(sph_body), complex_mode_(complex_mode),
          dv_e_edge_imag_(particles_->template getVariableByName<Vecd>(names.e_edge_recon_imag)),
          dv_j_edge_imag_(particles_->template getVariableByName<Vecd>(names.j_edge_recon_imag)),
          dv_q_edge_imag_(particles_->template getVariableByName<Real>(names.joule_heat_edge_recon_imag)),
          dv_e_edge_real_(particles_->template getVariableByName<Vecd>(names.e_edge_recon_real)),
          dv_j_edge_real_(particles_->template getVariableByName<Vecd>(names.j_edge_recon_real)),
          dv_q_complex_(particles_->template getVariableByName<Real>(names.joule_heat_edge_recon_complex)),
          dv_e_imag_(particles_->template getVariableByName<Vecd>(names.e_imag)),
          dv_j_imag_(particles_->template getVariableByName<Vecd>(names.j_imag)),
          dv_e_real_(particles_->template getVariableByName<Vecd>(names.e_real)),
          dv_j_real_(particles_->template getVariableByName<Vecd>(names.j_real)),
          dv_joule_heat_(particles_->template getVariableByName<Real>(names.joule_heat))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : complex_mode_(encloser.complex_mode_), e_edge_imag_(encloser.dv_e_edge_imag_->DelegatedData(ex_policy)),
              j_edge_imag_(encloser.dv_j_edge_imag_->DelegatedData(ex_policy)),
              q_edge_imag_(encloser.dv_q_edge_imag_->DelegatedData(ex_policy)),
              e_edge_real_(encloser.dv_e_edge_real_->DelegatedData(ex_policy)),
              j_edge_real_(encloser.dv_j_edge_real_->DelegatedData(ex_policy)),
              q_complex_(encloser.dv_q_complex_->DelegatedData(ex_policy)),
              e_imag_(encloser.dv_e_imag_->DelegatedData(ex_policy)),
              j_imag_(encloser.dv_j_imag_->DelegatedData(ex_policy)),
              e_real_(encloser.dv_e_real_->DelegatedData(ex_policy)),
              j_real_(encloser.dv_j_real_->DelegatedData(ex_policy)),
              joule_heat_(encloser.dv_joule_heat_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            e_imag_[index_i] = e_edge_imag_[index_i];
            j_imag_[index_i] = j_edge_imag_[index_i];
            if (complex_mode_)
            {
                e_real_[index_i] = e_edge_real_[index_i];
                j_real_[index_i] = j_edge_real_[index_i];
                joule_heat_[index_i] = q_complex_[index_i];
            }
            else
            {
                joule_heat_[index_i] = q_edge_imag_[index_i];
            }
        }

      protected:
        bool complex_mode_;
        Vecd *e_edge_imag_;
        Vecd *j_edge_imag_;
        Real *q_edge_imag_;
        Vecd *e_edge_real_;
        Vecd *j_edge_real_;
        Real *q_complex_;
        Vecd *e_imag_;
        Vecd *j_imag_;
        Vecd *e_real_;
        Vecd *j_real_;
        Real *joule_heat_;
    };

  protected:
    bool complex_mode_;
    DiscreteVariable<Vecd> *dv_e_edge_imag_;
    DiscreteVariable<Vecd> *dv_j_edge_imag_;
    DiscreteVariable<Real> *dv_q_edge_imag_;
    DiscreteVariable<Vecd> *dv_e_edge_real_;
    DiscreteVariable<Vecd> *dv_j_edge_real_;
    DiscreteVariable<Real> *dv_q_complex_;
    DiscreteVariable<Vecd> *dv_e_imag_;
    DiscreteVariable<Vecd> *dv_j_imag_;
    DiscreteVariable<Vecd> *dv_e_real_;
    DiscreteVariable<Vecd> *dv_j_real_;
    DiscreteVariable<Real> *dv_joule_heat_;
};

template <class ExecutionPolicy>
inline void syncOphelieEdgeReconToPrimaryEJQ(SolidBody &glass_body, const OphelieGlassFieldNames &names,
                                             const OphelieParameters &params)
{
    StateDynamics<ExecutionPolicy, CopyOphelieEdgeReconToPrimaryEJQCK> copy_primary(
        glass_body, names, params.edge_flux_complex_);
    copy_primary.exec();
}

inline void logOphelieEdgeFluxProductionFieldPolicy(const OphelieParameters &params)
{
    if (!ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        return;
    }
    std::cout << "[ophelie] edge-flux mode: primary E/J/Q = EEdgeRecon/JEdgeRecon/JouleHeatEdgeRecon; "
              << "particle divJ_L2_red is diagnostic-only and is not a production continuity gate.";
    if (params.output_particle_gradient_diagnostics_)
    {
        std::cout << " Particle-gradient fields saved to EImagParticleDiag/JImagParticleDiag/JouleHeatParticleDiag.";
    }
    std::cout << std::endl;
}

/** Pair conductance C_ij: reuse legacy-pairwise LHS weight (sigma_ij * pairwiseNegativeLaplaceWeight). */
inline Real computeOphelieEdgeFluxPairConductance(const Real sigma_i, const Real sigma_j, const Real dW_ijV_j,
                                                  const Real distance, const Real distance_sq,
                                                  const Real pair_weight_regularization,
                                                  const Real reference_smoothing_length)
{
    const Real sigma_ij = harmonicMean(sigma_i, sigma_j);
    return sigma_ij * pairwiseNegativeLaplaceWeight(dW_ijV_j, distance, distance_sq, pair_weight_regularization,
                                                    reference_smoothing_length);
}

inline Real computeOphelieEdgeFluxEdgeDrop(const Real phi_i, const Real phi_j, const Vecd &a_i, const Vecd &a_j,
                                           const Vecd &xj_minus_xi, const Real omega, const Real a_sign = Real(1.0))
{
    const Vecd a_avg = Real(0.5) * (a_i + a_j);
    return (phi_j - phi_i) + a_sign * omega * a_avg.dot(xj_minus_xi);
}

struct OphelieEdgeFluxPowerMetrics
{
    /** Graph/Laplace discrete energy Σ 0.25*C_ij*edge_drop² — diagnostic only. */
    Real p_graph_edge = 0.0;
    /** Physical Joule power from edge-reconstructed E/J. */
    Real p_total_recon = 0.0;
    Real p_graph_over_recon = 0.0;
    Real joule_heat_edge_max = 0.0;
    Real joule_heat_edge_mean = 0.0;
    Real joule_heat_edge_recon_max = 0.0;
    Real joule_heat_edge_recon_mean = 0.0;
    Real joule_heat_edge_max_over_mean = 0.0;
    Real recon_fallback_fraction = 0.0;
};

struct OphelieEdgeFluxEdgeDropMetrics
{
    Real edge_drop_l2 = 0.0;
    Real edge_drop_linf = 0.0;
    Real edge_drop_mean_abs = 0.0;
    size_t pair_count = 0;
};

template <typename... RelationTypes>
class ComputeOphelieEdgeFluxPhiRhsFromASrcCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieEdgeFluxPhiRhsFromASrcCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieEdgeFluxPhiRhsFromASrcCK(RelationType<Parameters...> &relation, const OphelieGlassFieldNames &names,
                                           const OphelieEdgeFluxComponent &component, Real omega,
                                           Real pair_weight_regularization,
                                           const OphelieEdgeFluxBoundaryClosureFlags &boundary_closure = {})
        : BaseInteraction(relation), omega_(omega), pair_weight_regularization_(pair_weight_regularization),
          a_sign_(component.a_sign), boundary_closure_(boundary_closure),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_active_a_(this->particles_->template getVariableByName<Vecd>(component.active_a_field)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_phi_rhs_(this->particles_->template getVariableByName<Real>(component.phi_rhs_field)),
          dv_signed_distance_(boundary_closure.phi_rhs_ghost
                                  ? this->particles_->template getVariableByName<Real>("SignedDistance")
                                  : nullptr),
          dv_normal_direction_(boundary_closure.phi_rhs_ghost
                                   ? this->particles_->template getVariableByName<Vecd>("NormalDirection")
                                   : nullptr)
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              active_a_(encloser.dv_active_a_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              phi_rhs_(encloser.dv_phi_rhs_->DelegatedData(ex_policy)),
              signed_distance_(encloser.dv_signed_distance_ != nullptr
                                   ? encloser.dv_signed_distance_->DelegatedData(ex_policy)
                                   : nullptr),
              normal_direction_(encloser.dv_normal_direction_ != nullptr
                                    ? encloser.dv_normal_direction_->DelegatedData(ex_policy)
                                    : nullptr),
              boundary_closure_(encloser.boundary_closure_), omega_(encloser.omega_),
              a_sign_(encloser.a_sign_), pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real sigma_i = sigma_[index_i];
            const Vecd a_i = active_a_[index_i];
            Real rhs_i = 0.0;
            double max_neighbor_conductance = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Vecd xj_minus_xi = -r_ij_vec;
                const Real distance = r_ij_vec.norm();
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real c_ij = computeOphelieEdgeFluxPairConductance(
                    sigma_i, sigma_[index_j], dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                    reference_smoothing_length_);
                max_neighbor_conductance = std::max(max_neighbor_conductance, static_cast<double>(std::abs(c_ij)));
                const Vecd a_avg = Real(0.5) * (a_i + active_a_[index_j]);
                rhs_i += c_ij * a_sign_ * omega_ * a_avg.dot(xj_minus_xi);
            }
            if (boundary_closure_.phi_rhs_ghost && signed_distance_ != nullptr && normal_direction_ != nullptr &&
                max_neighbor_conductance > TinyReal)
            {
                Vecd unit_normal = Vecd::Zero();
                if (ophelieBoundaryShellParticle(signed_distance_[index_i], boundary_closure_.boundary_width_m) &&
                    ophelieTryUnitBoundaryNormal(normal_direction_, index_i, unit_normal))
                {
                    const Real ghost_h = boundary_closure_.boundary_width_m;
                    const Vecd r_ig = ghost_h * unit_normal;
                    const Real c_ig = ophelieBoundaryGhostConductanceFromMaxNeighbor(
                        static_cast<Real>(max_neighbor_conductance));
                    rhs_i += c_ig * a_sign_ * omega_ * a_i.dot(r_ig);
                }
            }
            phi_rhs_[index_i] += rhs_i;
        }

      protected:
        Real *Vol_;
        Vecd *active_a_;
        Real *sigma_;
        Real *phi_rhs_;
        Real *signed_distance_;
        Vecd *normal_direction_;
        OphelieEdgeFluxBoundaryClosureFlags boundary_closure_;
        Real omega_;
        Real a_sign_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    Real pair_weight_regularization_;
    Real a_sign_;
    OphelieEdgeFluxBoundaryClosureFlags boundary_closure_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Vecd> *dv_active_a_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_phi_rhs_;
    DiscreteVariable<Real> *dv_signed_distance_;
    DiscreteVariable<Vecd> *dv_normal_direction_;
};

template <typename... Parameters>
class ComputeOphelieEdgeFluxPhiRhsFromASrcCK<Inner<Parameters...>>
    : public ComputeOphelieEdgeFluxPhiRhsFromASrcCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieEdgeFluxPhiRhsFromASrcCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names,
                                           const OphelieEdgeFluxComponent &component, Real omega,
                                           Real pair_weight_regularization,
                                           const OphelieEdgeFluxBoundaryClosureFlags &boundary_closure = {})
        : ComputeOphelieEdgeFluxPhiRhsFromASrcCK<Base, Inner<Parameters...>>(inner_relation, names, component, omega,
                                                                               pair_weight_regularization,
                                                                               boundary_closure)
    {
    }
};

template <typename... RelationTypes>
class ComputeOphelieEdgeFluxResidualCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieEdgeFluxResidualCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieEdgeFluxResidualCK(RelationType<Parameters...> &relation, const OphelieGlassFieldNames &names,
                                     const OphelieEdgeFluxComponent &component, Real omega,
                                     Real pair_weight_regularization)
        : BaseInteraction(relation), omega_(omega), pair_weight_regularization_(pair_weight_regularization),
          a_sign_(component.a_sign), reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_phi_(this->particles_->template getVariableByName<Real>(component.phi_field)),
          dv_active_a_(this->particles_->template getVariableByName<Vecd>(component.active_a_field)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_edge_residual_(this->particles_->template getVariableByName<Real>(component.edge_residual_field))
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
              edge_residual_(encloser.dv_edge_residual_->DelegatedData(ex_policy)), omega_(encloser.omega_),
              a_sign_(encloser.a_sign_), pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_[index_i];
            const Real sigma_i = sigma_[index_i];
            const Vecd a_i = active_a_[index_i];
            Real residual_i = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Vecd xj_minus_xi = -r_ij_vec;
                const Real distance = r_ij_vec.norm();
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real c_ij = computeOphelieEdgeFluxPairConductance(
                    sigma_i, sigma_[index_j], dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                    reference_smoothing_length_);
                const Real edge_drop = computeOphelieEdgeFluxEdgeDrop(phi_i, phi_[index_j], a_i, active_a_[index_j],
                                                                      xj_minus_xi, omega_, a_sign_);
                const Real q_ij = -c_ij * edge_drop;
                residual_i += q_ij;
            }
            edge_residual_[index_i] = residual_i;
        }

      protected:
        Real *Vol_;
        Real *phi_;
        Vecd *active_a_;
        Real *sigma_;
        Real *edge_residual_;
        Real omega_;
        Real a_sign_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    Real pair_weight_regularization_;
    Real a_sign_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Vecd> *dv_active_a_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_edge_residual_;
};

template <typename... Parameters>
class ComputeOphelieEdgeFluxResidualCK<Inner<Parameters...>>
    : public ComputeOphelieEdgeFluxResidualCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieEdgeFluxResidualCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names,
                                     const OphelieEdgeFluxComponent &component, Real omega,
                                     Real pair_weight_regularization)
        : ComputeOphelieEdgeFluxResidualCK<Base, Inner<Parameters...>>(inner_relation, names, component, omega,
                                                                       pair_weight_regularization)
    {
    }
};

class ComputeOphelieEdgeFluxReconstructedJouleHeatCK : public LocalDynamics
{
  public:
    ComputeOphelieEdgeFluxReconstructedJouleHeatCK(SPHBody &sph_body, const std::string &sigma_field,
                                                    const std::string &e_recon_field,
                                                    const std::string &joule_heat_recon_field)
        : LocalDynamics(sph_body),
          dv_sigma_(particles_->template getVariableByName<Real>(sigma_field)),
          dv_e_edge_recon_(particles_->template getVariableByName<Vecd>(e_recon_field)),
          dv_joule_heat_edge_recon_(particles_->template getVariableByName<Real>(joule_heat_recon_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              e_edge_recon_(encloser.dv_e_edge_recon_->DelegatedData(ex_policy)),
              joule_heat_edge_recon_(encloser.dv_joule_heat_edge_recon_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real sigma_i = sigma_[index_i];
            joule_heat_edge_recon_[index_i] = Real(0.5) * sigma_i * e_edge_recon_[index_i].squaredNorm();
        }

      protected:
        Real *sigma_;
        Vecd *e_edge_recon_;
        Real *joule_heat_edge_recon_;
    };

  protected:
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Vecd> *dv_e_edge_recon_;
    DiscreteVariable<Real> *dv_joule_heat_edge_recon_;
};

template <typename... RelationTypes>
class ComputeOphelieEdgeFluxJouleHeatCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieEdgeFluxJouleHeatCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieEdgeFluxJouleHeatCK(RelationType<Parameters...> &relation, const OphelieGlassFieldNames &names,
                                      Real omega, Real pair_weight_regularization, const std::string &a_real_field,
                                      Real a_sign = Real(1.0))
        : BaseInteraction(relation), omega_(omega), pair_weight_regularization_(pair_weight_regularization),
          a_sign_(a_sign), reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_phi_imag_(this->particles_->template getVariableByName<Real>(names.phi_imag)),
          dv_a_src_real_(this->particles_->template getVariableByName<Vecd>(a_real_field)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_joule_heat_edge_(this->particles_->template getVariableByName<Real>(names.joule_heat_edge)),
          dv_power_edge_particle_(this->particles_->template getVariableByName<Real>(names.power_edge_particle))
    {
    }

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType, typename... Args>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, Args &&...args)
            : BaseInteraction::InteractKernel(ex_policy, encloser, std::forward<Args>(args)...),
              Vol_(encloser.dv_Vol_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
              a_src_real_(encloser.dv_a_src_real_->DelegatedData(ex_policy)),
              sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              joule_heat_edge_(encloser.dv_joule_heat_edge_->DelegatedData(ex_policy)),
              power_edge_particle_(encloser.dv_power_edge_particle_->DelegatedData(ex_policy)), omega_(encloser.omega_),
              a_sign_(encloser.a_sign_), pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_imag_[index_i];
            const Real sigma_i = sigma_[index_i];
            const Vecd a_i = a_src_real_[index_i];
            const Real vol_i = Vol_[index_i];
            double power_i = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Vecd xj_minus_xi = -r_ij_vec;
                const Real distance = r_ij_vec.norm();
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real c_ij = computeOphelieEdgeFluxPairConductance(
                    sigma_i, sigma_[index_j], dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                    reference_smoothing_length_);
                const Real edge_drop = computeOphelieEdgeFluxEdgeDrop(phi_i, phi_imag_[index_j], a_i,
                                                                      a_src_real_[index_j], xj_minus_xi, omega_,
                                                                      a_sign_);
                power_i += Real(0.25) * static_cast<double>(c_ij) * static_cast<double>(edge_drop) *
                           static_cast<double>(edge_drop);
            }
            power_edge_particle_[index_i] = static_cast<Real>(power_i);
            joule_heat_edge_[index_i] = static_cast<Real>(power_i / (static_cast<double>(vol_i) + TinyReal));
        }

      protected:
        Real *Vol_;
        Real *phi_imag_;
        Vecd *a_src_real_;
        Real *sigma_;
        Real *joule_heat_edge_;
        Real *power_edge_particle_;
        Real omega_;
        Real a_sign_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    Real pair_weight_regularization_;
    Real a_sign_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_a_src_real_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_joule_heat_edge_;
    DiscreteVariable<Real> *dv_power_edge_particle_;
};

template <typename... Parameters>
class ComputeOphelieEdgeFluxJouleHeatCK<Inner<Parameters...>>
    : public ComputeOphelieEdgeFluxJouleHeatCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieEdgeFluxJouleHeatCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names,
                                      Real omega, Real pair_weight_regularization, const std::string &a_real_field,
                                      Real a_sign = Real(1.0))
        : ComputeOphelieEdgeFluxJouleHeatCK<Base, Inner<Parameters...>>(inner_relation, names, omega,
                                                                        pair_weight_regularization, a_real_field,
                                                                        a_sign)
    {
    }
};

template <typename... RelationTypes>
class ReconstructOphelieEdgeFluxElectricCurrentCK;

template <template <typename...> class RelationType, typename... Parameters>
class ReconstructOphelieEdgeFluxElectricCurrentCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ReconstructOphelieEdgeFluxElectricCurrentCK(RelationType<Parameters...> &relation, const OphelieGlassFieldNames &names,
                                                const OphelieEdgeFluxComponent &component, Real omega,
                                                Real pair_weight_regularization, Real recon_condition_threshold,
                                                Real solver_local_rhs_scale = Real(1),
                                                const OphelieEdgeFluxBoundaryClosureFlags &boundary_closure = {})
        : BaseInteraction(relation), omega_(omega), pair_weight_regularization_(pair_weight_regularization),
          a_sign_(component.a_sign), recon_condition_threshold_(recon_condition_threshold),
          solver_local_rhs_scale_(solver_local_rhs_scale), boundary_closure_(boundary_closure),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_phi_(this->particles_->template getVariableByName<Real>(component.phi_field)),
          dv_active_a_(this->particles_->template getVariableByName<Vecd>(component.active_a_field)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_e_edge_recon_(this->particles_->template getVariableByName<Vecd>(component.e_recon_field)),
          dv_j_edge_recon_(this->particles_->template getVariableByName<Vecd>(component.j_recon_field)),
          dv_edge_recon_condition_(this->particles_->template getVariableByName<Real>(names.edge_recon_condition)),
          dv_edge_recon_fallback_(this->particles_->template getVariableByName<Real>(names.edge_recon_fallback)),
          dv_signed_distance_((boundary_closure.ghost_edge_recon || boundary_closure.missing_moment_recon)
                                  ? this->particles_->template getVariableByName<Real>("SignedDistance")
                                  : nullptr),
          dv_normal_direction_((boundary_closure.ghost_edge_recon || boundary_closure.missing_moment_recon)
                                   ? this->particles_->template getVariableByName<Vecd>("NormalDirection")
                                   : nullptr)
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
              e_edge_recon_(encloser.dv_e_edge_recon_->DelegatedData(ex_policy)),
              j_edge_recon_(encloser.dv_j_edge_recon_->DelegatedData(ex_policy)),
              edge_recon_condition_(encloser.dv_edge_recon_condition_->DelegatedData(ex_policy)),
              edge_recon_fallback_(encloser.dv_edge_recon_fallback_->DelegatedData(ex_policy)),
              signed_distance_(encloser.dv_signed_distance_ != nullptr
                                   ? encloser.dv_signed_distance_->DelegatedData(ex_policy)
                                   : nullptr),
              normal_direction_(encloser.dv_normal_direction_ != nullptr
                                    ? encloser.dv_normal_direction_->DelegatedData(ex_policy)
                                    : nullptr),
              boundary_closure_(encloser.boundary_closure_), omega_(encloser.omega_),
              a_sign_(encloser.a_sign_), pair_weight_regularization_(encloser.pair_weight_regularization_),
              recon_condition_threshold_(encloser.recon_condition_threshold_),
              solver_local_rhs_scale_(encloser.solver_local_rhs_scale_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_[index_i];
            const Real sigma_i = sigma_[index_i];
            const Vecd a_i = active_a_[index_i];
            Eigen::Matrix<double, Dimensions, Dimensions> m_acc =
                Eigen::Matrix<double, Dimensions, Dimensions>::Zero();
            Eigen::Matrix<double, Dimensions, 1> b_acc = Eigen::Matrix<double, Dimensions, 1>::Zero();
            double max_neighbor_conductance = 0.0;
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
                max_neighbor_conductance = std::max(max_neighbor_conductance, static_cast<double>(std::abs(c_ij)));
                const Real edge_drop = computeOphelieEdgeFluxEdgeDrop(phi_i, phi_[index_j], a_i, active_a_[index_j],
                                                                      xj_minus_xi, omega_, a_sign_);
                const Eigen::Matrix<double, Dimensions, 1> e_hat =
                    r_ij_vec.template cast<double>() / static_cast<double>(distance);
                const double directional_e = static_cast<double>(edge_drop) / static_cast<double>(distance) *
                                             static_cast<double>(solver_local_rhs_scale_);
                const double w_ij = static_cast<double>(std::abs(c_ij));
                m_acc += w_ij * e_hat * e_hat.transpose();
                b_acc += w_ij * directional_e * e_hat;
            }
            const bool recon_boundary_active =
                (boundary_closure_.ghost_edge_recon || boundary_closure_.missing_moment_recon) &&
                signed_distance_ != nullptr && normal_direction_ != nullptr &&
                boundary_closure_.boundary_width_m > TinyReal && max_neighbor_conductance > TinyReal;
            if (recon_boundary_active)
            {
                Vecd unit_normal = Vecd::Zero();
                if (ophelieBoundaryShellParticle(signed_distance_[index_i], boundary_closure_.boundary_width_m) &&
                    ophelieTryUnitBoundaryNormal(normal_direction_, index_i, unit_normal))
                {
                    const Real ghost_conductance =
                        ophelieBoundaryGhostConductanceFromMaxNeighbor(static_cast<Real>(max_neighbor_conductance));
                    if (boundary_closure_.ghost_edge_recon)
                    {
                        accumulateOphelieNoFluxGhostEdgeLsContribution(m_acc, b_acc, unit_normal, ghost_conductance);
                    }
                    if (boundary_closure_.missing_moment_recon)
                    {
                        const Real moment_weight = ophelieBoundaryMissingMomentWeight(
                            ghost_conductance, boundary_closure_.boundary_width_m);
                        accumulateOphelieMissingMomentLsCorrection(m_acc, unit_normal, moment_weight);
                    }
                }
            }
            const double condition_proxy = m_acc.norm();
            edge_recon_condition_[index_i] = static_cast<Real>(condition_proxy);
            if (condition_proxy > static_cast<double>(recon_condition_threshold_))
            {
                const double tikhonov_eps = static_cast<double>(SqrtEps) * condition_proxy * condition_proxy;
                const Eigen::Matrix<double, Dimensions, 1> e_d = inverseTikhonov(m_acc, tikhonov_eps) * b_acc;
                Vecd e_i = e_d.cast<Real>();
                if (solver_local_rhs_scale_ < Real(1) - TinyReal)
                {
                    e_i /= solver_local_rhs_scale_;
                }
                if (!e_i.allFinite())
                {
                    e_i = Vecd::Zero();
                    e_edge_recon_[index_i] = e_i;
                    j_edge_recon_[index_i] = Vecd::Zero();
                    edge_recon_fallback_[index_i] = 1.0;
                }
                else
                {
                    e_edge_recon_[index_i] = e_i;
                    j_edge_recon_[index_i] = sigma_i * e_i;
                    edge_recon_fallback_[index_i] = 0.0;
                }
            }
            else
            {
                e_edge_recon_[index_i] = Vecd::Zero();
                j_edge_recon_[index_i] = Vecd::Zero();
                edge_recon_fallback_[index_i] = 1.0;
            }
        }

      protected:
        Real *Vol_;
        Real *phi_;
        Vecd *active_a_;
        Real *sigma_;
        Vecd *e_edge_recon_;
        Vecd *j_edge_recon_;
        Real *edge_recon_condition_;
        Real *edge_recon_fallback_;
        Real *signed_distance_;
        Vecd *normal_direction_;
        OphelieEdgeFluxBoundaryClosureFlags boundary_closure_;
        Real omega_;
        Real a_sign_;
        Real pair_weight_regularization_;
        Real recon_condition_threshold_;
        Real solver_local_rhs_scale_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    Real pair_weight_regularization_;
    Real a_sign_;
    Real recon_condition_threshold_;
    Real solver_local_rhs_scale_;
    OphelieEdgeFluxBoundaryClosureFlags boundary_closure_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Vecd> *dv_active_a_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Vecd> *dv_e_edge_recon_;
    DiscreteVariable<Vecd> *dv_j_edge_recon_;
    DiscreteVariable<Real> *dv_edge_recon_condition_;
    DiscreteVariable<Real> *dv_edge_recon_fallback_;
    DiscreteVariable<Real> *dv_signed_distance_;
    DiscreteVariable<Vecd> *dv_normal_direction_;
};

template <typename... Parameters>
class ReconstructOphelieEdgeFluxElectricCurrentCK<Inner<Parameters...>>
    : public ReconstructOphelieEdgeFluxElectricCurrentCK<Base, Inner<Parameters...>>
{
  public:
    ReconstructOphelieEdgeFluxElectricCurrentCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names,
                                                const OphelieEdgeFluxComponent &component, Real omega,
                                                Real pair_weight_regularization, Real recon_condition_threshold,
                                                Real solver_local_rhs_scale = Real(1),
                                                const OphelieEdgeFluxBoundaryClosureFlags &boundary_closure = {})
        : ReconstructOphelieEdgeFluxElectricCurrentCK<Base, Inner<Parameters...>>(
              inner_relation, names, component, omega, pair_weight_regularization, recon_condition_threshold,
              solver_local_rhs_scale, boundary_closure)
    {
    }
};

inline OphelieEdgeFluxResidualMetrics computeHostEdgeFluxResidualMetricsFromField(
    BaseParticles &particles, const std::string &edge_res_field, size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, edge_res_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *edge_res = particles.getVariableDataByName<Real>(edge_res_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OphelieEdgeFluxResidualMetrics metrics;
    Real sum_l2 = 0.0;
    Real total_vol = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real value = edge_res[i];
        const Real v = vol[i];
        sum_l2 += v * value * value;
        metrics.edge_res_linf = std::max(metrics.edge_res_linf, std::abs(value));
        metrics.edge_res_volume_integral += v * value;
        total_vol += v;
    }
    metrics.edge_res_l2 = std::sqrt(sum_l2 / (total_vol + TinyReal));
    metrics.global_conservation = metrics.edge_res_volume_integral;
    metrics.a_flux_sign = 1;
    return metrics;
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxResidualMetrics evaluateOphelieEdgeFluxResidualForComponent(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
    const OphelieEdgeFluxComponent &component, const OphelieParameters &params)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieEdgeFluxResidualCK<Inner<>>> compute_edge_res(
        inner, names, component, params.omega(), params.pair_weight_regularization_);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_edge_res.exec();
    return computeHostEdgeFluxResidualMetricsFromField(glass_body.getBaseParticles(), component.edge_residual_field,
                                                       glass_body.getBaseParticles().TotalRealParticles());
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxResidualMetrics evaluateOphelieEdgeFluxResidual(SolidBody &glass_body, Inner<> &inner,
                                                                      const OphelieGlassFieldNames &names,
                                                                      const OphelieParameters &params)
{
    return evaluateOphelieEdgeFluxResidualForComponent<ExecutionPolicy>(
        glass_body, inner, names, makeOphelieEdgeFluxImagComponent(names, params), params);
}

inline OphelieEdgeFluxPowerMetrics computeHostEdgeFluxPowerMetrics(BaseParticles &particles,
                                                                   const OphelieGlassFieldNames &names,
                                                                   size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, names.joule_heat_edge);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_imag);
    syncVariableToHost<Real>(particles, names.power_edge_particle);
    syncVariableToHost<Real>(particles, names.edge_recon_fallback);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Real *joule_edge = particles.getVariableDataByName<Real>(names.joule_heat_edge);
    const Real *joule_edge_recon = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_imag);
    const Real *power_particle = particles.getVariableDataByName<Real>(names.power_edge_particle);
    const Real *fallback = particles.getVariableDataByName<Real>(names.edge_recon_fallback);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OphelieEdgeFluxPowerMetrics metrics;
    Real joule_sum = 0.0;
    Real joule_recon_sum = 0.0;
    Real vol_sum = 0.0;
    size_t fallback_count = 0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        metrics.p_graph_edge += power_particle[i];
        metrics.joule_heat_edge_max = std::max(metrics.joule_heat_edge_max, joule_edge[i]);
        metrics.joule_heat_edge_recon_max = std::max(metrics.joule_heat_edge_recon_max, joule_edge_recon[i]);
        joule_sum += joule_edge[i] * vol[i];
        joule_recon_sum += joule_edge_recon[i] * vol[i];
        vol_sum += vol[i];
        if (fallback[i] > Real(0.5))
        {
            ++fallback_count;
        }
    }
    metrics.joule_heat_edge_mean = joule_sum / (vol_sum + TinyReal);
    metrics.joule_heat_edge_recon_mean = joule_recon_sum / (vol_sum + TinyReal);
    metrics.joule_heat_edge_max_over_mean =
        metrics.joule_heat_edge_max / (metrics.joule_heat_edge_mean + TinyReal);
    metrics.recon_fallback_fraction = static_cast<Real>(fallback_count) / static_cast<Real>(total_real_particles + 1);
    return metrics;
}

inline Real hostEdgeFluxReconPower(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                   size_t total_real_particles, const OphelieParameters &params);

/** Graph vs recon power breakdown (French gate expects p_graph/p_recon in [0.5, 2]). */
struct OphelieEdgeFluxPowerAuditDetail
{
    Real p_graph_sum = 0.0;
    Real p_joule_edge_vol = 0.0;
    Real p_recon_w = 0.0;
    Real p_graph_over_recon = 0.0;
    Real joule_recon_vol_w = 0.0;
};

inline OphelieEdgeFluxPowerAuditDetail hostOphelieEdgeFluxPowerAuditDetail(BaseParticles &particles,
                                                                         const OphelieGlassFieldNames &names,
                                                                         size_t n, const OphelieParameters &params)
{
    const OphelieEdgeFluxPowerMetrics metrics = computeHostEdgeFluxPowerMetrics(particles, names, n);
    OphelieEdgeFluxPowerAuditDetail detail;
    detail.p_graph_sum = metrics.p_graph_edge;
    syncVariableToHost<Real>(particles, names.joule_heat_edge);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *joule_edge = particles.getVariableDataByName<Real>(names.joule_heat_edge);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t i = 0; i < n; ++i)
    {
        detail.p_joule_edge_vol += joule_edge[i] * vol[i];
    }
    detail.p_recon_w = hostEdgeFluxReconPower(particles, names, n, params);
    detail.p_graph_over_recon = detail.p_graph_sum / (detail.p_recon_w + TinyReal);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_imag);
    const Real *joule_recon = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_imag);
    for (size_t i = 0; i < n; ++i)
    {
        detail.joule_recon_vol_w += joule_recon[i] * vol[i];
    }
    return detail;
}

inline void printOphelieEdgeFluxPowerAuditDetail(const OphelieEdgeFluxPowerAuditDetail &detail)
{
    std::cout << "[ophelie] edge-flux power audit (post-restore): p_graph_sum=" << detail.p_graph_sum
              << " p_joule_edge_vol=" << detail.p_joule_edge_vol << " p_recon_W=" << detail.p_recon_w
              << " joule_recon_vol_W=" << detail.joule_recon_vol_w
              << " p_graph/p_recon=" << detail.p_graph_over_recon
              << " (P_graph_edge_diagnostic vs P_recon_physical; graph is diagnostic-only)" << std::endl;
}

inline Real hostEdgeFluxReconPower(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                   size_t total_real_particles, const OphelieParameters &params)
{
    if (params.edge_flux_complex_)
    {
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
        syncVariableToHost<Real>(particles, "VolumetricMeasure");
        return hostVolWeightedSum(particles, names.joule_heat_edge_recon_complex, total_real_particles);
    }
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
    syncVariableToHost<Vecd>(particles, names.e_edge_recon_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Vecd *j_edge = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
    const Vecd *e_edge = particles.getVariableDataByName<Vecd>(names.e_edge_recon_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    Real power = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        power += Real(0.5) * (j_edge[i].dot(e_edge[i])) * vol[i];
    }
    return power;
}

class ComputeOphelieEdgeFluxComplexJouleHeatCK : public LocalDynamics
{
  public:
    ComputeOphelieEdgeFluxComplexJouleHeatCK(SPHBody &sph_body, const OphelieGlassFieldNames &names)
        : LocalDynamics(sph_body),
          dv_sigma_(particles_->template getVariableByName<Real>(names.sigma)),
          dv_e_real_(particles_->template getVariableByName<Vecd>(names.e_edge_recon_real)),
          dv_e_imag_(particles_->template getVariableByName<Vecd>(names.e_edge_recon_imag)),
          dv_q_complex_(particles_->template getVariableByName<Real>(names.joule_heat_edge_recon_complex))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)),
              e_real_(encloser.dv_e_real_->DelegatedData(ex_policy)),
              e_imag_(encloser.dv_e_imag_->DelegatedData(ex_policy)),
              q_complex_(encloser.dv_q_complex_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real sigma_i = sigma_[index_i];
            q_complex_[index_i] =
                Real(0.5) * sigma_i * (e_real_[index_i].squaredNorm() + e_imag_[index_i].squaredNorm());
        }

      protected:
        Real *sigma_;
        Vecd *e_real_;
        Vecd *e_imag_;
        Real *q_complex_;
    };

  protected:
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Vecd> *dv_e_real_;
    DiscreteVariable<Vecd> *dv_e_imag_;
    DiscreteVariable<Real> *dv_q_complex_;
};

inline OphelieEdgeFluxBoundaryClosureFlags resolveOphelieEdgeFluxBoundaryClosure(const OphelieParameters &params,
                                                                                 Real boundary_width_m,
                                                                                 BaseParticles &particles)
{
    OphelieEdgeFluxBoundaryClosureFlags flags =
        ophelieEdgeFluxBoundaryClosureFlagsFromMode(params.edge_recon_boundary_mode_, boundary_width_m);
    if (!ophelieEdgeFluxBoundaryClosureActive(flags))
    {
        return flags;
    }
    const bool have_normals = particles.getVariableDataByName<Vecd>("NormalDirection") != nullptr &&
                              particles.getVariableDataByName<Real>("SignedDistance") != nullptr;
    if (!have_normals)
    {
        std::cerr << "[ophelie] boundary closure requested (mode="
                  << ophelieEdgeReconBoundaryModeName(params.edge_recon_boundary_mode_)
                  << ") but NormalDirection/SignedDistance missing; falling back to standard operators" << std::endl;
        return {};
    }
    return flags;
}

inline void logOphelieEdgeFluxBoundaryClosure(const OphelieEdgeFluxBoundaryClosureFlags &flags,
                                            const std::string &chain_label)
{
    if (!ophelieEdgeFluxBoundaryClosureActive(flags))
    {
        return;
    }
    std::cout << "[ophelie] P5 boundary closure chain=" << chain_label << " width_m=" << flags.boundary_width_m
              << " ghost_edge=" << (flags.ghost_edge_recon ? 1 : 0)
              << " missing_moment=" << (flags.missing_moment_recon ? 1 : 0)
              << " phi_rhs_ghost=" << (flags.phi_rhs_ghost ? 1 : 0) << std::endl;
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxPowerMetrics execOphelieEdgeFluxPostPhiPipelineForComponent(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
    const OphelieEdgeFluxComponent &component, const OphelieParameters &params)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    const Real solver_local_rhs_scale = ophelieEdgeFluxEffectiveSolverLocalRhsScale(params);
    const Real boundary_width_m =
        params.edge_recon_boundary_width_factor_ *
        glass_body.getBaseParticles().getSPHAdaptation().ReferenceSmoothingLength();
    BaseParticles &particles_ref = glass_body.getBaseParticles();
    const OphelieEdgeFluxBoundaryClosureFlags boundary_closure =
        resolveOphelieEdgeFluxBoundaryClosure(params, boundary_width_m, particles_ref);
    logOphelieEdgeFluxBoundaryClosure(boundary_closure, component.chain_label);
    if (params.edge_flux_normalization_mode_ == OphelieEdgeFluxNormalizationMode::SolverLocal &&
        solver_local_rhs_scale < Real(1) - TinyReal)
    {
        std::cout << "[ophelie] edge recon solver-local rhs_scale=" << solver_local_rhs_scale << " chain="
                  << component.chain_label << std::endl;
    }
    InteractionDynamicsCK<ExecutionPolicy, ReconstructOphelieEdgeFluxElectricCurrentCK<Inner<>>> reconstruct_ej(
        inner, names, component, params.omega(), params.pair_weight_regularization_,
        params.edge_recon_condition_threshold_, solver_local_rhs_scale, boundary_closure);
    StateDynamics<ExecutionPolicy, ComputeOphelieEdgeFluxReconstructedJouleHeatCK> compute_joule_recon(
        glass_body, names.sigma, component.e_recon_field, component.joule_heat_recon_field);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    if (component.chain_label == std::string("imag"))
    {
        InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieEdgeFluxJouleHeatCK<Inner<>>> compute_joule(
            inner, names, params.omega(), params.pair_weight_regularization_, component.active_a_field,
            component.a_sign);
        compute_joule.exec();
    }
    reconstruct_ej.exec();
    compute_joule_recon.exec();
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    OphelieEdgeFluxPowerMetrics metrics = computeHostEdgeFluxPowerMetrics(particles, names, n);
    metrics.p_total_recon = hostEdgeFluxReconPower(particles, names, n, params);
    metrics.p_graph_over_recon = metrics.p_graph_edge / (metrics.p_total_recon + TinyReal);
    return metrics;
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxPowerMetrics execOphelieEdgeFluxPostPhiPipeline(SolidBody &glass_body, Inner<> &inner,
                                                                      const OphelieGlassFieldNames &names,
                                                                      const OphelieParameters &params,
                                                                      bool skip_real_phi_solve = false)
{
    const OphelieEdgeFluxComponent imag_component = makeOphelieEdgeFluxImagComponent(names, params);
    OphelieEdgeFluxPowerMetrics metrics =
        execOphelieEdgeFluxPostPhiPipelineForComponent<ExecutionPolicy>(glass_body, inner, names, imag_component, params);
    if (!params.edge_flux_complex_ || skip_real_phi_solve)
    {
        if (params.edge_flux_complex_ && skip_real_phi_solve)
        {
            StateDynamics<ExecutionPolicy, ComputeOphelieEdgeFluxComplexJouleHeatCK> compute_q_complex(glass_body,
                                                                                                         names);
            compute_q_complex.exec();
            BaseParticles &particles = glass_body.getBaseParticles();
            const size_t n = particles.TotalRealParticles();
            metrics.p_total_recon = hostEdgeFluxReconPower(particles, names, n, params);
            metrics.p_graph_over_recon = metrics.p_graph_edge / (metrics.p_total_recon + TinyReal);
        }
        return metrics;
    }
    const OphelieEdgeFluxComponent real_component = makeOphelieEdgeFluxRealComponent(names, params);
    (void)execOphelieEdgeFluxPostPhiPipelineForComponent<ExecutionPolicy>(glass_body, inner, names, real_component,
                                                                          params);
    StateDynamics<ExecutionPolicy, ComputeOphelieEdgeFluxComplexJouleHeatCK> compute_q_complex(glass_body, names);
    compute_q_complex.exec();
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    metrics.p_total_recon = hostEdgeFluxReconPower(particles, names, n, params);
    metrics.p_graph_over_recon = metrics.p_graph_edge / (metrics.p_total_recon + TinyReal);
    return metrics;
}

template <typename... RelationTypes>
class ComputeOphelieEdgeFluxEdgeDropPairStatsCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieEdgeFluxEdgeDropPairStatsCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieEdgeFluxEdgeDropPairStatsCK(RelationType<Parameters...> &relation,
                                              const OphelieGlassFieldNames &names, const std::string &phi_field,
                                              const std::string &active_a_field, Real omega,
                                              Real pair_weight_regularization, Real a_sign)
        : BaseInteraction(relation), omega_(omega), pair_weight_regularization_(pair_weight_regularization),
          a_sign_(a_sign), reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_phi_(this->particles_->template getVariableByName<Real>(phi_field)),
          dv_active_a_(this->particles_->template getVariableByName<Vecd>(active_a_field)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_edge_drop_abs_max_(this->particles_->template getVariableByName<Real>(names.edge_drop_abs_max)),
          dv_edge_drop_sq_mean_(this->particles_->template getVariableByName<Real>(names.edge_drop_sq_mean))
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
              edge_drop_abs_max_(encloser.dv_edge_drop_abs_max_->DelegatedData(ex_policy)),
              edge_drop_sq_mean_(encloser.dv_edge_drop_sq_mean_->DelegatedData(ex_policy)), omega_(encloser.omega_),
              a_sign_(encloser.a_sign_), pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_[index_i];
            const Real sigma_i = sigma_[index_i];
            const Vecd a_i = active_a_[index_i];
            Real abs_max = 0.0;
            Real sq_sum = 0.0;
            Real neighbor_count = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Vecd xj_minus_xi = -r_ij_vec;
                const Real edge_drop =
                    computeOphelieEdgeFluxEdgeDrop(phi_i, phi_[index_j], a_i, active_a_[index_j], xj_minus_xi, omega_,
                                                   a_sign_);
                abs_max = std::max(abs_max, std::abs(edge_drop));
                sq_sum += edge_drop * edge_drop;
                neighbor_count += 1.0;
                (void)sigma_i;
                (void)pair_weight_regularization_;
                (void)reference_smoothing_length_;
            }
            edge_drop_abs_max_[index_i] = abs_max;
            edge_drop_sq_mean_[index_i] = sq_sum / (neighbor_count + TinyReal);
        }

      protected:
        Real *Vol_;
        Real *phi_;
        Vecd *active_a_;
        Real *sigma_;
        Real *edge_drop_abs_max_;
        Real *edge_drop_sq_mean_;
        Real omega_;
        Real a_sign_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    Real pair_weight_regularization_;
    Real a_sign_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Vecd> *dv_active_a_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_edge_drop_abs_max_;
    DiscreteVariable<Real> *dv_edge_drop_sq_mean_;
};

template <typename... Parameters>
class ComputeOphelieEdgeFluxEdgeDropPairStatsCK<Inner<Parameters...>>
    : public ComputeOphelieEdgeFluxEdgeDropPairStatsCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieEdgeFluxEdgeDropPairStatsCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names,
                                              const OphelieEdgeFluxEdgeDropComponent &component, Real omega,
                                              Real pair_weight_regularization)
        : ComputeOphelieEdgeFluxEdgeDropPairStatsCK<Base, Inner<Parameters...>>(
              inner_relation, names, component.phi_field, component.active_a_field, omega, pair_weight_regularization,
              component.a_sign)
    {
    }
};

inline OphelieEdgeFluxEdgeDropMetrics computeHostEdgeFluxEdgeDropMetrics(BaseParticles &particles,
                                                                         const OphelieGlassFieldNames &names,
                                                                         size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, names.edge_drop_abs_max);
    syncVariableToHost<Real>(particles, names.edge_drop_sq_mean);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    const Real *abs_max = particles.getVariableDataByName<Real>(names.edge_drop_abs_max);
    const Real *sq_mean = particles.getVariableDataByName<Real>(names.edge_drop_sq_mean);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");

    OphelieEdgeFluxEdgeDropMetrics metrics;
    Real sum_l2 = 0.0;
    Real sum_abs = 0.0;
    Real total_vol = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        metrics.edge_drop_linf = std::max(metrics.edge_drop_linf, abs_max[i]);
        sum_l2 += vol[i] * sq_mean[i];
        sum_abs += vol[i] * std::sqrt(sq_mean[i]);
        total_vol += vol[i];
        ++metrics.pair_count;
    }
    metrics.edge_drop_l2 = std::sqrt(sum_l2 / (total_vol + TinyReal));
    metrics.edge_drop_mean_abs = sum_abs / (total_vol + TinyReal);
    return metrics;
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxEdgeDropMetrics evaluateOphelieEdgeFluxEdgeDropMetricsForComponent(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    const OphelieEdgeFluxEdgeDropComponent &component)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieEdgeFluxEdgeDropPairStatsCK<Inner<>>> compute_stats(
        inner, names, component, params.omega(), params.pair_weight_regularization_);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_stats.exec();
    return computeHostEdgeFluxEdgeDropMetrics(glass_body.getBaseParticles(), names,
                                              glass_body.getBaseParticles().TotalRealParticles());
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxEdgeDropMetrics evaluateOphelieEdgeFluxImagEdgeDropMetrics(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params)
{
    return evaluateOphelieEdgeFluxEdgeDropMetricsForComponent<ExecutionPolicy>(
        glass_body, inner, names, params, makeOphelieEdgeFluxImagEdgeDropComponent(names, params));
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxEdgeDropMetrics evaluateOphelieEdgeFluxRealEdgeDropMetrics(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params)
{
    return evaluateOphelieEdgeFluxEdgeDropMetricsForComponent<ExecutionPolicy>(
        glass_body, inner, names, params, makeOphelieEdgeFluxRealEdgeDropComponent(names, params));
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxEdgeDropMetrics evaluateOphelieEdgeFluxEdgeDropMetrics(SolidBody &glass_body, Inner<> &inner,
                                                                             const OphelieGlassFieldNames &names,
                                                                             const OphelieParameters &params)
{
    return evaluateOphelieEdgeFluxImagEdgeDropMetrics<ExecutionPolicy>(glass_body, inner, names, params);
}

struct OphelieEdgeFluxQAntisymMetrics
{
    Real q_antisym_l1 = 0.0;
    Real q_antisym_l2 = 0.0;
    Real q_antisym_linf = 0.0;
    Real q_antisym_rel_l2 = 0.0;
    size_t q_nonfinite_count = 0;
    size_t pair_count = 0;
};

struct OphelieEdgeFluxQSpatialMetrics
{
    Real q_min = 0.0;
    Real q_max = 0.0;
    Real q_mean = 0.0;
    Real q_max_over_mean = 0.0;
    size_t q_nonfinite_count = 0;
    size_t q_negative_count = 0;
    Real q_outer_mean = 0.0;
    Real q_center_mean = 0.0;
    Real q_outer_over_center = 0.0;
    bool soft_gate_passed = false;
};

template <typename... RelationTypes>
class ComputeOphelieEdgeFluxQAntisymmetryCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieEdgeFluxQAntisymmetryCK<Base, RelationType<Parameters...>> : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieEdgeFluxQAntisymmetryCK(RelationType<Parameters...> &relation, const OphelieGlassFieldNames &names,
                                          Real omega, Real pair_weight_regularization, const std::string &phi_field,
                                          const std::string &active_a_field, Real a_sign)
        : BaseInteraction(relation), omega_(omega), pair_weight_regularization_(pair_weight_regularization),
          a_sign_(a_sign), reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_phi_(this->particles_->template getVariableByName<Real>(phi_field)),
          dv_active_a_(this->particles_->template getVariableByName<Vecd>(active_a_field)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_q_antisym_max_(this->particles_->template getVariableByName<Real>(names.edge_q_antisym_max)),
          dv_q_antisym_sq_sum_(this->particles_->template getVariableByName<Real>(names.edge_q_antisym_sq_sum)),
          dv_q_scale_sq_sum_(this->particles_->template getVariableByName<Real>(names.edge_q_scale_sq_sum)),
          dv_q_neighbor_count_(this->particles_->template getVariableByName<Real>(names.edge_q_neighbor_count)),
          dv_q_nonfinite_count_(this->particles_->template getVariableByName<Real>(names.edge_q_nonfinite_count))
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
              q_antisym_max_(encloser.dv_q_antisym_max_->DelegatedData(ex_policy)),
              q_antisym_sq_sum_(encloser.dv_q_antisym_sq_sum_->DelegatedData(ex_policy)),
              q_scale_sq_sum_(encloser.dv_q_scale_sq_sum_->DelegatedData(ex_policy)),
              q_neighbor_count_(encloser.dv_q_neighbor_count_->DelegatedData(ex_policy)),
              q_nonfinite_count_(encloser.dv_q_nonfinite_count_->DelegatedData(ex_policy)), omega_(encloser.omega_),
              a_sign_(encloser.a_sign_), pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_[index_i];
            const Real sigma_i = sigma_[index_i];
            const Vecd a_i = active_a_[index_i];
            Real antisym_max = 0.0;
            Real antisym_sq_sum = 0.0;
            Real q_scale_sq_sum = 0.0;
            Real neighbor_count = 0.0;
            Real nonfinite_count = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Vecd xj_minus_xi = -r_ij_vec;
                const Real distance = r_ij_vec.norm();
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real dW_jiV_i = this->dW_ij(index_j, index_i) * Vol_[index_i];
                const Real c_ij = computeOphelieEdgeFluxPairConductance(
                    sigma_i, sigma_[index_j], dW_ijV_j, distance, distance_sq, pair_weight_regularization_,
                    reference_smoothing_length_);
                const Real c_ji = computeOphelieEdgeFluxPairConductance(
                    sigma_[index_j], sigma_i, dW_jiV_i, distance, distance_sq, pair_weight_regularization_,
                    reference_smoothing_length_);
                const Real edge_drop_ij = computeOphelieEdgeFluxEdgeDrop(phi_i, phi_[index_j], a_i, active_a_[index_j],
                                                                         xj_minus_xi, omega_, a_sign_);
                const Real ji_a_sign = a_sign_ < Real(0) ? a_sign_ : Real(1.0);
                const Real edge_drop_ji = computeOphelieEdgeFluxEdgeDrop(phi_[index_j], phi_i, active_a_[index_j], a_i,
                                                                       r_ij_vec, omega_, ji_a_sign);
                const Real q_ij = -c_ij * edge_drop_ij;
                const Real q_ji = -c_ji * edge_drop_ji;
                if (!std::isfinite(q_ij) || !std::isfinite(q_ji))
                {
                    nonfinite_count += 1.0;
                }
                const Real antisym = q_ij + q_ji;
                antisym_max = std::max(antisym_max, std::abs(antisym));
                antisym_sq_sum += antisym * antisym;
                q_scale_sq_sum += q_ij * q_ij;
                neighbor_count += 1.0;
            }
            q_antisym_max_[index_i] = antisym_max;
            q_antisym_sq_sum_[index_i] = antisym_sq_sum;
            q_scale_sq_sum_[index_i] = q_scale_sq_sum;
            q_neighbor_count_[index_i] = neighbor_count;
            q_nonfinite_count_[index_i] = nonfinite_count;
        }

      protected:
        Real *Vol_;
        Real *phi_;
        Vecd *active_a_;
        Real *sigma_;
        Real *q_antisym_max_;
        Real *q_antisym_sq_sum_;
        Real *q_scale_sq_sum_;
        Real *q_neighbor_count_;
        Real *q_nonfinite_count_;
        Real omega_;
        Real a_sign_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
    };

  protected:
    Real omega_;
    Real pair_weight_regularization_;
    Real a_sign_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_phi_;
    DiscreteVariable<Vecd> *dv_active_a_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_q_antisym_max_;
    DiscreteVariable<Real> *dv_q_antisym_sq_sum_;
    DiscreteVariable<Real> *dv_q_scale_sq_sum_;
    DiscreteVariable<Real> *dv_q_neighbor_count_;
    DiscreteVariable<Real> *dv_q_nonfinite_count_;
};

template <typename... Parameters>
class ComputeOphelieEdgeFluxQAntisymmetryCK<Inner<Parameters...>>
    : public ComputeOphelieEdgeFluxQAntisymmetryCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieEdgeFluxQAntisymmetryCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names,
                                          Real omega, Real pair_weight_regularization, const std::string &phi_field,
                                          const std::string &active_a_field, Real a_sign)
        : ComputeOphelieEdgeFluxQAntisymmetryCK<Base, Inner<Parameters...>>(
              inner_relation, names, omega, pair_weight_regularization, phi_field, active_a_field, a_sign)
    {
    }
};

inline OphelieEdgeFluxQAntisymMetrics computeHostEdgeFluxQAntisymMetrics(BaseParticles &particles,
                                                                         const OphelieGlassFieldNames &names,
                                                                         size_t total_real_particles)
{
    syncVariableToHost<Real>(particles, names.edge_q_antisym_max);
    syncVariableToHost<Real>(particles, names.edge_q_antisym_sq_sum);
    syncVariableToHost<Real>(particles, names.edge_q_scale_sq_sum);
    syncVariableToHost<Real>(particles, names.edge_q_neighbor_count);
    syncVariableToHost<Real>(particles, names.edge_q_nonfinite_count);

    const Real *antisym_max = particles.getVariableDataByName<Real>(names.edge_q_antisym_max);
    const Real *antisym_sq_sum = particles.getVariableDataByName<Real>(names.edge_q_antisym_sq_sum);
    const Real *q_scale_sq_sum = particles.getVariableDataByName<Real>(names.edge_q_scale_sq_sum);
    const Real *neighbor_count = particles.getVariableDataByName<Real>(names.edge_q_neighbor_count);
    const Real *nonfinite_count = particles.getVariableDataByName<Real>(names.edge_q_nonfinite_count);

    OphelieEdgeFluxQAntisymMetrics metrics;
    Real total_antisym_abs = 0.0;
    Real total_antisym_sq = 0.0;
    Real total_q_scale_sq = 0.0;
    Real total_pairs = 0.0;
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real n_neighbors = neighbor_count[i];
        if (n_neighbors <= TinyReal)
        {
            continue;
        }
        metrics.q_antisym_linf = std::max(metrics.q_antisym_linf, antisym_max[i]);
        total_antisym_abs += std::sqrt(antisym_sq_sum[i] / n_neighbors) * n_neighbors;
        total_antisym_sq += antisym_sq_sum[i];
        total_q_scale_sq += q_scale_sq_sum[i];
        total_pairs += n_neighbors;
        metrics.q_nonfinite_count += static_cast<size_t>(nonfinite_count[i] + Real(0.5));
        metrics.pair_count += static_cast<size_t>(n_neighbors + Real(0.5));
    }
    metrics.q_antisym_l1 = total_antisym_abs / (total_pairs + TinyReal);
    metrics.q_antisym_l2 = std::sqrt(total_antisym_sq / (total_pairs + TinyReal));
    metrics.q_antisym_rel_l2 = metrics.q_antisym_l2 / (std::sqrt(total_q_scale_sq / (total_pairs + TinyReal)) + TinyReal);
    return metrics;
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxQAntisymMetrics evaluateOphelieEdgeFluxQAntisymmetryForComponent(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names,
    const OphelieEdgeFluxComponent &component, const OphelieParameters &params)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieEdgeFluxQAntisymmetryCK<Inner<>>> compute_q_antisym(
        inner, names, params.omega(), params.pair_weight_regularization_, component.phi_field,
        component.active_a_field, component.a_sign);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_q_antisym.exec();
    return computeHostEdgeFluxQAntisymMetrics(glass_body.getBaseParticles(), names,
                                              glass_body.getBaseParticles().TotalRealParticles());
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxQAntisymMetrics evaluateOphelieEdgeFluxQAntisymmetry(SolidBody &glass_body, Inner<> &inner,
                                                                           const OphelieGlassFieldNames &names,
                                                                           const OphelieParameters &params)
{
    return evaluateOphelieEdgeFluxQAntisymmetryForComponent<ExecutionPolicy>(
        glass_body, inner, names, makeOphelieEdgeFluxImagComponent(names, params), params);
}

inline OphelieEdgeFluxQSpatialMetrics computeHostEdgeFluxQSpatialMetrics(
    BaseParticles &particles, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    size_t total_real_particles, const Vecd &center, Real outer_radius, Real center_radius,
    Real q_max_over_mean_max, Real outer_over_center_min)
{
    const std::string &q_field = params.edge_flux_complex_ ? names.joule_heat_edge_recon_complex
                                                           : names.joule_heat_edge_recon_imag;
    syncVariableToHost<Real>(particles, q_field);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, "Position");
    const Real *q = particles.getVariableDataByName<Real>(q_field);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");

    OphelieEdgeFluxQSpatialMetrics metrics;
    Real q_sum = 0.0;
    Real vol_sum = 0.0;
    Real outer_sum = 0.0;
    Real outer_vol = 0.0;
    Real center_sum = 0.0;
    Real center_vol = 0.0;
    metrics.q_min = q[0];
    metrics.q_max = q[0];
    for (size_t i = 0; i < total_real_particles; ++i)
    {
        const Real qi = q[i];
        if (!std::isfinite(qi))
        {
            ++metrics.q_nonfinite_count;
            continue;
        }
        if (qi < 0.0)
        {
            ++metrics.q_negative_count;
        }
        metrics.q_min = std::min(metrics.q_min, qi);
        metrics.q_max = std::max(metrics.q_max, qi);
        q_sum += qi * vol[i];
        vol_sum += vol[i];
        const Real r = (pos[i] - center).head<2>().norm();
        if (r >= outer_radius)
        {
            outer_sum += qi * vol[i];
            outer_vol += vol[i];
        }
        if (r <= center_radius)
        {
            center_sum += qi * vol[i];
            center_vol += vol[i];
        }
    }
    metrics.q_mean = q_sum / (vol_sum + TinyReal);
    metrics.q_max_over_mean = metrics.q_max / (metrics.q_mean + TinyReal);
    metrics.q_outer_mean = outer_sum / (outer_vol + TinyReal);
    metrics.q_center_mean = center_sum / (center_vol + TinyReal);
    metrics.q_outer_over_center = metrics.q_outer_mean / (metrics.q_center_mean + TinyReal);
    metrics.soft_gate_passed = metrics.q_nonfinite_count == 0 && metrics.q_negative_count == 0 &&
                               metrics.q_max > metrics.q_mean && metrics.q_mean > 0.0 &&
                               metrics.q_max_over_mean > 1.0 && metrics.q_max_over_mean < q_max_over_mean_max &&
                               metrics.q_outer_over_center > outer_over_center_min;
    return metrics;
}

inline void logOphelieEdgeFluxQAntisymMetrics(const OphelieEdgeFluxQAntisymMetrics &metrics,
                                              const char *chain_label = nullptr)
{
    std::cout << "[ophelie] edge_flux_q_antisym";
    if (chain_label != nullptr)
    {
        std::cout << " chain=" << chain_label;
    }
    std::cout << ": q_antisym_l1=" << metrics.q_antisym_l1
              << " q_antisym_l2=" << metrics.q_antisym_l2 << " q_antisym_linf=" << metrics.q_antisym_linf
              << " q_antisym_rel_l2=" << metrics.q_antisym_rel_l2 << " q_nonfinite_count=" << metrics.q_nonfinite_count
              << " pair_count=" << metrics.pair_count << std::endl;
}

inline void logOphelieEdgeFluxQSpatialMetrics(const OphelieEdgeFluxQSpatialMetrics &metrics)
{
    std::cout << "[ophelie] edge_flux_q_spatial: Q_min=" << metrics.q_min << " Q_max=" << metrics.q_max
              << " Q_mean=" << metrics.q_mean << " Q_max/mean=" << metrics.q_max_over_mean
              << " Q_nonfinite=" << metrics.q_nonfinite_count << " Q_negative=" << metrics.q_negative_count
              << " Q_outer_mean=" << metrics.q_outer_mean << " Q_center_mean=" << metrics.q_center_mean
              << " Q_outer/center=" << metrics.q_outer_over_center
              << " soft_gate=" << (metrics.soft_gate_passed ? 1 : 0) << std::endl;
}

inline void logOphelieEdgeFluxPowerMetrics(const OphelieEdgeFluxPowerMetrics &metrics)
{
    std::cout << "[ophelie] edge_flux_power: P_graph_edge=" << metrics.p_graph_edge
              << " P_total_recon=" << metrics.p_total_recon << " P_graph_over_recon=" << metrics.p_graph_over_recon
              << " Q_edge_graph_max=" << metrics.joule_heat_edge_max
              << " Q_edge_recon_max=" << metrics.joule_heat_edge_recon_max
              << " Q_edge_graph_mean=" << metrics.joule_heat_edge_mean
              << " Q_edge_recon_mean=" << metrics.joule_heat_edge_recon_mean
              << " Q_edge_graph_max/mean=" << metrics.joule_heat_edge_max_over_mean
              << " recon_fallback_frac=" << metrics.recon_fallback_fraction << std::endl;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_H
