#ifndef ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_DIAGNOSTICS_H
#define ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_DIAGNOSTICS_H

#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_cli.h"
#include "electromagnetic_ophelie_parameters.h"
#include "interaction_ck.h"
#include "update_body_relation.h"

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

template <typename... RelationTypes>
class ComputeOphelieLegacyEdgeFluxResidualCK;

template <template <typename...> class RelationType, typename... Parameters>
class ComputeOphelieLegacyEdgeFluxResidualCK<Base, RelationType<Parameters...>>
    : public Interaction<RelationType<Parameters...>>
{
    using BaseInteraction = Interaction<RelationType<Parameters...>>;

  public:
    ComputeOphelieLegacyEdgeFluxResidualCK(RelationType<Parameters...> &relation, const OphelieGlassFieldNames &names,
                                           Real omega, Real pair_weight_regularization, int a_flux_sign)
        : BaseInteraction(relation), omega_(omega), pair_weight_regularization_(pair_weight_regularization),
          a_flux_sign_(a_flux_sign),
          reference_smoothing_length_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          dv_Vol_(this->particles_->template getVariableByName<Real>("VolumetricMeasure")),
          dv_phi_imag_(this->particles_->template getVariableByName<Real>(names.phi_imag)),
          dv_a_src_real_(this->particles_->template getVariableByName<Vecd>(names.a_src_real)),
          dv_sigma_(this->particles_->template getVariableByName<Real>(names.sigma)),
          dv_edge_res_(this->particles_->template getVariableByName<Real>(names.div_j_imag))
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
              edge_res_(encloser.dv_edge_res_->DelegatedData(ex_policy)), omega_(encloser.omega_),
              pair_weight_regularization_(encloser.pair_weight_regularization_),
              reference_smoothing_length_(encloser.reference_smoothing_length_), a_flux_sign_(encloser.a_flux_sign_)
        {
        }

        void interact(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Real phi_i = phi_imag_[index_i];
            const Real sigma_i = sigma_[index_i];
            const Vecd a_i = a_src_real_[index_i];
            Real laplace_i = 0.0;
            Real a_rhs_i = 0.0;
            for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
            {
                const UnsignedInt index_j = this->neighbor_index_[n];
                const Vecd r_ij_vec = this->vec_r_ij(index_i, index_j);
                const Real distance = r_ij_vec.norm();
                const Real distance_sq = r_ij_vec.squaredNorm();
                const Real dW_ijV_j = this->dW_ij(index_i, index_j) * Vol_[index_j];
                const Real sigma_ij = harmonicMean(sigma_i, sigma_[index_j]);
                const Real pair_weight = sigma_ij * pairwiseNegativeLaplaceWeight(
                                                              dW_ijV_j, distance, distance_sq,
                                                              pair_weight_regularization_, reference_smoothing_length_);
                laplace_i += pair_weight * (phi_i - phi_imag_[index_j]);
                const Vecd g_ij = pairwiseGradientWeightUncorrected(dW_ijV_j, this->e_ij(index_i, index_j));
                a_rhs_i -= omega_ * sigma_ij * g_ij.dot(a_i - a_src_real_[index_j]);
            }
            edge_res_[index_i] = laplace_i - static_cast<Real>(a_flux_sign_) * a_rhs_i;
        }

      protected:
        Real *Vol_;
        Real *phi_imag_;
        Vecd *a_src_real_;
        Real *sigma_;
        Real *edge_res_;
        Real omega_;
        Real pair_weight_regularization_;
        Real reference_smoothing_length_;
        int a_flux_sign_;
    };

  protected:
    Real omega_;
    Real pair_weight_regularization_;
    int a_flux_sign_;
    Real reference_smoothing_length_;
    DiscreteVariable<Real> *dv_Vol_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_a_src_real_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Real> *dv_edge_res_;
};

template <typename... Parameters>
class ComputeOphelieLegacyEdgeFluxResidualCK<Inner<Parameters...>>
    : public ComputeOphelieLegacyEdgeFluxResidualCK<Base, Inner<Parameters...>>
{
  public:
    ComputeOphelieLegacyEdgeFluxResidualCK(Inner<Parameters...> &inner_relation, const OphelieGlassFieldNames &names,
                                           Real omega, Real pair_weight_regularization, int a_flux_sign)
        : ComputeOphelieLegacyEdgeFluxResidualCK<Base, Inner<Parameters...>>(inner_relation, names, omega,
                                                                             pair_weight_regularization, a_flux_sign)
    {
    }
};

template <class ExecutionPolicy>
inline OphelieEdgeFluxResidualMetrics evaluateOphelieLegacyEdgeFluxResidual(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params,
    int a_flux_sign)
{
    UpdateCellLinkedList<ExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<ExecutionPolicy, Inner<>> update_inner_relation(inner);
    InteractionDynamicsCK<ExecutionPolicy, ComputeOphelieLegacyEdgeFluxResidualCK<Inner<>>> compute_edge_res(
        inner, names, params.omega(), params.pair_weight_regularization_, a_flux_sign);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_edge_res.exec();
    OphelieEdgeFluxResidualMetrics metrics = computeHostEdgeFluxResidualMetricsFromField(
        glass_body.getBaseParticles(), names.div_j_imag, glass_body.getBaseParticles().TotalRealParticles());
    metrics.a_flux_sign = a_flux_sign;
    return metrics;
}

template <class ExecutionPolicy>
inline OphelieEdgeFluxResidualMetrics evaluateOphelieEdgeFluxResidualBestSignAtCurrentPhi(
    SolidBody &glass_body, Inner<> &inner, const OphelieGlassFieldNames &names, const OphelieParameters &params)
{
    if (ophelieUseEdgeFluxElectromotiveRhs(params))
    {
        return evaluateOphelieEdgeFluxResidual<ExecutionPolicy>(glass_body, inner, names, params);
    }
    const OphelieEdgeFluxResidualMetrics plus =
        evaluateOphelieLegacyEdgeFluxResidual<ExecutionPolicy>(glass_body, inner, names, params, 1);
    const OphelieEdgeFluxResidualMetrics minus =
        evaluateOphelieLegacyEdgeFluxResidual<ExecutionPolicy>(glass_body, inner, names, params, -1);
    return plus.edge_res_l2 <= minus.edge_res_l2 ? plus : minus;
}

inline void logOphelieEdgeFluxResidualMetrics(const char *stage, const OphelieEdgeFluxResidualMetrics &metrics)
{
    std::cout << "[ophelie] edge_flux_diag stage=" << stage << " a_flux_sign=" << metrics.a_flux_sign
              << " edge_res_l2=" << metrics.edge_res_l2 << " edge_res_linf=" << metrics.edge_res_linf
              << " edge_res_volume_integral=" << metrics.edge_res_volume_integral
              << " global_conservation=" << metrics.global_conservation << std::endl;
}

inline void logOphelieEdgeFluxDiagnosticReport(const OphelieEdgeFluxDiagnosticReport &report)
{
    logOphelieEdgeFluxResidualMetrics("level0", report.level0);
    logOphelieEdgeFluxResidualMetrics("post_phi", report.post_phi);
    std::cout << "[ophelie] edge_flux_diag edge_res_red_l2=" << report.edge_res_red_l2 << std::endl;
}

inline void appendOphelieEdgeFluxDiagnosticCsv(const std::string &path, const std::string &case_label, Real dp,
                                               const OphelieParameters &params,
                                               const OphelieEdgeFluxDiagnosticReport &report)
{
    namespace fs = std::filesystem;
    const bool write_header = path.empty() || !fs::exists(path) || fs::file_size(path) == 0;
    std::ofstream csv(path, std::ios::app);
    if (!csv)
    {
        return;
    }
    if (write_header)
    {
        csv << "case_label,dp,projection_operator,current_form,edge_res_l2_level0,edge_res_linf_level0,"
               "edge_res_integral_level0,edge_res_l2_post_phi,edge_res_linf_post_phi,edge_res_integral_post_phi,"
               "edge_res_red_l2,a_flux_sign_level0,a_flux_sign_post_phi\n";
    }
    csv << case_label << "," << dp << "," << phiProjectionOperatorKindName(inferOpheliePhiProjectionOperatorKind(params))
        << "," << ophelieCurrentFormKindName(params.ophelie_current_form_) << "," << report.level0.edge_res_l2 << ","
        << report.level0.edge_res_linf << "," << report.level0.edge_res_volume_integral << ","
        << report.post_phi.edge_res_l2 << "," << report.post_phi.edge_res_linf << ","
        << report.post_phi.edge_res_volume_integral << "," << report.edge_res_red_l2 << ","
        << report.level0.a_flux_sign << "," << report.post_phi.a_flux_sign << "\n";
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_EDGE_FLUX_DIAGNOSTICS_H
