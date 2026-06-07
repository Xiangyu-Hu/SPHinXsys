#ifndef APHI_ASSEMBLE_LHS_DEBUG_CK_HPP
#define APHI_ASSEMBLE_LHS_DEBUG_CK_HPP

#include "electromagnetic_dynamics/diagnostics/aphi_assemble_lhs_debug_ck.h"
#include "electromagnetic_dynamics/aphi_block_zero_ck.hpp"
#include "electromagnetic_dynamics/aphi_div_sigma_a_coupling_ck.hpp"
#include "electromagnetic_dynamics/aphi_grad_phi_coupling_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_gradient_divergence_debug_ck.hpp"
#include "electromagnetic_dynamics/aphi_laplace_ck.hpp"
#include "electromagnetic_dynamics/aphi_phi_gauge_penalty_ck.hpp"
#include "electromagnetic_dynamics/aphi_reaction_ck.hpp"

namespace SPH
{
namespace electromagnetics
{

template <class ExecutionPolicy>
inline AphiAssembleLhsDebugDynamicsBundle<ExecutionPolicy>::AphiAssembleLhsDebugDynamicsBundle(
    SPHBody &sph_body, Inner<> &inner_relation, const AphiVariableNames &variable_names,
    const AphiLhsAssemblyOptions &options)
    : options_(options), names_(variable_names),
      zero_lhs_(sph_body, variable_names.lhs),
      reaction_(sph_body, options.omega, variable_names),
      laplace_phi_real_(DynamicsArgs(inner_relation, variable_names.solution.phi_real, variable_names.material.sigma,
                                     variable_names.lhs.phi_real)),
      laplace_phi_imag_(DynamicsArgs(inner_relation, variable_names.solution.phi_imag, variable_names.material.sigma,
                                     variable_names.lhs.phi_imag)),
      laplace_a_real_(DynamicsArgs(inner_relation, variable_names.solution.a_real, variable_names.material.nu,
                                   variable_names.lhs.a_real)),
      laplace_a_imag_(DynamicsArgs(inner_relation, variable_names.solution.a_imag, variable_names.material.nu,
                                   variable_names.lhs.a_imag)),
      grad_phi_pairwise_(DynamicsArgs(inner_relation, variable_names)),
      div_sigma_a_pairwise_(DynamicsArgs(inner_relation, options.omega, variable_names)),
      linear_correction_matrix_(DynamicsArgs(inner_relation, 0.0)),
      phi_real_gradient_(DynamicsArgs(inner_relation, variable_names.solution.phi_real)),
      phi_imag_gradient_(DynamicsArgs(inner_relation, variable_names.solution.phi_imag)),
      a_real_gradient_(DynamicsArgs(inner_relation, variable_names.solution.a_real)),
      a_imag_gradient_(DynamicsArgs(inner_relation, variable_names.solution.a_imag)),
      div_a_real_diagnostic_(sph_body, variable_names.solution.a_real + "Gradient",
                             aphiDefaultADivergencePenaltyScratchNames().div_a_real),
      div_a_imag_diagnostic_(sph_body, variable_names.solution.a_imag + "Gradient",
                             aphiDefaultADivergencePenaltyScratchNames().div_a_imag),
      grad_phi_diagnostic_(sph_body, variable_names),
      div_sigma_a_diagnostic_(sph_body, options.omega, variable_names, variable_names.diagnostic),
      phi_gauge_penalty_(sph_body, variable_names.solution, variable_names.lhs, options.phi_gauge_penalty),
      div_a_real_gradient_(DynamicsArgs(inner_relation, aphiDefaultADivergencePenaltyScratchNames().div_a_real)),
      div_a_imag_gradient_(DynamicsArgs(inner_relation, aphiDefaultADivergencePenaltyScratchNames().div_a_imag)),
      div_a_real_pairwise_(DynamicsArgs(inner_relation, variable_names.solution.a_real,
                                        aphiDefaultADivergencePenaltyScratchNames().div_a_real)),
      div_a_imag_pairwise_(DynamicsArgs(inner_relation, variable_names.solution.a_imag,
                                        aphiDefaultADivergencePenaltyScratchNames().div_a_imag)),
      div_a_real_gradient_pairwise_(DynamicsArgs(inner_relation, aphiDefaultADivergencePenaltyScratchNames().div_a_real,
                                                 aphiDefaultADivergencePenaltyScratchNames().grad_div_a_real)),
      div_a_imag_gradient_pairwise_(DynamicsArgs(inner_relation, aphiDefaultADivergencePenaltyScratchNames().div_a_imag,
                                                 aphiDefaultADivergencePenaltyScratchNames().grad_div_a_imag)),
      a_divergence_penalty_(sph_body, aphiDefaultADivergencePenaltyScratchNames().grad_div_a_real,
                            aphiDefaultADivergencePenaltyScratchNames().grad_div_a_imag, variable_names.lhs,
                            options.a_divergence_penalty)
{
}

template <class ExecutionPolicy>
inline void AphiAssembleLhsDebugDynamicsBundle<ExecutionPolicy>::exec()
{
    zero_lhs_.exec();

    if (options_.terms.laplace_phi)
    {
        laplace_phi_real_.exec();
        laplace_phi_imag_.exec();
    }

    if (options_.terms.laplace_a)
    {
        laplace_a_real_.exec();
        laplace_a_imag_.exec();
    }

    if (options_.terms.reaction)
    {
        reaction_.exec();
    }

    const bool need_b_correction =
        (options_.terms.grad_phi_coupling &&
         options_.grad_phi_mode == AphiGradPhiCouplingMode::LinearGradientWithBCorrection) ||
        (options_.terms.div_sigma_a_coupling &&
         options_.div_mode == AphiDivCouplingMode::LinearGradientTraceWithBCorrection);

    if (need_b_correction)
    {
        linear_correction_matrix_.exec();
    }

    if (options_.terms.grad_phi_coupling)
    {
        if (options_.grad_phi_mode == AphiGradPhiCouplingMode::PairwiseUncorrected)
        {
            grad_phi_pairwise_.exec();
        }
        else
        {
            phi_real_gradient_.exec();
            phi_imag_gradient_.exec();
            grad_phi_diagnostic_.exec();
        }
    }

    if (options_.terms.div_sigma_a_coupling)
    {
        if (options_.div_mode == AphiDivCouplingMode::PairwiseUncorrected)
        {
            div_sigma_a_pairwise_.exec();
        }
        else
        {
            a_real_gradient_.exec();
            a_imag_gradient_.exec();
            div_a_real_diagnostic_.exec();
            div_a_imag_diagnostic_.exec();
            div_sigma_a_diagnostic_.exec();
        }
    }

    if (options_.use_phi_gauge_penalty)
    {
        phi_gauge_penalty_.exec();
    }

    if (options_.use_a_divergence_penalty)
    {
        if (options_.a_divergence_penalty_mode == AphiADivergencePenaltyMode::PairwiseUncorrected)
        {
            div_a_real_pairwise_.exec();
            div_a_imag_pairwise_.exec();
            div_a_real_gradient_pairwise_.exec();
            div_a_imag_gradient_pairwise_.exec();
        }
        else
        {
            linear_correction_matrix_.exec();
            a_real_gradient_.exec();
            a_imag_gradient_.exec();
            div_a_real_diagnostic_.exec();
            div_a_imag_diagnostic_.exec();
            div_a_real_gradient_.exec();
            div_a_imag_gradient_.exec();
        }
        a_divergence_penalty_.exec();
    }
}

template <class ExecutionPolicy>
inline AphiAssembleLhsDebugContactDynamicsBundle<ExecutionPolicy>::AphiAssembleLhsDebugContactDynamicsBundle(
    SPHBody &sph_body, Inner<> &inner_relation, Contact<> &contact_relation,
    const AphiVariableNames &variable_names, const AphiLhsAssemblyOptions &options)
    : options_(options), zero_lhs_(sph_body, variable_names.lhs),
      reaction_(sph_body, options.omega, variable_names),
      laplace_phi_real_(DynamicsArgs(inner_relation, variable_names.solution.phi_real, variable_names.material.sigma,
                                     variable_names.lhs.phi_real),
                        DynamicsArgs(contact_relation, variable_names.solution.phi_real, variable_names.material.sigma,
                                     variable_names.lhs.phi_real)),
      laplace_phi_imag_(DynamicsArgs(inner_relation, variable_names.solution.phi_imag, variable_names.material.sigma,
                                     variable_names.lhs.phi_imag),
                        DynamicsArgs(contact_relation, variable_names.solution.phi_imag, variable_names.material.sigma,
                                     variable_names.lhs.phi_imag)),
      laplace_a_real_(DynamicsArgs(inner_relation, variable_names.solution.a_real, variable_names.material.nu,
                                   variable_names.lhs.a_real),
                      DynamicsArgs(contact_relation, variable_names.solution.a_real, variable_names.material.nu,
                                   variable_names.lhs.a_real)),
      laplace_a_imag_(DynamicsArgs(inner_relation, variable_names.solution.a_imag, variable_names.material.nu,
                                   variable_names.lhs.a_imag),
                      DynamicsArgs(contact_relation, variable_names.solution.a_imag, variable_names.material.nu,
                                   variable_names.lhs.a_imag)),
      grad_phi_pairwise_(DynamicsArgs(inner_relation, variable_names),
                         DynamicsArgs(contact_relation, variable_names)),
      div_sigma_a_pairwise_(DynamicsArgs(inner_relation, options.omega, variable_names),
                            DynamicsArgs(contact_relation, options.omega, variable_names)),
      phi_gauge_penalty_(sph_body, variable_names.solution, variable_names.lhs, options.phi_gauge_penalty)
{
}

template <class ExecutionPolicy>
inline void AphiAssembleLhsDebugContactDynamicsBundle<ExecutionPolicy>::exec()
{
    zero_lhs_.exec();

    if (options_.terms.laplace_phi)
    {
        laplace_phi_real_.exec();
        laplace_phi_imag_.exec();
    }

    if (options_.terms.laplace_a)
    {
        laplace_a_real_.exec();
        laplace_a_imag_.exec();
    }

    if (options_.terms.reaction)
    {
        reaction_.exec();
    }

    if (options_.terms.grad_phi_coupling &&
        options_.grad_phi_mode == AphiGradPhiCouplingMode::PairwiseUncorrected)
    {
        grad_phi_pairwise_.exec();
    }

    if (options_.terms.div_sigma_a_coupling && options_.div_mode == AphiDivCouplingMode::PairwiseUncorrected)
    {
        div_sigma_a_pairwise_.exec();
    }

    if (options_.use_phi_gauge_penalty)
    {
        phi_gauge_penalty_.exec();
    }
}

} // namespace electromagnetics
} // namespace SPH

#endif // APHI_ASSEMBLE_LHS_DEBUG_CK_HPP
