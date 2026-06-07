/**
 * Stage 10.12-C1: monolithic Inner vs split two-body Inner+Contact pairwise A-penalty apply equivalence.
 */
#include "sphinxsys.h"
#include "electromagnetic_dynamics/diagnostics/aphi_contact_a_divergence_penalty_diagnostic_helpers.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

int main(int ac, char *av[])
{
    const Real max_apply_diff = 5.0e-5;
    const Real max_lhs_apply_diff = 8.0e-5;
    const AphiContactADivergencePenaltyEquivalenceMetrics metrics = runContactADivergencePenaltyEquivalenceMetrics(ac, av);
    printContactADivergencePenaltyEquivalenceMetrics("test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic",
                                                     metrics);

    const bool passed = metrics.div_a_missing == 0 && metrics.grad_div_a_stencil_safe_missing == 0 &&
                        metrics.lhs_a_stencil_safe_missing == 0 && metrics.div_a_matched > 0 &&
                        metrics.grad_div_a_stencil_safe_matched > 0 && metrics.lhs_a_stencil_safe_matched > 0 &&
                        metrics.div_a_max_abs_diff < max_apply_diff &&
                        metrics.grad_div_a_stencil_safe_max_abs_diff < max_apply_diff &&
                        metrics.lhs_a_stencil_safe_max_abs_diff < max_lhs_apply_diff;

    std::cout << "test_3d_aphi_ck_contact_a_divergence_penalty_pairwise_operator_diagnostic passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
