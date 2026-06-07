#ifndef ALL_ELECTROMAGNETIC_DYNAMICS_CK_H
#define ALL_ELECTROMAGNETIC_DYNAMICS_CK_H

#include "electromagnetic_dynamics/diagnostics/aphi_assemble_lhs_debug_ck.hpp"
#include "electromagnetic_dynamics/aphi_block_jacobi_preconditioner_ck.hpp"
#include "electromagnetic_dynamics/aphi_gmres_solver_ck.hpp"
#include "electromagnetic_dynamics/aphi_gmres_workspace_ck.hpp"
#include "electromagnetic_dynamics/aphi_joule_heating_ck.hpp"
#include "electromagnetic_dynamics/aphi_matrix_free_solve_ck.hpp"
#include "electromagnetic_dynamics/benchmark/aphi_benchmark_case_ck.hpp"
#include "electromagnetic_dynamics/aphi_block_zero_ck.hpp"
#include "electromagnetic_dynamics/aphi_block_vector_ops_ck.hpp"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_div_sigma_a_coupling_ck.hpp"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"
#include "electromagnetic_dynamics/aphi_grad_phi_coupling_ck.hpp"
#include "electromagnetic_dynamics/diagnostics/aphi_gradient_divergence_debug_ck.hpp"
#include "electromagnetic_dynamics/aphi_laplace_ck.hpp"
#include "electromagnetic_dynamics/aphi_matrix_free_operator_ck.hpp"
#include "electromagnetic_dynamics/aphi_pairwise_div_a_ck.hpp"
#include "electromagnetic_dynamics/aphi_a_divergence_penalty_ck.hpp"
#include "electromagnetic_dynamics/aphi_contact_pairwise_a_divergence_penalty_pipeline.h"
#include "electromagnetic_dynamics/aphi_phi_gauge_penalty_ck.hpp"
#include "electromagnetic_dynamics/aphi_reaction_ck.hpp"
#include "electromagnetic_dynamics/aphi_variables_ck.hpp"

#endif // ALL_ELECTROMAGNETIC_DYNAMICS_CK_H
