/**
 * @file electromagnetic_ophelie_analytic_cylinder_em_gates.h
 * @brief Umbrella for analytic-cylinder / EM round-acceptance gates.
 *
 * Individual gate headers live under diagnostics/gates/ (Doc199–Doc218 probes).
 * Production φ operators are NOT included here — see:
 *   electromagnetic_ophelie_phi_variational.h
 *   electromagnetic_ophelie_phi_p1_projection.h
 *   electromagnetic_ophelie_phi_d2_local_poly.h
 */
#ifndef ELECTROMAGNETIC_OPHELIE_ANALYTIC_CYLINDER_EM_GATES_H
#define ELECTROMAGNETIC_OPHELIE_ANALYTIC_CYLINDER_EM_GATES_H

#include "electromagnetic_ophelie_doc199_e2_grad_a_zero_j.h"
#include "electromagnetic_ophelie_doc199_e3_gauge_transform.h"
#include "electromagnetic_ophelie_doc199_edge_continuity_audit.h"
#include "electromagnetic_ophelie_doc199_p1_power_decomposition.h"
#include "electromagnetic_ophelie_doc201_a1_phase_rotate.h"
#include "electromagnetic_ophelie_doc201_unified_candidate.h"
#include "electromagnetic_ophelie_doc212_ec0_box_insulation.h"
#include "electromagnetic_ophelie_doc213_cylinder_volume_renorm.h"
#include "electromagnetic_ophelie_doc213_ec1_thin_plate.h"
#include "electromagnetic_ophelie_doc214_e0_constant_a_invariance.h"
#include "electromagnetic_ophelie_doc215_e1d_zone_residual.h"
#include "electromagnetic_ophelie_doc216_e1_operator_export.h"
#include "electromagnetic_ophelie_doc218_d2_manufactured.h"

#endif // ELECTROMAGNETIC_OPHELIE_ANALYTIC_CYLINDER_EM_GATES_H
