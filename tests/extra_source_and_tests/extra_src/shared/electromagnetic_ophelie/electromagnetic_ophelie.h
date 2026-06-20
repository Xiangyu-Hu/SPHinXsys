#ifndef ELECTROMAGNETIC_OPHELIE_H
#define ELECTROMAGNETIC_OPHELIE_H

// --- Core infrastructure (main directory) ---
#include "electromagnetic_ophelie_parameters.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_register_fields.h"
#include "electromagnetic_ophelie_cli.h"
#include "electromagnetic_ophelie_diagnostics.h"

// --- Field sources & geometry ---
#include "electromagnetic_ophelie_biot_savart.h"
#include "electromagnetic_ophelie_source.h"
#include "electromagnetic_ophelie_current_moment_source.h"
#include "electromagnetic_ophelie_multiloop_source.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"

// --- Phi solver (div-grad fallback + edge-flux production) ---
#include "electromagnetic_ophelie_laplace.h"
#include "electromagnetic_ophelie_phi_boundary.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_gmres.h"

// --- Edge-flux current formulation (Stage 1–2 production) ---
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_edge_flux_diagnostics.h"

// --- Postprocess ---
#include "electromagnetic_ophelie_postprocess.h"

// --- French literature pipeline ---
#include "electromagnetic_ophelie_french_literature.h"

// --- legacy/div_grad (particle-gradient fallback) ---
#include "electromagnetic_ophelie_phi_gradient.h"

// Optional diagnostics / benchmarks — include explicitly when needed:
//   diagnostics/electromagnetic_ophelie_aind_diagnostic.h
//   diagnostics/electromagnetic_ophelie_phi_boundary_diagnostics.h
//   diagnostics/electromagnetic_ophelie_vector_divergence_diagnostics.h
//   team7/electromagnetic_ophelie_team7_geometry.h
//   team7/electromagnetic_ophelie_team7_native_geometry.h
//   stage2/electromagnetic_ophelie_self_induction.h

#endif // ELECTROMAGNETIC_OPHELIE_H
