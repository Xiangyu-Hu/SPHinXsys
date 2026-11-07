/**
 * @file test_3d_ophelie_analytic_cylinder_skin.cpp
 * @brief Doc 23 §14.1 — long cylinder + prescribed uniform B₀ + Picard vs Bessel analytic.
 *
 * Excitation: A_coil = ½ B₀ (−y,x,0) (uniform axial B). Self-induction Picard updates A_ind.
 * Analytic: infinite cylinder H_z = H₀ J₀(αr)/J₀(αR). Compare mid-plane |J_θ|, Q shape, P/H.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_analytic_cylinder_self_induction.h"
#include "electromagnetic_ophelie_analytic_induction_weak_feedback.h"
#include "electromagnetic_ophelie_doc199_e2_grad_a_zero_j.h"
#include "electromagnetic_ophelie_doc199_e3_gauge_transform.h"
#include "electromagnetic_ophelie_doc199_edge_continuity_audit.h"
#include "electromagnetic_ophelie_doc212_ec0_box_insulation.h"
#include "electromagnetic_ophelie_doc213_ec1_thin_plate.h"
#include "electromagnetic_ophelie_doc199_p1_power_decomposition.h"
#include "electromagnetic_ophelie_doc201_a1_phase_rotate.h"
#include "electromagnetic_ophelie_doc201_unified_candidate.h"
#include "electromagnetic_ophelie_e1_phi_rhs_audit.h"
#include "electromagnetic_ophelie_e1_prescribed_a_consistency.h"
#include "electromagnetic_ophelie_edge_flux.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "electromagnetic_ophelie_french_skin_depth_diagnostic.h"
#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_phi_boundary.h"
#include "electromagnetic_ophelie_phi_component.h"
#include "electromagnetic_ophelie_self_induction.h"
#include "electromagnetic_ophelie_team7_boundary_normal.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

#if SPHINXSYS_USE_SYCL
inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape(par_ck).correctLevelSetSign().cleanLevelSet();
}
#else
inline LevelSetShape &defineOphelieSolidLevelSet(SolidBody &body)
{
    return body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
}
#endif

struct SkinLocalCli
{
    Real dp = 0.025;
    Real radius = 0.10;
    Real half_height = 0.30; // H=0.6 → H/R=6 (long)
    Real b0 = 0.01;
    Real sigma = 16.0;
    Real frequency_hz = 300000.0;
    size_t max_iter = 15;
    Real relax = 0.3;
    Real j_tol = 0.05;
    bool no_self_induction = false; // Doc197: pure A_coil → one edge-flux (no Biot / no Picard)
    Real lattice_rotate_deg = 0.0;  // Doc197: rotate Lattice xy about axis before CLL (anisotropy swap test)
    // Doc198 E0: unset → auto (no-self→prescribed_A_no_self, else full_skin)
    bool reference_mode_set = false;
    OphelieCylinderReferenceMode reference_mode = OphelieCylinderReferenceMode::FullSkin;
    std::string dump_radial_csv;            // Doc197: mid-plane radial J/Q vs analytic
    std::string dump_midplane_particle_csv; // Doc197: per-particle mid-plane for lattice audit
    // Doc198 E1: 0=off; 2=phi0 recon-only; 3=full phi+recon (implies --no-self + prescribed ref)
    int e1_layer = 0;
    std::string dump_e1_csv;
    bool e1_phi_rhs_audit = false; // Doc198: LegacyFlux vs DivSigmaA RHS region audit
    std::string dump_e1_phi_rhs_csv;
    // Doc199 E1: linear residual vs edge-current continuity (+ CSV dumps)
    bool doc199_e1_continuity = false;
    std::string dump_doc199_particle_csv;
    std::string dump_doc199_edge_csv;
    // Doc199 E2: grad-A zero-J manufactured cancel (phase1 inject φ; phase2 solve NOT_RUN until Neumann)
    int e2_phase = 0; // 0=off, 1=inject, 2=solve (blocked)
    OphelieDoc199E2ChiKind e2_chi = OphelieDoc199E2ChiKind::LinearX;
    std::string dump_e2_particle_csv;
    // Doc199 E3: gauge transform on rotational A0 (inject analytic φ1 or production solve)
    int e3_mode = 0; // 0=off, 1=inject, 2=solve
    OphelieDoc199E2ChiKind e3_chi = OphelieDoc199E2ChiKind::LinearX;
    std::string dump_e3_particle_csv;
    // Doc199 P1: power decomposition (prescribed A, no self)
    bool doc199_p1_power = false;
    std::string dump_p1_particle_csv;
    /** Doc200 T3: smooth strictly-positive σ(x) on top of e2-phase=2 solve. */
    bool doc200_t3_var_sigma = false;
    /** Doc201 A1: solve(A) then solve(iA); compare J/Q phase isomorphism (no self). */
    bool doc201_a1_phase_rotate = false;
    std::string dump_a1_csv;
    /** Doc201: apply frozen unified candidate (legacy-flux both + norm off). Named path only. */
    bool doc201_unified_candidate = false;
    /** Doc201 U3: solve once on total A0+∇χ (forbid Doc200 T2 split). */
    bool doc201_u3_unsplit = false;
    /** Doc201 U5: phase-rotate on A0+∇χ (nonzero φ). */
    bool doc201_u5_phase_rotate = false;
    /** Doc212 E-C0: thick-box insulation Neumann (early return). */
    bool doc212_ec0_box_insulation = false;
    bool doc213_ec1_thin_plate = false;
    /** Doc212 diagnostic: force imag DivΣA (shared complex edge OFF). Not production. */
    bool doc212_imag_div_sigma = false;
    bool doc212_shared_div_sigma = false;
    bool doc212_shared_pairwise = false;
    bool doc212_shared_legacy_flux = false;
};

inline void applyLocalCli(int ac, char *av[], SkinLocalCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            cli.dp = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--radius=", 9) == 0)
        {
            cli.radius = static_cast<Real>(std::atof(av[i] + 9));
        }
        else if (std::strncmp(av[i], "--half-height=", 14) == 0)
        {
            cli.half_height = static_cast<Real>(std::atof(av[i] + 14));
        }
        else if (std::strncmp(av[i], "--b0=", 5) == 0)
        {
            cli.b0 = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--sigma=", 8) == 0)
        {
            cli.sigma = static_cast<Real>(std::atof(av[i] + 8));
        }
        else if (std::strncmp(av[i], "--frequency=", 12) == 0)
        {
            cli.frequency_hz = static_cast<Real>(std::atof(av[i] + 12));
        }
        else if (std::strncmp(av[i], "--self-induction-max-iter=", 26) == 0)
        {
            cli.max_iter = static_cast<size_t>(std::atoi(av[i] + 26));
        }
        else if (std::strncmp(av[i], "--self-induction-relax=", 23) == 0)
        {
            cli.relax = static_cast<Real>(std::atof(av[i] + 23));
        }
        else if (std::strncmp(av[i], "--self-induction-tol=", 21) == 0)
        {
            cli.j_tol = static_cast<Real>(std::atof(av[i] + 21));
        }
        else if (std::strncmp(av[i], "--dump-radial-csv=", 18) == 0)
        {
            cli.dump_radial_csv = av[i] + 18;
        }
        else if (std::strncmp(av[i], "--dump-midplane-particles-csv=", 30) == 0)
        {
            cli.dump_midplane_particle_csv = av[i] + 30;
        }
        else if (std::strcmp(av[i], "--no-self-induction") == 0)
        {
            cli.no_self_induction = true;
        }
        else if (std::strncmp(av[i], "--lattice-rotate-deg=", 21) == 0)
        {
            cli.lattice_rotate_deg = static_cast<Real>(std::atof(av[i] + 21));
        }
        else if (std::strncmp(av[i], "--reference-mode=", 17) == 0)
        {
            const char *mode = av[i] + 17;
            cli.reference_mode_set = true;
            if (std::strcmp(mode, "prescribed_A_no_self") == 0 || std::strcmp(mode, "no_self") == 0)
            {
                cli.reference_mode = OphelieCylinderReferenceMode::PrescribedANoSelf;
            }
            else if (std::strcmp(mode, "full_skin") == 0 || std::strcmp(mode, "skin") == 0)
            {
                cli.reference_mode = OphelieCylinderReferenceMode::FullSkin;
            }
            else
            {
                std::cerr << "[ophelie][cyl-analytic] unknown --reference-mode=" << mode
                          << " (use full_skin|prescribed_A_no_self)" << std::endl;
                std::exit(2);
            }
        }
        else if (std::strncmp(av[i], "--e1-layer=", 11) == 0)
        {
            cli.e1_layer = std::atoi(av[i] + 11);
            if (cli.e1_layer != 2 && cli.e1_layer != 3)
            {
                std::cerr << "[ophelie][e1] --e1-layer must be 2 (phi0 recon) or 3 (phi+recon)" << std::endl;
                std::exit(2);
            }
        }
        else if (std::strncmp(av[i], "--dump-e1-csv=", 14) == 0)
        {
            cli.dump_e1_csv = av[i] + 14;
        }
        else if (std::strcmp(av[i], "--e1-phi-rhs-audit") == 0)
        {
            cli.e1_phi_rhs_audit = true;
        }
        else if (std::strncmp(av[i], "--dump-e1-phi-rhs-csv=", 22) == 0)
        {
            cli.dump_e1_phi_rhs_csv = av[i] + 22;
            cli.e1_phi_rhs_audit = true;
        }
        else if (std::strcmp(av[i], "--doc199-e1-continuity") == 0)
        {
            cli.doc199_e1_continuity = true;
        }
        else if (std::strncmp(av[i], "--dump-doc199-particle-csv=", 27) == 0)
        {
            cli.dump_doc199_particle_csv = av[i] + 27;
            cli.doc199_e1_continuity = true;
        }
        else if (std::strncmp(av[i], "--dump-doc199-edge-csv=", 23) == 0)
        {
            cli.dump_doc199_edge_csv = av[i] + 23;
            cli.doc199_e1_continuity = true;
        }
        else if (std::strncmp(av[i], "--doc199-e2-phase=", 18) == 0)
        {
            cli.e2_phase = std::atoi(av[i] + 18);
            if (cli.e2_phase != 1 && cli.e2_phase != 2)
            {
                std::cerr << "[ophelie][doc199][e2] --doc199-e2-phase must be 1 (inject) or 2 (solve)" << std::endl;
                std::exit(2);
            }
        }
        else if (std::strncmp(av[i], "--e2-chi=", 9) == 0)
        {
            const char *k = av[i] + 9;
            if (std::strcmp(k, "linear") == 0 || std::strcmp(k, "linear_x") == 0)
                cli.e2_chi = OphelieDoc199E2ChiKind::LinearX;
            else if (std::strcmp(k, "quadratic") == 0 || std::strcmp(k, "quadratic_z") == 0)
                cli.e2_chi = OphelieDoc199E2ChiKind::QuadraticZ;
            else
            {
                std::cerr << "[ophelie][doc199][e2] unknown --e2-chi=" << k << " (linear|quadratic)" << std::endl;
                std::exit(2);
            }
        }
        else if (std::strncmp(av[i], "--dump-e2-particle-csv=", 23) == 0)
        {
            cli.dump_e2_particle_csv = av[i] + 23;
        }
        else if (std::strncmp(av[i], "--doc199-e3-mode=", 17) == 0)
        {
            const char *m = av[i] + 17;
            if (std::strcmp(m, "inject") == 0 || std::strcmp(m, "1") == 0)
                cli.e3_mode = 1;
            else if (std::strcmp(m, "solve") == 0 || std::strcmp(m, "2") == 0)
                cli.e3_mode = 2;
            else
            {
                std::cerr << "[ophelie][doc199][e3] --doc199-e3-mode must be inject|solve" << std::endl;
                std::exit(2);
            }
        }
        else if (std::strncmp(av[i], "--e3-chi=", 9) == 0)
        {
            const char *k = av[i] + 9;
            if (std::strcmp(k, "linear") == 0 || std::strcmp(k, "linear_x") == 0)
                cli.e3_chi = OphelieDoc199E2ChiKind::LinearX;
            else if (std::strcmp(k, "quadratic") == 0 || std::strcmp(k, "quadratic_z") == 0)
                cli.e3_chi = OphelieDoc199E2ChiKind::QuadraticZ;
            else
            {
                std::cerr << "[ophelie][doc199][e3] unknown --e3-chi=" << k << " (linear|quadratic)" << std::endl;
                std::exit(2);
            }
        }
        else if (std::strncmp(av[i], "--dump-e3-particle-csv=", 23) == 0)
        {
            cli.dump_e3_particle_csv = av[i] + 23;
        }
        else if (std::strcmp(av[i], "--doc199-p1-power") == 0)
        {
            cli.doc199_p1_power = true;
        }
        else if (std::strncmp(av[i], "--dump-p1-particle-csv=", 23) == 0)
        {
            cli.dump_p1_particle_csv = av[i] + 23;
            cli.doc199_p1_power = true;
        }
        else if (std::strcmp(av[i], "--doc200-t3-var-sigma") == 0)
        {
            cli.doc200_t3_var_sigma = true;
            cli.e2_phase = 2; // T3 reuses T1 solve path
        }
        else if (std::strcmp(av[i], "--doc201-a1-phase-rotate") == 0)
        {
            cli.doc201_a1_phase_rotate = true;
        }
        else if (std::strncmp(av[i], "--dump-a1-csv=", 14) == 0)
        {
            // Path only — do not auto-enable A1 (would steal Doc201 U5 which reuses the dump).
            cli.dump_a1_csv = av[i] + 14;
        }
        else if (std::strcmp(av[i], "--doc201-unified-candidate") == 0)
        {
            cli.doc201_unified_candidate = true;
        }
        else if (std::strcmp(av[i], "--doc201-u3-unsplit") == 0)
        {
            cli.doc201_u3_unsplit = true;
            cli.doc201_unified_candidate = true;
        }
        else if (std::strcmp(av[i], "--doc201-u5-phase-rotate") == 0)
        {
            cli.doc201_u5_phase_rotate = true;
            cli.doc201_unified_candidate = true;
        }
        else if (std::strcmp(av[i], "--doc212-ec0-box-insulation") == 0)
        {
            cli.doc212_ec0_box_insulation = true;
        }
        else if (std::strcmp(av[i], "--doc213-ec1-thin-plate") == 0)
        {
            cli.doc213_ec1_thin_plate = true;
        }
        else if (std::strcmp(av[i], "--doc212-imag-rhs=div-sigma-a") == 0)
        {
            cli.doc212_imag_div_sigma = true;
        }
        else if (std::strcmp(av[i], "--doc212-complex-rhs=div-sigma-a") == 0)
        {
            cli.doc212_shared_div_sigma = true;
        }
        else if (std::strcmp(av[i], "--doc212-complex-rhs=legacy-pairwise") == 0)
        {
            cli.doc212_shared_pairwise = true;
        }
        else if (std::strcmp(av[i], "--doc212-complex-rhs=legacy-flux") == 0)
        {
            cli.doc212_shared_legacy_flux = true;
        }
    }
}

} // namespace

int main(int ac, char *av[])
{
    SkinLocalCli local;
    applyLocalCli(ac, av, local);
    if (local.e1_layer > 0 || local.e2_phase == 1 || local.e2_phase == 2 || local.e3_mode > 0 ||
        local.doc199_p1_power || local.doc200_t3_var_sigma || local.doc201_a1_phase_rotate ||
        local.doc201_u3_unsplit || local.doc201_u5_phase_rotate)
    {
        local.no_self_induction = true;
        if (!local.reference_mode_set)
        {
            local.reference_mode = OphelieCylinderReferenceMode::PrescribedANoSelf;
            local.reference_mode_set = true;
        }
    }
    if (!local.reference_mode_set)
    {
        local.reference_mode = local.no_self_induction ? OphelieCylinderReferenceMode::PrescribedANoSelf
                                                       : OphelieCylinderReferenceMode::FullSkin;
    }

    std::string id_msg;
    if (!verifyOphelieCylinderAnalyticWeakLimit(id_msg))
    {
        std::cerr << "cylinder analytic identity failed: " << id_msg << std::endl;
        return 1;
    }
    std::cout << "[ophelie][cyl-analytic] host_identity=" << id_msg << " passed=1" << std::endl;

    if (local.doc212_ec0_box_insulation)
    {
        OphelieDoc212Ec0BoxCli ec0;
        ec0.enabled = true;
        if (local.dp > TinyReal)
        {
            ec0.dp = local.dp;
        }
        return runOphelieDoc212Ec0BoxInsulation<MainExecutionPolicy>(ec0);
    }
    if (local.doc213_ec1_thin_plate)
    {
        OphelieDoc213Ec1ThinCli ec1;
        ec1.enabled = true;
        if (local.dp > TinyReal)
        {
            ec1.dp = local.dp;
            ec1.half[2] = Real(1.75) * local.dp; // shell=2dp overlaps both ±z; keep some interior
        }
        return runOphelieDoc213Ec1ThinPlate<MainExecutionPolicy>(ec1);
    }

    OphelieParameters params;
    params.sigma_glass_ = local.sigma;
    params.frequency_ = local.frequency_hz;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.enable_phi_correction_ = true;
    params.enable_self_induction_ = !local.no_self_induction;
    params.self_induction_max_iterations_ = local.max_iter;
    params.self_induction_relaxation_factor_ = local.relax;
    params.self_induction_j_tolerance_ = local.j_tol;
    params.aind_one_way_feedback_resolve_ = false;
    params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::SolverLocal;

    OphelieTestCliOptions cli_options;
    (void)filterOphelieTestCommandLine(ac, av, params, cli_options);
    if (local.no_self_induction)
    {
        params.enable_self_induction_ = false; // re-assert after global CLI filter
    }
    // Doc201 named unified candidate (after CLI so it wins over edge-flux DivΣA default).
    if (local.doc201_unified_candidate)
    {
        applyOphelieDoc201UnifiedCandidate(params);
        logOphelieDoc201UnifiedCandidate();
    }
    // Doc200 T1 / e2-phase=2: same pairing as Doc201 unified when not already applied.
    // Does NOT change production default (still DivΣA + boundary None).
    if (local.e2_phase == 2 && !local.doc201_unified_candidate)
    {
        params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
        params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::LegacyPairwise;
        params.phi_boundary_mode_ = OpheliePhiBoundaryMode::None;
        params.phi_rhs_project_zero_mean_ = true;
        params.phi_boundary_lhs_grad_neumann_ = false;
        params.phi_boundary_grad_neumann_projection_ = false;
        params.phi_gauge_penalty_ = Real(1.0e-8);
        params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::Off;
        params.edge_flux_solver_local_rhs_scale_ = Real(1);
        std::cout << "[ophelie][doc200][t1] pairing=legacy-pairwise+legacy-flux "
                  << "phi_boundary_mode=none phi_gauge_penalty=" << params.phi_gauge_penalty_
                  << " (opt-in T1; production DivΣA+Neumann shell remains FAIL pending redesign)"
                  << std::endl;
    }
    // After finalize (Doc213 default = shared DivΣA): optional RHS overrides.
    if (local.doc212_shared_div_sigma)
    {
        params.doc212_shared_complex_div_sigma_rhs_ = true;
        params.doc212_shared_complex_pairwise_rhs_ = false;
        params.doc212_shared_complex_edge_rhs_ = false;
        params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
        std::cout << "[ophelie][doc213] complex RHS=div-sigma-a (shared imag+real; matches production)"
                  << std::endl;
    }
    else if (local.doc212_shared_pairwise)
    {
        params.doc212_shared_complex_pairwise_rhs_ = true;
        params.doc212_shared_complex_div_sigma_rhs_ = false;
        params.doc212_shared_complex_edge_rhs_ = false;
        params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
        std::cout << "[ophelie][doc213] DIAGNOSTIC complex RHS=legacy-pairwise (shared imag+real; not production)"
                  << std::endl;
    }
    else if (local.doc212_shared_legacy_flux)
    {
        params.doc212_shared_complex_edge_rhs_ = true;
        params.doc212_shared_complex_div_sigma_rhs_ = false;
        params.doc212_shared_complex_pairwise_rhs_ = false;
        params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
        std::cout << "[ophelie][doc213] ROLLBACK complex RHS=legacy-flux shared EMF (Doc212 pin; not Doc213 default)"
                  << std::endl;
    }
    else if (local.doc212_imag_div_sigma)
    {
        params.doc212_shared_complex_edge_rhs_ = false;
        params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::DivSigmaA;
        std::cout << "[ophelie][doc212] DIAGNOSTIC imag RHS=div-sigma-a (shared complex edge OFF)" << std::endl;
    }
    logOphelieFinalParams(params, cli_options);

    const Vecd center(0.0, 0.0, local.half_height);
    const Real pad = Real(2) * local.dp;
    const BoundingBoxd bounds(Vecd(-local.radius - pad, -local.radius - pad, -pad),
                              Vecd(local.radius + pad, local.radius + pad, Real(2) * local.half_height + pad));
    SPHSystem sph_system(bounds, local.dp);
    SolidBody glass_body(sph_system, makeShared<OphelieFrenchReducedGlassCylinderShape>(
                                         "GlassBody", center, local.radius, local.half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(glass_body);
    glass_body.generateParticles<BaseParticles, Lattice>();
    if (std::abs(local.lattice_rotate_deg) > TinyReal)
    {
        BaseParticles &lat_particles = glass_body.getBaseParticles();
        const size_t n_lat = lat_particles.TotalRealParticles();
        syncVariableToHost<Vecd>(lat_particles, "Position");
        Vecd *lat_pos = lat_particles.getVariableDataByName<Vecd>("Position");
        const Real th = local.lattice_rotate_deg * Pi / Real(180);
        const Real cth = std::cos(th);
        const Real sth = std::sin(th);
        for (size_t i = 0; i < n_lat; ++i)
        {
            const Real dx = lat_pos[i][0] - center[0];
            const Real dy = lat_pos[i][1] - center[1];
            lat_pos[i][0] = center[0] + cth * dx - sth * dy;
            lat_pos[i][1] = center[1] + sth * dx + cth * dy;
        }
        syncVariableToDevice<Vecd>(lat_particles, "Position");
        std::cout << "[ophelie][cyl-analytic] lattice_rotate_deg=" << local.lattice_rotate_deg
                  << " (xy about cylinder axis; anisotropy swap test)" << std::endl;
    }
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    // Doc212: no-flux edge-recon boundary modes require NormalDirection/SignedDistance.
    // Cylinder already has LevelSet shape; reuse TEAM7 NormalFromBodyShapeCK helper.
    computeTeam7SolidBodyNormalsFromShape(sph_system, glass_body);
    {
        BaseParticles &pn = glass_body.getBaseParticles();
        syncVariableToHost<Real>(pn, "SignedDistance");
        syncVariableToHost<Vecd>(pn, "NormalDirection");
        const Real *sd = pn.getVariableDataByName<Real>("SignedDistance");
        const Vecd *nd = pn.getVariableDataByName<Vecd>("NormalDirection");
        size_t n_shell = 0;
        const Real w = Real(2) * local.dp;
        for (size_t i = 0; i < pn.TotalRealParticles(); ++i)
        {
            if (std::abs(sd[i]) <= w && nd[i].norm() > TinyReal)
                ++n_shell;
        }
        std::cout << "[ophelie][cyl-analytic] boundary_normals ok shell(|phi|<=" << double(w)
                  << ") n=" << n_shell << "/" << pn.TotalRealParticles() << std::endl;
    }

    OphelieGlassFieldNames names;
    RegisterOphelieGlassFields register_glass(glass_body, names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, names, params.sigma_glass_);
    assign_sigma.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    syncVariableToHost<Vecd>(particles, "Position");
    Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    if (local.doc200_t3_var_sigma)
    {
        // Smooth strictly-positive σ: σ0*(1.5 + 0.5*(x-cx)/R) ∈ [σ0, 2σ0] on the cylinder.
        syncVariableToHost<Real>(particles, names.sigma);
        Real *sigma_p = particles.getVariableDataByName<Real>(names.sigma);
        Real smin = params.sigma_glass_ * Real(10);
        Real smax = Real(0);
        for (size_t i = 0; i < n; ++i)
        {
            const Real xi = (pos[i][0] - center[0]) / std::max(local.radius, Real(TinyReal));
            const Real s = params.sigma_glass_ * (Real(1.5) + Real(0.5) * xi);
            sigma_p[i] = s;
            smin = std::min(smin, s);
            smax = std::max(smax, s);
        }
        syncVariableToDevice<Real>(particles, names.sigma);
        std::cout << "[ophelie][doc200][t3] var_sigma smooth-x sigma_range=[" << smin << "," << smax
                  << "] (strictly >0; J=0 still required)" << std::endl;
    }
    Vecd *a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
    Vecd *a_coil_imag = particles.getVariableDataByName<Vecd>(names.a_coil_imag);
    Vecd *b_coil_real = particles.getVariableDataByName<Vecd>(names.b_coil_real);
    Vecd *b_coil_imag = particles.getVariableDataByName<Vecd>(names.b_coil_imag);
    OphelieDoc199E2ChiField e2_chi_field;
    e2_chi_field.kind = local.e2_chi;
    e2_chi_field.center = center;
    e2_chi_field.g0 = ophelieDoc199E2GradChiScale(local.b0, local.radius);
    if (local.e2_phase == 1 || local.e2_phase == 2)
    {
        // A=∇χ (real), B=curl A = 0; scale |∇χ|~B0 R/2 vs rotational rim |A|.
        for (size_t i = 0; i < n; ++i)
        {
            a_coil_real[i] = ophelieDoc199E2GradChiAt(pos[i], e2_chi_field);
            a_coil_imag[i] = Vecd::Zero();
            b_coil_real[i] = Vecd::Zero();
            b_coil_imag[i] = Vecd::Zero();
        }
        std::cout << "[ophelie][doc199][e2] prescribe A=grad chi kind=" << ophelieDoc199E2ChiKindName(local.e2_chi)
                  << " |gradchi|_char=" << e2_chi_field.g0 << " (match B0*R/2) phase=" << local.e2_phase
                  << std::endl;
    }
    else
    {
        for (size_t i = 0; i < n; ++i)
        {
            a_coil_real[i] = ophelieRotationalAUniformBz(pos[i], center, local.b0);
            a_coil_imag[i] = Vecd::Zero();
            b_coil_real[i] = Vecd(0.0, 0.0, local.b0);
            b_coil_imag[i] = Vecd::Zero();
        }
    }
    syncVariableToDevice<Vecd>(particles, names.a_coil_real);
    syncVariableToDevice<Vecd>(particles, names.a_coil_imag);
    syncVariableToDevice<Vecd>(particles, names.b_coil_real);
    syncVariableToDevice<Vecd>(particles, names.b_coil_imag);

    hostZeroVecdField(particles, names.a_ind_real, n);
    hostZeroVecdField(particles, names.a_ind_imag, n);
    hostZeroVecdField(particles, names.b_ind_real, n);
    hostZeroVecdField(particles, names.b_ind_imag, n);
    syncVariableToDevice<Vecd>(particles, names.a_ind_real);
    syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
    syncVariableToDevice<Vecd>(particles, names.b_ind_real);
    syncVariableToDevice<Vecd>(particles, names.b_ind_imag);

    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine(glass_body, names);
    combine.exec();

    // U3/U5 before A1: --dump-a1-csv must not divert U5 into the A1-only path.
    if (local.doc201_u3_unsplit || local.doc201_u5_phase_rotate)
    {
        // Doc201 U3/U5: total A0+∇χ through ONE complex solve (no Doc200 T2 gradχ split).
        OphelieDoc199E2ChiField chi_field;
        chi_field.kind = local.e3_chi;
        chi_field.center = center;
        chi_field.g0 = ophelieDoc199E2GradChiScale(local.b0, local.radius);
        params.use_a_total_for_edge_flux_ = true;
        params.enable_self_induction_ = false;

        OphelieComplexEdgeFluxSolveReport base_report;
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, nullptr, 0.0, &base_report);
        syncGlassElectromagneticFieldsToHost(particles, names);
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
        syncVariableToHost<Real>(particles, names.phi_real);
        syncVariableToHost<Real>(particles, names.phi_imag);
        const Vecd *jr0p = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
        const Vecd *ji0p = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
        const Real *pr0p = particles.getVariableDataByName<Real>(names.phi_real);
        const Real *pi0p = particles.getVariableDataByName<Real>(names.phi_imag);
        const Real *q0p = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
        std::vector<Vecd> jr0(n), ji0(n);
        std::vector<Real> pr0(n), pi0(n), q0(n);
        for (size_t i = 0; i < n; ++i)
        {
            jr0[i] = jr0p[i];
            ji0[i] = ji0p[i];
            pr0[i] = pr0p[i];
            pi0[i] = pi0p[i];
            q0[i] = q0p[i];
        }
        std::cout << "[ophelie][doc201][u3] baseline A0 done phi_eq="
                  << ophelieSelfInductionPicardPhiEqResVol(base_report)
                  << " |gradchi|_char=" << chi_field.g0 << " chi=" << ophelieDoc199E2ChiKindName(chi_field.kind)
                  << std::endl;

        syncVariableToHost<Vecd>(particles, "Position");
        pos = particles.getVariableDataByName<Vecd>("Position");
        a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
        a_coil_imag = particles.getVariableDataByName<Vecd>(names.a_coil_imag);
        b_coil_real = particles.getVariableDataByName<Vecd>(names.b_coil_real);
        b_coil_imag = particles.getVariableDataByName<Vecd>(names.b_coil_imag);
        for (size_t i = 0; i < n; ++i)
        {
            const Vecd a0 = ophelieRotationalAUniformBz(pos[i], center, local.b0);
            const Vecd g = ophelieDoc199E2GradChiAt(pos[i], chi_field);
            a_coil_real[i] = a0 + g;
            a_coil_imag[i] = Vecd::Zero();
            b_coil_real[i] = Vecd(0.0, 0.0, local.b0);
            b_coil_imag[i] = Vecd::Zero();
        }
        syncVariableToDevice<Vecd>(particles, names.a_coil_real);
        syncVariableToDevice<Vecd>(particles, names.a_coil_imag);
        syncVariableToDevice<Vecd>(particles, names.b_coil_real);
        syncVariableToDevice<Vecd>(particles, names.b_coil_imag);
        hostZeroVecdField(particles, names.a_ind_real, n);
        hostZeroVecdField(particles, names.a_ind_imag, n);
        syncVariableToDevice<Vecd>(particles, names.a_ind_real);
        syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
        combine.exec();

        OphelieComplexEdgeFluxSolveReport u_report;
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, nullptr, 0.0, &u_report);
        std::cout << "[ophelie][doc201][u3] unsplit solve(A0+gradchi) done phi_eq="
                  << ophelieSelfInductionPicardPhiEqResVol(u_report)
                  << " phi_imag_rel=" << u_report.phi_imag_solver_rel_residual << std::endl;

        if (local.doc201_u3_unsplit)
        {
            const OphelieDoc199E3GaugeReport e3_rep = hostOphelieDoc199E3GaugeAudit(
                particles, names, chi_field, params.omega(), local.b0, params.sigma_glass_, local.radius,
                local.half_height, jr0, ji0, q0, "u3_unsplit", local.dump_e3_particle_csv);
            // Doc201 U3 gates: J/Q vs baseline A0 only. J0_vs_φ0-analytic is informational —
            // legacy-flux may leave discrete φ≠0 so φ=0 analytic is not the U3 contract.
            const bool passed = e3_rep.gate_dj && e3_rep.gate_dp;
            std::cout << "test_3d_ophelie_analytic_cylinder_skin doc201_u3 chi=" << e3_rep.chi_name << " n=" << n
                      << " dJ=" << e3_rep.j1_vs_j0_rel << " dP=" << e3_rep.p_rel_change
                      << " J0_an_info=" << e3_rep.j0_vs_analytic_rel
                      << " gate_j0_analytic_info=" << (e3_rep.gate_j0_analytic ? 1 : 0)
                      << " passed=" << (passed ? 1 : 0) << std::endl;
            return passed ? 0 : 1;
        }

        // U5: phase-rotate the nonzero-φ field A0+∇χ
        syncGlassElectromagneticFieldsToHost(particles, names);
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
        syncVariableToHost<Real>(particles, names.phi_real);
        syncVariableToHost<Real>(particles, names.phi_imag);
        const Vecd *jr_u = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
        const Vecd *ji_u = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
        const Real *pr_u = particles.getVariableDataByName<Real>(names.phi_real);
        const Real *pi_u = particles.getVariableDataByName<Real>(names.phi_imag);
        const Real *q_u = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
        for (size_t i = 0; i < n; ++i)
        {
            jr0[i] = jr_u[i];
            ji0[i] = ji_u[i];
            pr0[i] = pr_u[i];
            pi0[i] = pi_u[i];
            q0[i] = q_u[i];
        }

        hostOphelieMultiplyCoilFieldsByI(particles, names, n);
        hostZeroVecdField(particles, names.a_ind_real, n);
        hostZeroVecdField(particles, names.a_ind_imag, n);
        syncVariableToDevice<Vecd>(particles, names.a_ind_real);
        syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
        combine.exec();

        OphelieComplexEdgeFluxSolveReport rot_report;
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, nullptr, 0.0, &rot_report);
        std::cout << "[ophelie][doc201][u5] rotated solve(i(A0+gradchi)) done phi_eq="
                  << ophelieSelfInductionPicardPhiEqResVol(rot_report) << std::endl;

        const OphelieDoc201A1PhaseRotateReport u5 =
            hostOphelieDoc201A1ComparePhaseRotate(particles, names, jr0, ji0, pr0, pi0, q0, Real(0.01), Real(0.02),
                                                  local.dump_a1_csv);
        std::cout << "test_3d_ophelie_analytic_cylinder_skin doc201_u5 n=" << n << " dp=" << local.dp
                  << " J_rel=" << u5.j_rel_l2 << " Q_rel=" << u5.q_rel << " phi_rel=" << u5.phi_rel_l2
                  << " passed=" << (u5.passed ? 1 : 0) << std::endl;
        return u5.passed ? 0 : 1;
    }

    if (local.doc201_a1_phase_rotate)
    {
        // Doc201 A1: same production complex entry twice — A then iA. No solution swap.
        params.use_a_total_for_edge_flux_ = true;
        params.enable_self_induction_ = false;
        std::cout << "[ophelie][doc201][a1] path=production_edge_flux_complex "
                  << "phi_rhs="
                  << (params.phi_rhs_operator_kind_ == OpheliePhiRhsOperatorKind::LegacyFlux ? "legacy-flux"
                                                                                            : "div-sigma-a")
                  << " (Doc212: complex entry uses shared imag+real edge EMF; "
                     "DivΣA flag alone no longer splits chains)"
                  << std::endl;

        OphelieComplexEdgeFluxSolveReport base_report;
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, nullptr, 0.0, &base_report);
        syncGlassElectromagneticFieldsToHost(particles, names);
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
        syncVariableToHost<Real>(particles, names.phi_real);
        syncVariableToHost<Real>(particles, names.phi_imag);
        const Vecd *jr0p = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
        const Vecd *ji0p = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
        const Real *pr0p = particles.getVariableDataByName<Real>(names.phi_real);
        const Real *pi0p = particles.getVariableDataByName<Real>(names.phi_imag);
        const Real *q0p = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
        std::vector<Vecd> jr0(n), ji0(n);
        std::vector<Real> pr0(n), pi0(n), q0(n);
        for (size_t i = 0; i < n; ++i)
        {
            jr0[i] = jr0p[i];
            ji0[i] = ji0p[i];
            pr0[i] = pr0p[i];
            pi0[i] = pi0p[i];
            q0[i] = q0p[i];
        }
        std::cout << "[ophelie][doc201][a1] baseline solve(A) done phi_eq="
                  << ophelieSelfInductionPicardPhiEqResVol(base_report)
                  << " phi_imag_rel=" << base_report.phi_imag_solver_rel_residual << std::endl;

        hostOphelieMultiplyCoilFieldsByI(particles, names, n);
        hostZeroVecdField(particles, names.a_ind_real, n);
        hostZeroVecdField(particles, names.a_ind_imag, n);
        syncVariableToDevice<Vecd>(particles, names.a_ind_real);
        syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
        combine.exec();

        OphelieComplexEdgeFluxSolveReport rot_report;
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, nullptr, 0.0, &rot_report);
        std::cout << "[ophelie][doc201][a1] rotated solve(iA) done phi_eq="
                  << ophelieSelfInductionPicardPhiEqResVol(rot_report)
                  << " phi_imag_rel=" << rot_report.phi_imag_solver_rel_residual << std::endl;

        const OphelieDoc201A1PhaseRotateReport a1 =
            hostOphelieDoc201A1ComparePhaseRotate(particles, names, jr0, ji0, pr0, pi0, q0, Real(0.01), Real(0.02),
                                                  local.dump_a1_csv);
        std::cout << "test_3d_ophelie_analytic_cylinder_skin doc201_a1 n=" << n << " dp=" << local.dp
                  << " J_rel=" << a1.j_rel_l2 << " Q_rel=" << a1.q_rel << " phi_rel=" << a1.phi_rel_l2
                  << " passed=" << (a1.passed ? 1 : 0) << std::endl;
        return a1.passed ? 0 : 1;
    }

    if (local.e3_mode > 0)
    {
        // Doc199 E3: baseline A0 rotational → snapshot J0; then A1=A0+∇χ.
        OphelieDoc199E2ChiField e3_chi_field;
        e3_chi_field.kind = local.e3_chi;
        e3_chi_field.center = center;
        e3_chi_field.g0 = ophelieDoc199E2GradChiScale(local.b0, local.radius);
        params.use_a_total_for_edge_flux_ = true;
        params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::Off;
        params.edge_flux_solver_local_rhs_scale_ = Real(1);

        OphelieComplexEdgeFluxSolveReport base_report;
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, nullptr, 0.0, &base_report);
        syncGlassElectromagneticFieldsToHost(particles, names);
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
        const Vecd *jr0p = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
        const Vecd *ji0p = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);
        const Real *q0p = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
        std::vector<Vecd> jr0(n), ji0(n);
        std::vector<Real> q0(n);
        for (size_t i = 0; i < n; ++i)
        {
            jr0[i] = jr0p[i];
            ji0[i] = ji0p[i];
            q0[i] = q0p[i];
        }
        std::cout << "[ophelie][doc199][e3] baseline A0 rotational done phi_eq="
                  << ophelieSelfInductionPicardPhiEqResVol(base_report) << " chi="
                  << ophelieDoc199E2ChiKindName(local.e3_chi) << " |gradchi|_char=" << e3_chi_field.g0
                  << std::endl;

        syncVariableToHost<Vecd>(particles, "Position");
        pos = particles.getVariableDataByName<Vecd>("Position");
        a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
        a_coil_imag = particles.getVariableDataByName<Vecd>(names.a_coil_imag);
        b_coil_real = particles.getVariableDataByName<Vecd>(names.b_coil_real);
        b_coil_imag = particles.getVariableDataByName<Vecd>(names.b_coil_imag);
        for (size_t i = 0; i < n; ++i)
        {
            const Vecd a0 = ophelieRotationalAUniformBz(pos[i], center, local.b0);
            const Vecd g = ophelieDoc199E2GradChiAt(pos[i], e3_chi_field);
            a_coil_real[i] = a0 + g;
            a_coil_imag[i] = Vecd::Zero();
            b_coil_real[i] = Vecd(0.0, 0.0, local.b0); // curl(A0+∇χ)=curl A0
            b_coil_imag[i] = Vecd::Zero();
        }
        syncVariableToDevice<Vecd>(particles, names.a_coil_real);
        syncVariableToDevice<Vecd>(particles, names.a_coil_imag);
        syncVariableToDevice<Vecd>(particles, names.b_coil_real);
        syncVariableToDevice<Vecd>(particles, names.b_coil_imag);
        hostZeroVecdField(particles, names.a_ind_real, n);
        hostZeroVecdField(particles, names.a_ind_imag, n);
        syncVariableToDevice<Vecd>(particles, names.a_ind_real);
        syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
        combine.exec();

        const Real omega_e3 = params.omega();
        const char *mode_name = (local.e3_mode == 1) ? "inject" : "solve";
        if (local.e3_mode == 1)
        {
            // φ1 = -ωχ (baseline φ0≈0 for rotational A / DivΣA)
            Real *phi_i = particles.getVariableDataByName<Real>(names.phi_imag);
            Real *phi_r = particles.getVariableDataByName<Real>(names.phi_real);
            for (size_t i = 0; i < n; ++i)
            {
                phi_i[i] = -omega_e3 * ophelieDoc199E2ChiAt(pos[i], e3_chi_field);
                phi_r[i] = Real(0);
            }
            syncVariableToDevice<Real>(particles, names.phi_imag);
            syncVariableToDevice<Real>(particles, names.phi_real);
            (void)execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, names, params,
                                                                          /*skip_real_phi_solve=*/false);
            syncOphelieEdgeReconToPrimaryEJQ<MainExecutionPolicy>(glass_body, names, params);
        }
        else
        {
            // Doc200 T2: actual solve (not inject). Full A1=A0+∇χ with production DivΣA
            // collapses φ→0 (R1–R2); full A1 with legacy-flux reintroduces A0 EMF pollution.
            // Compatible path: solve φ on ∇χ only with T1 pairing (legacy-flux), then recon
            // with full A1 + that φ. Continuously E1=E0; discrete matches inject gauge test
            // but φ comes from PCG, not analytic assignment.
            const OpheliePhiRhsOperatorKind saved_rhs = params.phi_rhs_operator_kind_;
            const Real saved_gauge = params.phi_gauge_penalty_;
            params.phi_rhs_operator_kind_ = OpheliePhiRhsOperatorKind::LegacyFlux;
            params.phi_lhs_operator_kind_ = OpheliePhiLhsOperatorKind::LegacyPairwise;
            params.phi_boundary_mode_ = OpheliePhiBoundaryMode::None;
            params.phi_rhs_project_zero_mean_ = true;
            params.phi_gauge_penalty_ = Real(1.0e-8);
            params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::Off;
            params.edge_flux_solver_local_rhs_scale_ = Real(1);

            syncVariableToHost<Vecd>(particles, "Position");
            pos = particles.getVariableDataByName<Vecd>("Position");
            a_coil_real = particles.getVariableDataByName<Vecd>(names.a_coil_real);
            a_coil_imag = particles.getVariableDataByName<Vecd>(names.a_coil_imag);
            for (size_t i = 0; i < n; ++i)
            {
                a_coil_real[i] = ophelieDoc199E2GradChiAt(pos[i], e3_chi_field);
                a_coil_imag[i] = Vecd::Zero();
            }
            syncVariableToDevice<Vecd>(particles, names.a_coil_real);
            syncVariableToDevice<Vecd>(particles, names.a_coil_imag);
            hostZeroVecdField(particles, names.a_ind_real, n);
            hostZeroVecdField(particles, names.a_ind_imag, n);
            syncVariableToDevice<Vecd>(particles, names.a_ind_real);
            syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
            combine.exec();

            OphelieComplexEdgeFluxSolveReport gauge_report;
            (void)execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
                glass_body, *glass_inner, names, params, nullptr, 0.0, &gauge_report);
            std::cout << "[ophelie][doc200][t2] solved phi on grad-chi (legacy-flux) phi_eq="
                      << ophelieSelfInductionPicardPhiEqResVol(gauge_report)
                      << " phi_solver_rel=" << gauge_report.phi_imag_solver_rel_residual << std::endl;

            // Restore A1=A0+∇χ; keep solved φ; recon only.
            for (size_t i = 0; i < n; ++i)
            {
                const Vecd a0 = ophelieRotationalAUniformBz(pos[i], center, local.b0);
                const Vecd g = ophelieDoc199E2GradChiAt(pos[i], e3_chi_field);
                a_coil_real[i] = a0 + g;
                a_coil_imag[i] = Vecd::Zero();
            }
            syncVariableToDevice<Vecd>(particles, names.a_coil_real);
            syncVariableToDevice<Vecd>(particles, names.a_coil_imag);
            hostZeroVecdField(particles, names.a_ind_real, n);
            hostZeroVecdField(particles, names.a_ind_imag, n);
            syncVariableToDevice<Vecd>(particles, names.a_ind_real);
            syncVariableToDevice<Vecd>(particles, names.a_ind_imag);
            combine.exec();
            (void)execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, names, params,
                                                                          /*skip_real_phi_solve=*/false);
            syncOphelieEdgeReconToPrimaryEJQ<MainExecutionPolicy>(glass_body, names, params);
            params.phi_rhs_operator_kind_ = saved_rhs;
            params.phi_gauge_penalty_ = saved_gauge;
            std::cout << "[ophelie][doc200][t2] recon with A1=A0+gradchi + solved phi (not inject)" << std::endl;
        }

        syncGlassElectromagneticFieldsToHost(particles, names);
        syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
        syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);
        const OphelieDoc199E3GaugeReport e3_rep = hostOphelieDoc199E3GaugeAudit(
            particles, names, e3_chi_field, omega_e3, local.b0, params.sigma_glass_, local.radius,
            local.half_height, jr0, ji0, q0, mode_name, local.dump_e3_particle_csv);
        const bool passed = e3_rep.gate_dj && e3_rep.gate_dp && e3_rep.gate_j0_analytic;
        std::cout << "test_3d_ophelie_analytic_cylinder_skin doc199_e3 mode=" << mode_name
                  << " chi=" << e3_rep.chi_name << " n=" << n << " dJ=" << e3_rep.j1_vs_j0_rel
                  << " dP=" << e3_rep.p_rel_change << " J0_an=" << e3_rep.j0_vs_analytic_rel
                  << " passed=" << (passed ? 1 : 0) << std::endl;
        return passed ? 0 : 1;
    }

    if (local.e1_phi_rhs_audit || local.e1_layer == 3)
    {
        // Auto-run RHS audit on L3; optional explicit flag for L2/other.
        params.use_a_total_for_edge_flux_ = true;
        (void)execOphelieE1PrescribedAPhiRhsAudit<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, center, local.radius, local.half_height, local.b0,
            local.dump_e1_phi_rhs_csv);
    }

    Real phi_solver_rel = 0.0;
    Real phi_eq = 0.0;
    Real relaxed_step = 0.0;
    size_t iters = 0;
    bool picard_ok = false;
    Real fixed_point = 0.0;
    const char *e1_tag = "off";
    if (local.e2_phase == 1)
    {
        // Doc199 E2 phase1: inject φ_imag=-ωχ, φ_real=0; recon only (no phi solve).
        params.use_a_total_for_edge_flux_ = true;
        params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::Off;
        params.edge_flux_solver_local_rhs_scale_ = Real(1);
        const Real omega_e2 = params.omega();
        Real *phi_i = particles.getVariableDataByName<Real>(names.phi_imag);
        Real *phi_r = particles.getVariableDataByName<Real>(names.phi_real);
        for (size_t i = 0; i < n; ++i)
        {
            phi_i[i] = -omega_e2 * ophelieDoc199E2ChiAt(pos[i], e2_chi_field);
            phi_r[i] = Real(0);
        }
        syncVariableToDevice<Real>(particles, names.phi_imag);
        syncVariableToDevice<Real>(particles, names.phi_real);
        (void)execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, names, params,
                                                                      /*skip_real_phi_solve=*/false);
        syncOphelieEdgeReconToPrimaryEJQ<MainExecutionPolicy>(glass_body, names, params);
        iters = 1;
        picard_ok = true;
        fixed_point = 0.0;
        phi_eq = 0.0;
        e1_tag = "E2_phase1_inject";
        std::cout << "[ophelie][doc199][e2] phase=1 mode=inject_analytic_phi_recon_only" << std::endl;
    }
    else if (local.e1_layer == 2)
    {
        // Doc198 E1-L2: φ≡0, skip phi solve; reconstruct E/J from edge drops of prescribed A only.
        params.use_a_total_for_edge_flux_ = true;
        params.edge_flux_normalization_mode_ = OphelieEdgeFluxNormalizationMode::Off;
        params.edge_flux_solver_local_rhs_scale_ = Real(1);
        StateDynamics<MainExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi_i(glass_body, names.phi_imag);
        StateDynamics<MainExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi_r(glass_body, names.phi_real);
        zero_phi_i.exec();
        zero_phi_r.exec();
        (void)execOphelieEdgeFluxPostPhiPipeline<MainExecutionPolicy>(glass_body, *glass_inner, names, params,
                                                                      /*skip_real_phi_solve=*/false);
        syncOphelieEdgeReconToPrimaryEJQ<MainExecutionPolicy>(glass_body, names, params);
        iters = 1;
        picard_ok = true;
        fixed_point = 0.0;
        phi_eq = 0.0;
        e1_tag = "L2_phi0_recon";
        std::cout << "[ophelie][e1] layer=2 mode=phi0_recon_only (no phi solve)" << std::endl;
    }
    else if (local.no_self_induction)
    {
        // Pure A_coil → one complex edge-flux (A_ind remains zero; no final Biot shot).
        // Doc198 E1-L3 when --e1-layer=3.
        params.use_a_total_for_edge_flux_ = true;
        OpheliePhiBoundaryGeometryContext boundary_geom;
        const OpheliePhiBoundaryGeometryContext *boundary_geom_ptr = nullptr;
        if (params.phi_boundary_mode_ == OpheliePhiBoundaryMode::OneSidedNeumann)
        {
            boundary_geom.normal_source = OpheliePhiBoundaryNormalSource::AnalyticCylinder;
            boundary_geom.french.glass_center = center;
            boundary_geom.french.glass_radius = local.radius;
            boundary_geom.french.glass_half_height = local.half_height;
            boundary_geom.french.dp = local.dp;
            const size_t n_bnd =
                setupOpheliePhiBoundaryParticleFields(particles, names, params, boundary_geom, local.dp);
            boundary_geom_ptr = &boundary_geom;
            std::cout << "[ophelie][doc200][p2] one-sided-neumann n_boundary=" << n_bnd
                      << " (opt-in CLI; production default remains None)" << std::endl;
        }
        OphelieComplexEdgeFluxSolveReport solve_report;
        (void)execOphelieComplexEdgeFluxSolveReconAndPower<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, boundary_geom_ptr, local.dp, &solve_report);
        phi_solver_rel = solve_report.phi_imag_solver_rel_residual;
        phi_eq = ophelieSelfInductionPicardPhiEqResVol(solve_report);
        iters = 1;
        picard_ok = true; // N/A; treat as "single-shot ok" for logging
        fixed_point = 0.0;
        e1_tag = (local.e1_layer == 3) ? "L3_phi_solve_recon" : "no_self";
        std::cout << "[ophelie][cyl-analytic] mode=no-self-induction (A_coil only, one edge-flux)"
                  << " e1=" << e1_tag << std::endl;
    }
    else
    {
        fixed_point = runOphelieSelfInductionWithPhiSolve<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, phi_solver_rel, iters, phi_eq, picard_ok, &relaxed_step);
    }

    syncGlassElectromagneticFieldsToHost(particles, names);
    syncVariableToHost<Real>(particles, names.joule_heat_edge_recon_complex);
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Real>(particles, "VolumetricMeasure");
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_real);
    syncVariableToHost<Vecd>(particles, names.j_edge_recon_imag);

    pos = particles.getVariableDataByName<Vecd>("Position");
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const Real *q = particles.getVariableDataByName<Real>(names.joule_heat_edge_recon_complex);
    const Vecd *jr = particles.getVariableDataByName<Vecd>(names.j_edge_recon_real);
    const Vecd *ji = particles.getVariableDataByName<Vecd>(names.j_edge_recon_imag);

    const Real omega = params.omega();

    if (local.e2_phase == 1 || local.e2_phase == 2)
    {
        const OphelieDoc199E2ZeroJReport e2_rep = hostOphelieDoc199E2ZeroJAudit(
            particles, names, e2_chi_field, omega, local.radius, local.half_height, local.dump_e2_particle_csv,
            local.e2_phase);
        const bool a_ok = e2_rep.a_vs_gradchi_rel < Real(1e-12);
        const bool passed = local.e2_phase == 1
                                ? (e2_rep.gate_j_all && a_ok)
                                : (e2_rep.gate_j_all && e2_rep.gate_phi_match && a_ok);
        std::cout << "test_3d_ophelie_analytic_cylinder_skin doc199_e2 phase=" << local.e2_phase
                  << " chi=" << e2_rep.chi_name << " n=" << n << " dp=" << local.dp
                  << " J_rel_all=" << e2_rep.j_rel_l2_all << " J_rel_int=" << e2_rep.j_rel_l2_interior
                  << " phi_i_rel=" << e2_rep.phi_imag_vs_minus_omega_chi_rel
                  << " phi_r_rel=" << e2_rep.phi_real_vs_zero_rel << " a_ok=" << (a_ok ? 1 : 0)
                  << " phi_eq=" << phi_eq << " phi_solver_rel=" << phi_solver_rel
                  << " passed=" << (passed ? 1 : 0) << std::endl;
        // Doc201 A4: allow e1 continuity on the same e2 solve (was early-return INVALID_RUN).
        if (local.doc199_e1_continuity || !local.dump_doc199_particle_csv.empty() ||
            !local.dump_doc199_edge_csv.empty())
        {
            (void)execOphelieDoc199EdgeContinuityAudit<MainExecutionPolicy>(
                glass_body, *glass_inner, names, params, center, local.radius, local.half_height,
                local.dump_doc199_particle_csv, local.dump_doc199_edge_csv);
            std::cout << "[ophelie][doc199] final_params_echo phi_rhs_operator="
                      << (params.phi_rhs_operator_kind_ == OpheliePhiRhsOperatorKind::LegacyFlux ? "legacy-flux"
                                                                                                : "div-sigma-a")
                      << " phi_lhs_operator="
                      << (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::LegacyPairwise
                              ? "legacy-pairwise"
                              : "div-sigma-grad")
                      << std::endl;
        }
        return passed ? 0 : 1;
    }

    if (local.doc199_p1_power)
    {
        const OphelieDoc199P1PowerReport p1 = hostOphelieDoc199P1PowerDecomposition(
            particles, names, params, center, local.radius, local.half_height, local.b0,
            local.dump_p1_particle_csv);
        // Primary: field match on same particles. Geometry/total reported separately (may fail without field fail).
        const bool passed = p1.field_rel <= Real(0.01) && p1.j_rel_l2_vs_ref <= Real(0.05);
        std::cout << "test_3d_ophelie_analytic_cylinder_skin doc199_p1 n=" << n
                  << " field_rel=" << p1.field_rel << " geometry_rel=" << p1.geometry_rel
                  << " total_vs_cont_rel=" << p1.total_vs_cont_rel << " J_rel=" << p1.j_rel_l2_vs_ref
                  << " passed=" << (passed ? 1 : 0) << std::endl;
        return passed ? 0 : 1;
    }

    if (local.doc199_e1_continuity || !local.dump_doc199_particle_csv.empty() || !local.dump_doc199_edge_csv.empty())
    {
        (void)execOphelieDoc199EdgeContinuityAudit<MainExecutionPolicy>(
            glass_body, *glass_inner, names, params, center, local.radius, local.half_height,
            local.dump_doc199_particle_csv, local.dump_doc199_edge_csv);
        std::cout << "[ophelie][doc199] final_params_echo phi_rhs_operator="
                  << (params.phi_rhs_operator_kind_ == OpheliePhiRhsOperatorKind::LegacyFlux ? "legacy-flux"
                                                                                            : "div-sigma-a")
                  << " phi_lhs_operator="
                  << (params.phi_lhs_operator_kind_ == OpheliePhiLhsOperatorKind::LegacyPairwise ? "legacy-pairwise"
                                                                                                  : "div-sigma-grad")
                  << std::endl;
    }

    const Real z_band = Real(0.15) * local.half_height; // mid-plane slab
    std::cout << "[ophelie][cyl-analytic] compare_reference_mode="
              << ophelieCylinderReferenceModeName(local.reference_mode)
              << " solver_no_self=" << (local.no_self_induction ? 1 : 0) << std::endl;
    if (local.no_self_induction && local.reference_mode == OphelieCylinderReferenceMode::FullSkin)
    {
        std::cout << "[ophelie][cyl-analytic] WARN Doc198: no-self solver vs full_skin ref is approximate only"
                  << std::endl;
    }
    if (!local.no_self_induction && local.reference_mode == OphelieCylinderReferenceMode::PrescribedANoSelf)
    {
        std::cout << "[ophelie][cyl-analytic] WARN Doc198: Picard self-induction vs prescribed_A_no_self ref"
                  << std::endl;
    }
    const OphelieCylinderMidplaneCompareReport cmp = hostOphelieCylinderMidplaneCompare(
        pos, vol, jr, ji, q, n, center, local.radius, local.half_height, z_band, local.b0, omega, params.mu0_,
        params.sigma_glass_, 12, local.reference_mode);
    logOphelieCylinderMidplaneCompare(cmp, local.no_self_induction ? "no_self" : "picard");
    if (!local.dump_radial_csv.empty())
    {
        hostOphelieCylinderDumpMidplaneRadialCsv(pos, vol, jr, ji, q, n, center, local.radius, z_band, local.b0,
                                                 omega, params.mu0_, params.sigma_glass_, local.dump_radial_csv, 12,
                                                 local.reference_mode);
    }
    if (!local.dump_midplane_particle_csv.empty())
    {
        hostOphelieCylinderDumpMidplaneParticleCsv(pos, vol, jr, ji, n, center, local.radius, z_band, local.b0,
                                                   omega, params.mu0_, params.sigma_glass_,
                                                   local.dump_midplane_particle_csv, local.reference_mode);
    }

    if (local.e1_layer > 0)
    {
        const OphelieE1PrescribedAReport e1 = hostOphelieE1PrescribedAConsistencyAudit(
            particles, names, center, local.radius, local.half_height, local.b0, omega, params.sigma_glass_, e1_tag);
        logOphelieE1PrescribedAReport(e1);
        if (!local.dump_e1_csv.empty())
        {
            hostOphelieE1DumpParticleCsv(particles, names, center, local.radius, local.half_height, local.b0, omega,
                                         params.sigma_glass_, local.dump_e1_csv);
        }
        // Pre-registered Doc198 gates (informational in stdout; do not widen).
        const Real l2_int = ophelieE1RelL2(e1.interior);
        const Real l2_all = ophelieE1RelL2(e1.all);
        const bool e1_int_ok = l2_int <= Real(0.01);
        const bool e1_all_ok = l2_all <= Real(0.05);
        std::cout << "[ophelie][e1] gates interior_L2_J<=1%=" << (e1_int_ok ? 1 : 0) << " (" << l2_int << ")"
                  << " all_L2_J<=5%=" << (e1_all_ok ? 1 : 0) << " (" << l2_all << ")"
                  << " A_ok=" << (e1.a_prescription_rel_l2 < Real(1e-12) ? 1 : 0)
                  << " affine_edge_id=" << e1.edge_affine_identity_max << std::endl;
        if (local.e1_layer == 3)
        {
            UpdateCellLinkedList<MainExecutionPolicy, RealBody> upd_cll(glass_body);
            UpdateRelation<MainExecutionPolicy, Inner<>> upd_rel(*glass_inner);
            upd_cll.exec();
            upd_rel.exec();
            InteractionDynamicsCK<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad(
                *glass_inner, names);
            compute_grad.exec();
            hostOphelieE1LogGradPhiVsOmegaA(particles, names, center, local.radius, local.half_height, omega);
        }
    }

    const Real a_coil_n = hostComplexVecdPairVolWeightedNorm(particles, names.a_coil_real, names.a_coil_imag, n);
    const Real a_ind_n = hostComplexVecdPairVolWeightedNorm(particles, names.a_ind_real, names.a_ind_imag, n);
    const Real a_ratio = a_ind_n / (a_coil_n + TinyReal);

    auto sol = makeOphelieCylinderAnalyticSolution(local.radius, local.b0, omega, params.mu0_, params.sigma_glass_);
    finalizeOphelieCylinderAnalyticPower(sol);

    // Gates: shape first; absolute P tolerant (finite height + coarse dp)
    const bool corr_j_ok = cmp.corr_abs_Jtheta > Real(0.85);
    const bool corr_q_ok = cmp.corr_Q_shape > Real(0.85);
    const bool power_ok = cmp.rel_P < Real(0.40);
    const bool j_ok = fixed_point < local.j_tol;
    const bool phi_ok = phi_eq < ophelieSelfInductionPicardPhiEqResTolerance(params);
    const bool converged = picard_ok && j_ok && phi_ok;
    const bool finite_ok = n > 0 && std::isfinite(cmp.P_sph_W) && cmp.n_bins >= 4;
    const bool passed = finite_ok && converged && corr_j_ok && corr_q_ok && power_ok;

    std::cout << "test_3d_ophelie_analytic_cylinder_skin n=" << n << " dp=" << local.dp << " R=" << local.radius
              << " H=" << (Real(2) * local.half_height) << " B0=" << local.b0 << " sigma=" << params.sigma_glass_
              << " f_hz=" << params.frequency_ << " delta/R=" << sol.delta_over_R
              << " picard_iters=" << iters << " fixed_point=" << fixed_point << " phi_eq=" << phi_eq
              << " picard_converged=" << (picard_ok ? 1 : 0) << " A_ind/A_coil=" << a_ratio
              << " corr_J=" << cmp.corr_abs_Jtheta << " corr_Q=" << cmp.corr_Q_shape << " rel_P=" << cmp.rel_P
              << " corr_j_ok=" << (corr_j_ok ? 1 : 0) << " corr_q_ok=" << (corr_q_ok ? 1 : 0)
              << " power_ok=" << (power_ok ? 1 : 0) << " converged=" << (converged ? 1 : 0)
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
