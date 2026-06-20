/**
 * @file test_3d_ophelie_french_em_dp_scan.cpp
 * @brief French reduced geometry: edge-flux EM (and optional Picard) vs particle spacing dp.
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_french_reduced_geometry.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

struct EmDpScanCli
{
    bool scan = false;
    bool picard = false;
    Real single_dp = 0.08;
    size_t picard_max_iter = 8;
    std::string csv_path = "./output/ophelie_french_em_dp_scan.csv";
};

struct EmDpScanRow
{
    Real dp = 0.0;
    size_t n = 0;
    Real phi_eq_res_vol = 0.0;
    Real joule_power_w = 0.0;
    Real div_j_l2_red = 0.0;
    bool picard_mode = false;
    bool picard_converged = false;
    size_t picard_iters = 0;
    Real picard_j_rel = 0.0;
    Real a_ind_over_a_coil = 0.0;
};

inline void applyEmDpScanCli(int ac, char *av[], EmDpScanCli &cli)
{
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--em-dp-scan") == 0 || std::strcmp(av[i], "--em-dp-scan=1") == 0)
        {
            cli.scan = true;
        }
        else if (std::strcmp(av[i], "--em-dp-scan-picard") == 0 || std::strcmp(av[i], "--em-dp-scan-picard=1") == 0)
        {
            cli.scan = true;
            cli.picard = true;
        }
        else if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            cli.single_dp = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strncmp(av[i], "--self-induction-max-iter=", 26) == 0)
        {
            cli.picard_max_iter = static_cast<size_t>(std::atoi(av[i] + 26));
        }
        else if (std::strncmp(av[i], "--em-dp-scan-csv=", 17) == 0)
        {
            cli.csv_path = std::string(av[i] + 17);
        }
    }
}

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

inline Real divJL2Reduction(const OphelieFrenchEmSolveResult &em)
{
    const Real l2_level0 = em.div_j_level0.div_j_weighted_l2;
    const Real l2_phi = em.div_j_phi.div_j_weighted_l2 > 0.0 ? em.div_j_phi.div_j_weighted_l2 : l2_level0;
    return l2_level0 / (l2_phi + TinyReal);
}

inline EmDpScanRow runFrenchEmLatticeAtDp(Real dp, OphelieParameters &params, OphelieFrenchReducedCaseParams &french,
                                        bool picard_mode, size_t picard_max_iter)
{
    french.dp = dp;
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);
    applyOphelieCoilCurrentScale(french, params);

    const BoundingBoxd bounds = frenchReducedDomainBounds(french, 2.0 * french.dp);
    SPHSystem sph_system(bounds, french.dp);
    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(glass_body);
    glass_body.generateParticles<BaseParticles, Lattice>();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    EmDpScanRow row;
    row.dp = dp;
    row.n = glass_body.getBaseParticles().TotalRealParticles();
    row.picard_mode = picard_mode;

    if (picard_mode)
    {
        params.enable_self_induction_ = true;
        params.self_induction_max_iterations_ = picard_max_iter;
        const OphelieFrenchSelfInductionPicardResult picard = runFrenchReducedSelfInductionPicard<MainExecutionPolicy>(
            glass_body, *glass_inner, glass_names, params, french);
        row.phi_eq_res_vol = picard.phi_eq_res_vol;
        row.joule_power_w = picard.joule_power_w;
        row.picard_converged = picard.picard_converged;
        row.picard_iters = picard.self_induction_iterations;
        row.picard_j_rel = picard.final_j_rel_change;
        row.a_ind_over_a_coil = picard.a_ind_over_a_coil;
    }
    else
    {
        const OphelieFrenchEmSolveResult em =
            runFrenchReducedEmPipeline<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params, french);
        row.phi_eq_res_vol = em.phi_eq_res_vol;
        row.joule_power_w = em.joule_power_recon_edge > TinyReal ? em.joule_power_recon_edge : em.joule_power_raw;
        row.div_j_l2_red = divJL2Reduction(em);
    }
    return row;
}

inline void appendEmDpScanCsv(const std::string &path, const EmDpScanRow &row)
{
    const bool write_header = !std::ifstream(path).good();
    std::ofstream out(path, std::ios::app);
    if (write_header)
    {
        out << "dp,n,phi_eq_res_vol,joule_power_w,div_j_l2_red,picard_mode,picard_converged,picard_iters,picard_j_rel,"
               "a_ind_over_a_coil\n";
    }
    out << row.dp << "," << row.n << "," << row.phi_eq_res_vol << "," << row.joule_power_w << "," << row.div_j_l2_red
        << "," << (row.picard_mode ? 1 : 0) << "," << (row.picard_converged ? 1 : 0) << "," << row.picard_iters << ","
        << row.picard_j_rel << "," << row.a_ind_over_a_coil << "\n";
}

inline bool passEmDpScanSmoke(const EmDpScanRow &row)
{
    return row.n > 0 && std::isfinite(row.phi_eq_res_vol) && std::isfinite(row.joule_power_w) &&
           row.joule_power_w > TinyReal && row.phi_eq_res_vol < Real(0.05);
}

inline bool passEmDpScanRefinement(const StdVec<EmDpScanRow> &rows)
{
    if (rows.empty())
    {
        return false;
    }
    const EmDpScanRow &finest = rows.back();
    return passEmDpScanSmoke(finest) && finest.phi_eq_res_vol < Real(0.01);
}

inline bool passPicardDpScanGate(const EmDpScanRow &row, const OphelieParameters &params)
{
    return passEmDpScanSmoke(row) && row.picard_converged &&
           ophelieSelfInductionPicardConverged(row.picard_j_rel, row.phi_eq_res_vol, params) &&
           row.a_ind_over_a_coil > TinyReal;
}

} // namespace

int main(int ac, char *av[])
{
    EmDpScanCli cli;
    applyEmDpScanCli(ac, av, cli);

    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);

    OphelieTestCliOptions cli_options;
    (void)filterOphelieTestCommandLine(ac, av, params, cli_options);

    params.enable_phi_correction_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.self_induction_j_tolerance_ = 0.05;
    params.self_induction_phi_eq_res_tolerance_ = 0.01;
    params.self_induction_relaxation_factor_ = 0.3;
    syncFrenchReducedToParameters(french, params);

    StdVec<Real> dp_list;
    if (cli.scan)
    {
        dp_list = {Real(0.08), Real(0.06), Real(0.04)};
    }
    else
    {
        dp_list = {cli.single_dp};
    }

    StdVec<EmDpScanRow> rows;
    rows.reserve(dp_list.size());
    for (Real dp : dp_list)
    {
        EmDpScanRow row = runFrenchEmLatticeAtDp(dp, params, french, cli.picard, cli.picard_max_iter);
        rows.push_back(row);
        appendEmDpScanCsv(cli.csv_path, row);
        std::cout << "em_dp_scan dp=" << row.dp << " n=" << row.n << " picard=" << (cli.picard ? 1 : 0)
                  << " phi_eq_res_vol=" << row.phi_eq_res_vol << " P_joule_W=" << row.joule_power_w
                  << " div_j_l2_red=" << row.div_j_l2_red << " picard_converged=" << (row.picard_converged ? 1 : 0)
                  << " picard_iters=" << row.picard_iters << " J_rel=" << row.picard_j_rel
                  << " A_ind_over_A_coil=" << row.a_ind_over_a_coil << std::endl;
    }

    bool smoke_passed = true;
    for (const EmDpScanRow &row : rows)
    {
        smoke_passed = smoke_passed && passEmDpScanSmoke(row);
    }

    bool refinement_passed = true;
    if (cli.scan && !cli.picard)
    {
        refinement_passed = passEmDpScanRefinement(rows);
    }

    bool picard_passed = true;
    if (cli.picard)
    {
        picard_passed = passPicardDpScanGate(rows.back(), params);
    }

    const bool passed = smoke_passed && refinement_passed && picard_passed;
    std::cout << "test_3d_ophelie_french_em_dp_scan mode=" << (cli.picard ? "picard" : (cli.scan ? "em_scan" : "em_smoke"))
              << " dp_count=" << rows.size() << " csv=" << cli.csv_path << " smoke_passed=" << (smoke_passed ? 1 : 0)
              << " refinement_passed=" << (refinement_passed ? 1 : 0)
              << " picard_passed=" << (picard_passed ? 1 : 0) << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
