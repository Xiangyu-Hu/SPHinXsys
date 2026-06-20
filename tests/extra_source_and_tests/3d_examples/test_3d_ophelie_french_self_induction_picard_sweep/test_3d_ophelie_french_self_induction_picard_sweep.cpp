/**
 * @file test_3d_ophelie_french_self_induction_picard_sweep.cpp
 * @brief French reload: complex Picard convergence sweep over relax × max_iter (diagnostic).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_aind_diagnostic.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

struct SweepRow
{
    Real relax = 0.0;
    size_t max_iter = 0;
    size_t iters_used = 0;
    Real final_j_rel = 0.0;
    Real phi_eq_res_vol = 0.0;
    bool picard_converged = false;
};

inline bool reloadXmlExists(const std::string &folder)
{
    return fs::exists(fs::path(folder) / "Reload.xml");
}

inline std::string resolveDefaultFrenchReloadFolder()
{
    const StdVec<std::string> candidates = {"./reload", "../../../../../reload"};
    for (const std::string &candidate : candidates)
    {
        if (reloadXmlExists(candidate))
        {
            return candidate;
        }
    }
    return "./reload";
}

inline OphelieFrenchSelfInductionPicardResult runPicardSweepCase(
    SolidBody &glass_body, Inner<> &glass_inner, const OphelieGlassFieldNames &glass_names, OphelieParameters &params,
    const OphelieFrenchReducedCaseParams &french, Real relax, size_t max_iter)
{
    params.self_induction_relaxation_factor_ = relax;
    params.self_induction_max_iterations_ = max_iter;
    return runFrenchReducedSelfInductionPicard<MainExecutionPolicy>(glass_body, glass_inner, glass_names, params,
                                                                    french);
}

} // namespace

int main(int ac, char *av[])
{
    OphelieParameters params;
    OphelieFrenchReducedCaseParams french;
    applyFrenchReducedDefaults(params, french);

    bool use_reload = false;
    std::string reload_dir;
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--reload=1") == 0 || std::strcmp(av[i], "--reload") == 0)
        {
            use_reload = true;
            french.dp = 0.02;
        }
        else if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            reload_dir = std::string(av[i] + 13);
            use_reload = true;
        }
    }
    if (use_reload && reload_dir.empty())
    {
        reload_dir = resolveDefaultFrenchReloadFolder();
    }

    OphelieTestCliOptions cli_options;
    (void)filterOphelieTestCommandLine(ac, av, params, cli_options);
    if (!cli_options.reload_dir.empty())
    {
        reload_dir = cli_options.reload_dir;
        use_reload = true;
    }

    params.enable_phi_correction_ = true;
    params.enable_self_induction_ = true;
    params.ophelie_current_form_ = OphelieCurrentFormKind::EdgeFlux;
    params.edge_flux_complex_ = true;
    params.self_induction_j_tolerance_ = 0.05;
    params.self_induction_phi_eq_res_tolerance_ = 0.01;
    syncFrenchReducedToParameters(french, params);

    if (!use_reload || !reloadXmlExists(reload_dir))
    {
        std::cerr << "test_3d_ophelie_french_self_induction_picard_sweep: Reload.xml required under \"" << reload_dir
                  << "\"\n";
        return 1;
    }

    const BoundingBoxd bounds = frenchReducedDomainBounds(french, 3.0 * french.dp);
    SPHSystem sph_system(bounds, french.dp);
    sph_system.setReloadParticles(true);
    IO::getEnvironment().resetReloadFolder(reload_dir, true);

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape().correctLevelSetSign().cleanLevelSet();
    glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;
    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    const StdVec<Real> relax_values = {Real(0.15), Real(0.3), Real(0.45)};
    const StdVec<size_t> max_iter_values = {6, 8, 12};
    StdVec<SweepRow> rows;
    rows.reserve(relax_values.size() * max_iter_values.size());

    for (Real relax : relax_values)
    {
        for (size_t max_iter : max_iter_values)
        {
            const OphelieFrenchSelfInductionPicardResult result =
                runPicardSweepCase(glass_body, *glass_inner, glass_names, params, french, relax, max_iter);
            SweepRow row;
            row.relax = relax;
            row.max_iter = max_iter;
            row.iters_used = result.self_induction_iterations;
            row.final_j_rel = result.final_j_rel_change;
            row.phi_eq_res_vol = result.phi_eq_res_vol;
            row.picard_converged = result.picard_converged;
            rows.push_back(row);
            std::cout << "picard_sweep relax=" << relax << " max_iter=" << max_iter << " iters_used=" << row.iters_used
                      << " J_rel=" << row.final_j_rel << " phi_eq_res_vol=" << row.phi_eq_res_vol
                      << " picard_converged=" << (row.picard_converged ? 1 : 0) << std::endl;
        }
    }

    size_t converged_count = 0;
    bool reference_converged = false;
    for (const SweepRow &row : rows)
    {
        if (row.picard_converged)
        {
            ++converged_count;
        }
        if (std::abs(row.relax - Real(0.3)) < TinyReal && row.max_iter == 8)
        {
            reference_converged = row.picard_converged;
        }
    }

    auto sufficientPicardBudget = [](const SweepRow &row) {
        return (row.relax >= Real(0.3) && row.max_iter >= 8) ||
               (row.relax < Real(0.3) && row.max_iter >= 12);
    };
    size_t budget_converged = 0;
    size_t budget_total = 0;
    for (const SweepRow &row : rows)
    {
        if (!sufficientPicardBudget(row))
        {
            continue;
        }
        ++budget_total;
        if (row.picard_converged)
        {
            ++budget_converged;
        }
    }
    const bool passed = reference_converged && budget_converged == budget_total;
    std::cout << "test_3d_ophelie_french_self_induction_picard_sweep n_rows=" << rows.size()
              << " converged_count=" << converged_count << " budget_converged=" << budget_converged << "/"
              << budget_total << " reference_relax0.3_max8=" << (reference_converged ? 1 : 0)
              << " sweep_passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
