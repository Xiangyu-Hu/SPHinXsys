/**
 * @file test_3d_ophelie_phi_device_vector_ops.cpp
 * @brief P1: host vs SYCL device volume-weighted dot/norm/axpy (GMRES Krylov building blocks).
 */
#include "electromagnetic_ophelie.h"
#include "electromagnetic_ophelie_french_literature.h"
#include "electromagnetic_ophelie_phi_device_vector_ops.h"
#include "electromagnetic_ophelie_phi_gmres.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cmath>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <string>

namespace fs = std::filesystem;

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{

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

} // namespace

int main(int ac, char *av[])
{
    bool use_reload = false;
    bool gmres_parity = false;
    bool krylov_arnoldi = false;
    std::string reload_dir;
    Real dp = 0.08;
    for (int i = 1; i < ac; ++i)
    {
        if (std::strcmp(av[i], "--reload=1") == 0 || std::strcmp(av[i], "--reload") == 0)
        {
            use_reload = true;
            dp = 0.02;
        }
        else if (std::strncmp(av[i], "--reload-dir=", 13) == 0)
        {
            reload_dir = std::string(av[i] + 13);
        }
        else if (std::strncmp(av[i], "--dp=", 5) == 0)
        {
            dp = static_cast<Real>(std::atof(av[i] + 5));
        }
        else if (std::strcmp(av[i], "--gmres-parity") == 0)
        {
            gmres_parity = true;
        }
        else if (std::strcmp(av[i], "--krylov-arnoldi") == 0)
        {
            krylov_arnoldi = true;
        }
    }
    if (use_reload)
    {
        if (reload_dir.empty())
        {
            reload_dir = resolveDefaultFrenchReloadFolder();
        }
        if (!reloadXmlExists(reload_dir))
        {
            std::cerr << "test_3d_ophelie_phi_device_vector_ops: Reload.xml not found under \"" << reload_dir << "\"\n";
            return 1;
        }
    }

    OphelieFrenchReducedCaseParams french;
    OphelieParameters params;
    OphelieTestCliOptions cli_options;
    OphelieFrenchLiteratureProfile literature_profile;
    french.dp = dp;
    applyFrenchReducedDefaults(params, french);
    refreshFrenchReducedCoilStack(french);
    syncFrenchReducedToParameters(french, params);
    applyFrenchLiteratureMode(params, cli_options, literature_profile);

    const BoundingBoxd bounds = frenchReducedDomainBounds(french, use_reload ? 3.0 * french.dp : 2.0 * french.dp);
    SPHSystem sph_system(bounds, french.dp);
    sph_system.setReloadParticles(use_reload);
    if (use_reload)
    {
        IO::getEnvironment().resetReloadFolder(reload_dir, true);
    }

    SolidBody glass_body(sph_system,
                         makeShared<OphelieFrenchReducedGlassCylinderShape>("GlassBody", french.glass_center,
                                                                              french.glass_radius, french.glass_half_height));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    (void)defineOphelieSolidLevelSet(glass_body);

    if (use_reload)
    {
        glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());
    }
    else
    {
        glass_body.generateParticles<BaseParticles, Lattice>();
    }

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;

    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    applyMultiloopFilamentBiotToGlass(glass_body.getBaseParticles(), glass_names, french.coil, params.mu0_,
                                      params.softening_length_);
    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_a(glass_body, glass_names);
    combine_a.exec();
    setupOpheliePhiImagRhsFromASrc<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    StdVec<Real> lhs(n);
    StdVec<Real> rhs(n);
    hostReadScalarField(particles, glass_names.phi_rhs_imag, lhs.data(), n);
    hostReadScalarField(particles, glass_names.phi_laplace_diag, rhs.data(), n);

    OpheliePhiDeviceVectorOpsCheck check =
        checkHostDeviceVolWeightedVectorOps<MainExecutionPolicy>(glass_body, glass_names, particles, lhs.data(),
                                                                 rhs.data(), n);

    StdVec<Real> rand_lhs(n);
    StdVec<Real> rand_rhs(n);
    for (size_t i = 0; i < n; ++i)
    {
        rand_lhs[i] = std::sin(static_cast<Real>(i) * Real(0.13));
        rand_rhs[i] = std::cos(static_cast<Real>(i) * Real(0.07));
    }
    const OpheliePhiDeviceVectorOpsCheck rand_check =
        checkHostDeviceVolWeightedVectorOps<MainExecutionPolicy>(glass_body, glass_names, particles, rand_lhs.data(),
                                                                 rand_rhs.data(), n);

    bool passed = check.passed && rand_check.passed;

    if (krylov_arnoldi || gmres_parity)
    {
        StdVec<Real> v0(n);
        StdVec<Real> v1(n);
        StdVec<Real> w(n);
        for (size_t i = 0; i < n; ++i)
        {
            v0[i] = std::sin(static_cast<Real>(i) * Real(0.11));
            v1[i] = std::cos(static_cast<Real>(i) * Real(0.07));
            w[i] = std::sin(static_cast<Real>(i) * Real(0.05)) + Real(0.25) * v0[i];
        }
        const StdVec<const Real *> basis_host = {v0.data(), v1.data()};
        const OpheliePhiDeviceKrylovArnoldiCheck arnoldi_check =
            checkHostDeviceKrylovArnoldiStep<MainExecutionPolicy>(glass_body, glass_names, basis_host, w.data(), 1, n);
        passed = passed && arnoldi_check.passed;
        std::cout << "krylov_arnoldi host_norm=" << arnoldi_check.host_subdiagonal_norm
                  << " device_norm=" << arnoldi_check.device_subdiagonal_norm
                  << " norm_rel_err=" << arnoldi_check.norm_rel_err
                  << " krylov_arnoldi_passed=" << (arnoldi_check.passed ? 1 : 0) << std::endl;
    }

    std::cout << "test_3d_ophelie_phi_device_vector_ops n=" << n << " dp=" << french.dp
              << " reload=" << (use_reload ? 1 : 0) << " host_dot=" << check.host_dot << " device_dot=" << check.device_dot
              << " dot_rel_err=" << check.dot_rel_err << " host_norm=" << check.host_norm
              << " device_norm=" << check.device_norm << " norm_rel_err=" << check.norm_rel_err
              << " axpy_max_abs_err=" << check.axpy_max_abs_err << " rand_dot_rel_err=" << rand_check.dot_rel_err
              << " rand_norm_rel_err=" << rand_check.norm_rel_err << " vector_ops_passed=" << (passed ? 1 : 0)
              << std::endl;

    if (!gmres_parity)
    {
        return passed ? 0 : 1;
    }

    params.phi_gmres_max_outer_iterations_ = use_reload ? 40 : 25;
    params.phi_gmres_restart_dimension_ = use_reload ? 40 : 20;
    applyOpheliePhiImagLhsOperator<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);
    setupOpheliePhiImagRhsFromASrc<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    OphelieParameters params_host = params;
    params_host.phi_gmres_use_device_vector_ops_ = false;
    params_host.phi_gmres_use_device_krylov_storage_ = false;
    StateDynamics<MainExecutionPolicy, ZeroOphelieScalarFieldCK> zero_phi(glass_body, glass_names.phi_imag);
    zero_phi.exec();
    const Real eq_res_host = solvePhiImagGMRES<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params_host);

    zero_phi.exec();
    OphelieParameters params_device = params;
    params_device.phi_gmres_use_device_vector_ops_ = true;
    params_device.phi_gmres_use_device_krylov_storage_ = true;
    const Real eq_res_device = solvePhiImagGMRES<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params_device);

    const Real gmres_parity_err = std::abs(eq_res_host - eq_res_device) / (std::abs(eq_res_host) + Real(1));
    const bool gmres_parity_ok = gmres_parity_err < Real(1e-4);

    OphelieParameters params_device_ops_only = params;
    params_device_ops_only.phi_gmres_use_device_vector_ops_ = true;
    params_device_ops_only.phi_gmres_use_device_krylov_storage_ = false;
    zero_phi.exec();
    const Real eq_res_device_ops_only =
        solvePhiImagGMRES<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params_device_ops_only);
    const Real device_ops_only_err =
        std::abs(eq_res_host - eq_res_device_ops_only) / (std::abs(eq_res_host) + Real(1));
    const bool device_ops_only_ok = device_ops_only_err < Real(1e-4);

    passed = passed && gmres_parity_ok && device_ops_only_ok;
    std::cout << "gmres_parity eq_res_host=" << eq_res_host << " eq_res_device_krylov=" << eq_res_device
              << " rel_err=" << gmres_parity_err << " gmres_parity_passed=" << (gmres_parity_ok ? 1 : 0)
              << " eq_res_device_ops_only=" << eq_res_device_ops_only << " device_ops_only_err=" << device_ops_only_err
              << " device_ops_only_passed=" << (device_ops_only_ok ? 1 : 0) << " passed=" << (passed ? 1 : 0)
              << std::endl;
    return passed ? 0 : 1;
}
