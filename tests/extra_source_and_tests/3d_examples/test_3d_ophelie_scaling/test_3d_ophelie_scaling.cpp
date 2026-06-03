/**
 * @file test_3d_ophelie_scaling.cpp
 * @brief Before power normalization: J0 x2 => A/B/E/J ~2x, raw Joule integral ~4x.
 */
#include "electromagnetic_ophelie.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

#ifndef OPHELIE_TEST_RELOAD_DIR
#define OPHELIE_TEST_RELOAD_DIR "./reload"
#endif

namespace
{
struct Level0FieldPeaks
{
    Real max_a_src = 0.0;
    Real max_b_src = 0.0;
    Real max_e_imag = 0.0;
    Real max_j_imag = 0.0;
    Real joule_power_raw = 0.0;
};

inline Level0FieldPeaks runLevel0Peaks(Real j0)
{
    OphelieParameters params;
    params.coil_j0_override_ = j0;
    params.sigma_glass_ = 1.0e4;
    params.current_amplitude_ = 1.0;
    params.number_of_turns_ = 1.0;
    const Vecd coil_center(0.0, 0.0, 0.5);

    const Real dp = 0.2;
    const BoundingBoxd system_bounds(Vecd(-0.2, -0.2, -0.2), Vecd(0.4, 0.4, 0.8));

    SPHSystem sph_system(system_bounds, dp);
    sph_system.setReloadParticles(true);
    sph_system.setRunParticleRelaxation(false);
    IO::getEnvironment().resetReloadFolder(OPHELIE_TEST_RELOAD_DIR, true);

    SolidBody glass_body(sph_system, makeShared<ComplexShape>("GlassBody"));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.generateParticles<BaseParticles, Reload>(glass_body.Name());

    SolidBody coil_body(sph_system, makeShared<ComplexShape>("CoilSourceBody"));
    coil_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    coil_body.defineMatterMaterial<Solid>();
    coil_body.generateParticles<BaseParticles, Reload>(coil_body.Name());

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieCoilFieldNames coil_names;
    OphelieGlassFieldNames glass_names;
    RegisterOphelieCoilFields register_coil(coil_body, coil_names);
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_coil;
    (void)register_glass;

    StateDynamics<MainExecutionPolicy, InitializeOphelieCoilSourceCK> initialize_coil(coil_body, coil_names, params,
                                                                                      coil_center);
    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    StateDynamics<MainExecutionPolicy, ComputeOphelieCoilToGlassBiotSavartCK> compute_biot(
        glass_body, coil_body, glass_names, coil_names, params);
    StateDynamics<MainExecutionPolicy, CombineOphelieCoilAndInducedVectorPotentialCK> combine_a(glass_body, glass_names);
    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> level0(glass_body, glass_names, params);

    syncCoilSourceFieldsToDevice(coil_body.getBaseParticles(), coil_names);
    syncGlassElectromagneticFieldsToDevice(glass_body.getBaseParticles(), glass_names);
    assign_sigma.exec();
    initialize_coil.exec();
    compute_biot.exec();
    combine_a.exec();
    level0.exec();

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();
    syncGlassElectromagneticFieldsToHost(particles, glass_names);

    Level0FieldPeaks peaks;
    peaks.max_a_src = hostVecdFieldMax(particles, glass_names.a_src_real, n);
    peaks.max_b_src = hostVecdFieldMax(particles, glass_names.b_src_real, n);
    peaks.max_e_imag = hostVecdFieldMax(particles, glass_names.e_imag, n);
    peaks.max_j_imag = hostVecdFieldMax(particles, glass_names.j_imag, n);
    peaks.joule_power_raw = hostVolWeightedSum(particles, glass_names.joule_heat, n);
    return peaks;
}

inline bool ratioNear(Real ratio, Real expected, Real rtol = 0.1)
{
    return ratio > (1.0 - rtol) * expected && ratio < (1.0 + rtol) * expected;
}
} // namespace

int main(int ac, char *av[])
{
    (void)ac;
    (void)av;

    const Real j0_base = 1.0e6;
    const Level0FieldPeaks base = runLevel0Peaks(j0_base);
    const Level0FieldPeaks doubled = runLevel0Peaks(2.0 * j0_base);

    const Real ratio_a = doubled.max_a_src / (base.max_a_src + TinyReal);
    const Real ratio_b = doubled.max_b_src / (base.max_b_src + TinyReal);
    const Real ratio_e = doubled.max_e_imag / (base.max_e_imag + TinyReal);
    const Real ratio_j = doubled.max_j_imag / (base.max_j_imag + TinyReal);
    const Real ratio_p = doubled.joule_power_raw / (base.joule_power_raw + TinyReal);

    const bool fields_ok = ratioNear(ratio_a, 2.0) && ratioNear(ratio_b, 2.0) && ratioNear(ratio_e, 2.0) &&
                           ratioNear(ratio_j, 2.0);
    const bool power_ok = ratioNear(ratio_p, 4.0);
    const bool passed = fields_ok && power_ok && base.max_b_src > TinyReal;

    std::cout << "test_3d_ophelie_scaling"
              << " ratio_A=" << ratio_a << " ratio_B=" << ratio_b << " ratio_E=" << ratio_e << " ratio_J=" << ratio_j
              << " ratio_Praw=" << ratio_p << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
