/**
 * @file test_3d_ophelie_biot_savart_direction.cpp
 * @brief Minimal Biot-Savart B-direction check: J at (+R,0) pointing +y gives Bz>0 at origin.
 */
#include "electromagnetic_ophelie.h"
#include "io_environment.h"
#include "sphinxsys.h"

#include <cmath>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

#ifndef OPHELIE_TEST_RELOAD_DIR
#define OPHELIE_TEST_RELOAD_DIR "./reload"
#endif

int main(int ac, char *av[])
{
    (void)ac;
    (void)av;

    OphelieParameters params;
    params.current_amplitude_ = 1.0;
    params.number_of_turns_ = 1.0;
    params.coil_j0_override_ = 1.0e6;
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

    syncCoilSourceFieldsToDevice(coil_body.getBaseParticles(), coil_names);
    syncGlassElectromagneticFieldsToDevice(glass_body.getBaseParticles(), glass_names);
    assign_sigma.exec();
    initialize_coil.exec();
    compute_biot.exec();
    combine_a.exec();

    syncVariableToHost<Vecd>(glass_body.getBaseParticles(), glass_names.b_src_real);
    const Vecd *b_src = glass_body.getBaseParticles().getVariableDataByName<Vecd>(glass_names.b_src_real);
    const Vecd b0 = b_src[0];
    const Real bz = b0[2];
    const Real b_norm = b0.norm();

    const bool passed = bz > 0.0 && b_norm > TinyReal;

    std::cout << "test_3d_ophelie_biot_savart_direction B_src=" << b0.transpose() << " Bz=" << bz
              << " |B|=" << b_norm << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
