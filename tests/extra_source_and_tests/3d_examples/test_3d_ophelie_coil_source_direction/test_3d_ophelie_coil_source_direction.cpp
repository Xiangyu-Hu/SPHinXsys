/**
 * @file test_3d_ophelie_coil_source_direction.cpp
 * @brief Check azimuthal JSrcReal orientation and linear J0 scaling.
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

int main(int ac, char *av[])
{
    (void)ac;
    (void)av;

    const Vecd coil_center(0.6, 0.5, 0.5);
    const Real dp = 0.2;
    const BoundingBoxd system_bounds(Vecd(0.0, 0.0, 0.0), Vecd(1.0, 1.0, 1.0));

    auto run_with_j0 = [&](Real j0, Real &max_j_norm) -> void
    {
        OphelieParameters params;
        params.coil_j0_override_ = j0;

        SPHSystem sph_system(system_bounds, dp);
        sph_system.setReloadParticles(true);
        sph_system.setRunParticleRelaxation(false);
        IO::getEnvironment().resetReloadFolder(OPHELIE_TEST_RELOAD_DIR, true);

        SolidBody coil_body(sph_system, makeShared<ComplexShape>("CoilSourceBody"));
        coil_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
        coil_body.defineMatterMaterial<Solid>();
        coil_body.generateParticles<BaseParticles, Reload>(coil_body.Name());
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();

        OphelieCoilFieldNames coil_names;
        RegisterOphelieCoilFields register_coil(coil_body, coil_names);
        (void)register_coil;

        StateDynamics<MainExecutionPolicy, InitializeOphelieCoilSourceCK> initialize_coil(coil_body, coil_names, params,
                                                                                          coil_center);
        syncCoilSourceFieldsToDevice(coil_body.getBaseParticles(), coil_names);
        initialize_coil.exec();

        max_j_norm = hostVecdFieldMax(coil_body.getBaseParticles(), coil_names.j_src_real,
                                      coil_body.getBaseParticles().TotalRealParticles());
    };

    OphelieParameters params;
    params.coil_j0_override_ = 1.0e3;
    params.coil_center_ = coil_center;

    SPHSystem sph_system(system_bounds, dp);
    sph_system.setReloadParticles(true);
    sph_system.setRunParticleRelaxation(false);
    IO::getEnvironment().resetReloadFolder(OPHELIE_TEST_RELOAD_DIR, true);

    SolidBody coil_body(sph_system, makeShared<ComplexShape>("CoilSourceBody"));
    coil_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    coil_body.defineMatterMaterial<Solid>();
    coil_body.generateParticles<BaseParticles, Reload>(coil_body.Name());
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieCoilFieldNames coil_names;
    RegisterOphelieCoilFields register_coil(coil_body, coil_names);
    (void)register_coil;

    StateDynamics<MainExecutionPolicy, InitializeOphelieCoilSourceCK> initialize_coil(coil_body, coil_names, params,
                                                                                      coil_center);
    syncCoilSourceFieldsToDevice(coil_body.getBaseParticles(), coil_names);
    initialize_coil.exec();

    BaseParticles &particles = coil_body.getBaseParticles();
    syncVariableToHost<Vecd>(particles, coil_names.j_src_real);
    const Vecd *j_src = particles.getVariableDataByName<Vecd>(coil_names.j_src_real);

    const Vecd j_x_plus = j_src[0];
    const Vecd j_y_plus = j_src[1];
    const Vecd j_x_minus = j_src[2];

    Real max_j_base = 0.0;
    run_with_j0(1.0e3, max_j_base);
    Real max_j_doubled = 0.0;
    run_with_j0(2.0e3, max_j_doubled);

    const bool direction_ok = j_x_plus[1] > TinyReal && j_y_plus[0] < -TinyReal && j_x_minus[1] < -TinyReal;
    const bool scaling_ok = max_j_doubled > 1.9 * max_j_base && max_j_doubled < 2.1 * max_j_base;
    const bool passed = direction_ok && scaling_ok;

    std::cout << "test_3d_ophelie_coil_source_direction"
              << " J(x+)=" << j_x_plus.transpose() << " J(y+)=" << j_y_plus.transpose()
              << " J(x-)=" << j_x_minus.transpose() << " maxJ_base=" << max_j_base << " maxJ_2x=" << max_j_doubled
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
