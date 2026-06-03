/**
 * @file test_3d_ophelie_level0_uniform_A.cpp
 * @brief Level-0 postprocess: uniform ASrcReal gives EImag = -omega*A, JImag = sigma*E, Q = 0.5*sigma*|E|^2.
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

namespace
{
inline Real relativeError(Real measured, Real expected)
{
    return std::abs(measured - expected) / (std::abs(expected) + TinyReal);
}
} // namespace

int main(int ac, char *av[])
{
    (void)ac;
    (void)av;

    OphelieParameters params;
    params.frequency_ = 50.0;
    params.sigma_glass_ = 1.0e4;
    const Real a0 = 2.0;
    const Real omega = params.omega();
    const Vecd a_uniform(a0, 0.0, 0.0);
    const Vecd e_imag_expected(-omega * a0, 0.0, 0.0);
    const Vecd j_imag_expected = params.sigma_glass_ * e_imag_expected;
    const Real q_expected = Real(0.5) * params.sigma_glass_ * e_imag_expected.squaredNorm();

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

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;

    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    Vecd *a_src_real = particles.getVariableDataByName<Vecd>(glass_names.a_src_real);
    Vecd *a_src_imag = particles.getVariableDataByName<Vecd>(glass_names.a_src_imag);
    for (size_t i = 0; i < n; ++i)
    {
        a_src_real[i] = a_uniform;
        a_src_imag[i] = ZeroData<Vecd>::value;
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_src_real);
    syncVariableToDevice<Vecd>(particles, glass_names.a_src_imag);

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> level0(glass_body, glass_names, params);
    level0.exec();
    syncGlassElectromagneticFieldsToHost(particles, glass_names);

    Real max_e_err = 0.0;
    Real max_j_err = 0.0;
    Real max_q_err = 0.0;
    const Vecd *e_imag = particles.getVariableDataByName<Vecd>(glass_names.e_imag);
    const Vecd *j_imag = particles.getVariableDataByName<Vecd>(glass_names.j_imag);
    const Real *q = particles.getVariableDataByName<Real>(glass_names.joule_heat);
    for (size_t i = 0; i < n; ++i)
    {
        max_e_err = std::max(max_e_err, (e_imag[i] - e_imag_expected).norm() / (e_imag_expected.norm() + TinyReal));
        max_j_err = std::max(max_j_err, (j_imag[i] - j_imag_expected).norm() / (j_imag_expected.norm() + TinyReal));
        max_q_err = std::max(max_q_err, relativeError(q[i], q_expected));
    }

    const bool passed = n > 0 && max_e_err < 1.0e-6 && max_j_err < 1.0e-6 && max_q_err < 1.0e-6;
    std::cout << "test_3d_ophelie_level0_uniform_A n=" << n << " omega=" << omega << " sigma=" << params.sigma_glass_
              << " A0=" << a0 << " rel_err_E=" << max_e_err << " rel_err_J=" << max_j_err << " rel_err_Q=" << max_q_err
              << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
