/**
 * @file test_3d_ophelie_phi_rhs_constant_A.cpp
 * @brief Uniform ASrcReal -> PhiImag RHS = -div(omega sigma A) ~ 0.
 */
#include "electromagnetic_ophelie.h"
#include "sphinxsys.h"

#include <cmath>
#include <iostream>

using namespace SPH;
using namespace SPH::electromagnetics::ophelie;
using MainExecutionPolicy = execution::MainExecutionPolicy;

namespace
{
class OphelieTestGlassBoxShape : public ComplexShape
{
  public:
    OphelieTestGlassBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};
} // namespace

int main(int, char *[])
{
    OphelieParameters params;
    params.frequency_ = 50.0;
    params.sigma_glass_ = 1.0e4;
    const Vecd a_uniform(1.5, 0.0, 0.0);
    const Real rhs_scale = params.omega() * params.sigma_glass_ * a_uniform.norm();

    const Real dp = 0.08;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.16, 0.16, 0.16);
    const BoundingBoxd system_bounds(center - halfsize - Vecd(dp, dp, dp), center + halfsize + Vecd(dp, dp, dp));

    SPHSystem sph_system(system_bounds, dp);
    SolidBody glass_body(sph_system, makeShared<OphelieTestGlassBoxShape>("GlassBody", center, halfsize));
    glass_body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    glass_body.defineMatterMaterial<Solid>();
    glass_body.defineBodyLevelSetShape();
    glass_body.generateParticles<BaseParticles, Lattice>();

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    OphelieGlassFieldNames glass_names;
    RegisterOphelieGlassFields register_glass(glass_body, glass_names);
    (void)register_glass;

    UniquePtr<Inner<>> glass_inner = makeUnique<Inner<>>(glass_body);
    BaseParticles &particles = glass_body.getBaseParticles();
    const size_t n = particles.TotalRealParticles();

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    Vecd *a_src_real = particles.getVariableDataByName<Vecd>(glass_names.a_src_real);
    for (size_t i = 0; i < n; ++i)
    {
        a_src_real[i] = a_uniform;
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_src_real);

    setupOpheliePhiImagRhsProblem<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, params);

    syncVariableToHost<Real>(particles, glass_names.phi_rhs_imag);
    const Real *rhs = particles.getVariableDataByName<Real>(glass_names.phi_rhs_imag);
    Real max_rhs = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        max_rhs = std::max(max_rhs, std::abs(rhs[i]));
    }
    const Real rhs_rel = max_rhs / (rhs_scale + TinyReal);

    const bool passed = n > 0 && rhs_rel < 0.05;
    std::cout << "test_3d_ophelie_phi_rhs_constant_A n=" << n << " rhs_scale=" << rhs_scale
              << " max_rhs=" << max_rhs << " rhs_rel=" << rhs_rel << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
