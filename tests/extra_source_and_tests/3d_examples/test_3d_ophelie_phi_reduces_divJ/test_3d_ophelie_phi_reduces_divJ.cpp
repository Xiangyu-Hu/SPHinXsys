/**
 * @file test_3d_ophelie_phi_reduces_divJ.cpp
 * @brief Level0 divJ != 0 for linear A(x); with gradPhi = -omega A, E/J path yields lower divJ.
 */
#include "electromagnetic_ophelie.h"
#include "sphinxsys.h"

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
    params.enable_power_scaling_ = false;
    const Real a_grad = 1.2;
    const Real omega = params.omega();

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
    const Real div_j_length = halfsize[0];

    StateDynamics<MainExecutionPolicy, AssignOphelieGlassSigmaCK> assign_sigma(glass_body, glass_names, params.sigma_glass_);
    assign_sigma.exec();

    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Vecd *a_src_real = particles.getVariableDataByName<Vecd>(glass_names.a_src_real);
    Vecd *grad_phi_imag = particles.getVariableDataByName<Vecd>(glass_names.grad_phi_imag);
    for (size_t i = 0; i < n; ++i)
    {
        a_src_real[i] = Vecd(a_grad * pos[i][0], 0.0, 0.0);
    }
    syncVariableToDevice<Vecd>(particles, glass_names.a_src_real);

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQFromASrcNoPhiCK> level0(glass_body, glass_names, params);
    level0.exec();
    const OphelieDivJMetrics div_j_level0 =
        computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_j_length);

    for (size_t i = 0; i < n; ++i)
    {
        grad_phi_imag[i] = -omega * a_src_real[i];
    }
    syncVariableToDevice<Vecd>(particles, glass_names.grad_phi_imag);

    StateDynamics<MainExecutionPolicy, ComputeOphelieEJQWithPhiCK> compute_ejq_with_phi(glass_body, glass_names, params);
    compute_ejq_with_phi.exec();

    const OphelieDivJMetrics div_j_phi =
        computeOphelieDivJImag<MainExecutionPolicy>(glass_body, *glass_inner, glass_names, div_j_length);

    syncGlassElectromagneticFieldsToHost(particles, glass_names);
    const Real max_j_imag_phi = hostVecdFieldMax(particles, glass_names.j_imag, n);

    const Real div_j_reduction = div_j_level0.div_j_rel / (div_j_phi.div_j_rel + TinyReal);
    const bool passed = n > 0 && div_j_level0.div_j_rel > 0.1 &&
                        max_j_imag_phi < 0.01 &&
                        div_j_phi.div_j_weighted_l2 < 0.1 * div_j_level0.div_j_weighted_l2;

    std::cout << "test_3d_ophelie_phi_reduces_divJ n=" << n << " divJ_L0=" << div_j_level0.div_j_rel
              << " divJ_phi=" << div_j_phi.div_j_rel << " divJ_L2_L0=" << div_j_level0.div_j_weighted_l2
              << " divJ_L2_phi=" << div_j_phi.div_j_weighted_l2 << " divJ_reduction=" << div_j_reduction
              << " max_JImag_phi=" << max_j_imag_phi << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
