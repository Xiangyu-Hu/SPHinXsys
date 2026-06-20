/**
 * @file test_3d_ophelie_phi_gradient_linear.cpp
 * @brief PhiImag = x on a box lattice; ComputeOphelieScalarPhiGradientCK should yield GradPhi_x ≈ +1
 *        on interior particles (consistent with divJ using (f_j - f_i)·g_ij).
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

struct InteriorGradStats
{
    size_t n_interior = 0;
    Real vol_sum = 0.0;
    Real mean_grad_x = 0.0;
    Real mean_grad_y = 0.0;
    Real mean_grad_z = 0.0;
    Real max_abs_grad_yz = 0.0;
};

inline InteriorGradStats computeInteriorGradStats(BaseParticles &particles, const OphelieGlassFieldNames &names,
                                                  const Vecd &box_center, const Vecd &box_halfsize, Real margin)
{
    syncVariableToHost<Vecd>(particles, "Position");
    syncVariableToHost<Vecd>(particles, names.grad_phi_imag);
    syncVariableToHost<Real>(particles, "VolumetricMeasure");

    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *grad_phi = particles.getVariableDataByName<Vecd>(names.grad_phi_imag);
    const Real *vol = particles.getVariableDataByName<Real>("VolumetricMeasure");
    const size_t n = particles.TotalRealParticles();

    InteriorGradStats stats;
    for (size_t i = 0; i < n; ++i)
    {
        const Vecd r = pos[i] - box_center;
        if (std::abs(r[0]) > box_halfsize[0] - margin || std::abs(r[1]) > box_halfsize[1] - margin ||
            std::abs(r[2]) > box_halfsize[2] - margin)
        {
            continue;
        }
        stats.n_interior += 1;
        stats.vol_sum += vol[i];
        stats.mean_grad_x += vol[i] * grad_phi[i][0];
        stats.mean_grad_y += vol[i] * grad_phi[i][1];
        stats.mean_grad_z += vol[i] * grad_phi[i][2];
        stats.max_abs_grad_yz =
            std::max(stats.max_abs_grad_yz, std::max(std::abs(grad_phi[i][1]), std::abs(grad_phi[i][2])));
    }
    if (stats.vol_sum > TinyReal)
    {
        stats.mean_grad_x /= stats.vol_sum;
        stats.mean_grad_y /= stats.vol_sum;
        stats.mean_grad_z /= stats.vol_sum;
    }
    return stats;
}

} // namespace

int main(int, char *[])
{
    const Real dp = 0.08;
    const Vecd center(0.0, 0.0, 0.5);
    const Vecd halfsize(0.24, 0.24, 0.24);
    const Real interior_margin = 2.0 * dp;
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

    syncVariableToHost<Vecd>(particles, "Position");
    const Vecd *pos = particles.getVariableDataByName<Vecd>("Position");
    Real *phi_imag = particles.getVariableDataByName<Real>(glass_names.phi_imag);
    for (size_t i = 0; i < n; ++i)
    {
        phi_imag[i] = pos[i][0];
    }
    syncVariableToDevice<Real>(particles, glass_names.phi_imag);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(glass_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(*glass_inner);
    InteractionDynamicsCK<MainExecutionPolicy, ComputeOphelieScalarPhiGradientCK<Inner<>>> compute_grad_phi(
        *glass_inner, glass_names);
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    compute_grad_phi.exec();

    const InteriorGradStats stats =
        computeInteriorGradStats(particles, glass_names, center, halfsize, interior_margin);
    const Real expected_grad_x = 1.0;
    const Real rel_err_x = std::abs(stats.mean_grad_x - expected_grad_x) / expected_grad_x;
    const bool sign_ok = stats.mean_grad_x > 0.0;
    const bool passed = stats.n_interior > 0 && sign_ok && rel_err_x < 0.10 && stats.max_abs_grad_yz < 0.15;

    std::cout << "test_3d_ophelie_phi_gradient_linear n=" << n << " n_interior=" << stats.n_interior
              << " margin=" << interior_margin << " mean_GradPhi_x=" << stats.mean_grad_x
              << " mean_GradPhi_y=" << stats.mean_grad_y << " mean_GradPhi_z=" << stats.mean_grad_z
              << " rel_err_x=" << rel_err_x << " sign_ok=" << (sign_ok ? 1 : 0)
              << " max_abs_grad_yz=" << stats.max_abs_grad_yz;
    if (!sign_ok)
    {
        std::cout << " NOTE=GradPhi_x_negative_expect_flip_(phi_j-phi_i)";
    }
    std::cout << " passed=" << (passed ? 1 : 0) << std::endl;
    return passed ? 0 : 1;
}
