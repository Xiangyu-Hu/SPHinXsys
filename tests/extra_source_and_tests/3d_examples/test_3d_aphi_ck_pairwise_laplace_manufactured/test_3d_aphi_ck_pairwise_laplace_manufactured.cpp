#include "sphinxsys.h"
#include "electromagnetic_dynamics/all_electromagnetic_dynamics_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;

namespace
{

using MainExecutionPolicy = execution::MainExecutionPolicy;

class PairwiseLaplaceBoxShape : public ComplexShape
{
  public:
    PairwiseLaplaceBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

class AssignQuadraticAphiFieldsCK : public LocalDynamics
{
  public:
    explicit AssignQuadraticAphiFieldsCK(SPHBody &sph_body, const AphiVariableNames &names)
        : LocalDynamics(sph_body),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_phi_real_(particles_->template getVariableByName<Real>(names.solution.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(names.solution.phi_imag)),
          dv_a_real_(particles_->template getVariableByName<Vecd>(names.solution.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(names.solution.a_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            const Real x = position[0];
            const Real y = position[1];
            const Real z = position[2];
            const Real quadratic_sum = x * x + y * y + z * z;

            phi_real_[index_i] = quadratic_sum;
            phi_imag_[index_i] = quadratic_sum;
            a_real_[index_i] = Vecd(x * x, y * y, z * z);
            a_imag_[index_i] = Vecd(z * z, x * x, y * y);
        }

      protected:
        Vecd *position_;
        Real *phi_real_;
        Real *phi_imag_;
        Vecd *a_real_;
        Vecd *a_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
};

Real distanceToBoundary(const Vecd &position, Real length, Real height, Real width)
{
    return std::min({position[0], length - position[0], position[1], height - position[1], position[2], width - position[2]});
}

struct Summary
{
    size_t total_particles = 0;
    size_t core_particles = 0;
    Real scalar_laplace_core_max_error = 0.0;
    Real vector_laplace_core_max_error = 0.0;
};

} // namespace

int main(int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.5 * dp_0;
    const Real validation_threshold = 0.25;

    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width, body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);

    const Vecd center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    SolidBody body(sph_system, makeShared<PairwiseLaplaceBoxShape>("AphiPairwiseLaplaceBody", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    Inner<> inner_ck(body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(body, 1.0, 1.0, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_material_properties(body, 1.0, 1.0, names.material);
    StateDynamics<MainExecutionPolicy, AssignQuadraticAphiFieldsCK> assign_fields(body, names);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner_ck);

    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Real>>> laplace_phi_real(
        DynamicsArgs(inner_ck, names.solution.phi_real, names.material.sigma, names.diagnostic.laplace.phi_real));
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Real>>> laplace_phi_imag(
        DynamicsArgs(inner_ck, names.solution.phi_imag, names.material.sigma, names.diagnostic.laplace.phi_imag));
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Vecd>>> laplace_a_real(
        DynamicsArgs(inner_ck, names.solution.a_real, names.material.nu, names.diagnostic.laplace.a_real));
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Vecd>>> laplace_a_imag(
        DynamicsArgs(inner_ck, names.solution.a_imag, names.material.nu, names.diagnostic.laplace.a_imag));

    initialize_aphi_variables.exec();
    set_material_properties.exec();
    assign_fields.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    laplace_phi_real.exec();
    laplace_phi_imag.exec();
    laplace_a_real.exec();
    laplace_a_imag.exec();

    BaseParticles &particles = body.getBaseParticles();
    syncVariableToHost<Real>(particles, names.diagnostic.laplace.phi_real);
    syncVariableToHost<Real>(particles, names.diagnostic.laplace.phi_imag);
    syncVariableToHost<Vecd>(particles, names.diagnostic.laplace.a_real);
    syncVariableToHost<Vecd>(particles, names.diagnostic.laplace.a_imag);

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *lap_phi_real = particles.getVariableDataByName<Real>(names.diagnostic.laplace.phi_real);
    const Real *lap_phi_imag = particles.getVariableDataByName<Real>(names.diagnostic.laplace.phi_imag);
    const Vecd *lap_a_real = particles.getVariableDataByName<Vecd>(names.diagnostic.laplace.a_real);
    const Vecd *lap_a_imag = particles.getVariableDataByName<Vecd>(names.diagnostic.laplace.a_imag);

    const Real scalar_exact = -6.0;
    const Vecd vector_exact(-2.0, -2.0, -2.0);

    Summary summary;
    summary.total_particles = total_real_particles;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (distanceToBoundary(positions[i], body_length, body_height, body_width) <= core_shell)
        {
            continue;
        }

        summary.core_particles += 1;
        summary.scalar_laplace_core_max_error =
            std::max(summary.scalar_laplace_core_max_error, std::abs(lap_phi_real[i] - scalar_exact));
        summary.scalar_laplace_core_max_error =
            std::max(summary.scalar_laplace_core_max_error, std::abs(lap_phi_imag[i] - scalar_exact));
        summary.vector_laplace_core_max_error =
            std::max(summary.vector_laplace_core_max_error, (lap_a_real[i] - vector_exact).norm());
        summary.vector_laplace_core_max_error =
            std::max(summary.vector_laplace_core_max_error, (lap_a_imag[i] - vector_exact).norm());
    }

    const bool passed = summary.core_particles > 0 &&
                        summary.scalar_laplace_core_max_error < validation_threshold &&
                        summary.vector_laplace_core_max_error < validation_threshold;

    std::cout << "test_3d_aphi_ck_pairwise_laplace_manufactured"
              << " total_particles=" << summary.total_particles
              << " core_particles=" << summary.core_particles
              << " scalar_laplace_core_max_error=" << summary.scalar_laplace_core_max_error
              << " vector_laplace_core_max_error=" << summary.vector_laplace_core_max_error
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
