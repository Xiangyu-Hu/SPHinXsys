/**
 * Stage 2 diagnostic test: LinearGradient + LinearCorrectionMatrix.
 * DIAGNOSTIC_B_CORRECTED / DEBUG_REFERENCE_ONLY — not the matrix-free OPERATOR default path.
 */
#include "sphinxsys.h"
#include "electromagnetic_dynamics/all_electromagnetic_dynamics_ck.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>

using namespace SPH;
using namespace SPH::electromagnetics;

namespace
{

using MainExecutionPolicy = execution::MainExecutionPolicy;

class GradientDivergenceBoxShape : public ComplexShape
{
  public:
    GradientDivergenceBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize)
        : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

class AssignManufacturedAphiFieldsCK : public LocalDynamics
{
  public:
    explicit AssignManufacturedAphiFieldsCK(SPHBody &sph_body, const AphiVariableNames &names)
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

            phi_real_[index_i] = x + y + z;
            phi_imag_[index_i] = 0.5 * x - 0.25 * y + 0.75 * z;
            a_real_[index_i] = Vecd(x, y, z);
            a_imag_[index_i] = Vecd(y, z, x);
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
    Real phi_real_gradient_core_max_error = 0.0;
    Real phi_imag_gradient_core_max_error = 0.0;
    Real div_a_real_core_max_error = 0.0;
    Real div_a_imag_core_max_error = 0.0;
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
    SolidBody body(sph_system, makeShared<GradientDivergenceBoxShape>("AphiGradientDivergenceBody", center, halfsize));
    body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    body.defineMaterial<Solid>();
    body.defineBodyLevelSetShape();
    body.generateParticles<BaseParticles, Lattice>();

    InnerRelation inner_relation(body);
    Inner<> inner_ck(body);

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_aphi_variables(body, 1.0, 1.0, names);
    StateDynamics<MainExecutionPolicy, AssignManufacturedAphiFieldsCK> assign_fields(body, names);

    initialize_aphi_variables.exec();
    assign_fields.exec();

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner_ck);

    // DIAGNOSTIC_B_CORRECTED / DEBUG_REFERENCE_ONLY — not the matrix-free OPERATOR default path.
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrix<Inner<WithUpdate>>> linear_correction_matrix(
        DynamicsArgs(inner_ck, 0.0));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>> phi_real_gradient(
        DynamicsArgs(inner_ck, names.solution.phi_real));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Real>>> phi_imag_gradient(
        DynamicsArgs(inner_ck, names.solution.phi_imag));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_real_gradient(
        DynamicsArgs(inner_ck, names.solution.a_real));
    InteractionDynamicsCK<MainExecutionPolicy, LinearGradient<Inner<Vecd>>> a_imag_gradient(
        DynamicsArgs(inner_ck, names.solution.a_imag));
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_real(body,
                                                                                  names.solution.a_real + "Gradient",
                                                                                  names.diagnostic.div_a_real);
    StateDynamics<MainExecutionPolicy, AphiVectorGradientDivergenceCK> div_a_imag(body,
                                                                                  names.solution.a_imag + "Gradient",
                                                                                  names.diagnostic.div_a_imag);

    update_cell_linked_list.exec();
    update_inner_relation.exec();
    linear_correction_matrix.exec();
    phi_real_gradient.exec();
    phi_imag_gradient.exec();
    a_real_gradient.exec();
    a_imag_gradient.exec();
    div_a_real.exec();
    div_a_imag.exec();

    BaseParticles &particles = body.getBaseParticles();
    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *phi_real_gradient_data = particles.getVariableDataByName<Vecd>(names.solution.phi_real + "Gradient");
    const Vecd *phi_imag_gradient_data = particles.getVariableDataByName<Vecd>(names.solution.phi_imag + "Gradient");
    const Real *div_a_real_data = particles.getVariableDataByName<Real>(names.diagnostic.div_a_real);
    const Real *div_a_imag_data = particles.getVariableDataByName<Real>(names.diagnostic.div_a_imag);

    const Vecd phi_real_gradient_exact(1.0, 1.0, 1.0);
    const Vecd phi_imag_gradient_exact(0.5, -0.25, 0.75);
    const Real div_a_real_exact = 3.0;
    const Real div_a_imag_exact = 0.0;

    Summary summary;
    summary.total_particles = total_real_particles;

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (distanceToBoundary(positions[i], body_length, body_height, body_width) <= core_shell)
        {
            continue;
        }

        summary.core_particles += 1;
        summary.phi_real_gradient_core_max_error =
            std::max(summary.phi_real_gradient_core_max_error, (phi_real_gradient_data[i] - phi_real_gradient_exact).norm());
        summary.phi_imag_gradient_core_max_error =
            std::max(summary.phi_imag_gradient_core_max_error, (phi_imag_gradient_data[i] - phi_imag_gradient_exact).norm());
        summary.div_a_real_core_max_error =
            std::max(summary.div_a_real_core_max_error, std::abs(div_a_real_data[i] - div_a_real_exact));
        summary.div_a_imag_core_max_error =
            std::max(summary.div_a_imag_core_max_error, std::abs(div_a_imag_data[i] - div_a_imag_exact));
    }

    const bool passed = summary.core_particles > 0 &&
                        summary.phi_real_gradient_core_max_error < validation_threshold &&
                        summary.phi_imag_gradient_core_max_error < validation_threshold &&
                        summary.div_a_real_core_max_error < validation_threshold &&
                        summary.div_a_imag_core_max_error < validation_threshold;

    std::cout << "test_3d_aphi_ck_gradient_divergence_manufactured"
              << " total_particles=" << summary.total_particles
              << " core_particles=" << summary.core_particles
              << " phi_real_gradient_core_max_error=" << summary.phi_real_gradient_core_max_error
              << " phi_imag_gradient_core_max_error=" << summary.phi_imag_gradient_core_max_error
              << " div_a_real_core_max_error=" << summary.div_a_real_core_max_error
              << " div_a_imag_core_max_error=" << summary.div_a_imag_core_max_error
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
