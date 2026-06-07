/**
 * Stage 10-contact-1: monolithic Inner vs split two-body Inner+Contact pairwise Laplace equivalence.
 * Apply-level only (no GMRES).
 */
#include "sphinxsys.h"
#include "electromagnetic_dynamics/all_electromagnetic_dynamics_ck.h"
#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

using namespace SPH;
using namespace SPH::electromagnetics;
using namespace SPH::electromagnetics::test;
using namespace SPH::electromagnetics::benchmark;

namespace
{

using MainExecutionPolicy = execution::MainExecutionPolicy;

class AphiBoxShape : public ComplexShape
{
  public:
    AphiBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

class AssignQuadraticPhiRealCK : public LocalDynamics
{
  public:
    explicit AssignQuadraticPhiRealCK(SPHBody &sph_body, const std::string &phi_real_name)
        : LocalDynamics(sph_body),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_phi_real_(particles_->template getVariableByName<Real>(phi_real_name))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            phi_real_[index_i] = position.squaredNorm();
        }

      protected:
        Vecd *position_;
        Real *phi_real_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Real> *dv_phi_real_;
};

class AssignConstantSigmaCK : public LocalDynamics
{
  public:
    AssignConstantSigmaCK(SPHBody &sph_body, Real sigma, const AphiMaterialNames &material_names)
        : LocalDynamics(sph_body),
          dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma)),
          sigma_(sigma)
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), sigma_value_(encloser.sigma_)
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            sigma_[index_i] = sigma_value_;
        }

      protected:
        Real *sigma_;
        Real sigma_value_;
    };

  protected:
    DiscreteVariable<Real> *dv_sigma_;
    Real sigma_;
};

struct PositionKey
{
    Vecd position;
    bool operator<(const PositionKey &other) const
    {
        for (int d = 0; d < 3; ++d)
        {
            if (position[d] + TinyReal < other.position[d])
                return true;
            if (position[d] > other.position[d] + TinyReal)
                return false;
        }
        return false;
    }
};

using LaplaceByPosition = std::map<PositionKey, Real>;

inline PositionKey makePositionKey(const Vecd &position)
{
    PositionKey key;
    key.position = position;
    return key;
}

inline void runMonolithicLaplace(LaplaceByPosition &laplace_by_position, int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.5 * dp_0;
    const Real x_interface = 0.5;
    const Real sigma_left = 10.0;
    const Real sigma_right = 1.0;

    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width,
                                    body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);

    const Vecd center(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd halfsize(0.5 * body_length, 0.5 * body_height, 0.5 * body_width);
    SolidBody body(sph_system, makeShared<AphiBoxShape>("MonolithicBody", center, halfsize));
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
    StateDynamics<MainExecutionPolicy, AssignPiecewiseSigmaHalfSpaceCK> assign_sigma(
        body, x_interface, sigma_left, sigma_right, 1.0, names.material);
    StateDynamics<MainExecutionPolicy, AssignQuadraticPhiRealCK> assign_phi(body, names.solution.phi_real);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner_ck);

    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Real>>> laplace_phi_real(
        DynamicsArgs(inner_ck, names.solution.phi_real, names.material.sigma, names.diagnostic.laplace.phi_real));

    initialize_aphi_variables.exec();
    set_material_properties.exec();
    assign_sigma.exec();
    assign_phi.exec();
    update_cell_linked_list.exec();
    update_inner_relation.exec();
    laplace_phi_real.exec();

    BaseParticles &particles = body.getBaseParticles();
    syncVariableToHost<Real>(particles, names.diagnostic.laplace.phi_real);
    syncVariableToHost<Vecd>(particles, "Position");

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *laplace = particles.getVariableDataByName<Real>(names.diagnostic.laplace.phi_real);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (std::abs(positions[i][0] - x_interface) <= dp_0 + TinyReal)
        {
            continue;
        }
        laplace_by_position[makePositionKey(positions[i])] = laplace[i];
    }
}

inline void runSplitContactLaplace(LaplaceByPosition &laplace_by_position, int ac, char *av[])
{
    const Real dp_0 = 0.1;
    const Real body_length = 1.0;
    const Real body_height = 1.0;
    const Real body_width = 1.0;
    const Real boundary_width = 3.0 * dp_0;
    const Real core_shell = 2.5 * dp_0;
    const Real x_interface = 0.5;
    const Real sigma_left = 10.0;
    const Real sigma_right = 1.0;

    BoundingBoxd system_bounds(Vecd(-boundary_width, -boundary_width, -boundary_width),
                               Vecd(body_length + boundary_width, body_height + boundary_width,
                                    body_width + boundary_width));

    SPHSystem sph_system(system_bounds, dp_0);
    sph_system.handleCommandlineOptions(ac, av);

    const Vecd left_center(0.25 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd right_center(0.75 * body_length, 0.5 * body_height, 0.5 * body_width);
    const Vecd halfsize(0.25 * body_length, 0.5 * body_height, 0.5 * body_width);

    SolidBody left_body(sph_system, makeShared<AphiBoxShape>("LeftBody", left_center, halfsize));
    SolidBody right_body(sph_system, makeShared<AphiBoxShape>("RightBody", right_center, halfsize));
    for (auto *body_ptr : {&left_body, &right_body})
    {
        body_ptr->defineAdaptation<SPHAdaptation>(1.15, 1.0);
        body_ptr->defineMaterial<Solid>();
        body_ptr->defineBodyLevelSetShape();
        body_ptr->generateParticles<BaseParticles, Lattice>();
    }

    Inner<> left_inner(left_body);
    Inner<> right_inner(right_body);
    Contact<> left_to_right(left_body, {&right_body});
    Contact<> right_to_left(right_body, {&left_body});

    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();

    AphiVariableNames names;
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_left(left_body, 1.0, 1.0, names);
    StateDynamics<MainExecutionPolicy, InitializeAphiVariablesCK> initialize_right(right_body, 1.0, 1.0, names);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_left_material(left_body, 1.0, 1.0, names.material);
    StateDynamics<MainExecutionPolicy, SetAphiMaterialPropertiesCK> set_right_material(right_body, 1.0, 1.0, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantSigmaCK> assign_left_sigma(left_body, sigma_left, names.material);
    StateDynamics<MainExecutionPolicy, AssignConstantSigmaCK> assign_right_sigma(right_body, sigma_right, names.material);
    StateDynamics<MainExecutionPolicy, AssignQuadraticPhiRealCK> assign_left_phi(left_body, names.solution.phi_real);
    StateDynamics<MainExecutionPolicy, AssignQuadraticPhiRealCK> assign_right_phi(right_body, names.solution.phi_real);

    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_left_cell_linked_list(left_body);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_right_cell_linked_list(right_body);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_left_inner(left_inner);
    UpdateRelation<MainExecutionPolicy, Inner<>> update_right_inner(right_inner);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_left_contact(left_to_right);
    UpdateRelation<MainExecutionPolicy, Contact<>> update_right_contact(right_to_left);

    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Real>, Contact<Real>>> left_laplace(
        DynamicsArgs(left_inner, names.solution.phi_real, names.material.sigma, names.diagnostic.laplace.phi_real),
        DynamicsArgs(left_to_right, names.solution.phi_real, names.material.sigma, names.diagnostic.laplace.phi_real));
    InteractionDynamicsCK<MainExecutionPolicy, AphiPairwiseLaplaceCK<Inner<Real>, Contact<Real>>> right_laplace(
        DynamicsArgs(right_inner, names.solution.phi_real, names.material.sigma, names.diagnostic.laplace.phi_real),
        DynamicsArgs(right_to_left, names.solution.phi_real, names.material.sigma, names.diagnostic.laplace.phi_real));

    initialize_left.exec();
    initialize_right.exec();
    set_left_material.exec();
    set_right_material.exec();
    assign_left_sigma.exec();
    assign_right_sigma.exec();
    assign_left_phi.exec();
    assign_right_phi.exec();
    update_left_cell_linked_list.exec();
    update_right_cell_linked_list.exec();
    update_left_inner.exec();
    update_right_inner.exec();
    update_left_contact.exec();
    update_right_contact.exec();
    left_laplace.exec();
    right_laplace.exec();

    for (auto *body_ptr : {&left_body, &right_body})
    {
        BaseParticles &particles = body_ptr->getBaseParticles();
        syncVariableToHost<Real>(particles, names.diagnostic.laplace.phi_real);
        syncVariableToHost<Vecd>(particles, "Position");

        const size_t total_real_particles = particles.TotalRealParticles();
        const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
        const Real *laplace = particles.getVariableDataByName<Real>(names.diagnostic.laplace.phi_real);

        for (size_t i = 0; i != total_real_particles; ++i)
        {
            if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
            {
                continue;
            }
            if (std::abs(positions[i][0] - x_interface) <= dp_0 + TinyReal)
            {
                continue;
            }
            laplace_by_position[makePositionKey(positions[i])] = laplace[i];
        }
    }
}

struct Summary
{
    size_t matched_particles = 0;
    size_t missing_particles = 0;
    Real max_abs_diff = 0.0;
};

} // namespace

int main(int ac, char *av[])
{
    LaplaceByPosition monolithic_laplace;
    LaplaceByPosition contact_laplace;

    runMonolithicLaplace(monolithic_laplace, ac, av);
    runSplitContactLaplace(contact_laplace, ac, av);

    Summary summary;
    // float build: discrete pairwise accumulation tolerance; interface band excluded above
    const Real validation_threshold = 3.0e-5;

    for (const auto &entry : monolithic_laplace)
    {
        const auto it = contact_laplace.find(entry.first);
        if (it == contact_laplace.end())
        {
            summary.missing_particles += 1;
            continue;
        }
        summary.matched_particles += 1;
        summary.max_abs_diff = std::max(summary.max_abs_diff, std::abs(entry.second - it->second));
    }

    const bool passed = summary.missing_particles == 0 && summary.matched_particles > 0 &&
                        summary.max_abs_diff < validation_threshold;

    std::cout << "test_3d_aphi_ck_contact_pairwise_laplace_equivalence"
              << " matched_particles=" << summary.matched_particles
              << " missing_particles=" << summary.missing_particles
              << " max_abs_diff=" << summary.max_abs_diff
              << " passed=" << (passed ? 1 : 0) << std::endl;

    return passed ? 0 : 1;
}
