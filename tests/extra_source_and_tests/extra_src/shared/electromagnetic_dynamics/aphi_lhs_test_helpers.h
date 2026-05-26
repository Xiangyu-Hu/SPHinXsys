#ifndef APHI_LHS_TEST_HELPERS_H
#define APHI_LHS_TEST_HELPERS_H

#include "sphinxsys.h"
#include "electromagnetic_dynamics/all_electromagnetic_dynamics_ck.h"

#include <algorithm>
#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

using MainExecutionPolicy = execution::MainExecutionPolicy;

class AphiLhsBoxShape : public ComplexShape
{
  public:
    AphiLhsBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

inline Real distanceToBoundary(const Vecd &position, Real length, Real height, Real width)
{
    return std::min({position[0], length - position[0], position[1], height - position[1], position[2], width - position[2]});
}

inline bool isCoreParticle(const Vecd &position, Real length, Real height, Real width, Real core_shell)
{
    return distanceToBoundary(position, length, height, width) > core_shell;
}

struct AphiLhsTestBody
{
    SPHSystem sph_system;
    SolidBody body;
    Inner<> inner_ck;

    AphiLhsTestBody(Real dp_0, Real body_length, Real body_height, Real body_width, Real boundary_width)
        : sph_system(BoundingBoxd(Vecd(-boundary_width, -boundary_width, -boundary_width),
                                  Vecd(body_length + boundary_width, body_height + boundary_width,
                                       body_width + boundary_width)),
                     dp_0),
          body(sph_system,
               makeShared<AphiLhsBoxShape>("AphiLhsTestBody",
                                           Vecd(0.5 * body_length, 0.5 * body_height, 0.5 * body_width),
                                           Vecd(0.5 * body_length, 0.5 * body_height, 0.5 * body_width))),
          inner_ck(body)
    {
        body.defineAdaptation<SPHAdaptation>(1.15, 1.0);
        body.defineMaterial<Solid>();
        body.defineBodyLevelSetShape();
        body.generateParticles<BaseParticles, Lattice>();
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
    }

    void updateRelations()
    {
        UpdateCellLinkedList<MainExecutionPolicy, RealBody> update_cell_linked_list(body);
        UpdateRelation<MainExecutionPolicy, Inner<>> update_inner_relation(inner_ck);
        update_cell_linked_list.exec();
        update_inner_relation.exec();
    }
};

inline Real blockMaxAbsResidual(BaseParticles &particles, const AphiVariableNames &names, size_t total_real_particles)
{
    const Vecd *res_a_real = particles.getVariableDataByName<Vecd>(names.residual.a_real);
    const Vecd *res_a_imag = particles.getVariableDataByName<Vecd>(names.residual.a_imag);
    const Real *res_phi_real = particles.getVariableDataByName<Real>(names.residual.phi_real);
    const Real *res_phi_imag = particles.getVariableDataByName<Real>(names.residual.phi_imag);
    Real max_error = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        max_error = std::max(max_error, res_a_real[i].norm());
        max_error = std::max(max_error, res_a_imag[i].norm());
        max_error = std::max(max_error, std::abs(res_phi_real[i]));
        max_error = std::max(max_error, std::abs(res_phi_imag[i]));
    }
    return max_error;
}

inline Real coreBlockMaxAbsResidual(BaseParticles &particles, const AphiVariableNames &names, const Vecd *positions,
                                    size_t total_real_particles, Real body_length, Real body_height, Real body_width,
                                    Real core_shell)
{
    const Vecd *res_a_real = particles.getVariableDataByName<Vecd>(names.residual.a_real);
    const Vecd *res_a_imag = particles.getVariableDataByName<Vecd>(names.residual.a_imag);
    const Real *res_phi_real = particles.getVariableDataByName<Real>(names.residual.phi_real);
    const Real *res_phi_imag = particles.getVariableDataByName<Real>(names.residual.phi_imag);
    Real max_error = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        max_error = std::max(max_error, res_a_real[i].norm());
        max_error = std::max(max_error, res_a_imag[i].norm());
        max_error = std::max(max_error, std::abs(res_phi_real[i]));
        max_error = std::max(max_error, std::abs(res_phi_imag[i]));
    }
    return max_error;
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_LHS_TEST_HELPERS_H
