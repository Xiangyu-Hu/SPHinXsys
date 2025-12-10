#include "sphinxsys.h"
#include <gtest/gtest.h>

using namespace SPH;

TEST(test_meshes, split_for)
{
    Real length = 10;
    Real dp = 1;

    auto shape = makeShared<GeometricShapeBox>(Transform(), 0.5 * length * Vec3d::Ones(), "Shape");

    BoundingBoxd bb_system = shape->getBounds();

    SPHSystem system(bb_system, dp);

    SolidBody body(system, shape);
    body.defineMaterial<Solid>();
    body.generateParticles<BaseParticles, Lattice>();
    auto &particles = body.getBaseParticles();
    const auto pos = particles.registerStateVariableData<Vec3d>("Position");
    auto quantity = particles.registerStateVariableData<Real>("Quantity", [&](size_t i) -> Real
                                                              { return pos[i].norm(); });

    InnerRelation inner(body);

    body.updateCellLinkedList();
    inner.updateConfiguration();

    auto interaction = [&](size_t index_i)
    {
        const auto &configuration = inner.inner_configuration_;
        const auto &neighborhood = configuration[index_i];
        for (size_t n = 0; n < neighborhood.current_size_; n++)
        {
            size_t q_ij = quantity[index_i] - quantity[neighborhood.j_[n]];
            quantity[index_i] += 0.5 * q_ij;
            quantity[neighborhood.j_[n]] -= 0.5 * q_ij;
        }
    };

    auto &cell_linked_list = *dynamic_cast<CellLinkedList *>(&body.getCellLinkedList());

    // run the interaction in sequenced policy
    cell_linked_list.particle_for_split(execution::SequencedPolicy(), interaction);
    // record the result
    auto q_seq = quantity;

    // reset the value
    for (size_t i = 0; i < particles.TotalRealParticles(); i++)
    {
        quantity[i] = pos[i].norm();
    }
    // run the interaction in parallel policy
    cell_linked_list.particle_for_split(execution::ParallelPolicy(), interaction);

    // check the result
    for (size_t i = 0; i < particles.TotalRealParticles(); i++)
    {
        ASSERT_EQ(q_seq[i], quantity[i]);
    }
}