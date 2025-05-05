#include "device_copyable_variable.h"
#include "loop_range.h"
#include "particle_iterators.h"
#include "particle_iterators_sycl.h"
#include "reduce_functors.h"
#include "simtk_wrapper.h"
#include "sphinxsys_variable_sycl.hpp"

#include <gtest/gtest.h>

using namespace SPH;

TEST(particle_reduce, test_sycl)
{
    StdVec<SimTKVec3> torques(10000);
    StdVec<SimTKVec3> forces(10000);
    for (UnsignedInt i = 0; i != torques.size(); ++i)
    {
        torques[i] = SimTKVec3(rand_uniform(-1.0, 1.0),
                               rand_uniform(-1.0, 1.0),
                               rand_uniform(-1.0, 1.0));
        forces[i] = SimTKVec3(rand_uniform(-1.0, 1.0),
                              rand_uniform(-1.0, 1.0),
                              rand_uniform(-1.0, 1.0));
    }

    SimTK::SpatialVec sum = particle_reduce(SequencedPolicy{}, IndexRange(0, torques.size()),
                                            ZeroData<SimTK::SpatialVec>::value, ReduceSum<SimTK::SpatialVec>(),
                                            [&](size_t i)
                                            {
                                                SimTKVec3 a = SimTK::cross(torques[i], forces[i]);
                                                return SimTK::SpatialVec(a, forces[i]);
                                            });

    DiscreteVariable<SimTKVec3> dv_torque("Torque", torques.size());
    DiscreteVariable<SimTKVec3> dv_force("Force", forces.size());

    SimTKVec3 *torque = dv_torque.Data();
    SimTKVec3 *force = dv_force.Data();
    for (UnsignedInt i = 0; i != dv_torque.getDataSize(); ++i)
    {
        torque[i] = torques[i];
        force[i] = forces[i];
    }
    SingularVariable<UnsignedInt> sv_total_particles("TotalParticles", dv_torque.getDataSize());

    SimTKVec3 *torque_ck = dv_torque.DelegatedData(ParallelPolicy{});
    SimTKVec3 *force_ck = dv_force.DelegatedData(ParallelPolicy{});
    SimTK::SpatialVec sum_ck = particle_reduce<ReduceSum<SimTK::SpatialVec>>(
        LoopRangeCK<ParallelPolicy, SPHBody>(&sv_total_particles),
        ReduceReference<ReduceSum<SimTK::SpatialVec>>::value,
        [=](size_t i)
        {
            SimTKVec3 a = SimTK::cross(torque_ck[i], force_ck[i]);
            return SimTK::SpatialVec(a, force_ck[i]);
        });

    SimTKVec3 *torque_sycl = dv_torque.DelegatedData(ParallelDevicePolicy{});
    SimTKVec3 *force_sycl = dv_force.DelegatedData(ParallelDevicePolicy{});
    SimTK::SpatialVec sum_sycl = particle_reduce<ReduceSum<SimTK::SpatialVec>>(
        LoopRangeCK<ParallelDevicePolicy, SPHBody>(&sv_total_particles),
        ReduceReference<ReduceSum<SimTK::SpatialVec>>::value,
        [=](size_t i)
        {
            SimTKVec3 a = SimTK::cross(torque_sycl[i], force_sycl[i]);
            return SimTK::SpatialVec(a, force_sycl[i]);
        });
    EXPECT_EQ(sum, sum_ck);
    EXPECT_EQ(sum, sum_sycl);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
