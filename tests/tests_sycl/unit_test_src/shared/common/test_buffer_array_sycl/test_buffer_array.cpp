#include "loop_range.h"
#include "particle_iterators.h"
#include "particle_iterators_sycl.h"
#include "reduce_functors.h"
#include "simtk_wrapper.h"
#include "sphinxsys_buffer_array_sycl.hpp"
#include "sphinxsys_variable_sycl.hpp"

#include <gtest/gtest.h>

using namespace SPH;

TEST(buffer_array, test_sycl)
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

    DiscreteVariable<SimTKVec3> dv_torque("Torque", torques.size(), [&](size_t i)
                                          { return torques[i]; });
    DiscreteVariable<SimTKVec3> dv_force("Force", forces.size(), [&](size_t i)
                                         { return forces[i]; });

    SimTKVec3 *torque = dv_torque.DelegatedData(ParallelPolicy{});
    SimTKVec3 *force = dv_force.DelegatedData(ParallelPolicy{});
    SimTK::SpatialVec partial_sum =
        particle_reduce(ParallelPolicy{}, IndexRange(5000, 5050),
                        ZeroData<SimTK::SpatialVec>::value, ReduceSum<SimTK::SpatialVec>(),
                        [&](size_t i)
                        {
                            SimTKVec3 a = SimTK::cross(torque[i], force[i]);
                            return SimTK::SpatialVec(a, force[i]);
                        });

    StdVec<DiscreteVariable<SimTKVec3> *> variables = {&dv_torque, &dv_force};

    DiscreteVariableArray<SimTKVec3> variable_array(variables);
    BufferArray<SimTKVec3> buffer_array(variables, 100);

    AllocatedDataArray<SimTKVec3> variable_data_array =
        variable_array.DelegatedDataArray(ParallelPolicy{});
    AllocatedDataArray<SimTKVec3> buffer_data_array =
        buffer_array.DelegatedDataArray(ParallelPolicy{});

    AllocationDataArrayPair<SimTKVec3>
        variable_buffer_allocation_pair(variable_data_array, buffer_data_array);

    AllocationDataArrayPairSet<SimTKVec3>
        variable_buffer_allocation_pair_set(variable_buffer_allocation_pair, 2);

    CopyAllocationDataArrayPairSet copy_variable_to_buffer;

    particle_for(ParallelPolicy{}, IndexRange(5000, 5050),
                 [&](size_t i)
                 {
                     copy_variable_to_buffer(variable_buffer_allocation_pair_set, i, i - 5000);
                 });

    SimTKVec3 *buffer_torque = buffer_data_array[0];
    SimTKVec3 *buffer_force = buffer_data_array[1];
    SimTK::SpatialVec partial_sum_ck =
        particle_reduce(ParallelPolicy{}, IndexRange(0, 50),
                        ZeroData<SimTK::SpatialVec>::value, ReduceSum<SimTK::SpatialVec>(),
                        [&](size_t i)
                        {
                            SimTKVec3 a = SimTK::cross(buffer_torque[i], buffer_force[i]);
                            return SimTK::SpatialVec(a, buffer_force[i]);
                        });

    AllocatedDataArray<SimTKVec3> buffer_data_array_sycl =
        buffer_array.DelegatedDataArray(ParallelDevicePolicy{});

    AllocationDataArrayPair<SimTKVec3>
        variable_buffer_allocation_host(variable_data_array, buffer_data_array_sycl);

    AllocationDataArrayPairSet<SimTKVec3>
        variable_buffer_allocation_set_host(variable_buffer_allocation_host, 2);

    particle_for(ParallelPolicy{}, IndexRange(5000, 5050), // copy to buffer on host
                 [&](size_t i)
                 {
                     copy_variable_to_buffer(variable_buffer_allocation_set_host, i, i - 5000);
                 });

    SimTKVec3 *buffer_torque_sycl = buffer_data_array_sycl[0];
    SimTKVec3 *buffer_force_sycl = buffer_data_array_sycl[1];

    SingularVariable<UnsignedInt> sv_total_particles("TotalParticles", 50);
    SimTK::SpatialVec partial_sum_sycl =
        particle_reduce<ReduceSum<SimTK::SpatialVec>>( // summation on device
            LoopRangeCK<ParallelDevicePolicy, SPHBody>(&sv_total_particles),
            ReduceReference<ReduceSum<SimTK::SpatialVec>>::value,
            [=](size_t i)
            {
                SimTKVec3 a = SimTK::cross(buffer_torque_sycl[i], buffer_force_sycl[i]);
                return SimTK::SpatialVec(a, buffer_force_sycl[i]);
            });

    EXPECT_EQ(partial_sum, partial_sum_ck);
    EXPECT_EQ(partial_sum, partial_sum_sycl);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
