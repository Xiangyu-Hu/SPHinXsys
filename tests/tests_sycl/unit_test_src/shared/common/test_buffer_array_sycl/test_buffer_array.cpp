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
    DiscreteVariable<UnsignedInt> dv_copy_indexes(
        "BufferIndexes", 100,
        [&](size_t i)
        { return int(ceil(rand_uniform(0.0, 1.0) * int(torques.size()))); });

    SingularVariable<UnsignedInt> sv_buffer_particles("TotalBufferParticles", dv_copy_indexes.getDataSize());

    SimTKVec3 *torque = dv_torque.DelegatedData(ParallelPolicy{});
    SimTKVec3 *force = dv_force.DelegatedData(ParallelPolicy{});
    UnsignedInt *copy_indexes = dv_copy_indexes.DelegatedData(ParallelPolicy{});
    SimTK::SpatialVec partial_sum =
        particle_reduce(ParallelPolicy{}, IndexRange(0, sv_buffer_particles.getValue()),
                        ZeroData<SimTK::SpatialVec>::value, ReduceSum<SimTK::SpatialVec>(),
                        [&](size_t i)
                        {
                            UnsignedInt index = copy_indexes[i];
                            SimTKVec3 a = SimTK::cross(torque[index], force[index]);
                            return SimTK::SpatialVec(a, force[index]);
                        });

    StdVec<DiscreteVariable<SimTKVec3> *> variables = {&dv_torque, &dv_force};

    DiscreteVariableArray<SimTKVec3> variable_array(variables);
    BufferArray<SimTKVec3> buffer_array(variables, 200);

    AllocatedDataArray<SimTKVec3> variable_data_array =
        variable_array.DelegatedDataArray(ParallelPolicy{});
    AllocatedDataArray<SimTKVec3> buffer_data_array =
        buffer_array.DelegatedDataArray(ParallelPolicy{});

    AllocatedDataArrayPair<SimTKVec3> variable_buffer_allocation_pair(variable_data_array, buffer_data_array);

    SingularVariable<AllocatedDataArrayPairSet<SimTKVec3>>
        sv_variable_buffer_allocation_pair_set(
            "VariableBufferAllocationPairSet",
            AllocatedDataArrayPairSet<SimTKVec3>(variable_buffer_allocation_pair, 2));

    SingularVariable<CopyAllocatedDataArrayPairSet>
        sv_copy_variable_to_buffer("CopyVariableToBuffer", CopyAllocatedDataArrayPairSet());

    AllocatedDataArrayPairSet<SimTKVec3> *variable_buffer_allocation_pair_set =
        sv_variable_buffer_allocation_pair_set.DelegatedData(ParallelPolicy{});

    CopyAllocatedDataArrayPairSet *copy_variable_to_buffer =
        sv_copy_variable_to_buffer.DelegatedData(ParallelPolicy{});

    particle_for(LoopRangeCK<ParallelPolicy, SPHBody>(&sv_buffer_particles),
                 [=](size_t i)
                 {
                     UnsignedInt index = copy_indexes[i];
                     (*copy_variable_to_buffer)(*variable_buffer_allocation_pair_set, index, i);
                 });

    SimTKVec3 *buffer_torque = buffer_data_array[0];
    SimTKVec3 *buffer_force = buffer_data_array[1];
    SimTK::SpatialVec partial_sum_ck = particle_reduce<ReduceSum<SimTK::SpatialVec>>(
        LoopRangeCK<ParallelPolicy, SPHBody>(&sv_buffer_particles),
        ReduceReference<ReduceSum<SimTK::SpatialVec>>::value,
        [=](size_t i)
        {
            SimTKVec3 a = SimTK::cross(buffer_torque[i], buffer_force[i]);
            return SimTK::SpatialVec(a, buffer_force[i]);
        });

    AllocatedDataArray<SimTKVec3> variable_data_array_sycl =
        variable_array.DelegatedDataArray(ParallelDevicePolicy{});
    AllocatedDataArray<SimTKVec3> buffer_data_array_sycl =
        buffer_array.DelegatedDataArray(ParallelDevicePolicy{});

    AllocatedDataArrayPair<SimTKVec3> variable_buffer_allocation_pair_sycl(variable_data_array_sycl, buffer_data_array_sycl);

    SingularVariable<AllocatedDataArrayPairSet<SimTKVec3>>
        sv_variable_buffer_allocation_pair_set_sycl(
            "VariableBufferAllocationPairSet",
            AllocatedDataArrayPairSet<SimTKVec3>(variable_buffer_allocation_pair_sycl, 2));

    AllocatedDataArrayPairSet<SimTKVec3> *variable_buffer_allocation_pair_set_sycl =
        sv_variable_buffer_allocation_pair_set_sycl.DelegatedData(ParallelDevicePolicy{});

    CopyAllocatedDataArrayPairSet *copy_variable_to_buffer_sycl =
        sv_copy_variable_to_buffer.DelegatedData(ParallelDevicePolicy{});

    UnsignedInt *copy_indexes_sycl = dv_copy_indexes.DelegatedData(ParallelDevicePolicy{});
    particle_for(LoopRangeCK<ParallelDevicePolicy, SPHBody>(&sv_buffer_particles),
                 [=](size_t i)
                 {
                     UnsignedInt index = copy_indexes_sycl[i];
                     (*copy_variable_to_buffer_sycl)(*variable_buffer_allocation_pair_set_sycl, index, i);
                 });

    SimTKVec3 *buffer_torque_sycl = buffer_data_array_sycl[0];
    SimTKVec3 *buffer_force_sycl = buffer_data_array_sycl[1];

    SimTK::SpatialVec partial_sum_sycl =
        particle_reduce<ReduceSum<SimTK::SpatialVec>>( // summation on device
            LoopRangeCK<ParallelDevicePolicy, SPHBody>(&sv_buffer_particles),
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
