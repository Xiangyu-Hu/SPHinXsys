#include "loop_range.h"
#include "particle_iterators.h"
#include "particle_iterators_sycl.h"
#include "reduce_functors.h"
#include "simtk_wrapper.h"
#include "sphinxsys_buffer_array_sycl.hpp"

#include <gtest/gtest.h>

using namespace SPH;

struct CopyVariableToBufferArray
{
    template <typename DataType>
    void operator()(VariableBufferAllocationPair<AllocatedDataArray<DataType>> &variable_buffer_allocation_pair,
                    size_t variable_data_index, size_t buffer_data_index)
    {
        auto &variable_allocation = variable_buffer_allocation_pair.first.first;
        auto &buffer_allocation = variable_buffer_allocation_pair.first.second;
        for (size_t i = 0; i < variable_buffer_allocation_pair.second; ++i)
        {
            buffer_allocation[i][buffer_data_index] = variable_allocation[i][variable_data_index];
        }
    }
};

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

    VariableBufferAllocation<AllocatedDataArray<SimTKVec3>>
        variable_buffer_allocation(variable_data_array, buffer_data_array);

    VariableBufferAllocationPair<AllocatedDataArray<SimTKVec3>>
        variable_buffer_allocation_pair(variable_buffer_allocation, 2);

    CopyVariableToBufferArray copy_variable_to_buffer;

    particle_for(ParallelPolicy{}, IndexRange(5000, 5050),
                 [&](size_t i)
                 {
                     copy_variable_to_buffer(variable_buffer_allocation_pair, i, i - 5000);
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

    EXPECT_EQ(partial_sum, partial_sum_ck);
    // EXPECT_EQ(partial_sum, partial_sum_sycl);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
