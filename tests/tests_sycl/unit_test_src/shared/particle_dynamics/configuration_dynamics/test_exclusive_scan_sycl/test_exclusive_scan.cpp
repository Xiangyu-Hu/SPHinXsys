#include "base_configuration_dynamics.h"
#include "base_configuration_dynamics_sycl.h"
#include "particle_functors.h"
#include "particle_iterators_ck.h"
#include "particle_iterators_sycl.h"

#include <gtest/gtest.h>
using namespace SPH;

StdVec<UnsignedInt> cell_size_list{3, 2, 3, 5, 0, 1, 3, 2, 5, 1, 11}; //the last entry of the input is not used
StdVec<UnsignedInt> solution{0, 3, 5, 8, 13, 13, 14, 17, 19, 24, 25};

UnsignedInt list_size = cell_size_list.size();
StdVec<UnsignedInt> result(list_size, 0);
StdVec<UnsignedInt> sycl_result = result;

TEST(exclusive_scan, test_sycl)
{
    UnsignedInt *list_data = cell_size_list.data();
    UnsignedInt sum = exclusive_scan(SequencedPolicy{}, list_data, result.data(), list_size, PlusUnsignedInt<SequencedPolicy>::type());

    UnsignedInt *device_list_data = allocateDeviceOnly<UnsignedInt>(list_size);
    UnsignedInt *device_result = allocateDeviceOnly<UnsignedInt>(list_size);
    copyToDevice(list_data, device_list_data, list_size);
    UnsignedInt sycl_sum = exclusive_scan(ParallelDevicePolicy{}, device_list_data, device_result, list_size, PlusUnsignedInt<ParallelDevicePolicy>::type());
    copyFromDevice(sycl_result.data(), device_result, list_size);

    freeDeviceData(device_list_data);
    freeDeviceData(device_result);

    EXPECT_EQ(solution, result);
    EXPECT_EQ(solution, sycl_result);
    EXPECT_EQ(sum, sycl_sum);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
