#include "particle_functors.h"
#include "update_cell_linked_list_sycl.hpp"

#include <gtest/gtest.h>
using namespace SPH;

StdVec<UnsignedInt> cell_size_list{3, 2, 3, 5, 0, 1, 3, 2, 5, 1, 0};
UnsignedInt list_size = cell_size_list.size();
StdVec<UnsignedInt> result(list_size, 0);
StdVec<UnsignedInt> sycl_result = result;

StdVec<UnsignedInt> input_result{0, 3, 5, 8, 13, 13, 14, 17, 19, 24, 25};

TEST(exclusive_scan, test_sycl)
{
    UnsignedInt *list_data = cell_size_list.data();
    exclusive_scan(SequencedPolicy{}, list_data, list_data + list_size, result.data(), std::plus<UnsignedInt>());

    UnsignedInt *device_list_data = allocateDeviceOnly<UnsignedInt>(list_size);
    UnsignedInt *device_result = allocateDeviceOnly<UnsignedInt>(list_size);
    copyToDevice(list_data, device_list_data, list_size);
    exclusive_scan(ParallelDevicePolicy{}, device_list_data, device_list_data + list_size, device_result, sycl::plus<UnsignedInt>());

    copyFromDevice(sycl_result.data(), device_result, list_size);
    freeDeviceData(device_list_data);
    freeDeviceData(device_result);

    EXPECT_EQ(result, sycl_result);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
