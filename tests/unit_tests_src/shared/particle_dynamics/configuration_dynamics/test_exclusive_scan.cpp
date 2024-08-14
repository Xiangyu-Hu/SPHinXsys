#include "particle_functors.h"
#include "update_cell_linked_list_sycl.hpp"

#include <gtest/gtest.h>
using namespace SPH;

StdVec<UnsignedInt> cell_size_list{3, 2, 3, 5, 0, 1, 3, 2, 5, 1, 0};
StdVec<UnsignedInt> result(cell_size_list.size(), 0);
StdVec<UnsignedInt> tbb_result = result;

StdVec<UnsignedInt> input_result{0, 3, 5, 8, 13, 13, 14, 17, 19, 24, 25};

TEST(exclusive_scan, test_tbb_sycl)
{
    UnsignedInt *list_tbb = cell_size_list.data();

    std::exclusive_scan(list_tbb, list_tbb + 11, result.data(), 0, ReduceSum<UnsignedInt>());
    SPH::exclusive_scan(list_tbb, list_tbb + 11, tbb_result.data(), ReduceSum<UnsignedInt>());

    EXPECT_EQ(result, tbb_result);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
