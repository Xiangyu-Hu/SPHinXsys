#include "particle_functors.hpp"
#include "update_cell_linked_list_sycl.hpp"

#include <gtest/gtest.h>
using namespace SPH;

StdVec<UnsignedInt> cell_size_list_tbb{0, 2, 3, 5, 0, 1, 3, 2, 5, 0, 21};
//StdVec<UnsignedInt> cell_size_list_sycl{0, 2, 3, 5, 0, 1, 3, 2, 5, 0, 21};

StdVec<UnsignedInt> result{0, 2, 5, 10, 10, 11, 14, 16, 21, 21, 21};

TEST(exclusive_scan, test_tbb_sycl)
{
    UnsignedInt *list_tbb = cell_size_list_tbb.data();
    //    UnsignedInt *list_sycl = cell_size_list_sycl.data();

    exclusive_scan(list_tbb, list_tbb + 11, list_tbb, ReduceSum<UnsignedInt>());

    EXPECT_EQ(cell_size_list_tbb, result);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
