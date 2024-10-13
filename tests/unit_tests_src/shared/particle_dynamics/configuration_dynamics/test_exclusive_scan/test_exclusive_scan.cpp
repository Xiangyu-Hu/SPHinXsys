#include "particle_functors.h"
#include "base_configuration_dynamics.h"
#include "particle_iterators_ck.h"

#include <gtest/gtest.h>
using namespace SPH;

StdVec<UnsignedInt> cell_size_list{3, 2, 3, 5, 0, 1, 3, 2, 5, 1, 0};
StdVec<UnsignedInt> solution{0, 3, 5, 8, 13, 13, 14, 17, 19, 24, 25};

UnsignedInt list_size = cell_size_list.size();
StdVec<UnsignedInt> result(list_size, 0);
StdVec<UnsignedInt> tbb_result = result;

TEST(exclusive_scan, test_tbb)
{
    UnsignedInt *list_data = cell_size_list.data();

    UnsignedInt sum = exclusive_scan(SequencedPolicy{}, list_data, result.data(), list_size, PlusUnsignedInt<SequencedPolicy>::type());
    UnsignedInt sum_tbb = exclusive_scan(ParallelPolicy{}, list_data, tbb_result.data(), list_size, PlusUnsignedInt<ParallelPolicy>::type());

    EXPECT_EQ(solution, result);
    EXPECT_EQ(solution, tbb_result);
    EXPECT_EQ(sum, sum_tbb);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
