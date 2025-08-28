#include "scalar_functions.h"
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

TEST(test_scalar_functions, test_WENO)
{
    Real v1 = 1.0;
    Real v2 = 2.0;
    Real v3 = 3.0;
    Real v4 = 4.0;

    Real Left_state = getLeftStateInWeno(v1, v2, v3, v4);
    Real Right_state = getRightStateInWeno(v4, v3, v2, v1);
    EXPECT_EQ(2.5, Left_state);
    EXPECT_EQ(2.5, Right_state);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
