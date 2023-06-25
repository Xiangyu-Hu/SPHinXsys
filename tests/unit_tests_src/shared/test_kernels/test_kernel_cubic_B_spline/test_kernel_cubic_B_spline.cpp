#include "kernel_cubic_B_spline.h"
#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

TEST(test_Kernel_Cubic_B_Spline, test_W)
{
    Real q = 0.5;
    KernelCubicBSpline B_spline(1.0);
    Real W_1D = B_spline.W_1D(q);

    EXPECT_EQ(23.0 / 32.0, W_1D);
}

TEST(test_Kernel_Cubic_B_Spline, test_dW)
{
    Real q = 0.5;
    KernelCubicBSpline B_spline(1.0);
    Real dW_1D = B_spline.dW_1D(q);

    EXPECT_EQ(-15.0 / 16.0, dW_1D);
}

TEST(test_Kernel_Cubic_B_Spline, test_d2W)
{
    Real q = 0.5;
    KernelCubicBSpline B_spline(1.0);
    Real d2W_1D = B_spline.d2W_1D(q);

    EXPECT_EQ(-3.0 / 4.0, d2W_1D);
}
int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
