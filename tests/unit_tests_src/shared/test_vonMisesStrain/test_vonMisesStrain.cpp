#include "solid_particles.h"
#include <gtest/gtest.h>

using namespace SPH;

TEST(VonMisesStrainCalc, vonMisesStrain)
{
    Real tolerance = 1e-6;

    Matd strain_tensor_1{
        {1, 2, 3},
        {2, 5, 6},
        {3, 6, 9},
    };

    Real epsilonxx_1 = strain_tensor_1(0, 0);
    Real epsilonyy_1 = strain_tensor_1(1, 1);
    Real epsilonzz_1 = strain_tensor_1(2, 2);
    Real epsilonxy_1 = strain_tensor_1(0, 1);
    Real epsilonxz_1 = strain_tensor_1(0, 2);
    Real epsilonyz_1 = strain_tensor_1(1, 2);

    Real von_Mises_strain_1 = sqrt((1.0 / 3.0) * (pow(epsilonxx_1 - epsilonyy_1, 2) + pow(epsilonyy_1 - epsilonzz_1, 2) + pow(epsilonzz_1 - epsilonxx_1, 2)) + 2.0 * (pow(epsilonxy_1, 2) + pow(epsilonyz_1, 2) + pow(epsilonxz_1, 2)));
    Real von_Mises_strain_ref_1 = 11.4017543;

    Matd strain_tensor_2{
        {-5, 9, 1},
        {9, 4, -7},
        {1, -7, -12},
    };

    Real epsilonxx_2 = strain_tensor_2(0, 0);
    Real epsilonyy_2 = strain_tensor_2(1, 1);
    Real epsilonzz_2 = strain_tensor_2(2, 2);
    Real epsilonxy_2 = strain_tensor_2(0, 1);
    Real epsilonxz_2 = strain_tensor_2(0, 2);
    Real epsilonyz_2 = strain_tensor_2(1, 2);

    Real von_Mises_strain_2 = sqrt((1.0 / 3.0) * (pow(epsilonxx_2 - epsilonyy_2, 2) + pow(epsilonyy_2 - epsilonzz_2, 2) + pow(epsilonzz_2 - epsilonxx_2, 2)) + 2.0 * (pow(epsilonxy_2, 2) + pow(epsilonyz_2, 2) + pow(epsilonxz_2, 2)));
    Real von_Mises_strain_ref_2 = 19.7652894;

    EXPECT_NEAR(von_Mises_strain_1, von_Mises_strain_ref_1, tolerance);
    EXPECT_NEAR(von_Mises_strain_2, von_Mises_strain_ref_2, tolerance);
}
//=================================================================================================//
int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}