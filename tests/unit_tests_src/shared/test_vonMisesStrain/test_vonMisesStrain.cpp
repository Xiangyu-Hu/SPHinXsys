#include <gtest/gtest.h>
#include "solid_particles.h"

using namespace SPH;

TEST(VonMisesStrainCalc, vonMisesStrain)
{
    Real tolerance = 1e-6;
   
    Matd strain_tensor_1 = Matd( Vecd(1, 2, 3), Vecd(2, 5, 6), Vecd(3, 6, 9));
    
    Real epsilonxx_1 = strain_tensor_1(0, 0);
	Real epsilonyy_1 = strain_tensor_1(1, 1);
	Real epsilonzz_1 = strain_tensor_1(2, 2);
	Real epsilonxy_1 = strain_tensor_1(0, 1);
	Real epsilonxz_1 = strain_tensor_1(0, 2);
	Real epsilonyz_1 = strain_tensor_1(1, 2);
    
    Real von_Mises_strain_1 = sqrt( (1.0 / 3.0) * (std::pow(epsilonxx_1 - epsilonyy_1, 2.0) + std::pow(epsilonyy_1 - epsilonzz_1, 2.0) + std::pow(epsilonzz_1 - epsilonxx_1, 2.0))
		 + 2.0 * (std::pow(epsilonxy_1, 2.0) + std::pow(epsilonyz_1, 2.0) + std::pow(epsilonxz_1, 2.0)));
    Real von_Mises_strain_ref_1 = 11.4017543;
   
    Matd strain_tensor_2 = Matd( Vecd(-5, 9, 1), Vecd(9, 4, -7), Vecd(1, -7, -12));
    
    Real epsilonxx_2 = strain_tensor_2(0, 0);
	Real epsilonyy_2 = strain_tensor_2(1, 1);
	Real epsilonzz_2 = strain_tensor_2(2, 2);
	Real epsilonxy_2 = strain_tensor_2(0, 1);
	Real epsilonxz_2 = strain_tensor_2(0, 2);
	Real epsilonyz_2 = strain_tensor_2(1, 2);

    Real von_Mises_strain_2 = sqrt( (1.0 / 3.0) * (std::pow(epsilonxx_2 - epsilonyy_2, 2.0) + std::pow(epsilonyy_2 - epsilonzz_2, 2.0) + std::pow(epsilonzz_2 - epsilonxx_2, 2.0))
		 + 2.0 * (std::pow(epsilonxy_2, 2.0) + std::pow(epsilonyz_2, 2.0) + std::pow(epsilonxz_2, 2.0)));
    Real von_Mises_strain_ref_2 = 19.7652894;

	EXPECT_NEAR(von_Mises_strain_1, von_Mises_strain_ref_1, tolerance);
    EXPECT_NEAR(von_Mises_strain_2, von_Mises_strain_ref_2, tolerance);

}
//=================================================================================================//
int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}