#include <gtest/gtest.h>
#include "base_kernel_includes_nonisotropic.h"
#include "anisotropic_kernel.h"
#include "anisotropic_kernel.hpp"
#include "kernel_wenland_c2_anisotropic.h"
#include "sphinxsys.h"

using namespace SPH;

TEST(test_anisotropic_kernel, test_Laplacian)
{
	int y_num = 10; // particle number in y direction
	Real PH = 1.0; // domain size PH , PL
	Real ratio_ = 4.0; //  dp_x /dp_y
	Real PL = ratio_* PH;

	Vec2d scaling_vector = Vec2d(1.0, 1.0 / ratio_);    // For tensor
	Real resolution_ref = PH / Real(y_num);             //resolution in y direction
	Real resolution_ref_large = ratio_ * resolution_ref;       //resolution in x direction
	Real V_ = resolution_ref * resolution_ref_large;            // Particle volume 
	Vec2d center = Vec2d(resolution_ref_large * 5.0, resolution_ref * 5.0);   // Particle i location 

	int x_num = PL / resolution_ref_large;     // Particle number in x direction , the same as particle number in y direction
 
 	AnisotropicKernel<Anisotropic::KernelWendlandC2>  
 	    wendland(1.15 * resolution_ref_large, scaling_vector,  Vec2d(0.0, 0.0));
 	 

	Mat2d transform_tensor_ = wendland.getCoordinateTransformationTensorG(scaling_vector,  Vec2d(0.0, 0.0));  // tensor

    Mat2d Tensor_D = transform_tensor_ * transform_tensor_.transpose();
	Real trace_D = Tensor_D.trace();

	std::cout<< transform_tensor_<< std::endl;


	Real sum = 0.0;
	Vec2d first_order_rate = Vec2d(0.0, 0.0);
	Real second_order_rate =  0.0;

	 for (int i = 0; i < (x_num + 1); i++)
	{
		for (int j = 0; j < (y_num + 1); j++)
		{
			Real x = i * resolution_ref_large;
			Real y = j * resolution_ref;
			Vec2d displacement =  center - Vec2d(x, y) ;
			Real distance_ = displacement.norm();

			Real  sarutration_first_order =  y + x ;
			Real  sarutration_first_order_center =  center[1] + center[0];

			Real  sarutration_x = x * x  ;
			Real  sarutration_center_x = center[0] *  center[0];
		 
		 // if withincutoffradius
			if (wendland.checkIfWithinCutOffRadius(displacement)) 
			{
				Vec2d  eij_dwij_V = wendland.e(distance_, displacement)* wendland.dW(distance_, displacement) * V_;
              	sum += wendland.W(distance_, displacement)* V_;
				first_order_rate -= (sarutration_first_order_center - sarutration_first_order)* eij_dwij_V;



                Vec2d   isotropic_displacement = transform_tensor_ * displacement;
                Vec2d   isotropic_eij = isotropic_displacement / (isotropic_displacement.norm()+ TinyReal);
				
				// eij cross product, it seems cross() function can only be used in three dimensions
                Mat2d  eij_tensor =   Mat2d ({{isotropic_eij[0] * isotropic_eij[0], isotropic_eij[0] * isotropic_eij[1]}, 
				 							{isotropic_eij[1] * isotropic_eij[0], isotropic_eij[1] * isotropic_eij[1]}});

	 
	 			Real weight_ =  (Tensor_D *  eij_tensor).trace()* (2.0 + 2.0) - trace_D ;
              
			    // tensor.det has already been added in factor_dw_2d in  dW(distance_, displacement), so there is no tensor.det
				second_order_rate +=  (sarutration_center_x  - sarutration_x) 
									 /(isotropic_displacement.norm()+ TinyReal)  * weight_
					 				* V_ * wendland.dW(distance_, displacement); 
			
			}		
		}

	}

	EXPECT_EQ(1.0, sum); 
	EXPECT_EQ(1.0, first_order_rate[0]);
	EXPECT_EQ(1.0, first_order_rate[1]);
	EXPECT_EQ(2.0, second_order_rate);
}
 
int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
