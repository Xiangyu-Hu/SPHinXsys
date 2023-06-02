#include "anisotropic_kernel.hpp"
#include "base_kernel_includes_nonisotropic.h"
#include "kernel_wenland_c2_anisotropic.h"
#include "sphinxsys.h"
#include <gtest/gtest.h>

using namespace SPH;

Real getLinearProfile(const Vecd &input)
{
    return input.dot(Vecd::Ones());
}

Mat2d A{
    {0.0, 0.0},   // First row
    {0.0, 100.0}, // Second row
};

Real getQuadraticProfile(const Vecd &input)
{
    return (A * input * input.transpose()).trace();
}

TEST(test_anisotropic_kernel, test_Laplacian)
{
    Real PH = 1.0;                            // domain size in y direction
    Real ratio = 1.0;                         // dp_x / dp_y
    Real PL = ratio * PH;                     // domain size in x direction
    int y_num = 50;                           // particle number in y direction
    Real resolution_y = PH / Real(y_num);     // resolution in y direction
    Real resolution_x = ratio * resolution_y; // resolution in x direction
    int x_num = PL / resolution_x;            // Particle number in x direction , the same as particle number in y direction

    Vecd scaling_vector(1.0, 1.0 / ratio); // in x and y directions
    AnisotropicKernel<Anisotropic::KernelWendlandC2>
        wendland(1.15 * resolution_x, scaling_vector, Vecd(0.0, 0.0)); // no rotation introduced

    Mat2d transform_tensor = wendland.getCoordinateTransformationTensorG(scaling_vector, Vecd(0.0, 0.0)); // tensor

    Mat2d tensor_D = transform_tensor * transform_tensor.transpose();

    std::cout << transform_tensor << std::endl;

    Real predicted_kernel_integral = 0.0;
    Vecd predicted_gradient = Vecd(0.0, 0.0);
    Real predicted_laplacian = 0.0;

    Vecd pos_i = Vecd(resolution_x * 5.0, resolution_y * 5.0); // Particle i location
    Real V = resolution_y * resolution_x;                      // Particle volume
    for (int i = 0; i < (x_num + 1); i++)
    {
        for (int j = 0; j < (y_num + 1); j++)
        {
            Vecd pos_j(i * resolution_x, j * resolution_y);
            Vecd displacement = pos_i - pos_j;
            Real distance_ = displacement.norm();

            Real linear_profile = getLinearProfile(pos_j);
            Real linear_profile_pos_i = getLinearProfile(pos_i);

            Real parabolic_profile = getQuadraticProfile(pos_j);
            Real parabolic_profile_pos_i = getQuadraticProfile(pos_i);

            // if within cutoff radius
            if (wendland.checkIfWithinCutOffRadius(displacement))
            {
                predicted_kernel_integral += wendland.W(distance_, displacement) * V;
                Vecd eij_dwij_V = wendland.e(distance_, displacement) * wendland.dW(distance_, displacement) * V;
                predicted_gradient -= (linear_profile_pos_i - linear_profile) * eij_dwij_V;

                Vecd isotropic_displacement = transform_tensor * displacement;
                Vecd isotropic_eij = isotropic_displacement / (isotropic_displacement.norm() + TinyReal);

                Matd eij_tensor = isotropic_eij * isotropic_eij.transpose();
                Real weight_ = (tensor_D * ((Real(Dimensions) + 2.0) * eij_tensor - Mat2d::Identity())).trace();

                // tensor.det has already been added in factor_dw_2d in  dW(distance_, displacement), so there is no tensor.det
                predicted_laplacian += (parabolic_profile_pos_i - parabolic_profile) / (isotropic_displacement.norm() + TinyReal) *
                                       weight_ * V * wendland.dW(distance_, displacement);
            }
        }
    }

    EXPECT_EQ(1.0, predicted_kernel_integral);
    EXPECT_EQ(1.0, predicted_gradient[0]);
    EXPECT_EQ(1.0, predicted_gradient[1]);
    EXPECT_EQ(2.0 * A.trace(), predicted_laplacian);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
