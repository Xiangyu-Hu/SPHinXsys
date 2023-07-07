#include "vector_functions.h"
#include <gtest/gtest.h>

using namespace SPH;

TEST(small_vectors, CalculateBiDotProduct)
{
    Mat3d matrix_1{
        {6.5, 7.8, 0.0},
        {1.0, 2.0, 3.0},
        {10.0, 10.0, 17.3},
    };
    Mat3d matrix_2{
        {4.6, 4.9, 9.0},
        {1.0, 5.0, 7.0},
        {0.0, 0.0, 0.2},
    };
    EXPECT_EQ(CalculateBiDotProduct(matrix_1, matrix_2), 103.58);
}

TEST(small_vectors, getCosineOfAngleBetweenTwoVectors)
{
    Vec3d vector_1 = Vec3d(2, -4, 0);
    Vec3d vector_2 = Vec3d(3, 2, 5);

    Vec3d vector_3 = Vec3d(1, 0, 0);
    Vec3d vector_4 = Vec3d(0, 0, 1);

    Vec3d vector_5 = Vec3d(1, 0, 0);
    Vec3d vector_6 = Vec3d(-1, 0, 0);

    Real cos_teta_1 = getCosineOfAngleBetweenTwoVectors(vector_1, vector_2);
    Real cos_teta_ref_1 = -2 / (sqrt(20) * sqrt(38));

    Real cos_teta_2 = getCosineOfAngleBetweenTwoVectors(vector_3, vector_4);
    Real cos_teta_ref_2 = 0;

    Real cos_teta_3 = getCosineOfAngleBetweenTwoVectors(vector_5, vector_6);
    Real cos_teta_ref_3 = -1;

    EXPECT_EQ(cos_teta_1, cos_teta_ref_1);
    EXPECT_EQ(cos_teta_2, cos_teta_ref_2);
    EXPECT_EQ(cos_teta_3, cos_teta_ref_3);
}
//=================================================================================================//
TEST(small_vectors, getVectorProjectionOfVector)
{
    Vec3d vector_1 = Vec3d(1, 2, 3);
    Vec3d vector_2 = Vec3d(2, -1, 4);

    // angle exactly 90°
    Vec3d vector_3 = Vec3d(0, 0, 0);
    Vec3d vector_4 = Vec3d(0, 1, 0);

    // vectors in opposite direction
    Vec3d vector_5 = Vec3d(1, 1, 0);
    Vec3d vector_6 = Vec3d(-1, -1, 0);

    Vec3d vector_7 = Vec3d(5.44009282066326e-15, -0.000058164135562623, 2.22856150706995e-07);
    Vec3d vector_8 = Vec3d(-9.6283166736093e-11, 0.999698258170137, -0.00394363294839056);

    Vec3d proj_vector_1 = getVectorProjectionOfVector(vector_1, vector_2);
    Vec3d proj_vector_1_ref = Vec3d(1.142857, -0.571428, 2.285714);

    Vec3d proj_vector_2 = getVectorProjectionOfVector(vector_3, vector_4);
    Vec3d proj_vector_2_ref = Vec3d(0, 0, 0);

    Vec3d proj_vector_3 = getVectorProjectionOfVector(vector_5, vector_6);
    Vec3d proj_vector_3_ref = Vec3d(1, 1, 0);

    // values from spring normal to surface
    Vec3d proj_vector_4 = getVectorProjectionOfVector(vector_7, vector_8);
    Vec3d proj_vector_4_ref = Vec3d(5.60191499112963e-15, -0.000058164109562364, 2.29447132681609e-07);

    for (size_t i = 0; i < 3; ++i)
    {
        EXPECT_NEAR(proj_vector_1[i], proj_vector_1_ref[i], 1.0e-6);
        EXPECT_NEAR(proj_vector_2[i], proj_vector_2_ref[i], 1.0e-6);
        EXPECT_NEAR(proj_vector_3[i], proj_vector_3_ref[i], 1.0e-6);
        EXPECT_NEAR(proj_vector_4[i], proj_vector_4_ref[i], 1.0e-8);
    }
}
//=================================================================================================//
TEST(TransformationMatrix, getTransformationMatrix)
{
    for (int i = 0; i < 10000; ++i)
    {
        // Generate 3D unit vectors with a spherically symmetric probability distribution
        // according to the link: https://mathworld.wolfram.com/SpherePointPicking.html
        Real u = (Real)rand() / (RAND_MAX);
        Real v = (Real)rand() / (RAND_MAX);
        Real theta = 2.0 * Pi * u;
        Real cos_phi = 2.0 * v - 1.0;
        Real sin_phi = sqrt(1.0 - cos_phi * cos_phi);
        Vec3d unit_vector = Vec3d(sin_phi * cos(theta), sin_phi * sin(theta), cos_phi);
        Mat3d transformation_matrix = getTransformationMatrix(unit_vector);
        Vec3d transformed_unit_vector = transformation_matrix.transpose() * (transformation_matrix * unit_vector);

        EXPECT_NEAR(unit_vector[0], transformed_unit_vector[0], 1.0e-14);
        EXPECT_NEAR(unit_vector[1], transformed_unit_vector[1], 1.0e-14);
        EXPECT_NEAR(unit_vector[2], transformed_unit_vector[2], 1.0e-14);
    }
}
//=================================================================================================//
TEST(small_vectors, getPrincipalValuesFromMatrix)
{
    Mat3d stress_tensor_1{
        {50, 30, 20},
        {30, -20, -10},
        {20, -10, 10},
    };
    Vec3d principal_stress_1 = getPrincipalValuesFromMatrix(stress_tensor_1);
    Vec3d principal_stress_1_ref = {65.527, 11.531, -37.058};
    for (size_t i = 0; i < 3; i++)
        EXPECT_NEAR(principal_stress_1[i], principal_stress_1_ref[i], 1e-3);

    Mat3d stress_tensor_2{
        {12.5, 0, 0},
        {0, -50.4, 0},
        {0, 0, 15.3},
    };
    Vec3d principal_stress_2 = getPrincipalValuesFromMatrix(stress_tensor_2);
    Vec3d principal_stress_2_ref = {15.3, 12.5, -50.4};
    for (size_t i = 0; i < 3; i++)
        EXPECT_EQ(principal_stress_2[i], principal_stress_2_ref[i]);
}

//=================================================================================================//
//=================================================================================================//
int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}