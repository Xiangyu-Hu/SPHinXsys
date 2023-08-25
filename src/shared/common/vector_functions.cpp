#include "vector_functions.h"
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
Vec2d FirstAxisVector(const Vec2d &zero_vector)
{
    return Vec2d(1.0, 0.0);
}
//=================================================================================================//
Vec3d FirstAxisVector(const Vec3d &zero_vector)
{
    return Vec3d(1.0, 0.0, 0.0);
};
//=================================================================================================//
Vec3d upgradeToVec3d(const Real &input)
{
    return Vec3d(input, 0.0, 0.0);
}
//=================================================================================================//
Vec3d upgradeToVec3d(const Vec2d &input)
{
    return Vec3d(input[0], input[1], 0.0);
}
//=================================================================================================//
Vec3d upgradeToVec3d(const Vec3d &input)
{
    return input;
}
//=================================================================================================//
Mat3d upgradeToMat3d(const Mat2d &input)
{
    Mat3d output = Mat3d::Zero();
    output.block<2, 2>(0, 0) = input;
    return output;
}
//=================================================================================================//
Mat3d upgradeToMat3d(const Mat3d &input)
{
    return input;
}
//=================================================================================================//
void degradeToVecd(const Vec3d &input, Vec2d &output)
{
    output[0] = input[0];
    output[1] = input[1];
}
//=================================================================================================//
void degradeToMatd(const Mat3d &input, Mat2d &output)
{
    for (int i = 0; i != 2; i++)
        for (int j = 0; j != 2; j++)
            output(i, j) = input(i, j);
}
//=================================================================================================//
Mat2d getInverse(const Mat2d &A)
{
    Mat2d minv = Mat2d::Zero();
    Real det = A(0, 0) * A(1, 1) - A(0, 1) * A(1, 0);
    Real invdet = 1.0 / det;
    minv(0, 0) = A(1, 1) * invdet;
    minv(0, 1) = -A(0, 1) * invdet;
    minv(1, 0) = -A(1, 0) * invdet;
    minv(1, 1) = A(0, 0) * invdet;
    return minv;
}
//=================================================================================================//
Mat3d getInverse(const Mat3d &A)
{
    Real det = A(0, 0) * (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) -
               A(0, 1) * (A(1, 0) * A(2, 2) - A(1, 2) * A(2, 0)) +
               A(0, 2) * (A(1, 0) * A(2, 1) - A(1, 1) * A(2, 0));

    Real invdet = 1 / det;
    Mat3d minv = Mat3d::Zero();
    minv(0, 0) = (A(1, 1) * A(2, 2) - A(2, 1) * A(1, 2)) * invdet;
    minv(0, 1) = (A(0, 2) * A(2, 1) - A(0, 1) * A(2, 2)) * invdet;
    minv(0, 2) = (A(0, 1) * A(1, 2) - A(0, 2) * A(1, 1)) * invdet;
    minv(1, 0) = (A(1, 2) * A(2, 0) - A(1, 0) * A(2, 2)) * invdet;
    minv(1, 1) = (A(0, 0) * A(2, 2) - A(0, 2) * A(2, 0)) * invdet;
    minv(1, 2) = (A(1, 0) * A(0, 2) - A(0, 0) * A(1, 2)) * invdet;
    minv(2, 0) = (A(1, 0) * A(2, 1) - A(2, 0) * A(1, 1)) * invdet;
    minv(2, 1) = (A(2, 0) * A(0, 1) - A(0, 0) * A(2, 1)) * invdet;
    minv(2, 2) = (A(0, 0) * A(1, 1) - A(1, 0) * A(0, 1)) * invdet;

    return minv;
}
//=================================================================================================//
Mat2d getAverageValue(const Mat2d &A, const Mat2d &B)
{
    Mat2d C = Mat2d::Identity();
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            C(i, j) = 2.0 * A(i, j) * B(i, j) / (A(i, j) + B(i, j) + TinyReal);
        }
    }
    return C;
}
//=================================================================================================//
Mat3d getAverageValue(const Mat3d &A, const Mat3d &B)
{
    Mat3d C = Mat3d::Identity();
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            C(i, j) = 2.0 * A(i, j) * B(i, j) / (A(i, j) + B(i, j) + TinyReal);
        }
    }
    return C;
}
//=================================================================================================//
Mat2d inverseCholeskyDecomposition(const Mat2d &A)
{
    Mat2d lower = Mat2d::Zero();
    int n = 2;
    /** Decomposing a matrix into Lower Triangular. */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < (i + 1); j++)
        {
            Real sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum += lower(i, k) * lower(j, k);
            }
            if (i == j)
            {
                lower(i, j) = sqrt(A(i, i) - sum);
            }
            else
            {
                lower(i, j) = (1.0 / lower(j, j) * (A(i, j) - sum));
            }
        }
    }

    return lower.inverse();
}
//=================================================================================================//
Mat3d inverseCholeskyDecomposition(const Mat3d &A)
{
    Mat3d lower = Mat3d::Zero();
    int n = 3;
    /** Decomposing a matrix into Lower Triangular. */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < (i + 1); j++)
        {
            Real sum = 0;
            for (int k = 0; k < j; k++)
            {
                sum += lower(i, k) * lower(j, k);
            }
            if (i == j)
            {
                lower(i, j) = sqrt(A(i, i) - sum);
            }
            else
            {
                lower(i, j) = (1.0 / lower(j, j) * (A(i, j) - sum));
            }
        }
    }

    return lower.inverse();
}
//=================================================================================================//
Mat2d getTransformationMatrix(const Vec2d &direction_of_y)
{
    Mat2d transformation_matrix = Mat2d::Zero();
    transformation_matrix(0, 0) = direction_of_y[1];
    transformation_matrix(0, 1) = -direction_of_y[0];
    transformation_matrix(1, 0) = direction_of_y[0];
    transformation_matrix(1, 1) = direction_of_y[1];

    return transformation_matrix;
}
//=================================================================================================//
Mat3d getTransformationMatrix(const Vec3d &direction_of_z)
{
    Mat3d transformation_matrix = Mat3d::Zero();
    Real temp = 1.0 + direction_of_z[2];
    Real fraction = temp / (temp * temp + Eps);
    transformation_matrix(0, 0) = direction_of_z[2] + pow(direction_of_z[1], 2) * fraction;
    transformation_matrix(0, 1) = -direction_of_z[0] * direction_of_z[1] * fraction;
    transformation_matrix(0, 2) = -direction_of_z[0];
    transformation_matrix(1, 0) = transformation_matrix(0, 1);
    transformation_matrix(1, 1) = direction_of_z[2] + pow(direction_of_z[0], 2) * fraction;
    transformation_matrix(1, 2) = -direction_of_z[1];
    transformation_matrix(2, 0) = direction_of_z[0];
    transformation_matrix(2, 1) = direction_of_z[1];
    transformation_matrix(2, 2) = direction_of_z[2];

    return transformation_matrix;
}
//=================================================================================================//
//=================================================================================================//
Mat3d getTransformationMatrix(const Vec3d &direction_of_z, const Vec3d &direction_of_y)
{
    Mat3d transformation_matrix = Mat3d::Zero();
    Vec3d direction_of_x = direction_of_y.cross(direction_of_z);
    transformation_matrix.row(0) = direction_of_x.transpose();
    transformation_matrix.row(1) = direction_of_y.transpose();
    transformation_matrix.row(2) = direction_of_z.transpose();
     
    return transformation_matrix;
}

//=================================================================================================//

Mat2d getDiagonal(const Mat2d &A)
{
    Mat2d diag = Mat2d::Identity();
    diag(0, 0) = A(0, 0);
    diag(1, 1) = A(1, 1);

    return diag;
}
Mat3d getDiagonal(const Mat3d &A)
{
    Mat3d diag = Mat3d::Identity();
    diag(0, 0) = A(0, 0);
    diag(1, 1) = A(1, 1);
    diag(2, 2) = A(2, 2);

    return diag;
}
//=================================================================================================//
Real CalculateBiDotProduct(Mat2d Matrix1, Mat2d Matrix2)
{
    Real product = 0;
    for (int i = 0; i < 2; i++)
    {
        for (int j = 0; j < 2; j++)
        {
            product += Matrix1(i, j) * Matrix2(i, j);
        }
    }
    return product;
}
Real CalculateBiDotProduct(Mat3d Matrix1, Mat3d Matrix2)
{
    Real product = 0;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            product += Matrix1(i, j) * Matrix2(i, j);
        }
    }
    return product;
}
//=================================================================================================//
Real getCosineOfAngleBetweenTwoVectors(const Vec2d &vector_1, const Vec2d &vector_2)
{
    // returns the cosine of the angle between two vectors
    Real dot_product_1 = 0.0;
    for (int i = 0; i < vector_1.size(); i++)
    {
        dot_product_1 += vector_1[i] * vector_2[i];
    }
    Real cos_theta = dot_product_1 / (vector_1.norm() * vector_2.norm());

    return cos_theta;
}
//=================================================================================================//
Real getCosineOfAngleBetweenTwoVectors(const Vec3d &vector_1, const Vec3d &vector_2)
{
    // returns the cosine of the angle between two vectors
    Real dot_product_1 = 0.0;
    for (int i = 0; i < vector_1.size(); i++)
    {
        dot_product_1 += vector_1[i] * vector_2[i];
    }
    Real cos_theta = dot_product_1 / (vector_1.norm() * vector_2.norm());

    return cos_theta;
}
//=================================================================================================//
Vec2d getVectorProjectionOfVector(const Vec2d &vector_1, const Vec2d &vector_2)
{
    // get the projection of the vector_1 on vector 2, which is parallel to the vector_2, meaning it is the vector_2 * scalar
    Real dot_product_1 = 0.0;
    Real dot_product_2 = vector_2.squaredNorm();
    for (int i = 0; i < vector_1.size(); i++)
    {
        dot_product_1 += vector_1[i] * vector_2[i];
    }
    // get scalar, which to multiply n_0 with
    Real lambda = dot_product_1 / dot_product_2;
    Vec2d proj_vector_1 = lambda * vector_2;

    return proj_vector_1;
}
//=================================================================================================//
Vec3d getVectorProjectionOfVector(const Vec3d &vector_1, const Vec3d &vector_2)
{
    // get the projection of the vector_1 on vector 2, which is parallel to the vector_2, meaning it is the vector_2 * scalar
    Real dot_product_1 = 0.0;
    Real dot_product_2 = vector_2.squaredNorm();
    for (int i = 0; i < vector_1.size(); i++)
    {
        dot_product_1 += vector_1[i] * vector_2[i];
    }
    // get scalar, which to multiply n_0 with
    Real lambda = dot_product_1 / dot_product_2;
    Vec3d proj_vector_1 = lambda * vector_2;

    return proj_vector_1;
}
//=================================================================================================//
Real getVonMisesStressFromMatrix(const Mat2d &sigma)
{
    Real sigmaxx = sigma(0, 0);
    Real sigmayy = sigma(1, 1);
    Real sigmaxy = sigma(0, 1);

    return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy - sigmaxx * sigmayy + 3.0 * sigmaxy * sigmaxy);
}
//=================================================================================================//
Real getVonMisesStressFromMatrix(const Mat3d &sigma)
{
    Real sigmaxx = sigma(0, 0);
    Real sigmayy = sigma(1, 1);
    Real sigmazz = sigma(2, 2);
    Real sigmaxy = sigma(0, 1);
    Real sigmaxz = sigma(0, 2);
    Real sigmayz = sigma(1, 2);

    return sqrt(sigmaxx * sigmaxx + sigmayy * sigmayy + sigmazz * sigmazz -
                sigmaxx * sigmayy - sigmaxx * sigmazz - sigmayy * sigmazz +
                3.0 * (sigmaxy * sigmaxy + sigmaxz * sigmaxz + sigmayz * sigmayz));
}
//=================================================================================================//
Vec2d getPrincipalValuesFromMatrix(const Mat2d &A)
{
    Eigen::EigenSolver<EigMat> ces(A, /* computeEigenvectors = */ false);
    auto eigen_values = ces.eigenvalues();

    std::vector<Real> sorted_values = {
        Real(eigen_values(0).real()),
        Real(eigen_values(1).real())};
    // first sort into ascending order, and then reverse them
    std::sort(sorted_values.begin(), sorted_values.end());
    std::reverse(sorted_values.begin(), sorted_values.end());

    return {sorted_values[0], sorted_values[1]};
}
//=================================================================================================//
Vec3d getPrincipalValuesFromMatrix(const Mat3d &A)
{
    Eigen::EigenSolver<EigMat> ces(A, /* computeEigenvectors = */ false);
    auto eigen_values = ces.eigenvalues();

    std::vector<Real> sorted_values = {
        Real(eigen_values(0).real()),
        Real(eigen_values(1).real()),
        Real(eigen_values(2).real())};
    // first sort into ascending order, and then reverse them
    std::sort(sorted_values.begin(), sorted_values.end());
    std::reverse(sorted_values.begin(), sorted_values.end());

    return {sorted_values[0], sorted_values[1], sorted_values[2]};
}
//=================================================================================================//
Real getCrossProduct(const Vec2d &vector_1, const Vec2d &vector_2)
{
    return vector_1[1] * vector_2[0] - vector_1[0] * vector_2[1];
}
//=================================================================================================//
Vec3d getCrossProduct(const Vec3d &vector_1, const Vec3d &vector_2)
{ // Eigen cross product for only have 3D vector
    return vector_1.cross(vector_2);
}
//=================================================================================================//
} // namespace SPH
