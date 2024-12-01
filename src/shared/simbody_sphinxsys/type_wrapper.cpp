#include "type_wrapper.h"
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
SimTKVec2 EigenToSimTK(const Vec2d &eigen_vector)
{
    return SimTKVec2((double)eigen_vector[0], (double)eigen_vector[1]);
}
//=================================================================================================//
SimTKVec3 EigenToSimTK(const Vec3d &eigen_vector)
{
    return SimTKVec3((double)eigen_vector[0], (double)eigen_vector[1], (double)eigen_vector[2]);
}
//=================================================================================================//
Vec2d SimTKToEigen(const SimTKVec2 &simTK_vector)
{
    return Vec2d((Real)simTK_vector[0], (Real)simTK_vector[1]);
}
//=================================================================================================//
Vec3d SimTKToEigen(const SimTKVec3 &simTK_vector)
{
    return Vec3d((Real)simTK_vector[0], (Real)simTK_vector[1], (Real)simTK_vector[2]);
}
SimTKMat22 EigenToSimTK(const Mat2d &eigen_matrix)
{
    return SimTKMat22((double)eigen_matrix(0, 0), (double)eigen_matrix(0, 1),
                      (double)eigen_matrix(1, 0), (double)eigen_matrix(1, 1));
}
//=================================================================================================//
SimTKMat33 EigenToSimTK(const Mat3d &eigen_matrix)
{
    return SimTKMat33((double)eigen_matrix(0, 0), (double)eigen_matrix(0, 1), (double)eigen_matrix(0, 2),
                      (double)eigen_matrix(1, 0), (double)eigen_matrix(1, 1), (double)eigen_matrix(1, 2),
                      (double)eigen_matrix(2, 0), (double)eigen_matrix(2, 1), (double)eigen_matrix(2, 2));
}
Mat2d SimTKToEigen(const SimTKMat22 &simTK_matrix)
{
    return Mat2d{
        {(Real)simTK_matrix(0, 0), (Real)simTK_matrix(0, 1)},
        {(Real)simTK_matrix(1, 0), (Real)simTK_matrix(1, 1)}};
}
//=================================================================================================//
Mat3d SimTKToEigen(const SimTKMat33 &simTK_matrix)
{
    return Mat3d{
        {(Real)simTK_matrix(0, 0), (Real)simTK_matrix(0, 1), (Real)simTK_matrix(0, 2)},
        {(Real)simTK_matrix(1, 0), (Real)simTK_matrix(1, 1), (Real)simTK_matrix(1, 2)},
        {(Real)simTK_matrix(2, 0), (Real)simTK_matrix(2, 1), (Real)simTK_matrix(2, 2)}};
}
//=================================================================================================//
} // namespace SPH
