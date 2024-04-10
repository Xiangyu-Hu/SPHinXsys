/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	vector_functions.h
 * @brief 	Basic functions for vector type data.
 * @author	Chi Zhang and Xiangyu Hu
 */
#ifndef VECTOR_FUNCTIONS_H
#define VECTOR_FUNCTIONS_H

#include "data_type.h"
#include "execution_event.h"
#include "execution_queue.hpp"

namespace SPH
{
Vec2d FirstAxisVector(const Vec2d &zero_vector);
Vec3d FirstAxisVector(const Vec3d &zero_vector);

Vec3d upgradeToVec3d(const Real &input);
Vec3d upgradeToVec3d(const Vec2d &input);
Vec3d upgradeToVec3d(const Vec3d &input);
Mat3d upgradeToMat3d(const Mat2d &input);
Mat3d upgradeToMat3d(const Mat3d &input);

Vecd degradeToVecd(const Vec3d &input);
Matd degradeToMatd(const Mat3d &input);

Mat2d getInverse(const Mat2d &A);
Mat3d getInverse(const Mat3d &A);
Mat2d getAverageValue(const Mat2d &A, const Mat2d &B);
Mat3d getAverageValue(const Mat3d &A, const Mat3d &B);
Mat2d inverseCholeskyDecomposition(const Mat2d &A);
Mat3d inverseCholeskyDecomposition(const Mat3d &A);
Mat2d getDiagonal(const Mat2d &A);
Mat3d getDiagonal(const Mat3d &A);

/** Real dot product between two matrices, resulting in a scalar value (sum of products of element-wise) */
Real CalculateBiDotProduct(Mat2d Matrix1, Mat2d Matrix2); // calculate Real dot
Real CalculateBiDotProduct(Mat3d Matrix1, Mat3d Matrix2); // calculate Real dot

/** get transformation matrix. */
Mat2d getTransformationMatrix(const Vec2d &direction_of_y);
Mat3d getTransformationMatrix(const Vec3d &direction_of_z);
Mat3d getTransformationMatrix(const Vec3d &direction_of_z, const Vec3d &direction_of_y);
/** get angle between two vectors. */
Real getCosineOfAngleBetweenTwoVectors(const Vec2d &vector_1, const Vec2d &vector_2);
Real getCosineOfAngleBetweenTwoVectors(const Vec3d &vector_1, const Vec3d &vector_2);

/** get orthogonal projection of a vector. */
Vec2d getVectorProjectionOfVector(const Vec2d &vector_1, const Vec2d &vector_2);
Vec3d getVectorProjectionOfVector(const Vec3d &vector_1, const Vec3d &vector_2);

/** von Mises stress from stress matrix */
Real getVonMisesStressFromMatrix(const Mat2d &sigma);
Real getVonMisesStressFromMatrix(const Mat3d &sigma);

/** principal strain or stress from strain or stress matrix */
Vec2d getPrincipalValuesFromMatrix(const Mat2d &A);
Vec3d getPrincipalValuesFromMatrix(const Mat3d &A);

/** get transformation matrix. */
Real getCrossProduct(const Vec2d &vector_1, const Vec2d &vector_2);
Vec3d getCrossProduct(const Vec3d &vector_1, const Vec3d &vector_2);


/** convert host Vecd to device Vecd */
inline DeviceVec2d hostToDeviceVecd(const Vec2d& host) { return {host[0], host[1]}; }
inline DeviceVec3d hostToDeviceVecd(const Vec3d& host) { return {host[0], host[1], host[2]}; }
inline DeviceArray2i hostToDeviceArrayi(const Array2i& host) { return {host[0], host[1]}; }
inline DeviceArray3i hostToDeviceArrayi(const Array3i& host) { return {host[0], host[1], host[2]}; }

/** convert device Vecd to host Vecd */
inline Vec2d deviceToHostVecd(const DeviceVec2d& device) { return {device[0], device[1]}; }
inline Vec3d deviceToHostVecd(const DeviceVec3d& device) { return {device[0], device[1], device[2]}; }

/** Initialize Vecd of zeros for host or device */
template<class V, class Enable = std::true_type> inline V VecdZero();
template<> inline DeviceVec2d VecdZero<DeviceVec2d, is_device_type_different_from_host<DeviceVec2d>>()
{
    return DeviceVec2d{0};
}
template<> inline DeviceVec3d VecdZero<DeviceVec3d, is_device_type_different_from_host<DeviceVec3d>>()
{
    return DeviceVec3d{0};
}
template<> inline Vec2d VecdZero() { return Vec2d::Zero(); }
template<> inline Vec3d VecdZero() { return Vec3d::Zero(); }

/* specialization for specific vector operations */
template<class RealType, int Dimension>
inline RealType VecdDot(const sycl::vec<RealType,Dimension>& v1, const sycl::vec<RealType,Dimension>& v2) {
    return sycl::dot(v1, v2);
}
template<class RealType, int Dimension>
inline RealType VecdDot(const Eigen::Matrix<RealType,Dimension,1>& v1, const Eigen::Matrix<RealType,Dimension,1>& v2) {
    return v1.dot(v2);
}
template<class RealType, int Dimension>
inline RealType VecdNorm(const sycl::vec<RealType,Dimension>& vec) {
    return sycl::length(vec);
}
template<class RealType, int Dimension>
inline RealType VecdNorm(const Eigen::Matrix<RealType,Dimension,1>& vec) {
    return vec.norm();
}
template<class RealType, int Dimension>
inline RealType VecdSquareNorm(const sycl::vec<RealType,Dimension>& vec) {
    return sycl::dot(vec, vec);
}
template<class RealType, int Dimension>
inline RealType VecdSquareNorm(const Eigen::Matrix<RealType,Dimension,1>& vec) {
    return vec.squaredNorm();
}
template<class RealType, int Dimension>
inline sycl::vec<RealType,Dimension> VecdMax(const sycl::vec<RealType,Dimension>& v1, const sycl::vec<RealType,Dimension>& v2) {
    return sycl::max(v1, v2);
}
template<class Type, int Dimension>
inline Eigen::Array<Type,Dimension,1> VecdMax(const Eigen::Array<Type,Dimension,1>& v1, const Eigen::Array<Type,Dimension,1>& v2) {
    return v1.max(v2);
}
template<class RealType, int Dimension>
inline sycl::vec<RealType,Dimension> VecdMin(const sycl::vec<RealType,Dimension>& v1, const sycl::vec<RealType,Dimension>& v2) {
    return sycl::min(v1, v2);
}
template<class Type, int Dimension>
inline Eigen::Array<Type,Dimension,1> VecdMin(const Eigen::Array<Type,Dimension,1>& v1, const Eigen::Array<Type,Dimension,1>& v2) {
    return v1.min(v2);
}
template<class VecType, std::size_t ...Index>
inline auto VecdFoldingProd_impl(const VecType& vec, std::index_sequence<Index...>) {
    return ( vec[Index] * ... );
}
template<class Type, int Dimension>
inline Type VecdFoldProduct(const sycl::vec<Type,Dimension>& vec) {
    return VecdFoldingProd_impl(vec, std::make_index_sequence<Dimension>());
}
template<class Type, int Dimension>
inline Type VecdFoldProduct(const Eigen::Array<Type,Dimension,1>& vec) {
    return VecdFoldingProd_impl(vec, std::make_index_sequence<Dimension>());
}

/* SYCL memory transfer utilities */
template<class T>
inline T* allocateDeviceData(std::size_t size) {
    return sycl::malloc_device<T>(size, execution::executionQueue.getQueue());
}

template<class T>
inline T* allocateSharedData(std::size_t size) {
    return sycl::malloc_shared<T>(size, execution::executionQueue.getQueue());
}

template<class T>
inline void freeDeviceData(T* device_mem) {
    sycl::free(device_mem, execution::executionQueue.getQueue());
}

template<class T>
inline execution::ExecutionEvent copyDataToDevice(const T* host, T* device, std::size_t size) {
    return execution::executionQueue.getQueue().memcpy(device, host, size*sizeof(T));
}

template<class HostType, class DeviceType, class TransformFunc>
execution::ExecutionEvent transformAndCopyDataToDevice(const HostType* host, DeviceType* device, std::size_t size, TransformFunc&& transformation) {
    auto* transformed_host_copy = new DeviceType[size];
    for(size_t i = 0; i < size; ++i)
        transformed_host_copy[i] = transformation(host[i]);
    return std::move(copyDataToDevice(transformed_host_copy, device, size).then([=](){ delete[] transformed_host_copy; }));
}

inline execution::ExecutionEvent copyDataToDevice(const Vec2d* host, DeviceVec2d* device, std::size_t size) {
    return transformAndCopyDataToDevice(host, device, size, [](auto vec) { return hostToDeviceVecd(vec); });
}

inline execution::ExecutionEvent copyDataToDevice(const Vec3d* host, DeviceVec3d* device, std::size_t size) {
    return transformAndCopyDataToDevice(host, device, size, [](auto vec) { return hostToDeviceVecd(vec); });
}

inline execution::ExecutionEvent copyDataToDevice(const Real* host, DeviceReal* device, std::size_t size) {
    return transformAndCopyDataToDevice(host, device, size, [](Real val) { return static_cast<DeviceReal>(val); });
}

template<class T>
inline execution::ExecutionEvent copyDataToDevice(const T& value, T* device, std::size_t size) {
    return execution::executionQueue.getQueue().fill(device, value, size);
}

template<class T>
inline execution::ExecutionEvent copyDataFromDevice(T* host, const T* device, std::size_t size) {
    return execution::executionQueue.getQueue().memcpy(host, device, size*sizeof(T));
}

template<class HostType, class DeviceType, class TransformFunc>
execution::ExecutionEvent transformAndCopyDataFromDevice(HostType* host, const DeviceType* device, std::size_t size, TransformFunc&& transformation) {
    auto* device_copy = new DeviceType[size];
    return std::move(copyDataFromDevice(device_copy, device, size).then([=](){
                                                                            for(size_t i = 0; i < size; ++i)
                                                                                host[i] = transformation(device_copy[i]);
                                                                            delete[] device_copy;
                                                                        }));
}


inline execution::ExecutionEvent copyDataFromDevice(Vec2d* host, const DeviceVec2d* device, std::size_t size) {
    return transformAndCopyDataFromDevice(host, device, size, [](auto vec) { return deviceToHostVecd(vec); });
}

inline execution::ExecutionEvent copyDataFromDevice(Vec3d* host, const DeviceVec3d* device, std::size_t size) {
    return transformAndCopyDataFromDevice(host, device, size, [](auto vec) { return deviceToHostVecd(vec); });
}

inline execution::ExecutionEvent copyDataFromDevice(Real* host, const DeviceReal* device, std::size_t size) {
    return transformAndCopyDataFromDevice(host, device, size, [](DeviceReal val) { return static_cast<Real>(val); });
}
} // namespace SPH
#endif // VECTOR_FUNCTIONS_H
