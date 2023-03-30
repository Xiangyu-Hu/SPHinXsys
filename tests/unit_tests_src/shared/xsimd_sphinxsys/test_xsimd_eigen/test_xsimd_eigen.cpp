#include "xsimd_eigen.h"

#include <gtest/gtest.h>

using namespace SPH;

size_t vec_size = 10000000;

Real generateRandom()
{
	return (((double)rand() / (RAND_MAX)) - 0.5) * 2.0;
}

TEST(test_XsimdScalar, test_BasicOperations)
{
	StdLargeVec<Real> a, b;
	StdLargeVec<size_t> indexes;
	a.resize(vec_size);
	b.resize(vec_size);
	Real c(generateRandom());
	indexes.resize(vec_size);
	for (size_t i = 0; i < vec_size; ++i)
	{
		a[i] = generateRandom();
		b[i] = generateRandom();
		indexes[i] = i;
	}

	TickCount start_RealX = TickCount::now();

	RealX x_sum(0.0);
	RealX x_c(c);
	size_t floored_vec_size = vec_size - vec_size % XsimdSize;
	for (size_t i = 0; i < floored_vec_size; i += XsimdSize)
	{
		x_sum += (loadRealX(&a[i]) + gatherRealX<XsimdSize>(b, &indexes[i]) - x_c) / 2.0;
	}

	Real sum = reduceRealX(x_sum);
	for (size_t i = floored_vec_size; i < vec_size; ++i)
	{
		sum += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	TimeInterval internal_RealX = TickCount::now() - start_RealX;
	std::cout << "Total wall time for RealX computation: " << internal_RealX.seconds() << " seconds." << std::endl;

	TickCount start_Real = TickCount::now();

	Real sum_ref(0);
	for (size_t i = 0; i < vec_size; ++i)
	{
		sum_ref += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	TimeInterval internal_Real = TickCount::now() - start_Real;
	std::cout << "Total wall time for Real computation: " << internal_Real.seconds() << " seconds." << std::endl;

	EXPECT_FLOAT_EQ(sum_ref, sum);
}


TEST(test_XsimdVec2d, test_Vec2dOperations)
{
	StdLargeVec<Vec2d> a, b;
	StdLargeVec<size_t> indexes;
	a.resize(vec_size);
	b.resize(vec_size);
	Vec2d c(generateRandom(), generateRandom());
	indexes.resize(vec_size);
	for (size_t i = 0; i < vec_size; ++i)
	{
		a[i] = Vec2d(generateRandom(), generateRandom());
		b[i] = Vec2d(generateRandom(), generateRandom());
		indexes[i] = i;
	}

	Vec2dX x_sum = Vec2dX::Zero();
	size_t floored_vec_size = vec_size - vec_size % XsimdSize;
	for (size_t i = 0; i < floored_vec_size; i += XsimdSize)
	{
		x_sum += (loadVecdX<XsimdSize>(&a[i]) + gatherVecdX<XsimdSize>(b, &indexes[i]) - assignVecdX(c)) / 2.0;
	}

	Vec2d sum = reduceVecdX(x_sum);
	for (size_t i = floored_vec_size; i < vec_size; ++i)
	{
		sum += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	Vec2d sum_ref = Vec2d::Zero();
	for (size_t i = 0; i < vec_size; ++i)
	{
		sum_ref += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	EXPECT_FLOAT_EQ(sum_ref.norm(), sum.norm());
}

TEST(test_XsimdVec3d, test_Vec3dOperations)
{
	StdLargeVec<Vec3d> a, b;
	StdLargeVec<size_t> indexes;
	a.resize(vec_size);
	b.resize(vec_size);
	Vec3d c(generateRandom(), generateRandom(), generateRandom());
	indexes.resize(vec_size);
	for (size_t i = 0; i < vec_size; ++i)
	{
		a[i] = Vec3d(generateRandom(), generateRandom(), generateRandom());
		b[i] = Vec3d(generateRandom(), generateRandom(), generateRandom());
		indexes[i] = i;
	}

	TickCount start_Vec3dX = TickCount::now();

	Vec3dX x_sum = Vec3dX::Zero();
	size_t floored_vec_size = vec_size - vec_size % XsimdSize;
	for (size_t i = 0; i < floored_vec_size; i += XsimdSize)
	{
		x_sum += (loadVecdX<XsimdSize>(&a[i]) + gatherVecdX<XsimdSize>(b, &indexes[i]) - assignVecdX(c)) / 2.0;
	}

	Vec3d sum = reduceVecdX(x_sum);
	for (size_t i = floored_vec_size; i < vec_size; ++i)
	{
		sum += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	TimeInterval internal_Vec3dX = TickCount::now() - start_Vec3dX;
	std::cout << "Total wall time for Vec3dX computation: " << internal_Vec3dX.seconds() << " seconds." << std::endl;


	TickCount start_Vec3d = TickCount::now();

	Vec3d sum_ref = Vec3d::Zero();
	for (size_t i = 0; i < vec_size; ++i)
	{
		sum_ref += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	TimeInterval internal_Vec3d = TickCount::now() - start_Vec3d;
	std::cout << "Total wall time for Vec3d computation: " << internal_Vec3d.seconds() << " seconds." << std::endl;

	EXPECT_FLOAT_EQ(sum_ref.norm(), sum.norm());
}

TEST(test_XsimdMat2d, test_Mat2dOperations)
{
	StdLargeVec<Mat2d> a, b;
	StdLargeVec<size_t> indexes;
	a.resize(vec_size);
	b.resize(vec_size);
	Mat2d c{{generateRandom(), generateRandom()}, {generateRandom(), generateRandom()}};
	indexes.resize(vec_size);
	for (size_t i = 0; i < vec_size; ++i)
	{
		a[i] = Mat2d{{generateRandom(), generateRandom()}, {generateRandom(), generateRandom()}};
		b[i] = Mat2d{{generateRandom(), generateRandom()}, {generateRandom(), generateRandom()}};
		indexes[i] = i;
	}

	Mat2dX x_sum = Mat2dX::Zero();
	size_t floored_vec_size = vec_size - vec_size % XsimdSize;
	for (size_t i = 0; i < floored_vec_size; i += XsimdSize)
	{
		x_sum += (loadMatdX<XsimdSize>(&a[i]) + gatherMatdX<XsimdSize>(b, &indexes[i]) - assignMatdX(c)) / 2.0;
	}

	Mat2d sum = reduceMatdX(x_sum);
	for (size_t i = floored_vec_size; i < vec_size; ++i)
	{
		sum += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	Mat2d sum_ref = Mat2d::Zero();
	for (size_t i = 0; i < vec_size; ++i)
	{
		sum_ref += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	EXPECT_FLOAT_EQ(sum_ref.norm(), sum.norm());
}

TEST(test_XsimdMat3d, test_Mat3dOperations)
{
	StdLargeVec<Mat3d> a, b;
	StdLargeVec<size_t> indexes;
	a.resize(vec_size);
	b.resize(vec_size);
	Mat3d c{{generateRandom(), generateRandom(), generateRandom()},
			{generateRandom(), generateRandom(), generateRandom()},
			{generateRandom(), generateRandom(), generateRandom()}};
	indexes.resize(vec_size);
	for (size_t i = 0; i < vec_size; ++i)
	{
		a[i] = Mat3d{{generateRandom(), generateRandom(), generateRandom()},
					 {generateRandom(), generateRandom(), generateRandom()},
					 {generateRandom(), generateRandom(), generateRandom()}};
		b[i] = Mat3d{{generateRandom(), generateRandom(), generateRandom()},
					 {generateRandom(), generateRandom(), generateRandom()},
					 {generateRandom(), generateRandom(), generateRandom()}};
		indexes[i] = i;
	}

	TickCount start_Mat3dX = TickCount::now();

	Mat3dX x_sum = Mat3dX::Zero();
	size_t floored_vec_size = vec_size - vec_size % XsimdSize;
	for (size_t i = 0; i < floored_vec_size; i += XsimdSize)
	{
		x_sum += (loadMatdX<XsimdSize>(&a[i]) + gatherMatdX<XsimdSize>(b, &indexes[i]) - assignMatdX(c)) / 2.0;
	}

	Mat3d sum = reduceMatdX(x_sum);
	for (size_t i = floored_vec_size; i < vec_size; ++i)
	{
		sum += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	TimeInterval internal_Mat3dX = TickCount::now() - start_Mat3dX;
	std::cout << "Total wall time for Mat3dX computation: " << internal_Mat3dX.seconds() << " seconds." << std::endl;


	TickCount start_Mat3d = TickCount::now();

	Mat3d sum_ref = Mat3d::Zero();
	for (size_t i = 0; i < vec_size; ++i)
	{
		sum_ref += (a[i] + b[indexes[i]] - c) / 2.0;
	}

	TimeInterval internal_Mat3d = TickCount::now() - start_Mat3d;
	std::cout << "Total wall time for Mat3d computation: " << internal_Mat3d.seconds() << " seconds." << std::endl;

	EXPECT_FLOAT_EQ(sum_ref.norm(), sum.norm());
}

int main(int argc, char *argv[])
{
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
