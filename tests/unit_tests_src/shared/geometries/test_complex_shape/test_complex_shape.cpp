#include <gtest/gtest.h>
#include "geometric_shape.h"
#include "complex_shape.h"

using namespace SPH;

TEST(test_ComplexShape, test_findNormalDirection)
{
	Vec3d halfsize_inner(1.0, 0.5, 0.25);
	Vec3d halfsize_outer(1.1, 0.6, 0.35);
	SimTK::Transform transfrom(Vec3d(1.0, 0.5, 0.25));

	ComplexShape body_shape;
	body_shape.add<GeometricShapeBrick>(halfsize_outer, transfrom);
	body_shape.substract<GeometricShapeBrick>(halfsize_inner, transfrom);
	Vec3d point(1.0, 0.5, -0.095);
	Vec3d normal = body_shape.findNormalDirection(point);

	EXPECT_EQ(Vec3d(0.0, 0.0, -1.0), normal);
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
