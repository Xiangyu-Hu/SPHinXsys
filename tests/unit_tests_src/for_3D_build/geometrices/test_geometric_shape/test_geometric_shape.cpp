#include <gtest/gtest.h>
#include "geometric_shape.h"

using namespace SPH;

TEST(test_GeometricShapeBrick, test_findBounds)
{
	Vec3d halfsize(1.0, 0.5, 0.25);
	SimTK::Transform transfrom(Vec3d(1.0, 0.5, 0.25));

	GeometricShapeBrick brick(halfsize, transfrom);
	BoundingBox box = brick.findBounds();

	EXPECT_EQ(BoundingBox(Vec3d(0.0, 0.0, 0.0), Vec3d(2.0, 1.0, 0.5)), box);
}

int main(int argc, char* argv[])
{	
	testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
