#include "geometric_shape.h"
#include "transform_shape.h"
#include <gtest/gtest.h>

using namespace SPH;

TEST(test_GeometricShapeBox, test_findBounds)
{
    Vec3d halfsize(1.0, 0.5, 0.25);
    Transform transform(Vec3d(1.0, 0.5, 0.25));

    TransformShape<GeometricShapeBox> brick(transform, halfsize);
    BoundingBox box = brick.getBounds();

    EXPECT_EQ(BoundingBox(Vec3d(0.0, 0.0, 0.0), Vec3d(2.0, 1.0, 0.5)), box);
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
