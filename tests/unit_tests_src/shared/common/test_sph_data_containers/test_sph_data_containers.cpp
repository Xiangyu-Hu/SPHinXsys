#include "sphinxsys_containers.h"
#include <gtest/gtest.h>

using namespace SPH;

TEST(sph_data_containers, getIntersect)
{
    BoundingBoxd bb_1(Vec3d(-10, -20, -30), Vec3d(10, 20, 30));
    BoundingBoxd bb_2(Vec3d(-5, -10, -10), Vec3d(20, 10, 40));
    BoundingBoxd bb_ref(Vec3d(-5, -10, -10), Vec3d(10, 10, 30));
    EXPECT_EQ(bb_ref, bb_1.getIntersect(bb_2));
}

//=================================================================================================//
//=================================================================================================//
int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}