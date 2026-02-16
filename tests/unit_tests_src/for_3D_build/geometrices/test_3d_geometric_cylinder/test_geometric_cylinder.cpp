#include "geometric_shape.h"
#include "transform_geometry.h"
#include <gtest/gtest.h>

using namespace SPH;

auto tolerance = []()
{ return 100 * Eps; };

// Test GeometricShapeCylinder basic construction and bounds
TEST(test_GeometricShapeCylinder, test_findBounds)
{
    Vec3d center(0.0, 0.0, 0.0);
    Vec3d axis(0.0, 0.0, 1.0);  // z-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);
    BoundingBoxd box = cylinder.getBounds();

    // For a cylinder along z-axis: x,y should be [-radius, radius], z should be [-halflength, halflength]
    EXPECT_LE(std::abs(box.lower_[0] - (-radius)), tolerance());
    EXPECT_LE(std::abs(box.lower_[1] - (-radius)), tolerance());
    EXPECT_LE(std::abs(box.lower_[2] - (-halflength)), tolerance());
    EXPECT_LE(std::abs(box.upper_[0] - radius), tolerance());
    EXPECT_LE(std::abs(box.upper_[1] - radius), tolerance());
    EXPECT_LE(std::abs(box.upper_[2] - halflength), tolerance());
}

// Test containment - point inside cylinder
TEST(test_GeometricShapeCylinder, test_contain_inside)
{
    Vec3d center(0.0, 0.0, 0.0);
    Vec3d axis(0.0, 0.0, 1.0);  // z-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);

    // Point clearly inside
    Vec3d point_inside(0.5, 0.0, 1.0);
    EXPECT_TRUE(cylinder.checkContain(point_inside));

    // Point at center
    Vec3d point_center(0.0, 0.0, 0.0);
    EXPECT_TRUE(cylinder.checkContain(point_center));
}

// Test containment - point outside cylinder
TEST(test_GeometricShapeCylinder, test_contain_outside)
{
    Vec3d center(0.0, 0.0, 0.0);
    Vec3d axis(0.0, 0.0, 1.0);  // z-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);

    // Point outside radially
    Vec3d point_outside_radial(2.0, 0.0, 0.0);
    EXPECT_FALSE(cylinder.checkContain(point_outside_radial));

    // Point outside axially
    Vec3d point_outside_axial(0.0, 0.0, 3.0);
    EXPECT_FALSE(cylinder.checkContain(point_outside_axial));

    // Point outside both
    Vec3d point_outside_both(2.0, 0.0, 3.0);
    EXPECT_FALSE(cylinder.checkContain(point_outside_both));
}

// Test containment - point on boundary
TEST(test_GeometricShapeCylinder, test_contain_boundary)
{
    Vec3d center(0.0, 0.0, 0.0);
    Vec3d axis(0.0, 0.0, 1.0);  // z-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);

    // Point on cylindrical surface
    Vec3d point_on_surface(1.0, 0.0, 1.0);
    EXPECT_TRUE(cylinder.checkContain(point_on_surface));

    // Point on end cap
    Vec3d point_on_cap(0.5, 0.0, 2.0);
    EXPECT_TRUE(cylinder.checkContain(point_on_cap));
}

// Test closest point - for point inside cylinder
TEST(test_GeometricShapeCylinder, test_closest_point_inside)
{
    Vec3d center(0.0, 0.0, 0.0);
    Vec3d axis(0.0, 0.0, 1.0);  // z-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);

    // Point inside, closer to side
    Vec3d point_inside(0.9, 0.0, 0.0);
    Vec3d closest = cylinder.findClosestPoint(point_inside);
    // Should project to cylindrical surface
    EXPECT_LE(std::abs(closest[0] - 1.0), tolerance());
    EXPECT_LE(std::abs(closest[1] - 0.0), tolerance());
    EXPECT_LE(std::abs(closest[2] - 0.0), tolerance());
}

// Test closest point - for point outside radially
TEST(test_GeometricShapeCylinder, test_closest_point_outside_radial)
{
    Vec3d center(0.0, 0.0, 0.0);
    Vec3d axis(0.0, 0.0, 1.0);  // z-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);

    // Point outside radially
    Vec3d point_outside(2.0, 0.0, 1.0);
    Vec3d closest = cylinder.findClosestPoint(point_outside);
    // Should project to cylindrical surface at (1.0, 0.0, 1.0)
    EXPECT_LE(std::abs(closest[0] - 1.0), tolerance());
    EXPECT_LE(std::abs(closest[1] - 0.0), tolerance());
    EXPECT_LE(std::abs(closest[2] - 1.0), tolerance());
}

// Test closest point - for point outside axially
TEST(test_GeometricShapeCylinder, test_closest_point_outside_axial)
{
    Vec3d center(0.0, 0.0, 0.0);
    Vec3d axis(0.0, 0.0, 1.0);  // z-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);

    // Point outside axially
    Vec3d point_outside(0.5, 0.0, 3.0);
    Vec3d closest = cylinder.findClosestPoint(point_outside);
    // Should project to top cap at (0.5, 0.0, 2.0)
    EXPECT_LE(std::abs(closest[0] - 0.5), tolerance());
    EXPECT_LE(std::abs(closest[1] - 0.0), tolerance());
    EXPECT_LE(std::abs(closest[2] - 2.0), tolerance());
}

// Test with translated cylinder
TEST(test_GeometricShapeCylinder, test_translated_cylinder)
{
    Vec3d center(5.0, 3.0, 2.0);
    Vec3d axis(1.0, 0.0, 0.0);  // x-axis
    Real radius = 0.5;
    Real halflength = 1.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);

    // Point inside translated cylinder
    Vec3d point_inside(5.0, 3.0, 2.0);
    EXPECT_TRUE(cylinder.checkContain(point_inside));

    // Point outside
    Vec3d point_outside(5.0, 5.0, 2.0);
    EXPECT_FALSE(cylinder.checkContain(point_outside));
}

// Test with different axis orientations
TEST(test_GeometricShapeCylinder, test_different_axis)
{
    Vec3d center(0.0, 0.0, 0.0);
    Vec3d axis(1.0, 1.0, 0.0);  // diagonal in xy-plane (will be normalized)
    Real radius = 1.0;
    Real halflength = 1.0;

    GeometricShapeCylinder cylinder(center, axis, radius, halflength);

    // Point at center should be inside
    Vec3d point_center(0.0, 0.0, 0.0);
    EXPECT_TRUE(cylinder.checkContain(point_center));

    // Point along axis should be inside
    Vec3d normalized_axis = axis.normalized();
    Vec3d point_on_axis = 0.5 * normalized_axis;
    EXPECT_TRUE(cylinder.checkContain(point_on_axis));
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
