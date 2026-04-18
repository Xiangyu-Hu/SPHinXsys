#include "geometric_shape.h"
#include "transform_geometry.h"
#include <gtest/gtest.h>

using namespace SPH;

auto tolerance = []()
{ return 100 * Eps; };

/** Helper: build a Transform that maps UnitX to the given axis direction with translation center */
Transform makeTransform(const Vec2d &center, const Vec2d &axis)
{
    Real angle = std::atan2(axis[1], axis[0]);
    return Transform(Rotation2d(angle), center);
}

// Test GeometricShapeCylinder basic construction and bounds in 2D
// In 2D, a cylinder is like a rectangle (infinite cylinder perpendicular to the plane)
TEST(test_GeometricShapeCylinder, test_findBounds)
{
    Vec2d axis(0.0, 1.0); // y-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(makeTransform(Vec2d::Zero(), axis), radius, halflength);
    BoundingBoxd box = cylinder.getBounds();

    // For a cylinder along y-axis in 2D: x should be [-radius, radius], y should be [-halflength, halflength]
    EXPECT_LE(std::abs(box.lower_[0] - (-radius)), tolerance());
    EXPECT_LE(std::abs(box.lower_[1] - (-halflength)), tolerance());
    EXPECT_LE(std::abs(box.upper_[0] - radius), tolerance());
    EXPECT_LE(std::abs(box.upper_[1] - halflength), tolerance());
}

// Test containment - point inside cylinder
TEST(test_GeometricShapeCylinder, test_contain_inside)
{
    Vec2d axis(0.0, 1.0); // y-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(makeTransform(Vec2d::Zero(), axis), radius, halflength);

    // Point clearly inside
    Vec2d point_inside(0.5, 1.0);
    EXPECT_TRUE(cylinder.checkContain(point_inside));

    // Point at center
    Vec2d point_center(0.0, 0.0);
    EXPECT_TRUE(cylinder.checkContain(point_center));
}

// Test containment - point outside cylinder
TEST(test_GeometricShapeCylinder, test_contain_outside)
{
    Vec2d axis(0.0, 1.0); // y-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(makeTransform(Vec2d::Zero(), axis), radius, halflength);

    // Point outside radially
    Vec2d point_outside_radial(2.0, 0.0);
    EXPECT_FALSE(cylinder.checkContain(point_outside_radial));

    // Point outside axially
    Vec2d point_outside_axial(0.0, 3.0);
    EXPECT_FALSE(cylinder.checkContain(point_outside_axial));

    // Point outside both
    Vec2d point_outside_both(2.0, 3.0);
    EXPECT_FALSE(cylinder.checkContain(point_outside_both));
}

// Test containment - point on boundary
TEST(test_GeometricShapeCylinder, test_contain_boundary)
{
    Vec2d axis(0.0, 1.0); // y-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(makeTransform(Vec2d::Zero(), axis), radius, halflength);

    // Point on cylindrical surface
    Vec2d point_on_surface(1.0, 1.0);
    EXPECT_TRUE(cylinder.checkContain(point_on_surface));

    // Point on end cap
    Vec2d point_on_cap(0.5, 2.0);
    EXPECT_TRUE(cylinder.checkContain(point_on_cap));
}

// Test closest point - for point inside cylinder
TEST(test_GeometricShapeCylinder, test_closest_point_inside)
{
    Vec2d axis(0.0, 1.0); // y-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(makeTransform(Vec2d::Zero(), axis), radius, halflength);

    // Point inside, closer to side
    Vec2d point_inside(0.9, 0.0);
    Vec2d closest = cylinder.findClosestPoint(point_inside);
    // Should project to cylindrical surface
    EXPECT_LE(std::abs(closest[0] - 1.0), tolerance());
    EXPECT_LE(std::abs(closest[1] - 0.0), tolerance());
}

// Test closest point - for point outside radially
TEST(test_GeometricShapeCylinder, test_closest_point_outside_radial)
{
    Vec2d axis(0.0, 1.0); // y-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(makeTransform(Vec2d::Zero(), axis), radius, halflength);

    // Point outside radially
    Vec2d point_outside(2.0, 1.0);
    Vec2d closest = cylinder.findClosestPoint(point_outside);
    // Should project to cylindrical surface at (1.0, 1.0)
    EXPECT_LE(std::abs(closest[0] - 1.0), tolerance());
    EXPECT_LE(std::abs(closest[1] - 1.0), tolerance());
}

// Test closest point - for point outside axially
TEST(test_GeometricShapeCylinder, test_closest_point_outside_axial)
{
    Vec2d axis(0.0, 1.0); // y-axis
    Real radius = 1.0;
    Real halflength = 2.0;

    GeometricShapeCylinder cylinder(makeTransform(Vec2d::Zero(), axis), radius, halflength);

    // Point outside axially
    Vec2d point_outside(0.5, 3.0);
    Vec2d closest = cylinder.findClosestPoint(point_outside);
    // Should project to top cap at (0.5, 2.0)
    EXPECT_LE(std::abs(closest[0] - 0.5), tolerance());
    EXPECT_LE(std::abs(closest[1] - 2.0), tolerance());
}

// Test with translated cylinder (axis along x, so no rotation needed)
TEST(test_GeometricShapeCylinder, test_translated_cylinder)
{
    Vec2d center(5.0, 3.0);
    Real radius = 0.5;
    Real halflength = 1.0;

    // Axis = UnitX: use translation-only Transform
    GeometricShapeCylinder cylinder(Transform(center), radius, halflength);

    // Point inside translated cylinder
    Vec2d point_inside(5.0, 3.0);
    EXPECT_TRUE(cylinder.checkContain(point_inside));

    // Point outside
    Vec2d point_outside(5.0, 5.0);
    EXPECT_FALSE(cylinder.checkContain(point_outside));
}

// Test with diagonal axis orientation
TEST(test_GeometricShapeCylinder, test_diagonal_axis)
{
    Vec2d axis(1.0, 1.0); // diagonal
    Real radius = 1.0;
    Real halflength = 1.0;

    GeometricShapeCylinder cylinder(makeTransform(Vec2d::Zero(), axis), radius, halflength);

    // Point at center should be inside
    Vec2d point_center(0.0, 0.0);
    EXPECT_TRUE(cylinder.checkContain(point_center));

    // Point along axis should be inside
    Vec2d normalized_axis = axis.normalized();
    Vec2d point_on_axis = 0.5 * normalized_axis;
    EXPECT_TRUE(cylinder.checkContain(point_on_axis));
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
