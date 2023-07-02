#include "base_data_type.h"
#include "complex_shape.h"
#include "geometric_shape.h"
#include "transform_shape.h"
#include <gtest/gtest.h>

using namespace SPH;

Real DL = 5.366;                    /**< Tank length. */
Real DH = 5.366;                    /**< Tank height. */
Real particle_spacing_ref = 0.025;  /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Extending width for boundary conditions. */
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;

class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
Vec2d test_point(0.1, -0.025);
auto tolerance = []()
{ return 100 * Eps; }; // Invalid initialization order with static libraries on GCC 9.4 requires function

TEST(test_GeometricShapeBox, test_closest_point)
{
    TransformShape<GeometricShapeBox> inner_wall_box(Transform(inner_wall_translation), inner_wall_halfsize);

    EXPECT_LE((inner_wall_box.findClosestPoint(test_point) - Vec2d(0.1, 0.0)).cwiseAbs().maxCoeff(), tolerance());
}

TEST(test_Complex_GeometricShapeBox, test_contain)
{
    WallBoundary wall_boundary("WallBoundary");

    EXPECT_EQ(wall_boundary.checkContain(test_point), true);
}

TEST(test_Complex_GeometricShapeBox, test_closest_point)
{
    WallBoundary wall_boundary("WallBoundary");

    EXPECT_LE((wall_boundary.findClosestPoint(test_point) - Vec2d(0.1, 0.0)).cwiseAbs().maxCoeff(), tolerance());
}

TEST(test_Complex_GeometricShapeBox, test_normal_direction)
{
    WallBoundary wall_boundary("WallBoundary");

    EXPECT_LE((wall_boundary.findNormalDirection(test_point) - Vec2d(0.0, 1.0)).cwiseAbs().maxCoeff(), tolerance());
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
