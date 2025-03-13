#include "sphinxsys.h"
#include <gtest/gtest.h>
using namespace SPH;

TEST(material_verification, sanity_check)
{
    constexpr double tolerance = std::numeric_limits<float>::epsilon();
    const auto I = SPH::Mat3d::Identity();
    constexpr double E = 1e4;
    constexpr double nu = 0.49;
    const SPH::Mat3d R = Eigen::AngleAxisd(M_PI / 6, SPH::Vec3d{1, 1, 0}.normalized()).toRotationMatrix();
    constexpr double stretch = 2.0;
    const SPH::Mat3d U = SPH::Vec3d{stretch, 1 / sqrt(stretch), 1 / sqrt(stretch)}.asDiagonal();
    double J = U.determinant();
    SPH::Mat3d F = R * U;
    SPH::Mat3d B = F * F.transpose();
    SPH::Mat3d e = 0.5 * (I - B.inverse());

    auto sanity_check = [&](auto &&material)
    {
        auto S = material.StressPK2(F, 0);
        auto s = material.StressCauchy(e, 0);
        return (s.isApprox(F * S * F.transpose() / J, tolerance));
    };
    EXPECT_TRUE(sanity_check(SPH::LinearElasticSolid{1, E, nu}));
    EXPECT_TRUE(sanity_check(SPH::SaintVenantKirchhoffSolid{1, E, nu}));
    EXPECT_TRUE(sanity_check(SPH::NeoHookeanSolidIncompressible{1, E, nu}));
    EXPECT_TRUE(sanity_check(SPH::NeoHookeanSolid{1, E, nu}));
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
