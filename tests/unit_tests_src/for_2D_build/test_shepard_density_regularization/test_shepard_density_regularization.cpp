#include "sphinxsys.h"
#include <gtest/gtest.h>

#include <algorithm>
#include <limits>
#include <memory>

using namespace SPH;

namespace
{
class DensityTestShape : public ComplexShape
{
  public:
    DensityTestShape(const std::string &name, const Vecd &center, const Vecd &half_size)
        : ComplexShape(name)
    {
        add<GeometricShapeBox>(Transform(center), half_size);
    }
};

/** A truncated fluid patch above a wall covering only part of its lower edge. */
class ShepardDensityTest : public ::testing::Test
{
  protected:
    ShepardDensityTest()
        : system_(BoundingBoxd(Vecd(-0.5, -0.5), Vecd(1.5, 1.0)), 0.1),
          fluid_(system_, makeShared<DensityTestShape>("Fluid", Vecd(0.4, 0.2), Vecd(0.4, 0.2))),
          wall_(system_, makeShared<DensityTestShape>("Wall", Vecd(0.2, -0.1), Vecd(0.2, 0.1)))
    {
        fluid_.defineMatterMaterial<WeaklyCompressibleFluid>(1000.0, 10.0);
        fluid_.generateParticles<BaseParticles, Lattice>();
        wall_.defineMatterMaterial<Solid>(1200.0);
        wall_.generateParticles<BaseParticles, Lattice>();
        inner_ = std::make_unique<InnerRelation>(fluid_);
        contact_ = std::make_unique<ContactRelation>(fluid_, RealBodyVector{&wall_});
        system_.initializeSystemCellLinkedLists();
        system_.initializeSystemConfigurations();

        BaseParticles &particles = fluid_.getBaseParticles();
        count_ = particles.TotalRealParticles();
        rho_ = particles.getVariableDataByName<Real>("Density");
        mass_ = particles.getVariableDataByName<Real>("Mass");
        volume_ = particles.getVariableDataByName<Real>("VolumetricMeasure");
    }

    void setNonuniformField()
    {
        for (size_t i = 0; i != count_; ++i)
        {
            rho_[i] = 800.0 + 100.0 * (i % 5);
            mass_[i] *= 0.7 + 0.2 * (i % 4);
            volume_[i] = mass_[i] / rho_[i];
        }
    }

    void expectMixedWallSupport()
    {
        size_t with_wall = 0;
        for (size_t i = 0; i != count_; ++i)
            with_wall += contact_->contact_configuration_[0][i].current_size_ != 0;
        EXPECT_GT(with_wall, 0u);
        EXPECT_LT(with_wall, count_);
    }

    SPHSystem system_;
    FluidBody fluid_;
    SolidBody wall_;
    std::unique_ptr<InnerRelation> inner_;
    std::unique_ptr<ContactRelation> contact_;
    size_t count_;
    Real *rho_, *mass_, *volume_;
    // Round-off allowance for a weighted sum of a few dozen terms, not a CFD tolerance.
    const Real density_tolerance_ = 128.0 * std::numeric_limits<Real>::epsilon() * 1200.0;
};

TEST_F(ShepardDensityTest, PreservesConstantFieldAcrossTruncatedWallSupport)
{
    expectMixedWallSupport();
    const Real constant_density = 1037.0; // Deliberately different from reference density.
    for (size_t i = 0; i != count_; ++i)
    {
        rho_[i] = constant_density;
        mass_[i] *= 0.7 + 0.2 * (i % 4);
        volume_[i] = mass_[i] / rho_[i];
    }
    InteractionWithUpdate<fluid_dynamics::ShepardDensityRegularizationWithWall>
        regularize(*inner_, *contact_);
    regularize.exec();
    for (size_t i = 0; i != count_; ++i)
        EXPECT_NEAR(rho_[i], constant_density, density_tolerance_);
}

TEST_F(ShepardDensityTest, SmoothsWithoutWallContactAndStaysWithinDensityBounds)
{
    setNonuniformField();
    const StdVec<Real> before(rho_, rho_ + count_);
    ContactRelation no_wall(fluid_, RealBodyVector{});
    no_wall.updateConfiguration();
    InteractionWithUpdate<fluid_dynamics::ShepardDensityRegularizationWithWall>
        regularize(*inner_, no_wall);
    regularize.exec();
    Real maximum_change = 0.0;
    for (size_t i = 0; i != count_; ++i)
    {
        EXPECT_GE(rho_[i], 800.0 - density_tolerance_);
        EXPECT_LE(rho_[i], 1200.0 + density_tolerance_);
        maximum_change = std::max(maximum_change, std::abs(rho_[i] - before[i]));
    }
    EXPECT_GT(maximum_change, density_tolerance_);
}

TEST_F(ShepardDensityTest, PreservesParticleMassAndUpdatesPositiveConsistentVolume)
{
    setNonuniformField();
    const StdVec<Real> before_mass(mass_, mass_ + count_);
    InteractionWithUpdate<fluid_dynamics::ShepardDensityRegularizationWithWall>
        regularize(*inner_, *contact_);
    regularize.exec();
    for (size_t i = 0; i != count_; ++i)
    {
        EXPECT_GE(rho_[i], 800.0 - density_tolerance_);
        EXPECT_LE(rho_[i], 1200.0 + density_tolerance_);
        EXPECT_EQ(mass_[i], before_mass[i]);
        EXPECT_GT(volume_[i], 0.0);
        EXPECT_NEAR(rho_[i] * volume_[i], mass_[i],
                    8.0 * std::numeric_limits<Real>::epsilon() * mass_[i]);
    }
}

TEST_F(ShepardDensityTest, PreservesRawDensitySummationForBoundaryConsumers)
{
    setNonuniformField();
    // Compute the actual legacy geometry sum without applying its density update.
    // Boundary velocity correction needs this raw field after regularization.
    const StdVec<Real> before_rho(rho_, rho_ + count_);
    const StdVec<Real> before_volume(volume_, volume_ + count_);
    fluid_dynamics::DensitySummationComplex raw_sum(*inner_, *contact_);
    for (size_t i = 0; i != count_; ++i)
        raw_sum.interaction(i);
    for (size_t i = 0; i != count_; ++i)
    {
        ASSERT_EQ(rho_[i], before_rho[i]);
        ASSERT_EQ(volume_[i], before_volume[i]);
    }
    Real *raw = fluid_.getBaseParticles().getVariableDataByName<Real>("DensitySummation");
    const StdVec<Real> before_raw(raw, raw + count_);
    InteractionWithUpdate<fluid_dynamics::ShepardDensityRegularizationWithWall>
        regularize(*inner_, *contact_);
    regularize.exec();
    Real maximum_difference = 0.0;
    for (size_t i = 0; i != count_; ++i)
    {
        EXPECT_EQ(raw[i], before_raw[i]);
        maximum_difference = std::max(maximum_difference, std::abs(raw[i] - rho_[i]));
    }
    // Avoid a vacuous preservation check: the raw and normalized fields differ.
    EXPECT_GT(maximum_difference, density_tolerance_);
}

TEST_F(ShepardDensityTest, UsesWallReferenceVolumeRatherThanCurrentWallCompression)
{
    setNonuniformField();
    const StdVec<Real> before_rho(rho_, rho_ + count_);
    const StdVec<Real> before_volume(volume_, volume_ + count_);
    InteractionWithUpdate<fluid_dynamics::ShepardDensityRegularizationWithWall>
        regularize(*inner_, *contact_);
    regularize.exec();
    const StdVec<Real> result(rho_, rho_ + count_);
    std::copy(before_rho.begin(), before_rho.end(), rho_);
    std::copy(before_volume.begin(), before_volume.end(), volume_);

    BaseParticles &wall_particles = wall_.getBaseParticles();
    Real *wall_rho = wall_particles.getVariableDataByName<Real>("Density");
    Real *wall_volume = wall_particles.getVariableDataByName<Real>("VolumetricMeasure");
    for (size_t i = 0; i != wall_particles.TotalRealParticles(); ++i)
    {
        // Vary current compression without changing mass, reference material or geometry.
        wall_rho[i] *= 0.5;
        wall_volume[i] *= 2.0;
    }
    regularize.exec();
    for (size_t i = 0; i != count_; ++i)
        EXPECT_NEAR(rho_[i], result[i], density_tolerance_);
}

TEST_F(ShepardDensityTest, IncludesBothWallsWithDistinctReferenceDensities)
{
    ASSERT_GE(count_, 2u);
    SolidBody second_wall(system_, makeShared<DensityTestShape>(
                                      "SecondWall", Vecd(0.6, -0.1), Vecd(0.2, 0.1)));
    second_wall.defineMatterMaterial<Solid>(2400.0);
    second_wall.generateParticles<BaseParticles, Lattice>();
    ContactRelation two_walls(fluid_, RealBodyVector{&wall_, &second_wall});
    system_.initializeSystemCellLinkedLists();
    system_.initializeSystemConfigurations();

    // An exactly calculable four-point quadrature: target, one fluid neighbor,
    // and one particle in each wall. Equal kernel weights cancel. Their volume
    // weights are 2, 1, 3 and 5; extended densities are 1000, 3000, 1000, 1000.
    // This fixture deliberately controls weights; the other tests use geometric
    // neighbor searches. It detects omitted contacts, incorrect wall rho0, and
    // accidental use of current wall volume without mirroring the production sum.
    rho_[0] = 1000.0;
    volume_[0] = 2.0;
    mass_[0] = 2000.0;
    rho_[1] = 3000.0;
    volume_[1] = 1.0;
    mass_[1] = 3000.0;
    BaseParticles &first_particles = wall_.getBaseParticles();
    BaseParticles &second_particles = second_wall.getBaseParticles();
    ASSERT_GT(first_particles.TotalRealParticles(), 0u);
    ASSERT_GT(second_particles.TotalRealParticles(), 0u);
    first_particles.getVariableDataByName<Real>("Mass")[0] = 3600.0;
    second_particles.getVariableDataByName<Real>("Mass")[0] = 12000.0;
    // Reference volumes are 3600/1200=3 and 12000/2400=5. Current volumes
    // intentionally differ; they must not become wall quadrature weights.
    first_particles.getVariableDataByName<Real>("Density")[0] = 1800.0;
    first_particles.getVariableDataByName<Real>("VolumetricMeasure")[0] = 2.0;
    second_particles.getVariableDataByName<Real>("Density")[0] = 3000.0;
    second_particles.getVariableDataByName<Real>("VolumetricMeasure")[0] = 4.0;

    const Real equal_weight = fluid_.getSPHAdaptation().getKernel()->W0(ZeroVecd);
    const auto set_single_neighbor = [equal_weight](Neighborhood &neighborhood, size_t index)
    {
        neighborhood.current_size_ = neighborhood.allocated_size_ = 1;
        neighborhood.j_.assign(1, index);
        neighborhood.W_ij_.assign(1, equal_weight);
        neighborhood.dW_ij_.assign(1, 0.0);
        neighborhood.r_ij_.assign(1, 0.0);
        neighborhood.e_ij_.assign(1, Vecd::Zero());
    };
    set_single_neighbor(inner_->inner_configuration_[0], 1);
    set_single_neighbor(two_walls.contact_configuration_[0][0], 0);
    set_single_neighbor(two_walls.contact_configuration_[1][0], 0);
    fluid_dynamics::ShepardDensityRegularizationWithWall regularize(*inner_, two_walls);
    Real *reconstructed = fluid_.getBaseParticles().getVariableDataByName<Real>("ShepardDensity");

    regularize.interaction(0);
    EXPECT_NEAR(reconstructed[0], 13000.0 / 11.0, density_tolerance_);
    // Removing either wall must change the result for a nonuniform fluid field.
    two_walls.contact_configuration_[1][0].current_size_ = 0;
    regularize.interaction(0);
    EXPECT_NEAR(reconstructed[0], 8000.0 / 6.0, density_tolerance_);
    two_walls.contact_configuration_[0][0].current_size_ = 0;
    two_walls.contact_configuration_[1][0].current_size_ = 1;
    regularize.interaction(0);
    EXPECT_NEAR(reconstructed[0], 10000.0 / 8.0, density_tolerance_);
    two_walls.contact_configuration_[1][0].current_size_ = 0;
    regularize.interaction(0);
    EXPECT_NEAR(reconstructed[0], 5000.0 / 3.0, density_tolerance_);
}

TEST_F(ShepardDensityTest, SeparatesInteractionFromUpdateAndParticleTraversalOrder)
{
    setNonuniformField();
    const StdVec<Real> before_rho(rho_, rho_ + count_);
    const StdVec<Real> before_volume(volume_, volume_ + count_);
    fluid_dynamics::ShepardDensityRegularizationWithWall staged(*inner_, *contact_);
    Real *reconstructed = fluid_.getBaseParticles().getVariableDataByName<Real>("ShepardDensity");
    for (size_t i = 0; i != count_; ++i)
        staged.interaction(i);
    const StdVec<Real> forward_result(reconstructed, reconstructed + count_);
    for (size_t i = 0; i != count_; ++i)
    {
        EXPECT_EQ(rho_[i], before_rho[i]);
        EXPECT_EQ(volume_[i], before_volume[i]);
    }
    for (size_t i = count_; i != 0; --i)
        staged.interaction(i - 1);
    for (size_t i = 0; i != count_; ++i)
    {
        EXPECT_NEAR(reconstructed[i], forward_result[i], density_tolerance_);
        staged.update(i);
        EXPECT_EQ(rho_[i], reconstructed[i]);
    }
    std::copy(before_rho.begin(), before_rho.end(), rho_);
    std::copy(before_volume.begin(), before_volume.end(), volume_);
    InteractionWithUpdate<fluid_dynamics::ShepardDensityRegularizationWithWall>
        regularize(*inner_, *contact_);
    regularize.exec();
    for (size_t i = 0; i != count_; ++i)
        EXPECT_NEAR(rho_[i], forward_result[i], density_tolerance_);
}
} // namespace
