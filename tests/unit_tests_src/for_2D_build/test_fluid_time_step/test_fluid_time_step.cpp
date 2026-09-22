#include "sphinxsys.h"
#include <gtest/gtest.h>

using namespace SPH;
using namespace SPH::fluid_dynamics;

namespace
{
class FluidTimeStepTest : public ::testing::Test
{
  protected:
    SPHSystem system_;
    FluidBody fluid_;
    SolidBody wall_;
    std::unique_ptr<ContactRelation> contact_;
    Real *mass_;
    Vecd *velocity_, *force_, *force_prior_, *wall_acceleration_;

    FluidTimeStepTest()
        : system_(BoundingBoxd(Vecd(-1.0, -1.0), Vecd(1.0, 1.0)), 0.1),
          fluid_(system_, makeShared<GeometricShapeBox>(
                              BoundingBoxd(Vecd(-0.4, -0.2), Vecd(0.4, 0.2)), "Fluid")),
          wall_(system_, makeShared<GeometricShapeBox>(
                             BoundingBoxd(Vecd(0.4, -0.2), Vecd(0.6, 0.2)), "Wall")) {}

    void SetUp() override
    {
        fluid_.defineMatterMaterial<WeaklyCompressibleFluid>(1.0, 2.0);
        fluid_.addMaterialProperty<Viscosity>(1.0);
        fluid_.generateParticles<BaseParticles, Lattice>();
        wall_.defineMatterMaterial<Solid>();
        wall_.generateParticles<BaseParticles, Lattice>();
        BaseParticles &particles = fluid_.getBaseParticles();
        particles.registerStateVariable<Real>("Pressure");
        particles.registerSingleVariable<Real>("SurfaceTensionCoef", 1.0);
        mass_ = particles.getVariableDataByName<Real>("Mass");
        velocity_ = particles.registerStateVariableData<Vecd>("Velocity");
        force_ = particles.registerStateVariableData<Vecd>("Force");
        force_prior_ = particles.registerStateVariableData<Vecd>("ForcePrior");
        wall_acceleration_ = wall_.getBaseParticles().registerStateVariableData<Vecd>("Acceleration");
        for (size_t i = 0; i < particles.TotalRealParticles(); ++i)
            velocity_[i] = Vecd::Zero();
        setFluidAcceleration(Vecd::Zero());
        setWallAcceleration(Vecd::Zero());
        contact_ = std::make_unique<ContactRelation>(fluid_, RealBodyVector{&wall_});
        system_.initializeSystemCellLinkedLists();
        system_.initializeSystemConfigurations();
    }

    void setFluidAcceleration(const Vecd &acceleration)
    {
        for (size_t i = 0; i < fluid_.getBaseParticles().TotalRealParticles(); ++i)
        {
            force_[i] = mass_[i] * acceleration;
            force_prior_[i] = Vecd::Zero();
        }
    }

    void setWallAcceleration(const Vecd &acceleration)
    {
        for (size_t i = 0; i < wall_.getBaseParticles().TotalRealParticles(); ++i)
            wall_acceleration_[i] = acceleration;
    }
};

TEST_F(FluidTimeStepTest, DefaultOuterAccelerationCriterionIsPreserved)
{
    ReduceDynamics<AdvectionTimeStep> legacy(fluid_, 1.0);
    ReduceDynamics<AdvectionTimeStepWithoutAcceleration> velocity_only(fluid_, 1.0);
    const Real initial_step = legacy.exec();
    EXPECT_DOUBLE_EQ(initial_step, velocity_only.exec());

    setFluidAcceleration(Vecd(10000.0, 0.0));
    const Real accelerated_step = legacy.exec();
    EXPECT_LT(accelerated_step, initial_step);
    EXPECT_DOUBLE_EQ(initial_step, velocity_only.exec());
    setFluidAcceleration(Vecd(40000.0, 0.0));
    EXPECT_NEAR(legacy.exec(), accelerated_step / 2.0, 1.0e-12);
}

TEST_F(FluidTimeStepTest, AcousticAccelerationIsOptInAndUsesNetForcePerMass)
{
    ReduceDynamics<AcousticTimeStep> legacy(fluid_);
    ReduceDynamics<AcousticTimeStepWithAcceleration> accelerated(fluid_);
    const Real initial_step = legacy.exec();
    EXPECT_DOUBLE_EQ(initial_step, accelerated.exec());

    setFluidAcceleration(Vecd(10000.0, 0.0));
    const Real accelerated_step = accelerated.exec();
    EXPECT_LT(accelerated_step, initial_step);
    EXPECT_DOUBLE_EQ(initial_step, legacy.exec());
    for (size_t i = 0; i < fluid_.getBaseParticles().TotalRealParticles(); ++i)
    {
        // A common mass/force rescaling must not change acceleration.
        mass_[i] *= 3.0;
        force_[i] *= 3.0;
    }
    EXPECT_NEAR(accelerated.exec(), accelerated_step, 1.0e-12);
    for (size_t i = 0; i < fluid_.getBaseParticles().TotalRealParticles(); ++i)
        force_prior_[i] = -force_[i];
    EXPECT_DOUBLE_EQ(initial_step, accelerated.exec());
}

TEST_F(FluidTimeStepTest, ViscousBoundSurvivesTheExplicitAdvectionOptIn)
{
    ReduceDynamics<AdvectionTimeStepWithoutAcceleration> velocity_only(fluid_, 1.0);
    ReduceDynamics<AdvectionViscousTimeStep> legacy_viscous(fluid_, 1.0);
    ReduceDynamics<AdvectionViscousTimeStepWithoutAcceleration> opt_in_viscous(fluid_, 1.0);
    const Real viscous_step = legacy_viscous.exec();
    EXPECT_LT(viscous_step, velocity_only.exec());
    EXPECT_DOUBLE_EQ(viscous_step, opt_in_viscous.exec());
    setFluidAcceleration(Vecd(10000.0, 0.0));
    EXPECT_LT(legacy_viscous.exec(), viscous_step);
    EXPECT_DOUBLE_EQ(viscous_step, opt_in_viscous.exec());
}

TEST_F(FluidTimeStepTest, SurfaceTensionDefaultDoesNotInheritTheOptIn)
{
    ReduceDynamics<SurfaceTensionTimeStep> capillary(fluid_);
    const Real initial_step = capillary.exec();
    setFluidAcceleration(Vecd(10000.0, 0.0));
    EXPECT_DOUBLE_EQ(initial_step, capillary.exec());
}

TEST_F(FluidTimeStepTest, WallBoundUsesApproachingRelativeAccelerationOnly)
{
    ReduceDynamics<WallAccelerationTimeStep> wall_step(*contact_);
    WallAccelerationTimeStep local_bound(*contact_);
    setWallAcceleration(Vecd(-100.0, 0.0));
    const Real approaching_step = wall_step.exec();
    // The wall Riemann extension uses ForcePrior, not the fluid pressure force.
    setFluidAcceleration(Vecd(10000.0, 0.0));
    EXPECT_DOUBLE_EQ(wall_step.exec(), approaching_step);
    setFluidAcceleration(Vecd::Zero());
    size_t contacted = 0;
    size_t noncontacted = 0;
    for (size_t i = 0; i < fluid_.getBaseParticles().TotalRealParticles(); ++i)
    {
        if (contact_->contact_configuration_[0][i].current_size_ > 0)
        {
            ++contacted;
            EXPECT_GT(local_bound.reduce(i), 0.0);
        }
        else
        {
            ++noncontacted;
            EXPECT_DOUBLE_EQ(local_bound.reduce(i), 0.0);
        }
    }
    ASSERT_GT(contacted, 0u);
    ASSERT_GT(noncontacted, 0u);
    setWallAcceleration(Vecd(-400.0, 0.0));
    EXPECT_NEAR(wall_step.exec(), approaching_step / 2.0, 1.0e-12);

    // Equal common acceleration leaves no wall-relative pressure acceleration.
    for (size_t i = 0; i < fluid_.getBaseParticles().TotalRealParticles(); ++i)
        force_prior_[i] = mass_[i] * Vecd(-400.0, 0.0);
    for (size_t i = 0; i < fluid_.getBaseParticles().TotalRealParticles(); ++i)
        EXPECT_DOUBLE_EQ(local_bound.reduce(i), 0.0);

    setFluidAcceleration(Vecd::Zero());
    setWallAcceleration(Vecd(100.0, 0.0));
    for (size_t i = 0; i < fluid_.getBaseParticles().TotalRealParticles(); ++i)
        EXPECT_DOUBLE_EQ(local_bound.reduce(i), 0.0);
}

TEST_F(FluidTimeStepTest, WallBoundIncludesEveryContactBody)
{
    SolidBody second_wall(system_, makeShared<GeometricShapeBox>(
                                       BoundingBoxd(Vecd(-0.6, -0.2), Vecd(-0.4, 0.2)), "SecondWall"));
    second_wall.defineMatterMaterial<Solid>();
    second_wall.generateParticles<BaseParticles, Lattice>();
    Vecd *second_acceleration =
        second_wall.getBaseParticles().registerStateVariableData<Vecd>("Acceleration");
    for (size_t i = 0; i < second_wall.getBaseParticles().TotalRealParticles(); ++i)
        second_acceleration[i] = Vecd(100.0, 0.0);
    // Right-hand wall recedes; only the second, left-hand wall approaches.
    setWallAcceleration(Vecd(100.0, 0.0));
    ContactRelation both_walls(fluid_, {&wall_, &second_wall});
    ContactRelation second_only(fluid_, {&second_wall});
    system_.initializeSystemCellLinkedLists();
    system_.initializeSystemConfigurations();
    ReduceDynamics<WallAccelerationTimeStep> first_wall_step(*contact_);
    ReduceDynamics<WallAccelerationTimeStep> both_walls_step(both_walls);
    ReduceDynamics<WallAccelerationTimeStep> second_wall_step(second_only);
    EXPECT_LT(both_walls_step.exec(), first_wall_step.exec());
    EXPECT_DOUBLE_EQ(both_walls_step.exec(), second_wall_step.exec());
}
} // namespace
