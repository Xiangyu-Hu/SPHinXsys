#include "sphinxsys.h"
#include "fluid_mixture.h"
#include <gtest/gtest.h>

using namespace SPH;

namespace
{
class TestWeaklyCompressibleMixture : public WeaklyCompressibleMixture
{
  public:
    using WeaklyCompressibleMixture::WeaklyCompressibleMixture;

    void initializeLocalParameters(BaseParticles *base_particles) override
    {
        MatterMaterial::initializeLocalParameters(base_particles);
        WeaklyCompressibleMixture::initializeLocalParameters(base_particles);
    }

    Real getPressure(Real rho) override
    {
        return c0_ * c0_ * (rho - ReferenceDensity());
    }

    Real DensityFromPressure(Real p) override
    {
        return ReferenceDensity() + p / (c0_ * c0_);
    }

    Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override
    {
        return c0_;
    }
};

class MixtureBlock : public ComplexShape
{
  public:
    explicit MixtureBlock(const std::string &shape_name)
        : ComplexShape(shape_name)
    {
        Vecd half_size = ZeroData<Vecd>::value;
        for (int d = 0; d < Dimensions; ++d)
        {
            half_size[d] = 0.5;
        }
        Transform translate_to_center(half_size);
        add<GeometricShapeBox>(Transform(translate_to_center), half_size);
    }
};
} // namespace

TEST(ReferenceDesnityUpdate, UpdatesReferenceDensityAndMassFromMassFractions)
{
    Vecd lower_bound = ZeroData<Vecd>::value;
    Vecd upper_bound = ZeroData<Vecd>::value;
    for (int d = 0; d < Dimensions; ++d)
    {
        lower_bound[d] = -0.1;
        upper_bound[d] = 1.1;
    }

    BoundingBoxd system_domain_bounds(lower_bound, upper_bound);
    SPHSystem sph_system(system_domain_bounds, 0.25);

    FluidBody mixture_body(sph_system, makeShared<MixtureBlock>("MixtureBody"));

    StdVec<std::pair<std::string, Real>> species_data{{"A", 1.0}, {"B", 2.0}};
    mixture_body.defineMatterMaterial<TestWeaklyCompressibleMixture>(species_data, 10.0);
    mixture_body.generateParticles<BaseParticles, Lattice>();

    BaseParticles &particles = mixture_body.getBaseParticles();
    Real *rho0 = particles.getVariableDataByName<Real>("ReferenceDensity");
    Real *mass = particles.getVariableDataByName<Real>("Mass");
    MultiEntryView<Real> mass_fraction = particles.getVariableByName<Real>("MassFraction")->getMultiEntryView();

    const UnsignedInt total_particles = particles.TotalRealParticles();
    ASSERT_GT(total_particles, 0u);

    constexpr Real old_rho0 = 1.0;
    constexpr Real old_mass = 2.0;
    constexpr Real y_a = 0.25;
    constexpr Real y_b = 0.75;
    constexpr Real expected_rho0 = 1.0 / (y_a / 1.0 + y_b / 2.0);
    constexpr Real expected_mass = expected_rho0 * old_mass / old_rho0;

    for (UnsignedInt i = 0; i < total_particles; ++i)
    {
        rho0[i] = old_rho0;
        mass[i] = old_mass;
        mass_fraction[i][0] = y_a;
        mass_fraction[i][1] = y_b;
    }

    fluid_dynamics::ReferenceDesnityUpdate reference_density_update(mixture_body);
    fluid_dynamics::ReferenceDesnityUpdate::UpdateKernel update_kernel(execution::ParallelPolicy{}, reference_density_update);
    for (UnsignedInt i = 0; i < total_particles; ++i)
    {
        update_kernel.update(i);
    }

    for (UnsignedInt i = 0; i < total_particles; ++i)
    {
        EXPECT_NEAR(rho0[i], expected_rho0, 1.0e-12);
        EXPECT_NEAR(mass[i], expected_mass, 1.0e-12);
    }
}

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
