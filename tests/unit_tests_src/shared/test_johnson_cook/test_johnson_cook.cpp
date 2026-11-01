#include "sphinxsys.h"
#include <Eigen/Geometry>
#include <gtest/gtest.h>
#include <limits>

using namespace SPH;

static JohnsonCookParameters parameters()
{
    return {792e6, 510e6, Real(0.26), Real(0.014), 1, Real(1.03), 298, 1793, 477,
            Real(0.05), Real(3.44), Real(-2.12), Real(0.002), Real(0.61), true};
}

static Real tolerance(Real scale)
{
    return 128 * std::numeric_limits<Real>::epsilon() * std::max(Real(1), std::abs(scale));
}

TEST(Material, RejectsInvalidParameters)
{
    auto p = parameters();
    p.melting_temperature = p.reference_temperature;
    EXPECT_THROW(JohnsonCookSolid(7830, 200e9, 0.29, p), std::invalid_argument);
    p = parameters();
    p.specific_heat = 0;
    EXPECT_THROW(JohnsonCookSolid(7830, 200e9, 0.29, p), std::invalid_argument);
    EXPECT_THROW(JohnsonCookSolid(7830, 200e9, 0.5, parameters()), std::invalid_argument);
}

TEST(Material, ElasticShearAndDensityPressure)
{
    JohnsonCookSolid material(7830, 200e9, 0.29, parameters());
    auto state = material.InitialState();
    Mat3d increment = Mat3d::Zero();
    increment(0, 1) = increment(1, 0) = 1e-6;
    material.Integrate(increment, Real(1.0001), 1e-6, state);
    const Real shear = 200e9 / (1 + Real(0.29)) * Real(1e-6);
    const Real mean = -200e9 / (3 * (1 - 2 * Real(0.29))) * (1 / Real(1.0001) - 1);
    EXPECT_NEAR(state.stress(0, 1), shear, tolerance(shear));
    EXPECT_NEAR(state.stress.trace() / 3, mean, tolerance(mean));
    EXPECT_EQ(state.plastic_strain, 0);
    EXPECT_EQ(state.temperature, parameters().reference_temperature);
}

TEST(Material, LinearHardeningReturnAndAdiabaticWork)
{
    auto p = parameters();
    p.hardening_exponent = 1;
    p.strain_rate_coefficient = 0;
    p.failure_enabled = false;
    JohnsonCookSolid material(7830, 200e9, 0.29, p);
    auto state = material.InitialState();
    Mat3d increment = Mat3d::Zero();
    increment(2, 2) = Real(0.02);
    material.Integrate(increment, 1, 1e-6, state);
    const Real shear_modulus = 200e9 / (2 * (1 + Real(0.29)));
    const Real exact_increment = (2 * shear_modulus * Real(0.02) - p.yield_stress) / (3 * shear_modulus + p.hardening_coefficient);
    const Real exact_work = (p.yield_stress + p.hardening_coefficient * exact_increment) * exact_increment;
    EXPECT_NEAR(state.plastic_strain, exact_increment, tolerance(exact_increment));
    EXPECT_NEAR(state.temperature, p.reference_temperature + exact_work / (7830 * p.specific_heat), tolerance(state.temperature));
    EXPECT_NEAR(state.plastic_work_specific, exact_work / 7830, tolerance(state.plastic_work_specific));
}

TEST(Material, ThermalSofteningRemovesShearStrengthAtMelting)
{
    auto p = parameters();
    p.failure_enabled = false;
    JohnsonCookSolid material(7830, 200e9, 0.29, p);
    auto state = material.InitialState();
    state.temperature = p.melting_temperature;
    Mat3d increment = Mat3d::Zero();
    increment(0, 1) = increment(1, 0) = Real(0.01);
    material.Integrate(increment, 1, 1e-6, state);
    EXPECT_NEAR(state.stress.norm(), 0, tolerance(1));
    EXPECT_GT(state.plastic_strain, 0);
}

TEST(Material, FailureUsesTensionPositiveTriaxialityAndPreviousTemperature)
{
    auto p = parameters();
    JohnsonCookSolid material(7830, 200e9, 0.29, p);
    EXPECT_LT(material.FractureStrain(Real(1.0 / 3.0), 1, 298), material.FractureStrain(Real(-1.0 / 3.0), 1, 298));
    auto state = material.InitialState();
    Mat3d increment = Mat3d::Zero();
    increment(2, 2) = Real(0.02);
    const Real dt = 1e-6;
    material.Integrate(increment, 1, dt, state);
    const Real rate = Real(2.0 / 3.0) * increment(2, 2) / dt;
    const Real failure_strain = material.FractureStrain(0, rate, p.reference_temperature);
    EXPECT_NEAR(state.damage, state.plastic_strain / failure_strain, tolerance(state.damage));
    EXPECT_GT(state.temperature, p.reference_temperature);
}

TEST(Material, CompleteFailureRemovesHydrostaticAndDeviatoricStress)
{
    auto p = parameters();
    p.failure_D1 = Real(1e-8);
    p.failure_D2 = 0;
    JohnsonCookSolid material(7830, 200e9, 0.29, p);
    auto state = material.InitialState();
    Mat3d increment = Mat3d::Zero();
    increment(2, 2) = Real(0.02);
    material.Integrate(increment, Real(0.99), 1e-6, state);
    EXPECT_FALSE(state.active);
    EXPECT_EQ(state.damage, 1);
    EXPECT_EQ(state.stress.norm(), 0);
    const Real plastic_strain = state.plastic_strain;
    const Real temperature = state.temperature;
    material.Integrate(increment, 1, 1e-6, state);
    EXPECT_EQ(state.plastic_strain, plastic_strain);
    EXPECT_EQ(state.temperature, temperature);
    EXPECT_EQ(state.stress.norm(), 0);
}

TEST(Material, ConstitutiveResponseIsRotationCovariant)
{
    auto p = parameters();
    p.failure_enabled = false;
    JohnsonCookSolid material(7830, 200e9, 0.29, p);
    auto reference = material.InitialState();
    Mat3d increment;
    increment << Real(0.01), Real(0.003), 0, Real(0.003), Real(-0.004), 0, 0, 0, Real(-0.006);
    const Mat3d rotation = Eigen::AngleAxis<Real>(Real(0.7), Vec3d(1, 2, 3).normalized()).toRotationMatrix();
    auto rotated = reference;
    material.Integrate(increment, 1, 1e-6, reference);
    material.Integrate(rotation * increment * rotation.transpose(), 1, 1e-6, rotated);
    EXPECT_NEAR((rotated.stress - rotation * reference.stress * rotation.transpose()).norm(), 0, tolerance(reference.stress.norm()));
    EXPECT_NEAR(rotated.plastic_strain, reference.plastic_strain, tolerance(reference.plastic_strain));
    EXPECT_NEAR(rotated.temperature, reference.temperature, tolerance(reference.temperature));
}
