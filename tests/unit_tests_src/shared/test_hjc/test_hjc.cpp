// SPDX-License-Identifier: Apache-2.0
#include "sphinxsys.h"
#include <gtest/gtest.h>
#include <limits>

using namespace SPH;

static HJCParameters parameters()
{
    return {2.417e10, .29, 2.06, .0013, .866, 1.19e8, 8.2e6, 1., .01, 5.,
            4e7, .00124, 1.2e9, .011, .04, 1., 1.287e10, 1.631e10, 6.495e10};
}

// Keep the double-precision limits; account for rounding at the tested scale in float.
static Real tolerance(Real double_limit, Real scale, Real roundoff_factor = 8)
{
    return sizeof(Real) == sizeof(float)
               ? double_limit + roundoff_factor * std::numeric_limits<Real>::epsilon() * std::abs(scale)
               : double_limit;
}

TEST(Material, RejectsInvalidParameters)
{
    auto p = parameters();
    p.minimum_fracture_strain = 0;
    EXPECT_THROW(HJCSolid(2700, p), std::invalid_argument);
    p = parameters();
    p.K2 = std::numeric_limits<Real>::quiet_NaN();
    EXPECT_THROW(HJCSolid(2700, p), std::invalid_argument);
    EXPECT_THROW(HJCSolid(-1, parameters()), std::invalid_argument);
}

TEST(Material, VirginElasticResponse)
{
    HJCSolid material(2700, parameters());
    HJCState state;
    Mat3d de = Mat3d::Zero();
    de(0, 1) = de(1, 0) = 1e-6;
    material.Integrate(de, 1e-5, 1e-6, state);
    const Real shear_stress = 2 * parameters().shear_modulus * Real(1e-6);
    const Real pressure = parameters().crush_pressure / parameters().crush_strain * Real(1e-5);
    EXPECT_NEAR(state.stress(0, 1), shear_stress, tolerance(1e-7, shear_stress));
    EXPECT_NEAR(-state.stress.trace() / 3, pressure, tolerance(1e-7, pressure));
    EXPECT_EQ(state.damage, 0);
    EXPECT_EQ(state.plastic_strain, 0);
    EXPECT_EQ(state.plastic_volume, 0);
}

TEST(Material, MinimumFractureStrainIsAFloor)
{
    HJCSolid material(2700, parameters());
    EXPECT_EQ(material.FractureStrain(0), parameters().minimum_fracture_strain);
    HJCState state;
    Mat3d de = Mat3d::Zero();
    de(0, 1) = de(1, 0) = .001;
    material.Integrate(de, 0, 1e-6, state);
    EXPECT_GT(state.plastic_strain, 0);
    EXPECT_NEAR(state.damage, state.plastic_strain / parameters().minimum_fracture_strain,
                tolerance(1e-14, state.damage));
    for (int i = 0; i < 100; ++i)
        material.Integrate(de, 0, 1e-6, state);
    EXPECT_DOUBLE_EQ(state.damage, 1);
    EXPECT_NEAR(state.stress.norm(), 0, 1e-8);
    Real ep = state.plastic_strain;
    material.Integrate(de, 0, 1e-6, state);
    EXPECT_GT(state.plastic_strain, ep); // Plastic flow continues at zero cohesion.
}

TEST(Material, TensileStrengthFallsToZeroAtCutoff)
{
    HJCSolid material(2700, parameters());
    Real T = parameters().tensile_strength;
    Real s0 = material.YieldStrength(0, .3, 1);
    EXPECT_LT(material.YieldStrength(-T * .2, .3, 1), s0);
    EXPECT_NEAR(material.YieldStrength(-T * .7, .3, 1), 0, 1e-8);
    HJCState state;
    state.damage = .3;
    material.Integrate(Mat3d::Identity() * .001, -.003, 1e-6, state);
    EXPECT_NEAR(state.stress.trace() / 3, T * .7, 1e-8);
}

TEST(Material, ShearSofteningApproachesClosedForm)
{
    auto p = parameters();
    p.C = 0;
    HJCSolid material(2700, p);
    HJCState state;
    const Real gamma = .005;
    Mat3d de = Mat3d::Zero();
    const int steps = 10000;
    de(0, 1) = de(1, 0) = gamma / (2 * steps);
    for (int i = 0; i < steps; ++i)
        material.Integrate(de, 0, 1e-8, state);
    // Solve q = 3G*(gamma/sqrt(3)-ep) = A*fc*(1-ep/EFMIN).
    Real q0 = p.A * p.compressive_strength;
    Real exact_ep = (3 * p.shear_modulus * gamma / std::sqrt(3.) - q0) /
                    (3 * p.shear_modulus - q0 / p.minimum_fracture_strain);
    // Single-precision accumulated-roundoff budget: 4*sqrt(steps)*epsilon*gamma.
    const Real roundoff_factor = 4 * std::sqrt(Real(steps));
    EXPECT_NEAR(state.plastic_strain, exact_ep, tolerance(2e-8, gamma, roundoff_factor));
    EXPECT_NEAR(state.damage, exact_ep / p.minimum_fracture_strain,
                tolerance(2e-6, gamma / p.minimum_fracture_strain, roundoff_factor));
}

TEST(Material, RateSensitivityAndStrengthCap)
{
    HJCSolid material(2700, parameters());
    Real base = material.YieldStrength(1e8, .2, 1);
    EXPECT_NEAR(material.YieldStrength(1e8, .2, 1000) / base,
                1 + parameters().C * std::log(1000.), tolerance(1e-14, 1));
    EXPECT_EQ(material.YieldStrength(1e8, .2, .01), base);
    EXPECT_EQ(material.YieldStrength(1e14, 0, 1e6), parameters().compressive_strength * 5);
}

TEST(Material, CrushingUnloadingAndLockContinuity)
{
    auto p = parameters();
    HJCSolid material(2700, p);
    HJCState state;
    EXPECT_NEAR(material.Pressure(p.crush_strain, state), p.crush_pressure, tolerance(.01, p.crush_pressure));
    Real mu = (p.crush_strain + material.LockCompression()) / 2;
    const Real transition_pressure = (p.crush_pressure + p.lock_pressure) / 2;
    EXPECT_NEAR(material.Pressure(mu, state), transition_pressure, tolerance(.01, transition_pressure));
    EXPECT_NEAR(state.plastic_volume, p.lock_strain / 2, tolerance(1e-14, p.lock_strain / 2));
    EXPECT_NEAR(material.Pressure(state.plastic_volume, state), 0, .01);
    EXPECT_NEAR(material.Pressure(mu, state), transition_pressure, tolerance(.01, transition_pressure));
    EXPECT_NEAR(material.Pressure(material.LockCompression(), state), p.lock_pressure, tolerance(.01, p.lock_pressure));
    const Real below_lock = sizeof(Real) == sizeof(float)
                                ? std::nextafter(material.LockCompression(), Real(0))
                                : material.LockCompression() - 1e-10;
    EXPECT_NEAR(material.Pressure(below_lock, state), p.lock_pressure, tolerance(10, p.lock_pressure));
    EXPECT_EQ(state.plastic_volume, p.lock_strain);
    EXPECT_NEAR(material.Pressure(p.lock_strain, state), 0, .01);
    Real eta = .2;
    Real dense_mu = p.lock_strain + (1 + p.lock_strain) * eta;
    const Real dense_pressure = p.K1 * eta + p.K2 * eta * eta + p.K3 * eta * eta * eta;
    EXPECT_NEAR(material.Pressure(dense_mu, state), dense_pressure, tolerance(.01, dense_pressure));
    EXPECT_EQ(state.plastic_volume, p.lock_strain);
}

TEST(Material, CrushingAloneAccumulatesDamage)
{
    HJCSolid material(2700, parameters());
    HJCState state;
    Real mu = .01;
    material.Integrate(Mat3d::Identity() * (-std::log1p(mu) / 3), mu, 1e-6, state);
    Real pressure = -state.stress.trace() / 3;
    EXPECT_GT(state.plastic_volume, 0);
    EXPECT_EQ(state.plastic_strain, 0);
    EXPECT_NEAR(state.damage, state.plastic_volume / material.FractureStrain(pressure), 1e-13);
}

static HJCState cycle(int steps)
{
    HJCSolid material(2700, parameters());
    HJCState state;
    Mat3d previous = Mat3d::Zero();
    Real previous_damage = 0, previous_volume = 0;
    for (int i = 1; i <= steps; ++i)
    {
        Real t = Real(i) / steps;
        Real h = std::pow(std::sin(Pi * t), 2);
        Mat3d strain = -.01 * h * Mat3d::Identity();
        strain(0, 1) = strain(1, 0) = .008 * h;
        material.Integrate(strain - previous, std::expm1(-strain.trace()), 1e-4 / steps, state);
        EXPECT_TRUE(state.stress.allFinite());
        EXPECT_GE(state.damage, previous_damage);
        EXPECT_LE(state.damage, 1);
        EXPECT_GE(state.plastic_volume, previous_volume);
        previous = strain;
        previous_damage = state.damage;
        previous_volume = state.plastic_volume;
    }
    return state;
}

TEST(Material, MixedLoadingConvergesWithStepRefinement)
{
    auto coarse = cycle(500), fine = cycle(1000), reference = cycle(4000);
    EXPECT_LT((fine.stress - reference.stress).norm(), (coarse.stress - reference.stress).norm());
    EXPECT_LT(std::abs(fine.damage - reference.damage), .005);
}

TEST(Material, ParticleUpdateIsObjectiveAndStoresHistory)
{
    SPHSystem system(BoundingBoxd(Vecd::Constant(-1), Vecd::Constant(1)), .5);
    SolidBody body(system, makeShared<GeometricShapeBox>(BoundingBoxd(Vecd::Constant(-.5), Vecd::Constant(.5)), "Cube"));
    auto &material = body.defineMatterMaterial<HJCSolid>(2700, parameters());
    auto &particles = body.generateParticles<BaseParticles, Lattice>();
    Matd F = Matd::Identity();
    F(0, 1) = .0001;
    Matd initial = material.UpdateStress(F, 0, 1e-6);
    Matd rotation = Eigen::AngleAxis<Real>(.6, Vec3d::UnitZ()).toRotationMatrix();
    Matd rotated = material.UpdateStress(rotation * F, 0, 1e-6);
    // Rounding in a near-identity incremental stretch is amplified by the elastic moduli.
    EXPECT_LT((rotated - rotation * initial * rotation.transpose()).norm(),
              tolerance(.01, std::max(material.BulkModulus(), parameters().shear_modulus)));
    EXPECT_EQ(particles.getVariableDataByName<Real>("HJCDamage")[0], 0);
    EXPECT_THROW(material.UpdateStress(-F, 0, 1e-6), std::domain_error);
    EXPECT_THROW(material.StressPK1(F, 0), std::logic_error);
}
