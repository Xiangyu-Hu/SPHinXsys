// SPDX-License-Identifier: Apache-2.0
#include "sphinxsys.h"
#include "hjc_solid.hpp"
#include <gtest/gtest.h>
#include <array>
#include <limits>

using namespace SPH;

namespace
{
HJCParameters parameters()
{
    return {Real(2.417e10), Real(.29), Real(2.06), Real(.0013), Real(.866), Real(1.19e8), Real(8.2e6),
            Real(1), Real(.01), Real(5), Real(4e7), Real(.00124), Real(1.2e9), Real(.011),
            Real(.04), Real(1), Real(1.287e10), Real(1.631e10), Real(6.495e10)};
}

/** Run the actual material kernel through the same execution policy as CK examples. */
class PrescribedDeformation : public LocalDynamics
{
  public:
    explicit PrescribedDeformation(SPHBody &body)
        : LocalDynamics(body), material_(DynamicCast<HJCSolid>(this, body.getMatterMaterial())),
          deformation_(particles_->registerStateVariable<Matd>("PrescribedDeformation", IdentityMatrix<Matd>::value)),
          time_step_(particles_->registerStateVariable<Real>("PrescribedTimeStep", Real(1e-6))) {}

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &policy, EncloserType &encloser)
            : material_(policy, encloser.material_),
              deformation_(encloser.deformation_->DelegatedData(policy)),
              time_step_(encloser.time_step_->DelegatedData(policy)) {}

        void update(UnsignedInt i, Real = Real(0))
        {
            material_.UpdateStress(deformation_[i], i, time_step_[i]);
        }

      private:
        HJCSolid::ConstituteKernel material_;
        Matd *deformation_;
        Real *time_step_;
    };

  protected:
    HJCSolid &material_;
    DiscreteVariable<Matd> *deformation_;
    DiscreteVariable<Real> *time_step_;
};

constexpr std::array<const char *, 8> scalar_names{
    "HJCDamage", "HJCPlasticStrain", "HJCPlasticVolume", "HJCMaximumCompression",
    "Pressure", "VonMisesStress", "HJCAcousticModulus", "Density"};

struct Snapshot
{
    Mat3d stress;
    Matd previous;
    std::array<Real, scalar_names.size()> values;
    int status;
};

template <class Policy>
std::vector<Snapshot> snapshots(BaseParticles &particles, const Policy &policy)
{
    auto *stress = particles.getVariableByName<Mat3d>("StressCauchy");
    auto *previous = particles.getVariableByName<Matd>("HJCPreviousDeformation");
    auto *status = particles.getVariableByName<int>("HJCIntegrationStatus");
    stress->prepareForOutput(policy);
    previous->prepareForOutput(policy);
    status->prepareForOutput(policy);
    std::array<Real *, scalar_names.size()> scalars;
    for (size_t k = 0; k != scalar_names.size(); ++k)
    {
        auto *variable = particles.getVariableByName<Real>(scalar_names[k]);
        variable->prepareForOutput(policy);
        scalars[k] = variable->Data();
    }
    std::vector<Snapshot> result(particles.TotalRealParticles());
    for (size_t i = 0; i != result.size(); ++i)
    {
        result[i].stress = stress->Data()[i];
        result[i].previous = previous->Data()[i];
        result[i].status = status->Data()[i];
        for (size_t k = 0; k != scalar_names.size(); ++k)
            result[i].values[k] = scalars[k][i];
    }
    return result;
}

class HJCKernelTest : public ::testing::Test
{
  protected:
    SPHSystem system_{BoundingBoxd(Vecd::Constant(-1), Vecd::Constant(1)), Real(.5)};
    SolidBody device_body_{system_, makeShared<GeometricShapeBox>(
        BoundingBoxd(Vecd::Constant(-.5), Vecd::Constant(.5)), "KernelPoints")};
    SolidBody reference_body_{system_, makeShared<GeometricShapeBox>(
        BoundingBoxd(Vecd::Constant(-.5), Vecd::Constant(.5)), "ReferencePoints")};
    HJCSolid &device_material_ = device_body_.defineMatterMaterial<HJCSolid>(2700, parameters());
    HJCSolid &reference_material_ = reference_body_.defineMatterMaterial<HJCSolid>(2700, parameters());
    BaseParticles &device_particles_ = device_body_.generateParticles<BaseParticles, Lattice>();
    BaseParticles &reference_particles_ = reference_body_.generateParticles<BaseParticles, Lattice>();
    StateDynamics<MainExecutionPolicy, PrescribedDeformation> update_{device_body_};

    void SetUp() override
    {
        ASSERT_GE(device_particles_.TotalRealParticles(), 8u);
        ASSERT_EQ(device_particles_.TotalRealParticles(), reference_particles_.TotalRealParticles());
    }

    template <class Deformation, class TimeStep>
    void advance(Deformation deformation, TimeStep time_step, bool update_reference = true)
    {
        auto *input = device_particles_.getVariableByName<Matd>("PrescribedDeformation");
        auto *dt = device_particles_.getVariableByName<Real>("PrescribedTimeStep");
        for (size_t i = 0; i != device_particles_.TotalRealParticles(); ++i)
        {
            input->Data()[i] = deformation(i);
            dt->Data()[i] = time_step(i);
            if (update_reference)
                reference_material_.UpdateStress(input->Data()[i], i, dt->Data()[i]);
        }
        input->finalizeLoadIn(par_ck);
        dt->finalizeLoadIn(par_ck);
        update_.exec();
    }

    void compare_reference()
    {
        const auto actual = snapshots(device_particles_, par_ck);
        const auto expected = snapshots(reference_particles_, par_host);
        // A stress tolerance scaled by fc remains meaningful near zero stress.
        const Real relative = sizeof(Real) == sizeof(float) ? Real(2e-3) : Real(1e-9);
        for (size_t i = 0; i != actual.size(); ++i)
        {
            SCOPED_TRACE(i);
            ASSERT_EQ(actual[i].status, static_cast<int>(HJCIntegrationStatus::success));
            ASSERT_TRUE(actual[i].stress.allFinite());
            EXPECT_LE((actual[i].stress - expected[i].stress).norm(),
                      relative * std::max(parameters().compressive_strength, expected[i].stress.norm()));
            EXPECT_EQ((actual[i].previous - expected[i].previous).norm(), Real(0));
            for (size_t k = 0; k != scalar_names.size(); ++k)
            {
                SCOPED_TRACE(scalar_names[k]);
                EXPECT_TRUE(std::isfinite(actual[i].values[k]));
                Real floor = k < 4 ? Real(.01) : k == 7 ? Real(2700) : parameters().compressive_strength;
                EXPECT_NEAR(actual[i].values[k], expected[i].values[k],
                            relative * std::max(floor, std::abs(expected[i].values[k])));
            }
        }
    }
};

TEST_F(HJCKernelTest, LoadingUnloadingAndRateDependence)
{
    std::vector<Snapshot> previous = snapshots(device_particles_, par_ck);
    // Different points exercise dense crushing, shear damage, tension and two rates.
    for (int step = 1; step <= 48; ++step)
    {
        Real fraction = Real(step <= 24 ? step : 48 - step) / Real(24);
        advance([&](size_t i)
        {
            Matd F = Matd::Identity();
            if (i % 4 == 0)
                F *= std::pow(Real(1) + Real(.15) * fraction, Real(-1) / Real(3));
            else if (i % 4 == 1 || i % 4 == 3)
                F(0, 1) = Real(.008) * fraction;
            else
                F *= std::exp(Real(.001) * fraction);
            return F;
        }, [](size_t i) { return i % 4 == 3 ? Real(1e-8) : Real(1e-4); });
        compare_reference();
        auto current = snapshots(device_particles_, par_ck);
        for (size_t i = 0; i != current.size(); ++i)
        {
            EXPECT_GE(current[i].values[0], previous[i].values[0]);
            EXPECT_LE(current[i].values[0], Real(1));
            EXPECT_GE(current[i].values[1], previous[i].values[1]);
            EXPECT_GE(current[i].values[2], previous[i].values[2]);
        }
        if (step == 24)
        {
            EXPECT_GT(current[0].values[4], parameters().lock_pressure);
            EXPECT_NEAR(current[0].values[2], parameters().lock_strain, Real(1e-7));
            EXPECT_GT(current[1].values[0], Real(0));
            EXPECT_GT(current[3].values[5], current[1].values[5]);
            EXPECT_LE(-current[2].values[4], parameters().tensile_strength *
                      (Real(1) - current[2].values[0]) * Real(1.001));
        }
        previous = current;
    }
    EXPECT_GT(previous[0].values[2], Real(0));
    EXPECT_GT(previous[1].values[0], Real(0));
}

TEST_F(HJCKernelTest, RotationPreservesStressAndUndamagedHistory)
{
    auto initial_deformation = [](size_t i)
    {
        Matd F = Real(.9999) * Matd::Identity();
        if (i % 2)
        {
            F(0, 0) = Real(1.0001);
            F(1, 1) = Real(1.0001001);
            F(2, 2) = Real(.9998);
        }
        return F;
    };
    advance(initial_deformation, [](size_t) { return Real(1e-5); });
    const auto initial = snapshots(device_particles_, par_ck);
    Matd rotation = Matd::Identity();
    for (int step = 1; step <= 32; ++step)
    {
        rotation = Eigen::AngleAxis<Real>(Real(step) * Real(.06),
            Vec3d(Real(1), Real(2), Real(3)).normalized()).toRotationMatrix();
        advance([&](size_t i) { return Matd(rotation * initial_deformation(i)); },
                [](size_t) { return Real(1e-5); });
    }
    compare_reference();
    const auto actual = snapshots(device_particles_, par_ck);
    const Real tolerance = Real(128) * std::numeric_limits<Real>::epsilon() * device_material_.BulkModulus();
    for (size_t i = 0; i != actual.size(); ++i)
    {
        SCOPED_TRACE(i);
        EXPECT_LE((actual[i].stress - rotation * initial[i].stress * rotation.transpose()).norm(), tolerance);
        EXPECT_EQ(actual[i].values[0], Real(0));
        EXPECT_EQ(actual[i].values[1], Real(0));
        EXPECT_EQ(actual[i].values[2], Real(0));
    }
}

TEST_F(HJCKernelTest, InvalidUpdatesLeaveHistoryUnchanged)
{
    advance([](size_t)
    {
        Matd F = Real(.99) * Matd::Identity();
        F(0, 1) = Real(.004);
        return F;
    }, [](size_t) { return Real(1e-6); });
    auto *previous = device_particles_.getVariableByName<Matd>("HJCPreviousDeformation");
    previous->prepareForOutput(par_ck);
    previous->Data()[6] = Matd::Zero();
    previous->finalizeLoadIn(par_ck);
    const auto before = snapshots(device_particles_, par_ck);
    advance([](size_t i)
    {
        Matd F = Real(.98) * Matd::Identity();
        if (i == 0) F(0, 0) *= Real(-1);
        if (i == 1) F.setZero();
        if (i == 2) F(0, 0) = std::numeric_limits<Real>::quiet_NaN();
        return F;
    }, [](size_t i)
    {
        if (i == 3) return Real(-1e-6);
        if (i == 4) return std::numeric_limits<Real>::infinity();
        if (i == 5) return Real(0);
        return Real(1e-6);
    }, false);
    const auto after = snapshots(device_particles_, par_ck);
    for (size_t i = 0; i != 7; ++i)
    {
        SCOPED_TRACE(i);
        HJCIntegrationStatus expected = i < 3 ? HJCIntegrationStatus::invalid_deformation :
            i < 5 ? HJCIntegrationStatus::invalid_increment :
            i == 5 ? HJCIntegrationStatus::success : HJCIntegrationStatus::singular_increment;
        EXPECT_EQ(after[i].status, static_cast<int>(expected));
        EXPECT_EQ((after[i].stress - before[i].stress).norm(), Real(0));
        EXPECT_EQ((after[i].previous - before[i].previous).norm(), Real(0));
        EXPECT_EQ(after[i].values, before[i].values);
    }
    EXPECT_EQ(after[7].status, static_cast<int>(HJCIntegrationStatus::success));
    EXPECT_GT((after[7].stress - before[7].stress).norm(), Real(0));
}
} // namespace
