#include "johnson_cook_solid.hpp"
#include "sphinxsys.h"
#include <Eigen/Geometry>
#include <gtest/gtest.h>
#include <memory>
#include <type_traits>

using namespace SPH;

namespace
{
constexpr size_t point_count = 32;

JohnsonCookParameters Parameters()
{
    return {792e6, 510e6, Real(0.26), Real(0.014), 1, Real(1.03), 298, 1793, 477,
            Real(0.05), Real(3.44), Real(-2.12), Real(0.002), Real(0.61), true};
}

/** Real USM state on a SYCL build; the same checks also run on a native host. */
template <class T>
class SharedArray
{
  public:
    SharedArray()
    {
#if SPHINXSYS_USE_SYCL
        data_ = sycl::malloc_shared<T>(point_count, execution::execution_instance.getQueue());
        if (!data_)
            throw std::bad_alloc();
        std::uninitialized_value_construct_n(data_, point_count);
#else
        data_ = new T[point_count];
#endif
    }
    ~SharedArray()
    {
#if SPHINXSYS_USE_SYCL
        std::destroy_n(data_, point_count);
        sycl::free(data_, execution::execution_instance.getQueue());
#else
        delete[] data_;
#endif
    }
    SharedArray(const SharedArray &) = delete;
    SharedArray &operator=(const SharedArray &) = delete;
    T *data() { return data_; }

  private:
    T *data_;
};

template <class Function>
void Execute(Function function)
{
#if SPHINXSYS_USE_SYCL
    auto &queue = execution::execution_instance.getQueue();
    queue.parallel_for(sycl::range<1>(point_count), [=](sycl::id<1> index)
                       { function(size_t(index[0])); });
    queue.wait_and_throw();
#else
    for (size_t i = 0; i != point_count; ++i)
        function(i);
#endif
}

Matd Deformation(int step, size_t i)
{
    const Real fraction = Real(step <= 24 ? step : 48 - step) / 24;
    Matd F = Matd::Identity();
    F(0, 1) = (i % 2 ? Real(0.12) : Real(0.02)) * fraction;
    F(2, 2) = 1 + Real(0.015) * fraction;
    return F;
}

void Compare(const JohnsonCookState &actual, const JohnsonCookState &expected)
{
    const Real relative = sizeof(Real) == sizeof(float) ? Real(5e-4) : Real(1e-10);
    EXPECT_EQ(actual.active, expected.active);
    EXPECT_LE((actual.stress - expected.stress).norm(), relative * std::max(Real(792e6), expected.stress.norm()));
    EXPECT_NEAR(actual.plastic_strain, expected.plastic_strain, relative * std::max(Real(0.01), expected.plastic_strain));
    EXPECT_NEAR(actual.temperature, expected.temperature, relative * expected.temperature);
    EXPECT_NEAR(actual.damage, expected.damage, relative * std::max(Real(0.01), expected.damage));
    EXPECT_NEAR(actual.plastic_work_specific, expected.plastic_work_specific,
                relative * std::max(Real(1), expected.plastic_work_specific));
}
} // namespace

TEST(JohnsonCookDevice, MixedLoadingUnloadingAndRateHistory)
{
    static_assert(std::is_trivially_copyable<JohnsonCookConstitutiveModel>::value, "Device model must be a value type");
    JohnsonCookSolid material(7830, 200e9, Real(0.29), Parameters());
    const auto model = material.ConstitutiveModel();
    SharedArray<JohnsonCookState> storage;
    auto *actual = storage.data();
    for (size_t i = 0; i != point_count; ++i)
        actual[i] = material.InitialState();
    Execute([=](size_t i)
            {
        Matd previous = Matd::Identity();
        for (int step = 1; step <= 48; ++step)
        {
            const Matd F = Deformation(step, i);
            model.UpdateStress(F, previous, i % 4 < 2 ? Real(1e-6) : Real(1e-4), actual[i]);
            previous = F;
        } });
    for (size_t i = 0; i != point_count; ++i)
    {
        auto expected = material.InitialState();
        Matd previous = Matd::Identity();
        for (int step = 1; step <= 48; ++step)
        {
            const Matd F = Deformation(step, i);
            ASSERT_EQ(model.UpdateStress(F, previous, i % 4 < 2 ? Real(1e-6) : Real(1e-4), expected), JohnsonCookIntegrationStatus::success);
            previous = F;
        }
        Compare(actual[i], expected);
    }
    EXPECT_GT(actual[1].plastic_strain, Real(0));
    EXPECT_GT(actual[1].temperature, Parameters().reference_temperature);
    EXPECT_GT(actual[1].damage, Real(0));
}

TEST(JohnsonCookDevice, ThermalSofteningCompleteFailureAndFrozenHistory)
{
    auto p = Parameters();
    p.failure_D1 = Real(1e-5);
    p.failure_D2 = 0;
    JohnsonCookSolid material(7830, 200e9, Real(0.29), p);
    const auto model = material.ConstitutiveModel();
    SharedArray<JohnsonCookState> storage;
    SharedArray<JohnsonCookState> frozen_storage;
    auto *actual = storage.data();
    auto *frozen = frozen_storage.data();
    for (size_t i = 0; i != point_count; ++i)
    {
        actual[i] = material.InitialState();
        if (i % 2)
            actual[i].temperature = p.melting_temperature;
    }
    Execute([=](size_t i)
            {
        Mat3d increment = Mat3d::Zero();
        increment(0, 1) = increment(1, 0) = Real(0.02);
        model.Integrate(increment, i % 2 ? Real(1) : Real(0.99), Real(1e-6), actual[i]);
        frozen[i] = actual[i];
        model.Integrate(increment, Real(1.01), Real(1e-6), actual[i]); });
    for (size_t i = 0; i != point_count; ++i)
    {
        auto expected = material.InitialState();
        if (i % 2)
            expected.temperature = p.melting_temperature;
        Mat3d increment = Mat3d::Zero();
        increment(0, 1) = increment(1, 0) = Real(0.02);
        material.Integrate(increment, i % 2 ? Real(1) : Real(0.99), Real(1e-6), expected);
        material.Integrate(increment, Real(1.01), Real(1e-6), expected);
        Compare(actual[i], expected);
        EXPECT_FALSE(actual[i].active);
        EXPECT_EQ(actual[i].damage, 1);
        EXPECT_EQ(actual[i].stress.norm(), 0);
        EXPECT_FALSE(frozen[i].active);
        EXPECT_EQ(actual[i].temperature, frozen[i].temperature);
        EXPECT_EQ(actual[i].plastic_strain, frozen[i].plastic_strain);
    }
}

TEST(JohnsonCookDevice, InvalidDeformationReportsStatusWithoutChangingHistory)
{
    JohnsonCookSolid material(7830, 200e9, Real(0.29), Parameters());
    const auto model = material.ConstitutiveModel();
    SharedArray<JohnsonCookState> storage;
    SharedArray<int> statuses;
    auto *states = storage.data();
    auto *status = statuses.data();
    for (size_t i = 0; i != point_count; ++i)
        states[i] = material.InitialState();
    Execute([=](size_t i)
            {
        Matd F = Matd::Identity();
        if (i % 2)
            F(0, 0) = -1;
        status[i] = static_cast<int>(model.UpdateStress(F, Matd::Identity(), i % 2 ? Real(1e-6) : Real(-1), states[i])); });
    for (size_t i = 0; i != point_count; ++i)
    {
        EXPECT_NE(status[i], 0);
        EXPECT_EQ(states[i].stress.norm(), 0);
        EXPECT_EQ(states[i].plastic_strain, 0);
        EXPECT_EQ(states[i].damage, 0);
        EXPECT_EQ(states[i].temperature, Parameters().reference_temperature);
    }
}
