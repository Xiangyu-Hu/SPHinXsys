#ifndef NECKING_CASE_H
#define NECKING_CASE_H

#include "necking_geometry.h"
#include <vector>

namespace SPH
{
struct NeckingOptions
{
    Real spacing = 5e-5;
    Real end_time = 150e-6;
    int frames = 60;
    int threads = 8;
    bool state_recording = true;
    bool regression_test = false;
};

inline NeckingOptions ParseNeckingOptions(int argc, char **argv, bool *use_host = nullptr)
{
    NeckingOptions options;
    for (int i = 1; i < argc; ++i)
    {
        const std::string argument = argv[i];
        const auto value = [&]()
        { return std::stod(argument.substr(argument.find('=') + 1)); };
        if (argument == "--regression-test")
            options.regression_test = true;
        else if (argument == "--host" && use_host)
            *use_host = true;
        else if (argument.rfind("--spacing=", 0) == 0)
            options.spacing = value();
        else if (argument.rfind("--end-time=", 0) == 0)
            options.end_time = value();
        else if (argument.rfind("--frames=", 0) == 0)
            options.frames = int(value());
        else if (argument.rfind("--threads=", 0) == 0)
            options.threads = int(value());
        else if (argument == "--state_recording=false")
            options.state_recording = false;
        else if (argument == "--state_recording=true")
            options.state_recording = true;
        else
            throw std::invalid_argument("Unknown option: " + argument);
    }
    if (!(std::isfinite(options.spacing) && options.spacing > 0 && options.spacing <= 2e-4 &&
          std::isfinite(options.end_time) && options.end_time > 0 && options.frames > 0 && options.threads > 0))
        throw std::invalid_argument("Invalid necking options");
    if (options.regression_test &&
        (options.spacing != Real(0.0002) || options.end_time != Real(0.00008) || options.frames != 4))
        throw std::invalid_argument("Regression requires spacing=0.0002, end-time=0.00008 and frames=4");
    return options;
}

inline JohnsonCookParameters SteelParameters()
{
    return {792e6, 510e6, Real(0.26), Real(0.014), 1,
            Real(1.03), 298, 1793, 477, Real(0.05), Real(3.44), Real(-2.12), Real(0.002), Real(0.61), true};
}

/** Shared cached loading state; the same time is used by every constrained particle. */
struct GripMotion
{
    Real displacement = 0;
    Real velocity = 0;
    Real acceleration = 0;
    bool constrain_position = false;

    void SetTime(Real time, bool position)
    {
        constexpr Real speed = 10;
        constexpr Real ramp_time = 20e-6;
        constrain_position = position;
        if (time < ramp_time)
        {
            const Real phase = Pi * time / ramp_time;
            displacement = speed * (time / 2 - ramp_time * std::sin(phase) / (2 * Pi));
            velocity = speed * Real(0.5) * (1 - std::cos(phase));
            acceleration = speed * Pi * std::sin(phase) / (2 * ramp_time);
        }
        else
        {
            displacement = speed * (time - ramp_time / 2);
            velocity = speed;
            acceleration = 0;
        }
    }
};

class PrescribedNeckingGrips : public LocalDynamics
{
  public:
    PrescribedNeckingGrips(SPHBody &body, GripMotion &motion)
        : LocalDynamics(body), motion_(motion),
          position_(particles_->getVariableDataByName<Vecd>("Position")),
          velocity_(particles_->getVariableDataByName<Vecd>("Velocity")),
          initial_position_(particles_->getVariableDataByName<Vecd>("InitialPosition")) {}

    void update(size_t index_i, Real = 0)
    {
        if (std::abs(initial_position_[index_i][2]) < NeckingRod::grip_position)
            return;
        const Real side = initial_position_[index_i][2] > 0 ? Real(1) : Real(-1);
        velocity_[index_i] = Vecd(0, 0, side * motion_.velocity);
        if (motion_.constrain_position)
        {
            position_[index_i] = initial_position_[index_i];
            position_[index_i][2] += side * motion_.displacement;
        }
    }

  protected:
    GripMotion &motion_;
    Vecd *position_, *velocity_, *initial_position_;
};

/** Boundary reaction includes the prescribed grip inertia during the ramp. */
class GripReactions
{
  public:
    explicit GripReactions(BaseParticles &particles)
        : force_(particles.getVariableDataByName<Vecd>("Force")), right_mass_(0), left_mass_(0)
    {
        const Vecd *position = particles.getVariableDataByName<Vecd>("InitialPosition");
        const Real *mass = particles.getVariableDataByName<Real>("Mass");
        for (size_t i = 0; i != particles.TotalRealParticles(); ++i)
        {
            if (position[i][2] >= NeckingRod::grip_position)
            {
                right_.push_back(i);
                right_mass_ += mass[i];
            }
            else if (position[i][2] <= -NeckingRod::grip_position)
            {
                left_.push_back(i);
                left_mass_ += mass[i];
            }
        }
    }
    Vec2d Measure(Real acceleration) const
    {
        Real right = 0, left = 0;
        for (size_t i : right_)
            right -= force_[i][2];
        for (size_t i : left_)
            left -= force_[i][2];
        return Vec2d(right + right_mass_ * acceleration, left - left_mass_ * acceleration);
    }

  protected:
    Vecd *force_;
    std::vector<size_t> right_, left_;
    Real right_mass_, left_mass_;
};
} // namespace SPH
#endif // NECKING_CASE_H
