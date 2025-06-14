#include "sph_solver.h"

#include "sph_system.hpp"

namespace SPH
{
//=================================================================================================//
TimeStepper::TimeStepper(SPHSystem &sph_system, Real end_time, Real start_time)
    : end_time_(end_time), global_dt_(0.0)
{
    sv_physical_time_ = sph_system.getSystemVariableByName<Real>("PhysicalTime");
    sv_physical_time_->setValue(start_time);
}
//=================================================================================================//
TimeStepper::TriggerByPhysicalTime::
    TriggerByPhysicalTime(TimeStepper &time_stepper, Real trigger_time)
    : sv_physical_time_(time_stepper.sv_physical_time_),
      trigger_time_(trigger_time) {}
//=================================================================================================//
bool TimeStepper::TriggerByPhysicalTime::operator()()
{
    if (sv_physical_time_->getValue() > trigger_time_)
    {
        return true;
    }
    return false;
}
//=================================================================================================//
TimeStepper::TriggerByInterval::TriggerByInterval(Real initial_interval)
    : present_time_(0.0), interval_(initial_interval) {}
//=================================================================================================//
bool TimeStepper::TriggerByInterval::operator()(BaseDynamics<Real> &interval_evaluator)
{
    if (present_time_ > interval_)
    {
        present_time_ = 0.0;
        interval_ = interval_evaluator.exec();
        return true;
    }
    return false;
}
//=================================================================================================//
bool TimeStepper::TriggerByInterval::operator()()
{
    if (present_time_ > interval_)
    {
        present_time_ = 0.0;
        return true;
    }
    return false;
}
//=================================================================================================//
Real TimeStepper::TriggerByInterval::getInterval() const
{
    return interval_;
}
//=================================================================================================//
void TimeStepper::TriggerByInterval::incrementPresentTime(Real dt)
{
    present_time_ += dt;
}
//=================================================================================================//
Real TimeStepper::incrementPhysicalTime(BaseDynamics<Real> &step_evaluator)
{
    global_dt_ = step_evaluator.exec();
    sv_physical_time_->incrementValue(global_dt_);
    for (auto &interval_executor : interval_executers_)
    {
        interval_executor->incrementPresentTime(global_dt_);
    }
    return global_dt_;
}
//=================================================================================================//
void TimeStepper::integrateMatchedTimeInterval(
    Real interval, BaseDynamics<Real> &step_evaluator, BaseDynamics<void> &integrator)
{
    Real integrated_time_ = 0.0;
    Real dt = step_evaluator.exec();

    while (interval - integrated_time_ > 1.5 * dt)
    {
        integrator.exec(dt);
        dt = step_evaluator.exec();
        integrated_time_ += dt;
    }
    // designed to avoid too small last step
    if (interval - integrated_time_ > dt)
    {
        Real final_dt = 0.5 * (interval - integrated_time_);
        integrator.exec(final_dt);
        integrator.exec(final_dt);
    }
    else
    {
        integrator.exec(interval - integrated_time_);
    }
}
//=================================================================================================//
TimeStepper::TriggerByInterval &TimeStepper::addTriggerByInterval(Real initial_interval)
{

    TriggerByInterval *interval_executor =
        execution_by_interval_keeper_.createPtr<TriggerByInterval>(initial_interval);
    interval_executers_.push_back(interval_executor);
    return *interval_executor;
}
//=================================================================================================//
TimeStepper::TriggerByPhysicalTime &TimeStepper::addTriggerByPhysicalTime(Real trigger_time)
{

    TriggerByPhysicalTime *executor =
        execution_by_physical_time_keeper_.createPtr<TriggerByPhysicalTime>(*this, trigger_time);
    physical_time_executers_.push_back(executor);
    return *executor;
}
//=================================================================================================//
bool TimeStepper::isEndTime()
{
    return (sv_physical_time_->getValue() >= end_time_);
}
//=================================================================================================//
void TimeStepper::setPhysicalTime(Real time)
{
    sv_physical_time_->setValue(time);
}
//=================================================================================================//
Real TimeStepper::getPhysicalTime()
{
    return sv_physical_time_->getValue();
}
//=================================================================================================//
Real TimeStepper::getGlobalTimeStepSize()
{
    return global_dt_;
}
//=================================================================================================//
} // namespace SPH
