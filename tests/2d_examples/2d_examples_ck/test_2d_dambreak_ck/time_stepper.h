/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	time_stepper.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Xiangyu Hu and  Fabien Pean
 */

#ifndef TIME_STEPPER_H
#define TIME_STEPPER_H

#include "sphinxsys_ck.h"

namespace SPH
{
class PhysicalTimeStepper
{
  public:
    PhysicalTimeStepper(SPHSystem &sph_system, Real end_time, Real start_time = 0.0)
        : end_time_(end_time), global_dt_(0.0)
    {
        sv_physical_time_ = sph_system.getSystemVariableByName<Real>("PhysicalTime");
        sv_physical_time_->setValue(start_time);
    };
    ~PhysicalTimeStepper() {};

    class ExecutionByInterval
    {
      public:
        ExecutionByInterval(PhysicalTimeStepper &physical_time_stepper)
            : physical_time_stepper_(physical_time_stepper),
              present_(0.0), interval_(-1.0) {};

        template <class Executor> // first function
        bool operator()(BaseDynamics<Real> &interval_evaluator, const Executor &executor)
        {
            present_ += physical_time_stepper_.getGlobalTimeStepSize();
            if (present_ > interval_)
            {
                executor();
                present_ = 0.0;
                interval_ = interval_evaluator.exec();
                return true;
            }
            return false;
        };

        template <class Executor> // second function
        bool operator()(Real fix_interval, const Executor &executor)
        {
            present_ += physical_time_stepper_.getGlobalTimeStepSize();
            if (present_ > interval_)
            {
                executor();
                present_ = 0.0;
                interval_ = fix_interval;
                return true;
            }
            return false;
        };

        template <class Executor> // depending on either previous function
        void operator()(bool do_execution, const Executor &executor)
        {
            if (do_execution)
                executor();
        };

        Real getInterval() const
        {
            return interval_;
        };

      private:
        PhysicalTimeStepper &physical_time_stepper_;
        Real present_, interval_;
    };

    template <class Integrator>
    void incrementPhysicalTimeWithIntegration(
        BaseDynamics<Real> &step_evaluator, const Integrator &integrator)
    {
        global_dt_ = step_evaluator.exec();
        integrator(global_dt_);
        sv_physical_time_->incrementValue(global_dt_);
    };

    template <class Integrator>
    void TimeIntervalIntegration( // designed to avoid too small last step
        Real interval, BaseDynamics<Real> &step_evaluator, const Integrator &integrator)
    {
        Real present = 0.0;
        Real dt = step_evaluator.exec();

        while (interval - present > 1.5 * dt)
        {
            integrator(dt);
            dt = step_evaluator.exec();
            present += dt;
        }

        if (interval - present > dt)
        {
            Real final_dt = 0.5 * (interval - present);
            integrator(final_dt);
            integrator(final_dt);
        }
        else
        {
            integrator(interval - present);
        }
    };

    ExecutionByInterval &addExecutionByInterval()
    {
        return *execution_by_interval_keeper_
                    .createPtr<ExecutionByInterval>(*this);
    };

    bool isEndTime()
    {
        return (sv_physical_time_->getValue() >= end_time_);
    };

    void setPhysicalTime(Real time)
    {
        sv_physical_time_->setValue(time);
    };

    Real getPhysicalTime()
    {
        return sv_physical_time_->getValue();
    };

    Real getGlobalTimeStepSize()
    {
        return global_dt_;
    };

  private:
    UniquePtrsKeeper<ExecutionByInterval> execution_by_interval_keeper_;
    Real end_time_, start_time_;
    Real global_dt_;
    SingularVariable<Real> *sv_physical_time_;
};
} // namespace SPH
#endif // TIME_STEPPER_H
