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
 * @file 	sph_solver.h
 * @brief 	Here we define the execution policy relevant to parallel computing.
 * @details This analog of the standard library on the same functions.
 * @author	Xiangyu Hu and  Fabien Pean
 */

#ifndef SPH_SOLVER_H
#define SPH_SOLVER_H

#include "particle_method_container.h"

namespace SPH
{
class TimeStepper
{
  public:
    TimeStepper(SPHSystem &sph_system, Real end_time, Real start_time = 0.0)
        : end_time_(end_time), global_dt_(0.0)
    {
        sv_physical_time_ = sph_system.getSystemVariableByName<Real>("PhysicalTime");
        sv_physical_time_->setValue(start_time);
    };
    ~TimeStepper() {};

    class TriggerByInterval
    {
      public:
        TriggerByInterval(Real initial_interval)
            : present_time_(0.0), interval_(initial_interval) {};

        bool operator()(BaseDynamics<Real> &interval_evaluator)
        {
            if (present_time_ > interval_)
            {
                present_time_ = 0.0;
                interval_ = interval_evaluator.exec();
                return true;
            }
            return false;
        };

        bool operator()()
        {
            if (present_time_ > interval_)
            {
                present_time_ = 0.0;
                return true;
            }
            return false;
        };

        Real getInterval() const
        {
            return interval_;
        };

        void incrementPresentTime(Real dt)
        {
            present_time_ += dt;
        };

      private:
        Real present_time_, interval_;
    };

    Real incrementPhysicalTime(BaseDynamics<Real> &step_evaluator)
    {
        global_dt_ = step_evaluator.exec();
        sv_physical_time_->incrementValue(global_dt_);
        for (auto &interval_executor : interval_executers_)
        {
            interval_executor->incrementPresentTime(global_dt_);
        }
        return global_dt_;
    };

    template <class Integrator>
    void integrateMatchedTimeInterval( // designed to avoid too small last step
        Real interval, BaseDynamics<Real> &step_evaluator, const Integrator &integrator)
    {
        Real integrated_time_ = 0.0;
        Real dt = step_evaluator.exec();

        while (interval - integrated_time_ > 1.5 * dt)
        {
            integrator(dt);
            dt = step_evaluator.exec();
            integrated_time_ += dt;
        }

        if (interval - integrated_time_ > dt)
        {
            Real final_dt = 0.5 * (interval - integrated_time_);
            integrator(final_dt);
            integrator(final_dt);
        }
        else
        {
            integrator(interval - integrated_time_);
        }
    };

    TriggerByInterval &addTriggerByInterval(Real initial_interval)
    {

        TriggerByInterval *interval_executor =
            execution_by_interval_keeper_.createPtr<TriggerByInterval>(initial_interval);
        interval_executers_.push_back(interval_executor);
        return *interval_executor;
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
    UniquePtrsKeeper<TriggerByInterval> execution_by_interval_keeper_;

  protected:
    StdVec<TriggerByInterval *> interval_executers_;
    Real end_time_, start_time_;
    Real global_dt_;
    SingularVariable<Real> *sv_physical_time_;
};

class SPHSolver
{
    UniquePtrsKeeper<BaseMethodContainer> methods_keeper_;
    UniquePtrKeeper<TimeStepper> time_stepper_keeper_;

  public:
    SPHSolver(SPHSystem &sph_system) : sph_system_(sph_system) {};
    virtual ~SPHSolver() {};

    template <typename ExecutePolicy>
    auto &addParticleMethodContainer(const ExecutePolicy &ex_policy)
    {
        return *methods_keeper_.createPtr<ParticleMethodContainer<ExecutePolicy>>(ex_policy);
    };

    auto &defineTimeStepper(Real end_time, Real start_time = 0.0)
    {
        return *time_stepper_keeper_.createPtr<TimeStepper>(sph_system_, end_time, start_time);
    };

  protected:
    SPHSystem &sph_system_;
};
} // namespace SPH
#endif // SPH_SOLVER_H
