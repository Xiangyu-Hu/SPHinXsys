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
    TimeStepper(SPHSystem &sph_system, Real end_time, Real start_time = 0.0);
    ~TimeStepper() {};

    class TriggerByPhysicalTime
    {
      public:
        TriggerByPhysicalTime(TimeStepper &time_stepper, Real trigger_time);
        bool operator()();

      private:
        SingularVariable<Real> *sv_physical_time_;
        Real trigger_time_;
    };

    class TriggerByInterval
    {
      public:
        TriggerByInterval(Real initial_interval);
        bool operator()(BaseDynamics<Real> &interval_evaluator);
        bool operator()();
        Real getInterval() const;
        void incrementPresentTime(Real dt);

      private:
        Real present_time_, interval_;
    };

    Real incrementPhysicalTime(BaseDynamics<Real> &step_evaluator);

    void integrateMatchedTimeInterval(
        Real interval, BaseDynamics<Real> &step_evaluator, BaseDynamics<void> &integrator);

    TriggerByInterval &addTriggerByInterval(Real initial_interval);
    TriggerByPhysicalTime &addTriggerByPhysicalTime(Real trigger_time);

    bool isEndTime();
    void setPhysicalTime(Real time);
    Real getPhysicalTime();
    Real getGlobalTimeStepSize();

  private:
    UniquePtrsKeeper<TriggerByInterval> execution_by_interval_keeper_;
    UniquePtrsKeeper<TriggerByPhysicalTime> execution_by_physical_time_keeper_;

  protected:
    StdVec<TriggerByInterval *> interval_executers_;
    StdVec<TriggerByPhysicalTime *> physical_time_executers_;
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
