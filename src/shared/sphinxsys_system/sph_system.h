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
 * @file sph_system.h
 * @brief The SPH_System managing objects in the system level.
 * @details Note that the system operation prefer these are application independent.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef SPH_SYSTEM_H
#define SPH_SYSTEM_H

#include "base_body.h"
#include "execution_policy.h"
#include "simulation_context.h"
#include "sphinxsys_containers.h"

namespace SPH
{
class SPHBody;
/**
 * @class SPHSystem
 * @brief The SPH system managing objects in the system level.
 */
class SPHSystem
{
    UniquePtrKeeper<SimulationContext> simulation_ctx_ptr_keeper_;

  public:
    SPHSystem() {};
    virtual ~SPHSystem() {};

    template <typename... Args>
    SimulationContext &initializeSimulationContext(Args &&...args)
    {
        SimulationContext *simulation_ctx =
            simulation_ctx_ptr_keeper_.createPtr<SimulationContext>(
                std::forward<Args>(args)...);
        simulation_ctx->registerContextVariable<Real>("PhysicalTime", 0.0);
        return *simulation_ctx;
    };

    SimulationContext &getSimulationContext();
    void addRealBody(SPHBody *body) { real_bodies_.push_back(body); };
    void addSPHBody(SPHBody *body) { sph_bodies_.push_back(body); };
    void addObservationBody(SPHBody *body) { observation_bodies_.push_back(body); };
    void addSolidBody(SolidBody *solid_body) { solid_bodies_.push_back(solid_body); };
    void initializeSystemCellLinkedLists();
    void initializeSystemConfigurations();
    Real getSmallestTimeStepAmongSolidBodies(Real CFL = 0.6);

  protected:
    SPHBodyVector sph_bodies_;         /**< All sph bodies. */
    SPHBodyVector real_bodies_;        /**< The bodies with inner particle configuration. */
    SPHBodyVector observation_bodies_; /**< The bodies without inner particle configuration. */
    SolidBodyVector solid_bodies_;     /**< The bodies with inner particle configuration and acoustic time steps . */
};
} // namespace SPH
#endif // SPH_SYSTEM_H
