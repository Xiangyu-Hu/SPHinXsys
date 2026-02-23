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
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file    sph_simulation.h
 * @brief   High-level user-facing API for setting up and running SPH simulations.
 *          Provides a fluent builder interface to configure domain, bodies, solver,
 *          and run the simulation with minimal boilerplate.
 * @author  Xiangyu Hu
 */

#ifndef SPH_SIMULATION_H
#define SPH_SIMULATION_H

#include "base_data_type_package.h"

#include <memory>
#include <string>
#include <vector>

namespace SPH
{
/**
 * @class FluidBlockBuilder
 * @brief Builder for configuring a fluid body in a 2D simulation.
 *        Supports a fluent interface: addFluidBlock("name").rectangle(L, H).material(rho, c)
 */
class FluidBlockBuilder
{
  public:
    explicit FluidBlockBuilder(const std::string &name);

    /** Define the fluid block as a rectangle starting at the origin. */
    FluidBlockBuilder &rectangle(Real length, Real height);
    /** Set the weakly-compressible fluid material parameters. */
    FluidBlockBuilder &material(Real rho0, Real c);

    const std::string &getName() const { return name_; }
    Real getLength() const { return length_; }
    Real getHeight() const { return height_; }
    Real getRho0() const { return rho0_; }
    Real getC() const { return c_; }

  private:
    std::string name_;
    Real length_{0.0};
    Real height_{0.0};
    Real rho0_{1.0};
    Real c_{10.0};
};

/**
 * @class WallBuilder
 * @brief Builder for configuring a solid wall body in a 2D simulation.
 *        Supports a fluent interface: addWall("name").hollowBox(DL, DH, BW)
 */
class WallBuilder
{
  public:
    explicit WallBuilder(const std::string &name);

    /** Define the wall as a hollow rectangular box (outer rectangle minus inner rectangle). */
    WallBuilder &hollowBox(Real domain_length, Real domain_height, Real wall_width);

    const std::string &getName() const { return name_; }
    Real getDomainLength() const { return DL_; }
    Real getDomainHeight() const { return DH_; }
    Real getWallWidth() const { return BW_; }

  private:
    std::string name_;
    Real DL_{0.0};
    Real DH_{0.0};
    Real BW_{0.0};
};

/**
 * @class SolverConfig
 * @brief Fluent configuration object for the SPH solver algorithm choices.
 *        Supports: useSolver().dualTimeStepping().freeSurfaceCorrection()
 */
class SolverConfig
{
  public:
    SolverConfig() = default;

    /** Enable dual time stepping (advection + acoustic sub-stepping). */
    SolverConfig &dualTimeStepping();
    /** Enable density summation with free-surface correction. */
    SolverConfig &freeSurfaceCorrection();

    bool isDualTimeStepping() const { return dual_time_stepping_; }
    bool isFreeSurfaceCorrection() const { return free_surface_correction_; }

  private:
    bool dual_time_stepping_{false};
    bool free_surface_correction_{false};
};

/**
 * @class SPHSimulation
 * @brief High-level facade for a 2D SPH simulation using the CK execution backend.
 *
 * Typical usage:
 * @code
 *   SPHSimulation sim;
 *   sim.createDomain(DL, DH, dp_ref);
 *   sim.addFluidBlock("Water").rectangle(LL, LH).material(rho0_f, c_f);
 *   sim.addWall("Tank").hollowBox(DL, DH, BW);
 *   sim.enableGravity(0.0, -gravity_g);
 *   sim.addObserver("PressureProbe", Vec2d(DL, 0.2));
 *   sim.useSolver().dualTimeStepping().freeSurfaceCorrection();
 *   sim.run(20.0);
 * @endcode
 */
class SPHSimulation
{
  public:
    SPHSimulation() = default;

    /** Set the 2D domain dimensions and reference particle spacing. */
    void createDomain(Real domain_length, Real domain_height, Real particle_spacing);

    /** Add a named fluid block; configure it with the returned builder. */
    FluidBlockBuilder &addFluidBlock(const std::string &name);

    /** Add a named solid wall; configure it with the returned builder. */
    WallBuilder &addWall(const std::string &name);

    /** Enable uniform gravitational acceleration. */
    void enableGravity(Real gx, Real gy);

    /** Add a point observer at the given 2D position. */
    void addObserver(const std::string &name, const Vec2d &position);

    /** Return the solver configuration object for fluent setup. */
    SolverConfig &useSolver();

    /** Build all SPH objects and run the simulation until end_time. */
    void run(Real end_time);

  private:
    Real DL_{0.0};
    Real DH_{0.0};
    Real dp_ref_{0.0};
    Vec2d gravity_{Vec2d::Zero()};
    bool gravity_enabled_{false};

    std::vector<std::unique_ptr<FluidBlockBuilder>> fluid_blocks_;
    std::vector<std::unique_ptr<WallBuilder>> walls_;

    struct ObserverEntry
    {
        std::string name;
        Vec2d position;
    };
    std::vector<ObserverEntry> observers_;

    std::unique_ptr<SolverConfig> solver_config_;
};

} // namespace SPH
#endif // SPH_SIMULATION_H
