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
 *          Works for both 2D and 3D simulations via the dimension-agnostic Vecd type.
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
 * @brief Builder for configuring a fluid body in a 2D or 3D simulation.
 *
 * Fluent interface example (2D):
 * @code
 *   sim.addFluidBlock("Water").block(Vec2d(LL, LH)).material(rho0_f, c_f);
 * @endcode
 * Fluent interface example (3D):
 * @code
 *   sim.addFluidBlock("Water").block(Vec3d(LL, LH, LW)).material(rho0_f, c_f);
 * @endcode
 */
class FluidBlockBuilder
{
  public:
    explicit FluidBlockBuilder(const std::string &name);

    /** Define the fluid block dimensions (starting at the coordinate origin).
     *  Use Vec2d for 2D or Vec3d for 3D builds. */
    FluidBlockBuilder &block(Vecd dimensions);
    /** Set the weakly-compressible fluid material parameters. */
    FluidBlockBuilder &material(Real rho0, Real c);

    const std::string &getName() const { return name_; }
    const Vecd &getDimensions() const { return dimensions_; }
    Real getRho0() const { return rho0_; }
    Real getC() const { return c_; }

  private:
    std::string name_;
    Vecd dimensions_{Vecd::Zero()};
    Real rho0_{1.0};
    Real c_{10.0};
};

/**
 * @class WallBuilder
 * @brief Builder for configuring a solid wall body in a 2D or 3D simulation.
 *
 * Fluent interface example (2D):
 * @code
 *   sim.addWall("Tank").hollowBox(Vec2d(DL, DH), BW);
 * @endcode
 * Fluent interface example (3D):
 * @code
 *   sim.addWall("Tank").hollowBox(Vec3d(DL, DH, DW), BW);
 * @endcode
 */
class WallBuilder
{
  public:
    explicit WallBuilder(const std::string &name);

    /** Define the wall as a hollow rectangular box aligned with the origin.
     *  @param domain_dimensions Inner domain dimensions (Vecd for 2D/3D).
     *  @param wall_width Thickness of the wall. */
    WallBuilder &hollowBox(Vecd domain_dimensions, Real wall_width);

    const std::string &getName() const { return name_; }
    const Vecd &getDomainDimensions() const { return domain_dims_; }
    Real getWallWidth() const { return BW_; }

  private:
    std::string name_;
    Vecd domain_dims_{Vecd::Zero()};
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
 * @brief High-level facade for a 2D or 3D SPH simulation using the CK execution backend.
 *
 * Typical 2D usage:
 * @code
 *   SPHSimulation sim;
 *   sim.createDomain(Vec2d(DL, DH), dp_ref);
 *   sim.addFluidBlock("Water").block(Vec2d(LL, LH)).material(rho0_f, c_f);
 *   sim.addWall("Tank").hollowBox(Vec2d(DL, DH), BW);
 *   sim.enableGravity(Vec2d(0.0, -gravity_g));
 *   sim.addObserver("PressureProbe", Vec2d(DL, 0.2));
 *   sim.useSolver().dualTimeStepping().freeSurfaceCorrection();
 *   sim.run(20.0);
 * @endcode
 *
 * Typical 3D usage:
 * @code
 *   SPHSimulation sim;
 *   sim.createDomain(Vec3d(DL, DH, DW), dp_ref);
 *   sim.addFluidBlock("Water").block(Vec3d(LL, LH, LW)).material(rho0_f, c_f);
 *   sim.addWall("Tank").hollowBox(Vec3d(DL, DH, DW), BW);
 *   sim.enableGravity(Vec3d(0.0, -gravity_g, 0.0));
 *   sim.addObserver("Probe", {Vec3d(DL, 0.2, 0.5*DW), Vec3d(DL, 0.1, 0.5*DW)});
 *   sim.useSolver().dualTimeStepping().freeSurfaceCorrection();
 *   sim.run(20.0);
 * @endcode
 */
class SPHSimulation
{
  public:
    SPHSimulation() = default;

    /** Set the domain dimensions and reference particle spacing.
     *  Use Vec2d for 2D or Vec3d for 3D builds. */
    void createDomain(Vecd domain_dimensions, Real particle_spacing);

    /** Add a named fluid block; configure it with the returned builder. */
    FluidBlockBuilder &addFluidBlock(const std::string &name);

    /** Add a named solid wall; configure it with the returned builder. */
    WallBuilder &addWall(const std::string &name);

    /** Enable uniform gravitational acceleration.
     *  Use Vec2d for 2D or Vec3d for 3D builds. */
    void enableGravity(Vecd gravity);

    /** Add a single-point observer at the given position. */
    void addObserver(const std::string &name, const Vecd &position);

    /** Add a multi-point observer at the given positions. */
    void addObserver(const std::string &name, const StdVec<Vecd> &positions);

    /** Return the solver configuration object for fluent setup. */
    SolverConfig &useSolver();

    /** Build all SPH objects and run the simulation until end_time. */
    void run(Real end_time);

  private:
    Vecd domain_dims_{Vecd::Zero()};
    Real dp_ref_{0.0};
    Vecd gravity_{Vecd::Zero()};
    bool gravity_enabled_{false};

    std::vector<std::unique_ptr<FluidBlockBuilder>> fluid_blocks_;
    std::vector<std::unique_ptr<WallBuilder>> walls_;

    struct ObserverEntry
    {
        std::string name;
        StdVec<Vecd> positions;
    };
    std::vector<ObserverEntry> observers_;

    std::unique_ptr<SolverConfig> solver_config_;
};

} // namespace SPH
#endif // SPH_SIMULATION_H
