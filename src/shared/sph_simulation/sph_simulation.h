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

#include <functional>
#include <memory>
#include <string>
#include <vector>

namespace SPH
{
class SPHSystem;

using VecdRef = Eigen::Ref<const Vecd>;

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
    FluidBlockBuilder &block(VecdRef dimensions);
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
    WallBuilder &hollowBox(VecdRef domain_dimensions, Real wall_width);

    const std::string &getName() const { return name_; }
    const Vecd &getDomainDimensions() const { return domain_dims_; }
    Real getWallWidth() const { return BW_; }

  private:
    std::string name_;
    Vecd domain_dims_{Vecd::Zero()};
    Real BW_{0.0};
};

/**
 * @enum SolidMaterialType
 * @brief Selects the elastic material model for a solid body.
 *        The choice also determines which first-half integration method is used:
 *        - SaintVenantKirchhoff → StructureIntegration1stHalfPK2 (PK2 stress)
 *        - NeoHookean           → StructureIntegration1stHalf with NoKernelCorrectionCK
 */
enum class SolidMaterialType
{
    SaintVenantKirchhoff,
    NeoHookean,
};

/**
 * @class SolidBlockBuilder
 * @brief Builder for configuring an elastic solid body in a 2D or 3D simulation.
 *
 * The body geometry is constructed by union of axis-aligned boxes.
 * An optional constrained (holder) region can further be narrowed by subtracting
 * a second box (useful for the clamped-base pattern in beam tests).
 *
 * Fluent interface example (2D oscillating beam):
 * @code
 *   auto &beam = sim.addSolidBlock("BeamBody")
 *       .addBox(base_halfsize, base_translation)
 *       .addBox(column_halfsize, column_translation)
 *       .materialSVK(rho0_s, Youngs_modulus, poisson)
 *       .constrainBox(base_halfsize, base_translation)
 *       .subtractFromConstraint(column_halfsize, column_translation)
 *       .withNumericalDamping();
 *   beam.initialVelocity(LinearProfile(beam.getReferenceSoundSpeed()));
 * @endcode
 *
 * Fluent interface example (3D twisting column):
 * @code
 *   sim.addSolidBlock("Column")
 *       .addBox(halfsize_column, translation_column)
 *       .addBox(halfsize_holder, translation_holder)
 *       .materialNeoHookean(rho0_s, Youngs_modulus, poisson)
 *       .constrainBox(halfsize_holder, translation_holder)
 *       .initialVelocity(VelocityProfile());
 * @endcode
 */
class SolidBlockBuilder
{
  public:
    explicit SolidBlockBuilder(const std::string &name);

    /** Add a box (halfsize + translation) to the solid body ComplexShape. */
    SolidBlockBuilder &addBox(VecdRef halfsize, VecdRef translation = Vecd::Zero());

    /** Set SaintVenantKirchhoff elastic material. */
    SolidBlockBuilder &materialSVK(Real rho0, Real youngs_modulus, Real poisson);

    /** Set NeoHookean elastic material. */
    SolidBlockBuilder &materialNeoHookean(Real rho0, Real youngs_modulus, Real poisson);

    /** Specify a box region that will be kinematically constrained (fixed). */
    SolidBlockBuilder &constrainBox(VecdRef halfsize, VecdRef translation = Vecd::Zero());

    /** Subtract a box from the constrained region (e.g. to create a clamped-base strip). */
    SolidBlockBuilder &subtractFromConstraint(VecdRef halfsize, VecdRef translation = Vecd::Zero());

    /** Set an initial particle velocity field via a spatial function.
     *  The callable must accept @c const Vecd & and return @c Vecd. */
    SolidBlockBuilder &initialVelocity(std::function<Vecd(const Vecd &)> vel_func);

    /** Enable numerical damping during acoustic integration. */
    SolidBlockBuilder &withNumericalDamping();

    /** Compute the reference (bulk) sound speed from material parameters.
     *  Useful for constructing velocity profiles that scale with c0. */
    Real getReferenceSoundSpeed() const;

    const std::string &getName() const { return name_; }
    const std::vector<std::pair<Vecd, Vecd>> &getBoxes() const { return boxes_; }
    Real getRho0() const { return rho0_; }
    Real getYoungsModulus() const { return youngs_modulus_; }
    Real getPoissonRatio() const { return poisson_; }
    SolidMaterialType getMaterialType() const { return mat_type_; }
    bool hasConstraint() const { return has_constraint_; }
    const Vecd &getConstraintHalfsize() const { return constraint_halfsize_; }
    const Vecd &getConstraintTranslation() const { return constraint_translation_; }
    bool hasConstraintSubtract() const { return has_constraint_subtract_; }
    const Vecd &getConstraintSubtractHalfsize() const { return constraint_subtract_halfsize_; }
    const Vecd &getConstraintSubtractTranslation() const { return constraint_subtract_translation_; }
    bool hasInitialVelocity() const { return has_initial_velocity_; }
    const std::function<Vecd(const Vecd &)> &getInitialVelocityFunc() const { return initial_velocity_func_; }
    bool hasNumericalDamping() const { return numerical_damping_; }

  private:
    std::string name_;
    std::vector<std::pair<Vecd, Vecd>> boxes_; /**< (halfsize, translation) pairs */
    Real rho0_{1.0e3};
    Real youngs_modulus_{1.0e6};
    Real poisson_{0.3};
    SolidMaterialType mat_type_{SolidMaterialType::SaintVenantKirchhoff};
    bool has_constraint_{false};
    Vecd constraint_halfsize_{Vecd::Zero()};
    Vecd constraint_translation_{Vecd::Zero()};
    bool has_constraint_subtract_{false};
    Vecd constraint_subtract_halfsize_{Vecd::Zero()};
    Vecd constraint_subtract_translation_{Vecd::Zero()};
    bool has_initial_velocity_{false};
    std::function<Vecd(const Vecd &)> initial_velocity_func_;
    bool numerical_damping_{false};
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
 * Typical 2D fluid usage:
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
 * Typical 2D solid usage:
 * @code
 *   SPHSimulation sim;
 *   sim.defineDomain(Vec2d(lower_x, lower_y), Vec2d(upper_x, upper_y), dp_ref);
 *   auto &beam = sim.addSolidBlock("BeamBody")
 *       .addBox(base_halfsize, base_translation)
 *       .addBox(column_halfsize, column_translation)
 *       .materialSVK(rho0_s, Youngs_modulus, poisson)
 *       .constrainBox(base_halfsize, base_translation)
 *       .subtractFromConstraint(column_halfsize, column_translation)
 *       .withNumericalDamping();
 *   beam.initialVelocity(MyProfile(beam.getReferenceSoundSpeed()));
 *   sim.addObserver("TipObserver", Vec2d(PL, 0.0));
 *   sim.run(1.0);
 * @endcode
 */
class SPHSimulation
{
  public:
    SPHSimulation() = default;
    ~SPHSimulation();

    /** Set the domain dimensions and reference particle spacing.
     *  Use Vec2d for 2D or Vec3d for 3D builds. */
    void defineDomain(VecdRef domain_dimensions, Real particle_spacing);

    /** Set explicit domain lower/upper bounds and reference particle spacing.
     *  Use this overload for solid-dynamics simulations whose domain is not
     *  origin-aligned (e.g. oscillating beam, twisting column).
     *  Use Vec2d for 2D or Vec3d for 3D builds. */
    void defineDomain(VecdRef lower_bound, VecdRef upper_bound, Real particle_spacing);

    /** Set the domain dimensions and reference particle spacing.
     *  Use Vec2d for 2D or Vec3d for 3D builds. */
    void createDomain(VecdRef domain_dimensions, Real particle_spacing);

    /** Add a named fluid block; configure it with the returned builder. */
    FluidBlockBuilder &addFluidBlock(const std::string &name);

    /** Add a named solid wall; configure it with the returned builder. */
    WallBuilder &addWall(const std::string &name);

    /** Add a named elastic solid body; configure it with the returned builder. */
    SolidBlockBuilder &addSolidBlock(const std::string &name);

    /** Enable uniform gravitational acceleration.
     *  Use Vec2d for 2D or Vec3d for 3D builds. */
    void enableGravity(VecdRef gravity);

    /** Add a single-point observer at the given position. */
    void addObserver(const std::string &name, VecdRef position);

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
    bool explicit_domain_bounds_{false};
    Vecd domain_lower_bound_{Vecd::Zero()};

    std::vector<std::unique_ptr<FluidBlockBuilder>> fluid_blocks_;
    std::vector<std::unique_ptr<WallBuilder>> walls_;
    std::vector<std::unique_ptr<SolidBlockBuilder>> solid_blocks_;

    struct ObserverEntry
    {
        std::string name;
        StdVec<Vecd> positions;
    };
    std::vector<ObserverEntry> observers_;

    std::unique_ptr<SolverConfig> solver_config_;
    std::unique_ptr<SPHSystem> sph_system_;

    /** Run a fluid-only simulation (internal dispatch from run()). */
    void runFluid(Real end_time);
    /** Run a solid-only simulation (internal dispatch from run()). */
    void runSolid(Real end_time);
};

} // namespace SPH
#endif // SPH_SIMULATION_H
