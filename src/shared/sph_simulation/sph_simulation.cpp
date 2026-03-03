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
 * @file    sph_simulation.cpp
 * @brief   Implementation of the high-level SPHSimulation facade for 2D and 3D simulations.
 * @author  Xiangyu Hu
 */

#include "sph_simulation.h"

#include "sphinxsys.h"

#include <iomanip>
#include <iostream>

namespace SPH
{
//=================================================================================================//
FluidBlockBuilder::FluidBlockBuilder(const std::string &name)
    : name_(name) {}
//=================================================================================================//
FluidBlockBuilder &FluidBlockBuilder::block(Vecd dimensions)
{
    dimensions_ = dimensions;
    return *this;
}
//=================================================================================================//
FluidBlockBuilder &FluidBlockBuilder::material(Real rho0, Real c)
{
    rho0_ = rho0;
    c_ = c;
    return *this;
}
//=================================================================================================//
WallBuilder::WallBuilder(const std::string &name)
    : name_(name) {}
//=================================================================================================//
WallBuilder &WallBuilder::hollowBox(Vecd domain_dimensions, Real wall_width)
{
    domain_dims_ = domain_dimensions;
    BW_ = wall_width;
    return *this;
}
//=================================================================================================//
SolverConfig &SolverConfig::dualTimeStepping()
{
    dual_time_stepping_ = true;
    return *this;
}
//=================================================================================================//
SolverConfig &SolverConfig::freeSurfaceCorrection()
{
    free_surface_correction_ = true;
    return *this;
}
//=================================================================================================//
void SPHSimulation::createDomain(Vecd domain_dimensions, Real particle_spacing)
{
    domain_dims_ = domain_dimensions;
    dp_ref_ = particle_spacing;
}
//=================================================================================================//
FluidBlockBuilder &SPHSimulation::addFluidBlock(const std::string &name)
{
    fluid_blocks_.push_back(std::make_unique<FluidBlockBuilder>(name));
    return *fluid_blocks_.back();
}
//=================================================================================================//
WallBuilder &SPHSimulation::addWall(const std::string &name)
{
    walls_.push_back(std::make_unique<WallBuilder>(name));
    return *walls_.back();
}
//=================================================================================================//
void SPHSimulation::enableGravity(Vecd gravity)
{
    gravity_ = gravity;
    gravity_enabled_ = true;
}
//=================================================================================================//
void SPHSimulation::addObserver(const std::string &name, const Vecd &position)
{
    observers_.push_back({name, {position}});
}
//=================================================================================================//
void SPHSimulation::addObserver(const std::string &name, const StdVec<Vecd> &positions)
{
    observers_.push_back({name, positions});
}
//=================================================================================================//
SolverConfig &SPHSimulation::useSolver()
{
    if (!solver_config_)
        solver_config_ = std::make_unique<SolverConfig>();
    return *solver_config_;
}
//=================================================================================================//
void SPHSimulation::run(Real end_time)
{
    //----------------------------------------------------------------------
    // Validate configuration
    //----------------------------------------------------------------------
    if (fluid_blocks_.empty())
    {
        std::cerr << "SPHSimulation::run: no fluid block defined.\n";
        return;
    }
    if (walls_.empty())
    {
        std::cerr << "SPHSimulation::run: no wall defined.\n";
        return;
    }

    //----------------------------------------------------------------------
    // Derive geometry parameters
    //----------------------------------------------------------------------
    const FluidBlockBuilder &fluid_cfg = *fluid_blocks_[0];
    const WallBuilder &wall_cfg = *walls_[0];
    const Real BW = wall_cfg.getWallWidth();

    //----------------------------------------------------------------------
    // Build the SPH system
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(-BW * Vecd::Ones(), domain_dims_ + BW * Vecd::Ones());
    SPHSystem sph_system(system_domain_bounds, dp_ref_);

    //----------------------------------------------------------------------
    // Create fluid body (rectangular block starting at the coordinate origin)
    //----------------------------------------------------------------------
    Vecd water_halfsize = 0.5 * fluid_cfg.getDimensions();
    GeometricShapeBox initial_water_block(
        Transform(water_halfsize), water_halfsize, fluid_cfg.getName());
    FluidBody water_block(sph_system, initial_water_block);
    water_block.defineMaterial<WeaklyCompressibleFluid>(fluid_cfg.getRho0(), fluid_cfg.getC());
    water_block.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    // Create wall body (hollow box aligned with the domain origin)
    //----------------------------------------------------------------------
    Vecd inner_halfsize = 0.5 * wall_cfg.getDomainDimensions();
    Vecd outer_halfsize = inner_halfsize + BW * Vecd::Ones();
    // Both inner and outer box are centered at inner_halfsize
    ComplexShape wall_complex_shape(wall_cfg.getName());
    wall_complex_shape.add<GeometricShapeBox>(Transform(inner_halfsize), outer_halfsize);
    wall_complex_shape.subtract<GeometricShapeBox>(Transform(inner_halfsize), inner_halfsize);
    SolidBody wall_boundary(sph_system, wall_complex_shape);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    //----------------------------------------------------------------------
    // Create observer bodies and contacts (kept alive for the entire run)
    //----------------------------------------------------------------------
    std::vector<std::unique_ptr<Contact<>>> observer_contacts;
    for (const auto &obs : observers_)
    {
        ObserverBody &obs_body = sph_system.addBody<ObserverBody>(obs.name);
        obs_body.generateParticles<ObserverParticles>(obs.positions);
        observer_contacts.push_back(std::make_unique<Contact<>>(obs_body, StdVec<RealBody *>{&water_block}));
    }

    //----------------------------------------------------------------------
    // Define body relations
    //----------------------------------------------------------------------
    Inner<> water_block_inner(water_block);
    Contact<> water_wall_contact(water_block, {&wall_boundary});

    //----------------------------------------------------------------------
    // Build solver and particle method containers
    //----------------------------------------------------------------------
    SPHSolver sph_solver(sph_system);
    auto &main_methods = sph_solver.addParticleMethodContainer(par_ck);
    auto &host_methods = sph_solver.addParticleMethodContainer(par_host);

    //----------------------------------------------------------------------
    // Cell linked list and relation dynamics
    //----------------------------------------------------------------------
    auto &water_cell_linked_list = main_methods.addCellLinkedListDynamics(water_block);
    auto &wall_cell_linked_list = main_methods.addCellLinkedListDynamics(wall_boundary);
    auto &water_block_update_complex_relation =
        main_methods.addRelationDynamics(water_block_inner, water_wall_contact);

    std::vector<BaseDynamics<void> *> observer_relation_dynamics;
    for (auto &obs_contact : observer_contacts)
    {
        observer_relation_dynamics.push_back(
            &main_methods.addRelationDynamics(*obs_contact));
    }

    auto &particle_sort = main_methods.addSortDynamics(water_block);

    //----------------------------------------------------------------------
    // Physical dynamics
    //----------------------------------------------------------------------
    auto &wall_boundary_normal_direction =
        host_methods.addStateDynamics<NormalFromBodyShapeCK>(wall_boundary);
    auto &water_advection_step_setup =
        main_methods.addStateDynamics<fluid_dynamics::AdvectionStepSetup>(water_block);
    auto &water_update_particle_position =
        main_methods.addStateDynamics<fluid_dynamics::UpdateParticlePosition>(water_block);

    Gravity gravity_force(gravity_enabled_ ? gravity_ : Vecd::Zero());
    auto &constant_gravity =
        main_methods.addStateDynamics<GravityForceCK<Gravity>>(water_block, gravity_force);

    auto &fluid_linear_correction_matrix =
        main_methods
            .addInteractionDynamics<LinearCorrectionMatrix, WithUpdate>(water_block_inner, 0.5)
            .addPostContactInteraction(water_wall_contact);

    auto &fluid_acoustic_step_1st_half =
        main_methods
            .addInteractionDynamics<
                fluid_dynamics::AcousticStep1stHalf, OneLevel,
                AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(
                water_wall_contact);

    auto &fluid_acoustic_step_2nd_half =
        main_methods
            .addInteractionDynamics<
                fluid_dynamics::AcousticStep2ndHalf, OneLevel,
                AcousticRiemannSolverCK, LinearCorrectionCK>(water_block_inner)
            .addPostContactInteraction<Wall, AcousticRiemannSolverCK, LinearCorrectionCK>(
                water_wall_contact);

    auto &fluid_density_regularization =
        main_methods
            .addInteractionDynamics<fluid_dynamics::DensitySummationCK>(water_block_inner)
            .addPostContactInteraction(water_wall_contact)
            .addPostStateDynamics<fluid_dynamics::DensityRegularization, FreeSurface>(
                water_block);

    //----------------------------------------------------------------------
    // Time step estimators
    //----------------------------------------------------------------------
    const Real U_ref = fluid_cfg.getC() / 10.0; // c_f = 10 * U_ref => U_ref = c_f / 10
    auto &fluid_advection_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AdvectionTimeStepCK>(water_block, U_ref);
    auto &fluid_acoustic_time_step =
        main_methods.addReduceDynamics<fluid_dynamics::AcousticTimeStepCK<>>(water_block);

    //----------------------------------------------------------------------
    // I/O
    //----------------------------------------------------------------------
    auto &body_state_recorder =
        main_methods.addBodyStateRecorder<BodyStatesRecordingToVtpCK>(sph_system);
    body_state_recorder.addToWrite<Vecd>(wall_boundary, "NormalDirection");
    body_state_recorder.addToWrite<Real>(water_block, "Density");

    //----------------------------------------------------------------------
    // Observer output (pressure at observation points)
    //----------------------------------------------------------------------
    std::vector<BaseIO *> observer_pressure_outputs;
    for (auto &obs_contact : observer_contacts)
    {
        auto &recorder = main_methods.addObserveRecorder<Real>("Pressure", *obs_contact);
        observer_pressure_outputs.push_back(&recorder);
    }

    //----------------------------------------------------------------------
    // Define time stepper
    //----------------------------------------------------------------------
    TimeStepper &time_stepper = sph_solver.defineTimeStepper(end_time);

    //----------------------------------------------------------------------
    // Setup advection-step trigger
    //----------------------------------------------------------------------
    auto &advection_step = time_stepper.addTriggerByInterval(fluid_advection_time_step.exec());
    size_t advection_steps = 1;
    const int screening_interval = 100;
    const int observation_interval = screening_interval * 2;
    auto &state_recording_trigger = time_stepper.addTriggerByInterval(0.1);

    //----------------------------------------------------------------------
    // Initialise (must run host dynamics first, then device)
    //----------------------------------------------------------------------
    wall_boundary_normal_direction.exec();
    constant_gravity.exec();

    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_block_update_complex_relation.exec();
    for (auto *rel : observer_relation_dynamics)
        rel->exec();

    fluid_density_regularization.exec();
    water_advection_step_setup.exec();
    fluid_linear_correction_matrix.exec();

    //----------------------------------------------------------------------
    // First output
    //----------------------------------------------------------------------
    body_state_recorder.writeToFile();
    for (auto *obs_out : observer_pressure_outputs)
        obs_out->writeToFile(advection_steps);

    //----------------------------------------------------------------------
    // Timing
    //----------------------------------------------------------------------
    TimeInterval interval_output;
    TimeInterval interval_advection_step;
    TimeInterval interval_acoustic_step;
    TimeInterval interval_updating_configuration;

    //----------------------------------------------------------------------
    // Time integration loop
    //----------------------------------------------------------------------
    TickCount t0 = TickCount::now();
    while (!time_stepper.isEndTime())
    {
        // Fast acoustic sub-stepping
        TickCount time_instance = TickCount::now();
        Real acoustic_dt = time_stepper.incrementPhysicalTime(fluid_acoustic_time_step);
        fluid_acoustic_step_1st_half.exec(acoustic_dt);
        fluid_acoustic_step_2nd_half.exec(acoustic_dt);
        interval_acoustic_step += TickCount::now() - time_instance;

        // Slower advection stepping
        if (advection_step(fluid_advection_time_step))
        {
            advection_steps++;
            water_update_particle_position.exec();

            time_instance = TickCount::now();
            if (advection_steps % screening_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9)
                          << "N=" << advection_steps
                          << "  Time = " << time_stepper.getPhysicalTime()
                          << "  advection_dt = " << advection_step.getInterval()
                          << "  acoustic_dt = " << time_stepper.getGlobalTimeStepSize()
                          << "\n";
            }

            if (advection_steps % observation_interval == 0)
            {
                for (auto *rel : observer_relation_dynamics)
                    rel->exec();
                for (auto *obs_out : observer_pressure_outputs)
                    obs_out->writeToFile(advection_steps);
            }

            if (state_recording_trigger())
                body_state_recorder.writeToFile();

            interval_output += TickCount::now() - time_instance;

            // Update configuration
            time_instance = TickCount::now();
            if (advection_steps % 100)
                particle_sort.exec();
            water_cell_linked_list.exec();
            water_block_update_complex_relation.exec();
            interval_updating_configuration += TickCount::now() - time_instance;

            // Update dynamics for next advection step
            time_instance = TickCount::now();
            fluid_density_regularization.exec();
            water_advection_step_setup.exec();
            fluid_linear_correction_matrix.exec();
            interval_advection_step += TickCount::now() - time_instance;
        }
    }

    //----------------------------------------------------------------------
    // Summary
    //----------------------------------------------------------------------
    TimeInterval tt = TickCount::now() - t0 - interval_output;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds.\n";
    std::cout << std::fixed << std::setprecision(9)
              << "interval_advection_step = " << interval_advection_step.seconds() << "\n"
              << "interval_acoustic_step = " << interval_acoustic_step.seconds() << "\n"
              << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";
}
//=================================================================================================//
} // namespace SPH
