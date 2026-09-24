/*
 * @file modified_T_shaped_pipe.cpp
 * @brief This is the benchmark test of multi -inlet and multi - outlet.
 * @details We consider a flow with one inlet and two outlets in a T - shaped pipe in 2D.
 * @author Xiangyu Hu,Shuoguo Zhang
 */
#include "test_3d_flipping_L_plate.h" // case file to setup the test case
#include "bidirectional_buffer.h"
#include "density_correction.h"
#include "density_correction.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"

result_data run_fsi2(size_t res_factor = 1, size_t res_ratio = 1);

int main(int ac, char *av[])
{
    std::vector<result_data> results;
    results.reserve(3);
    for (auto res_factor : {8})
    {
        auto data = run_fsi2(res_factor, 1);
        results.push_back(data);
    }
    for (const auto &data : results)
    {
        std::cout << "Fluid particle number: " << data.fluid_particle_number << std::endl;
        std::cout << "Shell particle number: " << data.shell_particle_number << std::endl;
        std::cout << "Computational time: " << data.computational_time << std::endl;
    }
}
//----------------------------------------------------------------------
result_data run_fsi2(size_t res_factor, size_t res_ratio)
{
    bool run_relaxation = true;
    bool reload_particles = false;

    const Real dp_f = plate_y / Real(res_factor);
    const Real dp_s = dp_f / Real(res_ratio);
    std::cout << "dp_f = " << dp_f << ", dp_s = " << dp_s << std::endl;
    std::cout << "dp_shell/thickness = " << dp_s / plate_thickness << std::endl;
    std::cout << "U_f = " << U_f << std::endl;
    const Real nu_f = mu_f / rho0_f;
    std::cout << "The kinematic viscosity is " << nu_f << std::endl;

    const Real dp_wall = dp_f;
    const Real wall_thickness = 4 * dp_f;

    // create shape
    // fluid shape
    std::cout << "Creating fluid shape ... " << std::endl;
    auto fluid_shape = makeShared<ComplexShape>("fluid");
    Vec3d fluid_halfsize = 0.5 * (domain_upper_bound - domain_lower_bound);
    Vec3d fluid_translation = 0.5 * (domain_upper_bound + domain_lower_bound);
    fluid_shape->add<GeometricShapeBox>(Transform(fluid_translation), fluid_halfsize);
    std::cout << "Creating fluid shape - done. " << std::endl;
    auto bbox_fluid = fluid_shape->getBounds();

    // wall shape
    std::cout << "Creating wall shape ... " << std::endl;
    auto wall_halfsize = Vec3d(fluid_halfsize.x() + wall_thickness, fluid_halfsize.y() + wall_thickness, 0.5 * wall_thickness);

    // create the lower wall
    auto lower_wall_shape = makeShared<ComplexShape>("lower_wall");
    auto wall_lower_translation = Vec3d(fluid_translation.x(), fluid_translation.y(), domain_lower_bound.z() - 0.5 * wall_thickness);
    lower_wall_shape->add<GeometricShapeBox>(Transform(wall_lower_translation), wall_halfsize);
    // subtract the plate base from the wall
    Vec3d plate_base_halfsize = 0.5 * Vec3d(dp_s, plate_y, wall_thickness);
    Vec3d plate_base_translation = plate_tip_pos + 0.5 * Vec3d(0, plate_y, -wall_thickness);
    lower_wall_shape->subtract<GeometricShapeBox>(Transform(plate_base_translation), plate_base_halfsize);
    auto bbox_lower_wall = lower_wall_shape->getBounds();

    // create the upper wall
    auto upper_wall_shape = makeShared<ComplexShape>("upper_wall");
    auto wall_upper_translation = Vec3d(fluid_translation.x(), fluid_translation.y(), domain_upper_bound.z() + 0.5 * wall_thickness);
    upper_wall_shape->add<GeometricShapeBox>(Transform(wall_upper_translation), wall_halfsize);
    std::cout << "Creating wall shape - done. " << std::endl;
    auto bbox_upper_wall = upper_wall_shape->getBounds();

    // shell
    StdVec<Vec3d> positions_shell;
    {
        Real z = bbox_lower_wall.lower_.z() + 0.5 * dp_s;
        while (z < plate_z)
        {
            Real y = 0.5 * dp_s;
            while (y < plate_y)
            {
                Vec3d pos = plate_tip_pos + Vec3d(0.0, y, z);
                positions_shell.emplace_back(pos);
                y += dp_s;
            }
            z += dp_s;
        }
        while (z < plate_height)
        {
            Real y = 0.5 * dp_s;
            while (y < plate_width)
            {
                Vec3d pos = plate_tip_pos + Vec3d(0.0, y, z);
                positions_shell.emplace_back(pos);
                y += dp_s;
            }
            z += dp_s;
        }
    };
    StdVec<Vec3d> normals_shell(positions_shell.size(), Vec3d::UnitX());
    //----------------------------------------------------------------------
    //	Build up SPHSystem and IO environment.
    //----------------------------------------------------------------------
    auto bbox = bbox_fluid.add(bbox_lower_wall).add(bbox_upper_wall);
    SPHSystem sph_system(bbox, dp_f);
    auto &io_environment = IO::getEnvironment();
    io_environment.resetOutputFolder("./output_" + std::to_string(res_factor) + "_" + std::to_string(res_ratio), true);
    sph_system.setRunParticleRelaxation(run_relaxation); // Tag for run particle relaxation for body-fitted distribution
    sph_system.setReloadParticles(reload_particles);     // Tag for computation with save particles distribution
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_body(sph_system, fluid_shape);
    fluid_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    fluid_body.addMaterialProperty<Viscosity>(mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    fluid_body.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_lower_boundary(sph_system, lower_wall_shape);
    wall_lower_boundary.defineAdaptation<SPHAdaptation>(1.15, dp_f / dp_wall);
    wall_lower_boundary.defineBodyLevelSetShape();
    wall_lower_boundary.defineMatterMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_lower_boundary.generateParticles<BaseParticles, Reload>(wall_lower_boundary.Name())
        : wall_lower_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_upper_boundary(sph_system, upper_wall_shape);
    wall_upper_boundary.defineAdaptation<SPHAdaptation>(1.15, dp_f / dp_wall);
    wall_upper_boundary.defineBodyLevelSetShape();
    wall_upper_boundary.defineMatterMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_upper_boundary.generateParticles<BaseParticles, Reload>(wall_upper_boundary.Name())
        : wall_upper_boundary.generateParticles<BaseParticles, Lattice>();

    // run relaxation
    InnerRelation wall_lower_boundary_inner(wall_lower_boundary);
    InnerRelation wall_upper_boundary_inner(wall_upper_boundary);
    if (sph_system.RunParticleRelaxation())
    {
        relax_solid(wall_lower_boundary, wall_lower_boundary_inner);
        relax_solid(wall_upper_boundary, wall_upper_boundary_inner);
    }

    ShellObject shell_object(sph_system, "L_plate", positions_shell, normals_shell, dp_s, plate_thickness);

    std::cout << "fluid number: " << fluid_body.getBaseParticles().TotalRealParticles() << std::endl;
    std::cout << "shell number: " << shell_object.body_.getBaseParticles().TotalRealParticles() << std::endl;
    exit(0);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation fluid_inner(fluid_body);
    ContactRelation fluid_lower_wall_contact(fluid_body, {&wall_lower_boundary});
    ContactRelation fluid_upper_wall_contact(fluid_body, {&wall_upper_boundary});
    ContactRelationFSI2 fluid_shell_contact(fluid_body, {&shell_object.body_});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation fluid_body_complex(fluid_inner, {&fluid_lower_wall_contact, &fluid_upper_wall_contact, &fluid_shell_contact});
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    // For typical fluid-structure interaction, we first define structure dynamics,
    // Then fluid dynamics and the corresponding coupling dynamics.
    // The coupling with multi-body dynamics will be introduced at last.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> wall_lower_boundary_normal_direction(wall_lower_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_upper_boundary_normal_direction(wall_upper_boundary);
    // Boundary condition for the shell structure
    auto clamp_id = [&, z_min = domain_lower_bound.z()]()
    {
        IndexVector ids;
        const auto *pos = shell_object.body_.getBaseParticles().getVariableDataByName<Vec3d>("Position");
        for (size_t i = 0; i < shell_object.body_.getBaseParticles().TotalRealParticles(); ++i)
        {
            if (pos[i].z() < z_min)
                ids.push_back(i);
        }
        return ids;
    }();
    BodyPartByParticle clamped_part(shell_object.body_);
    clamped_part.body_part_particles_ = clamp_id;
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> shell_constraint(clamped_part);
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    InteractionDynamics<ComplexInteraction<NablaWV<Inner<>, Contact<>, Contact<>, Contact<>>>> kernel_summation(fluid_inner, fluid_lower_wall_contact, fluid_upper_wall_contact, fluid_shell_contact);
    InteractionWithUpdate<ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>, Contact<>, Contact<>>>> boundary_indicator(fluid_inner, fluid_lower_wall_contact, fluid_upper_wall_contact, fluid_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> pressure_relaxation(fluid_inner, fluid_lower_wall_contact, fluid_upper_wall_contact, fluid_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<FreeSlipWall>, Contact<Wall>>, AcousticRiemannSolver>> density_relaxation(fluid_inner, fluid_lower_wall_contact, fluid_upper_wall_contact, fluid_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::TransportVelocityCorrection<Inner<SPHAdaptation, NoLimiter>, Contact<Boundary>, Contact<Boundary>, Contact<Boundary>>, NoKernelCorrection, BulkParticles>> transport_correction(fluid_inner, fluid_lower_wall_contact, fluid_upper_wall_contact, fluid_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::ViscousForce<Inner<>, Contact<Wall>, Contact<FreeSlipWall>, Contact<Wall>>, fluid_dynamics::FixedViscosity, NoKernelCorrection>> viscous_force(fluid_inner, fluid_lower_wall_contact, fluid_upper_wall_contact, fluid_shell_contact);
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(fluid_inner);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(fluid_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(fluid_body);

    //-----------Inflow and periodic condition --------//
    Real buffer_length = 3 * dp_f;
    auto bidirectional_buffer_halfsize = Vec3d(0.5 * buffer_length, 2 * fluid_halfsize.y(), 2 * fluid_halfsize.z());
    auto left_bidirectional_translation = Vec3d(bbox_fluid.lower_.x() + 0.5 * buffer_length, fluid_translation.y(), fluid_translation.z());
    auto right_bidirectional_translation = Vec3d(bbox_fluid.upper_.x() - 0.5 * buffer_length, fluid_translation.y(), fluid_translation.z());
    OrientedBoxByCell left_emitter(fluid_body, OrientedBox(xAxis, Transform(Vec3d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    auto right_quaternion = Eigen::Quaternion<Real>::FromTwoVectors(Vec3d::UnitX(), -Vec3d::UnitX());
    OrientedBoxByCell right_emitter(fluid_body, OrientedBox(xAxis, Transform(Rotation3d(right_quaternion), Vec3d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationPressureComplex<Inner<>, Contact<>, Contact<>, Contact<>>> update_fluid_density(fluid_inner, fluid_lower_wall_contact, fluid_upper_wall_contact, fluid_shell_contact);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);

    PeriodicAlongAxis periodic_along_y(fluid_body.getSPHBodyBounds(), yAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(fluid_body, periodic_along_y);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    ShellFluidAlgorithms<decltype(density_relaxation)> shell_fluid_algorithms(shell_object.body_, {&fluid_body});
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(fluid_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Vec3d>(fluid_body, "Velocity");
    write_real_body_states.addToWrite<Real>(fluid_body, "Density");
    write_real_body_states.addToWrite<Real>(fluid_body, "Pressure");
    write_real_body_states.addToWrite<int>(fluid_body, "Indicator");
    write_real_body_states.addToWrite<int>(fluid_body, "BufferIndicator");
    write_real_body_states.addToWrite<Vec3d>(shell_object.body_, "NormalDirection");
    write_real_body_states.addDerivedVariableRecording<SimpleDynamics<Displacement>>(shell_object.body_);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    ObserverBody beam_observer(sph_system, "BeamObserver");
    StdVec<Vec3d> beam_observation_location = {plate_tip_pos + Vec3d(0.0, plate_width, plate_height)};
    beam_observer.generateParticles<ObserverParticles>(beam_observation_location);
    ContactRelation beam_observer_contact(beam_observer, {&shell_object.body_});
    ObservedQuantityRecording<Vec3d> write_tip_displacement("Position", beam_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    // buffer and surface
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    left_bidirection_buffer.deletion.exec();
    right_bidirection_buffer.deletion.exec();
    /** computing surface normal direction for the wall. */
    wall_lower_boundary_normal_direction.exec();
    wall_upper_boundary_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    shell_object.algs_->corrected_configuration_.exec();
    InteractionDynamics<CorrectInterpolationKernelWeights>{beam_observer_contact}.exec();
    //----------------------------------------------------------------------
    // initial relaxation of fluid body
    //----------------------------------------------------------------------
    auto run_fluid_relaxation = [&]()
    {
        size_t relaxation_fluid_itr = 0;
        std::cout << "Fluid relaxation starting..." << std::endl;
        while (relaxation_fluid_itr < 100)
        {
            boundary_indicator.exec();

            transport_correction.exec();
            relaxation_fluid_itr++;

            // first do injection for all buffers
            left_bidirection_buffer.injection.exec();
            right_bidirection_buffer.injection.exec();
            // then do deletion for all buffers
            left_bidirection_buffer.deletion.exec();
            right_bidirection_buffer.deletion.exec();

            periodic_condition.bounding_.exec();
            fluid_body.updateCellLinkedList();
            periodic_condition.update_cell_linked_list_.exec();
            fluid_body_complex.updateConfiguration();

            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }
        std::cout << "Fluid relaxation finished !" << std::endl;
    };
    run_fluid_relaxation();
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real output_interval = end_time / 100.0; // output 50 frames per cycle
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    Real Dt_ref = get_fluid_advection_time_step_size.exec();
    Real dt_ref = get_fluid_time_step_size.exec();
    Real dt_s_ref = shell_object.algs_->shell_computing_time_step_size_.exec();
    std::cout << "The reference fluid advection time step size is " << Dt_ref << std::endl;
    std::cout << "The reference fluid time step size is " << dt_ref << std::endl;
    std::cout << "The reference shell time step size is " << dt_s_ref << std::endl;
    auto run_fsi = [&]()
    {
        while (physical_time < end_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < output_interval)
            {
                Real Dt = get_fluid_advection_time_step_size.exec();
                Real dt = 0;
                Real dt_s = 0;
                if (Dt < Dt_ref / 20.0)
                    throw std::runtime_error("Error: The fluid advection time step size is too small!");

                update_fluid_density.exec();
                viscous_force.exec();
                transport_correction.exec();

                if (physical_time > fsi_start_time)
                {
                    /** FSI for viscous force. */
                    shell_fluid_algorithms.viscous_force_from_fluid_.exec();
                    /** Update normal direction on elastic body.*/
                    shell_object.algs_->update_normal_.exec();
                }

                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    dt = get_fluid_time_step_size.exec();
                    if (dt < dt_ref / 20.0)
                        throw std::runtime_error("Error: The fluid time step size is too small!");
                    dt = SMIN(dt, Dt);
                    /** Fluid pressure relaxation */
                    pressure_relaxation.exec(dt);

                    // Boundary conditions
                    kernel_summation.exec();
                    inflow_velocity_condition.exec();

                    if (physical_time > fsi_start_time)
                    {
                        /** FSI for pressure force. */
                        shell_fluid_algorithms.pressure_force_from_fluid_.exec();
                    }

                    /** Fluid density relaxation */
                    density_relaxation.exec(dt);

                    /** Solid dynamics. */
                    if (physical_time > fsi_start_time)
                    {
                        Real dt_s_sum = 0.0;
                        shell_fluid_algorithms.average_velocity_and_acceleration_.initialize_displacement_.exec();
                        while (dt_s_sum < dt)
                        {
                            dt_s = shell_object.algs_->shell_computing_time_step_size_.exec();
                            if (dt_s < dt_s_ref / 20.0)
                                throw std::runtime_error("Error: The shell time step size is too small!");
                            dt_s = SMIN(dt_s, dt - dt_s_sum);
                            shell_object.algs_->stress_relaxation_first_half_.exec(dt_s);
                            shell_constraint.exec();
                            shell_object.algs_->stress_relaxation_second_half_.exec(dt_s);
                            dt_s_sum += dt_s;
                        }
                        shell_fluid_algorithms.average_velocity_and_acceleration_.update_averages_.exec(dt);
                    }

                    relaxation_time += dt;
                    integration_time += dt;
                    physical_time += dt;
                }

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << physical_time
                              << "	Dt = " << Dt << "	dt = " << dt
                              << "	dt_s = " << dt_s
                              << "\n";
                }
                number_of_iterations++;

                /** Water block configuration and periodic condition. */
                // first do injection for all buffers
                left_bidirection_buffer.injection.exec();
                right_bidirection_buffer.injection.exec();
                // then do deletion for all buffers
                left_bidirection_buffer.deletion.exec();
                right_bidirection_buffer.deletion.exec();

                periodic_condition.bounding_.exec();
                if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                {
                    particle_sorting.exec();
                }
                fluid_body.updateCellLinkedList();
                periodic_condition.update_cell_linked_list_.exec();
                if (physical_time > fsi_start_time)
                {
                    shell_object.body_.updateCellLinkedList();
                }
                fluid_body_complex.updateConfiguration();
                shell_fluid_algorithms.contact_relation_.updateConfiguration();

                boundary_indicator.exec();
                left_bidirection_buffer.tag_buffer_particles.exec();
                right_bidirection_buffer.tag_buffer_particles.exec();

                TickCount t2 = TickCount::now();
                write_tip_displacement.writeToFile();
                TickCount t3 = TickCount::now();
                interval += t3 - t2;
            }

            TickCount t2 = TickCount::now();
            /** write run-time observation into file */
            compute_vorticity.exec();
            shell_object.body_.setNewlyUpdated();
            write_real_body_states.writeToFile();
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
        return tt.seconds();
    };

    result_data result{.fluid_particle_number = fluid_body.getBaseParticles().TotalRealParticles(), .shell_particle_number = shell_object.body_.getBaseParticles().TotalRealParticles()};
    try
    {
        auto computational_time = run_fsi();
        result.computational_time = computational_time;
    }
    catch (const std::exception &e)
    {
        std::cerr << "An error occurred during the simulation: " << e.what() << std::endl;
        write_real_body_states.writeToFile();
    }
    return result;
}
