/*
 * @file modified_T_shaped_pipe.cpp
 * @brief This is the benchmark test of multi -inlet and multi - outlet.
 * @details We consider a flow with one inlet and two outlets in a T - shaped pipe in 2D.
 * @author Xiangyu Hu,Shuoguo Zhang
 */
#include "test_2d_aortic_valve.h" // case file to setup the test case
#include "bidirectional_buffer.h"
#include "density_correction.h"
#include "density_correction.hpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"

#include <gtest/gtest.h>

return_data run_fsi2();

int main(int argc, char *argv[])
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

TEST(fluid_shell_interaction, aortic_valve)
{
    auto data = run_fsi2();
    auto [error_x, error_y] = get_average_error(data, time_cycle, num_cycles * time_cycle);
    std::cout << "Average error in x: " << error_x * 100 << "%" << std::endl;
    std::cout << "Average error in y: " << error_y * 100 << "%" << std::endl;
    EXPECT_LT(error_x, 0.1); // 10% error tolerance
    EXPECT_LT(error_y, 0.1); // 10% error tolerance
}

//----------------------------------------------------------------------
return_data run_fsi2()
{
    std::cout << "dp_shell/thickness = " << dp_shell / shell_thickness << std::endl;
    std::cout << "U_f = " << U_f << std::endl;
    const Real Re = rho0_f * U_f * height / mu_f;
    std::cout << "The Reynolds number is " << Re << std::endl;

    // Read reference data from csv file
    auto get_linear_interpolator_from_filename = [&](const std::string &filename)
    {
        auto data = read_csv_data(filename);
        auto linear_interpolator = LinearInterpolator(data[0], data[1]);
        return linear_interpolator;
    };
    auto linear_interpolator_x = get_linear_interpolator_from_filename("input/disp_x_ref.csv");
    auto linear_interpolator_y = get_linear_interpolator_from_filename("input/disp_y_ref.csv");
    return_data data;

    // create shape
    // fluid shape
    std::cout << "Creating fluid shape ... " << std::endl;
    MultiPolygon fluid_shape_poly;
    Vec2d fluid_halfsize = 0.5 * Vec2d(length + 2 * buffer_length, height);
    Vec2d translation = 0.5 * length * Vec2d::UnitX();
    fluid_shape_poly.addBox(Transform(translation), fluid_halfsize, GeometricOps::add);
    MultiPolygonShape fluid_shape(fluid_shape_poly, "Fluid");
    std::cout << "Creating fluid shape - done. " << std::endl;

    // wall shape
    std::cout << "Creating wall shape ... " << std::endl;
    MultiPolygon wall_shape_poly;
    Vec2d wall_outer_halfsize = fluid_halfsize + Vec2d(dp_fluid, wall_thickness);
    Vec2d wall_inner_halfsize = fluid_halfsize + dp_fluid * Vec2d::UnitX();
    wall_shape_poly.addBox(Transform(translation), wall_outer_halfsize, GeometricOps::add);
    wall_shape_poly.addBox(Transform(translation), wall_inner_halfsize, GeometricOps::sub);
    // subtract beam base
    wall_shape_poly.addBox(Transform(Vec2d(shell_pos_x, 0)), Vec2d(0.5 * dp_shell, 0.5 * height + wall_thickness), GeometricOps::sub);
    MultiPolygonShape wall_shape(wall_shape_poly, "Wall");
    std::cout << "Creating wall shape - done. " << std::endl;

    // shell
    auto get_pos = [&](Real yi, Real yf)
    {
        StdVec<Vec2d> positions_shell;
        auto number_of_particles = size_t(std::floor(std::abs(yf - yi) / dp_shell));
        Real dy = (yf - yi) / Real(number_of_particles);
        for (size_t i = 0; i < number_of_particles; ++i)
        {
            Real y = yi + (Real(i) + 0.5) * dy;
            positions_shell.emplace_back(shell_pos_x, y);
        }
        return positions_shell;
    };
    //----------------------------------------------------------------------
    //	Build up SPHSystem and IO environment.
    //----------------------------------------------------------------------
    const auto bbox = wall_shape.getBounds();
    SPHSystem sph_system(bbox, dp_fluid);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody blood_body(sph_system, fluid_shape);
    blood_body.defineMatterMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    blood_body.addMaterialProperty<Viscosity>(mu_f);
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(0.5);
    blood_body.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);

    SolidBody wall_boundary(sph_system, wall_shape);
    wall_boundary.defineAdaptation<SPHAdaptation>(1.15, dp_fluid / dp_solid);
    wall_boundary.defineMatterMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    std::vector<std::unique_ptr<ShellObject>> shell_objects;
    { // shell 1
        auto positions = get_pos(0.5 * height + wall_thickness, 0.5 * height - shell_length);
        StdVec<Vec2d> normals(positions.size(), Vec2d::UnitX());
        shell_objects.emplace_back(std::make_unique<ShellObject>(sph_system, "shell1", positions, normals));
    }
    {
        // shell 2
        auto positions = get_pos(-0.5 * height - wall_thickness, -0.5 * height + shell_length);
        StdVec<Vec2d> normals(positions.size(), Vec2d::UnitX());
        shell_objects.emplace_back(std::make_unique<ShellObject>(sph_system, "shell2", positions, normals));
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //----------------------------------------------------------------------
    InnerRelation blood_inner(blood_body);
    ContactRelation blood_wall_contact(blood_body, {&wall_boundary});
    RealBodyVector shell_bodies;
    for (auto &shell_object : shell_objects)
        shell_bodies.emplace_back(&shell_object->body_);
    ContactRelationFSI2 blood_shell_contact(blood_body, shell_bodies);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // and only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation blood_body_complex(blood_inner, {&blood_wall_contact, &blood_shell_contact});
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
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);

    auto get_clamp_id = [&](SPHBody &shell_body)
    {
        IndexVector ids;
        const auto *pos = shell_body.getBaseParticles().getVariableDataByName<Vec2d>("Position");
        for (size_t i = 0; i < shell_body.getBaseParticles().TotalRealParticles(); ++i)
        {
            if (abs(pos[i].y()) > 0.5 * height)
                ids.push_back(i);
        }
        return ids;
    };
    std::vector<std::unique_ptr<BodyPartByParticle>> shell_bases;
    for (auto &shell_object : shell_objects)
    {
        shell_bases.emplace_back(std::make_unique<BodyPartByParticle>(shell_object->body_));
        shell_bases.back()->body_part_particles_ = get_clamp_id(shell_object->body_);
    }
    std::vector<std::unique_ptr<SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion>>> shell_constraints;
    for (auto &shell_base : shell_bases)
        shell_constraints.emplace_back(std::make_unique<SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion>>(*shell_base));
    auto clamping_bc = [&]()
    {
        for (const auto &shell_bodies_constraint : shell_constraints)
            shell_bodies_constraint->exec();
    };
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    InteractionDynamics<ComplexInteraction<NablaWV<Inner<>, Contact<>, Contact<>>>> kernel_summation(blood_inner, blood_wall_contact, blood_shell_contact);
    InteractionWithUpdate<ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>, Contact<>>>> boundary_indicator(blood_inner, blood_wall_contact, blood_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> pressure_relaxation(blood_inner, blood_wall_contact, blood_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver>> density_relaxation(blood_inner, blood_wall_contact, blood_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::TransportVelocityCorrection<Inner<SPHAdaptation, NoLimiter>, Contact<Boundary>, Contact<Boundary>>, NoKernelCorrection, BulkParticles>> transport_correction(blood_inner, blood_wall_contact, blood_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::ViscousForce<Inner<>, Contact<Wall>, Contact<Wall>>, fluid_dynamics::FixedViscosity, NoKernelCorrection>> viscous_force(blood_inner, blood_wall_contact, blood_shell_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(blood_body, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(blood_body);

    //-----------Inflow and periodic condition --------//
    Vec2d bidirectional_buffer_halfsize = 0.5 * Vec2d(buffer_length, 1.2 * height);
    Vec2d left_bidirectional_translation = -0.5 * buffer_length * Vec2d::UnitX();
    Vec2d right_bidirectional_translation = (length + 0.5 * buffer_length) * Vec2d::UnitX();
    OrientedBoxByCell left_emitter(blood_body, OrientedBox(xAxis, Transform(Vec2d(left_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<LeftInflowPressure> left_bidirection_buffer(left_emitter, in_outlet_particle_buffer);
    OrientedBoxByCell right_emitter(blood_body, OrientedBox(xAxis, Transform(Rotation2d(Pi), Vec2d(right_bidirectional_translation)), bidirectional_buffer_halfsize));
    fluid_dynamics::BidirectionalBuffer<RightInflowPressure> right_bidirection_buffer(right_emitter, in_outlet_particle_buffer);

    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationPressureComplex<Inner<>, Contact<>, Contact<>>> update_fluid_density(blood_inner, blood_wall_contact, blood_shell_contact);

    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_velocity_condition(left_emitter);
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<OutflowVelocity>> outflow_velocity_condition(right_emitter);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    std::vector<std::unique_ptr<ShellFluidAlgorithms<decltype(density_relaxation)>>> shell_fluid_algorithms;
    for (auto &object : shell_objects)
        shell_fluid_algorithms.emplace_back(std::make_unique<ShellFluidAlgorithms<decltype(density_relaxation)>>(object->body_, RealBodyVector{&blood_body}));
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(blood_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Vec2d>(blood_body, "Velocity");
    write_real_body_states.addToWrite<Real>(blood_body, "Density");
    write_real_body_states.addToWrite<Real>(blood_body, "Pressure");
    write_real_body_states.addToWrite<int>(blood_body, "Indicator");
    write_real_body_states.addToWrite<int>(blood_body, "BufferIndicator");
    write_real_body_states.addToWrite<Vec2d>(wall_boundary, "NormalDirection");
    for (const auto &object : shell_objects)
    {
        write_real_body_states.addToWrite<Vec2d>(object->body_, "NormalDirection");
        write_real_body_states.addDerivedVariableRecording<SimpleDynamics<Displacement>>(object->body_);
    }
    write_real_body_states.writeToFile();
    // Observer for the tip of the beam
    struct BeamObserver
    {
        ObserverBody body_;
        std::unique_ptr<ContactRelation> contact_relation_;
        std::unique_ptr<ObservedQuantityRecording<Vecd>> write_displacement_;

        BeamObserver(SPHSystem &sph_system, const std::string &name, const StdVec<Vec2d> &observation_locations, const RealBodyVector &contact_bodies)
            : body_(sph_system, name)
        {
            body_.generateParticles<ObserverParticles>(observation_locations);
            contact_relation_ = std::make_unique<ContactRelation>(body_, contact_bodies);
            write_displacement_ = std::make_unique<ObservedQuantityRecording<Vec2d>>("Displacement", *contact_relation_);
        }

        void correct_kernel()
        {
            InteractionDynamics<CorrectInterpolationKernelWeights>{*contact_relation_}.exec();
        }
    };
    std::vector<std::unique_ptr<BeamObserver>> beam_observers;
    beam_observers.emplace_back(std::make_unique<BeamObserver>(sph_system, "observer_1", StdVec<Vec2d>{Vec2d(shell_pos_x, 0.5 * height - shell_length)}, RealBodyVector{&shell_objects[0]->body_}));
    beam_observers.emplace_back(std::make_unique<BeamObserver>(sph_system, "observer_2", StdVec<Vec2d>{Vec2d(shell_pos_x, -0.5 * height + shell_length)}, RealBodyVector{&shell_objects[1]->body_}));
    auto write_beam_tip_displacement = [&]()
    {
        for (const auto &beam_observer : beam_observers)
            beam_observer->write_displacement_->writeToFile();
    };
    auto record_data = [&](Real time)
    {
        time -= time_flow_init;
        data.time_vec.emplace_back(time);
        const auto &upper_beam_disp = beam_observers[0]->write_displacement_->getObservedQuantity()[0];
        const auto &lower_beam_disp = beam_observers[1]->write_displacement_->getObservedQuantity()[0];
        Real disp_x = 0.5 * (upper_beam_disp[0] + lower_beam_disp[0]);
        Real disp_y = 0.5 * (upper_beam_disp[1] - lower_beam_disp[1]);
        data.disp_x_vec.emplace_back(disp_x / cm_to_m); // m to cm
        data.disp_y_vec.emplace_back(disp_y / cm_to_m); // m to cm
        data.disp_x_ref_vec.emplace_back(linear_interpolator_x(time));
        data.disp_y_ref_vec.emplace_back(linear_interpolator_y(time));
    };
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    // buffer and surface
    boundary_indicator.exec();
    left_bidirection_buffer.tag_buffer_particles.exec();
    right_bidirection_buffer.tag_buffer_particles.exec();
    left_bidirection_buffer.deletion.exec();
    right_bidirection_buffer.deletion.exec();
    /** computing surface normal direction for the wall. */
    wall_boundary_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    for (const auto &object : shell_objects)
        object->algs_->corrected_configuration_.exec();
    for (const auto &beam_observer : beam_observers)
        beam_observer->correct_kernel();
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

            blood_body.updateCellLinkedList();
            blood_body_complex.updateConfiguration();
            left_bidirection_buffer.tag_buffer_particles.exec();
            right_bidirection_buffer.tag_buffer_particles.exec();
        }
        std::cout << "Fluid relaxation finished !" << std::endl;
    };
    run_fluid_relaxation();
    write_real_body_states.writeToFile();
    write_beam_tip_displacement();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = time_flow_init + time_cycle * num_cycles;
    Real output_interval = time_cycle / 100.0; // output 50 frames per cycle
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
    Real dt_s_ref = shell_objects[0]->algs_->shell_computing_time_step_size_.exec();
    auto run_fsi = [&]()
    {
        while (physical_time < end_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < output_interval)
            {
                Real Dt = get_fluid_advection_time_step_size.exec();
                if (Dt < Dt_ref / 20.0)
                    throw std::runtime_error("Error: The fluid advection time step size is too small!");
                update_fluid_density.exec();
                viscous_force.exec();
                transport_correction.exec();

                if (physical_time > time_flow_init)
                {
                    /** FSI for viscous force. */
                    for (const auto &shell_fluid_algorithm : shell_fluid_algorithms)
                        shell_fluid_algorithm->viscous_force_from_fluid_.exec();
                    /** Update normal direction on elastic body.*/
                    for (const auto &object : shell_objects)
                        object->algs_->update_normal_.exec();
                }

                size_t inner_ite_dt = 0;
                size_t inner_ite_dt_s = 0;
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    Real dt = get_fluid_time_step_size.exec();
                    if (dt < dt_ref / 20.0)
                        throw std::runtime_error("Error: The fluid time step size is too small!");
                    dt = SMIN(dt, Dt);
                    /** Fluid pressure relaxation */
                    pressure_relaxation.exec(dt);

                    // Boundary conditions
                    kernel_summation.exec();
                    inflow_velocity_condition.exec();
                    outflow_velocity_condition.exec();

                    if (physical_time > time_flow_init)
                    {
                        /** FSI for pressure force. */
                        for (const auto &shell_fluid_algorithm : shell_fluid_algorithms)
                            shell_fluid_algorithm->pressure_force_from_fluid_.exec();
                    }

                    /** Fluid density relaxation */
                    density_relaxation.exec(dt);

                    /** Solid dynamics. */
                    if (physical_time > time_flow_init)
                    {
                        inner_ite_dt_s = 0;
                        Real dt_s_sum = 0.0;
                        for (const auto &shell_fluid_algorithm : shell_fluid_algorithms)
                            shell_fluid_algorithm->average_velocity_and_acceleration_.initialize_displacement_.exec();
                        while (dt_s_sum < dt)
                        {
                            Real dt_s = [&]()
                            {
                                Real dt_s_min = std::numeric_limits<Real>::max();
                                for (const auto &object : shell_objects)
                                    dt_s_min = SMIN(dt_s_min, object->algs_->shell_computing_time_step_size_.exec());
                                return dt_s_min;
                            }();
                            if (dt_s < dt_s_ref / 20.0)
                                throw std::runtime_error("Error: The shell time step size is too small!");
                            dt_s = SMIN(dt_s, dt - dt_s_sum);
                            for (const auto &object : shell_objects)
                                object->algs_->stress_relaxation_first_half_.exec(dt_s);
                            clamping_bc();
                            for (const auto &object : shell_objects)
                            {
                                object->algs_->shell_position_damping_.exec(dt_s);
                                object->algs_->shell_rotation_damping_.exec(dt_s);
                            }
                            clamping_bc();
                            for (const auto &object : shell_objects)
                                object->algs_->stress_relaxation_second_half_.exec(dt_s);
                            dt_s_sum += dt_s;
                            inner_ite_dt_s++;
                        }
                        for (const auto &shell_fluid_algorithm : shell_fluid_algorithms)
                            shell_fluid_algorithm->average_velocity_and_acceleration_.update_averages_.exec(dt);
                    }

                    relaxation_time += dt;
                    integration_time += dt;
                    physical_time += dt;
                    inner_ite_dt++;
                }

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << physical_time
                              << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt
                              << "	dt / dt_s = " << inner_ite_dt_s
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

                if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
                {
                    particle_sorting.exec();
                }
                blood_body.updateCellLinkedList();
                if (physical_time > time_flow_init)
                {
                    for (auto &object : shell_objects)
                        object->body_.updateCellLinkedList();
                }
                blood_body_complex.updateConfiguration();
                for (const auto &shell_fluid_algorithm : shell_fluid_algorithms)
                    shell_fluid_algorithm->contact_relation_.updateConfiguration();

                boundary_indicator.exec();
                left_bidirection_buffer.tag_buffer_particles.exec();
                right_bidirection_buffer.tag_buffer_particles.exec();
            }

            TickCount t2 = TickCount::now();
            /** write run-time observation into file */
            if (physical_time < time_flow_init)
            {
                for (const auto &object : shell_objects)
                    object->body_.setNewlyUpdated();
            }
            write_real_body_states.writeToFile();
            write_beam_tip_displacement();
            record_data(physical_time);
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    };

    try
    {
        run_fsi();
    }
    catch (const std::exception &e)
    {
        std::cerr << "An error occurred during the simulation: " << e.what() << std::endl;
        write_real_body_states.writeToFile();
    }
    write_data_to_csv("./return_data.csv", data);
    return data;
}
