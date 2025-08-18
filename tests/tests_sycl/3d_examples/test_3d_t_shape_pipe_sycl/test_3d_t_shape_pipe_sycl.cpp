/**
 * @file dambreak.cpp
 * @brief 3D dambreak example using computing kernels.
 * @author Xiangyu Hu
 */
#include "sphinxsys_sycl.h"
using namespace SPH;

struct Parameters
{
    Real inlet_pressure = 1333.0; // 10mmHg
    int number_of_particles = 10;
    Real t_ref = 0.5;
    std::string fluid_file_path = "./input/full_fluid_raw.stl";
    std::string wall_file_path = "./input/wall.stl";
    // Time and output parameters
    Real end_time = 2.5;
    Real output_dt = end_time / 100.0;
    // Material parameters
    const Real rho0_f = 1060.0; /**< Density of Blood. */
    const Real mu_f = 3.5e-3;
    Real U_max = 2.5;        // From experience
    Real c_f = 10.0 * U_max; // Speed of sound
};

struct BoundaryParameter
{
    std::string name{};
    Vec3d center;
    Vec3d normal;
    Real diameter;
    bool is_outlet = false;
    Real pressure = std::numeric_limits<Real>::quiet_NaN();
    Real L_emitter = std::numeric_limits<Real>::quiet_NaN();
};

/**
 * @brief Compute the rotation matrix that aligns the X-axis to the specified
 * normal vector.
 *
 * This function returns a Rotation3d (which can be interpreted as an AngleAxis)
 * that rotates the unit X-axis to align with the provided normalized target
 * vector.
 *
 * @param normal The target normal vector.
 * @param tolerance A small tolerance value for comparing near-zero differences.
 * @return A Rotation3d representing the rotation.
 */
inline Rotation3d compute_rotation(const Vec3d &normal,
                                   Real tolerance = Real(1e-6))
{
    // 'source' is the X-axis.
    Vec3d source(Real(1), Real(0), Real(0));
    Vec3d target = normal.normalized();

    if ((source - target).norm() < tolerance)
    {
        return Rotation3d(Real(0), Vec3d::UnitX());
    }
    if ((source + target).norm() < tolerance)
    {
        Vec3d orthogonal = Vec3d(Real(0), Real(1), Real(0)).cross(source);
        if (orthogonal.norm() < tolerance)
            orthogonal = Vec3d(Real(0), Real(0), Real(1)).cross(source);
        orthogonal.normalize();
        return Rotation3d(static_cast<Real>(M_PI), orthogonal);
    }
    Eigen::Quaternion<Real> quaternion =
        Eigen::Quaternion<Real>::FromTwoVectors(source, target);
    return Rotation3d(quaternion);
}

struct BoundaryPressurePrescribed
{
    typedef WeaklyCompressibleFluid Fluid;
    Real pressure_;
    Real t_ref_;

    explicit BoundaryPressurePrescribed(Real pressure, Real t_ref)
        : pressure_(pressure), t_ref_(t_ref) {};
    Real getPressure(const Real &input_pressure, Real time) { return time < t_ref_ ? pressure_ * time / t_ref_ : pressure_; };
    Real getAxisVelocity(const Vecd &input_position, const Real &input_axis_velocity, Real time)
    {
        return input_axis_velocity;
    };
};

template <typename ExecutionPolicy>
struct PressureBC
{
    Vec3d normal;
    Vec3d center;
    Rotation3d rot;
    Vec3d buffer_halfsize;
    AlignedBox alignedbox;
    AlignedBoxByCell alignedbox_by_cell;
    fluid_dynamics::BidirectionalBoundaryCK<ExecutionPolicy, LinearCorrectionCK, BoundaryPressurePrescribed> boundary_condition;

    PressureBC(FluidBody &fluid_body, const BoundaryParameter &params, ParticleBuffer<Base> &in_outlet_particle_buffer, Real t_ref)
        : normal(params.normal.normalized()),
          center(params.center + 0.5 * params.L_emitter * normal),
          rot(compute_rotation(normal)),
          buffer_halfsize(params.L_emitter * 0.5,
                          params.diameter * 0.505,
                          params.diameter * 0.505),
          alignedbox(xAxis, Transform(rot, center), buffer_halfsize),
          alignedbox_by_cell(fluid_body, alignedbox),
          boundary_condition(alignedbox_by_cell, in_outlet_particle_buffer,
                             params.pressure, t_ref) {}
};

void run_t_shape_pipe(int ac, char *av[], Parameters &params);

int main(int ac, char *av[])
{
    Parameters params;
    run_t_shape_pipe(ac, av, params);
    return 0;
}

void run_t_shape_pipe(int ac, char *av[], Parameters &params)
{
    constexpr Real mmHg_to_Pa = 133.3224;
    // --- Define the main execution policy for this case  ---
    using MainExecutionPolicy = execution::ParallelDevicePolicy;
    // --- Section 1: Initialization of Scale and Simulation Parameters ---
    const Real scale = 0.001;
    // Process boundaries
    Real outlet_pressure = params.inlet_pressure - 0.1 * mmHg_to_Pa;
    std::vector<BoundaryParameter> boundaries;
    {
        BoundaryParameter inlet_boundary{
            .name = "inlet",
            .center = -10 * scale * Vec3d::UnitX(),
            .normal = Vec3d::UnitX(),
            .diameter = Real(1.85 * 2.0 * scale),
            .is_outlet = false,
            .pressure = params.inlet_pressure};
        boundaries.emplace_back(inlet_boundary);
    }
    {
        BoundaryParameter upper_outlet_boundary{
            .name = "upper_outlet",
            .center = scale * Vec3d(5, 30, 0),
            .normal = -Vec3d::UnitY(),
            .diameter = Real(1.45 * 2.0 * scale),
            .is_outlet = true,
            .pressure = outlet_pressure};
        boundaries.emplace_back(upper_outlet_boundary);
    }
    {
        BoundaryParameter lower_outlet_boundary{
            .name = "lower_outlet",
            .center = scale * Vec3d(5, -30, 0),
            .normal = Vec3d::UnitY(),
            .diameter = Real(1.21 * 2.0 * scale),
            .is_outlet = true,
            .pressure = outlet_pressure};
        boundaries.emplace_back(lower_outlet_boundary);
    }
    // Set additional simulation parameters
    Real min_diameter = std::min_element(boundaries.begin(), boundaries.end(),
                                         [](const BoundaryParameter &a, const BoundaryParameter &b)
                                         {
                                             return a.diameter < b.diameter;
                                         })
                            ->diameter;
    Real min_vessels_diameter = min_diameter;
    Real resolution_ref = min_diameter / params.number_of_particles;
    std::cout << "simulation_param.min_vessels_diameter = " << min_vessels_diameter << std::endl;
    std::cout << "simulation_param.resolution_ref = " << resolution_ref << std::endl;

    // --- Section 5: Set Simulation Resolution and Boundary Emitter Length ---
    Real boundary_length = 3.0 * resolution_ref;
    for (auto &BoundaryParameter : boundaries)
        BoundaryParameter.L_emitter = boundary_length;

    // --- Section 6: Build SPH System and IO Environment Setup ---
    auto water_block_shape = makeShared<ComplexShape>("WaterBlock");
    water_block_shape->add<TriangleMeshShapeSTL>(params.fluid_file_path.c_str(), Vec3d::Zero(), scale, "OuterBoundary");
    // water_block_shape->subtract<TriangleMeshShapeSTL>(params.wall_file_path.c_str(), Vec3d::Zero(), scale);
    auto wall_boundary_shape = makeShared<TriangleMeshShapeSTL>(params.wall_file_path.c_str(), Vec3d::Zero(), scale, "WallBoundary");
    auto system_bounds = wall_boundary_shape->getBounds();
    std::cout << "Domain lower bounds: " << system_bounds.first_.transpose() << std::endl;
    std::cout << "Domain upper bounds: " << system_bounds.second_.transpose() << std::endl;

    SPHSystem sph_system(system_bounds, resolution_ref);
    auto &io_env = sph_system.getIOEnvironment();
    std::filesystem::path fluid_file_path(params.fluid_file_path);
    io_env.output_folder_ += std::to_string(params.number_of_particles) + "_" +
                             fluid_file_path.parent_path().filename().string();
    if (std::filesystem::exists(io_env.output_folder_))
        std::filesystem::remove_all(io_env.output_folder_);
    std::filesystem::create_directories(io_env.output_folder_);

    // --- Section 7: Create Fluid and Solid Bodies ---
    FluidBody water_block(sph_system, water_block_shape);
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(params.rho0_f, params.c_f), params.mu_f);
    water_block.defineComponentLevelSetShape("OuterBoundary");
    ParticleBuffer<ReserveSizeFactor> in_outlet_particle_buffer(10.);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(in_outlet_particle_buffer);
    std::cout << "total water particles:" << water_block.getBaseParticles().TotalRealParticles() << std::endl;

    SolidBody wall_boundary(sph_system, wall_boundary_shape);
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    { // relaxation of wall boundary
      // Define inner relations
        InnerRelation wall_block_inner(wall_boundary);
        InnerRelation water_block_inner(water_block);
        ContactRelation water_contact(water_block, {&wall_boundary});

        // Randomize particle positions for relaxation
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_body_particles(
            wall_boundary);
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(
            water_block);

        /** A  Physics relaxation step. */
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(
            wall_block_inner);
        RelaxationStepLevelSetCorrectionComplex relaxation_step_complex(
            DynamicsArgs(water_block_inner, std::string("OuterBoundary")),
            water_contact);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_body_particles.exec(0.25);
        random_water_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_complex.SurfaceBounding().exec();
        wall_boundary.updateCellLinkedList();
        water_block.updateCellLinkedList();
        //----------------------------------------------------------------------
        //	Particle relaxation time stepping start here.
        //----------------------------------------------------------------------
        size_t ite_p = 0;
        while (ite_p < 400)
        {
            relaxation_step_inner.exec();
            relaxation_step_complex.exec();
            ite_p += 1;
            if (ite_p % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(16)
                          << "Relaxation steps for the imported model N = " << ite_p
                          << "\n";
            }
        }
        std::cout << "The physics relaxation process of wall finish !" << std::endl;
    }

    // --- Section 8: Define Body Relations and Cell Linking ---
    Inner<> water_body_inner(water_block);
    Contact<> water_wall_contact(water_block, {&wall_boundary});
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> water_cell_linked_list(water_block);
    UpdateCellLinkedList<MainExecutionPolicy, RealBody> wall_cell_linked_list(wall_boundary);
    UpdateRelation<MainExecutionPolicy, Inner<>, Contact<>> water_body_update_relation(water_body_inner, water_wall_contact);
    ParticleSortCK<MainExecutionPolicy> particle_sort(water_block);
    StateDynamics<execution::ParallelPolicy, NormalFromBodyShapeCK> wall_normal_direction(wall_boundary);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::AdvectionStepSetup> water_advection_step_setup(water_block);
    StateDynamics<MainExecutionPolicy, fluid_dynamics::UpdateParticlePosition> water_advection_step_close(water_block);
    InteractionDynamicsCK<MainExecutionPolicy, LinearCorrectionMatrixComplex>
        fluid_linear_correction_matrix(DynamicsArgs(water_body_inner, 0.5),
                                       water_wall_contact);

    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::FreeSurfaceIndicationComplexSpatialTemporalCK>
        fluid_boundary_indicator(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<
        MainExecutionPolicy,
        fluid_dynamics::AcousticStep1stHalfWithWallRiemannCorrectionCK>
        fluid_acoustic_step_1st_half(water_body_inner, water_wall_contact);

    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::AcousticStep2ndHalfWithWallNoRiemannCK>
        fluid_acoustic_step_2nd_half(water_body_inner, water_wall_contact);
    // TODO: might need to switch to the other density regularization
    InteractionDynamicsCK<
        MainExecutionPolicy,
        fluid_dynamics::DensityRegularizationComplexInternalPressureBoundary>
        fluid_density_regularization(water_body_inner, water_wall_contact);
    InteractionDynamicsCK<MainExecutionPolicy, fluid_dynamics::TransportVelocityCorrectionWallNoCorrectionBulkParticlesCK>
        transport_correction_ck(water_body_inner, water_wall_contact);

    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AdvectionTimeStepCK>
        fluid_advection_time_step(water_block, params.U_max);
    ReduceDynamicsCK<MainExecutionPolicy, fluid_dynamics::AcousticTimeStepCK<>>
        fluid_acoustic_time_step(water_block);
    InteractionDynamicsCK<MainExecutionPolicy,
                          fluid_dynamics::ViscousForceWithWallCK>
        fluid_viscous_force(water_body_inner, water_wall_contact);

    const Real Dt_ref = fluid_advection_time_step.exec() / 20.0;
    const Real dt_ref = fluid_acoustic_time_step.exec() / 20.0;

    // --- Section 9: Boundary Conditions Setup ---
    std::vector<std::unique_ptr<PressureBC<MainExecutionPolicy>>> bidirectional_pressure_conditions;
    for (const auto &boundary : boundaries)
        bidirectional_pressure_conditions.emplace_back(
            std::make_unique<PressureBC<MainExecutionPolicy>>(
                water_block, boundary, in_outlet_particle_buffer, params.t_ref));
    StateDynamics<MainExecutionPolicy, fluid_dynamics::OutflowParticleDeletion> particle_deletion(water_block);

    // --- Section 12: Setup Recording for Body States and Observers ---
    BodyStatesRecordingToVtpCK<MainExecutionPolicy> body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");
    body_states_recording.addToWrite<int>(water_block, "Indicator");
    body_states_recording.addToWrite<Real>(water_block, "VolumetricMeasure");
    body_states_recording.addToWrite<Real>(water_block, "Density");
    body_states_recording.addToWrite<int>(water_block, "BufferIndicator");
    body_states_recording.addToWrite<Vecd>(water_block, "ZeroGradientResidue");

    // --- Section 13: Simulation Initialization and Particle Updates ---
    wall_normal_direction.exec();
    water_cell_linked_list.exec();
    wall_cell_linked_list.exec();
    water_body_update_relation.exec();
    fluid_boundary_indicator.exec();
    for (auto &bc : bidirectional_pressure_conditions)
        bc->boundary_condition.tagBufferParticles();
    // Relaxation of fluid particles
    {
        size_t relaxation_fluid_itr = 0;
        while (relaxation_fluid_itr < 100)
        {
            /** Integrate time (loop) until the next output time. */
            water_advection_step_setup.exec();
            transport_correction_ck.exec();
            water_advection_step_close.exec();
            if (relaxation_fluid_itr % 10 == 0)
            {
                std::cout << std::fixed << std::setprecision(16)
                          << "N=" << relaxation_fluid_itr << "\n";
            }
            relaxation_fluid_itr++;
            /** Update cell linked list and configuration. */
            if (relaxation_fluid_itr % 25 == 0 && relaxation_fluid_itr != 1)
            {
                particle_sort.exec();
            }
            water_cell_linked_list.exec();
            water_body_update_relation.exec();
            for (auto &bc : bidirectional_pressure_conditions)
                bc->boundary_condition.tagBufferParticles();
            fluid_boundary_indicator.exec();
        }
    }
    body_states_recording.writeToFile();

    // --- Section 17: Main Simulation Loop to Find Baseline ---
    {
        SingularVariable<Real> *sv_physical_time =
            sph_system.getSystemVariableByName<Real>("PhysicalTime");
        // ─────────────────────────────────────────────────────────────────────────────
        // ─── BEGIN INSERTION: BASELINE PID (PRESSURE-DRIVEN) SECTION
        // ─────────────────
        // ─────────────────────────────────────────────────────────────────────────────
        // 4) Begin baseline time‐stepping loop.
        size_t iterations = sph_system.RestartStep();
        const int screen_output_interval = 100;
        Real dt = 0.0;
        TickCount start_time = TickCount::now();
        TimeInterval interval;

        while (sv_physical_time->getValue() < params.end_time)
        {
            Real integration_time = 0.0;

            while (integration_time < params.output_dt)
            {
                // ─── SPH PHYSICS STEPS (unchanged)
                // ─────────────────────────────────────
                fluid_density_regularization.exec();
                water_advection_step_setup.exec();
                fluid_linear_correction_matrix.exec();
                fluid_viscous_force.exec();
                transport_correction_ck.exec();

                Real Dt = fluid_advection_time_step.exec();
                if (Dt < Dt_ref)
                {
                    body_states_recording.writeToFile();
                    std::cout << "Small Dt N=" << iterations << " Time = "
                              << sph_system.getSystemVariableByName<Real>("PhysicalTime")
                                     ->getValue()
                              << " Dt = " << Dt << " dt = " << dt << "\n";
                    return;
                }

                Real current_time =
                    sph_system.getSystemVariableByName<Real>("PhysicalTime")
                        ->getValue();

                Real relaxation_time = 0.0;

                // ─── INTEGRATE ONE ACOUSTIC‐STEP
                // ──────────────────────────────────────
                while (relaxation_time < Dt)
                {
                    dt = fluid_acoustic_time_step.exec();
                    if (fluid_acoustic_time_step.exec() < dt_ref)
                    {
                        body_states_recording.writeToFile();
                        std::cout << "Small dt N=" << iterations << " Time = "
                                  << sph_system
                                         .getSystemVariableByName<Real>("PhysicalTime")
                                         ->getValue()
                                  << " Dt = " << Dt << " dt = " << dt << "\n";
                        return;
                    }
                    dt = SMIN(dt, Dt);
                    fluid_acoustic_step_1st_half.exec(dt);

                    // APPLY PRESSURE on both inlet and outlet
                    for (auto &bc : bidirectional_pressure_conditions)
                        bc->boundary_condition.applyBoundaryCondition(dt);

                    fluid_acoustic_step_2nd_half.exec(dt);
                    relaxation_time += dt;
                    integration_time += dt;
                    sph_system.getSystemVariableByName<Real>("PhysicalTime")
                        ->incrementValue(dt);
                }

                water_advection_step_close.exec();
                iterations++;

                // ─── PARTICLE INJECTION / DELETION & CELL UPDATES
                // ───────────────────
                for (auto &bc : bidirectional_pressure_conditions)
                    bc->boundary_condition.injectParticles();
                for (auto &bc : bidirectional_pressure_conditions)
                    bc->boundary_condition.indicateOutFlowParticles();
                particle_deletion.exec();
                if (iterations % 100 == 0 && iterations != 1)
                {
                    particle_sort.exec();
                }
                water_cell_linked_list.exec();
                water_body_update_relation.exec();
                for (auto &bc : bidirectional_pressure_conditions)
                    bc->boundary_condition.tagBufferParticles();
                fluid_boundary_indicator.exec();

                // ─── SCREEN‐OUTPUT LOGGING
                // ──────────────────────────────────────────
                if (iterations % screen_output_interval == 0)
                {
                    std::cout << "N=" << iterations << " Time=" << current_time << " Dt=" << Dt
                              << " dt=" << dt << "\n";
                }
            } // end inner-time‐step loop

            // Write observers/bodies once at the end of each outer iteration
            TickCount t2 = TickCount::now();
            body_states_recording.writeToFile();
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
        } // end outer “PhysicalTime < end_time” loop

        TickCount t4 = TickCount::now();
        TimeInterval total_time = t4 - start_time - interval;
        std::cout << "Total wall time for computation: " << total_time.seconds()
                  << " seconds.\n";
    }
}