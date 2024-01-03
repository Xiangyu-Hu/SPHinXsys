/**
 * @file 	test_2d_fluid_around_balloon_shell.cpp
 * @brief 	Test on fluid-shell interaction when 2 shell particles are close to each other
 * @details This is a case to test fluid-shell interaction.
 * @author 	Weiyi Kong, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

// Make sure the initial normal direction points from shell to fluid
class CorrectNormalDirection
{
  private:
    BaseParticles *particles_;
    Vecd direction_point_;
    bool if_reverse_;

  public:
    CorrectNormalDirection(SPHBody &shell, Vecd direction_point, bool if_reverse = false)
        : particles_(&shell.getBaseParticles()), direction_point_(direction_point), if_reverse_(if_reverse){};
    inline void exec()
    {
        auto &n = *particles_->getVariableByName<Vecd>("NormalDirection");
        particle_for(
            par,
            particles_->total_real_particles_,
            [&, this](size_t index_i)
            {
                Vecd displacement = particles_->pos_[index_i] - direction_point_;
                if (if_reverse_)
                    displacement *= -1;
                if (n[index_i].dot(displacement) > 0)
                    n[index_i] *= -1;
            });
    }
};
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const bool if_fsi = true;
const Real scale = 1.0;
const Real DL = 7.5 * scale; /**< Reference length. */
const Real DH = 3.0 * scale; /**< Reference and the height of main channel. */
const Real DL_balloon = 4.0 * scale;
const Real radius_balloon = 0.5 * scale; // radius of the mid surface

const Real resolution_ref = 0.15 * scale;
const Real resolution_shell = resolution_ref * 0.5; // outer surface radius

const Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
const Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */

const Real thickness_balloon = resolution_shell;                            // thickness of the balloon
const Real radius_balloon_outer = radius_balloon + 0.5 * thickness_balloon; // radius of the outer surface
const Real radius_balloon_inner = radius_balloon - 0.5 * thickness_balloon; // radius of the outer surface
const Real level_set_refinement_ratio = resolution_shell / (0.1 * thickness_balloon);

const BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -0.5 * DH - BW), Vec2d(DL + BW, 0.5 * DH + BW));
const Vec2d balloon_center(0.5 * DL, 0);
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1.0; /**< Reference density of fluid. */
const Real U_f = 1.0;    /**< Characteristic velocity. */
const Real U_max = 1.5 * U_f * DH / (DH - 2 * radius_balloon_outer);
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
const Real c_f = 10.0 * U_max;
const Real Re = 100.0;                    /**< Reynolds number. */
const Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
const Real rho0_s = 10.0;    /**< Reference density.*/
const Real poisson = 0.4;    /**< Poisson ratio.*/
const Real Ae = 3.5 * 1.4e5; /**< Normalized Youngs Modulus. */
const Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * thickness_balloon;
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_shape;
    water_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));
    water_shape.push_back(Vecd(-DL_sponge, 0.5 * DH));
    water_shape.push_back(Vecd(DL, 0.5 * DH));
    water_shape.push_back(Vecd(DL, -0.5 * DH));
    water_shape.push_back(Vecd(-DL_sponge, -0.5 * DH));

    return water_shape;
}
/** create a wall outer shape */
std::vector<Vecd> createWallOuterShape()
{
    // geometry
    std::vector<Vecd> outer_shape;
    outer_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));
    outer_shape.push_back(Vecd(-DL_sponge - BW, 0.5 * DH + BW));
    outer_shape.push_back(Vecd(DL + BW, 0.5 * DH + BW));
    outer_shape.push_back(Vecd(DL + BW, -0.5 * DH - BW));
    outer_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH - BW));

    return outer_shape;
}
/** create a wall inner shape */
std::vector<Vecd> createWallInnerShape()
{
    // geometry
    std::vector<Vecd> inner_shape;
    inner_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));
    inner_shape.push_back(Vecd(-DL_sponge - BW, 0.5 * DH));
    inner_shape.push_back(Vecd(DL + BW, 0.5 * DH));
    inner_shape.push_back(Vecd(DL + BW, -0.5 * DH));
    inner_shape.push_back(Vecd(-DL_sponge - BW, -0.5 * DH));

    return inner_shape;
}
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
const Vec2d circle_center_1 = balloon_center - Vec2d(0.5 * DL_balloon, 0);
const Vec2d circle_center_2 = balloon_center + Vec2d(0.5 * DL_balloon, 0);
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        multi_polygon_.addABox(Transform(balloon_center), Vec2d(0.5 * DL_balloon, radius_balloon_outer), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_1, radius_balloon_outer, 100, ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_2, radius_balloon_outer, 100, ShapeBooleanOps::sub);
    }
};
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWallOuterShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createWallInnerShape(), ShapeBooleanOps::sub);
    }
};
class Shell : public MultiPolygonShape
{
  public:
    explicit Shell(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addABox(Transform(balloon_center), Vec2d(0.5 * DL_balloon, radius_balloon_outer), ShapeBooleanOps::add);
        multi_polygon_.addACircle(circle_center_1, radius_balloon_outer, 100, ShapeBooleanOps::add);
        multi_polygon_.addACircle(circle_center_2, radius_balloon_outer, 100, ShapeBooleanOps::add);
        multi_polygon_.addABox(Transform(balloon_center), Vec2d(0.5 * DL_balloon, radius_balloon_inner), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_1, radius_balloon_inner, 100, ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_2, radius_balloon_inner, 100, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(1.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Define the boundary geometry
//----------------------------------------------------------------------
class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };
    virtual ~BoundaryGeometry(){};

  private:
    void tagManually(size_t index_i)
    {
        if (std::abs(base_particles_.pos_[index_i][1]) < 2 * resolution_shell)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    SolidBody shell(sph_system, makeShared<Shell>("Shell"));
    shell.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        shell.generateParticles<ParticleGeneratorReload>(io_environment, shell.getName());
    }
    else
    {
        shell.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(io_environment);
        shell.generateParticles<ThickSurfaceParticleGeneratorLattice>(thickness_balloon);
    }

    if (!sph_system.RunParticleRelaxation() && !sph_system.ReloadParticles())
    {
        std::cout << "Error: This case requires reload shell particles for simulation!" << std::endl;
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_inner(water_block);
    InnerRelation shell_inner(shell);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelationToShell water_shell_contact(water_block, {&shell});
    ContactRelationFromShell shell_water_contact(shell, {&water_block});
    ComplexRelation water_block_complex(water_inner, {&water_wall_contact, &water_shell_contact});
    // inner relation to compute curvature
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell, water_block);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation shell_inner(shell);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for wall boundary.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> shell_random_particles(shell);
        relax_dynamics::ShellRelaxationStep relaxation_step_shell_inner(shell_inner);
        relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(shell_inner, thickness_balloon, cos(Pi / 2.0));
        shell.addBodyStateForRecording<int>("UpdatedIndicator");
        CorrectNormalDirection correct_normal_direction(shell, balloon_center, true);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
        MeshRecordingToPlt write_mesh_cell_linked_list(io_environment, shell.getCellLinkedList());
        ReloadParticleIO write_particle_reload_files(io_environment, {&shell});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        shell_random_particles.exec(0.25);

        relaxation_step_shell_inner.MidSurfaceBounding().exec();
        write_relaxed_particles.writeToFile(0);
        shell.updateCellLinkedList();
        write_mesh_cell_linked_list.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            for (int k = 0; k < 2; ++k)
                relaxation_step_shell_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_relaxed_particles.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        correct_normal_direction.exec();
        shell_normal_prediction.exec();
        correct_normal_direction.exec();
        shell.setNewlyUpdated();
        write_relaxed_particles.writeToFile(ite);
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Algorithm for fluid dynamics. */
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization(water_block);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeStream>, Contact<>, Contact<>>> update_fluid_density_by_summation(water_inner, water_wall_contact, water_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> fluid_pressure_relaxation(water_inner, water_wall_contact, water_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, NoRiemannSolver>> fluid_density_relaxation(water_inner, water_wall_contact, water_shell_contact);
    InteractionDynamics<ComplexInteraction<fluid_dynamics::ViscousAcceleration<Inner<>, Contact<Wall>, Contact<Wall>>>> viscous_acceleration(water_inner, water_wall_contact, water_shell_contact);
    InteractionWithUpdate<ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>, Contact<>>>> inlet_outlet_surface_particle_indicator(water_inner, water_wall_contact, water_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::TransportVelocityCorrection<Inner<SingleResolution>, Contact<Boundary>, Contact<Boundary>>, NoKernelCorrection, BulkParticles>> transport_velocity_correction(water_inner, water_wall_contact, water_shell_contact);
    /** Algorithm for in-/outlet. */
    Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
    Vec2d emitter_translation = Vec2d(-DL_sponge + 0.5 * BW, 0.0);
    BodyAlignedBoxByParticle emitter(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);

    Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
    Vec2d inlet_buffer_translation = Vec2d(-0.5 * DL_sponge, 0.0);
    BodyAlignedBoxByCell emitter_buffer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(emitter_buffer);

    Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.6 * DH);
    Vec2d disposer_translation = Vec2d(DL, 0.0);
    BodyAlignedBoxByCell disposer(
        water_block, makeShared<AlignedBoxShape>(Transform(Vec2d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, xAxis);
    /** Algorithm for solid dynamics. */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(shell_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_curvature(shell_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell);
    /** FSI */
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_shell(shell_water_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_shell_update(shell_water_contact, viscous_force_on_shell);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell);
    /** constraint and damping */
    BoundaryGeometry shell_boundary_geometry(shell, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> shell_constrain(shell_boundary_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_position_damping(0.2, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_rotation_damping(0.2, shell_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    shell.addBodyStateForRecording<Real>("TotalMeanCurvature");
    BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    shell_corrected_configuration.exec();
    shell_curvature.exec();
    water_block_complex.updateConfiguration();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 10;
    Real end_time = 5.0;
    Real output_interval = end_time / 100.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                           /**< Default acoustic time step sizes. */
    Real dt_s = 0.0;                         /**< Default acoustic time step sizes for solid. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_body_states.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            fluid_step_initialization.exec();
            Real Dt = fluid_advection_time_step.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_fluid_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_velocity_correction.exec();

            if (if_fsi)
                viscous_force_on_shell.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(fluid_acoustic_time_step.exec(), Dt - relaxation_time);
                fluid_pressure_relaxation.exec(dt);
                emitter_buffer_inflow_condition.exec();
                if (if_fsi)
                    fluid_force_on_shell_update.exec();
                fluid_density_relaxation.exec(dt);

                /** Solid dynamics time stepping. */
                if (if_fsi)
                {
                    Real dt_s_sum = 0.0;
                    average_velocity_and_acceleration.initialize_displacement_.exec();
                    while (dt_s_sum < dt)
                    {
                        dt_s = std::min(shell_time_step_size.exec(), dt - dt_s_sum);
                        shell_stress_relaxation_first.exec(dt_s);
                        shell_constrain.exec();
                        shell_position_damping.exec(dt_s);
                        shell_rotation_damping.exec(dt_s);
                        shell_constrain.exec();
                        shell_stress_relaxation_second.exec(dt_s);
                        dt_s_sum += dt_s;
                    }
                    average_velocity_and_acceleration.update_averages_.exec(dt);
                }

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt;
                if (if_fsi)
                    std::cout << "  dt_s = " << dt_s;
                std::cout << "\n";
            }
            number_of_iterations++;

            if (if_fsi)
                shell_update_normal.exec();

            /** inflow injection*/
            emitter_inflow_injection.exec();
            disposer_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            if (if_fsi)
            {
                shell.updateCellLinkedList();
                shell_water_contact.updateConfiguration();
                shell_curvature_inner.updateConfiguration();
                shell_curvature.exec();
            }
            water_block_complex.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        write_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;

    return 0;
}
