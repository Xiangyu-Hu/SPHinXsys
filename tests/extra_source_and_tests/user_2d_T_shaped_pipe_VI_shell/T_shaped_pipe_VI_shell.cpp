/**
 * @file 	T_shaped_pipe.cpp
 * @brief 	This is the benchmark test of multi-inlet and multi-outlet.
 * @details We consider a flow with one inlet and two outlets in a T-shaped pipe in 2D.
 * @author 	Xiangyu Hu, Shuoguo Zhang
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.0;                        /**< Reference length. */
Real DH = 3.0;                        /**< Reference and the height of main channel. */
Real DL1 = 0.7 * DL;                  /**< The length of the main channel. */
Real resolution_ref = 0.15;           /**< Initial reference particle spacing. */
Real resolution_shell = resolution_ref;
Real BW = resolution_ref * 4.0;
Real buffer_width = resolution_ref * 4.0;
Real thickness = resolution_shell * 1.0; 
Real DL_sponge = resolution_ref * 20; /**< Reference size of the emitter buffer to impose inflow condition. */
Real level_set_refinement_ratio = resolution_ref / (0.1 * thickness);
//-------------------------------------------------------
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
Real rho0_f = 1.0; /**< Reference density of fluid. */
Real U_f = 1.0;    /**< Characteristic velocity. */
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
Real c_f = 10.0 * U_f * SMAX(Real(1), DH / (Real(2.0) * (DL - DL1)));
Real Re = 100.0;                    /**< Reynolds number. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Material parameters of the shell
//----------------------------------------------------------------------
Real rho0_s = 1.2;           /** Normalized density. */
Real Youngs_modulus = 1.0e4; /** Normalized Youngs Modulus. */
Real poisson = 0.3;          /** Poisson ratio. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
/** the water block in T shape polygon. */
std::vector<Vecd> water_block_shape{
    Vecd(-DL_sponge, 0.0), Vecd(-DL_sponge, DH), Vecd(DL1, DH), Vecd(DL1, 2.0 * DH),
    Vecd(DL, 2.0 * DH), Vecd(DL, -DH), Vecd(DL1, -DH), Vecd(DL1, 0.0), Vecd(-DL_sponge, 0.0)};
/** the outer wall polygon. */
std::vector<Vecd> outer_wall_shape{
    Vecd(-DL_sponge - BW, -thickness), Vecd(-DL_sponge - BW, DH + thickness), Vecd(DL1 - thickness, DH + thickness), Vecd(DL1 - thickness, 2.0 * DH + BW),
    Vecd(DL + thickness, 2.0 * DH + BW), Vecd(DL + thickness, -DH - BW), Vecd(DL1 - thickness, -DH - BW), Vecd(DL1 - thickness, -thickness), Vecd(-DL_sponge - BW, -thickness)};
/** the inner wall polygon. */
std::vector<Vecd> inner_wall_shape{
    Vecd(-DL_sponge - BW, 0.0), Vecd(-DL_sponge - BW, DH), Vecd(DL1, DH), Vecd(DL1, 2.0 * DH + BW),
    Vecd(DL, 2.0 * DH + BW), Vecd(DL, -DH - BW), Vecd(DL1, -DH - BW), Vecd(DL1, 0.0), Vecd(-DL_sponge - BW, 0.0)};
//----------------------------------------------------------------------
//	Define case dependent body shapes.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};

class ShellShape : public MultiPolygonShape
{
  public:
    explicit ShellShape(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};

class WallBoundary;
template <>
class ParticleGenerator<SurfaceParticles, WallBoundary> : public ParticleGenerator<SurfaceParticles>
{
    Real resolution_shell_;
    Real shell_thickness_;

  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles,
                               Real resolution_shell, Real shell_thickness)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles),
          resolution_shell_(resolution_shell), shell_thickness_(shell_thickness){};
    void prepareGeometricData() override
    {
        auto particle_number_mid_surface_01 = int((DL1 + DL_sponge + BW) / resolution_shell_);
        //std::cout << " particle_number_mid_surface_01 = " << particle_number_mid_surface_01 << std::endl;
        for (int i = 0; i < particle_number_mid_surface_01 - 1; i++)
        {
            Real x = -DL_sponge - BW + (Real(i) + 0.5) * resolution_shell_;
            // upper wall
            Real y1 = DH + 0.5 * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_shell_);
            Vec2d normal_direction_1 = Vec2d(0, 1.0);
            addSurfaceProperties(normal_direction_1, shell_thickness_);
            // lower wall
            Real y2 = - 0.5 * resolution_shell_; // lower wall
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_shell_);
            Vec2d normal_direction_2 = Vec2d(0, -1.0);
            addSurfaceProperties(normal_direction_2, shell_thickness_);
        }

        addPositionAndVolumetricMeasure(Vecd(DL1 - 0.5 * resolution_shell_ , DH + 0.5 * resolution_shell_), resolution_shell_);
        addSurfaceProperties(Vec2d(-1.0, 1.0).normalized(), shell_thickness_);
        addPositionAndVolumetricMeasure(Vecd(DL1 - 0.5 * resolution_shell_ , - 0.5 * resolution_shell_), resolution_shell_);
        addSurfaceProperties(Vec2d(-1.0, -1.0).normalized(), shell_thickness_);

        auto particle_number_mid_surface_02 = int((DH + BW) / resolution_shell_);
        //std::cout << " particle_number_mid_surface_02 = " << particle_number_mid_surface_02 << std::endl;
        for (int i = 1; i < particle_number_mid_surface_02; i++)
        {
            // upper wall
            Real y1 = DH + (Real(i) + 0.5) * resolution_shell_;
            Real x = DL1 - 0.5 * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_shell_);
            Vec2d normal_direction = Vec2d(-1.0, 0);
            addSurfaceProperties(normal_direction, shell_thickness_);
            // lower wall
            Real y2 =  -(Real(i) + 0.5) * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_shell_);
            addSurfaceProperties(normal_direction, shell_thickness_);
        }

        auto particle_number_mid_surface_03 = int((1.5 * DH + BW) / resolution_shell_);
        //std::cout << " particle_number_mid_surface_03 = " << particle_number_mid_surface_03 << std::endl;
        for (int i = 0; i < particle_number_mid_surface_03; i++)
        {
            Real y1 = 0.5 * DH + (Real(i) + 0.5) * resolution_shell_;
            // upper wall
            Real x = DL + 0.5 * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y1), resolution_shell_);
            Vec2d normal_direction = Vec2d(1.0, 0);
            addSurfaceProperties(normal_direction, shell_thickness_);
            // lower wall
            Real y2 = 0.5 * DH - (Real(i) + 0.5) * resolution_shell_;
            addPositionAndVolumetricMeasure(Vecd(x, y2), resolution_shell_);
            addSurfaceProperties(normal_direction, shell_thickness_);
        }
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
        : u_ref_(U_f), t_ref_(2.0),
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
        if (base_particles_.ParticlePositions()[index_i][0] < -DL_sponge + buffer_width
            || base_particles_.ParticlePositions()[index_i][1] > 2.0 * DH - buffer_width
            || base_particles_.ParticlePositions()[index_i][1] < -DH + buffer_width)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//-----------------------------------------------------------------------------------------------------------
//	Main program starts here.
//-----------------------------------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -DH - BW), Vec2d(DL + BW, 2.0 * DH + BW));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.cd
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    ParticleBuffer<ReserveSizeFactor> inlet_particle_buffer(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, Lattice>(inlet_particle_buffer);
    
    SolidBody shell_body(sph_system, makeShared<ShellShape>("ShellBody"));
    shell_body.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    shell_body.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(sph_system);
    shell_body.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    shell_body.generateParticles<SurfaceParticles, WallBoundary>(resolution_shell, thickness);

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation shell_inner(shell_body);
    ContactRelationFromShellToFluid water_shell_contact(water_block, {&shell_body}, {false});
    ContactRelationFromFluidToShell shell_water_contact(shell_body, {&water_block}, {false});
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell_body, water_block);
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_shell_contact);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    // shell dynamics
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(shell_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell_body);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell_body);

    /** Exert constrain on shell. */
    BoundaryGeometry boundary_geometry(shell_body, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constrain_holder(boundary_geometry);

    // fluid dynamics
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_block_inner, water_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(water_block_inner, water_shell_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_shell_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_shell_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(water_block_inner, water_shell_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);

    // Left buffer
    Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH);
    Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize;
    BodyAlignedBoxByParticle emitter(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, inlet_particle_buffer);

    // Up buffer
    Vec2d inlet_flow_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
    Vec2d inlet_flow_buffer_translation = Vec2d(-DL_sponge, 0.0) + inlet_flow_buffer_halfsize;
    BodyAlignedBoxByCell inlet_flow_buffer(water_block, makeShared<AlignedBoxShape>(xAxis, Transform(Vec2d(inlet_flow_buffer_translation)), inlet_flow_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> inflow_condition(inlet_flow_buffer);
    Vec2d disposer_up_halfsize = Vec2d(0.3 * DH, 0.5 * BW);
    Vec2d disposer_up_translation = Vec2d(DL + 0.05 * DH, 2.0 * DH) - disposer_up_halfsize;
    BodyAlignedBoxByCell disposer_up(
        water_block, makeShared<AlignedBoxShape>(yAxis, Transform(Vec2d(disposer_up_translation)), disposer_up_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_up_outflow_deletion(disposer_up);

    // Down buffer
    Vec2d disposer_down_halfsize = disposer_up_halfsize;
    Vec2d disposer_down_translation = Vec2d(DL1 - 0.05 * DH, -DH) + disposer_down_halfsize;
    BodyAlignedBoxByCell disposer_down(
        water_block, makeShared<AlignedBoxShape>(yAxis, Transform(Rotation2d(Pi), Vec2d(disposer_down_translation)), disposer_down_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_down_outflow_deletion(disposer_down);
    
    
    // FSI
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_shell(shell_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_on_shell(shell_water_contact);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_body_states(sph_system);
    write_body_states.addToWrite<Real>(water_block, "Pressure"); // output for debug
    write_body_states.addToWrite<int>(water_block, "Indicator"); // output for debug
    write_body_states.addToWrite<Real>(water_block, "Density");
    write_body_states.addToWrite<Vecd>(shell_body, "NormalDirection");
    write_body_states.addToWrite<Vecd>(shell_body, "PressureForceFromFluid");
    write_body_states.addToWrite<Real>(shell_body, "Average1stPrincipleCurvature");
    write_body_states.addToWrite<Real>(shell_body, "Average2ndPrincipleCurvature");
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    shell_corrected_configuration.exec();
    shell_average_curvature.exec();
    constrain_holder.exec();
    water_block_complex.updateConfiguration();
    shell_water_contact.updateConfiguration();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 100;
    Real end_time = 100.0;
    Real output_interval = end_time / 200.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                           /**< Default acoustic time step sizes. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
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
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();

            /** FSI for viscous force. */
            viscous_force_on_shell.exec();

            /** Dynamics including pressure relaxation. */
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                pressure_relaxation.exec(dt);

                /** FSI for pressure force. */
                pressure_force_on_shell.exec();

                inflow_condition.exec();
                density_relaxation.exec(dt);

                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    dt_s = shell_time_step_size.exec();
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    shell_stress_relaxation_first.exec(dt_s);

                    constrain_holder.exec(dt_s);

                    shell_stress_relaxation_second.exec(dt_s);
                    dt_s_sum += dt_s;

                    //write_body_states.writeToFile();
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;

            /** inflow injection*/
            emitter_inflow_injection.exec();
            disposer_up_outflow_deletion.exec();
            disposer_down_outflow_deletion.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            shell_update_normal.exec();
            shell_body.updateCellLinkedList();
            shell_curvature_inner.updateConfiguration();
            shell_average_curvature.exec();
            shell_water_contact.updateConfiguration();
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
