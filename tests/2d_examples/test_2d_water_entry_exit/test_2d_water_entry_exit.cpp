/**
 * @file	water entry and exit.cpp
 * @brief	2D water entry and exit example with surface wetting considered.
 * @details	This is the one of FSI test cases, also one case for
 * 			understanding spatial temporal identification approach,
 *          especially when coupled with the wetting.
 * @author  Shuoguo Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real cylinder_radius = 0.055;                             /**< Cylinder radius. */
Real DL = 8.0 * cylinder_radius;                          /**< Water tank length. */
Real DH = 7.0 * cylinder_radius;                          /**< Water tank height. */
Real LL = DL;                                             /**< Water column length. */
Real LH = 3.0 * cylinder_radius;                          /**< Water column height. */
Real particle_spacing_ref = 2.0 * cylinder_radius / 40.0; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;                       /**< Thickness of tank wall. */
Vec2d cylinder_center(0.5 * DL, LH + 0.15);               /**< Location of the cylinder center. */
Vecd tethering_point(0.5 * DL, DH);                       /**< The tethering point. */
StdVec<Vecd> observer_location = {cylinder_center};       /**< Displacement observation point. */
StdVec<Vecd> wetting_observer_location =
    {cylinder_center - Vecd(0.0, cylinder_radius)}; /**< wetting observation point. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                       /**< Fluid density. */
Real rho0_s = 0.5;                       /**< Cylinder density. */
Real gravity_g = 9.81;                   /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;                 /**< Reference sound speed. */
Real mu_f = 8.9e-7;                      /**< Water dynamics viscosity. */
//----------------------------------------------------------------------
//	Wetting parameters
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 100.0 * pow(particle_spacing_ref, 2); /**< Wetting coefficient. */
Real fluid_moisture = 1.0;                                   /**< fluid moisture. */
Real cylinder_moisture = 0.0;                                /**< cylinder moisture. */
Real wall_moisture = 1.0;                                    /**< wall moisture. */
//----------------------------------------------------------------------
//	Definition for water body
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block;
    water_block.push_back(Vecd(0.0, 0.0));
    water_block.push_back(Vecd(0.0, LH));
    water_block.push_back(Vecd(LL, LH));
    water_block.push_back(Vecd(LL, 0.0));
    water_block.push_back(Vecd(0.0, 0.0));

    return water_block;
}
class WettingFluidBody : public MultiPolygonShape
{
  public:
    explicit WettingFluidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};

class WettingFluidBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingFluidBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = fluid_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Definition for wall body
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall;
    outer_wall.push_back(Vecd(-BW, -BW));
    outer_wall.push_back(Vecd(-BW, DH + BW));
    outer_wall.push_back(Vecd(DL + BW, DH + BW));
    outer_wall.push_back(Vecd(DL + BW, -BW));
    outer_wall.push_back(Vecd(-BW, -BW));

    return outer_wall;
}
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall;
    inner_wall.push_back(Vecd(0.0, 0.0));
    inner_wall.push_back(Vecd(0.0, DH));
    inner_wall.push_back(Vecd(DL, DH));
    inner_wall.push_back(Vecd(DL, 0.0));
    inner_wall.push_back(Vecd(0.0, 0.0));

    return inner_wall;
}
class WettingWallBody : public MultiPolygonShape
{
  public:
    explicit WettingWallBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
class WettingWallBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingWallBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = wall_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Definition for cylinder body
//----------------------------------------------------------------------
class WettingCylinderBody : public MultiPolygonShape
{
  public:
    explicit WettingCylinderBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addACircle(cylinder_center, cylinder_radius, 100, ShapeBooleanOps::add);
    }
};
class WettingCylinderBodyInitialCondition : public LocalDynamics
{
  public:
    explicit WettingCylinderBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = cylinder_moisture;
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};

//----------------------------------------------------------------------
//	The diffusion model of wetting
//----------------------------------------------------------------------
using CylinderFluidDiffusionDirichlet =
    DiffusionRelaxationRK2<DiffusionRelaxation<Dirichlet<KernelGradientContact>, IsotropicDiffusion>>;
//------------------------------------------------------------------------------
// Constrained part for Simbody
//------------------------------------------------------------------------------
MultiPolygon createSimbodyConstrainShape(SPHBody &sph_body)
{
    MultiPolygon multi_polygon;
    multi_polygon.addACircle(cylinder_center, cylinder_radius, 100, ShapeBooleanOps::add);
    return multi_polygon;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
    SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WettingFluidBody>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WettingWallBody>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody cylinder(sph_system, makeShared<WettingCylinderBody>("Cylinder"));
    cylinder.defineAdaptationRatios(1.15, 1.0);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineClosure<Solid, IsotropicDiffusion>(
        rho0_s, ConstructArgs(diffusion_species_name, diffusion_coeff));
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();

    ObserverBody cylinder_observer(sph_system, "CylinderObserver");
    cylinder_observer.generateParticles<ObserverParticles>(observer_location);

    ObserverBody wetting_observer(sph_system, "WettingObserver");
    wetting_observer.generateParticles<ObserverParticles>(wetting_observer_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation cylinder_inner(cylinder);
    ContactRelation water_block_contact(water_block, {&wall_boundary, &cylinder});
    ContactRelation cylinder_contact(cylinder, {&water_block});
    ContactRelation cylinder_observer_contact(cylinder_observer, {&cylinder});
    ContactRelation wetting_observer_contact(wetting_observer, {&cylinder});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        /** body topology only for particle relaxation */
        InnerRelation cylinder_inner(cylinder);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        /** Write the body state to Vtp file. */
        BodyStatesRecordingToVtp write_inserted_body_to_vtp(cylinder);
        /** Write the particle reload files. */
        ReloadParticleIO write_particle_reload_files(cylinder);
        /** A  Physics relaxation step. */
        RelaxationStepInner relaxation_step_inner(cylinder_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_inserted_body_to_vtp.writeToFile(0);
        //----------------------------------------------------------------------
        //	Relax particles of the insert body.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_inserted_body_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the fluid dynamics used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    GetDiffusionTimeStepSize get_thermal_time_step(cylinder);
    CylinderFluidDiffusionDirichlet cylinder_wetting(cylinder_contact);
    SimpleDynamics<WettingFluidBodyInitialCondition> wetting_water_initial_condition(water_block);
    SimpleDynamics<WettingWallBodyInitialCondition> wetting_wall_initial_condition(wall_boundary);
    SimpleDynamics<WettingCylinderBodyInitialCondition> wetting_cylinder_initial_condition(cylinder);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    InteractionWithUpdate<WettingCoupledSpatialTemporalFreeSurfaceIndicationComplex> free_stream_surface_indicator(water_block_inner, water_block_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> fluid_density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> fluid_density_by_summation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> fluid_acoustic_time_step(water_block);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(fluid_density_relaxation)>> pressure_force_from_fluid(cylinder_contact);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Building Simbody.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the fixed spot. */
    SimTK::Body::Rigid fixed_spot_info(SimTK::MassProperties(1.0, SimTKVec3(0), SimTK::UnitInertia(1)));
    SolidBodyPartForSimbody cylinder_constraint_area(cylinder, makeShared<MultiPolygonShape>(createSimbodyConstrainShape(cylinder), "cylinder"));
    /** Mass properties of the constrained spot. */
    SimTK::Body::Rigid tethered_spot_info(*cylinder_constraint_area.body_part_mass_properties_);
    /** Mobility of the fixed spot. */
    SimTK::MobilizedBody::Weld fixed_spot(matter.Ground(), SimTK::Transform(SimTKVec3(tethering_point[0], tethering_point[1], 0.0)),
                                          fixed_spot_info, SimTK::Transform(SimTKVec3(0)));
    /** Mobility of the tethered spot.
     * Set the mass center as the origin location of the planar mobilizer
     */
    Vecd displacement0 = cylinder_constraint_area.initial_mass_center_ - tethering_point;
    SimTK::MobilizedBody::Planar tethered_spot(fixed_spot,
                                               SimTK::Transform(SimTKVec3(displacement0[0], displacement0[1], 0.0)),
                                               tethered_spot_info, SimTK::Transform(SimTKVec3(0)));
    // discrete forces acting on the bodies
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTK::Vec3(0.0, Real(-9.81), 0.0), 0.0);
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    fixed_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.02));
    tethered_spot_info.addDecoration(SimTK::Transform(), SimTK::DecorativeSphere(0.4));
    SimTK::State state = MBsystem.realizeTopology();

    /** Time stepping method for multibody system.*/
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);
    //----------------------------------------------------------------------
    //	Coupling between SimBody and SPH..
    //----------------------------------------------------------------------
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_tethered_spot(cylinder_constraint_area, MBsystem, tethered_spot, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_tethered_spot(cylinder_constraint_area, MBsystem, tethered_spot, integ);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(sph_system);
    body_states_recording.addToWrite<Real>(water_block, "Pressure");          // output for debug
    body_states_recording.addToWrite<Real>(water_block, "Density");           // output for debug
    body_states_recording.addToWrite<int>(water_block, "Indicator");          // output for debug
    body_states_recording.addToWrite<Vecd>(wall_boundary, "NormalDirection"); // output for debug
    RestartIO restart_io(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_cylinder_displacement("Position", cylinder_observer_contact);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Real>> write_cylinder_wetting("Phi", wetting_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    cylinder_normal_direction.exec();
    wetting_water_initial_condition.exec();
    wetting_wall_initial_condition.exec();
    wetting_cylinder_initial_condition.exec();
    Real dt_thermal = get_thermal_time_step.exec();
    free_stream_surface_indicator.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    int observation_sample_interval = screen_output_interval * 2;
    int restart_output_interval = screen_output_interval * 10;
    Real end_time = 1.0;
    Real output_interval = end_time / 70.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    TimeInterval interval_computing_time_step;
    TimeInterval interval_computing_fluid_pressure_relaxation;
    TimeInterval interval_updating_configuration;
    TickCount time_instance;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    body_states_recording.writeToFile();
    write_cylinder_displacement.writeToFile(number_of_iterations);
    write_cylinder_wetting.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** outer loop for dual-time criteria time-stepping. */
            time_instance = TickCount::now();
            Real Dt = fluid_advection_time_step.exec();

            fluid_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();
            interval_computing_time_step += TickCount::now() - time_instance;

            time_instance = TickCount::now();
            Real relaxation_time = 0.0;
            Real dt = 0.0;
            viscous_force_from_fluid.exec();
            while (relaxation_time < Dt)
            {
                /** inner loop for dual-time criteria time-stepping.  */
                dt = SMIN(SMIN(dt_thermal, fluid_acoustic_time_step.exec()), Dt);
                fluid_pressure_relaxation.exec(dt);
                pressure_force_from_fluid.exec();
                fluid_density_relaxation.exec(dt);
                cylinder_wetting.exec(dt);

                integ.stepBy(dt);
                SimTK::State &state_for_update = integ.updAdvancedState();
                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, tethered_spot, force_on_tethered_spot.exec());
                constraint_tethered_spot.exec();

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }
            interval_computing_fluid_pressure_relaxation += TickCount::now() - time_instance;

            /** screen output, write body reduced values and restart files  */
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "\n";

                if (number_of_iterations % observation_sample_interval == 0 && number_of_iterations != sph_system.RestartStep())
                {
                    write_cylinder_displacement.writeToFile(number_of_iterations);
                    write_cylinder_wetting.writeToFile(number_of_iterations);
                }
                if (number_of_iterations % restart_output_interval == 0)
                    restart_io.writeToFile(number_of_iterations);
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            time_instance = TickCount::now();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            cylinder.updateCellLinkedList();
            water_block_inner.updateConfiguration();
            cylinder_inner.updateConfiguration();
            cylinder_contact.updateConfiguration();
            water_block_complex.updateConfiguration();
            free_stream_surface_indicator.exec();
            interval_updating_configuration += TickCount::now() - time_instance;
        }

        body_states_recording.writeToFile();
        TickCount t2 = TickCount::now();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds()
              << " seconds." << std::endl;
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_time_step ="
              << interval_computing_time_step.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_computing_fluid_pressure_relaxation = "
              << interval_computing_fluid_pressure_relaxation.seconds() << "\n";
    std::cout << std::fixed << std::setprecision(9) << "interval_updating_configuration = "
              << interval_updating_configuration.seconds() << "\n";

    if (sph_system.GenerateRegressionData())
    {
        write_cylinder_displacement.generateDataBase(1.0e-3);
        write_cylinder_wetting.generateDataBase(1.0e-3);
    }
    else
    {
        write_cylinder_displacement.testResult();
        write_cylinder_wetting.testResult();
    }

    return 0;
};
