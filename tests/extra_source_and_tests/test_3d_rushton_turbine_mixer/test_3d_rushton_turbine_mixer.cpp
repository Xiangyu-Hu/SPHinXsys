/**
 * @file	HVMixer.cpp
 * @brief	Non-Newtonian Rushton Turbine Mixer
 * @details	This case is based on a Rushton turbine impeller in a cylindrical tank fitted with vertical baffles using the Herschel Bulkley model
 * @author	Theodor Hennings
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// setup parameters
bool relaxation = false;
bool linearized_iteration = false;
bool reload = false;

// simulation setup
Real particle_spacing = 0.002;                                  // particle spacing (must be small enough so blades are resolved)
Real gravity_g = 9.81;                                          // gravity
Real end_time = 0.5;                                            // end time
Real RPS = 5;                                                   // revolutions per second
Real omega = RPS * 3.14 * 2;                                    // angular velocity
Real U_ref = 0.0735 * 0.5 * omega;                              // tip velocity
Real SOS = 10.0 * SMAX(U_ref, (Real)std::sqrt(2 * gravity_g * 0.09)); // numerical speed of sound

// material properties
Real rho = 1000.0; // reference density
Real K = 4.58;     // consistency index
Real n = 0.46;     // power index
Real tau_y = 18.9; // yield stress

Real min_shear_rate = 5e-2; // cutoff low shear rate
Real max_shear_rate = 1e+5; // cutoff high shear rate

// mesh geometry data
std::string full_path_to_shaft = "./input/shaft.stl";
std::string full_path_to_housing = "./input/housing.stl";
std::string full_path_to_fluid = "./input/fluid.stl";

Vecd translation(0.0, 0.0, 0.0);
Vecd housing_translation(0.0, 0.0, 0.0);
Real length_scale = 1.0;

//----------------------------------------------------------------------
//	Wall boundary geometry classes
//----------------------------------------------------------------------
class Mixer_Shaft : public ComplexShape
{
  public:
    explicit Mixer_Shaft(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_shaft, translation, length_scale);
    }
};
class Mixer_Housing : public ComplexShape
{
  public:
    explicit Mixer_Housing(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_housing, housing_translation, length_scale);
    }
};
class Fluid_Filling : public ComplexShape
{
  public:
    explicit Fluid_Filling(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_fluid, translation, length_scale, "OuterBoundary");
    }
};

// Setup output helper function
void output_setup()
{
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "XXXXXXXXXXXX Mixer Case XXXXXXXXXXXX" << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

    std::cout << "     particle_spacing= " << particle_spacing << std::endl;
    std::cout << "     end_time= " << end_time << std::endl;
    std::cout << "     K= " << K << std::endl;
    std::cout << "     n= " << n << std::endl;
    std::cout << "     tau_y= " << tau_y << std::endl;
    std::cout << "     min_shear_rate= " << min_shear_rate << std::endl;
    std::cout << "     max_shear_rate= " << max_shear_rate << std::endl;
    std::cout << "     rho= " << rho << std::endl;
    std::cout << "     u_ref= " << U_ref << std::endl;
    std::cout << "     SOS= " << SOS << std::endl;

    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
}

int main(int ac, char *av[])
{
    output_setup();
    //	Build up an SPHSystem
    BoundingBox system_domain_bounds(Vecd(-0.16, -0.16, -0.001), Vecd(0.16, 0.16, 0.32));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    sph_system.setRunParticleRelaxation(relaxation);
    sph_system.setReloadParticles(reload);

    //	Creating bodies with corresponding materials and particles
    FluidBody fluid(sph_system, makeShared<Fluid_Filling>("Fluid"));
    fluid.defineComponentLevelSetShape("OuterBoundary");
    fluid.defineParticlesAndMaterial<BaseParticles, HerschelBulkleyFluid>(rho, SOS, min_shear_rate, max_shear_rate, K, n, tau_y);
    sph_system.ReloadParticles()
        ? fluid.generateParticles<Reload>(fluid.getName())
        : fluid.generateParticles<Lattice>();

    SolidBody mixer_housing(sph_system, makeShared<Mixer_Housing>("Mixer_Housing"));
    mixer_housing.defineParticlesAndMaterial<SolidParticles, Solid>();
    sph_system.ReloadParticles()
        ? mixer_housing.generateParticles<Reload>(mixer_housing.getName())
        : mixer_housing.generateParticles<Lattice>();
    mixer_housing.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody mixer_shaft(sph_system, makeShared<Mixer_Shaft>("Mixer_Shaft"));
    mixer_shaft.defineAdaptationRatios(1.15, 2.0);
    mixer_shaft.defineBodyLevelSetShape();
    mixer_shaft.defineParticlesAndMaterial<SolidParticles, Solid>();
    sph_system.ReloadParticles()
        ? mixer_shaft.generateParticles<Reload>(mixer_shaft.getName())
        : mixer_shaft.generateParticles<Lattice>();
    mixer_shaft.addBodyStateForRecording<Vec3d>("NormalDirection");

    ObserverBody observer_body(sph_system, makeShared<Fluid_Filling>("ObserverBody"));
    observer_body.generateParticles<Lattice>();

    //	Define body relation map
    InnerRelation fluid_inner(fluid);
    InnerRelation shaft_inner(mixer_shaft);
    InnerRelation housing_inner(mixer_housing);
    ContactRelation fluid_observer_contact(observer_body, {&fluid});
    ContactRelation fluid_wall_contact(fluid, {&mixer_housing, &mixer_shaft});
    ContactRelation fluid_shaft_contact(fluid, {&mixer_shaft});
    ContactRelation shaft_fluid_contact(mixer_shaft, {&fluid});
    ContactRelation housing_fluid_contact(mixer_housing, {&fluid});

    ComplexRelation fluid_wall_complex(fluid_inner, fluid_wall_contact);

    //	Define the numerical methods used in the simulation
    Gravity gravity(Vec3d(0.0, 0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(fluid, gravity);

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(ConstructorArgs(fluid_inner, 0.3), fluid_wall_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> housing_normal_direction(mixer_housing);
    SimpleDynamics<NormalDirectionFromBodyShape> shaft_normal_direction(mixer_shaft);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWall<AcousticRiemannSolver, LinearGradientCorrection>> pressure_relaxation(fluid_inner, fluid_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall<AcousticRiemannSolver>> density_relaxation(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(fluid_inner, fluid_wall_contact);

    InteractionDynamics<fluid_dynamics::DistanceFromWall, SequencedPolicy> distance_to_wall(fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<LinearGradientCorrection>> vel_grad_calculation(fluid_inner, fluid_wall_contact);
    SimpleDynamics<fluid_dynamics::ShearRateDependentViscosity> shear_dependent_viscosity(fluid);
    InteractionWithUpdate<fluid_dynamics::NonNewtonianViscousForceWithWall<AngularConservative>> viscous_acceleration(fluid_inner, fluid_wall_contact);

    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_surface_indicator(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseTransportVelocityCorrectionComplex<SingleResolution, TruncatedLinear, NoKernelCorrection, BulkParticles>> transport_velocity_correction(fluid_inner, fluid_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(fluid, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_acoustic_time_step_size(fluid);
    ReduceDynamics<fluid_dynamics::SRDViscousTimeStepSize> get_viscous_time_step_size(fluid);

    ObservingAQuantity<Real> observing_viscosity(fluid_observer_contact, "VariableViscosity");
    SimpleDynamics<ParticleSnapshotAverage<Real>> average_viscosity(observer_body, "VariableViscosity");

    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_input_shaft_particles(mixer_shaft);
        SimpleDynamics<RandomizeParticlePosition> random_input_fluid_particles(fluid);
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner_shaft(shaft_inner);
        RelaxationStepLevelSetCorrectionComplex relaxation_step_complex_fluid(ConstructorArgs(fluid_inner, std::string("OuterBoundary")), fluid_shaft_contact);

        Real relax_time = 0.1;
        random_input_shaft_particles.exec(relax_time);
        random_input_fluid_particles.exec(relax_time);
        relaxation_step_inner_shaft.SurfaceBounding().exec();
        relaxation_step_complex_fluid.SurfaceBounding().exec();

        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 50;
        while (ite < relax_step)
        {
            // relaxation_step_inner_housing.exec();
            relaxation_step_inner_shaft.exec();
            relaxation_step_complex_fluid.exec();
            ite += 1;
        }
        std::cout << "Relaxation Complete" << std::endl;
    }

    //	Define the methods for I/O operations, observations
    fluid.addBodyStateForRecording<Real>("Pressure");
    fluid.addBodyStateForRecording<Real>("Density");
    fluid.addBodyStateForRecording<Real>("Mass");
    BodyStatesRecordingToVtp write_fluid_states(sph_system.real_bodies_);
    observer_body.addBodyStateForRecording<Real>("VariableViscosity");
    BodyStatesRecordingToVtp write_observation_states(observer_body);

    //----------------------------------------------------------------------
    //	Building Simbody.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    SimTK::GeneralForceSubsystem forces(MBsystem);
    SolidBodyPartForSimbody shaft_constraint_area(mixer_shaft, makeShared<TriangleMeshShapeSTL>(full_path_to_shaft, translation, length_scale));

    SimTK::Body::Rigid info(*shaft_constraint_area.body_part_mass_properties_);
    SimTK::MobilizedBody::Free mobBody(matter.updGround(), SimTK::Transform(), info, SimTK::Transform());

    SimTK::State state = MBsystem.realizeTopology();
    // Rigid Body Rotation
    mobBody.setOneU(state, 2, omega);

    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_tethered_spot(shaft_constraint_area, MBsystem, mobBody, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_rotation(shaft_constraint_area, MBsystem, mobBody, integ);

    //	Prepare the simulation
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    housing_normal_direction.exec();
    shaft_normal_direction.exec();
    distance_to_wall.exec();
    constant_gravity.exec();

    //	Setup for time-stepping control
    // size_t number_of_iterations = sph_system.RestartStep();
    int nmbr_of_outputs = 100;
    Real output_interval = end_time / nmbr_of_outputs;
    Real dt = 0;
    Real Dt = 0;
    Real Dt_visc = 0;
    Real Dt_adv = 0;
    Real Dt_aco = 0;
    int iteration = 0;
    int output_counter = 1;

    //	First output before the main loop.
    write_fluid_states.writeToFile(0);
    TickCount t1 = TickCount::now();
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        iteration++;
        TimeInterval tt;
        TickCount t2 = TickCount::now();
        tt = t2 - t1;
        Dt_adv = get_fluid_advection_time_step_size.exec();
        Dt_visc = get_viscous_time_step_size.exec();
        update_density_by_summation.exec();

        // Viscous Integration
        if (linearized_iteration == true && Dt_visc < Dt_adv && GlobalStaticVariables::physical_time_ < end_time * 0.001)
        {
            Real viscous_time = 0.0;
            corrected_configuration_fluid.exec();
            distance_to_wall.exec();
            free_surface_indicator.exec();
            vel_grad_calculation.exec();
            shear_dependent_viscosity.exec();

            // Viscous Substepping
            while (viscous_time < Dt_adv)
            {
                viscous_acceleration.exec(Dt_visc);
                viscous_time += Dt_visc;
                // Last substep
                if (viscous_time + Dt_visc > Dt_adv)
                {
                    Dt_visc = Dt_adv - viscous_time;
                }
            }
            transport_velocity_correction.exec(Dt);
        }
        else
        {
            Dt = SMIN(Dt_visc, Dt_adv);
            corrected_configuration_fluid.exec();
            distance_to_wall.exec();
            free_surface_indicator.exec(Dt);
            vel_grad_calculation.exec(Dt);
            shear_dependent_viscosity.exec(Dt);
            viscous_acceleration.exec(Dt);
            transport_velocity_correction.exec(Dt);
        }

        Real relaxation_time = 0.0;
        while (relaxation_time < Dt)
        {
            Dt_aco = get_acoustic_time_step_size.exec();
            dt = SMIN(Dt_aco, Dt);
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);
            relaxation_time += dt;
            GlobalStaticVariables::physical_time_ += dt;

            integ.stepBy(dt);
            constraint_rotation.exec();
        }

        if (iteration < 100 || output_counter * output_interval < GlobalStaticVariables::physical_time_)
        {
            std::cout << "Iteration: " << iteration << " | sim time in %: " << GlobalStaticVariables::physical_time_ / end_time * 100 << " | physical time in s: " << GlobalStaticVariables::physical_time_ << " | computation time in s: " << tt.seconds() << " | dt_adv: " << Dt_adv << " | dt_visc: " << Dt_visc << " | dt_aco: " << Dt_aco << "\r" << std::flush;
        }

        if (output_counter * output_interval < GlobalStaticVariables::physical_time_)
        {
            write_fluid_states.writeToFile();
            output_counter++;
        }
        fluid.updateCellLinkedListWithParticleSort(100);
        fluid_wall_complex.updateConfiguration();
    }
    TickCount t3 = TickCount::now();
    TimeInterval te;
    te = t3 - t1;
    std::cout << "Done with iterations: " << iteration << " | Total computation time in s: " << (t3 - t1).seconds() << std::endl;
    return 0;
}