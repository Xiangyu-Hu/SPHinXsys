/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: Lid Driven Cavity 2.5D                     *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// setup properties
Real particle_spacing = 0.01;
Real end_time = 50.0;
int number_of_outputs = 500;

// non-Newtonian properties
Real K = 1;     // consistency index
Real n = 0.5;   // power index
Real tau_y = 0; // yield stress

Real min_shear_rate = 1e-2; // cutoff low shear rate
Real max_shear_rate = 1e+3; // cutoff high shear rate

// material properties
Real rho = 1000.0;       // reference density
Real u_lid = 1.0;        // lid velocity
Real SOS = 10.0 * u_lid; // numerical speed of sound

// geometry data
Real height = 1;
Real width = 1;
Real boundary_width = particle_spacing * 4; // boundary width

//----------------------------------------------------------------------
//	Complex shapes for wall boundary
//----------------------------------------------------------------------
class Lid_Boundary : public ComplexShape
{
  public:
    explicit Lid_Boundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container(0.5 * width + boundary_width, 0.5 * boundary_width);
        Transform translate_to_origin(scaled_container);
        Vecd transform(-boundary_width, height);
        Transform translate_to_position(transform + scaled_container);
        add<GeometricShapeBox>(Transform(translate_to_position), scaled_container);
    }
};
class No_Slip_Boundary : public ComplexShape
{
  public:
    explicit No_Slip_Boundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container_outer(0.5 * width + boundary_width, 0.5 * height + boundary_width);
        Vecd scaled_container(0.5 * width, 0.5 * height);
        Transform translate_to_origin_outer(Vec2d(-boundary_width, -boundary_width) + scaled_container_outer);
        Transform translate_to_origin_inner(scaled_container);

        add<GeometricShapeBox>(Transform(translate_to_origin_outer), scaled_container_outer);
        subtract<GeometricShapeBox>(Transform(translate_to_origin_inner), scaled_container);
    }
};
class FluidFilling : public ComplexShape
{
  public:
    explicit FluidFilling(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd scaled_container(0.5 * width, 0.5 * height);
        Transform translate_to_origin(scaled_container);
        add<GeometricShapeBox>(Transform(translate_to_origin), scaled_container);
    }
};

class BoundaryVelocity : public BodyPartMotionConstraint
{
  public:
    BoundaryVelocity(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part) {}

    void update(size_t index_i, Real dt = 0.0)
    {
        Vec2d velocity{Vecd::Zero()};
        velocity[0] = u_lid;
        vel_[index_i] = velocity;
    };
};

class ChangingBoundaryVelocity : public fluid_dynamics::FluidInitialCondition
{
  public:
    ChangingBoundaryVelocity(SPHBody &sph_body)
        : fluid_dynamics::FluidInitialCondition(sph_body) {};
    void update(size_t index_i, Real dt)
    {
        if (vel_[index_i][0] < 1)
        {
            vel_[index_i][0] += u_lid / 1000;
        }
        else
        {
            vel_[index_i][0] = u_lid;
        }
    }
};

void output_setup()
{
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "XXXXXX Lid Driven Cavity Case XXXXXX" << std::endl;
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;

    std::cout << "     particle_spacing= " << particle_spacing << std::endl;
    std::cout << "     end_time= " << end_time << std::endl;
    std::cout << "     K= " << K << std::endl;
    std::cout << "     n= " << n << std::endl;
    std::cout << "     tau_y= " << tau_y << std::endl;
    std::cout << "     min_shear_rate= " << min_shear_rate << std::endl;
    std::cout << "     max_shear_rate= " << max_shear_rate << std::endl;
    std::cout << "     rho= " << rho << std::endl;
    std::cout << "     u_lid= " << u_lid << std::endl;
    std::cout << "     SOS= " << SOS << std::endl;

    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
}

int main(int ac, char *av[])
{
    output_setup();
    //	Build up an SPHSystem
    BoundingBoxd system_domain_bounds(Vecd(-boundary_width * 2, -boundary_width * 2),
                                     Vecd(width + boundary_width * 2, height + boundary_width * 2));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av);

    //	Creating bodies with corresponding materials and particles
    FluidBody fluid(sph_system, makeShared<FluidFilling>("FluidBody"));
    fluid.defineClosure<WeaklyCompressibleFluid, HerschelBulkleyViscosity>(
        ConstructArgs(rho, SOS), ConstructArgs(min_shear_rate, max_shear_rate, K, n, tau_y));
    fluid.generateParticles<BaseParticles, Lattice>();

    SolidBody no_slip_boundary(sph_system, makeShared<No_Slip_Boundary>("NoSlipWall"));
    no_slip_boundary.defineMaterial<Solid>();
    no_slip_boundary.generateParticles<BaseParticles, Lattice>();

    ObserverBody observer_body(sph_system, makeShared<FluidFilling>("ObserverBody"));
    observer_body.generateParticles<BaseParticles, Lattice>();

    //	Define body relation map
    InnerRelation fluid_inner(fluid);
    ContactRelation fluid_all_walls(fluid, {&no_slip_boundary});
    ContactRelation fluid_observer_contact(observer_body, {&fluid});
    ComplexRelation fluid_walls_complex(fluid_inner, fluid_all_walls);

    //	Define the numerical methods used in the simulation
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(no_slip_boundary);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(DynamicsArgs(fluid_inner, 0.3), fluid_all_walls);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWall<AcousticRiemannSolver, LinearGradientCorrection>> pressure_relaxation(fluid_inner, fluid_all_walls);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall<NoRiemannSolver>> density_relaxation(fluid_inner, fluid_all_walls);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_inner, fluid_all_walls);

    InteractionDynamics<fluid_dynamics::DistanceFromWall> distance_to_wall(fluid_all_walls);
    InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<LinearGradientCorrection>> vel_grad_calculation(fluid_inner, fluid_all_walls);
    SimpleDynamics<fluid_dynamics::ShearRateDependentViscosity> shear_dependent_viscosity(fluid);
    InteractionWithUpdate<fluid_dynamics::NonNewtonianViscousForceWithWall<AngularConservative>> viscous_acceleration(fluid_inner, fluid_all_walls);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityLimitedCorrectionComplex<AllParticles>> transport_velocity_correction(fluid_inner, fluid_all_walls);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_fluid_advection_time_step_size(fluid, u_lid);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_acoustic_time_step_size(fluid);
    ReduceDynamics<fluid_dynamics::SRDViscousTimeStepSize> get_viscous_time_step_size(fluid);

    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(fluid_inner);
    BodyRegionByParticle lid_boundary(no_slip_boundary, makeShared<Lid_Boundary>("LidWall"));
    SimpleDynamics<BoundaryVelocity> lid_velocity(lid_boundary);
    ObservingAQuantity<Real> observing_viscosity(fluid_observer_contact, "VariableViscosity");
    SimpleDynamics<ParticleSnapshotAverage<Real>> average_viscosity(observer_body, "VariableViscosity");

    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(fluid);
    //	Define the methods for I/O operations, observations
    BodyStatesRecordingToVtp write_fluid_states(sph_system);
    write_fluid_states.addToWrite<Real>(fluid, "Pressure");
    write_fluid_states.addToWrite<Vecd>(no_slip_boundary, "NormalDirection");
    BodyStatesRecordingToVtp write_observation_states(observer_body);
    write_observation_states.addToWrite<Real>(observer_body, "VariableViscosity");
    //	Prepare the simulation
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    lid_velocity.exec();
    wall_boundary_normal_direction.exec();

    //	Setup for time-stepping control
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real output_interval = end_time / number_of_outputs;
    Real dt = 0;
    Real Dt = 0;
    Real Dt_visc = 0;
    Real Dt_adv = 0;
    Real Dt_aco = 0;
    int iteration = 0;
    int output_counter = 0;

    //	First output before the main loop.
    write_fluid_states.writeToFile(0);
    TickCount t1 = TickCount::now();
    while (physical_time < end_time)
    {
        iteration++;
        TimeInterval tt;
        TickCount t2 = TickCount::now();
        tt = t2 - t1;
        Dt_adv = get_fluid_advection_time_step_size.exec();
        Dt_visc = get_viscous_time_step_size.exec();
        Dt = SMIN(Dt_visc, Dt_adv);

        update_density_by_summation.exec();
        corrected_configuration_fluid.exec();
        distance_to_wall.exec();
        vel_grad_calculation.exec();
        shear_dependent_viscosity.exec();
        viscous_acceleration.exec();
        transport_velocity_correction.exec();

        Real relaxation_time = 0.0;
        while (relaxation_time < Dt)
        {
            Dt_aco = get_acoustic_time_step_size.exec();
            dt = SMIN(Dt_aco, Dt);
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);
            relaxation_time += dt;
            physical_time += dt;
        }

        if (iteration % 100 == 0 && iteration != 1)
        {
            particle_sorting.exec();
        }
        fluid.updateCellLinkedList();
        fluid_walls_complex.updateConfiguration();

        if (iteration % 100 == 0 || output_counter * output_interval < physical_time)
        {
            std::cout << "Iteration: " << iteration << " | sim time in %: " << physical_time / end_time * 100
                      << " | physical time in s: " << physical_time
                      << " | computation time in s: " << tt.seconds() << " | dt_adv: " << Dt_adv << " | dt_visc: " << Dt_visc
                      << " | dt_aco: " << Dt_aco << "\n"
                      << std::flush;

            if (output_counter > number_of_outputs * 9 / 10)
            {
                fluid_observer_contact.updateConfiguration();
                observing_viscosity.exec();
                average_viscosity.exec();
            }
        }

        if (output_counter * output_interval < physical_time)
        {
            compute_vorticity.exec();
            write_fluid_states.writeToFile();
            write_observation_states.writeToFile();
            output_counter++;
        }
    }
    TickCount t3 = TickCount::now();
    TimeInterval te;
    te = t3 - t1;
    std::cout << "Done with iterations: " << iteration << " | Total computation time in s: " << (t3 - t1).seconds() << std::endl;
    return 0;
}