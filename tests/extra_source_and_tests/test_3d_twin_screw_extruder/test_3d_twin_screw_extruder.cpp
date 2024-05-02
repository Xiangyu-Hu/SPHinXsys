/**
 * @file	test_3d_twin_screw_extruder.cpp
 * @brief	Co-Rotating Twin Screw Extruder
 * @details	This case is based on a Co-Rotating Twin Screw Extruder with periodic boundaries using a power law model
 * @author	Theodor Hennings
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// setup data
bool relaxation = true;
bool linearized_iteration = true;
bool reload = false;

// simulation setup
Real particle_spacing = 0.002; // particle spacing
Real gravity_g = 0;            // gravity
Real end_time = 1;             // end time
Real RPS = 1;                  // revolutions per second
Real omega = RPS * 3.14 * 2;   // angular velocity
Real U_ref = 0.03 * omega;     // tip velocity
Real SOS = 10.0 * U_ref;       // numerical speed of sound (0.287 is the fluid column height)

// material properties
Real rho = 1000.0; // reference density
Real K = 1;        // consistency index
Real n = 1.25;     // power index
Real tau_y = 0;    // yield stress

Real min_shear_rate = 5e-2; // cutoff low shear rate
Real max_shear_rate = 1e+3; // cutoff high shear rate

// mesh geometry data
std::string full_path_to_left_screw = "./input/left_screw.stl";

std::string full_path_to_right_screw = "./input/right_screw.stl";

std::string full_path_to_fluid = "./input/fluid.stl";

std::string full_path_to_barrel = "./input/barrel.stl";

Vecd translation(0.0, 0.0, 0.0);
Real length_scale = 1.0;

//----------------------------------------------------------------------
//	Wall boundary geometry classes
//----------------------------------------------------------------------
class Left_Screw : public ComplexShape
{
  public:
    explicit Left_Screw(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_left_screw, translation, length_scale);
    }
};
class Right_Screw : public ComplexShape
{
  public:
    explicit Right_Screw(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_right_screw, translation, length_scale);
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
class Barrel : public ComplexShape
{
  public:
    explicit Barrel(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_barrel, translation, length_scale);
    }
};

void output_setup()
{
    std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
    std::cout << "XXXXXXXXXXX Extruder Case XXXXXXXXXX" << std::endl;
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
    BoundingBox system_domain_bounds(Vecd(-0.036, -0.046, -0.011), Vecd(0.036, 0.093, 0.072));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.setRunParticleRelaxation(relaxation);
    sph_system.setReloadParticles(reload);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();

    //	Creating bodies with corresponding materials and particles
    FluidBody fluid(sph_system, makeShared<Fluid_Filling>("Fluid"));
    fluid.defineComponentLevelSetShape("OuterBoundary");
    fluid.defineParticlesAndMaterial<BaseParticles, HerschelBulkleyFluid>(rho, SOS, min_shear_rate, max_shear_rate, K, n, tau_y);
    sph_system.ReloadParticles()
        ? fluid.generateParticles<Reload>(fluid.getName())
        : fluid.generateParticles<Lattice>();

    SolidBody barrel(sph_system, makeShared<Barrel>("Barrel"));
    barrel.defineParticlesAndMaterial<SolidParticles, Solid>();
    sph_system.ReloadParticles()
        ? barrel.generateParticles<Reload>(barrel.getName())
        : barrel.generateParticles<Lattice>();
    barrel.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody left_screw(sph_system, makeShared<Left_Screw>("Left_Screw"));
    left_screw.defineParticlesAndMaterial<SolidParticles, Solid>();
    sph_system.ReloadParticles()
        ? left_screw.generateParticles<Reload>(left_screw.getName())
        : left_screw.generateParticles<Lattice>();
    left_screw.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody right_screw(sph_system, makeShared<Right_Screw>("Right_Screw"));
    right_screw.defineParticlesAndMaterial<SolidParticles, Solid>();
    sph_system.ReloadParticles()
        ? right_screw.generateParticles<Reload>(right_screw.getName())
        : right_screw.generateParticles<Lattice>();
    right_screw.addBodyStateForRecording<Vec3d>("NormalDirection");

    //	Define body relation map
    InnerRelation fluid_inner(fluid);
    InnerRelation left_screw_inner(left_screw);
    InnerRelation right_screw_inner(right_screw);
    InnerRelation barrel_inner(barrel);
    ContactRelation fluid_wall_contact(fluid, {&barrel, &left_screw, &right_screw});
    ContactRelation left_screw_fluid_contact(left_screw, {&fluid});
    ContactRelation right_screw_fluid_contact(right_screw, {&fluid});
    ContactRelation barrel_fluid_contact(barrel, {&fluid});

    ComplexRelation fluid_wall_complex(fluid_inner, fluid_wall_contact);

    //	Define the numerical methods used in the simulation
    Gravity gravity(Vec3d(0.0, 0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(fluid, gravity);
    PeriodicAlongAxis periodic_along_z(fluid.getSPHBodyBounds(), zAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition_z(fluid, periodic_along_z);

    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> corrected_configuration_fluid(ConstructorArgs(fluid_inner, 0.3), fluid_wall_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> left_screw_normal_direction(left_screw);
    SimpleDynamics<NormalDirectionFromBodyShape> right_screw_normal_direction(right_screw);
    SimpleDynamics<NormalDirectionFromBodyShape> barrel_normal_direction(barrel);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWall<AcousticRiemannSolver, LinearGradientCorrection>> pressure_relaxation(fluid_inner, fluid_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall<AcousticRiemannSolver>> density_relaxation(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(fluid_inner, fluid_wall_contact);

    InteractionDynamics<fluid_dynamics::DistanceFromWall, SequencedPolicy> distance_to_wall(fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::VelocityGradientWithWall<LinearGradientCorrection>> vel_grad_calculation(fluid_inner, fluid_wall_contact);
    SimpleDynamics<fluid_dynamics::ShearRateDependentViscosity> shear_dependent_viscosity(fluid);
    InteractionWithUpdate<fluid_dynamics::NonNewtonianViscousForceWithWall<AngularConservative>> viscous_acceleration(fluid_inner, fluid_wall_contact);

    //InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> free_surface_indicator(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::BaseTransportVelocityCorrectionComplex<SingleResolution, TruncatedLinear, NoKernelCorrection, BulkParticles>> transport_velocity_correction(fluid_inner, fluid_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(fluid, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_acoustic_time_step_size(fluid);
    ReduceDynamics<fluid_dynamics::SRDViscousTimeStepSize> get_viscous_time_step_size(fluid);

    //	Define the methods for I/O operations, observations
    BodyStatesRecordingToVtp write_fluid_states(sph_system.real_bodies_);
    fluid.addBodyStateForRecording<Real>("Pressure");
    fluid.addBodyStateForRecording<Real>("Density");
    fluid.addBodyStateForRecording<Real>("Mass");

    //----------------------------------------------------------------------
    //	Building Simbody.
    //----------------------------------------------------------------------
    SimTK::MultibodySystem MBsystem;
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    SimTK::GeneralForceSubsystem forces(MBsystem);
    SolidBodyPartForSimbody left_screw_constraint_area(left_screw, makeShared<TriangleMeshShapeSTL>(full_path_to_left_screw, translation, length_scale));
    SolidBodyPartForSimbody right_screw_constraint_area(right_screw, makeShared<TriangleMeshShapeSTL>(full_path_to_right_screw, translation, length_scale));

    SimTK::Body::Rigid info_left(*left_screw_constraint_area.body_part_mass_properties_);
    SimTK::Body::Rigid info_right(*right_screw_constraint_area.body_part_mass_properties_);

    SimTK::Vec3 rightScrewOrigin(0, 0.05, 0);
    SimTK::Transform offsetTransform(SimTK::Rotation(), rightScrewOrigin);

    SimTK::MobilizedBody::Free mobBody_left(matter.updGround(), SimTK::Transform(), info_left, SimTK::Transform());
    SimTK::MobilizedBody::Free mobBody_right(matter.updGround(), offsetTransform, info_right, SimTK::Transform());

    SimTK::State state = MBsystem.realizeTopology();
    // Rigid Body Rotations
    mobBody_left.setOneU(state, 2, -omega);
    mobBody_right.setOneU(state, 2, -omega);

    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);

    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_tethered_spot_left(left_screw_constraint_area, MBsystem, mobBody_left, integ);
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_tethered_spot_right(right_screw_constraint_area, MBsystem, mobBody_right, integ);

    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_rotation_left(left_screw_constraint_area, MBsystem, mobBody_left, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_rotation_right(right_screw_constraint_area, MBsystem, mobBody_right, integ);

    //	Prepare the simulation
    sph_system.initializeSystemCellLinkedLists();
    periodic_condition_z.update_cell_linked_list_.exec();
    sph_system.initializeSystemConfigurations();
    constant_gravity.exec();
    barrel_normal_direction.exec();
    left_screw_normal_direction.exec();
    right_screw_normal_direction.exec();
    distance_to_wall.exec();

    //	Setup for time-stepping control
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
        update_density_by_summation.exec(Dt);

        // Viscous Integration
        if (linearized_iteration == true && Dt_visc < Dt_adv && GlobalStaticVariables::physical_time_ < end_time * 0.001)
        {
            Real viscous_time = 0.0;
            corrected_configuration_fluid.exec();
            distance_to_wall.exec();
            //free_surface_indicator.exec();
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
            //free_surface_indicator.exec(Dt);
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
            constraint_rotation_left.exec();
            constraint_rotation_right.exec();
        }

        if (iteration < 100 || output_counter * output_interval < GlobalStaticVariables::physical_time_)
        {
            std::cout << std::fixed << std::setprecision(2) << std::scientific << "Iteration: " << iteration << " | sim time in %: " << GlobalStaticVariables::physical_time_ / end_time * 100 << " | physical time in s: " << GlobalStaticVariables::physical_time_ << " | computation time in s: " << tt.seconds() << " | dt_adv: " << Dt_adv << " | dt_visc: " << Dt_visc << " | dt_aco: " << Dt_aco << "\r" << std::flush;
        }

        if (output_counter * output_interval < GlobalStaticVariables::physical_time_)
        {
            write_fluid_states.writeToFile();
            output_counter++;
        }
        periodic_condition_z.bounding_.exec();
        fluid.updateCellLinkedListWithParticleSort(100);
        periodic_condition_z.update_cell_linked_list_.exec();
        fluid_wall_complex.updateConfiguration();
        left_screw.updateCellLinkedList();
        right_screw.updateCellLinkedList();
    }
    TickCount t3 = TickCount::now();
    TimeInterval te;
    te = t3 - t1;
    std::cout << "Done with iterations: " << iteration << " | Total computation time in s: " << (t3 - t1).seconds() << std::endl;
    return 0;
}