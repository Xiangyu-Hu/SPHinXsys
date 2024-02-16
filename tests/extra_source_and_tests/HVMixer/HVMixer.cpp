/**
 * @file	HVMixer.cpp
 * @brief	Non-Newtonian Mixer
 * @details	This case is based on a Rushton turbine impeller in a cylindrical tank fitted with vertical baffles
 * @author	Theodor Hennings
 */
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// setup data
Real particle_spacing = 0.001;
Real gravity_g = 9.81;
Real end_time = 0.5;

// material properties
Real rho = 1000.0; // reference density
Real omega = 7 * 3.14 * 2;
Real U_ref = 0.0735 * 0.5 * omega;
Real SOS = 10.0 * U_ref; /**< Reference sound speed. */

// non-Newtonian properties
Real K = 4.58;     // consistency index
Real n = 0.46;     // power index
Real tau_y = 18.9; // yield stress

Real min_shear_rate = 1e-3; // cutoff low shear rate
Real max_shear_rate = 1e+5; // cutoff high shear rate

// mesh geometry data
// std::string full_path_to_shaft = "./input/Shaft_Fusion.stl";
std::string full_path_to_shaft = "./input/Shaft_Fusion_Large_Blade.stl";

std::string full_path_to_housing = "./input/Housing_Fusion_2.stl";

std::string full_path_to_fluid = "./input/Fluid_Reduced_Height.stl";

std::string full_path_to_refinement = "./input/Refinement.stl";

Vecd translation(0.0, 0.0, 0.0);
Vecd housing_translation(0.0, 0.0, 0.0);
Real length_scale = 1.0;

Eigen::AngleAxis<Real> angle(omega / 1000, Eigen::Vector3d::UnitZ());
Transform rotation_transform(angle);

//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
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
        if (true)
        {
            add<TriangleMeshShapeSTL>(full_path_to_housing, housing_translation, length_scale);
        }
        else
        {
            add<ExtrudeShape<TriangleMeshShapeSTL>>(4.0 * particle_spacing, full_path_to_housing, housing_translation, length_scale);
            subtract<TriangleMeshShapeSTL>(full_path_to_housing, housing_translation, length_scale);
        }
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
class Refinement : public ComplexShape
{
  public:
    explicit Refinement(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_refinement, translation, length_scale);
    }
};

int main(int ac, char *av[])
{
    //	Build up an SPHSystem
    BoundingBox system_domain_bounds(Vecd(-0.085, -0.085, -0.01), Vecd(0.085, 0.085, 0.22));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    sph_system.setRunParticleRelaxation(true);

    //	Creating bodies with corresponding materials and particles
    FluidBody fluid(sph_system, makeShared<Fluid_Filling>("Fluid"));
    fluid.defineComponentLevelSetShape("OuterBoundary");
    fluid.defineParticlesAndMaterial<BaseParticles, HerschelBulkleyFluid>(rho, SOS, min_shear_rate, max_shear_rate, K, n, tau_y);
    fluid.generateParticles<ParticleGeneratorLattice>();

    SolidBody mixer_housing(sph_system, makeShared<Mixer_Housing>("Mixer_Housing"));
    mixer_housing.defineParticlesAndMaterial<SolidParticles, Solid>();
    mixer_housing.generateParticles<ParticleGeneratorLattice>();
    mixer_housing.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody mixer_shaft(sph_system, makeShared<Mixer_Shaft>("Mixer_Shaft"));
    mixer_shaft.defineBodyLevelSetShape();
    mixer_shaft.defineParticlesAndMaterial<SolidParticles, Solid>();
    mixer_shaft.generateParticles<ParticleGeneratorLattice>();
    mixer_shaft.addBodyStateForRecording<Vec3d>("NormalDirection");

    //	Define body relation map
    InnerRelation fluid_inner(fluid);
    InnerRelation shaft_inner(mixer_shaft);
    InnerRelation housing_inner(mixer_housing);
    ContactRelation fluid_wall_contact(fluid, {&mixer_housing, &mixer_shaft});
    ContactRelation fluid_shaft_contact(fluid, {&mixer_shaft});
    ContactRelation shaft_fluid_contact(mixer_shaft, {&fluid});
    ContactRelation housing_fluid_contact(mixer_housing, {&fluid});

    ComplexRelation fluid_wall_complex(fluid_inner, fluid_wall_contact);

    //	Define the numerical methods used in the simulation
    Gravity gravity(Vec3d(0.0, 0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(fluid, gravity);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(fluid_inner, fluid_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(fluid_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(fluid_inner, fluid_wall_contact);

    InteractionDynamics<fluid_dynamics::VelocityGradientWithWall> vel_grad_calculation(fluid_inner, fluid_wall_contact);
    InteractionDynamics<fluid_dynamics::ShearRateDependentViscosity> shear_rate_calculation(fluid_inner);
    InteractionWithUpdate<fluid_dynamics::GeneralizedNewtonianViscousForceWithWall> viscous_acceleration(fluid_inner, fluid_wall_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(fluid, U_ref);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(fluid);
    ReduceDynamics<fluid_dynamics::SRDViscousTimeStepSize> get_viscous_time_step_size(fluid);

    SimpleDynamics<NormalDirectionFromBodyShape> housing_normal_direction(mixer_housing);
    SimpleDynamics<NormalDirectionFromBodyShape> shaft_normal_direction(mixer_shaft);

    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation.
        //----------------------------------------------------------------------
        // SimpleDynamics<SPH::relax_dynamics::RandomizeParticlePosition> random_input_housing_particles(mixer_housing);
        SimpleDynamics<SPH::relax_dynamics::RandomizeParticlePosition> random_input_shaft_particles(mixer_shaft);
        SimpleDynamics<SPH::relax_dynamics::RandomizeParticlePosition> random_input_fluid_particles(fluid);
        // relax_dynamics::RelaxationStepLevelSetCorrectionInner relaxation_step_inner_housing(housing_inner);
        relax_dynamics::RelaxationStepLevelSetCorrectionInner relaxation_step_inner_shaft(shaft_inner);
        relax_dynamics::RelaxationStepLevelSetCorrectionComplex relaxation_step_complex_fluid(ConstructorArgs(fluid_inner, "OuterBoundary"), fluid_shaft_contact);

        Real relax_time = 0.1;
        random_input_shaft_particles.exec(relax_time);
        // random_input_housing_particles.exec(relax_time);
        random_input_fluid_particles.exec(relax_time);
        relaxation_step_inner_shaft.SurfaceBounding().exec();
        // relaxation_step_inner_housing.SurfaceBounding().exec();
        relaxation_step_complex_fluid.SurfaceBounding().exec();
        //? mixer_shaft.updateCellLinkedList();
        //? mixer_housing.updateCellLinkedList();

        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 10;
        while (ite < relax_step)
        {
            // relaxation_step_inner_housing.exec();
            relaxation_step_inner_shaft.exec();
            relaxation_step_complex_fluid.exec();
            ite += 1;
        }
        std::cout << "The relaxation process of shaft & housing particles finished !" << std::endl;
    }

    //	Define the methods for I/O operations, observations
    BodyStatesRecordingToVtp write_fluid_states(sph_system.real_bodies_);
    fluid.addBodyStateForRecording<Real>("Pressure");
    fluid.addBodyStateForRecording<Real>("Density");

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
    constant_gravity.exec();
    housing_normal_direction.exec();
    shaft_normal_direction.exec();

    //	Setup for time-stepping control
    // size_t number_of_iterations = sph_system.RestartStep();
    int nmbr_of_outputs = 100;
    Real output_interval = end_time / nmbr_of_outputs;
    Real dt = 0.0;
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
        Real Dt = get_fluid_advection_time_step_size.exec();
        Real Dt_visc = get_viscous_time_step_size.exec();
        Dt = SMIN(Dt_visc, Dt);

        update_density_by_summation.exec(Dt);

        vel_grad_calculation.exec(Dt);
        shear_rate_calculation.exec(Dt);
        viscous_acceleration.exec(Dt);

        Real relaxation_time = 0.0;
        while (relaxation_time < Dt)
        {
            dt = SMIN(dt, Dt);
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);
            dt = get_fluid_time_step_size.exec();
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
        shaft_fluid_contact.updateConfiguration();
        housing_fluid_contact.updateConfiguration();
        mixer_shaft.updateCellLinkedList();
    }
    TickCount t3 = TickCount::now();
    TimeInterval te;
    te = t3 - t1;
    std::cout << "Done with iterations: " << iteration << " | Total computation time in s: " << (t3 - t1).seconds() << std::endl;
    return 0;
}