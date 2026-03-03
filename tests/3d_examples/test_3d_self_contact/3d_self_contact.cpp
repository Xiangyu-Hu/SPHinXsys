/**
 * @file 	3d_self_contact.cpp
 * @brief 	This is the test to check self contact for solid dynamics
 * @author 	Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string full_path_to_file = "./input/coil.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real half_width = 55.0;
Real global_resolution = half_width / 30.0;
Real BW = global_resolution * 4;
Vec3d domain_lower_bound(-half_width - BW, -half_width - 1.5 * BW, -BW);
Vec3d domain_upper_bound(half_width + BW, half_width + BW, 2.0 * half_width + BW);
// Domain bounds of the system.
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
Real rho0_s = 1.265;
Real poisson = 0.45;
Real Youngs_modulus = 5e4;
Real physical_viscosity = 200.0;
//----------------------------------------------------------------------
//	Body shapes used in the case.
//----------------------------------------------------------------------
class Coil : public ComplexShape
{
  public:
    explicit Coil(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_file, Vecd::Zero(), 1.0);
    }
};
class StationaryPlate : public ComplexShape
{
  public:
    explicit StationaryPlate(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_plate(half_width + BW, 0.5 * BW, half_width + BW);
        Vecd translation_plate(0.0, -half_width - 0.75 * BW, half_width);
        add<GeometricShapeBox>(Transform(translation_plate), halfsize_plate);
    }
};
//----------------------------------------------------------------------
//	The main program
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    // Tag for run particle relaxation for the initial body fitted distribution.
    sph_system.setRunParticleRelaxation(false);
    // Tag for reload initially relaxed particles.
    sph_system.setReloadParticles(false);
#ifdef BOOST_AVAILABLE
    // handle command line arguments
    sph_system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody coil(sph_system, makeShared<Coil>("Coil"));
    coil.defineBodyLevelSetShape()->writeLevelSet();
    coil.defineMaterial<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? coil.generateParticles<BaseParticles, Reload>(coil.getName())
        : coil.generateParticles<BaseParticles, Lattice>();

    SolidBody stationary_plate(sph_system, makeShared<StationaryPlate>("StationaryPlate"));
    stationary_plate.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    stationary_plate.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation coil_inner(coil);
    SelfSurfaceContactRelation coil_self_contact(coil);
    SurfaceContactRelation coil_contact(coil_self_contact, {&stationary_plate});
    //----------------------------------------------------------------------
    //	Define simple file input and outputs functions.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states(sph_system);
    //----------------------------------------------------------------------
    //	check whether run particle relaxation for body fitted particle distribution.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(coil);
        // Write the particle reload files.
        ReloadParticleIO write_particle_reload_files(coil);
        // A  Physics relaxation step.
        RelaxationStepInner relaxation_step_inner(coil_inner);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_states.writeToFile(0);
        //----------------------------------------------------------------------
        //	Particle relaxation loop.
        //----------------------------------------------------------------------
        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
                write_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        // Output particles position for reload.
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	This section define all numerical methods will be used in this case.
    //----------------------------------------------------------------------
    Gravity gravity(Vec3d(0.0, -1.0, 0.0));
    SimpleDynamics<GravityForce<Gravity>> coil_constant_gravity(coil, gravity);
    // Corrected configuration for reproducing rigid rotation.
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> corrected_configuration(coil_inner);

    // stress relaxation.
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> stress_relaxation_first_half(coil_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> stress_relaxation_second_half(coil_inner);
    // Algorithms for solid-solid contacts.
    InteractionDynamics<solid_dynamics::ContactFactorSummation> coil_update_contact_density(coil_contact);
    InteractionWithUpdate<solid_dynamics::ContactForceFromWall> coil_compute_solid_contact_forces(coil_contact);
    InteractionDynamics<solid_dynamics::SelfContactFactorSummation> coil_self_contact_density(coil_self_contact);
    InteractionWithUpdate<solid_dynamics::SelfContactForce> coil_self_contact_forces(coil_self_contact);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> computing_time_step_size(coil);

    // Damping the velocity field for quasi-static solution
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d, FixedDampingRate>>>
        coil_damping(0.2, coil_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<TotalKineticEnergy>> write_coil_kinetic_energy(coil);
    //----------------------------------------------------------------------
    //	From here the time stepping begins.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    corrected_configuration.exec();
    coil_constant_gravity.exec();
    write_states.writeToFile(0);
    // Setup time stepping control parameters.
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int ite = 0;
    Real end_time = 10.0;
    Real output_period = end_time / 100.0;
    Real dt = 0.0;
    // Statistics for computing time.
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_period)
        {
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << physical_time << "	dt: "
                          << dt << "\n";
                write_coil_kinetic_energy.writeToFile(ite);
            }
            // contact dynamics.
            coil_self_contact_density.exec();
            coil_self_contact_forces.exec();
            coil_update_contact_density.exec();
            coil_compute_solid_contact_forces.exec();
            // Stress relaxation and damping.
            stress_relaxation_first_half.exec(dt);
            coil_damping.exec(dt);
            stress_relaxation_second_half.exec(dt);

            ite++;
            dt = computing_time_step_size.exec();
            integration_time += dt;
            physical_time += dt;

            // update particle neighbor relations for contact dynamics
            coil.updateCellLinkedList();
            coil_self_contact.updateConfiguration();
            coil_contact.updateConfiguration();
        }
        TickCount t2 = TickCount::now();
        write_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_coil_kinetic_energy.generateDataBase(1.0e-3);
    }
    else
    {
        write_coil_kinetic_energy.testResult();
    }

    return 0;
}
