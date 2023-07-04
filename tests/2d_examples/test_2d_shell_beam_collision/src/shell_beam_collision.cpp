/**
 * @file 	shell_beam_collision.cpp
 * @brief 	A rigid shell rigid box hitting an elastic wall boundary
 * @details This is a case to test shell contact formulations in a reverse way (shell to elastic solid).
 * @author 	Massoud Rezavand, Virtonomy GmbH and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0;                        /**< box length. */
Real DH = 4.0;                        /**< box height. */
Real resolution_ref = 0.025;          /**< reference resolution. */
Real BW = resolution_ref * 4.;        /**< wall width for BCs. */
Real thickness = resolution_ref * 1.; /**< shell thickness. */
Real level_set_refinement_ratio = resolution_ref / (0.1 * thickness);
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
Vec2d circle_center(2.0, 2.0);
Real circle_radius = 0.5;
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1.0;               /** Normalized density. */
Real Youngs_modulus = 5e3;       /** Normalized Young's modulus. */
Real poisson = 0.45;             /** Poisson ratio. */
Real physical_viscosity = 200.0; /** physical damping, here we choose the same value as numerical viscosity. */
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class Beam : public MultiPolygonShape
{
  public:
    explicit Beam(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        std::vector<Vecd> outer_beam_shape;
        outer_beam_shape.push_back(Vecd(-BW, -BW));
        outer_beam_shape.push_back(Vecd(-BW, DH + BW));
        outer_beam_shape.push_back(Vecd(0.0, DH + BW));
        outer_beam_shape.push_back(Vecd(0.0, -BW));
        outer_beam_shape.push_back(Vecd(-BW, -BW));

        multi_polygon_.addAPolygon(outer_beam_shape, ShapeBooleanOps::add);
    }
};
class Shell : public ComplexShape
{
  public:
    explicit Shell(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBall>(circle_center, circle_radius + resolution_ref);
        subtract<GeometricShapeBall>(circle_center, circle_radius);
    }
};
//----------------------------------------------------------------------
//	define the beam base which will be constrained.
//----------------------------------------------------------------------
MultiPolygon createBeamConstrainShape()
{
    // a beam base shape
    std::vector<Vecd> bottom_beam_base_shape{
        Vecd(-1.5 * BW, -1.5 * BW), Vecd(-1.5 * BW, 0.5 * resolution_ref),
        Vecd(0.5 * resolution_ref, 0.5 * resolution_ref),
        Vecd(0.5 * resolution_ref, -1.5 * BW), Vecd(-1.5 * BW, -1.5 * BW)};

    std::vector<Vecd> top_beam_base_shape{
        Vecd(-1.5 * BW, DH - 0.5 * resolution_ref), Vecd(-1.5 * BW, DH + 1.5 * BW),
        Vecd(0.5 * resolution_ref, DH + 1.5 * BW),
        Vecd(0.5 * resolution_ref, DH - 0.5 * resolution_ref), Vecd(-1.5 * BW, DH - 0.5 * resolution_ref)};

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(bottom_beam_base_shape, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(top_beam_base_shape, ShapeBooleanOps::add);
    return multi_polygon;
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
    sph_system.setRunParticleRelaxation(true);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody shell(sph_system, makeShared<Shell>("Shell"));
    shell.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    // here dummy linear elastic solid is use because no solid dynamics in particle relaxation
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1.0, 1.0, 0.0);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        shell.generateParticles<ParticleGeneratorReload>(io_environment, shell.getName());
    }
    else
    {
        shell.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(io_environment);
        shell.generateParticles<ThickSurfaceParticleGeneratorLattice>(thickness);
    }

    if (!sph_system.RunParticleRelaxation() && !sph_system.ReloadParticles())
    {
        std::cout << "Error: This case requires reload shell particles for simulation!" << std::endl;
        return 0;
    }

    SolidBody beam(sph_system, makeShared<Beam>("Beam"));
    beam.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    beam.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation beam_inner(beam);
    SurfaceContactRelation shell_contact(shell, {&beam});
    SurfaceContactRelation beam_contact(beam, {&shell});
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
        relax_dynamics::ShellRelaxationStepInner
            relaxation_step_shell_inner(shell_inner, thickness, level_set_refinement_ratio);
        relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(shell_inner, thickness, cos(Pi / 3.75));
        shell.addBodyStateForRecording<int>("UpdatedIndicator");
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

        relaxation_step_shell_inner.mid_surface_bounding_.exec();
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
        shell_normal_prediction.exec();
        write_relaxed_particles.writeToFile(ite);
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Define external force.*/
    SimpleDynamics<TimeStepInitialization> beam_initialize_timestep(beam);
    InteractionWithUpdate<CorrectedConfigurationInner> beam_corrected_configuration(beam_inner);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> shell_get_time_step_size(beam);
    /** stress relaxation for the walls. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> beam_stress_relaxation_first_half(beam_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> beam_stress_relaxation_second_half(beam_inner);
    /** Algorithms for shell-solid contact. */
    InteractionDynamics<solid_dynamics::ContactDensitySummation> beam_shell_update_contact_density(beam_contact);
    InteractionDynamics<solid_dynamics::ContactForceFromWall> beam_compute_solid_contact_forces(beam_contact);
    InteractionDynamics<solid_dynamics::ContactForceToWall> shell_compute_solid_contact_forces(shell_contact);
    BodyRegionByParticle holder(beam, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> constraint_holder(holder);
    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        beam_damping(0.5, beam_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    /**
     * The multi body system from simbody.
     */
    SimTK::MultibodySystem MBsystem;
    /** The bodies or matter of the MBsystem. */
    SimTK::SimbodyMatterSubsystem matter(MBsystem);
    /** The forces of the MBsystem.*/
    SimTK::GeneralForceSubsystem forces(MBsystem);
    /** Mass properties of the rigid shell box. */
    SolidBodyPartForSimbody shell_multibody(shell, makeShared<Shell>("Shell"));
    SimTK::Body::Rigid rigid_info(*shell_multibody.body_part_mass_properties_);
    SimTK::MobilizedBody::Slider
        shellMBody(matter.Ground(), SimTK::Transform(SimTKVec3(0)), rigid_info, SimTK::Transform(SimTKVec3(0)));
    /** Gravity. */
    SimTK::Force::UniformGravity sim_gravity(forces, matter, SimTKVec3(Real(-150.), 0.0, 0.0));
    /** discrete forces acting on the bodies. */
    SimTK::Force::DiscreteForces force_on_bodies(forces, matter);
    /** Time stepping method for multibody system.*/
    SimTK::State state = MBsystem.realizeTopology();
    SimTK::RungeKuttaMersonIntegrator integ(MBsystem);
    integ.setAccuracy(1e-3);
    integ.setAllowInterpolation(false);
    integ.initialize(state);
    /** Coupling between SimBody and SPH.*/
    ReduceDynamics<solid_dynamics::TotalForceOnBodyPartForSimBody>
        force_on_shell(shell_multibody, MBsystem, shellMBody, integ);
    SimpleDynamics<solid_dynamics::ConstraintBodyPartBySimBody>
        constraint_shell(shell_multibody, MBsystem, shellMBody, integ);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    beam_corrected_configuration.exec();
    /** Initial states output. */
    body_states_recording.writeToFile(0);
    /** Main loop. */
    int ite = 0;
    Real T0 = 1.0;
    Real end_time = T0;
    Real output_interval = 0.01 * T0;
    Real dt = 0.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            beam_initialize_timestep.exec();
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
            }
            beam_shell_update_contact_density.exec();
            beam_compute_solid_contact_forces.exec();
            shell_compute_solid_contact_forces.exec();

            {
                SimTK::State &state_for_update = integ.updAdvancedState();
                force_on_bodies.clearAllBodyForces(state_for_update);
                force_on_bodies.setOneBodyForce(state_for_update, shellMBody, force_on_shell.exec());
                integ.stepBy(dt);
                constraint_shell.exec();
            }

            beam_stress_relaxation_first_half.exec(dt);
            constraint_holder.exec(dt);
            beam_damping.exec(dt);
            constraint_holder.exec(dt);
            beam_stress_relaxation_second_half.exec(dt);

            shell.updateCellLinkedList();
            shell_contact.updateConfiguration();
            beam.updateCellLinkedList();
            beam_contact.updateConfiguration();

            ite++;
            Real dt_free = shell_get_time_step_size.exec();
            dt = dt_free;
            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
        }
        TickCount t2 = TickCount::now();
        body_states_recording.writeToFile(ite);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
    return 0;
}
