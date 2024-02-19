/**
 * @file 	self_contact.cpp
 * @brief 	This is the case file for the test of dynamic self contact.
 * @author   Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;

//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
const Real scale = 0.001;
const Real radius_shell = 3.5 * scale;
const Real resolution_shell = 0.075 * scale;
const Real resolution_plate = 0.075 * scale;
const Real thickness_shell = 0.15 * scale; // thickness of the balloon
const Real half_length_plate = 2.5 * radius_shell;
const Real thickness_plate = 0.6 * scale; // thickness of the plate
const Real gravity_g = 1.0;
const BoundingBox system_domain_bounds(Vec2d(-half_length_plate, -resolution_plate), Vec2d(half_length_plate, resolution_plate));
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
// shell
const Real rho0_s = 1.0e3;         // reference density
const Real Youngs_modulus = 2.0e3; // reference Youngs modulus
const Real poisson = 0.45;         // Poisson ratio
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * thickness_shell * thickness_shell;
// plate
const Real Youngs_modulus_plate = 1e2; // reference Youngs modulus
//----------------------------------------------------------------------
//	Geometric elements used in the case.
//----------------------------------------------------------------------
class PlateParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit PlateParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        Real y = -0.5 * resolution_plate;
        Real x = -half_length_plate + 0.5 * resolution_plate;
        while (x < half_length_plate)
        {
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_plate);
            initializeSurfaceProperties(-Vecd(0, 1), thickness_plate);
            x += resolution_plate;
        }
    }
};

class ShellParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit ShellParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        Vecd shell_center(0, radius_shell + 0.5 * resolution_shell);
        int number_of_particles_circle = int(M_PI * radius_shell / resolution_shell) + 1;
        Real dtheta = M_PI / Real(number_of_particles_circle - 1);
        // circle
        Real theta = 0;
        for (int i = 0; i < number_of_particles_circle; i++)
        {
            Real x1 = radius_shell * cos(theta);
            Real y1 = -radius_shell * sin(theta);
            // left circle
            initializePositionAndVolumetricMeasure(shell_center + Vecd(x1, y1), resolution_shell);
            initializeSurfaceProperties(Vecd(cos(theta), -sin(theta)), thickness_shell);
            theta += dtheta;
        }
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
        if (base_particles_.pos_[index_i][0] <= -half_length_plate + 4 * resolution_plate ||
            base_particles_.pos_[index_i][0] >= half_length_plate - 4 * resolution_plate)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_shell);
// handle command line arguments
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av);
#endif
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody shell(sph_system, makeShared<DefaultShape>("Shell"));
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    shell.generateParticles<ShellParticleGenerator>();

    SolidBody plate(sph_system, makeShared<DefaultShape>("Plate"));
    plate.defineAdaptation<SPHAdaptation>(1.15, resolution_shell / resolution_plate);
    plate.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus_plate, poisson);
    plate.generateParticles<PlateParticleGenerator>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation shell_inner(shell);
    InnerRelation plate_inner(plate);
    SurfaceContactRelationFromShellToShell shell_plate_contact(shell, {&plate}, {false});
    SurfaceContactRelationFromShellToShell plate_shell_contact(plate, {&shell}, {true});
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell, plate);
    ShellInnerRelationWithContactKernel plate_curvature_inner(plate, shell);
    //-----------------------------------------------------------------------------
    // this section define all numerical methods will be used in this case
    //-----------------------------------------------------------------------------
    // shell
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(shell, gravity);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(shell_inner);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_position_damping(0.5, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_rotation_damping(0.5, shell_inner, "AngularVelocity", physical_viscosity);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell);
    // plate
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> plate_corrected_configuration(plate_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> plate_time_step_size(plate);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> plate_stress_relaxation_first(plate_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> plate_stress_relaxation_second(plate_inner);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        plate_position_damping(0.5, plate_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        plate_rotation_damping(0.5, plate_inner, "AngularVelocity", physical_viscosity);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> plate_average_curvature(plate_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> plate_update_normal(plate);
    BoundaryGeometry plate_boundary_geometry(plate, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> plate_constrain(plate_boundary_geometry);
    // contact
    InteractionDynamics<solid_dynamics::ContactDensitySummation> shell_update_contact_density(shell_plate_contact);
    InteractionDynamics<solid_dynamics::ContactDensitySummation> plate_update_contact_density(plate_shell_contact);
    InteractionWithUpdate<solid_dynamics::ShellContactForce> shell_compute_solid_contact_forces(shell_plate_contact);
    InteractionWithUpdate<solid_dynamics::ShellContactForce> plate_compute_solid_contact_forces(plate_shell_contact);
    //-----------------------------------------------------------------------------
    //	outputs
    //-----------------------------------------------------------------------------
    IOEnvironment io_environment(sph_system);
    shell.addBodyStateForRecording<Real>("RepulsionDensity");
    shell.addBodyStateForRecording<Vecd>("RepulsionForce");
    plate.addBodyStateForRecording<Real>("RepulsionDensity");
    plate.addBodyStateForRecording<Vecd>("RepulsionForce");
    shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    plate.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    BodyStatesRecordingToVtp body_states_recording(sph_system.real_bodies_);
    //-----------------------------------------------------------------------------
    //	Setup particle configuration and initial conditions
    //-----------------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    plate_corrected_configuration.exec();
    shell_corrected_configuration.exec();
    constant_gravity.exec();
    shell_average_curvature.exec();
    plate_average_curvature.exec();
    plate_shell_contact.updateConfiguration();
    shell_plate_contact.updateConfiguration();
    //-----------------------------------------------------------------------------
    // from here the time stepping begins
    //-----------------------------------------------------------------------------
    // starting time zero
    GlobalStaticVariables::physical_time_ = 0.0;
    body_states_recording.writeToFile(0);

    int ite = 0;
    Real T0 = 1.0;
    Real end_time = T0;
    Real output_interval = 0.01 * T0;
    Real Dt = 0.1 * output_interval;
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
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                if (ite % 100 == 0)
                {
                    std::cout << "N=" << ite << " Time: "
                              << GlobalStaticVariables::physical_time_ << "	dt: " << dt
                              << "\n";
                }
                shell_update_contact_density.exec();
                plate_update_contact_density.exec();
                shell_compute_solid_contact_forces.exec();
                plate_compute_solid_contact_forces.exec();

                dt = std::min(shell_time_step_size.exec(), plate_time_step_size.exec());

                shell_stress_relaxation_first.exec(dt);
                plate_stress_relaxation_first.exec(dt);
                plate_constrain.exec();
                shell_position_damping.exec(dt);
                shell_rotation_damping.exec(dt);
                plate_position_damping.exec(dt);
                plate_rotation_damping.exec(dt);
                plate_constrain.exec();
                shell_stress_relaxation_second.exec(dt);
                plate_stress_relaxation_second.exec(dt);

                shell.updateCellLinkedList();
                plate.updateCellLinkedList();
                shell_update_normal.exec();
                shell_curvature_inner.updateConfiguration();
                shell_average_curvature.exec();
                plate_update_normal.exec();
                plate_curvature_inner.updateConfiguration();
                plate_average_curvature.exec();
                shell_plate_contact.updateConfiguration();
                plate_shell_contact.updateConfiguration();

                ite++;

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
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
