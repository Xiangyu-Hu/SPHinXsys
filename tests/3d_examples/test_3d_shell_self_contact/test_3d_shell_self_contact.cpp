/**
 * @file 	shell_self_contact.cpp
 * @brief 	Self contact of a cylinder shell
 * @details This is a case to test shell self contact formulation.
 * @author 	Weiyi Kong, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real scale = 0.001;
const Real radius = 0.5 * scale; // outer surface radius
const Real length = 4.0 * scale;
const Real resolution_ref = 0.05 * scale;
const Real thickness = resolution_ref; // thickness of the balloon
const BoundingBox system_domain_bounds(Vec3d(-0.5 * length, -radius, -radius),
                                       Vec3d(0.5 * length, radius, radius));
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
const Real rho0_s = 1000.0;      /** Normalized density. */
const Real Youngs_modulus = 5e3; /** Normalized Young's modulus. */
const Real poisson = 0.45;       /** Poisson ratio. */
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * thickness;
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
const Real radius_mid_surface = radius - 0.5 * thickness; // mid surface radius
const int particle_number_mid_surface = int(2.0 * radius_mid_surface * Pi / resolution_ref);
const int particle_number_height = int(length / resolution_ref);

class ShellParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit ShellParticleGenerator(SPHBody &sph_body)
        : SurfaceParticleGenerator(sph_body){};
    void initializeGeometricVariables() override
    {
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            for (int j = 0; j < particle_number_height; j++)
            {
                Real theta = (i + 0.5) * 2 * Pi / (Real)particle_number_mid_surface;
                Real x = -0.5 * length + length * j / (Real)particle_number_height + 0.5 * resolution_ref;
                Real y = radius_mid_surface * cos(theta);
                Real z = radius_mid_surface * sin(theta);
                initializePositionAndVolumetricMeasure(Vec3d(x, y, z),
                                                       resolution_ref * resolution_ref);
                const Vec3d n_0 = Vec3d(0.0, y / radius_mid_surface, z / radius_mid_surface);
                initializeSurfaceProperties(n_0, thickness);
            }
        }
    }
};
//----------------------------------------------------------------------
//	Define boundary conditions
//----------------------------------------------------------------------
/** Define the displacement-control geometry. */
class UpperBoundaryGeometry : public BodyPartByParticle
{
  public:
    UpperBoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&UpperBoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][1] > radius - resolution_ref &&
            std::abs(base_particles_.pos_[index_i][0]) < radius &&
            std::abs(base_particles_.pos_[index_i][2]) < radius)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
// class LowerBoundaryGeometry : public BodyPartByParticle
// {
//   public:
//     LowerBoundaryGeometry(SPHBody &body, const std::string &body_part_name)
//         : BodyPartByParticle(body, body_part_name)
//     {
//         TaggingParticleMethod tagging_particle_method = std::bind(&LowerBoundaryGeometry::tagManually, this, _1);
//         tagParticles(tagging_particle_method);
//     };

//   private:
//     void tagManually(size_t index_i)
//     {
//         if (base_particles_.pos_[index_i][1] < -radius + resolution_ref &&
//             std::abs(base_particles_.pos_[index_i][0]) < 2 * resolution_ref &&
//             std::abs(base_particles_.pos_[index_i][2]) < 2 * resolution_ref)
//         {
//             body_part_particles_.push_back(index_i);
//         }
//     };
// };
class UpperDisplacement : public thin_structure_dynamics::ConstrainShellBodyRegion
{
  public:
    explicit UpperDisplacement(BodyPartByParticle &body_part)
        : ConstrainShellBodyRegion(body_part){};

  protected:
    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i] = Vec3d(0, -0.1, 0);
    };
};
// class LowerDisplacement : public thin_structure_dynamics::ConstrainShellBodyRegion
// {
//   public:
//     explicit LowerDisplacement(BodyPartByParticle &body_part)
//         : ConstrainShellBodyRegion(body_part), acc_prior_(particles_->acc_prior_){};

//   protected:
//     StdLargeVec<Vecd> &acc_prior_;
//     void update(size_t index_i, Real dt = 0.0)
//     {
//         acc_prior_[index_i] = Vec3d(0, 300, 0);
//     };
// };
class ConstrainAlongYAxis : public LocalDynamics, public thin_structure_dynamics::ShellDataSimple
{
  public:
    explicit ConstrainAlongYAxis(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          thin_structure_dynamics::ShellDataSimple(sph_body),
          vel_(particles_->vel_){};

  protected:
    StdLargeVec<Vecd> &vel_;
    double vel_y_avg_;
    void setupDynamics(Real dt = 0.0)
    {
        vel_y_avg_ = 0;
        SPH::particle_for(
            SPH::par,
            particles_->total_real_particles_,
            [&, this](size_t index_i)
            {
                vel_y_avg_ += vel_[index_i][1];
            });
        vel_y_avg_ /= particles_->total_real_particles_;
    }
    void update(size_t index_i, Real dt = 0.0)
    {
        vel_[index_i][1] -= vel_y_avg_;
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
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody shell(sph_system, makeShared<DefaultShape>("Shell"));
    shell.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    shell.generateParticles<ShellParticleGenerator>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation shell_inner(shell);
    ShellSelfContactRelation shell_self_contact(shell);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Corrected configuration. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    /** Normal update. */
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_normal_update(shell);
    /** Shell curvature. */
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> shell_curvature(shell_inner);
    /** Define external force.*/
    SimpleDynamics<TimeStepInitialization> shell_initialize_timestep(shell);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_get_time_step_size(shell);
    /** stress relaxation for the walls. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first_half(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second_half(shell_inner);
    /** Algorithms for shell self contact. */
    InteractionDynamics<solid_dynamics::ShellSelfContactDensitySummation> shell_self_contact_density(shell_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> shell_self_contact_forces(shell_self_contact);
    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        shell_position_damping(0.2, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        shell_rotation_damping(0.2, shell_inner, "AngularVelocity", physical_viscosity);
    /** Boundary condition*/
    UpperBoundaryGeometry upper_bc_geometry(shell, "UpperBcGeometry");
    SimpleDynamics<UpperDisplacement> upper_dis(upper_bc_geometry);
    // LowerBoundaryGeometry lower_bc_geometry(shell, "LowerBcGeometry");
    // SimpleDynamics<LowerDisplacement> lower_dis(lower_bc_geometry);
    SimpleDynamics<ConstrainAlongYAxis> constrain_y_axis(shell);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    shell_corrected_configuration.exec();
    shell_curvature.compute_initial_curvature();
    /** Initial states output. */
    shell.addBodyStateForRecording<Real>("1stPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("2ndPrincipleCurvature");
    body_states_recording.writeToFile(0);

    /** Main loop. */
    int ite = 0;
    Real end_time = 0.015;
    Real output_interval = end_time / 100;
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
            shell_initialize_timestep.exec();
            if (ite % 100 == 0)
            {
                std::cout << "N=" << ite << " Time: "
                          << GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
            }

            // lower_dis.exec();
            // upper_dis.exec();

            shell_self_contact_density.exec();
            shell_self_contact_forces.exec();

            dt = shell_get_time_step_size.exec();

            shell_stress_relaxation_first_half.exec(dt);
            upper_dis.exec();
            constrain_y_axis.exec();
            shell_position_damping.exec(dt);
            shell_rotation_damping.exec(dt);
            upper_dis.exec();
            constrain_y_axis.exec();
            shell_stress_relaxation_second_half.exec(dt);

            shell_normal_update.exec();
            shell_curvature.exec();
            shell.updateCellLinkedList();
            shell_self_contact.updateConfiguration();

            ite++;
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
