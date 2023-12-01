/**
 * @file 	shell_self_contact.cpp
 * @brief 	Self contact of a shell balloon and its contact with a rigid solid plate
 * @details This is a case to test shell self contact formulation.
 * @author 	Weiyi Kong, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real scale = 0.001;
const Real radius_balloon = 0.5 * scale; // outer surface radius
const Real DL_balloon = 4.0 * scale;
const Real resolution_ref = 0.05 * scale;
const Real thickness_balloon = resolution_ref; // thickness of the balloon
const Real level_set_refinement_ratio = resolution_ref / (0.1 * thickness_balloon);
const Real DL = DL_balloon + 2 * radius_balloon;
const Real DH = 2 * radius_balloon;
const BoundingBox system_domain_bounds(Vec2d(-0.5 * DL, 0), Vec2d(0.5 * DL, 2 * DH));
const Vec2d balloon_center(0, radius_balloon);
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real rho0_s = 1000.0;      /** Normalized density. */
Real Youngs_modulus = 5e3; /** Normalized Young's modulus. */
Real poisson = 0.45;       /** Poisson ratio. */
Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * thickness_balloon;
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
/** create a shell outer shape */
std::vector<Vecd> createShellOuterShape()
{
    // geometry
    std::vector<Vecd> outer_shape;
    outer_shape.push_back(Vecd(-0.5 * DL_balloon, 0.0));
    outer_shape.push_back(Vecd(-0.5 * DL_balloon, 2 * radius_balloon));
    outer_shape.push_back(Vecd(0.5 * DL_balloon, 2 * radius_balloon));
    outer_shape.push_back(Vecd(0.5 * DL_balloon, 0.0));
    outer_shape.push_back(Vecd(-0.5 * DL_balloon, 0.0));

    return outer_shape;
}
/** create a shell inner shape */
std::vector<Vecd> createShellInnerShape()
{
    // geometry
    std::vector<Vecd> inner_shape;
    inner_shape.push_back(Vecd(-0.5 * DL_balloon, thickness_balloon));
    inner_shape.push_back(Vecd(-0.5 * DL_balloon, 2 * radius_balloon - thickness_balloon));
    inner_shape.push_back(Vecd(0.5 * DL_balloon, 2 * radius_balloon - thickness_balloon));
    inner_shape.push_back(Vecd(0.5 * DL_balloon, thickness_balloon));
    inner_shape.push_back(Vecd(-0.5 * DL_balloon, thickness_balloon));

    return inner_shape;
}
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class Shell : public MultiPolygonShape
{
  public:
    explicit Shell(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        const Vec2d circle_center_1 = balloon_center - Vec2d(0.5 * DL_balloon, 0);
        const Vec2d circle_center_2 = balloon_center + Vec2d(0.5 * DL_balloon, 0);

        multi_polygon_.addAPolygon(createShellOuterShape(), ShapeBooleanOps::add);
        multi_polygon_.addACircle(circle_center_1, radius_balloon, 100, ShapeBooleanOps::add);
        multi_polygon_.addACircle(circle_center_2, radius_balloon, 100, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createShellInnerShape(), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_1, radius_balloon - thickness_balloon, 100, ShapeBooleanOps::sub);
        multi_polygon_.addACircle(circle_center_2, radius_balloon - thickness_balloon, 100, ShapeBooleanOps::sub);
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
        if (base_particles_.pos_[index_i][1] > 2 * radius_balloon - resolution_ref &&
            base_particles_.pos_[index_i][0] < 2 * resolution_ref &&
            base_particles_.pos_[index_i][0] > -2 * resolution_ref)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class LowerBoundaryGeometry : public BodyPartByParticle
{
  public:
    LowerBoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&LowerBoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][1] < resolution_ref &&
            base_particles_.pos_[index_i][0] < 2 * resolution_ref &&
            base_particles_.pos_[index_i][0] > -2 * resolution_ref)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class UpperDisplacement : public thin_structure_dynamics::ConstrainShellBodyRegion
{
  public:
    explicit UpperDisplacement(BodyPartByParticle &body_part)
        : ConstrainShellBodyRegion(body_part), acc_prior_(particles_->acc_prior_){};

  protected:
    StdLargeVec<Vecd> &acc_prior_;
    void update(size_t index_i, Real dt = 0.0)
    {
        acc_prior_[index_i] = Vec2d(0, -30);
    };
};
class LowerDisplacement : public thin_structure_dynamics::ConstrainShellBodyRegion
{
  public:
    explicit LowerDisplacement(BodyPartByParticle &body_part)
        : ConstrainShellBodyRegion(body_part), acc_prior_(particles_->acc_prior_){};

  protected:
    StdLargeVec<Vecd> &acc_prior_;
    void update(size_t index_i, Real dt = 0.0)
    {
        acc_prior_[index_i] = Vec2d(0, 30);
    };
};
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
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody shell(sph_system, makeShared<Shell>("Shell"));
    shell.defineAdaptation<SPHAdaptation>(1.15, 1.0);
    // here dummy linear elastic solid is use because no solid dynamics in particle relaxation
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
    {
        shell.generateParticles<ParticleGeneratorReload>(io_environment, shell.getName());
    }
    else
    {
        shell.defineBodyLevelSetShape(level_set_refinement_ratio)->writeLevelSet(io_environment);
        shell.generateParticles<ThickSurfaceParticleGeneratorLattice>(thickness_balloon);
    }

    if (!sph_system.RunParticleRelaxation() && !sph_system.ReloadParticles())
    {
        std::cout << "Error: This case requires reload shell particles for simulation!" << std::endl;
        return 0;
    }
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
        relax_dynamics::ShellRelaxationStep relaxation_step_shell_inner(shell_inner);
        relax_dynamics::ShellNormalDirectionPrediction shell_normal_prediction(shell_inner, thickness_balloon, cos(Pi / 3.75));
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

        relaxation_step_shell_inner.MidSurfaceBounding().exec();
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
    InteractionDynamics<solid_dynamics::SelfContactDensitySummation> shell_self_contact_density(shell_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> shell_self_contact_forces(shell_self_contact);
    /** Damping with the solid body*/
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_position_damping(0.2, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_rotation_damping(0.2, shell_inner, "AngularVelocity", physical_viscosity);
    /** Boundary condition*/
    UpperBoundaryGeometry upper_bc_geometry(shell, "UpperBcGeometry");
    SimpleDynamics<UpperDisplacement> upper_dis(upper_bc_geometry);
    LowerBoundaryGeometry lower_bc_geometry(shell, "LowerBcGeometry");
    SimpleDynamics<LowerDisplacement> lower_dis(lower_bc_geometry);
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
    shell_curvature.exec();
    /** Initial states output. */
    shell.addBodyStateForRecording<Real>("TotalMeanCurvature");
    body_states_recording.writeToFile(0);
    /** Main loop. */
    int ite = 0;
    Real end_time = 0.05;
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

            lower_dis.exec();
            upper_dis.exec();

            shell_self_contact_density.exec();
            shell_self_contact_forces.exec();

            dt = shell_get_time_step_size.exec();

            shell_stress_relaxation_first_half.exec(dt);
            constrain_y_axis.exec();
            shell_position_damping.exec(dt);
            shell_rotation_damping.exec(dt);
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
