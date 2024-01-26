#include "sphinxsys.h"
using namespace SPH;

class CheckKernelCompleteness
{
  private:
    BaseParticles *particles_;
    Kernel *kernel_;
    std::vector<SPH::BaseParticles *> contact_particles_;
    ParticleConfiguration *inner_configuration_;
    std::vector<ParticleConfiguration *> contact_configuration_;

    StdLargeVec<Real> W_ijV_j_ttl;
    StdLargeVec<Real> W_ijV_j_ttl_contact;
    StdLargeVec<Vecd> dW_ijV_je_ij_ttl;
    StdLargeVec<int> number_of_inner_neighbor;
    StdLargeVec<int> number_of_contact_neighbor;

  public:
    CheckKernelCompleteness(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : particles_(&inner_relation.base_particles_),
          kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
          inner_configuration_(&inner_relation.inner_configuration_)
    {
        for (size_t i = 0; i != contact_relation.contact_bodies_.size(); ++i)
        {
            contact_particles_.push_back(&contact_relation.contact_bodies_[i]->getBaseParticles());
            contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
        }
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl, "TotalKernel");
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl_contact, "TotalKernelContact");
        inner_relation.base_particles_.registerVariable(dW_ijV_je_ij_ttl, "TotalKernelGrad");
        inner_relation.base_particles_.registerVariable(number_of_inner_neighbor, "InnerNeighborNumber");
        inner_relation.base_particles_.registerVariable(number_of_contact_neighbor, "ContactNeighborNumber");
    }

    inline void exec()
    {
        particle_for(
            par,
            particles_->total_real_particles_,
            [&, this](size_t index_i)
            {
                int N_inner_number = 0;
                int N_contact_number = 0;
                Real W_ijV_j_ttl_i = particles_->Vol_[index_i] * kernel_->W(0, ZeroVecd);
                Vecd dW_ijV_je_ij_ttl_i = Vecd::Zero();
                const Neighborhood &inner_neighborhood = (*inner_configuration_)[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    W_ijV_j_ttl_i += inner_neighborhood.W_ij_[n] * particles_->Vol_[index_j];
                    dW_ijV_je_ij_ttl_i += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                    N_inner_number++;
                }

                double W_ijV_j_ttl_contact_i = 0;
                for (size_t k = 0; k < contact_configuration_.size(); ++k)
                {
                    const SPH::Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
                    for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = wall_neighborhood.j_[n];
                        W_ijV_j_ttl_contact_i += wall_neighborhood.W_ij_[n] * contact_particles_[k]->Vol_[index_j];
                        dW_ijV_je_ij_ttl_i += wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
                        N_contact_number++;
                    }
                }
                W_ijV_j_ttl[index_i] = W_ijV_j_ttl_i + W_ijV_j_ttl_contact_i;
                W_ijV_j_ttl_contact[index_i] = W_ijV_j_ttl_contact_i;
                dW_ijV_je_ij_ttl[index_i] = dW_ijV_je_ij_ttl_i;
                number_of_inner_neighbor[index_i] = N_inner_number;
                number_of_contact_neighbor[index_i] = N_contact_number;
            });
    }
};

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const bool if_fsi = true;
const Real unit_scale = 0.001;               // mm to m
const Real fluid_diameter = 12 * unit_scale; /**<Inlet and outlet diameter. */
const Real fluid_radius = 0.5 * fluid_diameter;
const Real full_length = fluid_diameter * 7.5;           /**<Total length og fluid. */
const Real shell_thickness = 0.6 * unit_scale;           /**<Balloon thickness. */
const Real resolution_fluid = fluid_diameter / Real(15); /**< Global reference resolution. */
const Real resolution_shell = 0.4 * unit_scale;          /**< Balloon resolution. */
const Real wall_thickness = 3.2 * unit_scale;            /**< Wall boundary thickness. */
const Real inflow_length = 20 * resolution_fluid;        /**< Inflow region. */

const Real balloon_fixed_length_distal = 2.15 * unit_scale;
const Real balloon_fixed_length_proximal = 1.5 * unit_scale;
const Real balloon_full_length = 45.5 * unit_scale;
const Vecd balloon_translation(fluid_diameter * 3.5, 0, 0);
const Real balloon_distal_end = -5.5 * unit_scale + balloon_translation.x();

const Vec3d emitter_halfsize(resolution_fluid * 2, fluid_radius, fluid_radius);
const Vec3d emitter_translation(resolution_fluid * 2, 0., 0.);
const Vec3d buffer_halfsize(inflow_length * 0.5, fluid_radius + 2 * resolution_fluid, fluid_radius + 2 * resolution_fluid);
const Vec3d buffer_translation(inflow_length * 0.5 - 2 * resolution_fluid, 0., 0.);
const Vec3d disposer_halfsize(resolution_fluid * 2, fluid_radius * 1.1, fluid_radius * 1.1);
const Vec3d disposer_translation(full_length - disposer_halfsize[0], 0., 0.);

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-wall_thickness, -fluid_radius - wall_thickness, -fluid_radius - wall_thickness),
                                 Vec3d(full_length + wall_thickness, fluid_radius + shell_thickness, fluid_radius + shell_thickness));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1056.0;                           // blood density
const Real mu_f = 3.5e-3;                             // blood viscosity
const Real Re = 200;                                  /**< Reynolds number. */
const Real U_f = Re * mu_f / rho0_f / fluid_diameter; /**< Average velocity at throat. */
const Real U_max = 4.0 * U_f;
const Real c_f = 10.0 * U_max; /**< Speed of sound. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
const Real rho0_s = 1250.0; // https://plastena.lt/en/rubber-silicone-rubber-vmq
const Real hardness = 50;   // Durometer hardnes: 50A
const Real youngs_modulus =
    std::pow(10, 0.0235 * hardness - 0.6403) * 1e3; // actual: 1e6, ref: https://doi.org/10.5254/1.3547752 eq. 12A
const Real poisson_ratio = 0.495;
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * youngs_modulus) * shell_thickness;
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
const std::string path_to_fluid_file = "./input/fluid.stl";
const std::string path_to_balloon_srf_file = "./input/balloon_outer_srf_dp_0_4.stl";
class FluidBlock : public ComplexShape
{
  public:
    explicit FluidBlock(const std::string &shape_name)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_fluid_file, Vecd::Zero(), unit_scale);
        subtract<TriangleMeshShapeSTL>(path_to_balloon_srf_file, balloon_translation, unit_scale);
    }
};
const std::string path_to_wall_file = "./input/wall_3_2.stl";
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_wall_file, Vecd::Zero(), unit_scale);
    }
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(1.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        // Real run_time = GlobalStaticVariables::physical_time_;
        //  Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        // Real u_ave = 0.5 * u_ref_ * (1 - cos(2.0 * Pi * run_time / t_ref_));
        Real radius2 = position[1] * position[1] + position[2] * position[2];
        target_velocity[0] = 2.0 * u_ref_ * SMAX(0.0, 1.0 - radius2 / fluid_radius / fluid_radius);

        return target_velocity;
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

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][0] < balloon_distal_end + balloon_fixed_length_distal ||
            base_particles_.pos_[index_i][0] > balloon_distal_end + balloon_full_length - balloon_fixed_length_proximal)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
class ForceBoundaryGeometry : public BodyPartByParticle
{
  public:
    ForceBoundaryGeometry(SPHBody &body, const std::string &body_part_name)
        : BodyPartByParticle(body, body_part_name)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&ForceBoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    void tagManually(size_t index_i)
    {
        if (base_particles_.pos_[index_i][0] < balloon_distal_end + balloon_fixed_length_distal + 3.5 * unit_scale ||
            base_particles_.pos_[index_i][0] > balloon_distal_end + balloon_full_length - balloon_fixed_length_proximal)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
const Real t_ref = 0.1;
const Real balloon_force = 5;
class BalloonForce : public thin_structure_dynamics::ConstrainShellBodyRegion
{
  public:
    explicit BalloonForce(BodyPartByParticle &body_part)
        : ConstrainShellBodyRegion(body_part),
          force_prior_(particles_->force_prior_),
          pos0_(particles_->pos0_),
          n_(particles_->n_){};

  protected:
    StdLargeVec<Vecd> &force_prior_;
    StdLargeVec<Vecd> &pos0_;
    StdLargeVec<Vecd> &n_;
    void update(size_t index_i, Real dt = 0.0)
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real x_ratio = (pos0_[index_i].x() - (balloon_distal_end + balloon_fixed_length_distal)) / (balloon_full_length - balloon_fixed_length_proximal - balloon_fixed_length_distal);
        Real force_avg = 0.5 * balloon_force * (1 - cos(Pi * (2.0 * run_time / t_ref + x_ratio)));
        force_prior_[index_i] = -force_avg * n_[index_i];
    };
};
int main(int ac, char *av[])
{
    std::cout << "U_max = " << U_max << std::endl;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_fluid);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setRunParticleRelaxation(false);
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidBlock>("fluid"));
    // fluid_block.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet(0);
    fluid_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    fluid_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("wall_3_2"));
    wall_boundary.defineBodyLevelSetShape();
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>(); // dummy material parameters
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? wall_boundary.generateParticles<ParticleGeneratorReload>(io_environment, wall_boundary.getName())
        : wall_boundary.generateParticles<ParticleGeneratorLattice>();

    SolidBody shell(sph_system, makeShared<DefaultShape>("balloon_0_4"));
    shell.defineAdaptation<SPHAdaptation>(1.15, resolution_fluid / resolution_shell);
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    shell.generateParticles<ParticleGeneratorReload>(io_environment, shell.getName());
    auto update_shell_thickness = [&]()
    {
        auto &thickness = *shell.getBaseParticles().getVariableByName<Real>("Thickness");
        particle_for(
            par,
            shell.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                thickness[index_i] = shell_thickness;
            });
    };
    update_shell_thickness();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation wall_boundary_inner(wall_boundary);
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> random_wall_boundary_particles(wall_boundary);
        relax_dynamics::RelaxationStepInner relaxation_step_inner(wall_boundary_inner);
        BodyStatesRecordingToVtp write_wall_boundary_to_vtp(io_environment, {&wall_boundary});
        ReloadParticleIO write_particle_reload_files(io_environment, {&wall_boundary});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_wall_boundary_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        write_wall_boundary_to_vtp.writeToFile(0);
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
                write_wall_boundary_to_vtp.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
        /** Output results. */
        write_particle_reload_files.writeToFile(0);
        return 0;
    }
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    // Must construct ShellCurvature before ShellContactRelation
    InnerRelation fluid_inner(fluid_block);
    InnerRelation shell_inner(shell);
    ContactRelation fluid_wall_contact(fluid_block, {&wall_boundary});
    ContactRelationToShell fluid_shell_contact(fluid_block, {&shell});
    ContactRelationFromShell shell_fluid_contact(shell, {&fluid_block});
    ComplexRelation fluid_block_complex(fluid_inner, {&fluid_wall_contact, &fluid_shell_contact});
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell, fluid_block);
    ShellSelfContactRelation shell_self_contact(shell);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Algorithm for fluid dynamics. */
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization(fluid_block);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(fluid_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(fluid_block);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeStream>, Contact<>, Contact<>>> update_fluid_density_by_summation(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> fluid_pressure_relaxation(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, NoRiemannSolver>> fluid_density_relaxation(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    InteractionDynamics<ComplexInteraction<fluid_dynamics::ViscousAcceleration<Inner<>, Contact<Wall>, Contact<Wall>>>> viscous_acceleration(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    InteractionWithUpdate<ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>, Contact<>>>> inlet_outlet_surface_particle_indicator(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::TransportVelocityCorrection<Inner<SingleResolution>, Contact<Boundary>, Contact<Boundary>>, NoKernelCorrection, BulkParticles>> transport_velocity_correction(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    /** Algorithm for in-/outlet. */
    BodyAlignedBoxByParticle emitter(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec3d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);
    BodyAlignedBoxByCell buffer(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec3d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(buffer);
    BodyAlignedBoxByCell disposer(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec3d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, xAxis);
    /** Algorithm for solid dynamics. */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(shell_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell);
    /** Algorithms for shell self contact. */
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> shell_curvature(shell_inner);
    InteractionDynamics<solid_dynamics::ShellSelfContactDensitySummation> shell_self_contact_density(shell_self_contact);
    InteractionDynamics<solid_dynamics::SelfContactForce> shell_self_contact_forces(shell_self_contact);
    auto reset_force_prior = [&]()
    {
        const auto &force_from_fluid = *shell.getBaseParticles().getVariableByName<Vecd>("AllForceFromFluid");
        particle_for(
            par,
            shell.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                shell.getBaseParticles().force_prior_[index_i] = force_from_fluid[index_i];
            });
    };
    auto update_shell_volume = [&]()
    {
        particle_for(
            par,
            shell.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                shell.getBaseParticles().Vol_[index_i] = shell.getBaseParticles().mass_[index_i] / shell.getBaseParticles().rho_[index_i];
            });
    };
    /** FSI */
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_shell(shell_fluid_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_shell_update(shell_fluid_contact, viscous_force_on_shell);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(shell);
    /** constraint and damping */
    BoundaryGeometry shell_boundary_geometry(shell, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> shell_constrain(shell_boundary_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        shell_position_damping(0.2, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec3d>>>
        shell_rotation_damping(0.2, shell_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    fluid_block.addBodyStateForRecording<Real>("Pressure");
    fluid_block.addBodyStateForRecording<int>("Indicator");
    shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("Average2ndPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("1stPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("2ndPrincipleCurvature");
    shell.addBodyStateForRecording<Vecd>("AllForceFromFluid");
    shell.addBodyStateForRecording<Vecd>("PriorForce");
    shell.addBodyStateForRecording<Real>("Density");
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    shell_corrected_configuration.exec();
    // shell_curvature.compute_initial_curvature();
    shell_average_curvature.exec();
    fluid_block_complex.updateConfiguration();
    shell_fluid_contact.updateConfiguration();
    //   Check dWijVjeij
    CheckKernelCompleteness check_kernel_completeness(fluid_inner, fluid_shell_contact);
    check_kernel_completeness.exec();
    fluid_block.addBodyStateForRecording<Real>("TotalKernel");
    fluid_block.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    fluid_block.addBodyStateForRecording<int>("InnerNeighborNumber");
    fluid_block.addBodyStateForRecording<int>("ContactNeighborNumber");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 10;
    Real end_time = 2.0;
    Real output_interval = end_time / 200.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                           /**< Default acoustic time step sizes. */
    Real dt_s = 0.0;                         /**< Default acoustic time step sizes for solid. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    const double Dt_ref = fluid_advection_time_step.exec();
    const double dt_ref = fluid_acoustic_time_step.exec();
    const double dt_s_ref = shell_time_step_size.exec();
    auto run_simulation = [&]()
    {
        std::cout << "Simulation starts here" << std::endl;
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < output_interval)
            {
                fluid_step_initialization.exec();
                Real Dt = fluid_advection_time_step.exec();
                if (Dt < Dt_ref / 20)
                {
                    std::cout << "Dt = " << Dt << ", Dt_ref = " << Dt_ref << std::endl;
                    std::cout << "Advective time step decreased too much!" << std::endl;
                    throw std::runtime_error("Advective time step decreased too much!");
                }

                inlet_outlet_surface_particle_indicator.exec();
                update_fluid_density_by_summation.exec();
                viscous_acceleration.exec();
                transport_velocity_correction.exec();

                if (if_fsi)
                    viscous_force_on_shell.exec();

                /** Dynamics including pressure relaxation. */
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    Real dt_temp = fluid_acoustic_time_step.exec();
                    if (dt_temp < dt_ref / 20)
                    {
                        std::cout << "dt = " << dt_temp << ", dt_ref = " << dt_ref << std::endl;
                        std::cout << "Acoustic time step decreased too much!" << std::endl;
                        throw std::runtime_error("Acoustic time step decreased too much!");
                    }
                    dt = SMIN(dt_temp, Dt - relaxation_time);
                    fluid_pressure_relaxation.exec(dt);
                    emitter_buffer_inflow_condition.exec();
                    if (if_fsi)
                        fluid_force_on_shell_update.exec();
                    fluid_density_relaxation.exec(dt);

                    /** Solid dynamics time stepping. */
                    if (if_fsi)
                    {
                        Real dt_s_sum = 0.0;
                        average_velocity_and_acceleration.initialize_displacement_.exec();
                        while (dt_s_sum < dt)
                        {
                            reset_force_prior();

                            shell_self_contact_density.exec();
                            shell_self_contact_forces.exec();

                            Real dt_s_temp = shell_time_step_size.exec();
                            if (dt_s_temp < dt_s_ref / 100)
                            {
                                std::cout << "dt_s = " << dt_s_temp << ", dt_s_ref = " << dt_s_ref << std::endl;
                                std::cout << "Shell time step decreased too much!" << std::endl;
                                throw std::runtime_error("Shell time step decreased too much!");
                            }
                            dt_s = std::min(dt_s_temp, dt - dt_s_sum);
                            shell_stress_relaxation_first.exec(dt_s);
                            shell_constrain.exec();
                            shell_position_damping.exec(dt_s);
                            shell_rotation_damping.exec(dt_s);
                            shell_constrain.exec();
                            shell_stress_relaxation_second.exec(dt_s);

                            update_shell_volume();
                            shell_update_normal.exec();
                            // shell_curvature.exec();
                            shell.updateCellLinkedList();
                            shell_self_contact.updateConfiguration();

                            dt_s_sum += dt_s;
                        }
                        average_velocity_and_acceleration.update_averages_.exec(dt);
                    }

                    relaxation_time += dt;
                    integration_time += dt;
                    GlobalStaticVariables::physical_time_ += dt;
                }

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << GlobalStaticVariables::physical_time_
                              << "	Dt = " << Dt << "	dt = " << dt;
                    if (if_fsi)
                        std::cout << "  dt_s = " << dt_s;
                    std::cout << "\n";
                }
                number_of_iterations++;

                /** inflow injection*/
                emitter_inflow_injection.exec();
                disposer_outflow_deletion.exec();

                /** Update cell linked list and configuration. */
                fluid_block.updateCellLinkedList();
                if (if_fsi)
                {
                    update_shell_volume();
                    shell_update_normal.exec();
                    shell.updateCellLinkedList();
                    shell_curvature_inner.updateConfiguration();
                    shell_average_curvature.exec();
                    shell_fluid_contact.updateConfiguration();
                }
                fluid_block_complex.updateConfiguration();

                // check_kernel_completeness.exec();
                // fluid_block.setNewlyUpdated();
                // shell.setNewlyUpdated();
                // write_real_body_states.writeToFile();
            }

            TickCount t2 = TickCount::now();
            check_kernel_completeness.exec();
            write_real_body_states.writeToFile();
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
            // exit(0);
        }
        TickCount t4 = TickCount::now();

        TimeInterval tt;
        tt = t4 - t1 - interval;
        std::cout << "Total wall time for computation: " << tt.seconds()
                  << " seconds." << std::endl;
    };

    try
    {
        run_simulation();
    }
    catch (const std::exception &e)
    {
        std::cout << "Error catched..." << std::endl;
        fluid_block.setNewlyUpdated();
        shell.setNewlyUpdated();
        write_real_body_states.writeToFile(1e8);
    }
    return 0;
}
