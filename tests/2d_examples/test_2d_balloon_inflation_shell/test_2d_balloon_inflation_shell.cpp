/**
 * @file 	test_2d_fluid_around_balloon_shell.cpp
 * @brief 	Test on fluid-shell interaction when 2 shell particles are close to each other
 * @details This is a case to test fluid-shell interaction.
 * @author 	Weiyi Kong, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

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
    CheckKernelCompleteness(BaseInnerRelation &inner_relation, std::vector<BaseContactRelation *> &contact_relations)
        : particles_(&inner_relation.base_particles_),
          kernel_(inner_relation.getSPHBody().sph_adaptation_->getKernel()),
          inner_configuration_(&inner_relation.inner_configuration_)
    {
        for (size_t n = 0; n < contact_relations.size(); n++)
        {
            auto &contact_relation = *contact_relations[n];
            for (size_t i = 0; i != contact_relation.contact_bodies_.size(); ++i)
            {
                contact_particles_.push_back(&contact_relation.contact_bodies_[i]->getBaseParticles());
                contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
            }
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
const Real scale = 1.0;
const Real radius_balloon = 2.0 * scale;     // initial radius of the balloon
const Real hole_diameter = 0.6245 * scale;   // diameter of the hole
const Real thickness_balloon = 2e-3 * scale; // thickness of the balloon

const Real resolution_ref = hole_diameter / Real(10);
const Real resolution_shell = 0.5 * resolution_ref;

const Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
const Real DL_sponge = resolution_ref * 10; /**< Reference size of the emitter buffer to impose inflow condition. */

const Real radius_balloon_outer = radius_balloon + 0.5 * resolution_shell; // radius of the outer surface
const Real radius_balloon_inner = radius_balloon - 0.5 * resolution_shell; // radius of the outer surface

const Real hole_radius = 0.5 * hole_diameter;
const Vecd balloon_center(sqrt(radius_balloon *radius_balloon - hole_radius * hole_radius), 0);
const BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -radius_balloon_outer),
                                       Vec2d(balloon_center.x() + radius_balloon_outer, radius_balloon_outer));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1.0;  /**< Reference density of air. */
const Real mu_f = 1.5e-5; /**< Dynamics viscosity. */
const Real Q_f = 0.45;
const Real U_f = Q_f / hole_diameter; /**< Characteristic velocity. */
const Real U_max = 1.5 * U_f;
/** Reference sound speed needs to consider the flow speed in the narrow channels. */
const Real c_f = 10.0 * U_max;
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
const Real rho0_s = 100.0; /**< Reference density.*/
const Real youngs_modulus = 1e3;
const Real poisson_ratio = 0.495;
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * youngs_modulus) * thickness_balloon;
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
/** create a Fluid block shape */
std::vector<Vecd> createFluidInletShape()
{
    // geometry
    std::vector<Vecd> Fluid_shape;
    Fluid_shape.push_back(Vecd(-DL_sponge, -hole_radius));
    Fluid_shape.push_back(Vecd(-DL_sponge, hole_radius));
    Fluid_shape.push_back(Vecd(0, hole_radius));
    Fluid_shape.push_back(Vecd(0, -hole_radius));
    Fluid_shape.push_back(Vecd(-DL_sponge, -hole_radius));

    return Fluid_shape;
}
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class FluidBlock : public MultiPolygonShape
{
  public:
    explicit FluidBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createFluidInletShape(), ShapeBooleanOps::add);
        // multi_polygon_.addACircle(balloon_center, radius_balloon_inner, 100, ShapeBooleanOps::add);
    }
};
class ShellParticleGenerator : public SurfaceParticleGenerator
{
  public:
    explicit ShellParticleGenerator(SPHBody &sph_body) : SurfaceParticleGenerator(sph_body){};
    void initializeGeometricVariables() override
    {
        Real alpha_0 = asin(hole_radius / radius_balloon);
        Real arc_length = 2 * (Pi - alpha_0) * radius_balloon;
        int particle_number_mid_surface = int(arc_length / resolution_shell) + 1;
        Real dtheta = 2 * (Pi - alpha_0) / Real(particle_number_mid_surface);
        Real theta = alpha_0;
        for (int i = 0; i < particle_number_mid_surface; i++)
        {
            Vecd pos = radius_balloon * Vecd(-cos(theta), sin(theta)) + balloon_center;
            Vecd normal(cos(theta), -sin(theta));
            initializePositionAndVolumetricMeasure(pos, resolution_shell);
            initializeSurfaceProperties(normal, thickness_balloon);
            theta += dtheta;
        }
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
        : u_ref_(U_f), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_avg = 0.5 * u_ref_ * (1 - cos(2.0 * Pi * run_time / t_ref_));
        return Vecd(u_avg, 0);
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
        if (base_particles_.pos_[index_i][0] <= 0.5 * resolution_shell)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    std::cout << "U_max = " << U_max << std::endl;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidBlock>("FluidBody"));
    fluid_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    fluid_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody shell(sph_system, makeShared<DefaultShape>("Shell"));
    shell.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    shell.generateParticles<ShellParticleGenerator>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation fluid_inner(fluid_block);
    InnerRelation shell_inner(shell);
    ContactRelationToShell fluid_shell_contact(fluid_block, {&shell});
    ContactRelationFromShell shell_fluid_contact(shell, {&fluid_block});
    ComplexRelation fluid_block_complex(fluid_inner, {&fluid_shell_contact});
    // inner relation to compute curvature
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell, fluid_block);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Algorithm for fluid dynamics. */
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization(fluid_block);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(fluid_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(fluid_block);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_fluid_density_by_summation(fluid_inner, fluid_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(fluid_inner, fluid_shell_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> fluid_density_relaxation(fluid_inner, fluid_shell_contact);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(fluid_inner, fluid_shell_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(fluid_inner, fluid_shell_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(fluid_inner, fluid_shell_contact);
    /** Algorithm for in-/outlet. */
    Vec2d emitter_halfsize = Vec2d(0.5 * BW, hole_radius);
    Vec2d emitter_translation = Vec2d(0.5 * BW - DL_sponge, 0.0);
    BodyAlignedBoxByParticle emitter(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec2d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 20, 0);

    Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, hole_radius);
    Vec2d inlet_buffer_translation = Vec2d(-0.5 * DL_sponge, 0.0);
    BodyAlignedBoxByCell emitter_buffer(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec2d(inlet_buffer_translation)), inlet_buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(emitter_buffer);

    /** Algorithm for solid dynamics. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> shell_time_step_size(shell);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> shell_corrected_configuration(shell_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> shell_stress_relaxation_first(shell_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> shell_stress_relaxation_second(shell_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> shell_update_normal(shell);
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
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> shell_constrain(shell_boundary_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_position_damping(0.2, shell_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        shell_rotation_damping(0.2, shell_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    fluid_block.addBodyStateForRecording<Real>("Pressure");
    fluid_block.addBodyStateForRecording<Real>("Density");
    fluid_block.addBodyStateForRecording<int>("Indicator");
    shell.addBodyStateForRecording<Real>("Thickness");
    shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    shell.addBodyStateForRecording<Vecd>("AllForceFromFluid");
    shell.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    shell_corrected_configuration.exec();
    shell_average_curvature.exec();
    fluid_block_complex.updateConfiguration();
    shell_fluid_contact.updateConfiguration();

    //   Check dWijVjeij
    std::vector<BaseContactRelation *> contact_relations = {&fluid_shell_contact};
    CheckKernelCompleteness check_kernel_completeness(fluid_inner, contact_relations);
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
    Real end_time = 16.0;
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
    write_body_states.writeToFile();
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

                    fluid_force_on_shell_update.exec();

                    fluid_density_relaxation.exec(dt);

                    /** Solid dynamics time stepping. */
                    Real dt_s_sum = 0.0;
                    average_velocity_and_acceleration.initialize_displacement_.exec();
                    while (dt_s_sum < dt)
                    {
                        Real dt_s_temp = 0.25 * shell_time_step_size.exec();
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

                        dt_s_sum += dt_s;
                    }
                    average_velocity_and_acceleration.update_averages_.exec(dt);

                    relaxation_time += dt;
                    integration_time += dt;
                    GlobalStaticVariables::physical_time_ += dt;
                }
                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << GlobalStaticVariables::physical_time_
                              << "	Dt = " << Dt << "	dt = " << dt
                              << "  dt_s = " << dt_s
                              << "\n";
                }
                number_of_iterations++;

                /** inflow injection*/
                emitter_inflow_injection.exec();

                /** Update cell linked list and configuration. */
                fluid_block.updateCellLinkedList();

                update_shell_volume();
                shell_update_normal.exec();
                shell.updateCellLinkedList();
                shell_curvature_inner.updateConfiguration();
                shell_average_curvature.exec();

                shell_fluid_contact.updateConfiguration();
                fluid_block_complex.updateConfiguration();
            }
            // exit(0);
            TickCount t2 = TickCount::now();
            check_kernel_completeness.exec();
            write_body_states.writeToFile();
            TickCount t3 = TickCount::now();
            interval += t3 - t2;
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
        write_body_states.writeToFile(1e8);
    }
    return 0;
}
