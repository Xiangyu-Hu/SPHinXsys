/**
 * @file 	test_2d_fluid_around_balloon_container.cpp
 * @brief 	Test on fluid-container interaction when 2 container particles are close to each other
 * @details This is a case to test fluid-container interaction.
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
const bool if_fsi = true;
const Real scale = 1.0;
const Real B = 4.87 * scale; /**< Width of the funnel entrance. */
const Real b = 1.3 * scale;  /**< Width of the funnel bottom. */
const Real h = 2.5 * scale;  /**< Height of the funnel. */
const Real H = 3.75 * scale; /**< Height of the container. */
const Real R = 2.25 * scale; /**< Radius of the container. */
const Real s = 0.2 * scale;  // Thickness of the container

const Real resolution_container = s / 4.0;
const Real resolution_ref = 2 * resolution_container;

const Real BW = resolution_ref * 4; /**< Rigid wall thickness. */

const BoundingBox system_domain_bounds(Vec2d(-B * 0.5 - BW, -R - s), Vec2d(B * 0.5 + BW, H + h + 2 * BW));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1000.0;                     /**< Reference density of fluid. */
const Real gravity_g = 9.81;                    /**< Gravity. */
const Real mu_f = 50;                           /**< Dynamics viscosity. */
const Real bulk_modulus_f = 1.75e7;             /**< Bulk modulus of fluid. */
const Real c_f = sqrt(bulk_modulus_f / rho0_f); /** Reference sound speed */
const Real U_f = c_f / 10.0;                    /**< Characteristic velocity. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
const Real rho0_s = 20.0;          /**< Reference density.*/
const Real youngs_modulus = 2.1e7; // actual: 1e6, ref: https://doi.org/10.5254/1.3547752 eq. 12A
const Real poisson_ratio = 0.3;
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * youngs_modulus) * s;
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_shape;
    water_shape.push_back(Vecd(-0.5 * B, H + h));
    water_shape.push_back(Vecd(0.5 * B, H + h));
    water_shape.push_back(Vecd(0.5 * b, H));
    water_shape.push_back(Vecd(-0.5 * b, H));
    water_shape.push_back(Vecd(-0.5 * B, H + h));

    return water_shape;
}
/** create a wall outer shape */
std::vector<Vecd> createWallLeftShape()
{
    // geometry
    std::vector<Vecd> left_shape;
    left_shape.push_back(Vecd(-0.5 * b, H));
    left_shape.push_back(Vecd(-R - s, H));
    left_shape.push_back(Vecd(-R - s, H + BW));
    left_shape.push_back(Vecd(-0.5 * b - 2 * BW, H + BW));
    left_shape.push_back(Vecd(-0.5 * B - BW, H + h));
    left_shape.push_back(Vecd(-0.5 * B, H + h));
    left_shape.push_back(Vecd(-0.5 * b, H));

    return left_shape;
}
/** create a wall inner shape */
std::vector<Vecd> createWallRightShape()
{
    // geometry
    std::vector<Vecd> right_shape;
    right_shape.push_back(Vecd(0.5 * b, H));
    right_shape.push_back(Vecd(0.5 * B, H + h));
    right_shape.push_back(Vecd(0.5 * B + BW, H + h));
    right_shape.push_back(Vecd(0.5 * b + 2 * BW, H + BW));
    right_shape.push_back(Vecd(R + s, H + BW));
    right_shape.push_back(Vecd(R + s, H));
    right_shape.push_back(Vecd(0.5 * b, H));

    return right_shape;
}
//----------------------------------------------------------------------
//	Define case dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};
class Funnel : public MultiPolygonShape
{
  public:
    explicit Funnel(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWallLeftShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createWallRightShape(), ShapeBooleanOps::add);
    }
};
class Container : public MultiPolygonShape
{
  public:
    explicit Container(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addABox(Transform(Vec2d(0, 0.5 * H)), Vec2d(R + s, 0.5 * H), ShapeBooleanOps::add);
        multi_polygon_.addACircle(Vec2d(0, 0), R + s, 100, ShapeBooleanOps::add);
        multi_polygon_.addABox(Transform(Vec2d(0, 0.5 * H)), Vec2d(R, 0.5 * H), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(Vec2d(0, 0), R, 100, ShapeBooleanOps::sub);
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
        if (std::abs(base_particles_.pos_[index_i][1]) > H - BW)
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
    std::cout << "U_f: " << U_f << std::endl;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for running particle relaxation for the initially body-fitted distribution */
    sph_system.setRunParticleRelaxation(false);
    /** Tag for starting with relaxed body-fitted particles distribution */
    sph_system.setReloadParticles(false);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody funnel(sph_system, makeShared<Funnel>("Funnel"));
    funnel.defineParticlesAndMaterial<SolidParticles, Solid>();
    funnel.generateParticles<ParticleGeneratorLattice>();

    SolidBody container(sph_system, makeShared<Container>("Container"));
    container.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_container);
    container.defineBodyLevelSetShape();
    container.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    if (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        container.generateParticles<ParticleGeneratorReload>(io_environment, container.getName());
    else
        container.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Define body relation map used for particle relaxation.
        //----------------------------------------------------------------------
        InnerRelation container_inner(container);
        //----------------------------------------------------------------------
        //	Define the methods for particle relaxation for wall boundary.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> container_random_particles(container);
        relax_dynamics::RelaxationStepInner relaxation_step_container_inner(container_inner);
        //----------------------------------------------------------------------
        //	Output for particle relaxation.
        //----------------------------------------------------------------------
        BodyStatesRecordingToVtp write_relaxed_particles(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_particle_reload_files(io_environment, {&container});
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        container_random_particles.exec(0.25);
        relaxation_step_container_inner.SurfaceBounding().exec();
        write_relaxed_particles.writeToFile(0);
        //----------------------------------------------------------------------
        //	From here iteration for particle relaxation begins.
        //----------------------------------------------------------------------
        int ite = 0;
        int relax_step = 1000;
        while (ite < relax_step)
        {
            relaxation_step_container_inner.exec();
            ite += 1;
            if (ite % 100 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite << "\n";
                write_relaxed_particles.writeToFile(ite);
            }
        }
        std::cout << "The physics relaxation process of ball particles finish !" << std::endl;
        write_particle_reload_files.writeToFile(0);
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
    InnerRelation water_inner(water_block);
    InnerRelation container_inner(container);
    ContactRelation water_solid_contact(water_block, {&funnel, &container});
    ContactRelation container_water_contact(container, {&water_block});
    ComplexRelation water_block_complex(water_inner, {&water_solid_contact});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Algorithm for fluid dynamics. */
    SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization(water_block, gravity_ptr);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_fluid_density_by_summation(water_inner, water_solid_contact);
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> fluid_pressure_relaxation(water_inner, water_solid_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> fluid_density_relaxation(water_inner, water_solid_contact);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_inner, water_solid_contact);
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex> inlet_outlet_surface_particle_indicator(water_inner, water_solid_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_velocity_correction(water_inner, water_solid_contact);
    /** Algorithm for solid dynamics. */
    SimpleDynamics<NormalDirectionFromBodyShape> funnel_normal_direction(funnel);
    SimpleDynamics<NormalDirectionFromBodyShape> container_normal_direction(container);
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> container_time_step_size(container);
    InteractionWithUpdate<KernelCorrectionMatrixInner> container_corrected_configuration(container_inner);
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> container_stress_relaxation_first(container_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> container_stress_relaxation_second(container_inner);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> container_update_normal(container);
    // auto update_container_volume = [&]()
    // {
    //     particle_for(
    //         par,
    //         container.getBaseParticles().total_real_particles_,
    //         [&](size_t index_i)
    //         {
    //             container.getBaseParticles().Vol_[index_i] = container.getBaseParticles().mass_[index_i] / container.getBaseParticles().rho_[index_i];
    //         });
    // };
    auto update_average_velocity = [&]()
    {
        auto avg_vel = *container.getBaseParticles().getVariableByName<Vecd>("AverageVelocity");
        auto avg_force = *container.getBaseParticles().getVariableByName<Vecd>("AverageForce");
        particle_for(
            par,
            container.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                avg_vel[index_i] = container.getBaseParticles().vel_[index_i];
                avg_force[index_i] = container.getBaseParticles().force_[index_i];
            });
    };
    /** FSI */
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_container(container_water_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_container_update(container_water_contact, viscous_force_on_container);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(container);
    /** constraint and damping */
    BoundaryGeometry container_boundary_geometry(container, "BoundaryGeometry");
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> container_constrain(container_boundary_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        container_position_damping(0.2, container_inner, "Velocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<Real>("Density");
    container.addBodyStateForRecording<Vecd>("AllForceFromFluid");
    container.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    funnel_normal_direction.exec();
    container_normal_direction.exec();
    container_corrected_configuration.exec();

    //   Check dWijVjeij
    std::vector<BaseContactRelation *> contact_relations = {&water_solid_contact};
    CheckKernelCompleteness check_kernel_completeness(water_inner, contact_relations);
    check_kernel_completeness.exec();
    water_block.addBodyStateForRecording<Real>("TotalKernel");
    water_block.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    water_block.addBodyStateForRecording<int>("InnerNeighborNumber");
    water_block.addBodyStateForRecording<int>("ContactNeighborNumber");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 10;
    Real end_time = 10;
    Real output_interval = end_time / 200.0; /**< Time stamps for output of body states. */
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
    // const double Dt_ref = fluid_advection_time_step.exec();
    const double dt_ref = fluid_acoustic_time_step.exec();
    const double dt_s_ref = container_time_step_size.exec();
    // const double dt_stable = 4e-6;
    // const double dt_s_stable = 8.6e-6;
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
                /** Dynamics including pressure relaxation. */
                Real dt_temp = fluid_acoustic_time_step.exec();
                if (dt_temp < dt_ref / 20)
                {
                    std::cout << "dt = " << dt_temp << ", dt_ref = " << dt_ref << std::endl;
                    std::cout << "Acoustic time step decreased too much!" << std::endl;
                    throw std::runtime_error("Acoustic time step decreased too much!");
                }
                Real dt_s_temp = container_time_step_size.exec();
                if (dt_s_temp < dt_s_ref / 100)
                {
                    std::cout << "dt_s = " << dt_s_temp << ", dt_s_ref = " << dt_s_ref << std::endl;
                    std::cout << "container time step decreased too much!" << std::endl;
                    throw std::runtime_error("container time step decreased too much!");
                }
                Real dt = std::min(dt_s_temp, dt_temp);
                // inlet_outlet_surface_particle_indicator.exec();
                update_fluid_density_by_summation.exec();
                viscous_acceleration.exec();
                // transport_velocity_correction.exec();
                if (if_fsi)
                    viscous_force_on_container.exec();
                container_update_normal.exec();
                fluid_pressure_relaxation.exec(dt);
                if (if_fsi)
                    fluid_force_on_container_update.exec();
                fluid_density_relaxation.exec(dt);
                /** Solid dynamics time stepping. */
                if (if_fsi)
                {
                    // average_velocity_and_acceleration.initialize_displacement_.exec();
                    container_stress_relaxation_first.exec(dt);
                    container_constrain.exec();
                    container_position_damping.exec(dt);
                    container_constrain.exec();
                    container_stress_relaxation_second.exec(dt);
                    // average_velocity_and_acceleration.update_averages_.exec(dt);
                    update_average_velocity();
                }
                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << GlobalStaticVariables::physical_time_
                              << "	dt = " << dt << "\n";
                }
                number_of_iterations++;

                /** Update cell linked list and configuration. */
                water_block.updateCellLinkedList();
                if (if_fsi)
                {
                    // update_container_volume();
                    container.updateCellLinkedList();
                    container_water_contact.updateConfiguration();
                }
                water_block_complex.updateConfiguration();

                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
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
        water_block.setNewlyUpdated();
        container.setNewlyUpdated();
        write_body_states.writeToFile(1e8);
    }
    return 0;
}
