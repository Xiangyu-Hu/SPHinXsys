/**
 * @file 	test_2d_fluid_filling_elastic_container_shell.cpp
 * @brief 	Test on fluid filling into an elastic shell container
 * @details This is a case to test fluid-shell interaction.
 * @author 	Weiyi Kong, Virtonomy GmbH
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

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
const Real rho0_f = 1000.0;                           /**< Reference density of fluid. */
const Real gravity_g = 9.81;                          /**< Gravity. */
const Real mu_f = 50;                                 /**< Dynamics viscosity. */
const Real bulk_modulus_f = 1.75e7;                   /**< Bulk modulus of fluid. */
const Real c_f = sqrt(bulk_modulus_f / rho0_f);       /** Reference sound speed */
const Real U_f = 2.0 * sqrt(gravity_g * (H + h + R)); /**< Characteristic velocity. */
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
class ContainerParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit ContainerParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        Real mid_surface_radius = R + 0.5 * resolution_container;
        int particle_circle = int(Pi * mid_surface_radius / resolution_container) + 1;
        Real dtheta = Pi / Real(particle_circle - 1);
        for (int i = 0; i < particle_circle; i++)
        {
            Real theta = dtheta * i;
            Real x = -mid_surface_radius * cos(theta);
            Real y = -mid_surface_radius * sin(theta);
            Vecd normal = Vecd(x, y).normalized();
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_container);
            initializeSurfaceProperties(normal, s);
        }
        Real y = resolution_container;
        while (y < H)
        {
            initializePositionAndVolumetricMeasure(Vecd(mid_surface_radius, y), resolution_container);
            initializeSurfaceProperties(Vecd(1, 0), s);
            initializePositionAndVolumetricMeasure(Vecd(-mid_surface_radius, y), resolution_container);
            initializeSurfaceProperties(Vecd(-1, 0), s);
            y += resolution_container;
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

    SolidBody container(sph_system, makeShared<DefaultShape>("Container"));
    container.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_container);
    container.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    container.generateParticles<ContainerParticleGenerator>();

    ObserverBody disp_observer(sph_system, "Observer");
    StdVec<Vecd> observation_location = {Vec2d(0, -R - 0.5 * resolution_container)};
    disp_observer.generateParticles<ParticleGeneratorObserver>(observation_location);
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
    ContactRelation water_funnel_contact(water_block, {&funnel});
    ContactRelationToShell water_container_contact(water_block, {&container}, false);
    ContactRelationFromShell container_water_contact(container, {&water_block}, false);
    ComplexRelation water_block_complex(water_inner, {&water_funnel_contact, &water_container_contact});
    ShellInnerRelationWithContactKernel container_curvature_inner(container, water_block);
    ContactRelation disp_observer_contact(disp_observer, {&container});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Algorithm for fluid dynamics. */
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(water_block, gravity);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeSurface>, Contact<>, Contact<>>> update_fluid_density_by_summation(water_inner, water_funnel_contact, water_container_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> fluid_pressure_relaxation(water_inner, water_funnel_contact, water_container_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, NoRiemannSolver>> fluid_density_relaxation(water_inner, water_funnel_contact, water_container_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::ViscousForce<Inner<>, Contact<Wall>, Contact<Wall>>, fluid_dynamics::FixedViscosity>> viscous_acceleration(water_inner, water_funnel_contact, water_container_contact);
    /** Algorithm for solid dynamics. */
    SimpleDynamics<NormalDirectionFromBodyShape> funnel_normal_direction(funnel);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> container_time_step_size(container);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> container_corrected_configuration(container_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> container_stress_relaxation_first(container_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> container_stress_relaxation_second(container_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> container_average_curvature(container_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> container_update_normal(container);
    auto update_container_volume = [&]()
    {
        particle_for(
            par,
            container.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                container.getBaseParticles().Vol_[index_i] = container.getBaseParticles().mass_[index_i] / container.getBaseParticles().rho_[index_i] / s;
            });
    };
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
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_container(container_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid> pressure_force_on_container(container_water_contact);
    /** constraint and damping */
    BoundaryGeometry container_boundary_geometry(container, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> container_constrain(container_boundary_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        container_position_damping(0.2, container_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        container_rotation_damping(0.2, container_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<int>("Indicator");
    water_block.addBodyStateForRecording<Real>("Pressure");
    water_block.addBodyStateForRecording<Real>("Density");
    container.addBodyStateForRecording<Vecd>("PressureForceFromFluid");
    container.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    container.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    BodyStatesRecordingToVtp write_body_states(sph_system.real_bodies_);
    ObservedQuantityRecording<Vecd> write_tip_displacement("Displacement", disp_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    funnel_normal_direction.exec();
    container_corrected_configuration.exec();
    container_average_curvature.exec();
    water_block_complex.updateConfiguration();
    container_water_contact.updateConfiguration();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 10;
    Real end_time = 10.0;
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
    write_tip_displacement.writeToFile(number_of_iterations);
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    const double Dt_ref = fluid_advection_time_step.exec();
    const double dt_ref = fluid_acoustic_time_step.exec();
    const double dt_s_ref = container_time_step_size.exec();

    auto run_simulation = [&]()
    {
        std::cout << "Simulation starts here" << std::endl;
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < output_interval)
            {
                /** Dynamics including pressure relaxation. */
                Real Dt = fluid_advection_time_step.exec();
                if (Dt < Dt_ref / 20)
                {
                    std::cout << "Dt = " << Dt << ", Dt_ref = " << Dt_ref << std::endl;
                    std::cout << "Advective time step decreased too much!" << std::endl;
                    throw std::runtime_error("Advective time step decreased too much!");
                }
                Real dt_temp = fluid_acoustic_time_step.exec();
                if (dt_temp < dt_ref / 20)
                {
                    std::cout << "dt = " << dt_temp << ", dt_ref = " << dt_ref << std::endl;
                    std::cout << "Acoustic time step decreased too much!" << std::endl;
                    throw std::runtime_error("Acoustic time step decreased too much!");
                }
                Real dt_s_temp = container_time_step_size.exec();
                if (dt_s_temp < dt_s_ref / 20)
                {
                    std::cout << "dt_s = " << dt_s_temp << ", dt_s_ref = " << dt_s_ref << std::endl;
                    std::cout << "container time step decreased too much!" << std::endl;
                    throw std::runtime_error("container time step decreased too much!");
                }
                Real dt = std::min({dt_s_temp, dt_temp, Dt});
                update_fluid_density_by_summation.exec();
                viscous_acceleration.exec();
                if (if_fsi)
                    viscous_force_on_container.exec();

                fluid_pressure_relaxation.exec(dt);
                if (if_fsi)
                    pressure_force_on_container.exec();
                fluid_density_relaxation.exec(dt);
                /** Solid dynamics time stepping. */
                if (if_fsi)
                {
                    container_stress_relaxation_first.exec(dt);
                    container_constrain.exec();
                    container_position_damping.exec(dt);
                    container_rotation_damping.exec(dt);
                    container_constrain.exec();
                    container_stress_relaxation_second.exec(dt);
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
                    update_container_volume();
                    container_update_normal.exec();
                    container.updateCellLinkedList();
                    container_curvature_inner.updateConfiguration();
                    container_average_curvature.exec();
                    container_water_contact.updateConfiguration();
                }
                water_block_complex.updateConfiguration();

                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            TickCount t2 = TickCount::now();
            write_body_states.writeToFile();
            write_tip_displacement.writeToFile(number_of_iterations);
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
