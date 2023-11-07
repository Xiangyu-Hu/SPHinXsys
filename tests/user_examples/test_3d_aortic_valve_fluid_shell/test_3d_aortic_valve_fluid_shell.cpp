#include "sphinxsys.h"
using namespace SPH;

class CheckKernelCompleteness
{
  private:
    BaseParticles *particles_;
    std::vector<SPH::BaseParticles *> contact_particles_;
    ParticleConfiguration *inner_configuration_;
    std::vector<ParticleConfiguration *> contact_configuration_;

    StdLargeVec<Real> W_ijV_j_ttl;
    StdLargeVec<Vecd> dW_ijV_je_ij_ttl;
    StdLargeVec<int> number_of_inner_neighbor;
    StdLargeVec<int> number_of_contact_neighbor;

  public:
    CheckKernelCompleteness(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : particles_(&inner_relation.base_particles_), inner_configuration_(&inner_relation.inner_configuration_)
    {
        for (size_t i = 0; i != contact_relation.contact_bodies_.size(); ++i)
        {
            contact_particles_.push_back(&contact_relation.contact_bodies_[i]->getBaseParticles());
            contact_configuration_.push_back(&contact_relation.contact_configuration_[i]);
        }
        inner_relation.base_particles_.registerVariable(W_ijV_j_ttl, "TotalKernel");
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
                Real W_ijV_j_ttl_i = 0;
                Vecd dW_ijV_je_ij_ttl_i = Vecd::Zero();
                const Neighborhood &inner_neighborhood = (*inner_configuration_)[index_i];
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    W_ijV_j_ttl_i += inner_neighborhood.W_ij_[n] * particles_->Vol_[index_j];
                    dW_ijV_je_ij_ttl_i += inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
                    N_inner_number++;
                }

                for (size_t k = 0; k < contact_configuration_.size(); ++k)
                {
                    const SPH::Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
                    for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
                    {
                        size_t index_j = wall_neighborhood.j_[n];
                        W_ijV_j_ttl_i += wall_neighborhood.W_ij_[n] * contact_particles_[k]->Vol_[index_j];
                        dW_ijV_je_ij_ttl_i += wall_neighborhood.dW_ijV_j_[n] * wall_neighborhood.e_ij_[n];
                        N_contact_number++;
                    }
                }
                W_ijV_j_ttl[index_i] = W_ijV_j_ttl_i;
                dW_ijV_je_ij_ttl[index_i] = dW_ijV_je_ij_ttl_i;
                number_of_inner_neighbor[index_i] = N_inner_number;
                number_of_contact_neighbor[index_i] = N_contact_number;
            });
    }
};

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real scale = 0.001;
const Real diameter = 25 * scale; /**<Inlet and outlet diameter. */
const Real radius = 0.5 * diameter;
const Real upstream_length = 1.5 * diameter;   /**<Length before valve. */
const Real downstream_length = 3.0 * diameter; /**<Length after valve. */
const Real wall_thickness = 1 * scale;
const Real valve_thickness = 0.5 * scale; /**<Thickness of valve. */
const Real resolution_fluid = 1 * scale;  /**< Global reference resolution. */
const Real resolution_wall = wall_thickness;
const Real resolution_valve = valve_thickness;
const Real wall_extension = resolution_fluid * 4.0; /**<Extension of wall. */
const Real inflow_length = 10 * resolution_fluid;   /**< Inflow region. */

/** Domain bounds of the system. */
const BoundingBox system_domain_bounds(Vec3d(-upstream_length - wall_extension, -radius - wall_thickness, -radius - wall_thickness),
                                       Vec3d(downstream_length + wall_extension, radius + wall_thickness, radius + wall_thickness));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1100.0;                      // blood density
const Real mu_f = 3.6e-3;                        // blood viscosity
const Real Q_max = 11 * 1e-3 / 60.0;             // 1L/min
const Real U_f = Q_max / (Pi * radius * radius); // average velocity at peak systole
const Real U_max = 6 * U_f;
const Real Re = rho0_f * U_f * diameter / mu_f; /**< Reynolds number. */
const Real c_f = 10.0 * U_max;                  /**< Speed of sound. */
//----------------------------------------------------------------------
//	Global parameters on the solid properties
//----------------------------------------------------------------------
const double rho0_s = 1000.0;
const double Youngs_modulus = 3.6e6 * std::pow(0.15 * scale / valve_thickness, 3); // use softer material for development
const double poisson = 0.45;
const Real shape_constant = 0.4;
const Real physical_viscosity = 20 * shape_constant / 4.0 * std::sqrt(rho0_s * Youngs_modulus) * valve_thickness;
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
const std::string path_to_fluid_file = "./input/fluid_block.stl";
const std::string path_to_leaflet_1_file = "./input/leaflet_1.stl";
const std::string path_to_leaflet_2_file = "./input/leaflet_2.stl";
const std::string path_to_leaflet_3_file = "./input/leaflet_3.stl";

class FluidShape : public ComplexShape
{
  public:
    explicit FluidShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_fluid_file, Vecd::Zero(), 1.0);
        subtract<TriangleMeshShapeSTL>(path_to_leaflet_1_file, Vecd::Zero(), 1.0);
        subtract<TriangleMeshShapeSTL>(path_to_leaflet_2_file, Vecd::Zero(), 1.0);
        subtract<TriangleMeshShapeSTL>(path_to_leaflet_3_file, Vecd::Zero(), 1.0);
    }
};
const Vec3d buffer_halfsize(0.5 * inflow_length, radius, radius);
const Vec3d buffer_translation(-upstream_length + 0.5 * inflow_length, 0.0, 0.0);

/** create the wall constrain shape. */
const std::string path_to_base_file = "./input/base.stl";
class FixedShape : public ComplexShape
{
  public:
    explicit FixedShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_base_file, Vecd::Zero(), 1.0);
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
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = 0.5 * u_ref_ * (1.0 - cos(2 * Pi * run_time / t_ref_));

        Real radius2 = position[1] * position[1] + position[2] * position[2];
        target_velocity[0] = 2.0 * u_ave * SMAX(0.0, 1.0 - radius2 / radius / radius);

        return target_velocity;
    }
};

int main(int ac, char *av[])
{
    std::cout << "U_f = " << U_f << std::endl;
    std::cout << "Inlet Re = " << Re << std::endl;
    std::cout << "Young's modulus after scaling = " << Youngs_modulus << std::endl;
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_fluid);
#ifdef BOOST_AVAILABLE
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
#endif
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidShape>("fluid"));
    fluid_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    fluid_block.generateParticles<ParticleGeneratorLattice>();
    std::cout << "Fluid object generation finished..." << std::endl;

    SolidBody wall_boundary(sph_system, makeShared<DefaultShape>("wall_boundary"));
    wall_boundary.defineAdaptation<SPH::SPHAdaptation>(1.15, resolution_fluid / resolution_wall);
    wall_boundary.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1.0, 1.0, 0.0); // dummy material parameters
    wall_boundary.generateParticles<ParticleGeneratorReload>(io_environment, wall_boundary.getName());
    std::cout << "Wall object generation finished..." << std::endl;

    SolidBody leaflet_1(sph_system, makeShared<DefaultShape>("leaflet_1"));
    leaflet_1.defineAdaptationRatios(1.15, resolution_fluid / resolution_valve);
    leaflet_1.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    leaflet_1.generateParticles<ParticleGeneratorReload>(io_environment, leaflet_1.getName());
    std::cout << "leaflet_1 object generation finished..." << std::endl;

    SolidBody leaflet_2(sph_system, makeShared<DefaultShape>("leaflet_2"));
    leaflet_2.defineAdaptationRatios(1.15, resolution_fluid / resolution_valve);
    leaflet_2.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    leaflet_2.generateParticles<ParticleGeneratorReload>(io_environment, leaflet_2.getName());
    std::cout << "leaflet_2 object generation finished..." << std::endl;

    SolidBody leaflet_3(sph_system, makeShared<DefaultShape>("leaflet_3"));
    leaflet_3.defineAdaptationRatios(1.15, resolution_fluid / resolution_valve);
    leaflet_3.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    leaflet_3.generateParticles<ParticleGeneratorReload>(io_environment, leaflet_3.getName());
    std::cout << "leaflet_3 object generation finished..." << std::endl;
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation wall_boundary_inner(wall_boundary);
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> wall_boundary_curvature(wall_boundary_inner);

    InnerRelation leaflet_1_inner(leaflet_1);
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> leaflet_1_curvature(leaflet_1_inner);

    InnerRelation leaflet_2_inner(leaflet_2);
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> leaflet_2_curvature(leaflet_2_inner);

    InnerRelation leaflet_3_inner(leaflet_3);
    SimpleDynamics<thin_structure_dynamics::ShellCurvature> leaflet_3_curvature(leaflet_3_inner);

    ComplexRelation fluid_block_complex(fluid_block, RealBodyVector{&wall_boundary, &leaflet_1, &leaflet_2, &leaflet_3});

    ContactRelation leaflet_1_fluid_contact(leaflet_1, {&fluid_block});
    ContactRelation leaflet_2_fluid_contact(leaflet_2, {&fluid_block});
    ContactRelation leaflet_3_fluid_contact(leaflet_3, {&fluid_block});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Initialize particle acceleration. */
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(fluid_block);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_block_complex);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(fluid_block, U_max);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(fluid_block);
    /** Pressure relaxation using verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(fluid_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(fluid_block_complex);
    /** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_correction(fluid_block_complex);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(fluid_block_complex);
    /** Inflow boundary condition. */
    BodyAlignedBoxByCell inflow_buffer(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec3d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
    /** Periodic BCs in x direction. */
    PeriodicConditionUsingCellLinkedList periodic_condition(fluid_block, fluid_block.getBodyShapeBounds(), xAxis);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    /** Corrected configuration for shell body. */
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> wall_boundary_corrected_configuration(wall_boundary_inner);
    // leaflet 1
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> leaflet_1_corrected_configuration(leaflet_1_inner);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::FluidViscousForceOnShell> viscous_force_on_leaflet_1(leaflet_1_fluid_contact);
    InteractionDynamics<solid_dynamics::FluidForceOnShellUpdateRiemann>
        fluid_force_on_leaflet_1_update(leaflet_1_fluid_contact, viscous_force_on_leaflet_1);
    /** Compute the average velocity of the insert body. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration_leaflet_1(leaflet_1);
    // leaflet 2
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> leaflet_2_corrected_configuration(leaflet_2_inner);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::FluidViscousForceOnShell> viscous_force_on_leaflet_2(leaflet_2_fluid_contact);
    InteractionDynamics<solid_dynamics::FluidForceOnShellUpdateRiemann>
        fluid_force_on_leaflet_2_update(leaflet_2_fluid_contact, viscous_force_on_leaflet_2);
    /** Compute the average velocity of the insert body. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration_leaflet_2(leaflet_2);
    // leaflet 1
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> leaflet_3_corrected_configuration(leaflet_3_inner);
    /** Compute the force exerted on solid body due to fluid pressure and viscosity. */
    InteractionDynamics<solid_dynamics::FluidViscousForceOnShell> viscous_force_on_leaflet_3(leaflet_3_fluid_contact);
    InteractionDynamics<solid_dynamics::FluidForceOnShellUpdateRiemann>
        fluid_force_on_leaflet_3_update(leaflet_3_fluid_contact, viscous_force_on_leaflet_3);
    /** Compute the average velocity of the insert body. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration_leaflet_3(leaflet_3);
    //----------------------------------------------------------------------
    //	Algorithms of solid dynamics.
    //----------------------------------------------------------------------
    /** Compute time step size of elastic solid. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> leaflet_1_computing_time_step_size(leaflet_1);
    /** Stress relaxation for the inserted body. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> leaflet_1_stress_relaxation_first_half(leaflet_1_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> leaflet_1_stress_relaxation_second_half(leaflet_1_inner);
    /** Constrain region of the inserted body. */
    BodyRegionByParticle leaflet_1_base(leaflet_1, makeShared<FixedShape>("valve_fixation_shape"));
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constraint_leaflet_1_base(leaflet_1_base);
    /** damping */
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        leaflet_1_position_damping(0.2, leaflet_1_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        leaflet_1_rotation_damping(0.2, leaflet_1_inner, "AngularVelocity", physical_viscosity);

    /** Compute time step size of elastic solid. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> leaflet_2_computing_time_step_size(leaflet_2);
    /** Stress relaxation for the inserted body. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> leaflet_2_stress_relaxation_first_half(leaflet_2_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> leaflet_2_stress_relaxation_second_half(leaflet_2_inner);
    /** Constrain region of the inserted body. */
    BodyRegionByParticle leaflet_2_base(leaflet_2, makeShared<FixedShape>("valve_fixation_shape"));
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constraint_leaflet_2_base(leaflet_2_base);
    /** damping */
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        leaflet_2_position_damping(0.2, leaflet_2_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        leaflet_2_rotation_damping(0.2, leaflet_2_inner, "AngularVelocity", physical_viscosity);

    /** Compute time step size of elastic solid. */
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> leaflet_3_computing_time_step_size(leaflet_3);
    /** Stress relaxation for the inserted body. */
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> leaflet_3_stress_relaxation_first_half(leaflet_3_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> leaflet_3_stress_relaxation_second_half(leaflet_3_inner);
    /** Constrain region of the inserted body. */
    BodyRegionByParticle leaflet_3_base(leaflet_3, makeShared<FixedShape>("valve_fixation_shape"));
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> constraint_leaflet_3_base(leaflet_3_base);
    /** damping */
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        leaflet_3_position_damping(0.2, leaflet_3_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingBySplittingInner<Vecd>>>
        leaflet_3_rotation_damping(0.2, leaflet_3_inner, "AngularVelocity", physical_viscosity);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    fluid_block.addBodyStateForRecording<Real>("Pressure");
    fluid_block.addBodyStateForRecording<Real>("Density");
    fluid_block.addBodyStateForRecording<Real>("VolumetricMeasure");
    fluid_block.addBodyStateForRecording<Real>("MassiveMeasure");
    wall_boundary.addBodyStateForRecording<Real>("MeanCurvature");
    leaflet_1.addBodyStateForRecording<Real>("MeanCurvature");
    leaflet_1.addBodyStateForRecording<Vecd>("ForceFromFluid");
    leaflet_1.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    leaflet_2.addBodyStateForRecording<Real>("MeanCurvature");
    leaflet_2.addBodyStateForRecording<Vecd>("ForceFromFluid");
    leaflet_2.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    leaflet_3.addBodyStateForRecording<Real>("MeanCurvature");
    leaflet_3.addBodyStateForRecording<Vecd>("ForceFromFluid");
    leaflet_3.addBodyStateForRecording<Vecd>("ViscousForceFromFluid");
    BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing linear reproducing configuration for the insert body. */
    wall_boundary_corrected_configuration.exec();
    wall_boundary_curvature.compute_initial_curvature();
    leaflet_1_corrected_configuration.exec();
    leaflet_1_curvature.compute_initial_curvature();
    leaflet_2_corrected_configuration.exec();
    leaflet_2_curvature.compute_initial_curvature();
    leaflet_3_corrected_configuration.exec();
    leaflet_3_curvature.compute_initial_curvature();
    fluid_block_complex.getContactRelation().updateConfiguration();

    //  Check dWijVjeij
    CheckKernelCompleteness check_kernel_completeness(fluid_block_complex.getInnerRelation(), fluid_block_complex.getContactRelation());
    check_kernel_completeness.exec();
    fluid_block.addBodyStateForRecording<Real>("TotalKernel");
    fluid_block.addBodyStateForRecording<Vecd>("TotalKernelGrad");
    fluid_block.addBodyStateForRecording<int>("InnerNeighborNumber");
    fluid_block.addBodyStateForRecording<int>("ContactNeighborNumber");
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 10;
    Real end_time = 1.0; // 1 cycle
    Real output_interval = end_time / 100.0;
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_correction.exec();

            /** FSI for viscous force. */
            viscous_force_on_leaflet_1.exec();
            viscous_force_on_leaflet_2.exec();
            viscous_force_on_leaflet_3.exec();
            /** Update normal direction on elastic body.*/
            size_t inner_ite_dt = 0;
            size_t inner_ite_dt_s = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** FSI for pressure force. */
                fluid_force_on_leaflet_1_update.exec();
                fluid_force_on_leaflet_2_update.exec();
                fluid_force_on_leaflet_3_update.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                /** Solid dynamics. */
                inner_ite_dt_s = 0;
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration_leaflet_1.initialize_displacement_.exec();
                average_velocity_and_acceleration_leaflet_2.initialize_displacement_.exec();
                average_velocity_and_acceleration_leaflet_3.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    Real dt_s = std::min({leaflet_1_computing_time_step_size.exec(),
                                          leaflet_2_computing_time_step_size.exec(),
                                          leaflet_3_computing_time_step_size.exec(),
                                          dt - dt_s_sum});
                    leaflet_1_stress_relaxation_first_half.exec(dt_s);
                    constraint_leaflet_1_base.exec();
                    leaflet_1_position_damping.exec(dt_s);
                    leaflet_1_rotation_damping.exec(dt_s);
                    constraint_leaflet_1_base.exec();
                    leaflet_1_stress_relaxation_second_half.exec(dt_s);

                    leaflet_2_stress_relaxation_first_half.exec(dt_s);
                    constraint_leaflet_2_base.exec();
                    leaflet_2_position_damping.exec(dt_s);
                    leaflet_2_rotation_damping.exec(dt_s);
                    constraint_leaflet_2_base.exec();
                    leaflet_2_stress_relaxation_second_half.exec(dt_s);

                    leaflet_3_stress_relaxation_first_half.exec(dt_s);
                    constraint_leaflet_3_base.exec();
                    leaflet_3_position_damping.exec(dt_s);
                    leaflet_3_rotation_damping.exec(dt_s);
                    constraint_leaflet_3_base.exec();
                    leaflet_3_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    inner_ite_dt_s++;
                }
                average_velocity_and_acceleration_leaflet_1.update_averages_.exec(dt);
                average_velocity_and_acceleration_leaflet_2.update_averages_.exec(dt);
                average_velocity_and_acceleration_leaflet_3.update_averages_.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
                parabolic_inflow.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "	dt / dt_s = " << inner_ite_dt_s << "\n";
            }
            number_of_iterations++;

            leaflet_1_curvature.exec();
            leaflet_2_curvature.exec();
            leaflet_3_curvature.exec();

            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();

            fluid_block.updateCellLinkedListWithParticleSort(100);
            leaflet_1.updateCellLinkedList();
            leaflet_2.updateCellLinkedList();
            leaflet_3.updateCellLinkedList();

            periodic_condition.update_cell_linked_list_.exec();

            /** one need update configuration after periodic condition. */
            fluid_block_complex.updateConfiguration();
            leaflet_1_fluid_contact.updateConfiguration();
            leaflet_2_fluid_contact.updateConfiguration();
            leaflet_3_fluid_contact.updateConfiguration();
        }

        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
