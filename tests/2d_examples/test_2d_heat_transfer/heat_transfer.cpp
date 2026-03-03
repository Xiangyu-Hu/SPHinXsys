/**
 * @file 	heat_transfer.cpp
 * @brief 	This is a test to validate heat transfer between a flow and channel walls.
 * @author 	Xiaojing Tang, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 2.0;                          /**< Channel length. */
Real DH = 0.4;                          /**< Channel height. */
Real global_resolution = DH / 25.0;        /**< Global reference resolution. */
Real DL_sponge = global_resolution * 20.0; /**< Sponge region to impose inflow condition. */
/** Boundary width, determined by specific layer of boundary particles. */
Real BW = global_resolution * 4.0; /** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
// temperature observer location
StdVec<Vecd> observation_location = {Vecd(0.0, DH * 0.5)};
//----------------------------------------------------------------------
//	Global parameters on the material properties
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 1.0e-3;
Real rho0_f = 1.0;                  /**< Density. */
Real U_f = 1.0;                     /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;              /**< Speed of sound. */
Real Re = 100.0;                    /**< Reynolds number100. */
Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Global parameters on the initial condition
//----------------------------------------------------------------------
Real phi_upper_wall = 20.0;
Real phi_lower_wall = 40.0;
Real phi_fluid_initial = 20.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
/** create a water block shape */
std::vector<Vecd> createShape()
{
    // geometry
    std::vector<Vecd> shape;
    shape.push_back(Vecd(0.0 - DL_sponge, 0.0));
    shape.push_back(Vecd(0.0 - DL_sponge, DH));
    shape.push_back(Vecd(DL, DH));
    shape.push_back(Vecd(DL, 0.0));
    shape.push_back(Vecd(0.0 - DL_sponge, 0.0));
    return shape;
}
/** create outer wall shape */
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
    outer_wall_shape.push_back(Vecd(-DL_sponge - BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
    outer_wall_shape.push_back(Vecd(DL + BW, -BW));
    outer_wall_shape.push_back(Vecd(-DL_sponge - BW, -BW));
    return outer_wall_shape;
}
/** create inner wall shape */
std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
    inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
    inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, DH));
    inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
    inner_wall_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

    return inner_wall_shape;
}
Vec2d buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH);
Vec2d buffer_translation = Vec2d(-DL_sponge, 0.0) + buffer_halfsize;
//----------------------------------------------------------------------
//	Case-dependent geometries
//----------------------------------------------------------------------
class ThermofluidBody : public MultiPolygonShape
{
  public:
    explicit ThermofluidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createShape(), ShapeBooleanOps::add);
    }
};
class ThermosolidBody : public MultiPolygonShape
{
  public:
    explicit ThermosolidBody(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Application dependent solid body initial condition
//----------------------------------------------------------------------
class ThermosolidBodyInitialCondition : public LocalDynamics
{
  public:
    explicit ThermosolidBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        if (-BW <= pos_[index_i][1] && pos_[index_i][1] <= 0.0)
        {
            phi_[index_i] = phi_lower_wall;
        }

        if (DH <= pos_[index_i][1] && pos_[index_i][1] <= DH + BW)
        {
            phi_[index_i] = phi_upper_wall;
        }
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Application dependent fluid body initial condition
//----------------------------------------------------------------------
class ThermofluidBodyInitialCondition : public LocalDynamics
{
  public:
    explicit ThermofluidBodyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        if (0 <= pos_[index_i][1] && pos_[index_i][1] <= DH)
        {
            phi_[index_i] = phi_fluid_initial;
        }
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
using ThermalRelaxationComplex = DiffusionBodyRelaxationComplex<
    IsotropicDiffusion, KernelGradientInner, KernelGradientContact, Dirichlet>;
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBox &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f), t_ref_(2.0),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity, Real current_time)
    {
        Vecd target_velocity = velocity;
        Real u_ave = current_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * current_time / t_ref_)) : u_ref_;
        if (aligned_box_.checkInBounds(position))
        {
            target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        }
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, global_resolution);
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody thermofluid_body(sph_system, makeShared<ThermofluidBody>("ThermofluidBody"));
    thermofluid_body.defineClosure<WeaklyCompressibleFluid, Viscosity, IsotropicDiffusion>(
        ConstructArgs(rho0_f, c_f), mu_f, ConstructArgs(diffusion_species_name, diffusion_coeff));
    thermofluid_body.generateParticles<BaseParticles, Lattice>();

    SolidBody thermosolid_body(sph_system, makeShared<ThermosolidBody>("ThermosolidBody"));
    thermosolid_body.defineMaterial<Solid>();
    thermosolid_body.generateParticles<BaseParticles, Lattice>();

    ObserverBody temperature_observer(sph_system, "FluidObserver");
    temperature_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation fluid_body_inner(thermofluid_body);
    InnerRelation solid_body_inner(thermosolid_body);
    ContactRelation fluid_body_contact(thermofluid_body, {&thermosolid_body});
    ContactRelation fluid_observer_contact(temperature_observer, {&thermofluid_body});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation fluid_body_complex(fluid_body_inner, fluid_body_contact);
    //----------------------------------------------------------------------
    // Define the numerical methods used in the simulation.
    // Note that there may be data dependence on the sequence of constructions.
    // Generally, the geometric models or simple objects without data dependencies,
    // such as gravity, should be initiated first.
    // Then the major physical particle dynamics model should be introduced.
    // Finally, the auxiliary models such as time step estimator, initial condition,
    // boundary condition and other constraints should be defined.
    // For typical fluid-structure interaction, we first define structure dynamics,
    // Then fluid dynamics and the corresponding coupling dynamics.
    // The coupling with multi-body dynamics will be introduced at last.
    //----------------------------------------------------------------------
    SimpleDynamics<NormalDirectionFromBodyShape> thermosolid_body_normal_direction(thermosolid_body);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(fluid_body_inner, fluid_body_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(fluid_body_inner, fluid_body_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_body_inner, fluid_body_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(fluid_body_inner, fluid_body_contact);
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(fluid_body_inner, fluid_body_contact);

    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step(thermofluid_body, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step(thermofluid_body);
    PeriodicAlongAxis periodic_along_x(thermofluid_body.getSPHBodyBounds(), xAxis);
    PeriodicConditionUsingCellLinkedList periodic_condition(thermofluid_body, periodic_along_x);
    AlignedBoxByCell inflow_buffer(thermofluid_body, AlignedBox(xAxis, Transform(Vec2d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);

    GetDiffusionTimeStepSize get_thermal_time_step(thermofluid_body);
    ThermalRelaxationComplex thermal_relaxation_complex(fluid_body_inner, fluid_body_contact);
    SimpleDynamics<ThermosolidBodyInitialCondition> thermosolid_condition(thermosolid_body);
    SimpleDynamics<ThermofluidBodyInitialCondition> thermofluid_initial_condition(thermofluid_body);

    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(fluid_body_inner);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(thermofluid_body);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Real>(thermofluid_body, "Phi");
    write_real_body_states.addToWrite<Real>(thermosolid_body, "Phi");
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>> write_fluid_phi("Phi", fluid_observer_contact);
    ObservedQuantityRecording<Vecd> write_fluid_velocity("Velocity", fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** computing surface normal direction for the wall. */
    thermosolid_body_normal_direction.exec();
    thermosolid_condition.exec();
    thermofluid_initial_condition.exec();
    Real dt_thermal = get_thermal_time_step.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    Real end_time = 10;
    Real output_interval = end_time / 100.0; /**< time stamps for output,WriteToFile*/
    int number_of_iterations = 0;
    int screen_output_interval = 40;
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
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step.exec();
            update_density_by_summation.exec();
            viscous_force.exec();
            transport_velocity_correction.exec();

            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(SMIN(dt_thermal, get_fluid_time_step.exec()), Dt);
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                thermal_relaxation_complex.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
                parabolic_inflow.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            thermofluid_body.updateCellLinkedList();
            periodic_condition.update_cell_linked_list_.exec();
            fluid_body_complex.updateConfiguration();
        }
        TickCount t2 = TickCount::now();
        /** write run-time observation into file */
        compute_vorticity.exec();
        fluid_observer_contact.updateConfiguration();
        write_real_body_states.writeToFile();
        write_fluid_phi.writeToFile(number_of_iterations);
        write_fluid_velocity.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }

    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_fluid_phi.generateDataBase(1.0e-3, 1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_fluid_phi.testResult();
    }

    return 0;
}
