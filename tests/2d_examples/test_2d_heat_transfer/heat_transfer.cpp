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
Real resolution_ref = DH / 25.0;        /**< Global reference resolution. */
Real DL_sponge = resolution_ref * 20.0; /**< Sponge region to impose inflow condition. */
/** Boundary width, determined by specific layer of boundary particles. */
Real BW = resolution_ref * 4.0; /** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - BW, -BW), Vec2d(DL + BW, DH + BW));
// temperature observer location
StdVec<Vecd> observation_location = {Vecd(0.0, DH * 0.5)};
//----------------------------------------------------------------------
//	Global parameters on the material properties
//----------------------------------------------------------------------
Real diffusion_coff = 1.0e-3;
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
//	Setup heat conduction material properties for diffusion fluid body
//----------------------------------------------------------------------
class ThermofluidBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
  public:
    ThermofluidBodyMaterial()
        : DiffusionReaction<WeaklyCompressibleFluid>({"Phi"}, SharedPtr<NoReaction>(), rho0_f, c_f, mu_f)
    {
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
    };
};
using DiffusionBaseParticles = DiffusionReactionParticles<BaseParticles, ThermofluidBodyMaterial>;
//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion solid body
//----------------------------------------------------------------------
class ThermosolidBodyMaterial : public DiffusionReaction<Solid>
{
  public:
    ThermosolidBodyMaterial() : DiffusionReaction<Solid>({"Phi"}, SharedPtr<NoReaction>())
    {
        // only default property is given, as no heat transfer within solid considered here.
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi");
    };
};
using DiffusionSolidParticles = DiffusionReactionParticles<SolidParticles, ThermosolidBodyMaterial>;
//----------------------------------------------------------------------
//	Application dependent solid body initial condition
//----------------------------------------------------------------------
class ThermosolidBodyInitialCondition
    : public DiffusionReactionInitialCondition<DiffusionSolidParticles>
{
  protected:
    size_t phi_;

  public:
    explicit ThermosolidBodyInitialCondition(SPHBody &sph_body)
        : DiffusionReactionInitialCondition<DiffusionSolidParticles>(sph_body)
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    };

    void update(size_t index_i, Real dt)
    {
        if (-BW <= pos_[index_i][1] && pos_[index_i][1] <= 0.0)
        {
            all_species_[phi_][index_i] = phi_lower_wall;
        }

        if (DH <= pos_[index_i][1] && pos_[index_i][1] <= DH + BW)
        {
            all_species_[phi_][index_i] = phi_upper_wall;
        }
    };
};
//----------------------------------------------------------------------
//	Application dependent fluid body initial condition
//----------------------------------------------------------------------
class ThermofluidBodyInitialCondition
    : public DiffusionReactionInitialCondition<DiffusionBaseParticles>
{
  protected:
    size_t phi_;

  public:
    explicit ThermofluidBodyInitialCondition(SPHBody &sph_body)
        : DiffusionReactionInitialCondition<DiffusionBaseParticles>(sph_body)
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    };

    void update(size_t index_i, Real dt)
    {
        if (0 <= pos_[index_i][1] && pos_[index_i][1] <= DH)
        {
            all_species_[phi_][index_i] = phi_fluid_initial;
        }
    };
};

using FluidDiffusionInner = DiffusionRelaxationInner<DiffusionBaseParticles>;
using FluidDiffusionDirichlet = DiffusionRelaxationDirichlet<DiffusionBaseParticles, DiffusionSolidParticles>;
//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
class ThermalRelaxationComplex
    : public DiffusionRelaxationRK2<
          ComplexInteraction<FluidDiffusionInner, FluidDiffusionDirichlet>>
{
  public:
    explicit ThermalRelaxationComplex(BaseInnerRelation &inner_relation, BaseContactRelation &body_contact_relation_Dirichlet)
        : DiffusionRelaxationRK2<ComplexInteraction<FluidDiffusionInner, FluidDiffusionDirichlet>>(inner_relation, body_contact_relation_Dirichlet){};
    virtual ~ThermalRelaxationComplex(){};
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
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        if (aligned_box_.checkInBounds(0, position))
        {
            target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        }
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, resolution_ref);
    GlobalStaticVariables::physical_time_ = 0.0;
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody thermofluid_body(system, makeShared<ThermofluidBody>("ThermofluidBody"));
    thermofluid_body.defineParticlesAndMaterial<DiffusionBaseParticles, ThermofluidBodyMaterial>();
    thermofluid_body.generateParticles<ParticleGeneratorLattice>();

    SolidBody thermosolid_body(system, makeShared<ThermosolidBody>("ThermosolidBody"));
    thermosolid_body.defineParticlesAndMaterial<DiffusionSolidParticles, ThermosolidBodyMaterial>();
    thermosolid_body.generateParticles<ParticleGeneratorLattice>();

    ObserverBody temperature_observer(system, "FluidObserver");
    temperature_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation fluid_body_inner(thermofluid_body);
    InnerRelation solid_body_inner(thermosolid_body);
    ContactRelation fluid_wall_contact_Dirichlet(thermofluid_body, {&thermosolid_body});
    ComplexRelation fluid_body_complex(fluid_body_inner, {&thermosolid_body});
    ContactRelation fluid_observer_contact(temperature_observer, {&thermofluid_body});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    PeriodicConditionUsingCellLinkedList periodic_condition(thermofluid_body, thermofluid_body.getBodyShapeBounds(), xAxis);
    SimpleDynamics<ThermosolidBodyInitialCondition> thermosolid_condition(thermosolid_body);
    SimpleDynamics<ThermofluidBodyInitialCondition> thermofluid_initial_condition(thermofluid_body);
    SimpleDynamics<NormalDirectionFromBodyShape> thermosolid_body_normal_direction(thermosolid_body);
    /** Initialize particle acceleration. */
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(thermofluid_body);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplex> update_density_by_summation(fluid_body_complex);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step(thermofluid_body, U_f);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step(thermofluid_body);
    /** Time step size calculation. */
    GetDiffusionTimeStepSize<DiffusionBaseParticles> get_thermal_time_step(thermofluid_body);
    /** Diffusion process between two diffusion bodies. */
    ThermalRelaxationComplex thermal_relaxation_complex(fluid_body_inner, fluid_wall_contact_Dirichlet);
    /** Pressure relaxation using verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(fluid_body_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(fluid_body_complex);
    /** Computing viscous acceleration. */
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(fluid_body_complex);
    /** Apply transport velocity formulation. */
    InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionComplex<AllParticles>> transport_velocity_correction(fluid_body_complex);
    /** Computing vorticity in the flow. */
    InteractionDynamics<fluid_dynamics::VorticityInner> compute_vorticity(fluid_body_inner);
    /** Inflow boundary condition. */
    BodyAlignedBoxByCell inflow_buffer(
        thermofluid_body, makeShared<AlignedBoxShape>(Transform(Vec2d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(io_environment, system.real_bodies_);
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>>
        write_fluid_phi("Phi", io_environment, fluid_observer_contact);
    ObservedQuantityRecording<Vecd>
        write_fluid_velocity("Velocity", io_environment, fluid_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    system.initializeSystemCellLinkedLists();
    /** periodic condition applied after the mesh cell linked list build up
     * but before the configuration build up. */
    periodic_condition.update_cell_linked_list_.exec();
    /** initialize configurations for all bodies. */
    system.initializeSystemConfigurations();
    /** computing surface normal direction for the wall. */
    thermosolid_body_normal_direction.exec();
    thermosolid_condition.exec();
    thermofluid_initial_condition.exec();
    Real dt_thermal = get_thermal_time_step.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
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
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
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
                GlobalStaticVariables::physical_time_ += dt;
                parabolic_inflow.exec();
                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** Water block configuration and periodic condition. */
            periodic_condition.bounding_.exec();
            thermofluid_body.updateCellLinkedListWithParticleSort(100);
            periodic_condition.update_cell_linked_list_.exec();
            fluid_body_complex.updateConfiguration();

            fluid_body_inner.updateConfiguration();
            fluid_wall_contact_Dirichlet.updateConfiguration();
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

    write_fluid_phi.testResult();

    return 0;
}
