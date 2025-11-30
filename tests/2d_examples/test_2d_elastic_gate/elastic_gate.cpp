/**
 * @file 	elastic_gate.cpp
 * @brief 	2D elastic gate deformation due to dam break force.
 * @details This is the one of the basic test cases, also the first case for
 * 			understanding SPH method for fluid-structure-interaction (FSI) simulation.
 * @author 	Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 500.0;                        /**< Tank length. */
Real DH = 200.1;                        /**< Tank height. */
Real Dam_L = 100.0;                     /**< Water block width. */
Real Dam_H = 140.0;                     /**< Water block height. */
Real Gate_width = 5.0;                  /**< Width of the gate. */
Real Base_bottom_position = 79.0;       /**< Position of gate base. (In Y direction) */
Real resolution_ref = Gate_width / 2.0; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4.0;         /**< Extending width for BCs. */
/** The offset that the rubber gate shifted above the tank. */
Real dp_s = 0.5 * resolution_ref;
Vec2d offset = Vec2d(0.0, Base_bottom_position - floor(Base_bottom_position / dp_s) * dp_s);
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/** Define the corner points of the water block geometry. */
Vec2d DamP_lb(DL - Dam_L, 0.0);   /**< Left bottom. */
Vec2d DamP_lt(DL - Dam_L, Dam_H); /**< Left top. */
Vec2d DamP_rt(DL, Dam_H);         /**< Right top. */
Vec2d DamP_rb(DL, 0.0);           /**< Right bottom. */
/** Define the corner points of the gate geometry. */
Vec2d GateP_lb(DL - Dam_L - Gate_width, 0.0);        /**< Left bottom. */
Vec2d GateP_lt(DL - Dam_L - Gate_width, Dam_H + BW); /**< Left top. */
Vec2d GateP_rt(DL - Dam_L, Dam_H + BW);              /**< Right top. */
Vec2d GateP_rb(DL - Dam_L, 0.0);                     /**< Right bottom. */
/** Define the corner points of the gate constrain. */
Vec2d ConstrainP_lb(DL - Dam_L - Gate_width, Base_bottom_position); /**< Left bottom. */
Vec2d ConstrainP_lt(DL - Dam_L - Gate_width, Dam_H + BW);           /**< Left top. */
Vec2d ConstrainP_rt(DL - Dam_L, Dam_H + BW);                        /**< Right top. */
Vec2d ConstrainP_rb(DL - Dam_L, Base_bottom_position);              /**< Right bottom. */
// observer location
StdVec<Vecd> observation_location = {GateP_lb};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                         /**< Reference density of fluid. */
Real gravity_g = 9.8e-3;                   /**< Value of gravity. */
Real U_f = 1.0;                            /**< Characteristic velocity. */
Real c_f = 20.0 * sqrt(140.0 * gravity_g); /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Material parameters of the elastic gate.
//----------------------------------------------------------------------
Real rho0_s = 1.1;   /**< Reference density of gate. */
Real poisson = 0.47; /**< Poisson ratio. */
Real Ae = 7.8e3;     /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae * rho0_f * U_f * U_f;
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(DamP_lb);
        water_block_shape.push_back(DamP_lt);
        water_block_shape.push_back(DamP_rt);
        water_block_shape.push_back(DamP_rb);
        water_block_shape.push_back(DamP_lb);
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Wall cases-dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(-BW, -BW));
        outer_wall_shape.push_back(Vecd(-BW, DH + BW));
        outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
        outer_wall_shape.push_back(Vecd(DL + BW, -BW));
        outer_wall_shape.push_back(Vecd(-BW, -BW));

        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(0.0, 0.0));
        inner_wall_shape.push_back(Vecd(0.0, DH));
        inner_wall_shape.push_back(Vecd(DL, DH));
        inner_wall_shape.push_back(Vecd(DL, 0.0));
        inner_wall_shape.push_back(Vecd(0.0, 0.0));

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	create a gate shape
//----------------------------------------------------------------------
MultiPolygon createGateShape()
{
    std::vector<Vecd> gate_shape;
    gate_shape.push_back(GateP_lb);
    gate_shape.push_back(GateP_lt);
    gate_shape.push_back(GateP_rt);
    gate_shape.push_back(GateP_rb);
    gate_shape.push_back(GateP_lb);

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(gate_shape, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
// Create the gate constrain shape
//----------------------------------------------------------------------
MultiPolygon createGateConstrainShape()
{
    // geometry
    std::vector<Vecd> gate_constraint_shape;
    gate_constraint_shape.push_back(ConstrainP_lb);
    gate_constraint_shape.push_back(ConstrainP_lt);
    gate_constraint_shape.push_back(ConstrainP_rt);
    gate_constraint_shape.push_back(ConstrainP_rb);
    gate_constraint_shape.push_back(ConstrainP_lb);

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(gate_constraint_shape, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineMaterial<WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody gate(sph_system, makeShared<MultiPolygonShape>(createGateShape(), "Gate"));
    gate.defineAdaptationRatios(1.15, 2.0);
    gate.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    gate.generateParticles<BaseParticles, Lattice>();

    ObserverBody gate_observer(sph_system, "Observer");
    gate_observer.defineAdaptationRatios(1.15, 2.0);
    gate_observer.generateParticles<ObserverParticles>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_block_contact(water_block, RealBodyVector{&wall_boundary, &gate});
    InnerRelation gate_inner(gate);
    ContactRelation gate_water_contact(gate, {&water_block});
    ContactRelation gate_observer_contact(gate_observer, {&gate});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, water_block_contact);
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
    SimpleDynamics<OffsetInitialPosition> gate_offset_position(gate, offset);
    InteractionWithUpdate<LinearGradientCorrectionMatrixInner> gate_corrected_configuration(gate_inner);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> gate_normal_direction(gate);

    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> gate_stress_relaxation_first_half(gate_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> gate_stress_relaxation_second_half(gate_inner);

    ReduceDynamics<solid_dynamics::AcousticTimeStep> gate_computing_time_step_size(gate);
    BodyRegionByParticle gate_constraint_part(gate, makeShared<MultiPolygonShape>(createGateConstrainShape()));
    SimpleDynamics<FixBodyPartConstraint> gate_constraint(gate_constraint_part);
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> gate_update_normal(gate);

    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);

    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::DensitySummationComplexFreeSurface> update_density_by_summation(water_block_inner, water_block_contact);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> fluid_pressure_force_on_gate(gate_water_contact);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(gate);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToPlt write_real_body_states_to_plt(sph_system);
    BodyStatesRecordingToVtp write_real_body_states_to_vtp(sph_system);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
        write_beam_tip_displacement("Position", gate_observer_contact);
    // TODO: observing position is not as good observing displacement.
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    gate_offset_position.exec();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    gate_normal_direction.exec();
    gate_corrected_configuration.exec();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    int number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 400.0;
    Real output_interval = end_time / 200.0;
    Real dt = 0.0;   /**< Default acoustic time step sizes. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states_to_vtp.writeToFile();
    write_beam_tip_displacement.writeToFile();
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            /** Update normal direction at elastic body surface. */
            gate_update_normal.exec();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                /** Fluid relaxation and force computation. */
                pressure_relaxation.exec(dt);
                fluid_pressure_force_on_gate.exec();
                density_relaxation.exec(dt);
                /** Solid dynamics time stepping. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    dt_s = gate_computing_time_step_size.exec();
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    gate_stress_relaxation_first_half.exec(dt_s);
                    gate_constraint.exec();
                    gate_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                }
                average_velocity_and_acceleration.update_averages_.exec(dt);
                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();
            gate.updateCellLinkedList();
            water_block_complex.updateConfiguration();
            gate_water_contact.updateConfiguration();
            /** Output the observed data. */
            write_beam_tip_displacement.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        write_real_body_states_to_vtp.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();
    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_beam_tip_displacement.generateDataBase(1.0e-3);
    }
    else if (sph_system.RestartStep() == 0)
    {
        write_beam_tip_displacement.testResult();
    }

    return 0;
}
