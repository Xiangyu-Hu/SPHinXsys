/**
 * @file 	hydrostatic_fsi.cpp
 * @brief 	structure deformation due to hydrostatic pressure under gravity.
 * @details This is the one of the basic test cases
 * for understanding SPH method for fluid-structure-interaction (FSI) simulation.
 * @author 	Yujie Zhu, Chi Zhang and Xiangyu Hu
 * @version 0.1
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 1.0;                                /**< Tank length. */
Real DH = 2.1;                                /**< Tank height. */
Real Dam_L = 1.0;                             /**< Water block width. */
Real Dam_H = 2.0;                             /**< Water block height. */
Real Gate_width = 0.05;                       /**< Width of the gate. */
Real particle_spacing_ref = Gate_width / 4.0; /**< Initial reference particle spacing. 8, 10, 12 */
Real BW = 4.0 * particle_spacing_ref;         /**< Extending width for BCs. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Define the corner point of water block geometry.
//----------------------------------------------------------------------
Vec2d DamP_lb(0.0, 0.0);     /**< Left bottom. */
Vec2d DamP_lt(0.0, Dam_H);   /**< Left top. */
Vec2d DamP_rt(Dam_L, Dam_H); /**< Right top. */
Vec2d DamP_rb(Dam_L, 0.0);   /**< Right bottom. */
//----------------------------------------------------------------------
//	Define the corner point of gate geometry.
//----------------------------------------------------------------------
Vec2d GateP_lb(-BW, -Gate_width);
Vec2d GateP_lt(-BW, 0.0);
Vec2d GateP_rt(Dam_L + BW, 0.0);
Vec2d GateP_rb(Dam_L + BW, -Gate_width);
//----------------------------------------------------------------------
//	Define the geometry for gate constrain.
//----------------------------------------------------------------------
Vec2d ConstrainLP_lb(-BW, -Gate_width);
Vec2d ConstrainLP_lt(-BW, 0.0);
Vec2d ConstrainLP_rt(0.0, 0.0);
Vec2d ConstrainLP_rb(0.0, -Gate_width);
Vec2d ConstrainRP_lb(Dam_L, -Gate_width);
Vec2d ConstrainRP_lt(Dam_L, 0.0);
Vec2d ConstrainRP_rt(Dam_L + BW, 0.0);
Vec2d ConstrainRP_rb(Dam_L + BW, -Gate_width);
// observer location
StdVec<Vecd> observation_location = {Vecd(0.5 * Dam_L, -0.5 * Gate_width)};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;  /**< Reference density of fluid. */
Real gravity_g = 9.81; /**< Value of gravity. */
Real U_max = 2.0 * sqrt(Dam_H * gravity_g);
;                                     /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;              /**< Reference sound speed. */
Real Re = 0.1;                        /**< Reynolds number. */
Real mu_f = rho0_f * U_max * DL / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Material properties of the elastic gate.
//----------------------------------------------------------------------
Real rho0_s = 2700.0; /**< Reference solid density. */
Real poisson = 0.34;  /**< Poisson ratio. */
Real Ae = 6.75e10;    /**< Normalized Youngs Modulus. */
Real Youngs_modulus = Ae;
//----------------------------------------------------------------------
//	Geometry definition.
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(DamP_lb);
    water_block_shape.push_back(DamP_lt);
    water_block_shape.push_back(DamP_rt);
    water_block_shape.push_back(DamP_rb);
    water_block_shape.push_back(DamP_lb);

    return water_block_shape;
}
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	wall body shape definition.
//----------------------------------------------------------------------
std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> outer_wall_shape;
    outer_wall_shape.push_back(Vecd(-BW, 0.0));
    outer_wall_shape.push_back(Vecd(-BW, DH));
    outer_wall_shape.push_back(Vecd(0.0, DH));
    outer_wall_shape.push_back(Vecd(0.0, 0.0));
    outer_wall_shape.push_back(Vecd(-BW, 0.0));

    return outer_wall_shape;
}

std::vector<Vecd> createInnerWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(DL, 0.0));
    inner_wall_shape.push_back(Vecd(DL, DH));
    inner_wall_shape.push_back(Vecd(DL + BW, DH));
    inner_wall_shape.push_back(Vecd(DL + BW, 0.0));
    inner_wall_shape.push_back(Vecd(DL, 0.0));

    return inner_wall_shape;
}

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	create a gate shape
//----------------------------------------------------------------------
std::vector<Vecd> createGateShape()
{
    std::vector<Vecd> gate_shape;
    gate_shape.push_back(GateP_lb);
    gate_shape.push_back(GateP_lt);
    gate_shape.push_back(GateP_rt);
    gate_shape.push_back(GateP_rb);
    gate_shape.push_back(GateP_lb);

    return gate_shape;
}
//----------------------------------------------------------------------
//	Define the gate body shape.
//----------------------------------------------------------------------
class Gate : public MultiPolygonShape
{
  public:
    explicit Gate(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createGateShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	create left Gate constrain shape
//----------------------------------------------------------------------
MultiPolygon createGateConstrainShape()
{
    // geometry
    std::vector<Vecd> gate_constraint_shape_left;
    gate_constraint_shape_left.push_back(ConstrainLP_lb);
    gate_constraint_shape_left.push_back(ConstrainLP_lt);
    gate_constraint_shape_left.push_back(ConstrainLP_rt);
    gate_constraint_shape_left.push_back(ConstrainLP_rb);
    gate_constraint_shape_left.push_back(ConstrainLP_lb);

    std::vector<Vecd> gate_constraint_shape_right;
    gate_constraint_shape_right.push_back(ConstrainRP_lb);
    gate_constraint_shape_right.push_back(ConstrainRP_lt);
    gate_constraint_shape_right.push_back(ConstrainRP_rt);
    gate_constraint_shape_right.push_back(ConstrainRP_rb);
    gate_constraint_shape_right.push_back(ConstrainRP_lb);

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(gate_constraint_shape_left, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(gate_constraint_shape_right, ShapeBooleanOps::add);
    return multi_polygon;
}
//----------------------------------------------------------------------
//	create right Gate constrain shape
//----------------------------------------------------------------------
std::vector<Vecd> createGateConstrainShapeRight()
{
    // geometry
    std::vector<Vecd> gate_constraint_shape;
    gate_constraint_shape.push_back(ConstrainRP_lb);
    gate_constraint_shape.push_back(ConstrainRP_lt);
    gate_constraint_shape.push_back(ConstrainRP_rt);
    gate_constraint_shape.push_back(ConstrainRP_rb);
    gate_constraint_shape.push_back(ConstrainRP_lb);

    return gate_constraint_shape;
}
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
    //----------------------------------------------------------------------
    //	Build up -- a SPHSystem
    //----------------------------------------------------------------------
    SPHSystem system(system_domain_bounds, particle_spacing_ref);
    /** Set the starting time to zero. */
    GlobalStaticVariables::physical_time_ = 0.0;
    IOEnvironment io_environment(system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    SolidBody gate(system, makeShared<Gate>("Gate"));
    gate.defineParticlesAndMaterial<ElasticSolidParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    gate.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Particle and body creation of gate observer.
    //----------------------------------------------------------------------
    ObserverBody gate_observer(system, "Observer");
    gate_observer.generateParticles<ObserverParticleGenerator>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation gate_inner(gate);
    ComplexRelation water_block_complex(water_block_inner, {&wall_boundary, &gate});
    ContactRelation gate_contact(gate, {&water_block});
    ContactRelation gate_observer_contact(gate_observer, {&gate});
    //----------------------------------------------------------------------
    //	Define all numerical methods which are used in this case.
    //----------------------------------------------------------------------
    /** Initialize particle acceleration. */ // TODO: this is missing for solid body.
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block, makeShared<Gravity>(Vecd(0.0, -gravity_g)));
    /** Evaluation of fluid density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> update_fluid_density(water_block_complex);
    /** Compute time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_max);
    /** Compute time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    /** Pressure relaxation using verlet time stepping. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> pressure_relaxation(water_block_complex);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWall> density_relaxation(water_block_complex);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseWithWall<Vec2d, DampingPairwiseInner>>>
        fluid_damping(0.2, water_block_complex, "Velocity", mu_f);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> gate_normal_direction(gate);
    /** Corrected configuration. */
    InteractionWithUpdate<CorrectedConfigurationInner> gate_corrected_configuration(gate_inner);
    /** Compute time step size of elastic solid. */
    ReduceDynamics<solid_dynamics::AcousticTimeStepSize> gate_computing_time_step_size(gate);
    /** Stress relaxation stepping for the elastic gate. */
    Dynamics1Level<solid_dynamics::Integration1stHalfPK2> gate_stress_relaxation_first_half(gate_inner);
    Dynamics1Level<solid_dynamics::Integration2ndHalf> gate_stress_relaxation_second_half(gate_inner);
    /**Constrain a solid body part.  */
    BodyRegionByParticle gate_constraint_part(gate, makeShared<MultiPolygonShape>(createGateConstrainShape()));
    SimpleDynamics<solid_dynamics::FixBodyPartConstraint> gate_constraint(gate_constraint_part);
    /** Update the norm of elastic gate. */
    SimpleDynamics<solid_dynamics::UpdateElasticNormalDirection> gate_update_normal(gate);
    /** Compute the average velocity of gate. */
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(gate);
    /** Compute the force exerted on elastic gate due to fluid pressure. */
    InteractionDynamics<solid_dynamics::PressureForceAccelerationFromFluid> fluid_pressure_force_on_gate(gate_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    /** Output body states for visualization. */
    BodyStatesRecordingToVtp write_real_body_states_to_vtp(io_environment, system.real_bodies_);
    /** Output the observed displacement of gate free end. */
    RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
        write_beam_tip_displacement("Position", io_environment, gate_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    system.initializeSystemConfigurations();
    /** computing surface normal direction for the wall. */
    wall_boundary_normal_direction.exec();
    /** computing surface normal direction for the insert body. */
    gate_normal_direction.exec();
    /** computing linear reproducing configuration for the insert body. */
    gate_corrected_configuration.exec();
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states_to_vtp.writeToFile(0);
    write_beam_tip_displacement.writeToFile(0);
    //----------------------------------------------------------------------
    //	Basic control parameters for time stepping.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.5; /**< End time. */
    Real output_interval = end_time / 50.0;
    Real dt = 0.0;   /**< Default acoustic time step sizes. */
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	Main loop of time stepping starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Acceleration due to viscous force and gravity. */
            initialize_a_fluid_step.exec();
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_fluid_density.exec();
            /** Update normal direction on elastic body. */
            gate_update_normal.exec();
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                dt = SMIN(get_fluid_time_step_size.exec(), Dt);
                fluid_damping.exec(dt);
                /** Fluid relaxation and force computation. */
                pressure_relaxation.exec(dt);
                fluid_pressure_force_on_gate.exec();
                density_relaxation.exec(dt);
                /** Solid dynamics time stepping. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    gate_stress_relaxation_first_half.exec(dt_s);
                    gate_constraint.exec();
                    gate_stress_relaxation_second_half.exec(dt_s);
                    dt_s_sum += dt_s;
                    dt_s = gate_computing_time_step_size.exec();
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
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedList(); // water particle motion is small
            water_block_complex.updateConfiguration();
            /** one need update configuration after periodic condition. */
            gate.updateCellLinkedList();
            gate_contact.updateConfiguration();

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

    write_beam_tip_displacement.testResult();

    return 0;
}
