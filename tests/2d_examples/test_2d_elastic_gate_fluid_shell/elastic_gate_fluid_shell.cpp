/**
 * @file 	elastic_gate.cpp
 * @brief 	2D elastic gate deformation due to dam break force.
 * @details This is the one of the basic test cases for
 * 			understanding SPH method for fluid-shell-interaction simulation.
 * @author 	Weiyi Kong
 */
#include "sphinxsys.h" //	SPHinXsys Library.
using namespace SPH;   //	Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real scale = 0.001;
Real DL = 500.0 * scale;                   /**< Tank length. */
Real DH = 200.1 * scale;                   /**< Tank height. */
Real Dam_L = 100.0 * scale;                /**< Water block width. */
Real Dam_H = 140.0 * scale;                /**< Water block height. */
Real Gate_width = 5.0 * scale;             /**< Width of the gate. */
Real Base_bottom_position = 79.0 * scale;  /**< Position of gate base. (In Y direction) */
Real resolution_gate = Gate_width / 4.0;   /**< Initial reference particle spacing of the gate. */
Real resolution_ref = 2 * resolution_gate; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4.0;            /**< Extending width for BCs. */
/** The offset that the rubber gate shifted above the tank. */
Vec2d offset = Vec2d(0.0, Base_bottom_position - floor(Base_bottom_position / resolution_gate) * resolution_gate);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
/** Define the corner points of the water block geometry. */
Vec2d DamP_lb(DL - Dam_L, 0.0);   /**< Left bottom. */
Vec2d DamP_lt(DL - Dam_L, Dam_H); /**< Left top. */
Vec2d DamP_rt(DL, Dam_H);         /**< Right top. */
Vec2d DamP_rb(DL, 0.0);           /**< Right bottom. */
/** Define the corner points of the gate constrain. */
Vec2d ConstrainP_lb(DL - Dam_L - Gate_width, Base_bottom_position); /**< Left bottom. */
Vec2d ConstrainP_lt(DL - Dam_L - Gate_width, Dam_H + BW);           /**< Left top. */
Vec2d ConstrainP_rt(DL - Dam_L, Dam_H + BW);                        /**< Right top. */
Vec2d ConstrainP_rb(DL - Dam_L, Base_bottom_position);              /**< Right bottom. */
// observer location
StdVec<Vecd> observation_location = {Vec2d(DL - Dam_L - resolution_gate, 0.0)};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                   /**< Reference density of fluid. */
Real gravity_g = 9.8;                   /**< Value of gravity. */
Real U_f = 2 * sqrt(Dam_H * gravity_g); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                  /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Material parameters of the elastic gate.
//----------------------------------------------------------------------
Real rho0_s = 1100;  /**< Reference density of gate. */
Real poisson = 0.47; /**< Poisson ratio. */
Real Youngs_modulus = 7.8e6;
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
class GateParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit GateParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        Real x = DL - Dam_L - 0.5 * resolution_gate;
        Real y = 0.5 * resolution_gate;
        while (y < Dam_H + BW)
        {
            initializePositionAndVolumetricMeasure(Vecd(x, y), resolution_gate);
            initializeSurfaceProperties(Vec2d(-1, 0), Gate_width);
            y += resolution_gate;
        }
    }
};
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
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();

    SolidBody gate(sph_system, makeShared<DefaultShape>("Gate"));
    gate.defineAdaptationRatios(1.15, resolution_ref / resolution_gate);
    gate.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, Youngs_modulus, poisson);
    gate.generateParticles<GateParticleGenerator>();

    ObserverBody gate_observer(sph_system, "Observer");
    gate_observer.defineAdaptationRatios(1.15, resolution_ref / resolution_gate);
    gate_observer.generateParticles<ParticleGeneratorObserver>(observation_location);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    ContactRelation water_wall_contact(water_block, {&wall_boundary});
    ContactRelationToShell water_gate_contact(water_block, {&gate}, {false});
    InnerRelation gate_inner(gate);
    ContactRelationFromShell gate_water_contact(gate, {&water_block}, {false});
    ShellInnerRelationWithContactKernel gate_curvature_inner(gate, water_block);
    ContactRelation gate_observer_contact(gate_observer, {&gate});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, {&water_wall_contact, &water_gate_contact});
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Algorithms of fluid dynamics.
    //----------------------------------------------------------------------
    Gravity gravity(Vecd(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(water_block, gravity);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> pressure_relaxation(water_block_inner, water_wall_contact, water_gate_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver>> density_relaxation(water_block_inner, water_wall_contact, water_gate_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeSurface>, Contact<>, Contact<>>> update_density_by_summation(water_block_inner, water_wall_contact, water_gate_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    //----------------------------------------------------------------------
    //	Algorithms of FSI.
    //----------------------------------------------------------------------
    SimpleDynamics<OffsetInitialPosition> gate_offset_position(gate, offset);
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> gate_corrected_configuration(gate_inner);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> fluid_pressure_force_on_gate(gate_water_contact);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(gate);
    //----------------------------------------------------------------------
    //	Algorithms of Elastic dynamics.
    //----------------------------------------------------------------------
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> gate_stress_relaxation_first_half(gate_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> gate_stress_relaxation_second_half(gate_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> gate_computing_time_step_size(gate);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> gate_average_curvature(gate_curvature_inner);

    BodyRegionByParticle gate_constraint_part(gate, makeShared<MultiPolygonShape>(createGateConstrainShape()));
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> gate_constraint(gate_constraint_part);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> gate_update_normal(gate);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<Real>("Pressure");
    gate.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    BodyStatesRecordingToVtp write_real_body_states_to_vtp(sph_system.real_bodies_);
    ObservedQuantityRecording<Vecd> write_beam_tip_displacement("Displacement", gate_observer_contact);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    gate_offset_position.exec();
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    gate_corrected_configuration.exec();
    gate_average_curvature.exec();
    water_block_complex.updateConfiguration();
    gate_water_contact.updateConfiguration();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    int number_of_iterations = 0;
    int screen_output_interval = 100;
    Real end_time = 0.4;
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
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        /** Integrate time (loop) until the next output time. */
        while (integration_time < output_interval)
        {
            /** Force Prior due to viscous force and gravity. */
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();

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
                    dt_s = 0.5 * gate_computing_time_step_size.exec();
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
                GlobalStaticVariables::physical_time_ += dt;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	dt = " << dt << "	dt_s = " << dt_s << "\n";
            }
            number_of_iterations++;

            /** Update normal direction at elastic body surface. */
            gate_update_normal.exec();

            /** Update cell linked list and configuration. */
            water_block.updateCellLinkedListWithParticleSort(100);
            gate.updateCellLinkedList();
            gate_curvature_inner.updateConfiguration();
            gate_average_curvature.exec();
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

    return 0;
}
