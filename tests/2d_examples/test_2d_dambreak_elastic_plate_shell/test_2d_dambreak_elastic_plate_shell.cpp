/**
 * @file dambreak.cpp
 * @brief 2D dambreak example.
 * @details This is the one of the basic test cases, also the first case for
 * understanding SPH method for free surface flow simulation.
 * @author Luhui Han, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

// time
const Real gate_moving_time = 0.1;
const Real end_time = 0.5;
// general parameters for geometry
const Real t = 4e-3;                                        // gate thickness
const Real plate_height = 0.09;                             // plate height
const Real LH = 0.2;                                        // liquid initial height
const Real LL = 0.2;                                        // liquid initial length
const Real DH = 0.4;                                        // tank height
const Real DL = 0.8;                                        // tank length
const Real resolution_shell = t / 4.0;                      // shell particle spacing
const Real resolution_ref = 2 * resolution_shell;           // system particle spacing
const Real BW = resolution_ref * 4;                         // boundary width
const Real plate_x_pos = DL - 0.2 + 0.5 * resolution_shell; // center x coordiniate of plate

const Real marker_h = 0.0875; // height of marker
const std::vector<Vec2d> observer_position_1 = {Vec2d(plate_x_pos, marker_h)};

const BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));

// for material properties of the fluid
const Real rho0_f = 997.0;
const Real mu_f = 8.93e-7 * rho0_f;
const Real gravity_g = 9.8;
const Real U_f = 2.0 * sqrt(gravity_g * LH);
const Real c_f = 10.0 * U_f;

// material properties of the plate
const Real rho0_s = 1161.54; /**< Reference density.*/
const Real youngs_modulus = 3.5e6;
const Real poisson_ratio = 0.49;
const Real physical_viscosity = 0.4 / 4.0 * std::sqrt(rho0_s * youngs_modulus) * t * t;

//	define the water block shape
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        Vec2d halfsize_water(0.5 * LL, 0.5 * LH);
        Transform translation_water(halfsize_water);
        multi_polygon_.addABox(translation_water, halfsize_water, ShapeBooleanOps::add);
    }
};
//	define the static solid wall boundary shape
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        Vec2d halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW);
        Vec2d halfsize_inner(0.5 * DL, 0.5 * DH);
        Transform translation_wall(halfsize_inner);
        multi_polygon_.addABox(translation_wall, halfsize_outer, ShapeBooleanOps::add);
        multi_polygon_.addABox(translation_wall, halfsize_inner, ShapeBooleanOps::sub);

        // subtract the plate from wall to avoid overlapping
        Vec2d halfsize_plate(0.5 * resolution_ref, 0.5 * BW);
        Transform translation_plate(halfsize_plate + Vec2d(plate_x_pos, -BW));
        multi_polygon_.addABox(translation_plate, halfsize_plate, ShapeBooleanOps::sub);
    }
};
//	define the rigid solid gate shape
class MovingGate : public MultiPolygonShape
{
  public:
    explicit MovingGate(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        Vec2d halfsize_gate(0.5 * BW, 0.5 * DH);
        Transform translation_gate(Vec2d(LL, 0) + halfsize_gate);
        multi_polygon_.addABox(translation_gate, halfsize_gate, ShapeBooleanOps::add);
    }
};
//	define the elastic plate shape
class PlateParticleGenerator : public ParticleGeneratorSurface
{
  public:
    explicit PlateParticleGenerator(SPHBody &sph_body) : ParticleGeneratorSurface(sph_body){};
    void initializeGeometricVariables() override
    {
        Real y = -BW + 0.5 * resolution_shell;
        while (y < plate_height)
        {
            initializePositionAndVolumetricMeasure(Vec2d(plate_x_pos, y), resolution_shell);
            initializeSurfaceProperties(Vec2d(1, 0), t);
            y += resolution_shell;
        }
    }
};
// Define boundary geometry

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
        if (base_particles_.pos_[index_i].y() <= 0)
        {
            body_part_particles_.push_back(index_i);
        }
    };
};

// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("Wall"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary.generateParticles<ParticleGeneratorLattice>();
    wall_boundary.addBodyStateForRecording<Vec2d>("NormalDirection");

    SolidBody gate(sph_system, makeShared<MovingGate>("Gate"));
    gate.defineParticlesAndMaterial<SolidParticles, Solid>();
    gate.generateParticles<ParticleGeneratorLattice>();
    gate.addBodyStateForRecording<Vec2d>("NormalDirection");

    SolidBody plate(sph_system, makeShared<DefaultShape>("Plate"));
    plate.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    plate.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    plate.generateParticles<PlateParticleGenerator>();

    ObserverBody disp_observer_1(sph_system, "Observer1");
    disp_observer_1.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    disp_observer_1.generateParticles<ParticleGeneratorObserver>(observer_position_1);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation plate_inner(plate);
    ContactRelation water_wall_contact(water_block, {&wall_boundary, &gate});
    ContactRelationToShell water_plate_contact(water_block, {&plate}, {false});
    ContactRelationFromShell plate_water_contact(plate, {&water_block}, {false});
    ShellInnerRelationWithContactKernel plate_curvature_inner(plate, water_block);
    ContactRelation disp_observer_contact_1(disp_observer_1, {&plate});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, {&water_wall_contact, &water_plate_contact});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    // fluid
    Gravity gravity(Vec2d(0.0, -gravity_g));
    SimpleDynamics<GravityForce> constant_gravity(water_block, gravity);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> pressure_relaxation(water_block_inner, water_wall_contact, water_plate_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver>> density_relaxation(water_block_inner, water_wall_contact, water_plate_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeSurface>, Contact<>, Contact<>>> update_density_by_summation(water_block_inner, water_wall_contact, water_plate_contact);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::ViscousForce<Inner<>, Contact<Wall>, Contact<Wall>>, fluid_dynamics::FixedViscosity>> viscous_acceleration(water_block_inner, water_wall_contact, water_plate_contact);
    // solid
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> gate_normal_direction(gate);
    auto update_gate_position = [&](Real run_time)
    {
        const auto &pos0 = *gate.getBaseParticles().getVariableByName<Vec2d>("InitialPosition");
        Real h_g = -285.115 * run_time * run_time * run_time + 72.305 * run_time * run_time + 0.1463 * run_time;
        particle_for(
            par,
            gate.getBaseParticles().total_real_particles_,
            [&](size_t index_i)
            {
                gate.getBaseParticles().pos_[index_i][1] = pos0[index_i][1] + h_g;
            });
    };
    // Shell
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> plate_time_step_size(plate);
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> plate_corrected_configuration(plate_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> plate_stress_relaxation_first(plate_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> plate_stress_relaxation_second(plate_inner);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> plate_average_curvature(plate_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> plate_update_normal(plate);
    /** constraint and damping */
    BoundaryGeometry plate_boundary_geometry(plate, "BoundaryGeometry");
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> plate_constraint(plate_boundary_geometry);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        plate_position_damping(0.2, plate_inner, "Velocity", physical_viscosity);
    DampingWithRandomChoice<InteractionSplit<DampingPairwiseInner<Vec2d>>>
        plate_rotation_damping(0.2, plate_inner, "AngularVelocity", physical_viscosity);
    // FSI
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_plate(plate_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluidRiemann> pressure_force_on_plate(plate_water_contact);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(plate);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    water_block.addBodyStateForRecording<Real>("Pressure");
    plate.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    plate.addBodyStateForRecording<Vec2d>("PressureForceFromFluid");
    BodyStatesRecordingToVtp write_water_block_states(sph_system.real_bodies_);
    ObservedQuantityRecording<Vec2d> write_displacement_1("Displacement", disp_observer_contact_1);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    gate_normal_direction.exec();
    plate_corrected_configuration.exec();
    plate_average_curvature.exec();
    water_block_complex.updateConfiguration();
    plate_water_contact.updateConfiguration();
    constant_gravity.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 10;
    Real output_interval = end_time / 200.0;
    Real dt = 0.0;   // default acoustic time step sizes
    Real dt_s = 0.0; /**< Default acoustic time step sizes for solid. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_water_block_states.writeToFile(0);
    write_displacement_1.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();

            viscous_force_on_plate.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                pressure_relaxation.exec(dt);

                pressure_force_on_plate.exec();
                density_relaxation.exec(dt);

                /** Solid dynamics time stepping. */
                Real dt_s_sum = 0.0;
                average_velocity_and_acceleration.initialize_displacement_.exec();
                while (dt_s_sum < dt)
                {
                    dt_s = 0.5 * plate_time_step_size.exec();
                    if (dt - dt_s_sum < dt_s)
                        dt_s = dt - dt_s_sum;
                    plate_stress_relaxation_first.exec(dt_s);
                    plate_constraint.exec();
                    plate_position_damping.exec(dt_s);
                    plate_rotation_damping.exec(dt_s);
                    plate_constraint.exec();
                    plate_stress_relaxation_second.exec(dt_s);
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

            water_block.updateCellLinkedListWithParticleSort(100);

            if (GlobalStaticVariables::physical_time_ < gate_moving_time)
            {
                update_gate_position(GlobalStaticVariables::physical_time_);
                gate.updateCellLinkedList();
            }

            /** Update normal direction at elastic body surface. */
            plate_update_normal.exec();
            plate.updateCellLinkedList();
            plate_curvature_inner.updateConfiguration();
            plate_average_curvature.exec();
            plate_water_contact.updateConfiguration();

            water_block_complex.updateConfiguration();

            write_displacement_1.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        if (GlobalStaticVariables::physical_time_ < gate_moving_time)
            gate.setNewlyUpdated();
        write_water_block_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
};