/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D dambreak example                        *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// time
const Real gate_moving_time = 0.1;
const Real contact_time = 0.25;
const Real end_time = 0.5;
// general parameters for geometry
const Real t = 4e-3;                                        // gate thickness
const Real plate_height = 0.09;                             // plate height
const Real plate_width = 0.1995;                            // plate width
const Real LH = 0.2;                                        // liquid initial height
const Real LL = 0.2;                                        // liquid initial length
const Real LW = 0.2;                                        // liquid initial width
const Real DH = 0.4;                                        // tank height
const Real DL = 0.8;                                        // tank length
const Real DW = 0.2;                                        // tank width
const Real resolution_shell = t;                            // shell particle spacing
const Real resolution_ref = 2 * resolution_shell;           // system particle spacing
const Real BW = resolution_ref * 4;                         // boundary width
const Real plate_x_pos = DL - 0.2 + 0.5 * resolution_shell; // center x coordinate of plate

const Real marker_h = 0.0875; // height of marker
const std::vector<Vec3d> observer_position_1 = {Vec3d(plate_x_pos, marker_h, (DW - plate_width) * 0.5)};
const std::vector<Vec3d> observer_position_2 = {Vec3d(plate_x_pos, marker_h, 0.5 * DW)};

const BoundingBoxd system_domain_bounds(Vec3d(-BW, -BW, -BW), Vec3d(DL + BW, DH + BW, DW + BW));

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

//	define the water block shape
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
        Transform translation_water(halfsize_water);
        add<GeometricShapeBox>(Transform(translation_water), halfsize_water);
    }
};
//	define the static solid wall boundary shape
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d halfsize_outer(0.5 * DL + BW, 0.5 * DH + BW, 0.5 * DW + BW);
        Vec3d halfsize_inner(0.5 * DL, 0.5 * DH, 0.5 * DW);
        Transform translation_wall(halfsize_inner);
        add<GeometricShapeBox>(Transform(translation_wall), halfsize_outer);
        subtract<GeometricShapeBox>(Transform(translation_wall), halfsize_inner);

        Vec3d halfsize_plate(0.5 * resolution_ref, 0.5 * plate_width, 0.5 * (plate_height + BW));
        Transform translation_plate(halfsize_plate + Vec3d(plate_x_pos, -BW, (DW - plate_width) * 0.5));
        subtract<GeometricShapeBox>(Transform(translation_plate), halfsize_plate);
    }
};
//	define the rigid solid gate shape
class MovingGate : public ComplexShape
{
  public:
    explicit MovingGate(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vec3d halfsize_gate(0.5 * BW, 0.5 * DH, 0.5 * DW);
        Transform translation_gate(Vec3d(LL, 0, 0) + halfsize_gate);
        add<GeometricShapeBox>(Transform(translation_gate), halfsize_gate);
    }
};
//	define the elastic plate shape
namespace SPH
{
class Plate;
template <>
class ParticleGenerator<SurfaceParticles, Plate> : public ParticleGenerator<SurfaceParticles>
{
  public:
    explicit ParticleGenerator(SPHBody &sph_body, SurfaceParticles &surface_particles)
        : ParticleGenerator<SurfaceParticles>(sph_body, surface_particles) {};
    void prepareGeometricData() override
    {
        Real y = -BW + 0.5 * resolution_shell;
        while (y < plate_height)
        {
            Real z = (DW - plate_width + resolution_shell) * 0.5;
            while (z < 0.5 * (DW + plate_width))
            {
                addPositionAndVolumetricMeasure(Vec3d(plate_x_pos, y, z), resolution_shell * resolution_shell);
                addSurfaceProperties(Vec3d(1, 0, 0), t);
                z += resolution_shell;
            }
            y += resolution_shell;
        }
    }
};
} // namespace SPH

// Define boundary geometry

class BoundaryGeometry : public BodyPartByParticle
{
  public:
    BoundaryGeometry(SPHBody &body) : BodyPartByParticle(body)
    {
        TaggingParticleMethod tagging_particle_method = std::bind(&BoundaryGeometry::tagManually, this, _1);
        tagParticles(tagging_particle_method);
    };

  private:
    bool tagManually(size_t index_i)
    {
        return base_particles_.ParticlePositions()[index_i].y() <= 0;
    };
};

class GateMotionConstraint : public MotionConstraint<SPHBody>
{
  public:
    GateMotionConstraint(SPHBody &body)
        : MotionConstraint<SPHBody>(body),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")) {};
    virtual ~GateMotionConstraint() {};
    void update(size_t index_i, Real dt)
    {
        Real run_time = *physical_time_;
        Real h_g = -285.115 * run_time * run_time * run_time + 72.305 * run_time * run_time + 0.1463 * run_time;
        pos_[index_i][1] = pos0_[index_i][1] + h_g;
    };

  protected:
    Real *physical_time_;
};

// the main program with commandline options
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up an SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    sph_system.setGenerateRegressionData(false);
    //----------------------------------------------------------------------
    //	Creating bodies with corresponding materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary.defineMaterial<Solid>();
    wall_boundary.generateParticles<BaseParticles, Lattice>();

    SolidBody gate(sph_system, makeShared<MovingGate>("Gate"));
    gate.defineMaterial<Solid>();
    gate.generateParticles<BaseParticles, Lattice>();

    SolidBody plate(sph_system, makeShared<DefaultShape>("Plate"));
    plate.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    plate.defineMaterial<SaintVenantKirchhoffSolid>(rho0_s, youngs_modulus, poisson_ratio);
    plate.generateParticles<SurfaceParticles, Plate>();

    ObserverBody disp_observer_1(sph_system, "Observer1");
    disp_observer_1.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    ObserverBody disp_observer_2(sph_system, "Observer2");
    disp_observer_2.defineAdaptation<SPHAdaptation>(1.15, resolution_ref / resolution_shell);
    disp_observer_1.generateParticles<ObserverParticles>(observer_position_1);
    disp_observer_2.generateParticles<ObserverParticles>(observer_position_2);
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
    // shell normal should point from fluid to shell
    // normal corrector set to false if shell normal is already pointing from fluid to shell
    ContactRelationFromShellToFluid water_plate_contact(water_block, {&plate}, {false});
    ContactRelationFromFluidToShell plate_water_contact(plate, {&water_block}, {false});
    ShellInnerRelationWithContactKernel plate_curvature_inner(plate, water_block);
    ContactRelation disp_observer_contact_1(disp_observer_1, {&plate});
    ContactRelation disp_observer_contact_2(disp_observer_2, {&plate});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block_inner, {&water_wall_contact, &water_plate_contact});
    //----------------------------------------------------------------------
    //	Define the numerical methods used in the simulation.
    //	Note that there may be data dependence on the sequence of constructions.
    //----------------------------------------------------------------------
    // solid
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<NormalDirectionFromBodyShape> gate_normal_direction(gate);
    SimpleDynamics<GateMotionConstraint> update_gate_position(gate);
    // Shell
    InteractionDynamics<thin_structure_dynamics::ShellCorrectConfiguration> plate_corrected_configuration(plate_inner);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationFirstHalf> plate_stress_relaxation_first(plate_inner, 3, true);
    Dynamics1Level<thin_structure_dynamics::ShellStressRelaxationSecondHalf> plate_stress_relaxation_second(plate_inner);
    ReduceDynamics<thin_structure_dynamics::ShellAcousticTimeStepSize> plate_time_step_size(plate);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> plate_average_curvature(plate_curvature_inner);
    SimpleDynamics<thin_structure_dynamics::UpdateShellNormalDirection> plate_update_normal(plate);
    /** constraint and damping */
    BoundaryGeometry plate_boundary_geometry(plate);
    SimpleDynamics<thin_structure_dynamics::ConstrainShellBodyRegion> plate_constraint(plate_boundary_geometry);
    // fluid
    Gravity gravity(Vec3d(0.0, -gravity_g, 0.0));
    SimpleDynamics<GravityForce<Gravity>> constant_gravity(water_block, gravity);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>>
        pressure_relaxation(water_block_inner, water_wall_contact, water_plate_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver>>
        density_relaxation(water_block_inner, water_wall_contact, water_plate_contact);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeSurface>, Contact<>, Contact<>>>
        update_density_by_summation(water_block_inner, water_wall_contact, water_plate_contact);
    ReduceDynamics<fluid_dynamics::AdvectionViscousTimeStep> get_fluid_advection_time_step_size(water_block, U_f);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::ViscousForce<Inner<>, Contact<Wall>, Contact<Wall>>,
                                             fluid_dynamics::FixedViscosity, NoKernelCorrection>>
        viscous_acceleration(water_block_inner, water_wall_contact, water_plate_contact);
    // FSI
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_on_plate(plate_water_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_on_plate(plate_water_contact);
    solid_dynamics::AverageVelocityAndAcceleration average_velocity_and_acceleration(plate);
    //----------------------------------------------------------------------
    //	Define the configuration related particles dynamics.
    //----------------------------------------------------------------------
    ParticleSorting particle_sorting(water_block);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations, observations
    //	and regression tests of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<Real>(water_block, "Pressure");
    write_real_body_states.addToWrite<Vec3d>(gate, "NormalDirection");
    write_real_body_states.addToWrite<Real>(plate, "Average1stPrincipleCurvature");
    write_real_body_states.addToWrite<Real>(plate, "Average2ndPrincipleCurvature");
    write_real_body_states.addToWrite<Vec3d>(plate, "PressureForceFromFluid");
    write_real_body_states.addToWrite<Vec3d>(wall_boundary, "NormalDirection");
    write_real_body_states.addDerivedVariableRecording<SimpleDynamics<Displacement>>(plate);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_displacement_1("Displacement", disp_observer_contact_1);
    RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>> write_displacement_2("Displacement", disp_observer_contact_2);
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
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
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
    write_real_body_states.writeToFile(0);
    write_displacement_1.writeToFile(0);
    write_displacement_2.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();

            if (physical_time > contact_time)
                viscous_force_on_plate.exec();

            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                pressure_relaxation.exec(dt);
                if (physical_time > contact_time)
                    pressure_force_on_plate.exec();
                density_relaxation.exec(dt);

                if (physical_time > contact_time)
                {
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
                        plate_stress_relaxation_second.exec(dt_s);
                        dt_s_sum += dt_s;
                    }
                    average_velocity_and_acceleration.update_averages_.exec(dt);
                }

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

            if (number_of_iterations % 100 == 0 && number_of_iterations != 1)
            {
                particle_sorting.exec();
            }
            water_block.updateCellLinkedList();

            if (physical_time < gate_moving_time)
            {
                update_gate_position.exec();
                gate.updateCellLinkedList();
            }
            if (physical_time > contact_time)
            {
                /** Update normal direction at elastic body surface. */
                plate_update_normal.exec();
                plate.updateCellLinkedList();
                plate_curvature_inner.updateConfiguration();
                plate_average_curvature.exec();
                plate_water_contact.updateConfiguration();
            }
            water_block_complex.updateConfiguration();

            write_displacement_1.writeToFile(number_of_iterations);
            write_displacement_2.writeToFile(number_of_iterations);
        }
        TickCount t2 = TickCount::now();
        if (physical_time <= gate_moving_time)
            gate.setNewlyUpdated();
        if (physical_time <= contact_time)
            plate.setNewlyUpdated();
        write_real_body_states.writeToFile();
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    if (sph_system.GenerateRegressionData())
    {
        write_displacement_1.generateDataBase(1.0e-3);
        write_displacement_2.generateDataBase(1.0e-3);
    }
    else
    {
        write_displacement_1.testResult();
        write_displacement_2.testResult();
    }

    return 0;
}