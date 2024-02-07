/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D flow around rigid plate                       *
 *-----------------------------------------------------------------------------*
 * Test of 3D fluid-rigid shell interaction                                    *
 * where shell particles are close to each other                               *
 *-----------------------------------------------------------------------------*/

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real unit_scale = 0.001;               // mm to m
const Real fluid_diameter = 12 * unit_scale; /**<Inlet and outlet diameter. */
const Real fluid_radius = 0.5 * fluid_diameter;
const Real full_length = fluid_diameter * 7.5;           /**<Total length og fluid. */
const Real shell_thickness = 0.4 * unit_scale;           /**<Balloon thickness. */
const Real resolution_fluid = fluid_diameter / Real(15); /**< Global reference resolution. */
const Real resolution_shell = shell_thickness;           /**< Balloon resolution. */
const Real wall_thickness = 4 * resolution_fluid;        /**< Wall boundary thickness. */
const Real inflow_length = 10 * resolution_fluid;        /**< Inflow region. */

const Vecd plate_translation(fluid_diameter * 3.5, 0, 0);

const Vec3d emitter_halfsize(resolution_fluid * 2, fluid_radius, fluid_radius);
const Vec3d emitter_translation(resolution_fluid * 2, 0., 0.);
const Vec3d buffer_halfsize(inflow_length * 0.5, fluid_radius, fluid_radius);
const Vec3d buffer_translation(inflow_length * 0.5, 0., 0.);
const Vec3d disposer_halfsize(resolution_fluid * 2, fluid_radius * 1.1, fluid_radius * 1.1);
const Vec3d disposer_translation(full_length - disposer_halfsize[0], 0., 0.);

/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(-wall_thickness, -fluid_radius - wall_thickness, -fluid_radius - wall_thickness),
                                 Vec3d(full_length + wall_thickness, fluid_radius + shell_thickness, fluid_radius + shell_thickness));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1056.0;                           // blood density
const Real mu_f = 3.5e-3;                             // blood viscosity
const Real Re = 100;                                  /**< Reynolds number. */
const Real U_f = Re * mu_f / rho0_f / fluid_diameter; /**< Average velocity at throat. */
const Real U_max = 2.0 * U_f;
const Real c_f = 10.0 * U_max; /**< Speed of sound. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
const std::string path_to_fluid_file = "./input/fluid.stl";
const std::string path_to_plate_srf_file = "./input/plate_outer_srf_0_8.stl";
class FluidBlock : public ComplexShape
{
  public:
    explicit FluidBlock(const std::string &shape_name)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_fluid_file, Vecd::Zero(), unit_scale);
        subtract<TriangleMeshShapeSTL>(path_to_plate_srf_file, plate_translation, unit_scale);
    }
};
const std::string path_to_wall_file = "./input/wall_3_2.stl";
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name)
        : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_wall_file, Vecd::Zero(), unit_scale);
    }
};
//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_;
    AlignedBoxShape &aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType &boundary_condition)
        : u_ref_(U_f),
          aligned_box_(boundary_condition.getAlignedBox()),
          halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd &position, Vecd &velocity)
    {
        Vecd target_velocity = Vecd::Zero();
        Real radius2 = position[1] * position[1] + position[2] * position[2];
        target_velocity[0] = 2.0 * u_ref_ * (1.0 - radius2 / halfsize_[1] / halfsize_[1]);

        return target_velocity;
    }
};

int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem with global controls.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_fluid);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidBlock>("fluid"));
    fluid_block.defineBodyLevelSetShape()->correctLevelSetSign()->cleanLevelSet(0);
    fluid_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    fluid_block.generateParticles<ParticleGeneratorLattice>();

    SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("wall_3_2"));
    wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>(); // dummy material parameters
    wall_boundary.generateParticles<ParticleGeneratorReload>(wall_boundary.getName());

    SolidBody shell(sph_system, makeShared<DefaultShape>("plate_0_4"));
    shell.defineAdaptation<SPHAdaptation>(1.15, resolution_fluid / resolution_shell);
    shell.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1.0, 1.0, 0.0);
    shell.generateParticles<ParticleGeneratorReload>(shell.getName());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    // Must construct ShellCurvature before ShellContactRelation
    InnerRelation fluid_inner(fluid_block);
    InnerRelation shell_inner(shell);
    ContactRelation fluid_wall_contact(fluid_block, {&wall_boundary});
    ContactRelationToShell fluid_shell_contact(fluid_block, {&shell}, {true});
    ComplexRelation fluid_block_complex(fluid_inner, {&fluid_wall_contact, &fluid_shell_contact});
    ShellInnerRelationWithContactKernel shell_curvature_inner(shell, fluid_block);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Algorithm for fluid dynamics. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(fluid_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(fluid_block);
    InteractionWithUpdate<fluid_dynamics::BaseDensitySummationComplex<Inner<FreeStream>, Contact<>, Contact<>>> update_fluid_density_by_summation(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration1stHalf<Inner<>, Contact<Wall>, Contact<Wall>>, AcousticRiemannSolver, NoKernelCorrection>> fluid_pressure_relaxation(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    Dynamics1Level<ComplexInteraction<fluid_dynamics::Integration2ndHalf<Inner<>, Contact<Wall>, Contact<Wall>>, NoRiemannSolver>> fluid_density_relaxation(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::ViscousForce<Inner<>, Contact<Wall>, Contact<Wall>>, fluid_dynamics::FixedViscosity>> viscous_acceleration(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    InteractionWithUpdate<ComplexInteraction<FreeSurfaceIndication<Inner<SpatialTemporal>, Contact<>, Contact<>>>> inlet_outlet_surface_particle_indicator(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    InteractionWithUpdate<ComplexInteraction<fluid_dynamics::TransportVelocityCorrection<Inner<SingleResolution>, Contact<Boundary>, Contact<Boundary>>, NoKernelCorrection, BulkParticles>> transport_velocity_correction(fluid_inner, fluid_wall_contact, fluid_shell_contact);
    /** Algorithm for in-/outlet. */
    BodyAlignedBoxByParticle emitter(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec3d(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);
    BodyAlignedBoxByCell buffer(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec3d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> emitter_buffer_inflow_condition(buffer);
    BodyAlignedBoxByCell disposer(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec3d(disposer_translation)), disposer_halfsize));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_outflow_deletion(disposer, xAxis);
    /** Algorithm for solid dynamics. */
    SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary);
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_average_curvature(shell_curvature_inner);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    fluid_block.addBodyStateForRecording<Real>("Pressure");
    fluid_block.addBodyStateForRecording<int>("Indicator");
    shell.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    shell.addBodyStateForRecording<Real>("Average2ndPrincipleCurvature");
    BodyStatesRecordingToVtp write_real_body_states(sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    wall_boundary_normal_direction.exec();
    shell_average_curvature.exec();
    fluid_block_complex.updateConfiguration();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = sph_system.RestartStep();
    int screen_output_interval = 10;
    Real end_time = 1.0;
    Real output_interval = end_time / 100.0; /**< Time stamps for output of body states. */
    Real dt = 0.0;                           /**< Default acoustic time step sizes. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile();
    //----------------------------------------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------------------------------------
    const double Dt_ref = fluid_advection_time_step.exec();
    const double dt_ref = fluid_acoustic_time_step.exec();
    auto run_simulation = [&]()
    {
        std::cout << "Simulation starts here" << std::endl;
        while (GlobalStaticVariables::physical_time_ < end_time)
        {
            Real integration_time = 0.0;
            /** Integrate time (loop) until the next output time. */
            while (integration_time < output_interval)
            {
                Real Dt = fluid_advection_time_step.exec();
                if (Dt < Dt_ref / 20)
                {
                    std::cout << "Dt = " << Dt << ", Dt_ref = " << Dt_ref << std::endl;
                    std::cout << "Advective time step decreased too much!" << std::endl;
                    throw std::runtime_error("Advective time step decreased too much!");
                }

                inlet_outlet_surface_particle_indicator.exec();
                update_fluid_density_by_summation.exec();
                viscous_acceleration.exec();
                transport_velocity_correction.exec();

                /** Dynamics including pressure relaxation. */
                Real relaxation_time = 0.0;
                while (relaxation_time < Dt)
                {
                    Real dt_temp = fluid_acoustic_time_step.exec();
                    if (dt_temp < dt_ref / 20)
                    {
                        std::cout << "dt = " << dt_temp << ", dt_ref = " << dt_ref << std::endl;
                        std::cout << "Acoustic time step decreased too much!" << std::endl;
                        throw std::runtime_error("Acoustic time step decreased too much!");
                    }
                    dt = SMIN(dt_temp, Dt - relaxation_time);
                    fluid_pressure_relaxation.exec(dt);
                    emitter_buffer_inflow_condition.exec();
                    fluid_density_relaxation.exec(dt);

                    relaxation_time += dt;
                    integration_time += dt;
                    GlobalStaticVariables::physical_time_ += dt;
                }

                if (number_of_iterations % screen_output_interval == 0)
                {
                    std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                              << GlobalStaticVariables::physical_time_
                              << "	Dt = " << Dt << "	dt = " << dt;
                    std::cout << "\n";
                }
                number_of_iterations++;

                /** inflow injection*/
                emitter_inflow_injection.exec();
                disposer_outflow_deletion.exec();

                /** Update cell linked list and configuration. */
                fluid_block.updateCellLinkedListWithParticleSort(100);
                fluid_block_complex.updateConfiguration();
            }

            TickCount t2 = TickCount::now();
            write_real_body_states.writeToFile();
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
        fluid_block.setNewlyUpdated();
        shell.setNewlyUpdated();
        write_real_body_states.writeToFile(1e8);
    }
    return 0;
}
