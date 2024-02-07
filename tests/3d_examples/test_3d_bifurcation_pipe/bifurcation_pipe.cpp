/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: 3D bifurcation pipe                       *
 *-----------------------------------------------------------------------------*
 * Test of 3D fluid-rigid shell interaction, where shell pipe has a sharp angle*
 *-----------------------------------------------------------------------------*/

#include "sphinxsys.h"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
const Real unit_scale = 1.0;                 // m
const Real diameter = 6.0 * unit_scale;      /**<Inlet diameter. */
const Real diameter_side = 4.0 * unit_scale; /**<Side pipe diameter. */
const Real angle = 10 * M_PI / 180.0;        /**<Angle between main and side pipe>*/
const Real radius = 0.5 * diameter;
const Real radius_side = 0.5 * diameter_side;
const Real wall_extension = 2.4 * unit_scale; /**<Extension of wall. */
const Real shell_thickness = 0.2 * unit_scale;
const Real resolution_ref = diameter / Real(15); /**< Global reference resolution. */
const Real resolution_shell = shell_thickness;
const Real inflow_length = 10 * resolution_ref; /**< Inflow region. */

const Vecd inlet_center(0, 0, 0);
const Vecd outlet_center_1(60, 0, 0);
const Vecd outlet_center_2(59.18, 0, 9.38);

/**
 * @brief Geometry parameters for boundary condition.
 */
const Vec3d emitter_halfsize(resolution_ref * 2, radius, radius);
const Vec3d emitter_translation = inlet_center + Vecd(resolution_ref * 2, 0., 0.);
const Vec3d buffer_halfsize(inflow_length * 0.5, radius, radius);
const Vec3d buffer_translation = inlet_center + Vecd(inflow_length * 0.5, 0., 0.);
// main pipe disposer
const Vec3d disposer_halfsize_1(resolution_ref * 2, radius * 1.1, radius * 1.1);
const Vec3d disposer_translation_1 = outlet_center_1 - Vec3d(resolution_ref * 2, 0, 0);
// side pipe disposer
const Vec3d disposer_halfsize_2(resolution_ref * 2, radius_side * 1.1, radius_side * 1.1);
const Vec3d disposer_translation_2 = outlet_center_2 - resolution_ref * 2 * Vec3d(cos(angle), 0, sin(angle));

/** Domain bounds of the system. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(inlet_center[0] - wall_extension, -radius - shell_thickness, -radius - shell_thickness),
                                 Vec3d(outlet_center_1[0] + wall_extension, radius + shell_thickness, 12.0));
//----------------------------------------------------------------------
//	Global parameters on the fluid properties
//----------------------------------------------------------------------
const Real rho0_f = 1.0;
const Real U_f = 1.0;
const Real Re = 100;                            /**< Reynolds number. */
const Real mu_f = rho0_f * U_f * diameter / Re; /**< Dynamics viscosity. */
const Real U_max = 3 * U_f;
const Real c_f = 10.0 * U_max; /**< Speed of sound. */
//----------------------------------------------------------------------
//	define geometry of SPH bodies
//----------------------------------------------------------------------
const std::string path_to_fluid_file = "./input/fluid.stl";
class FluidShape : public ComplexShape
{
  public:
    explicit FluidShape(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_fluid_file, Vecd::Zero(), 1.0);
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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    /** Tag for computation start with relaxed body fitted particles distribution. */
    sph_system.setReloadParticles(true);
    sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody fluid_block(sph_system, makeShared<FluidShape>("fluid_0_4"));
    fluid_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    // Using relaxed particle distribution if needed
    fluid_block.generateParticles<ParticleGeneratorReload>(fluid_block.getName());

    SolidBody wall_boundary(sph_system, makeShared<DefaultShape>("wall"));
    wall_boundary.defineAdaptationRatios(1.15, resolution_ref / resolution_shell);
    wall_boundary.defineParticlesAndMaterial<ShellParticles, SaintVenantKirchhoffSolid>(1.0, 1.0, 0.0); // dummy material parameters
    wall_boundary.generateParticles<ParticleGeneratorReload>(wall_boundary.getName());
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //----------------------------------------------------------------------
    // Must construct ShellCurvature before ShellContactRelation
    InnerRelation fluid_block_inner(fluid_block);
    InnerRelation wall_boundary_inner(wall_boundary);
    ContactRelationToShell fluid_wall_contact(fluid_block, {&wall_boundary}, {false});
    ComplexRelation fluid_block_complex(fluid_block_inner, fluid_wall_contact);
    ShellInnerRelationWithContactKernel shell_curvature_inner(wall_boundary, fluid_block);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** time-space method to detect surface particles. (Important for
     * DensitySummationFreeSurfaceComplex work correctly.)*/
    InteractionWithUpdate<SpatialTemporalFreeSurfaceIndicationComplex>
        inlet_outlet_surface_particle_indicator(fluid_block_inner, fluid_wall_contact);
    /** Evaluation of density by summation approach. */
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeStreamComplex> update_density_by_summation(fluid_block_inner, fluid_wall_contact);
    /** Time step size without considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(fluid_block, U_max);
    /** Time step size with considering sound wave speed. */
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(fluid_block);
    /** Pressure relaxation using verlet time stepping. */
    /** Here, we do not use Riemann solver for pressure as the flow is viscous. */
    Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann> pressure_relaxation(fluid_block_inner, fluid_wall_contact);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallNoRiemann> density_relaxation(fluid_block_inner, fluid_wall_contact);
    /** viscous acceleration and transport velocity correction can be combined because they are independent dynamics. */
    InteractionWithUpdate<fluid_dynamics::TransportVelocityCorrectionComplex<BulkParticles>> transport_correction(fluid_block_inner, fluid_wall_contact);
    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_acceleration(fluid_block_inner, fluid_wall_contact);
    /** Inflow boundary condition. */
    BodyAlignedBoxByParticle emitter(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vecd(emitter_translation)), emitter_halfsize));
    SimpleDynamics<fluid_dynamics::EmitterInflowInjection> emitter_inflow_injection(emitter, 10, 0);
    BodyAlignedBoxByCell inflow_buffer(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vec3d(buffer_translation)), buffer_halfsize));
    SimpleDynamics<fluid_dynamics::InflowVelocityCondition<InflowVelocity>> parabolic_inflow(inflow_buffer);

    BodyAlignedBoxByCell disposer_main(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Vecd(disposer_translation_1)), disposer_halfsize_1));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_main_outflow_deletion(disposer_main, xAxis);

    BodyAlignedBoxByCell disposer_side(
        fluid_block, makeShared<AlignedBoxShape>(Transform(Rotation3d(-angle, Vec3d::UnitY()), Vecd(disposer_translation_2)), disposer_halfsize_2));
    SimpleDynamics<fluid_dynamics::DisposerOutflowDeletion> disposer_side_outflow_deletion(disposer_side, xAxis);
    // wall curvature
    SimpleDynamics<thin_structure_dynamics::AverageShellCurvature> shell_curvature(shell_curvature_inner);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    fluid_block.addBodyStateForRecording<Real>("Pressure");
    fluid_block.addBodyStateForRecording<int>("Indicator");
    wall_boundary.addBodyStateForRecording<Real>("Average1stPrincipleCurvature");
    wall_boundary.addBodyStateForRecording<Real>("Average2ndPrincipleCurvature");
    BodyStatesRecordingToVtp write_real_body_states(sph_system.real_bodies_);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    /** initialize cell linked lists for all bodies. */
    sph_system.initializeSystemCellLinkedLists();
    /** initialize configurations for all bodies. */
    sph_system.initializeSystemConfigurations();
    /** initial curvature*/
    shell_curvature.exec();
    fluid_block_complex.updateConfiguration();
    //----------------------------------------------------------------------
    //	Setup computing and initial conditions.
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 10;
    Real end_time = 10.0;
    Real output_interval = end_time / 200.0;
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
            Real Dt = get_fluid_advection_time_step_size.exec();
            inlet_outlet_surface_particle_indicator.exec();
            update_density_by_summation.exec();
            viscous_acceleration.exec();
            transport_correction.exec();

            size_t inner_ite_dt = 0;
            Real relaxation_time = 0.0;
            while (relaxation_time < Dt)
            {
                Real dt = SMIN(get_fluid_time_step_size.exec(), Dt - relaxation_time);
                /** Fluid pressure relaxation */
                pressure_relaxation.exec(dt);
                /** velocity */
                parabolic_inflow.exec();
                /** Fluid density relaxation */
                density_relaxation.exec(dt);

                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;

                inner_ite_dt++;
            }

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << GlobalStaticVariables::physical_time_
                          << "	Dt = " << Dt << "	Dt / dt = " << inner_ite_dt << "\n";
            }
            number_of_iterations++;

            /** fluid block configuration and periodic condition. */
            emitter_inflow_injection.exec();
            disposer_main_outflow_deletion.exec();
            disposer_side_outflow_deletion.exec();

            fluid_block.updateCellLinkedListWithParticleSort(100);
            fluid_block_complex.updateConfiguration();
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
