/**
 * @file 	2d_eulerian_flow_around_cylinder_LG.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder coupling with the Laguerre Gauss kernel.
 * @details We consider a Eulerian flow passing by a cylinder in 2D.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 15.0;                        /**< Channel length. */
Real DH = 10.0;                        /**< Channel height. */
Real resolution_ref = 1.0 / 4.0;       /**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Vec2d cylinder_center(4.0, DH / 2.0);  /**< Location of the cylinder center. */
Real cylinder_radius = 1.0;            /**< Radius of the cylinder. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Define geometries and body shapes
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, DH + DH_sponge));
    water_block_shape.push_back(Vecd(DL, -DH_sponge));
    water_block_shape.push_back(Vecd(-DL_sponge, -DH_sponge));

    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon circle(cylinder_center, cylinder_radius, 100);
        subtract<MultiPolygonShape>(circle);
    }
};
class Cylinder : public MultiPolygonShape
{
  public:
    explicit Cylinder(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        multi_polygon_.addACircle(cylinder_center, cylinder_radius, 100, ShapeBooleanOps::add);
    }
};

class FarFieldBoundary : public fluid_dynamics::NonReflectiveBoundaryCorrection
{
  public:
    explicit FarFieldBoundary(BaseInnerRelation &inner_relation)
        : fluid_dynamics::NonReflectiveBoundaryCorrection(inner_relation)
    {
        rho_farfield_ = rho0_f;
        sound_speed_ = c_f;
        vel_farfield_ = Vecd(U_f, 0.0);
    };
    virtual ~FarFieldBoundary() {};
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    // Tag for run particle relaxation for the initial body fitted distribution.
    sph_system.setRunParticleRelaxation(false);
    // Tag for computation start with relaxed body fitted particles distribution.
    sph_system.setReloadParticles(false);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.getSPHAdaptation().resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
    water_block.defineComponentLevelSetShape("OuterBoundary");
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<BaseParticles, Reload>(water_block.getName())
        : water_block.generateParticles<BaseParticles, Lattice>();

    SolidBody cylinder(sph_system, makeShared<Cylinder>("Cylinder"));
    cylinder.defineAdaptationRatios(1.3, 2.0);
    cylinder.getSPHAdaptation().resetKernel<KernelTabulated<KernelLaguerreGauss>>(20);
    cylinder.defineBodyLevelSetShape();
    cylinder.defineMaterial<Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<BaseParticles, Reload>(cylinder.getName())
        : cylinder.generateParticles<BaseParticles, Lattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //	Note that the same relation should be defined only once.
    //----------------------------------------------------------------------
    InnerRelation water_block_inner(water_block);
    InnerRelation cylinder_inner(cylinder);
    ContactRelation water_block_contact(water_block, {&cylinder});
    ContactRelation cylinder_contact(cylinder, {&water_block});
    //----------------------------------------------------------------------
    // Combined relations built from basic relations
    // which is only used for update configuration.
    //----------------------------------------------------------------------
    ComplexRelation water_wall_complex(water_block_inner, water_block_contact);
    //----------------------------------------------------------------------
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        using namespace relax_dynamics;
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
        BodyStatesRecordingToVtp write_real_body_states(sph_system);
        ReloadParticleIO write_real_body_particle_reload_files({&cylinder, &water_block});
        RelaxationStepLevelSetCorrectionInner relaxation_step_inner(cylinder_inner);
        RelaxationStepLevelSetCorrectionComplex relaxation_step_complex(
            DynamicsArgs(water_block_inner, std::string("OuterBoundary")), water_block_contact);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        random_water_body_particles.exec(0.25);
        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_complex.SurfaceBounding().exec();
        write_real_body_states.writeToFile(0);

        int ite_p = 0;
        while (ite_p < 1000)
        {
            relaxation_step_inner.exec();
            relaxation_step_complex.exec();
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                write_real_body_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;

        write_real_body_particle_reload_files.writeToFile(0);

        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    InteractionWithUpdate<FreeSurfaceIndicationComplex> surface_indicator(water_block_inner, water_block_contact);
    InteractionDynamics<SmearedSurfaceIndication> smeared_surface(water_block_inner);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> cylinder_kernel_correction_matrix(cylinder_inner, cylinder_contact);
    InteractionWithUpdate<LinearGradientCorrectionMatrixComplex> water_block_kernel_correction_matrix(water_block_inner, water_block_contact);
    InteractionDynamics<KernelGradientCorrectionComplex> kernel_gradient_update(water_block_inner, water_block_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);

    InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfWithWallRiemann> pressure_relaxation(water_block_inner, water_block_contact);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfWithWallRiemann> density_relaxation(water_block_inner, water_block_contact);

    InteractionWithUpdate<fluid_dynamics::ViscousForceWithWall> viscous_force(water_block_inner, water_block_contact);
    SimpleDynamics<NormalDirectionFromBodyShape> water_block_normal_direction(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStep> get_fluid_time_step_size(water_block, 0.5);
    InteractionWithUpdate<FarFieldBoundary> variable_reset_in_boundary_condition(water_block_inner);
    //----------------------------------------------------------------------
    //	Compute the force exerted on solid body due to fluid pressure and viscosity
    //----------------------------------------------------------------------
    InteractionWithUpdate<solid_dynamics::ViscousForceFromFluid> viscous_force_from_fluid(cylinder_contact);
    InteractionWithUpdate<solid_dynamics::PressureForceFromFluid<decltype(density_relaxation)>> pressure_force_from_fluid(cylinder_contact);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_real_body_states(sph_system);
    write_real_body_states.addToWrite<int>(water_block, "Indicator");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<QuantitySummation<Vecd>>>
        write_total_viscous_force_from_fluid(cylinder, "ViscousForceFromFluid");
    ReducedQuantityRecording<QuantitySummation<Vecd>>
        write_total_pressure_force_from_fluid_body(cylinder, "PressureForceFromFluid");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    cylinder_normal_direction.exec();
    surface_indicator.exec();
    smeared_surface.exec();
    water_block_normal_direction.exec();
    variable_reset_in_boundary_condition.exec();
    cylinder_kernel_correction_matrix.exec();
    water_block_kernel_correction_matrix.exec();
    kernel_gradient_update.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 80.0;
    Real output_interval = 5.0; /**< time stamps for output. */
    //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_real_body_states.writeToFile(0);
    //----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real dt = get_fluid_time_step_size.exec();
            viscous_force.exec();
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);

            integration_time += dt;
            physical_time += dt;
            variable_reset_in_boundary_condition.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }

        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
        viscous_force_from_fluid.exec();
        pressure_force_from_fluid.exec();
        write_total_viscous_force_from_fluid.writeToFile(number_of_iterations);
        write_total_pressure_force_from_fluid_body.writeToFile(number_of_iterations);

        write_maximum_speed.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    write_total_viscous_force_from_fluid.testResult();

    return 0;
}
