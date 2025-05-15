/**
 * @file 	2d_eulerian_flow_around_elliptic_cylinder.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a elliptic cylinder.
 * @details We consider a Eulerian flow passing by a cylinder in 2D.
 * @author 	Bo Zhang and Xiangyu Hu
 */
#include "general_eulerian_fluid_dynamics.hpp" // eulerian classes for weakly compressible fluid only.
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 64;                          /**< Channel length. */
Real DH = 42;                          /**< Channel height. */
Real resolution_ref = 1.0 / 10.0;      /**< Initial reference particle spacing. */
Real DL_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Real DH_sponge = resolution_ref * 2.0; /**< Sponge region to impose inflow condition. */
Vec2d cylinder_center(21, 21);         /**< Location of the cylinder center. */
Real cylinder_radius = 1.0;            /**< Radius of the cylinder. */
Real aspect_ratio = 0.6;
Real shape_resolution = 500;
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Create elliptic shape
//----------------------------------------------------------------------
std::vector<Vecd> CreateEllipticShape(Real AR, Real sample_point)
{
    Real a = cylinder_radius;
    Real b = AR * cylinder_radius;

    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(1, 0) + cylinder_center);

    for (int n = 1; n < sample_point; ++n)
    {
        Real theta = 2 * PI * n / sample_point;
        Real x = a * cos(theta);
        Real y = b * sin(theta);
        pnts.push_back(Vecd(x, y) + cylinder_center);
    }
    pnts.push_back(Vecd(1, 0) + cylinder_center);
    return pnts;
}
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
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon ellptic_shape(CreateEllipticShape(aspect_ratio, shape_resolution));
        subtract<MultiPolygonShape>(ellptic_shape);
    }
};
class Cylinder : public MultiPolygonShape
{
public:
    explicit Cylinder(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> ellptic_shape = CreateEllipticShape(aspect_ratio, shape_resolution);
        multi_polygon_.addAPolygon(ellptic_shape, ShapeBooleanOps::add);
    }
};

class FarFieldBoundary : public fluid_dynamics::NonReflectiveBoundaryCorrection
{
public:
    explicit FarFieldBoundary(BaseInnerRelation& inner_relation)
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
    BoundingBox system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    // Tag for run particle relaxation for the initial body fitted distribution.
    sph_system.setRunParticleRelaxation(false);
    // Tag for computation start with relaxed body fitted particles distribution.
    sph_system.setReloadParticles(true);
    // Handle command line arguments and override the tags for particle relaxation and reload.
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineAdaptationRatios(0.95, 1.0);
    water_block.defineComponentLevelSetShape("OuterBoundary")->writeLevelSet(io_environment);
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f, mu_f);
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? water_block.generateParticles<ParticleGeneratorReload>(io_environment, water_block.getName())
        : water_block.generateParticles<ParticleGeneratorLattice>();
    water_block.addBodyStateForRecording<int>("Indicator");

    SolidBody cylinder(sph_system, makeShared<Cylinder>("Cylinder"));
    cylinder.defineAdaptationRatios(0.95, 2.0);
    cylinder.defineBodyLevelSetShape()->writeLevelSet(io_environment);
    cylinder.defineParticlesAndMaterial<SolidParticles, Solid>();
    (!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
        ? cylinder.generateParticles<ParticleGeneratorReload>(io_environment, cylinder.getName())
        : cylinder.generateParticles<ParticleGeneratorLattice>();
    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //	Note that the same relation should be defined only once.
    //----------------------------------------------------------------------
    ComplexRelation water_block_complex(water_block, {&cylinder});
    ComplexRelation water_block_complex_correction(water_block, {&cylinder});
    ContactRelation cylinder_contact(cylinder, {&water_block});
    InnerRelation cylinder_inner(cylinder);
    //----------------------------------------------------------------------
    //	Run particle relaxation for body-fitted distribution if chosen.
    //----------------------------------------------------------------------
    if (sph_system.RunParticleRelaxation())
    {
        InnerRelation water_block_inner(water_block);
        InnerRelation cylinder_inner(cylinder); // extra body topology only for particle relaxation
        //----------------------------------------------------------------------
        //	Methods used for particle relaxation.
        //----------------------------------------------------------------------
        SimpleDynamics<RandomizeParticlePosition> random_inserted_body_particles(cylinder);
        SimpleDynamics<RandomizeParticlePosition> random_water_body_particles(water_block);
        InteractionWithUpdate<KernelCorrectionMatrixInnerWithLevelSet> kernel_correction_inner(cylinder_inner);
        InteractionWithUpdate<KernelCorrectionMatrixComplexWithLevelSet> kernel_correction_complex(water_block_complex, "OuterBoundary");
        relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_step_inner(cylinder_inner, true);
        relax_dynamics::RelaxationStepComplexImplicit<CorrectionMatrixRelaxation> relaxation_step_complex(water_block_complex, "OuterBoundary", true);
        SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_water_block_kinetic_energy(water_block_inner);
        SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_cylinder_kinetic_energy(cylinder_inner);
        ReduceDynamics<Average<QuantitySummation<Real>>> calculate_water_block_average_kinetic_energy(water_block, "ParticleKineticEnergy");
        ReduceDynamics<Average<QuantitySummation<Real>>> calculate_cylinder_average_kinetic_energy(cylinder, "ParticleKineticEnergy");
        ReduceDynamics<QuantityMaximum<Real>> calculate_water_block_maximum_kinetic_energy(water_block, "ParticleKineticEnergy");
        ReduceDynamics<QuantityMaximum<Real>> calculate_cylinder_maximum_kinetic_energy(cylinder, "ParticleKineticEnergy");
        cylinder.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
        BodyStatesRecordingToVtp write_real_body_states(io_environment, sph_system.real_bodies_);
        ReloadParticleIO write_real_body_particle_reload_files(io_environment, sph_system.real_bodies_);
        //----------------------------------------------------------------------
        //	Particle relaxation starts here.
        //----------------------------------------------------------------------
        random_inserted_body_particles.exec(0.25);
        random_water_body_particles.exec(0.25);
        sph_system.initializeSystemCellLinkedLists();
        sph_system.initializeSystemConfigurations();
        relaxation_step_inner.SurfaceBounding().exec();
        relaxation_step_complex.SurfaceBounding().exec();
        write_real_body_states.writeToFile(0);

        Real water_block_average_energy = 100.0;
        Real water_block_maximum_energy = 100.0;
        Real cylinder_average_energy = 100.0;
        Real cylinder_maximum_energy = 100.0;
        Real last_water_block_maximum_energy = 100.0;
        Real dt = 1;

        int ite_p = 0;
        while (ite_p < 10000)
        {
            kernel_correction_inner.exec();
            relaxation_step_inner.exec();
            kernel_correction_complex.exec();
            relaxation_step_complex.exec(dt);
            ite_p += 1;
            if (ite_p % 200 == 0)
            {
                update_water_block_kinetic_energy.exec();
                update_cylinder_kinetic_energy.exec();
                water_block_maximum_energy = calculate_water_block_maximum_kinetic_energy.exec();
                water_block_average_energy = calculate_water_block_average_kinetic_energy.exec();
                cylinder_maximum_energy = calculate_cylinder_maximum_kinetic_energy.exec();
                cylinder_average_energy = calculate_cylinder_average_kinetic_energy.exec();
                std::cout << std::fixed << std::setprecision(9) << "Relaxation steps N = " << ite_p << "\n";
                std::cout << "WaterBody: " << " Maximum: " << water_block_maximum_energy << " Average: " << water_block_average_energy << std::endl;
                std::cout << "Cylinder: " <<  " Maximum: " << cylinder_maximum_energy << " Average: " << cylinder_average_energy << std::endl;

                if (water_block_maximum_energy > last_water_block_maximum_energy)
                {
                    dt = 1.0 * dt;
                }
                else if (water_block_maximum_energy < last_water_block_maximum_energy)
                {
                    dt = 1.0 * dt;
                }
                last_water_block_maximum_energy = water_block_maximum_energy;
                std::cout << "dt ratio is " << dt << std::endl;

                write_real_body_states.writeToFile(ite_p);
            }
        }
        std::cout << "The physics relaxation process finish !" << std::endl;

        std::cout << "Maximum residual: "
                  << " cylinder: " << calculate_cylinder_maximum_kinetic_energy.exec()
                  << " water_block: " << calculate_water_block_maximum_kinetic_energy.exec() << std::endl;

        std::cout << "Average residual: "
                  << " cylinder: " << calculate_cylinder_average_kinetic_energy.exec()
                  << " water_block: " << calculate_water_block_average_kinetic_energy.exec() << std::endl;

        write_real_body_particle_reload_files.writeToFile(0);

        return 0;
    }
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------    
    InteractionWithUpdate<fluid_dynamics::ICEIntegration1stHalfHLLERiemannWithWall> pressure_relaxation(water_block_complex_correction);
    InteractionWithUpdate<fluid_dynamics::ICEIntegration2ndHalfHLLERiemannWithWall> density_relaxation(water_block_complex_correction);
    InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_matrix_correction(water_block_complex_correction);
    InteractionWithUpdate<KernelCorrectionMatrixInnerWithLevelSet> kernel_correction_matrix_cylinder(cylinder_inner);
    InteractionWithUpdate<KernelCorrectionMatrixComplex> kernel_correction_matrix(water_block_complex);
    InteractionDynamics<KernelGradientCorrectionComplex> kernel_gradient_update(kernel_correction_matrix);
    SimpleDynamics<TimeStepInitialization> initialize_a_fluid_step(water_block);
    SimpleDynamics<NormalDirectionFromBodyShape> cylinder_normal_direction(cylinder);
    InteractionDynamics<fluid_dynamics::ViscousAccelerationWithWall> viscous_acceleration(water_block_complex);
    SimpleDynamics<NormalDirectionFromBodyShape> water_block_normal_direction(water_block);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(water_block, 0.5 / Dimensions);
    InteractionWithUpdate<FarFieldBoundary> variable_reset_in_boundary_condition(water_block_complex.getInnerRelation());
    InteractionWithUpdate<fluid_dynamics::FreeSurfaceIndicationComplex> surface_indicator(water_block_complex);
    InteractionDynamics<fluid_dynamics::SmearedSurfaceIndication> smeared_surface(water_block_complex.getInnerRelation());
    water_block.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
    //----------------------------------------------------------------------
    //	Compute the force exerted on solid body due to fluid pressure and viscosity
    //----------------------------------------------------------------------
    InteractionDynamics<solid_dynamics::ViscousForceFromFluid> viscous_force_on_solid(cylinder_contact);
    InteractionDynamics<solid_dynamics::AllForceAccelerationFromFluid> fluid_force_on_solid_update(cylinder_contact, viscous_force_on_solid);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToPlt write_real_body_states(io_environment, sph_system.real_bodies_);
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>>
        write_total_viscous_force_on_inserted_body(io_environment, viscous_force_on_solid, "TotalViscousForceOnSolid");
    ReducedQuantityRecording<ReduceDynamics<solid_dynamics::TotalForceFromFluid>>
        write_total_force_on_inserted_body(io_environment, fluid_force_on_solid_update, "TotalPressureForceOnSolid");
    ReducedQuantityRecording<ReduceDynamics<MaximumSpeed>> write_maximum_speed(io_environment, water_block);
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
    kernel_correction_matrix.exec();
    kernel_gradient_update.exec();
    kernel_correction_matrix_correction.exec();
    kernel_correction_matrix_cylinder.exec();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 200.0;
    Real output_interval = 1.0; /**< time stamps for output. */
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
    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            initialize_a_fluid_step.exec();
            Real dt = get_fluid_time_step_size.exec();
            viscous_acceleration.exec();
            pressure_relaxation.exec(dt);
            density_relaxation.exec(dt);

            integration_time += dt;
            GlobalStaticVariables::physical_time_ += dt;
            variable_reset_in_boundary_condition.exec();

            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                    << GlobalStaticVariables::physical_time_
                    << "	dt = " << dt << "\n";
            }
            number_of_iterations++;
        }

        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();

        write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
        write_total_force_on_inserted_body.writeToFile(number_of_iterations);

        write_maximum_speed.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    write_total_viscous_force_on_inserted_body.testResult();

    return 0;
}
