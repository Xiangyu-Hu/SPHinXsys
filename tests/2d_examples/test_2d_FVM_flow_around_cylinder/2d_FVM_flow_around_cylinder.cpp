/**
 * @file 	2d_FVM_flow_around_cylinder.cpp
 * @brief 	This is the test file for the weakly compressible viscous flow around a cylinder in FVM.
 * @details We consider a flow passing by a cylinder in 2D in FVM framework.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
#include "FVM_flow_around_cylinder_2d_def.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 50.0;                  /**< Channel length. */
Real DH = 30.0;                  /**< Channel height. */
Real resolution_ref = 1.0 / 5.0; /**< Initial reference particle spacing. */
Real DL_sponge = 2.0;            /**< Sponge region to impose inflow condition. */
Real DH_sponge = 2.0;            /**< Sponge region to impose inflow condition. */
Real cylinder_radius = 1.0;      /**< Radius of the cylinder. */
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1.0;                                       /**< Density. */
Real U_f = 1.0;                                          /**< freestream velocity. */
Real c_f = 10.0 * U_f;                                   /**< Speed of sound. */
Real Re = 100.0;                                         /**< Reynolds number. */
Real mu_f = rho0_f * U_f * (2.0 * cylinder_radius) / Re; /**< Dynamics viscosity. */
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string ansys_mesh_file_path = "./input/fluent_0.3.msh";
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
        MultiPolygon water_block(createWaterBlockShape());
        add<MultiPolygonShape>(water_block, "WaterBlock");
    }
};
//----------------------------------------------------------------------
//	Case dependent boundary condition
//----------------------------------------------------------------------
class FACBoundaryConditionSetup : public BoundaryConditionSetupInFVM
{
  public:
    FACBoundaryConditionSetup(BaseInnerRelationInFVM &inner_relation, GhostCreationFromMesh &ghost_creation)
        : BoundaryConditionSetupInFVM(inner_relation, ghost_creation),
          fluid_(DynamicCast<WeaklyCompressibleFluid>(this, particles_->getBaseMaterial())) {};
    virtual ~FACBoundaryConditionSetup() {};

    void applyNonSlipWallBoundary(size_t ghost_index, size_t index_i) override
    {
        vel_[ghost_index] = -vel_[index_i];
        p_[ghost_index] = p_[index_i];
        rho_[ghost_index] = rho_[index_i];
    }
    void applyFarFieldBoundary(size_t ghost_index) override
    {
        Vecd far_field_velocity(1.0, 0.0);
        Real far_field_density = 1.0;
        Real far_field_pressure = fluid_.getPressure(far_field_density);

        vel_[ghost_index] = far_field_velocity;
        p_[ghost_index] = far_field_pressure;
        rho_[ghost_index] = far_field_density;
    }

  protected:
    Fluid &fluid_;
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    // read data from ANSYS mesh.file
    ANSYSMesh ansys_mesh(ansys_mesh_file_path);
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    BoundingBoxd system_domain_bounds(Vec2d(-DL_sponge, -DH_sponge), Vec2d(DL, DH + DH_sponge));
    SPHSystem sph_system(system_domain_bounds, ansys_mesh.MinMeshEdge());
    sph_system.handleCommandlineOptions(ac, av);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBlock"));
    water_block.defineClosure<WeaklyCompressibleFluid, Viscosity>(ConstructArgs(rho0_f, c_f), mu_f);
    Ghost<ReserveSizeFactor> ghost_boundary(0.5);
    water_block.generateParticlesWithReserve<BaseParticles, UnstructuredMesh>(ghost_boundary, ansys_mesh);
    GhostCreationFromMesh ghost_creation(water_block, ansys_mesh, ghost_boundary);
    //----------------------------------------------------------------------
    //	Define body relation map.
    //----------------------------------------------------------------------
    InnerRelationInFVM water_block_inner(water_block, ansys_mesh);
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    /** Here we introduce the limiter in the Riemann solver and 0 means the no extra numerical dissipation.
    the value is larger, the numerical dissipation larger*/
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration1stHalfInnerRiemann> pressure_relaxation(water_block_inner, 200.0);
    InteractionWithUpdate<fluid_dynamics::EulerianIntegration2ndHalfInnerRiemann> density_relaxation(water_block_inner, 200.0);
    FACBoundaryConditionSetup boundary_condition_setup(water_block_inner, ghost_creation);
    ReduceDynamics<fluid_dynamics::WCAcousticTimeStepSizeInFVM> get_fluid_time_step_size(water_block, ansys_mesh.MinMeshEdge());
    InteractionWithUpdate<fluid_dynamics::ViscousForceInner> viscous_force(water_block_inner);
    //----------------------------------------------------------------------
    //	Compute the force exerted on solid body due to fluid pressure and viscosity
    //----------------------------------------------------------------------
    InteractionDynamics<fluid_dynamics::ViscousForceFromFluidInFVM> viscous_force_on_solid(water_block_inner, ghost_creation.each_boundary_type_contact_real_index_);
    InteractionDynamics<fluid_dynamics::PressureForceFromFluidInFVM<decltype(density_relaxation)>> pressure_force_on_solid(water_block_inner, ghost_creation.each_boundary_type_contact_real_index_);
    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToMeshVtu write_real_body_states(water_block, ansys_mesh);
    write_real_body_states.addToWrite<Real>(water_block, "Density");
    RegressionTestDynamicTimeWarping<ReducedQuantityRecording<QuantitySummation<Vecd>>> write_total_viscous_force_on_inserted_body(water_block, "ViscousForceOnSolid");
    ReducedQuantityRecording<QuantitySummation<Vecd>> write_total_pressure_force_on_inserted_body(water_block, "PressureForceOnSolid");
    ReducedQuantityRecording<MaximumSpeed> write_maximum_speed(water_block);
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    water_block_inner.updateConfiguration();
    //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    int screen_output_interval = 1000;
    Real end_time = 70.0;
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
            boundary_condition_setup.resetBoundaryConditions();
            viscous_force.exec();
            pressure_relaxation.exec(dt);
            boundary_condition_setup.resetBoundaryConditions();
            density_relaxation.exec(dt);

            integration_time += dt;
            physical_time += dt;
            if (number_of_iterations % screen_output_interval == 0)
            {
                std::cout << std::fixed << std::setprecision(9) << "N=" << number_of_iterations << "	Time = "
                          << physical_time
                          << "	dt = " << dt << "\n";
                viscous_force_on_solid.exec();
                pressure_force_on_solid.exec();
                write_total_viscous_force_on_inserted_body.writeToFile(number_of_iterations);
                write_total_pressure_force_on_inserted_body.writeToFile(number_of_iterations);
                write_maximum_speed.writeToFile(number_of_iterations);
            }
            number_of_iterations++;
        }
        TickCount t2 = TickCount::now();
        write_real_body_states.writeToFile();
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
