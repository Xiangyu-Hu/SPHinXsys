/*-----------------------------------------------------------------------------*
 *                       SPHinXsys: Lid Driven Cavity 2.5D                     *
 *-----------------------------------------------------------------------------*
 * This is the one of the basic test cases for efficient and accurate time     *
 * integration scheme investigation                                            *
 *-----------------------------------------------------------------------------*/
#include "sphinxsys.h" // SPHinXsys Library.
using namespace SPH;

// geometric data
Real particle_spacing = 0.01;

// material properties
Real rho = 1.0;          // reference density
Real u_lid = 1.0;        // lid velocity
Real SOS = 10.0 * u_lid; // numerical speed of sound

// non-Newtonian properties
Real K = 1.0;   // consistency index
Real n = 1.5;   // power index
Real tau_y = 0; // yield stress

Real min_shear_rate = 1e-3; // cutoff low shear rate
Real max_shear_rate = 1e+3; // cutoff high shear rate

// mesh geometry data
std::string path_to_lid_boundary = "./input/lid_boundary.stl";
std::string path_to_no_slip_boundary = "./input/no_slip_boundary.stl";
std::string path_to_slip_boundary = "./input/slip_boundary.stl";
std::string path_to_fluid = "./input/fluid.stl";

//----------------------------------------------------------------------
//	Complex shapes for wall boundary
//----------------------------------------------------------------------
class Lid_Boundary : public ComplexShape
{
    explicit Lid_Boundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_lid_boundary, Vecd(0.0, 0.0, 0.0), 1.0);
    }
};
class No_Slip_Boundary : public ComplexShape
{
    explicit No_Slip_Boundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_no_slip_boundary, Vecd(0.0, 0.0, 0.0), 1.0);
    }
};
class Slip_Boundary : public ComplexShape
{
    explicit Slip_Boundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_slip_boundary, Vecd(0.0, 0.0, 0.0), 1.0);
    }
};
class FluidFilling : public ComplexShape
{
    explicit FluidFilling(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(path_to_fluid, Vecd(0.0, 0.0, 0.0), 1.0);
    }
};

int main(int ac, char *av[])
{
    //	Build up an SPHSystem
    BoundingBox system_domain_bounds(Vecd(-0.3, -0.7, -0.7), Vecd(0.3, 0.7, 0.7));
    SPHSystem sph_system(system_domain_bounds, particle_spacing);
    sph_system.handleCommandlineOptions(ac, av)->setIOEnvironment();

    //	Creating bodies with corresponding materials and particles
    FluidBody fluid(sph_system, makeShared<FluidFilling>("FluidBody"));
    fluid.defineParticlesAndMaterial<BaseParticles, HerschelBulkleyFluid>(rho, SOS, K, n, tau_y, min_shear_rate, max_shear_rate);
    fluid.generateParticles<ParticleGeneratorLattice>();

    SolidBody no_slip_boundary(sph_system, makeShared<No_Slip_Boundary>("NoSlipWall"));
    no_slip_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    no_slip_boundary.generateParticles<ParticleGeneratorLattice>();
    no_slip_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody slip_boundary(sph_system, makeShared<Slip_Boundary>("SlipWall"));
    slip_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    slip_boundary.generateParticles<ParticleGeneratorLattice>();
    slip_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

    SolidBody lid_boundary(sph_system, makeShared<Lid_Boundary>("SlipWall"));
    lid_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
    lid_boundary.generateParticles<ParticleGeneratorLattice>();
    lid_boundary.addBodyStateForRecording<Vec3d>("NormalDirection");

    //	Define body relation map
    InnerRelation fluid_inner(fluid);
    ContactRelation fluid_all_walls(fluid, {&lid_boundary, &slip_boundary, &no_slip_boundary});
    ContactRelation fluid_no_slip(fluid, {&no_slip_boundary});
    ContactRelation fluid_slip(fluid, {&slip_boundary});
    ContactRelation fluid_lid(fluid, {&lid_boundary});

    ComplexRelation fluid_walls_complex(fluid_inner, fluid_all_walls)

        //	Define the numerical methods used in the simulation
        Dynamics1Level<fluid_dynamics::Integration1stHalfWithWallRiemann>
            pressure_relaxation(fluid_inner, fluid_all_walls);
    Dynamics1Level<fluid_dynamics::Integration2ndHalfWithWallRiemann> density_relaxation(fluid_inner, fluid_all_walls);

    InteractionDynamics<fluid_dynamics::VelocityGradientInner> vel_grad_calc_inner(fluid_inner);
    InteractionDynamics<fluid_dynamics::VelocityGradientContact> vel_grad_calc_contact(fluid_no_slip);
    InteractionDynamics<fluid_dynamics::ShearRateDependentViscosity> shear_rate_calculation(fluid_no_slip);
    InteractionWithUpdate<fluid_dynamics::ViscousShearRateDependent> viscous_acceleration(fluid_inner, fluid_no_slip);

    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> get_fluid_advection_time_step_size(fluid_inner, u_lid);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> get_fluid_time_step_size(fluid_inner);
    //! SimpleDynamics<NormalDirectionFromBodyShape> wall_boundary_normal_direction(wall_boundary); [needed?]

    //	Define the methods for I/O operations, observations
    BodyStatesRecordingToVtp write_fluid_states(sph_system.real_bodies_);

    //	Prepare the simulation
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    //! wall_boundary_normal_direction.exec(); [needed?]

    //	Setup for time-stepping control
    size_t number_of_iterations = sph_system.RestartStep();
    Real end_time = 1.0;
    Real output_interval = end_time / 100.0;
    Real dt = 0.0;

    //	First output before the main loop.
    write_fluid_states.writeToFile(0);

    while (GlobalStaticVariables::physical_time_ < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real Dt = get_fluid_advection_time_step_size.exec();
            update_density_by_summation.exec(Dt);

            vel_grad_calc_inner.exec(Dt);
            vel_grad_calc_contact.exec(Dt);
            shear_rate_calculation.exec(Dt);
            viscous_acceleration.exec(Dt);

            while (relaxation_time < Dt)
            {
                // TODO Write Viscous Time Step Criterion
                //  Real Dt = SMIN(get_fluid_advection_time_step_size.exec(), get_fluid_viscous_time_step_size.exec());
                pressure_relaxation.exec(dt);
                density_relaxation.exec(dt);
                dt = get_fluid_time_step_size.exec();
                relaxation_time += dt;
                integration_time += dt;
                GlobalStaticVariables::physical_time_ += dt;
            }
            fluid.updateCellLinkedListWithParticleSort(100);
            fluid_walls_complex.updateConfiguration();
        }
        write_water_block_states.writeToFile();
    }

    return;
}