#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.

template<class ReferenceFunc, class SYCLFunc>
void benchmark(ReferenceFunc&& refFunc, SYCLFunc&& syclFunc, const std::string& func_name, std::size_t iterations) {
    TickCount timer_start = TickCount::now();
    for(auto i = 0; i < iterations; ++i)
        syclFunc();
    TimeInterval sycl_interval = TickCount::now() - timer_start;
    std::cout << "SYCL " <<  func_name << " computation: " << sycl_interval.seconds() << " seconds" << std::endl;

    timer_start = TickCount::now();
    for(auto i = 0; i < iterations; ++i)
        refFunc();
    TimeInterval ref_interval = TickCount::now() - timer_start;
    std::cout << "Reference " <<  func_name << " computation: " << ref_interval.seconds() << " seconds" << std::endl;

    std::cout << "Speedup " << func_name << " (main loop): " << ref_interval.seconds()/sycl_interval.seconds() << 'x' << std::endl;
}

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;					/**< Water tank length. */
Real DH = 5.366;					/**< Water tank height. */
Real LL = 2.0;						/**< Water column length. */
Real LH = 1.0;						/**< Water column height. */
Real particle_spacing_ref = 0.005;	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4; /**< Thickness of tank wall. */
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_f = 1.0;						 /**< Reference density of fluid. */
Real gravity_g = 1.0;					 /**< Gravity. */
Real U_max = 2.0 * sqrt(gravity_g * LH); /**< Characteristic velocity. */
Real c_f = 10.0 * U_max;				 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin
Vec2d water_block_translation = water_block_halfsize;	// translation to global coordinates
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex shape for wall boundary, note that no partial overlap is allowed
//	for the shapes in a complex shape.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<TransformShape<GeometricShapeBox>>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
	SPHSystem sph_system(system_domain_bounds, particle_spacing_ref);
	sph_system.handleCommandlineOptions(ac, av);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
    FluidBody water_block(
            sph_system, makeShared<TransformShape<GeometricShapeBox>>(
                    Transform(water_block_translation), water_block_halfsize, "WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    FluidBody water_block_sycl(
            sph_system, makeShared<TransformShape<GeometricShapeBox>>(
                    Transform(water_block_translation), water_block_halfsize, "WaterBody"));
    water_block_sycl.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block_sycl.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary(sph_system, makeShared<WallBoundary>("WallBoundary"));
	wall_boundary.defineParticlesAndMaterial<SolidParticles, Solid>();
	wall_boundary.generateParticles<ParticleGeneratorLattice>();
	wall_boundary.addBodyStateForRecording<Vecd>("NormalDirection");

    SolidBody wall_boundary_sycl(sph_system, makeShared<WallBoundary>("WallBoundary"));
    wall_boundary_sycl.defineParticlesAndMaterial<SolidParticles, Solid>();
    wall_boundary_sycl.generateParticles<ParticleGeneratorLattice>();
    wall_boundary_sycl.addBodyStateForRecording<Vecd>("NormalDirection");
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation water_block_complex(water_block, {&wall_boundary});
	ComplexRelation water_block_complex_sycl(water_block_sycl, {&wall_boundary_sycl});

    TickCount timer_mem = TickCount::now();
    water_block_sycl.getBaseParticles().registerExtraDeviceMemory();
    TimeInterval tt_mem = TickCount::now() - timer_mem;

    water_block.getBaseParticles().registerExtraDeviceMemory();

    SharedPtr<Gravity> gravity_ptr = makeSharedDevice<Gravity>(Vecd(0.0, -gravity_g));

    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall, ParallelSYCLDevicePolicy> fluid_pressure_relaxation_sycl(water_block_complex_sycl);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex, ParallelSYCLDevicePolicy> fluid_density_by_summation_sycl(water_block_complex_sycl);
	SimpleDynamics<TimeStepInitialization, ParallelSYCLDevicePolicy> fluid_step_initialization_sycl(water_block_sycl, gravity_ptr);
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize, ParallelSYCLDevicePolicy> fluid_advection_time_step_sycl(water_block_sycl, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize, ParallelSYCLDevicePolicy> fluid_acoustic_time_step_sycl(water_block_sycl);

    Dynamics1Level<fluid_dynamics::Integration1stHalfRiemannWithWall> fluid_pressure_relaxation(water_block_complex);
    InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceComplex> fluid_density_by_summation(water_block_complex);
    SimpleDynamics<TimeStepInitialization> fluid_step_initialization(water_block, gravity_ptr);
    ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(water_block, U_max);
    ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(water_block);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();

    //----------------------------------------------------------------------
    //	Copy data to device
    //----------------------------------------------------------------------
    timer_mem = TickCount::now();
    water_block_sycl.getBaseParticles().copyToDeviceMemory();
    water_block_sycl.getBaseParticles().copyToExtraDeviceMemory();
    water_block_complex_sycl.getInnerRelation().copyInnerConfigurationToDevice();
    water_block_complex_sycl.getContactRelation().copyContactConfigurationToDevice();
    tt_mem += TickCount::now() - timer_mem;

    executionQueue.setWorkGroupSize(32);

    //----------------------------------------------------------------------
    //	Main loop benchmarks
    //----------------------------------------------------------------------
    std::cout << "Number of particles: " << water_block_sycl.getBaseParticles().total_real_particles_ << std::endl;

    std::size_t iterations = 2000;
    std::cout << "Number of iterations per test: " << iterations << std::endl;

    std::cout << "------------" << std::endl;
    std::cout << "SYCL memory operations: " << tt_mem.seconds()
              << " seconds." << std::endl;

    std::cout << "------------" << std::endl;
    benchmark([&](){
                  fluid_step_initialization.exec();
                  Real advection_dt = fluid_advection_time_step.exec();
                  fluid_density_by_summation.exec();
                  Real acoustic_dt = fluid_acoustic_time_step.exec();
                  fluid_pressure_relaxation.exec();
              },
              [&](){
                  fluid_step_initialization_sycl.exec();
                  Real advection_dt = fluid_advection_time_step_sycl.exec();
                  fluid_density_by_summation_sycl.exec();
                  Real acoustic_dt = fluid_acoustic_time_step_sycl.exec();
                  fluid_pressure_relaxation_sycl.exec();
              },
              "all methods", iterations);

    std::cout << "------------" << std::endl;
    benchmark([&](){ fluid_step_initialization.exec(); },
              [&](){ fluid_step_initialization_sycl.exec(); },
              "fluid_step_initialization", iterations);

    std::cout << "------------" << std::endl;
    benchmark([&](){ Real advection_dt = fluid_advection_time_step.exec(); },
              [&](){ Real advection_dt = fluid_advection_time_step_sycl.exec(); },
              "fluid_advection_time_step", iterations);

    std::cout << "------------" << std::endl;
    benchmark([&](){ fluid_density_by_summation.exec(); },
              [&](){ fluid_density_by_summation_sycl.exec(); },
              "fluid_density_by_summation", iterations);

    std::cout << "------------" << std::endl;
    benchmark([&](){ Real acoustic_dt = fluid_acoustic_time_step.exec(); },
              [&](){ Real acoustic_dt = fluid_acoustic_time_step_sycl.exec(); },
              "fluid_acoustic_time_step", iterations);

    std::cout << "------------" << std::endl;
    benchmark([&](){ fluid_pressure_relaxation.exec(); },
              [&](){ fluid_pressure_relaxation_sycl.exec(); },
              "fluid_pressure_relaxation", iterations);

	return 0;
};
