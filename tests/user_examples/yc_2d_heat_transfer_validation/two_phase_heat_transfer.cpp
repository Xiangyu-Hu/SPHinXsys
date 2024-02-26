/**
 * @file 	two_phase_heat_transfer.cpp
 * @brief 	Heat Transfer in Slabs
 */
#include "sphinxsys.h" //SPHinXsys Library.
#define PI (3.14159265358979323846)
#include "two_phase_heat_transfer_particles.h"
#include "extra_two_phase_heat_transfer.h"
#include "extra_two_phase_heat_transfer.hpp"
using namespace SPH;   // Namespace cite here.

//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_r = 1000.0;					 /**< Reference density of right material. */
Real rho0_l = 1.226;						 /**< Reference density of left material. */
Real c_p_r = 4.179;
Real c_p_l = 1.012;
Real k_r = 0.620;
Real k_l = 0.0254;
Real diffusion_coff_r = k_r / (c_p_r * rho0_r);
Real diffusion_coff_l = k_l / (c_p_l * rho0_l);
Real dp = 0.0125;	/**< Initial reference particle spacing. */

Real c_ = 1.0;				 /**< Reference sound speed. */
//----------------------------------------------------------------------
//	Geometric elements used in shape modeling.
//----------------------------------------------------------------------
std::vector<Vecd> createRightBlockShape()
{
	std::vector<Vecd> right_block_shape;
	right_block_shape.push_back(Vecd(40*dp, 0.0));
	right_block_shape.push_back(Vecd(40*dp, 40*dp));
	right_block_shape.push_back(Vecd(80*dp, 40*dp));
	right_block_shape.push_back(Vecd(80*dp, 0.0));
	right_block_shape.push_back(Vecd(40*dp, 0.0));
	return right_block_shape;
}

std::vector<Vecd> createLeftBlockShape()
{
	std::vector<Vecd> left_block_shape;
	left_block_shape.push_back(Vecd(0.0, 0.0));
	left_block_shape.push_back(Vecd(0.0, 40*dp));
	left_block_shape.push_back(Vecd(40*dp, 40*dp));
	left_block_shape.push_back(Vecd(40*dp, 0.0));
	left_block_shape.push_back(Vecd(0.0, 0.0));
	return left_block_shape;
}
//----------------------------------------------------------------------
//	cases-dependent geometric shape for right block.
//----------------------------------------------------------------------
class RightBlock : public MultiPolygonShape
{
public:
	explicit RightBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createRightBlockShape(), ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	cases-dependent geometric shape for left block.
//----------------------------------------------------------------------
class LeftBlock : public MultiPolygonShape
{
public:
	explicit LeftBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createLeftBlockShape(), ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion right body
//----------------------------------------------------------------------
class ThermoRightBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
  public:
    ThermoRightBodyMaterial()
        : DiffusionReaction<WeaklyCompressibleFluid>({"Phi"}, SharedPtr<NoReaction>(), rho0_r, c_)
    {
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff_r);
    };
};
using DiffusionRightParticles = DiffusionReactionParticles<BaseParticles, ThermoRightBodyMaterial>;
//class ThermoRightBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
//{
//public:
//	ThermoRightBodyMaterial()
//		: DiffusionReaction<WeaklyCompressibleFluid>({ "Phi" }, rho0_r,c_)
//	{
//		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff_r);
//	};
//};
//----------------------------------------------------------------------
//	Setup heat conduction material properties for diffusion left body
//----------------------------------------------------------------------
class ThermoLeftBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
{
  public:
    ThermoLeftBodyMaterial() : DiffusionReaction<WeaklyCompressibleFluid>({"Phi"}, SharedPtr<NoReaction>(),rho0_l,c_)
    {
        // only default property is given, as no heat transfer within solid considered here.
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi");
    };
};
using DiffusionLiftParticles = DiffusionReactionParticles<BaseParticles, ThermoLeftBodyMaterial>;

//class ThermoLeftBodyMaterial : public DiffusionReaction<WeaklyCompressibleFluid>
//{
//public:
//	ThermoLeftBodyMaterial() : DiffusionReaction<WeaklyCompressibleFluid>({ "Phi" }, rho0_l,c_)
//	{
//		// only default property is given, as no heat transfer within solid considered here.
//		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi",diffusion_coff_l);
//	};
//};
//----------------------------------------------------------------------
//	Application dependent right body initial condition
//----------------------------------------------------------------------
class ThermoRightBodyInitialCondition
	: public TwoPhaseDiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

public:
	explicit ThermoRightBodyInitialCondition(SPHBody& sph_body)
		: TwoPhaseDiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		species_n_[phi_][index_i] = 1.0;
		thermal_conductivity_[index_i] = k_r;
	};
};

class ThermoRightBodyInitialCondition
    : public DiffusionReactionInitialCondition<DiffusionRightParticles>
{
  protected:
    size_t phi_;

  public:
    explicit ThermoRightBodyInitialCondition(SPHBody &sph_body)
        : DiffusionReactionInitialCondition<DiffusionRightParticles>(sph_body)
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    };

    void update(size_t index_i, Real dt)
	{
		all_species_[phi_][index_i] = 1.0;
		thermal_conductivity_[index_i] = k_r;
    };
};
//----------------------------------------------------------------------
//	Application dependent left body initial condition
//----------------------------------------------------------------------
class ThermoLeftBodyInitialCondition
	: public TwoPhaseDiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>
{
protected:
	size_t phi_;

public:
	explicit ThermoLeftBodyInitialCondition(SPHBody &sph_body)
		: TwoPhaseDiffusionReactionInitialCondition<FluidParticles, WeaklyCompressibleFluid>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.SpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		species_n_[phi_][index_i] = 0.0;
		thermal_conductivity_[index_i] = k_l;
	};
};

//----------------------------------------------------------------------
//	Set thermal relaxation between different bodies
//----------------------------------------------------------------------
class ThermalRelaxationComplex
	: public TwoPhaseRelaxationOfAllDiffusionSpeciesRK2<
	TwoPhaseRelaxationOfAllDiffusionSpeciesComplex<
	FluidParticles, WeaklyCompressibleFluid, FluidParticles, WeaklyCompressibleFluid>>
{
public:
	explicit ThermalRelaxationComplex(ComplexRelation &body_complex_relation)
		: TwoPhaseRelaxationOfAllDiffusionSpeciesRK2(body_complex_relation) {};
	virtual ~ThermalRelaxationComplex() {};
};
//----------------------------------------------------------------------
//	An observer particle generator.
//----------------------------------------------------------------------
class TemperatureObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit TemperatureObserverParticleGenerator(SPHBody &sph_body)
		: ObserverParticleGenerator(sph_body)
	{
		size_t number_of_observation_points = 80;
		
		for(int i = 0;i< number_of_observation_points;i++)
		{
			positions_.push_back(Vecd(i*dp,20*dp));
		}
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
	BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(80*dp, 40*dp));
	SPHSystem sph_system(system_domain_bounds, dp);
	sph_system.handleCommandlineOptions(ac, av);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	FluidBody right_block(sph_system, makeShared<RightBlock>("RightBody"));
	right_block.defineAdaptation<SPHAdaptation>(1.0, 1.0);
	right_block.defineParticlesAndMaterial<TwoPhaseDiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermoRightBodyMaterial>();
	right_block.generateParticles<ParticleGeneratorLattice>();

	FluidBody left_block(sph_system, makeShared<LeftBlock>("LeftBody"));
	left_block.defineAdaptation<SPHAdaptation>(1.0, 1.0);
	left_block.defineParticlesAndMaterial<TwoPhaseDiffusionReactionParticles<FluidParticles, WeaklyCompressibleFluid>, ThermoLeftBodyMaterial>();
	left_block.generateParticles<ParticleGeneratorLattice>();

	ObserverBody temperature_observer(sph_system, "TemperatureObserver");
	temperature_observer.generateParticles<TemperatureObserverParticleGenerator>();

	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation right_left_complex(right_block, { &left_block });
	ComplexRelation left_right_complex(left_block, { &right_block });
	ContactRelation temperature_observer_contact(temperature_observer, RealBodyVector{ &right_block,&left_block });
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	/** Initialize particle acceleration. */
	SimpleDynamics<ThermoRightBodyInitialCondition> thermo_right_initial_condition(right_block);
	SimpleDynamics<ThermoLeftBodyInitialCondition> thermo_left_initial_condition(left_block);
	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_right(right_block);
	GetDiffusionTimeStepSize<FluidParticles, WeaklyCompressibleFluid> get_thermal_time_step_left(left_block);
	ThermalRelaxationComplex thermal_relaxation_complex_right(right_left_complex);
	ThermalRelaxationComplex thermal_relaxation_complex_left(left_right_complex);

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp body_states_recording(io_environment, sph_system.real_bodies_);
	RestartIO restart_io(io_environment, sph_system.real_bodies_);
	ObservedQuantityRecording<Real> write_temperature("Phi", io_environment, temperature_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	//----------------------------------------------------------------------
	//	Load restart file if necessary.
	//----------------------------------------------------------------------
	if (sph_system.RestartStep() != 0)
	{
		GlobalStaticVariables::physical_time_ = restart_io.readRestartFiles(sph_system.RestartStep());
		right_block.updateCellLinkedList();
		left_block.updateCellLinkedList();
		right_left_complex.updateConfiguration();
		left_right_complex.updateConfiguration();
	}
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	size_t number_of_iterations = sph_system.RestartStep();
	int screen_output_interval = 100;
	int observation_sample_interval = screen_output_interval * 2;
	int restart_output_interval = screen_output_interval * 10;
	int ite=0.0;
	Real end_time = 1.0;
	Real Output_Time = 0.1 * end_time;
	Real Observe_time =  0.05*Output_Time;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	tick_count::interval_t interval_computing_time_step;
	tick_count::interval_t interval_computing_fluid_pressure_relaxation;
	tick_count::interval_t interval_updating_configuration;
	tick_count time_instance;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	thermo_right_initial_condition.parallel_exec();
	thermo_left_initial_condition.parallel_exec();
	body_states_recording.writeToFile(0);
	write_temperature.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------

	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		/** Integrate time (loop) until the next output time. */
		while (integration_time < Output_Time)
		{
			time_instance = tick_count::now();
			interval_computing_time_step += tick_count::now() - time_instance;
			time_instance = tick_count::now();
			Real relaxation_time = 0.0;
			while (relaxation_time < Observe_time)
			{
				Real dt_thermal_right = get_thermal_time_step_right.parallel_exec();
				Real dt_thermal_left = get_thermal_time_step_left.parallel_exec();
				dt = SMIN(dt_thermal_right, dt_thermal_left);

				thermal_relaxation_complex_left.parallel_exec(dt);
				thermal_relaxation_complex_right.parallel_exec(dt);

				if (ite % 100== 0)
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}
				ite++;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}

			interval_computing_fluid_pressure_relaxation += tick_count::now() - time_instance;
			/** Update cell linked list and configuration. */
			right_block.updateCellLinkedListWithParticleSort(100);
			left_block.updateCellLinkedListWithParticleSort(100);

			right_left_complex.updateConfiguration();
			left_right_complex.updateConfiguration();
			temperature_observer_contact.updateConfiguration();
			write_temperature.writeToFile();
			time_instance = tick_count::now();
			interval_updating_configuration += tick_count::now() - time_instance;
		}

		tick_count t2 = tick_count::now();
		body_states_recording.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds()
		<< " seconds." << std::endl;
	return 0;
};
