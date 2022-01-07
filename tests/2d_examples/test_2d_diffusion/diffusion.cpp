/**
 * @file 	diffusion.cpp
 * @brief 	This is the first test to validate our anisotropic diffusion solver.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library
using namespace SPH;   //Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 2.0;
Real H = 0.4;
Real resolution_ref = H / 40.0;
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coff = 1.0e-4;
Real bias_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
StdVec<std::string> species_name_list{"Phi"};
//----------------------------------------------------------------------
//	Geometric shapes used in the case.
//----------------------------------------------------------------------
std::vector<Vecd> createShape()
{
	//geometry
	std::vector<Vecd> shape;
	shape.push_back(Vecd(0.0, 0.0));
	shape.push_back(Vecd(0.0, H));
	shape.push_back(Vecd(L, H));
	shape.push_back(Vecd(L, 0.0));
	shape.push_back(Vecd(0.0, 0.0));

	return shape;
}
//----------------------------------------------------------------------
//	Define SPH bodies.
//----------------------------------------------------------------------
class DiffusionBody : public SolidBody
{
public:
	DiffusionBody(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createShape(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	Setup diffusion material properties.
//----------------------------------------------------------------------
class DiffusionMaterial
	: public DiffusionReaction<SolidParticles, Solid>
{
public:
	DiffusionMaterial()
		: DiffusionReaction<SolidParticles, Solid>(species_name_list)
	{
		initializeAnDiffusion<DirectionalDiffusion>("Phi", "Phi", diffusion_coff, bias_coff, bias_direction);
	};
};
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
class DiffusionInitialCondition
	: public DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{

		if (pos_n_[index_i][0] >= 0.45 && pos_n_[index_i][0] <= 0.55)
		{
			species_n_[phi_][index_i] = 1.0;
		}
		if (pos_n_[index_i][0] >= 1.0)
		{
			species_n_[phi_][index_i] = exp(-2500.0 * ((pos_n_[index_i][0] - 1.5) * (pos_n_[index_i][0] - 1.5)));
		}
	};

public:
	explicit DiffusionInitialCondition(SolidBody &diffusion_body)
		: DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>(diffusion_body)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};
//----------------------------------------------------------------------
//	Specify diffusion relaxation method.
//----------------------------------------------------------------------
class DiffusionBodyRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2<SolidBody, SolidParticles, Solid,
												RelaxationOfAllDiffussionSpeciesInner<SolidBody, SolidParticles, Solid>,
												BodyRelationInner>
{
public:
	explicit DiffusionBodyRelaxation(BodyRelationInner &body_inner_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation){};
	virtual ~DiffusionBodyRelaxation(){};
};
//----------------------------------------------------------------------
//	An observer particle generator.
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		size_t number_of_observation_points = 11;
		Real range_of_measure = 0.9 * L;
		Real start_of_measure = 0.05 * L;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec2d point_coordinate(range_of_measure * (Real)i / (Real)(number_of_observation_points - 1) + start_of_measure, 0.5 * H);
			positions_volumes_.push_back(std::make_pair(point_coordinate, 0.0));
		}
	}
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** output environment. */
	In_Output in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	DiffusionBody diffusion_body(sph_system, "DiffusionBody");
	DiffusionReactionParticles<SolidParticles, Solid>
		diffusion_body_particles(diffusion_body, makeShared<DiffusionMaterial>());
	//----------------------------------------------------------------------
	//	Particle and body creation of fluid observers.
	//----------------------------------------------------------------------
	ObserverBody temperature_observer(sph_system, "TemperatureObserver");
	ObserverParticles temperature_observer_particles(temperature_observer, makeShared<ObserverParticleGenerator>());
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner diffusion_body_inner_relation(diffusion_body);
	BodyRelationContact temperature_observer_contact(temperature_observer, {&diffusion_body});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	DiffusionInitialCondition setup_diffusion_initial_condition(diffusion_body);
	/** Corrected configuration for diffusion body. */
	solid_dynamics::CorrectConfiguration correct_configuration(diffusion_body_inner_relation);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<SolidBody, SolidParticles, Solid> get_time_step_size(diffusion_body);
	/** Diffusion process for diffusion body. */
	DiffusionBodyRelaxation diffusion_relaxation(diffusion_body_inner_relation);
	/** Periodic BCs. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList periodic_condition_y(diffusion_body, yAxis);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(in_output, sph_system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<indexScalar, Real>>
		write_solid_temperature("Phi", in_output, temperature_observer_contact);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	periodic_condition_y.update_cell_linked_list_.parallel_exec();
	sph_system.initializeSystemConfigurations();
	correct_configuration.parallel_exec();
	setup_diffusion_initial_condition.exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 1.0;
	Real End_Time = T0;
	Real Output_Time = 0.1 * End_Time;
	Real Observe_time = 0.1 * Output_Time;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_states.writeToFile();
	write_solid_temperature.writeToFile();
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < Output_Time)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Observe_time)
			{
				if (ite % 1 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: "
							  << dt << "\n";
				}

				diffusion_relaxation.parallel_exec(dt);

				ite++;
				dt = get_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}

		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		write_solid_temperature.writeToFile(ite);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	write_solid_temperature.newResultTest();

	return 0;
}
