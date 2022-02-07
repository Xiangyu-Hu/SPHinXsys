/**
 * @file 	regression_test.cpp
 * @brief 	This is a test case based on diffusion, which can be used to 
            validate the generation of the converged database in a regression test. 
            It can be run successfully (using CMake's CTest) in Linux system installed with Python 3.
 * @author 	Bo Zhang and Xiangyu Hu
 */

#include "sphinxsys.h" //SPHinXsys Library
using namespace SPH;   //namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 0.2;
Real H = 0.2;
Real resolution_ref = H / 40.0;
Real BW = resolution_ref * 4;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Global parameters on material properties
//----------------------------------------------------------------------
Real diffusion_coff = 1.0e-3;
Real bias_coff = 0.0;
Real alpha = Pi / 4.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
Real initaltemperature = 0.0;
Real hightemperature = 1.0;
Real lowtemperature = 0.0;
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
std::vector<Vecd> createDiffusionDomain()
{
	//thermal solid domain geometry.
	std::vector<Vecd> diffusiondomain;
	diffusiondomain.push_back(Vecd(-BW, -BW));
	diffusiondomain.push_back(Vecd(-BW, H + BW));
	diffusiondomain.push_back(Vecd(L + BW, H + BW));
	diffusiondomain.push_back(Vecd(L + BW, -BW));
	diffusiondomain.push_back(Vecd(-BW, -BW));

	return diffusiondomain;
}

MultiPolygon createInnerDomain()
{
	//thermal solid inner domain geometry.
	std::vector<Vecd> innerdomain;
	innerdomain.push_back(Vecd(0.0, 0.0));
	innerdomain.push_back(Vecd(0.0, H));
	innerdomain.push_back(Vecd(L, H));
	innerdomain.push_back(Vecd(L, 0.0));
	innerdomain.push_back(Vecd(0.0, 0.0));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(innerdomain, ShapeBooleanOps::add);

	return multi_polygon;
}

MultiPolygon createLeftSideBoundary()
{
	//left isothermal boundary geometry.
	std::vector<Vecd> leftsideboundary;
	leftsideboundary.push_back(Vecd(-BW, -BW));
	leftsideboundary.push_back(Vecd(-BW, H + BW));
	leftsideboundary.push_back(Vecd(0.0, H));
	leftsideboundary.push_back(Vecd(0.0, 0.0));
	leftsideboundary.push_back(Vecd(-BW, -BW));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(leftsideboundary, ShapeBooleanOps::add);

	return multi_polygon;
}

MultiPolygon createOtherSideBoundary()
{
	//other side isothermal boundary geometry.
	std::vector<Vecd> othersideboundary;
	othersideboundary.push_back(Vecd(-BW, -BW));
	othersideboundary.push_back(Vecd(0.0, 0.0));
	othersideboundary.push_back(Vecd(L, 0.0));
	othersideboundary.push_back(Vecd(L, H));
	othersideboundary.push_back(Vecd(0.0, L));
	othersideboundary.push_back(Vecd(-BW, H + BW));
	othersideboundary.push_back(Vecd(L + BW, H + BW));
	othersideboundary.push_back(Vecd(L + BW, -BW));
	othersideboundary.push_back(Vecd(-BW, -BW));

	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(othersideboundary, ShapeBooleanOps::add);

	return multi_polygon;
}
//----------------------------------------------------------------------
//	Bodies with cases-dependent geometries (ComplexShape).
//----------------------------------------------------------------------
class DiffusionBody : public SolidBody
{
public:
	DiffusionBody(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(createDiffusionDomain(), ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	setup diffusion material properties.
//----------------------------------------------------------------------
class DiffusionMaterial
	: public DiffusionReaction<SolidParticles, Solid>
{
public:
	DiffusionMaterial()
		: DiffusionReaction<SolidParticles, Solid>({"Phi"})
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
		if (pos_n_[index_i][0] >= 0 && pos_n_[index_i][0] <= L && pos_n_[index_i][1] >= 0 && pos_n_[index_i][1] <= H)
		{
			species_n_[phi_][index_i] = initaltemperature;
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
//	Set left side boundary condition.
//----------------------------------------------------------------------
class LeftSideBoundaryCondition
	: public ConstrainDiffusionBodyRegion<SolidBody, SolidParticles, BodyPartByParticle, Solid>
{
protected:
	size_t phi_;
	void Update(size_t index_i, Real dt) override
	{
		species_n_[phi_][index_i] = hightemperature;
	};

public:
	LeftSideBoundaryCondition(SolidBody &diffusion_body, BodyPartByParticle &body_part)
		: ConstrainDiffusionBodyRegion<SolidBody, SolidParticles, BodyPartByParticle, Solid>(diffusion_body, body_part)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};
//----------------------------------------------------------------------
//	Set other side boundary condition.
//----------------------------------------------------------------------
class OtherSideBoundaryCondition
	: public ConstrainDiffusionBodyRegion<SolidBody, SolidParticles, BodyPartByParticle, Solid>
{
protected:
	size_t phi_;
	void Update(size_t index_i, Real dt) override
	{
		species_n_[phi_][index_i] = lowtemperature;
	};

public:
	OtherSideBoundaryCondition(SolidBody &diffusion_body, BodyPartByParticle &body_part)
		: ConstrainDiffusionBodyRegion<SolidBody, SolidParticles, BodyPartByParticle, Solid>(diffusion_body, body_part)
	{
		phi_ = material_->SpeciesIndexMap()["Phi"];
	};
};
//----------------------------------------------------------------------
//	Specify diffusion relaxation method.
//----------------------------------------------------------------------
class DiffusionBodyRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2<
		  SolidBody, SolidParticles, Solid,
		  RelaxationOfAllDiffussionSpeciesInner<SolidBody, SolidParticles, Solid>,
		  BodyRelationInner>
{
public:
	explicit DiffusionBodyRelaxation(BodyRelationInner &body_inner_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation){};
	virtual ~DiffusionBodyRelaxation(){};
};
//----------------------------------------------------------------------
//	an observer body to measure temperature at given positions.
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		/** A line of measuring points at the middle line. */
		size_t number_of_observation_points = 11;
		Real range_of_measure = L - BW;
		Real start_of_measure = BW;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec2d point_coordinate(0.5 * L, start_of_measure + range_of_measure * (Real)i / (Real)(number_of_observation_points - 1));
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
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	/** output environment. */
	In_Output in_output(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.cd
	//----------------------------------------------------------------------
	DiffusionBody diffusion_body(sph_system, "DiffusionBody");
	DiffusionReactionParticles<SolidParticles, Solid>
		diffusion_body_particles(diffusion_body, makeShared<DiffusionMaterial>());
	//----------------------------------------------------------------------
	//	Several parts of the diffusion body.
	//----------------------------------------------------------------------
	MultiPolygonShape inner_domain_shape(createInnerDomain());
	BodyRegionByParticle inner_domain(diffusion_body, "InnerDomian", inner_domain_shape);

	MultiPolygonShape left_boundary_shape(createLeftSideBoundary());
	BodyRegionByParticle left_boundary(diffusion_body, "LeftBoundaryDomain", left_boundary_shape);

	MultiPolygonShape other_boundary_shape(createOtherSideBoundary());
	BodyRegionByParticle other_boundary(diffusion_body, "OtherBoundaryDomain", other_boundary_shape);
	//----------------------------------------------------------------------
	//	Observer body
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
	/** Left and lower wall boundary conditions */
	LeftSideBoundaryCondition left_boundary_condition(diffusion_body, left_boundary);
	OtherSideBoundaryCondition other_boundary_condition(diffusion_body, other_boundary);
	/** Corrected configuration for diffusion body. */
	solid_dynamics::CorrectConfiguration correct_configuration(diffusion_body_inner_relation);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<SolidBody, SolidParticles, Solid> get_time_step_size(diffusion_body);
	/** Diffusion process for diffusion body. */
	DiffusionBodyRelaxation diffusion_relaxation(diffusion_body_inner_relation);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations of the simulation.
	//	Regression tests are also defined here.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(in_output, sph_system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<Real>>
		write_solid_temperature("Phi", in_output, temperature_observer_contact);
	RegressionTestDynamicTimeWarping<BodyReducedQuantityRecording<TotalAveragedParameterOnPartlyDiffusionBody<SolidBody, SolidParticles, Solid>>>
		write_solid_average_temperature_part(in_output, diffusion_body, inner_domain, "Phi");
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();
	correct_configuration.parallel_exec();
	setup_diffusion_initial_condition.parallel_exec();
	left_boundary_condition.parallel_exec();
	other_boundary_condition.parallel_exec();
	/** Output global basic parameters. */
	write_states.writeToFile(0);
	write_solid_temperature.writeToFile(0);
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 20.0;
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
				left_boundary_condition.parallel_exec();
				other_boundary_condition.parallel_exec();
				ite++;
				dt = get_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				if (ite % 100 == 0)
				{
					write_solid_temperature.writeToFile(ite);
					write_solid_average_temperature_part.writeToFile(ite);
					write_states.writeToFile(ite);
				}
			}
		}

		tick_count t2 = tick_count::now();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();
	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	//----------------------------------------------------------------------
	//	@ensemble_average_method.
	//	The first argument is the threshold of meanvalue convergence.
	//	The second argument is the threshold of variance convergence.
	//----------------------------------------------------------------------
	write_solid_temperature.generateDataBase(0.001, 0.001);
	//----------------------------------------------------------------------
	//	@dynamic_time_warping_method.
	//	The value is the threshold of dynamic time warping (dtw) distance.
	//----------------------------------------------------------------------
	write_solid_average_temperature_part.generateDataBase(0.001);

	return 0;
}