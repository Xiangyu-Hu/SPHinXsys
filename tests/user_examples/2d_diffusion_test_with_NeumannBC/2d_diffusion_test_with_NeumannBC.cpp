/**
 * @file 	2d_diffusion_test_with_NeumannBC.cpp
 * @brief 	2D diffusion test of diffusion problem with Neumann boundary condition.
 * @details This is a case to implement Neumann boundary condition.
 * @author 	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
#include "2d_diffusion_test_with_NeumannBC.h"

using namespace SPH; //Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 100.0;
Real BW = resolution_ref * 2.0;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coff = 1;
std::array<std::string, 1> species_name_list{ "Phi" };
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 0.0;
Real left_temperature = 300.0;
Real right_temperature = 350.0;
Real heat_flux = 2000.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
std::vector<Vecd> createThermalDomain()
{
	std::vector<Vecd> thermalDomainShape;
	thermalDomainShape.push_back(Vecd(0.0, 0.0));
	thermalDomainShape.push_back(Vecd(0.0, H));
	thermalDomainShape.push_back(Vecd(L, H));
	thermalDomainShape.push_back(Vecd(L, 0.0));
	thermalDomainShape.push_back(Vecd(0.0, 0.0));

	return thermalDomainShape;
}

std::vector<Vecd> left_temperature_region
{
	Vecd(0.3 * L, H), Vecd(0.3 * L, H + BW), Vecd(0.4 * L, H + BW),
	Vecd(0.4 * L, H), Vecd(0.3 * L, H)
};

std::vector<Vecd> right_temperature_region
{
	Vecd(0.6 * L, H), Vecd(0.6 * L, H + BW), Vecd(0.7 * L, H + BW),
	Vecd(0.7 * L, H), Vecd(0.6 * L, H)
};

std::vector<Vecd> heat_flux_region
{
	Vecd(0.45 * L, -BW), Vecd(0.45 * L, 0), Vecd(0.55 * L, 0),
	Vecd(0.55 * L, -BW), Vecd(0.45 * L, -BW)
};

//----------------------------------------------------------------------
//	Define SPH bodies. 
//----------------------------------------------------------------------
class DiffusionBody : public MultiPolygonShape
{
public:
	explicit DiffusionBody(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::add);
	}
};

class WallBoundaryDirichlet : public MultiPolygonShape
{
public:
	explicit WallBoundaryDirichlet(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(left_temperature_region, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(right_temperature_region, ShapeBooleanOps::add);
	}
};

class WallBoundaryNeumann : public MultiPolygonShape
{
public:
	explicit WallBoundaryNeumann(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(heat_flux_region, ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	Setup diffusion material properties. 
//----------------------------------------------------------------------
class DiffusionMaterial : public DiffusionReaction<Solid>
{
public:
	DiffusionMaterial() : DiffusionReaction<Solid>({ "Phi" }, SharedPtr<NoReaction>())
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
	}
};
using DiffusionParticlesWithBoundary = DiffusionReactionParticlesWithBoundary<SolidParticles, DiffusionMaterial>;
using WallParticles = DiffusionReactionParticlesWithBoundary<SolidParticles, DiffusionMaterial>;
//----------------------------------------------------------------------
//	Application dependent initial condition. 
//----------------------------------------------------------------------
class DiffusionInitialCondition
	: public DiffusionReactionInitialConditionWithBoundary<DiffusionParticlesWithBoundary>
{
protected:
	size_t phi_;

public:
	explicit DiffusionInitialCondition(SPHBody& sph_body)
		: DiffusionReactionInitialConditionWithBoundary<DiffusionParticlesWithBoundary>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		all_species_[phi_][index_i] = 400 + 50 * (double)rand() / RAND_MAX;
	};
};

class WallBoundaryInitialCondition
	: public DiffusionReactionInitialConditionWithBoundary<WallParticles>
{
protected:
	size_t phi_;

public:
	WallBoundaryInitialCondition(SolidBody& diffusion_body) :
		DiffusionReactionInitialConditionWithBoundary<WallParticles>(diffusion_body)
	{
		phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
	}

	void update(size_t index_i, Real dt)
	{
		all_species_[phi_][index_i] = -0.0;
		if (pos_[index_i][1] > H && pos_[index_i][0] > 0.3 * L && pos_[index_i][0] < 0.4 * L)
		{
			all_species_[phi_][index_i] = left_temperature;
		}
		if (pos_[index_i][1] > H && pos_[index_i][0] > 0.6 * L && pos_[index_i][0] < 0.7 * L)
		{
			all_species_[phi_][index_i] = right_temperature;
		}
		if (pos_[index_i][1] < 0 && pos_[index_i][0] > 0.45 * L && pos_[index_i][0] < 0.55 * L)
		{
			heat_flux_[index_i] = heat_flux;
		}
	}
};

//----------------------------------------------------------------------
//	Specify diffusion relaxation method. 
//----------------------------------------------------------------------
class DiffusionBodyRelaxationWithDirichlet
	: public RelaxationOfAllDiffusionSpeciesRK2<
	RelaxationOfAllDiffusionSpeciesWithDirichlet<DiffusionParticlesWithBoundary, WallParticles>>
{
public:
	explicit DiffusionBodyRelaxationWithDirichlet(ComplexRelation& body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation) {};
	virtual ~DiffusionBodyRelaxationWithDirichlet() {};
};

class DiffusionBodyRelaxationWithNeumann
	:public RelaxationOfAllDiffusionSpeciesRK2<
	RelaxationOfAllDiffusionSpeciesWithNeumann<DiffusionParticlesWithBoundary, WallParticles>>
{
public:
	explicit DiffusionBodyRelaxationWithNeumann(ComplexRelation& body_complex_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_complex_relation) {};
	virtual ~DiffusionBodyRelaxationWithNeumann() {};
};
//----------------------------------------------------------------------
//	An observer body to measure temperature at given positions. 
//----------------------------------------------------------------------
class TemperatureObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	TemperatureObserverParticleGenerator(SPHBody& sph_body) : ObserverParticleGenerator(sph_body)
	{
		/** A line of measuring points at the middle line. */
		size_t number_of_observation_points = 11;
		Real range_of_measure = L;
		Real start_of_measure = 0;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec2d point_coordinate(0.5 * L, range_of_measure * Real(i) /
				Real(number_of_observation_points - 1) + start_of_measure);
			positions_.push_back(point_coordinate);
		}
	}
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	//---------------------------------------------------------------------- 
	SolidBody diffusion_body(sph_system, makeShared<DiffusionBody>("DiffusionBody"));
	diffusion_body.defineParticlesAndMaterial<DiffusionParticlesWithBoundary, DiffusionMaterial>();
	diffusion_body.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary_Dirichlet(sph_system, makeShared<WallBoundaryDirichlet>("WallBoundaryDirichlet"));
	wall_boundary_Dirichlet.defineParticlesAndMaterial<WallParticles, DiffusionMaterial>();
	wall_boundary_Dirichlet.generateParticles<ParticleGeneratorLattice>();

	SolidBody wall_boundary_Neumann(sph_system, makeShared<WallBoundaryNeumann>("WallBoundaryNeumann"));
	wall_boundary_Neumann.defineParticlesAndMaterial<WallParticles, DiffusionMaterial>();
	wall_boundary_Neumann.generateParticles<ParticleGeneratorLattice>();

	//----------------------------------------------------------------------
	//	Particle and body creation of temperature observers.
	//----------------------------------------------------------------------
	ObserverBody temperature_observer(sph_system, "TemperatureObserver");
	temperature_observer.generateParticles<TemperatureObserverParticleGenerator>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ComplexRelation diffusion_body_complex_Dirichlet(diffusion_body, { &wall_boundary_Dirichlet });
	ComplexRelation wall_boundary_complex_Dirichlet(wall_boundary_Dirichlet, { &diffusion_body });

	ComplexRelation diffusion_body_complex_Neumann(diffusion_body, { &wall_boundary_Neumann });
	ComplexRelation wall_boundary_complex_Neumann(wall_boundary_Neumann, { &diffusion_body });

	ContactRelation temperature_observer_contact(temperature_observer, { &diffusion_body });
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body);
	SimpleDynamics<WallBoundaryInitialCondition> setup_boundary_condition_Dirichlet(wall_boundary_Dirichlet);
	SimpleDynamics<WallBoundaryInitialCondition> setup_boundary_condition_Neumann(wall_boundary_Neumann);
	GetDiffusionTimeStepSize<DiffusionParticlesWithBoundary> get_time_step_size(diffusion_body);

	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(io_environment, sph_system.real_bodies_);
	ObservedQuantityRecording<Real> write_solid_temperature("Phi", io_environment, temperature_observer_contact);
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	DiffusionBodyRelaxationWithDirichlet temperature_relaxation_Dirichlet(diffusion_body_complex_Dirichlet);
	DiffusionBodyRelaxationWithNeumann temperature_relaxation_Neumann(diffusion_body_complex_Neumann);

	InteractionDynamics<UpdateUnitVectorNormalToBoundary<DiffusionParticlesWithBoundary, WallParticles>> update_diffusion_body_normal_vector_Dirichlet(diffusion_body_complex_Dirichlet);
	InteractionDynamics<UpdateUnitVectorNormalToBoundary<DiffusionParticlesWithBoundary, WallParticles>> update_wall_boundary_normal_vector_Dirichlet(wall_boundary_complex_Dirichlet);

	InteractionDynamics<UpdateUnitVectorNormalToBoundary<DiffusionParticlesWithBoundary, WallParticles>> update_diffusion_body_normal_vector_Neumann(diffusion_body_complex_Neumann);
	InteractionDynamics<UpdateUnitVectorNormalToBoundary<DiffusionParticlesWithBoundary, WallParticles>> update_wall_boundary_normal_vector_Neumann(wall_boundary_complex_Neumann);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	sph_system.initializeSystemConfigurations();

	setup_diffusion_initial_condition.exec();

	setup_boundary_condition_Dirichlet.exec();
	setup_boundary_condition_Neumann.exec();

	update_diffusion_body_normal_vector_Dirichlet.exec();
	update_wall_boundary_normal_vector_Dirichlet.exec();

	update_diffusion_body_normal_vector_Neumann.exec();
	update_wall_boundary_normal_vector_Neumann.exec();

	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 1;
	Real End_Time = T0;
	Real Observe_time = 0.01 * End_Time;
	Real Output_Time = 0.1 * End_Time;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TickCount::interval_t interval;
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
				if (ite % 500 == 0)
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				temperature_relaxation_Dirichlet.exec(dt);
				temperature_relaxation_Neumann.exec(dt);

				ite++;
				dt = get_time_step_size.exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}

		TickCount t2 = TickCount::now();
		write_states.writeToFile();
		write_solid_temperature.writeToFile(ite);
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TickCount::interval_t tt;
	tt = t4 - t1 - interval;

	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	std::cout << "Total physical time for computation: " << GlobalStaticVariables::physical_time_ << " seconds." << std::endl;
	return 0;
}