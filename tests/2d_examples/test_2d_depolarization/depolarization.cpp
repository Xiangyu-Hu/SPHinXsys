/**
 * @file 	depolarization.cpp
 * @brief 	This is the first test to validate our PED-ODE solver for solving
 * 			electrophysiology mono-domain model closed by a physiology reaction.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 50.0;
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));
// observer location
StdVec<Vecd> observation_location = {Vecd(0.3, 0.7)};
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coff = 1.0;
Real bias_coff = 0.0;
Vec2d fiber_direction(1.0, 0.0);
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.15;
Real b = 0.0;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.04;
Real k_a = 0.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
class MuscleBlock : public MultiPolygonShape
{
public:
	explicit MuscleBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		std::vector<Vecd> shape;
		shape.push_back(Vecd(0.0, 0.0));
		shape.push_back(Vecd(0.0, H));
		shape.push_back(Vecd(L, H));
		shape.push_back(Vecd(L, 0.0));
		shape.push_back(Vecd(0.0, 0.0));
		multi_polygon_.addAPolygon(shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
class DepolarizationInitialCondition
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

public:
	explicit DepolarizationInitialCondition(SPHBody &sph_body)
		: electro_physiology::ElectroPhysiologyInitialCondition(sph_body)
	{
		voltage_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Voltage"];
	};

	void update(size_t index_i, Real dt)
	{
		all_species_[voltage_][index_i] = exp(-4.0 * ((pos_[index_i][0] - 1.0) * (pos_[index_i][0] - 1.0) + pos_[index_i][1] * pos_[index_i][1]));
	};
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	IOEnvironment io_environment(system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody muscle_body(system, makeShared<MuscleBlock>("MuscleBlock"));
	SharedPtr<AlievPanfilowModel> muscle_reaction_model_ptr = makeShared<AlievPanfilowModel>(k_a, c_m, k, a, b, mu_1, mu_2, epsilon);
	muscle_body.defineParticlesAndMaterial<ElectroPhysiologyParticles, MonoFieldElectroPhysiology>(
		muscle_reaction_model_ptr, TypeIdentity<DirectionalDiffusion>(), diffusion_coff, bias_coff, fiber_direction);
	muscle_body.generateParticles<ParticleGeneratorLattice>();

	ObserverBody voltage_observer(system, "VoltageObserver");
	voltage_observer.generateParticles<ObserverParticleGenerator>(observation_location);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation muscle_body_inner_relation(muscle_body);
	ContactRelation voltage_observer_contact_relation(voltage_observer, {&muscle_body});
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simulation.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	SimpleDynamics<DepolarizationInitialCondition> initialization(muscle_body);
	InteractionWithUpdate<CorrectedConfigurationInner> correct_configuration(muscle_body_inner_relation);
	electro_physiology::GetElectroPhysiologyTimeStepSize get_time_step_size(muscle_body);
	// Diffusion process for diffusion body.
	electro_physiology::ElectroPhysiologyDiffusionInnerRK2 diffusion_relaxation(muscle_body_inner_relation);
	// Solvers for ODE system or reactions
	electro_physiology::ElectroPhysiologyReactionRelaxationForward reaction_relaxation_forward(muscle_body);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward reaction_relaxation_backward(muscle_body);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>>
		write_recorded_voltage("Voltage", io_environment, voltage_observer_contact_relation);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	initialization.exec();
	correct_configuration.exec();
	//----------------------------------------------------------------------
	//	Initial states output.
	//----------------------------------------------------------------------
	write_states.writeToFile(0);
	write_recorded_voltage.writeToFile(0);
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 16.0;
	Real end_time = T0;
	Real output_interval = 0.5;		  /**< Time period for output */
	Real Dt = 0.01 * output_interval; /**< Time period for data observing */
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_interval)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				if (ite % 1000 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: "
							  << dt << "\n";
				}
				/**Strang splitting method. */
				reaction_relaxation_forward.exec(0.5 * dt);
				diffusion_relaxation.exec(dt);
				reaction_relaxation_backward.exec(0.5 * dt);

				ite++;
				dt = get_time_step_size.exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_recorded_voltage.writeToFile(ite);
		}

		TickCount t2 = TickCount::now();
		write_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();

	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	write_recorded_voltage.testResult();

	return 0;
}
