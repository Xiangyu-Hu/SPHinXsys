/**
 * @file 	diffusion.cpp
 * @brief 	This is the first test to validate our anisotropic diffusion solver.
 * @author 	Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library
using namespace SPH; //Namespace cite here
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
Real bias_diffusion_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
std::vector<Vecd> createShape()
{
	//geometry
	std::vector<Vecd> shape;
	shape.push_back(Vecd(0.0, 0.0));
	shape.push_back(Vecd(0.0,  H));
	shape.push_back(Vecd( L,  H));
	shape.push_back(Vecd( L, 0.0));
	shape.push_back(Vecd(0.0, 0.0));

	return shape;
}
//----------------------------------------------------------------------
//	Define SPH bodies. 
//----------------------------------------------------------------------
class DiffusionBody : public SolidBody
{
public:
	DiffusionBody(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		std::vector<Vecd> body_shape = createShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(body_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Setup diffusion material properties. 
//----------------------------------------------------------------------
class DiffusionBodyMaterial
	: public DiffusionReactionMaterial<SolidParticles, Solid>
{
public:
	DiffusionBodyMaterial()
		: DiffusionReactionMaterial<SolidParticles, Solid>()
	{
		insertASpecies("Phi");
		assignDerivedMaterialParameters();
		initializeDiffusion();
	}

	/** Initialize diffusion reaction material. */
	virtual void initializeDiffusion() override {
		DirectionalDiffusion* phi_diffusion
			= new DirectionalDiffusion(species_indexes_map_["Phi"], species_indexes_map_["Phi"],
				diffusion_coff, bias_diffusion_coff, bias_direction);
		species_diffusion_.push_back(phi_diffusion);
	};
};
//----------------------------------------------------------------------
//	Application dependent initial condition. 
//----------------------------------------------------------------------
class DiffusionBodyInitialCondition 
	: public DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>
{
protected:
	size_t phi_;

	void Update(size_t index_i, Real dt) override
	{

        if(pos_n_[index_i][0] >= 0.45 && pos_n_[index_i][0] <= 0.55)
		{
			species_n_[phi_][index_i] = 1.0;
		}
		if(pos_n_[index_i][0] >= 1.0)
		{
			species_n_[phi_][index_i] = exp(-2500.0 
				* ((pos_n_[index_i][0] - 1.5) * (pos_n_[index_i][0] - 1.5)));
		}
	};
public:
	DiffusionBodyInitialCondition(SolidBody* diffusion_body) : 
		DiffusionReactionInitialCondition<SolidBody, SolidParticles, Solid>(diffusion_body) {
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
	DiffusionBodyRelaxation(BodyRelationInner* body_inner_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation) {
	};
	virtual ~DiffusionBodyRelaxation() {};
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
	DiffusionBody *diffusion_body  =  new DiffusionBody(sph_system, "DiffusionBody");
	DiffusionBodyMaterial *diffusion_body_material = new DiffusionBodyMaterial();
	DiffusionReactionParticles<SolidParticles, Solid>	diffusion_body_particles(diffusion_body, diffusion_body_material);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner* diffusion_body_inner_relation = new BodyRelationInner(diffusion_body);
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simultion.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	DiffusionBodyInitialCondition setup_diffusion_initial_condition(diffusion_body);
	/** Corrected strong configuration for diffusion body. */	
	solid_dynamics::CorrectConfiguration 			correct_configuration(diffusion_body_inner_relation);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<SolidBody, SolidParticles, Solid> get_time_step_size(diffusion_body);
	/** Diffusion process for diffusion body. */
	DiffusionBodyRelaxation 			diffusion_relaxation(diffusion_body_inner_relation);
	/** Periodic BCs. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList	periodic_condition_y(diffusion_body, yAxis);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtu 				write_states(in_output, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	sph_system.initializeSystemCellLinkedLists();
	periodic_condition_y.update_cell_linked_list_.parallel_exec();
	sph_system.initializeSystemConfigurations();
	correct_configuration.parallel_exec();
	setup_diffusion_initial_condition.exec();
	/** Output global basic parameters. */
	write_states.writeToFile(0);

	int ite 				= 0;
	Real T0 				= 1.0;
	Real End_Time 			= T0;
	Real Output_Time 	    = 0.1 * End_Time;
	Real Observe_time 		= 0.1 * Output_Time;
	Real dt		 			= 0.0;
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

				ite++;
				dt = get_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}

		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
