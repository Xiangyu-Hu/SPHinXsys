/**
 * @file 	diffusion.cpp
 * @brief 	This is the first test to validate our anisotropic diffusion solver.
 * @author 	Chi Zhang and Xiangyu Hu
 */
/** SPHinXsys Library. */
#include "sphinxsys.h"
/** Namespace cite here. */
using namespace SPH;
/** Geometry parameter. */
Real L = 2.0; 	
Real H = 0.4;
/** Particle spacing. */
Real relosution_ref = H / 40.0;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));

/** Material properties. */
Real diffusion_coff = 1.0e-4;
Real bias_diffusion_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
/**
* @brief create a block shape
*/
std::vector<Vecd> CreatShape()
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
/** Define geometry and initial conditions of SPH bodies. */
class DiffusionBody : public SolidBody
{
public:
	DiffusionBody(SPHSystem &system, string body_name)
		: SolidBody(system, body_name)
	{
		std::vector<Vecd> body_shape = CreatShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(body_shape, ShapeBooleanOps::add);
	}
};
/**
 * Setup diffusion material properties
 */
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
/**
 * application dependent initial condition 
 */
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
/** Set diffusion relaxation. */
class DiffusionBodyRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2<SolidBody, SolidParticles, Solid,
	RelaxationOfAllDiffussionSpeciesInner<SolidBody, SolidParticles, Solid>, 
	InnerBodyRelation>
{
public:
	DiffusionBodyRelaxation(InnerBodyRelation* body_inner_relation)
		: RelaxationOfAllDiffusionSpeciesRK2(body_inner_relation) {
	};
	virtual ~DiffusionBodyRelaxation() {};
};

/** The main program. */
int main()
{
	/** Build up context -- a SPHSystem. */
	SPHSystem system(system_domain_bounds, relosution_ref);
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Configuration of materials, crate particle container and diffusion body. */
	DiffusionBody *diffusion_body  =  new DiffusionBody(system, "DiffusionBody");
	DiffusionBodyMaterial *diffusion_body_material = new DiffusionBodyMaterial();
	DiffusionReactionParticles<SolidParticles, Solid>	diffusion_body_particles(diffusion_body, diffusion_body_material);

	/** topology */
	InnerBodyRelation* diffusion_body_inner_relation = new InnerBodyRelation(diffusion_body);

	/**
	 * The main dynamics algorithm is defined start here.
	 */
	/** Case setup */
	DiffusionBodyInitialCondition setup_diffusion_initial_condition(diffusion_body);
	/** Corrected strong configuration for diffusion body. */	
	solid_dynamics::CorrectConfiguration 			correct_configuration(diffusion_body_inner_relation);
	/** Time step size calculation. */
	GetDiffusionTimeStepSize<SolidBody, SolidParticles, Solid> get_time_step_size(diffusion_body);
	/** Diffusion process for diffusion body. */
	DiffusionBodyRelaxation 			diffusion_relaxation(diffusion_body_inner_relation);
	/** Periodic BCs. */
	PeriodicConditionInAxisDirectionUsingCellLinkedList	periodic_condition_y(diffusion_body, 1);
	/**
	 * @brief simple input and outputs.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToVtu 				write_states(in_output, system.real_bodies_);

	/** Pre-simultion*/
	system.initializeSystemCellLinkedLists();
	periodic_condition_y.update_cell_linked_list_.parallel_exec();
	system.initializeSystemConfigurations();
	correct_configuration.parallel_exec();
	setup_diffusion_initial_condition.exec();
	/** Output global basic parameters. */
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite 				= 0;
	Real T0 				= 1.0;
	Real End_Time 			= T0;
	Real Output_Time 	    = 0.1 * End_Time;
	Real Observe_time 		= 0.1 * Output_Time;
	Real dt		 			= 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main loop starts here. */ 
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
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				diffusion_relaxation.parallel_exec(dt);

				ite++;
				dt = get_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				write_states.WriteToFile(GlobalStaticVariables::physical_time_);
			}
		}

		tick_count t2 = tick_count::now();
		write_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
