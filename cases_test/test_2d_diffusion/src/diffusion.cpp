/**
 * @file 	diffusion.cpp
 * @brief 	This is the first test to validate our anisotropic diffusion solver.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Chi Zhang
 */
/** SPHinXsys Library. */
#include "sphinxsys.h"
/** Namespace cite here. */
using namespace SPH;
/** Geometry parameter. */
Real L = 2.0; 	
Real H = 0.4;
/** Particle spacing. */
Real particle_spacing_ref = H / 40.0;
/** Material properties. */
Real diffusion_coff = 1.0e-4;
Real bias_diffusion_coff = 0.0;
Real alpha = Pi / 6.0;
Vec2d bias_direction(cos(alpha), sin(alpha));
/**
* @brief create a block shape
*/
std::vector<Point> CreatShape()
{
	//geometry
	std::vector<Point> shape;
	shape.push_back(Point(0.0, 0.0));
	shape.push_back(Point(0.0,  H));
	shape.push_back(Point( L,  H));
	shape.push_back(Point( L, 0.0));
	shape.push_back(Point(0.0, 0.0));

	return shape;
}
/** Define geometry and initial conditions of SPH bodies. */
class DiffusionBody : public SolidBody
{
public:
	DiffusionBody(SPHSystem &system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		std::vector<Point> body_shape = CreatShape();
		body_region_.add_geometry(new Geometry(body_shape), RegionBooleanOps::add);
		body_region_.done_modeling();
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
	: public DiffusionReactionSimple<SolidBody, SolidParticles, Solid>
{
protected:
	size_t phi_;

	void Update(size_t index_particle_i, Real dt) override
	{
		BaseParticleData& base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		DiffusionReactionData& diffusion_data_i = particles_->diffusion_reaction_data_[index_particle_i];

        if(0.45 <= base_particle_data_i.pos_n_[0] && base_particle_data_i.pos_n_[0] <= 0.55)
		{
			diffusion_data_i.species_n_[phi_] = 1.0;
		}
		if(base_particle_data_i.pos_n_[0] >= 1.0)
		{
			diffusion_data_i.species_n_[phi_] = exp(-2500.0 * ((base_particle_data_i.pos_n_[0] - 1.5)
					* (base_particle_data_i.pos_n_[0] - 1.5)));
		}
	};
public:
	DiffusionBodyInitialCondition(SolidBody* diffusion_body)
		: DiffusionReactionSimple<SolidBody, SolidParticles, Solid>(diffusion_body) {
		phi_ = material_->getSpeciesIndexMap()["Phi"];
	};
};
/** Set diffusion relaxation. */
class DiffusionBodyRelaxation
	: public RelaxationOfAllDifussionSpeciesRK2<SolidBody, SolidParticles, Solid>
{
public:
	DiffusionBodyRelaxation(SolidBody* body)
		: RelaxationOfAllDifussionSpeciesRK2<SolidBody, SolidParticles, Solid>(body) {
	};
	virtual ~DiffusionBodyRelaxation() {};
};

/** The main program. */
int main()
{
	/** Build up context -- a SPHSystem. */
	SPHSystem system(Vec2d(0.0, 0.0), Vec2d(L, H), particle_spacing_ref, 4);
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Configuration of materials, crate particle container and diffusion body. */
	DiffusionBody *diffusion_body  =  new DiffusionBody(system, "DiffusionBody", 0, ParticlesGeneratorOps::lattice);
	DiffusionBodyMaterial *diffusion_body_material = new DiffusionBodyMaterial();
	DiffusionReactionParticles<SolidParticles, Solid>	diffusion_body_particles(diffusion_body, diffusion_body_material);
	/** Set body contact map. */
	SPHBodyTopology body_topology = { { diffusion_body, {  } }};
	system.SetBodyTopology(&body_topology);

	/**
	 * The main dynamics algorithm is defined start here.
	 */
	/** Case setup */
	DiffusionBodyInitialCondition setup_diffusion_initial_condition(diffusion_body);
	/** Corrected strong configuration for diffusion body. */	
	solid_dynamics::CorrectConfiguration 			correct_configuration(diffusion_body);
	/** Time step size caclutation. */
	GetDiffusionTimeStepSize<SolidBody, SolidParticles, Solid> get_time_step_size(diffusion_body);
	/** Diffusion process for diffusion body. */
	DiffusionBodyRelaxation 			diffusion_relaxation(diffusion_body);
	/** Periodic BCs. */
	PeriodicConditionInAxisDirection 					periodic_condition_y(diffusion_body, 1);
	/**
	 * @brief simple input and outputs.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToVtu 				write_states(in_output, system.real_bodies_);

	/** Pre-simultion*/
	system.InitializeSystemCellLinkedLists();
	system.InitializeSystemConfigurations();
	setup_diffusion_initial_condition.exec();
	periodic_condition_y.parallel_exec();
	diffusion_body->BuildInnerConfiguration();
	correct_configuration.parallel_exec();
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
		Real integeral_time = 0.0;
		while (integeral_time < Output_Time) 
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
				integeral_time += dt;
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
