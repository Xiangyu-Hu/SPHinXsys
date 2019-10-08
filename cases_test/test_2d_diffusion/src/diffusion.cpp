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
/** Thickness of the surounding boundary. */
Real BW = 4.0 * particle_spacing_ref;
/** Material properties. */
Real rho_0 = 1.0; 	
Real a_0[4] = {1.0, 0.0, 0.0, 0.0};
Real b_0[4] = {1.0, 0.0, 0.0, 0.0};
Vec2d d_0(1.0e-4, 0.0);
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
 * application dependent initial condition 
 */
class DiffusionInitialCondition
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
public:
	DiffusionInitialCondition(SolidBody *diffusion)
		: electro_physiology::ElectroPhysiologyInitialCondition(diffusion) {};
protected:
	void Update(size_t index_particle_i, Real dt) override 
	{
		/** first set all particle at rest*/
		electro_physiology::ElectroPhysiologyInitialCondition::Update(index_particle_i, dt);

		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		MuscleParticleData &muscle_particle_data_i = particles_->muscle_body_data_[index_particle_i];

        if(0.45 <= base_particle_data_i.pos_n_[0] && base_particle_data_i.pos_n_[0] <= 0.55)
		{
			muscle_particle_data_i.voltage_n_ = 1.0;
		}
		if(base_particle_data_i.pos_n_[0] >= 1.0)
		{
			muscle_particle_data_i.voltage_n_ = exp(-2500 * ((base_particle_data_i.pos_n_[0] - 1.5) 
					* (base_particle_data_i.pos_n_[0] - 1.5)));
		}
	};
};
/** The main program. */
int main()
{
	/** Build up context -- a SPHSystem. */
	SPHSystem system(Vec2d(- BW, - BW), Vec2d(L + BW, H + BW), particle_spacing_ref);
	GlobalStaticVariables::physical_time_ = 0.0;
	/** Configuration of materials, crate particle container and diffusion body. */
	DiffusionBody *diffusion_body  =  new DiffusionBody(system, "DiffusionBody", 0, ParticlesGeneratorOps::lattice);
	Muscle 						material("Muscle", diffusion_body, a_0, b_0,d_0, rho_0, 1.0);
	MuscleParticles 			particles(diffusion_body);
	/** Set body contact map. */
	SPHBodyTopology body_topology = { { diffusion_body, {  } }};
	system.SetBodyTopology(&body_topology);
	/**
	 * @brief 	Simulation set up.
	 */
	system.SetupSPHSimulation();
	/**
	 * The main dynamics algorithm is defined start here.
	 */
	/** Case setup */
	DiffusionInitialCondition setup_diffusion_initial_condition(diffusion_body);
	/** Corrected strong configuration for diffusion body. */	
	electro_physiology::CorrectConfiguration 				correct_configuration(diffusion_body);
	/** Time step size caclutation. */
	electro_physiology::getDiffusionTimeStepSize 			get_time_step_size(diffusion_body);
	/** Diffusion process for diffusion body. */
	electro_physiology::DiffusionRelaxation 				diffusion_relaxation(diffusion_body);
	/**
	 * @brief simple input and outputs.
	 */
	In_Output 							in_output(system);
	WriteBodyStatesToPlt 				write_states(in_output, system.real_bodies_);
	/** Pre-simultion*/
	setup_diffusion_initial_condition.exec();
	correct_configuration.parallel_exec();
	/** Output global basic parameters. */
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite 		= 0;
	Real T0 		= 1.0;
	Real End_Time 	= T0;
	Real D_Time 	= 0.1 * End_Time;
	Real Dt 		= 0.001 * D_Time;
	Real dt		 	= 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main loop starts here. */ 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		while (integeral_time < D_Time) 
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) 
			{
				if (ite % 100 == 0) 
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
