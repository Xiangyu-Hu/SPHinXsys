/**
 * @file 	depolarization.cpp
 * @brief 	This is the first test to validate our PED-ODE solver for solving
 * 			electrophysiology monodomain model closed by a physiology reaction.
 * @author 	Chi Zhang and Xiangyu Hu
 */

/** SPHinXsys Library. */
#include "sphinxsys.h"
/** Namespace cite here. */
using namespace SPH;
/** Geometry parameter. */
Real L = 1.0; 	
Real H = 1.0;
/** Reference particle spacing. */
Real resolution_ref = H / 50.0;
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(0.0, 0.0), Vec2d(L, H));

/** Electrophysiology parameters. */
Real diffusion_coff_ = 1.0;
Real bias_diffusion_coff_ = 0.0;
Vec2d fiber_direction(1.0, 0.0);
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.04;
Real k_a = 0.0;
/** create a block shape */
std::vector<Vecd> CreatShape()
{
	std::vector<Vecd> shape;
	shape.push_back(Vecd(0.0, 0.0));
	shape.push_back(Vecd(0.0,  H));
	shape.push_back(Vecd(L ,  H));
	shape.push_back(Vecd(L, 0.0));
	shape.push_back(Vecd(0.0, 0.0));
	return shape;
}
/** 
 * Define geometry and initial conditions of SPH bodies. 
 * */
class MuscleBody : public SolidBody
{
public:
	MuscleBody(SPHSystem& system, std::string body_name) : SolidBody(system, body_name)
	{

		std::vector<Vecd> block_shape = CreatShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(block_shape, ShapeBooleanOps::add);
	}
};
/**
 * Voltage observer body definition.
 */
class VoltageObserver : public FictitiousBody
{
public:
	VoltageObserver(SPHSystem &system, std::string body_name) : FictitiousBody(system, body_name)
	{
		/** postion and volume. */
		body_input_points_volumes_.push_back(std::make_pair(Vecd(0.3, 0.7), 0.0));
	}
};

/**
 * Setup electro_physiology reaction properties
 */
class MuscleReactionModel : public AlievPanfilowModel
{
public:
	MuscleReactionModel() : AlievPanfilowModel()
	{
		/** Basic reaction parameters*/
		k_a_ = k_a;
		c_m_ = c_m;
		k_ = k;
		a_ = a;
		mu_1_ = mu_1;
		mu_2_ = mu_2;
		epsilon_ = epsilon;

		/** Compute the derived material parameters*/
		assignDerivedReactionParameters();
	}
};

/**
 * Setup material properties of myocardium
 */
class MyocardiumMuscle
 	: public MonoFieldElectroPhysiology
{
 public:
 	MyocardiumMuscle(ElectroPhysiologyReaction* electro_physiology_reaction)
		: MonoFieldElectroPhysiology(electro_physiology_reaction)
	{
		/** Basic material parameters*/
		diff_cf_ = diffusion_coff_;
		bias_diff_cf_ = bias_diffusion_coff_;
		bias_direction_ = fiber_direction;

		/** Compute the derived material parameters. */
		assignDerivedMaterialParameters();
		/** Create the vector of diffusion properties. */
		initializeDiffusion();
	}
};
 /**
 * application dependent initial condition 
 */
class DepolarizationInitialCondition
	: public electro_physiology::ElectroPhysiologyInitialCondition
{
protected:
	size_t voltage_;

	void Update(size_t index_i, Real dt) override
	{
		species_n_[voltage_][index_i] = exp(-4.0 * ((pos_n_[index_i][0] - 1.0)
				* (pos_n_[index_i][0] - 1.0) + pos_n_[index_i][1] *	pos_n_[index_i][1]));
	};
public:
	DepolarizationInitialCondition(SolidBody* muscle)
		: electro_physiology::ElectroPhysiologyInitialCondition(muscle) {
		voltage_ = material_->SpeciesIndexMap()["Voltage"];
	};
};

/** 
 * The main program. 
 */
int main()
{
	/** 
	 * Build up context -- a SPHSystem. 
	 */
	SPHSystem system(system_domain_bounds, resolution_ref);
		GlobalStaticVariables::physical_time_ = 0.0;
	/** output environment. */
	In_Output in_output(system);
	/**
	 * Configuration of materials, crate particle container and muscle body. 
	 */
	MuscleBody *muscle_body  =  new MuscleBody(system, "MuscleBody");
	MuscleReactionModel *muscle_reaction_model = new MuscleReactionModel();
	MyocardiumMuscle 	*myocardium_muscle = new MyocardiumMuscle(muscle_reaction_model);
	ElectroPhysiologyParticles 		myocardium_muscle_particles(muscle_body, myocardium_muscle);
	/**
	 * Particle and body creation of fluid observer.
	 */
	VoltageObserver *voltage_observer = new VoltageObserver(system, "VoltageObserver");
	BaseParticles 					observer_particles(voltage_observer);

	/** topology */
	InnerBodyRelation* muscle_body_inner_relation = new InnerBodyRelation(muscle_body);
	ContactBodyRelation* voltage_observer_contact_relation = new ContactBodyRelation(voltage_observer, { muscle_body });

	/**
	 * The main dynamics algorithm is defined start here.
	 */
	/**
	 * Initialization for electrophysiology computation. 
	 */	
	DepolarizationInitialCondition 								initialization(muscle_body);
	/** 
	 * Corrected strong configuration. 
	 */	
	solid_dynamics::CorrectConfiguration 					correct_configuration(muscle_body_inner_relation);
	/** 
	 * Time step size calculation. 
	 */
	electro_physiology::GetElectroPhysiologyTimeStepSize 		get_time_step_size(muscle_body);
	/** 
	 * Diffusion process for diffusion body. 
	 */
	electro_physiology::ElectroPhysiologyDiffusionRelaxationInner 		diffusion_relaxation(muscle_body_inner_relation);
	/** 
	 * Solvers for ODE system 
	 */
	electro_physiology::ElectroPhysiologyReactionRelaxationForward 		reaction_relaxation_forward(muscle_body);
	electro_physiology::ElectroPhysiologyReactionRelaxationBackward 	reaction_relaxation_backward(muscle_body);
	/** Outputs. */
	/**
	 * Simple input and outputs.
	 */
	WriteBodyStatesToVtu 				write_states(in_output, system.real_bodies_);
	WriteAnObservedQuantity<indexScalar, Real>
		write_recorded_voltage("Voltage", in_output, voltage_observer_contact_relation);

	/** 
	 * Pre-simultion. 
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	initialization.exec();
	correct_configuration.parallel_exec();
	/** 
	 * Output global basic parameters. 
	 */
	write_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_recorded_voltage.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite 		= 0;
	Real T0 		= 16.0;
	Real End_Time 	= T0;
	Real D_Time 	= 0.5; 				/**< Time period for output */
	Real Dt 		= 0.01 * D_Time;	/**< Time period for data observing */
	Real dt		 	= 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/** Main loop starts here. */ 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < D_Time) 
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
				reaction_relaxation_forward.parallel_exec(0.5 * dt);
				diffusion_relaxation.parallel_exec(dt);
				reaction_relaxation_backward.parallel_exec(0.5 * dt);

				ite++;
				dt = get_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
			write_recorded_voltage.WriteToFile(GlobalStaticVariables::physical_time_);
		}

		tick_count t2 = tick_count::now();
		write_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
