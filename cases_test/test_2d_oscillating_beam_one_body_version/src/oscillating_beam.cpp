/* ---------------------------------------------------------------------------*
*            SPHinXsys: 2D oscilation beam example-one body version           *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases, also the first case for            *
* understanding SPH method for solid similation.                              *
* In this case, the constriant of the beam is implemented with                *
* internal constrianed subregion.                                             *
* ----------------------------------------------------------------------------*/
/**
  * @brief 	SPHinXsys Library.
  */
#include "sphinxsys.h"
/**
 * @brief Namespace cite here.
 */
using namespace SPH;

//------------------------------------------------------------------------------
//global parameters for the case
//------------------------------------------------------------------------------

//for geometry
Real PL = 0.2; 	//beam lenght
Real PH = 0.02; //for thick plate; =0.01 for thin plate
Real SL = 0.06; //depth of the insert
//particle spacing, at least three particles
Real particle_spacing_ref = PH / 10.0;
Real BW = particle_spacing_ref * 4; 	//boundary width

//for material properties of the beam
Real rho0_s = 1.0e3; 			//reference density
Real Youngs_modulus = 2.0e6;	//reference Youngs modulus
Real poisson = 0.3975; 			//Poisson ratio

//for initial condition
Real initial_pressure = 0.0;
//for initial velocity
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl)*sinh(kl) - sin(kl)*cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * pi);

/**
* @brief define geometry and initial conditions of SPH bodies
*/
/**
* @brief create a beam base shape
*/
std::vector<Point> CreatBeamBaseShape()
{
	//geometry
	std::vector<Point> beam_base_shape;
	beam_base_shape.push_back(Point(-SL - BW, -PH / 2 - BW));
	beam_base_shape.push_back(Point(-SL - BW, PH / 2 + BW));
	beam_base_shape.push_back(Point(0.0, PH / 2 + BW));
	beam_base_shape.push_back(Point(0.0, -PH / 2 - BW));
	beam_base_shape.push_back(Point(-SL - BW, -PH / 2 - BW));
	
	return beam_base_shape;
}
/**
* @brief create a beam shape
*/
std::vector<Point> CreatBeamShape()
{
	std::vector<Point> beam_shape;
	beam_shape.push_back(Point(-SL, -PH / 2));
	beam_shape.push_back(Point(-SL, PH / 2));
	beam_shape.push_back(Point(PL, PH / 2));
	beam_shape.push_back(Point(PL, -PH / 2));
	beam_shape.push_back(Point(-SL, -PH / 2));

	return beam_shape;
}
/**
* @brief define the beam body
*/
class Beam : public SolidBody
{
public:
	Beam(SPHSystem &system, string body_name, ElasticSolid &material, ElasticSolidParticles 
						&elastic_particles, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, material, elastic_particles, refinement_level, op)
	{
		/** Geometry defination. */
		std::vector<Point> beam_base_shape = CreatBeamBaseShape();
		Geometry * beam_base_gemetry = new Geometry(beam_base_shape);
		body_region_.add_geometry(beam_base_gemetry, RegionBooleanOps::add);

		std::vector<Point> beam_shape = CreatBeamShape();
		Geometry * beam_gemetry = new Geometry(beam_shape);
		body_region_.add_geometry(beam_gemetry, RegionBooleanOps::add);

		/** Finish the region modeling. */
		body_region_.done_modeling();
	}
};
/**
 * application dependent initial condition 
 */
class BeamInitialCondition
	: public solid_dynamics::ElasticSolidDynamicsInitialCondition
{
public:
	BeamInitialCondition(SolidBody *beam)
		: solid_dynamics::ElasticSolidDynamicsInitialCondition(beam) {};
protected:
	void Update(size_t index_particle_i, Real dt) override {
		/** first set all particle at rest*/
		solid_dynamics::ElasticSolidDynamicsInitialCondition::Update(index_particle_i, dt);

		/** initial velocity profile */
		BaseParticleData &base_particle_data_i = particles_->base_particle_data_[index_particle_i];
		SolidParticleData &solid_body_data_i = particles_->solid_body_data_[index_particle_i];

		Real x = solid_body_data_i.pos_0_[0] / PL;
		if (x > 0.0) {
			base_particle_data_i.vel_n_[1]
				= vf * material_->c_0_*(M*(cos(kl*x) - cosh(kl*x)) - N * (sin(kl*x) - sinh(kl*x))) / Q;
		}
	};
};
/**
* @brief define the beam base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class BeamBase : public SolidBodyPart
{
public:
	BeamBase(SolidBody *solid_body, string constrianed_region_name)
		: SolidBodyPart(solid_body, constrianed_region_name)
	{
		/* Geometry defination */
		std::vector<Point> beam_base_shape = CreatBeamBaseShape();
		Geometry * beam_base_gemetry = new Geometry(beam_base_shape);
		soild_body_part_region_.add_geometry(beam_base_gemetry, RegionBooleanOps::add);
		std::vector<Point> beam_shape = CreatBeamShape();
		Geometry * beam_gemetry = new Geometry(beam_shape);
		soild_body_part_region_.add_geometry(beam_gemetry, RegionBooleanOps::sub);

		/** Finish the region modeling. */
		soild_body_part_region_.done_modeling();

		//tag the constrained particle
		TagBodyPartParticles();
	}
};

//define an observer body
class BeamObserver : public ObserverLagrangianBody
{
public:
	BeamObserver(SPHSystem &system, string body_name, ObserverParticles 
			   &observer_particles, int refinement_level, ParticlesGeneratorOps op)
		: ObserverLagrangianBody(system, body_name, observer_particles, refinement_level, op)
	{
		body_input_points_volumes_.push_back(make_pair(Point(PL, 0.0), 0.0));
	}
};

//------------------------------------------------------------------------------
//the main program
//------------------------------------------------------------------------------

int main()
{

	//build up context -- a SPHSystem
	SPHSystem system(Vec2d(-SL - BW, -PL / 2.0), 
		Vec2d(PL + 3.0*BW, PL / 2.0), particle_spacing_ref);

	//Configuration of soild materials
	ElasticSolid solid_material("ElasticSolid", rho0_s, Youngs_modulus, poisson, 0.0);

	//creat a particle cotainer the elastic body
	ElasticSolidParticles beam_particles("BeamBody");
	//the osillating beam
	Beam *beam_body = 
		new Beam(system, "BeamBody", solid_material, beam_particles, 0, ParticlesGeneratorOps::lattice);

	//create a observer particle container
	ObserverParticles observer_particles("BeamObserver");
	BeamObserver *beam_observer 
		= new BeamObserver(system, "BeamObserver", observer_particles, 1, ParticlesGeneratorOps::direct);

	//set body contact map
	//the contact map gives the data conntections between the bodies
	//basically the the range of bidies to build neighbor particle lists
	SPHBodyTopology body_topology 
		= { { beam_body, {} }, { beam_observer,{ beam_body} } };
	system.SetBodyTopology(&body_topology);

	//setting up the simulation
	system.SetupSPHSimulation();


	//-----------------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-----------------------------------------------------------------------------
	/** initial condition */
	BeamInitialCondition beam_initial_velocity(beam_body);
	//corrected strong configuration	
	solid_dynamics::CorrectConfiguration
		beam_corrected_configuration_in_strong_form(beam_body, {});

	//time step size caclutation
	solid_dynamics::GetAcousticTimeStepSize computing_time_step_size(beam_body);

	//stress relaxation for the beam
	solid_dynamics::StressRelaxationFirstStep
		stress_relaxation_first_step(beam_body);
	solid_dynamics::StressRelaxationSecondStep
		stress_relaxation_second_step(beam_body);

	/**
	 * @brief Constrain a solid body part
	 */
	solid_dynamics::ConstrainSolidBodyRegion 
		constrain_beam_base(beam_body, new BeamBase(beam_body, "BeamBase"));

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	WriteBodyStatesToVtu write_beam_states(in_output, system.real_bodies_);
	WriteObservedElasticDisplacement
		write_beam_tip_displacement(in_output, beam_observer, { beam_body });


	/**
	 * @brief Setup goematrics and initial conditions
	 */
	beam_initial_velocity.exec();
	beam_corrected_configuration_in_strong_form.parallel_exec();

	//-----------------------------------------------------------------------------
	//from here the time stepping begines
	//-----------------------------------------------------------------------------
	//starting time zero
	GlobalStaticVariables::physical_time_ = 0.0;
	write_beam_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);

	int ite = 0;
	Real T0 = 1.0;
	Real End_Time = T0;
	//time step size for oupt file
	Real D_Time = 0.01*T0;
	Real Dt = 0.01*D_Time;			//default advection time step sizes
	Real dt = 0.0; 					//default accoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	//computation loop starts 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integeral_time = 0.0;
		//integrate time (loop) until the next output time
		while (integeral_time < D_Time) {

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {

				if (ite % 100 == 0) {
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				stress_relaxation_first_step.parallel_exec(dt);
				constrain_beam_base.parallel_exec(dt);
				stress_relaxation_second_step.parallel_exec(dt);

				ite++;
				dt = computing_time_step_size.parallel_exec();
				relaxation_time += dt;
				integeral_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}

		write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);

		tick_count t2 = tick_count::now();
		write_beam_states.WriteToFile(GlobalStaticVariables::physical_time_);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	cout << "Total wall time for computation: " << tt.seconds() << " seconds." << endl;

	return 0;
}
