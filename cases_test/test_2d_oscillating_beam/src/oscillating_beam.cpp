/* ---------------------------------------------------------------------------*
*            SPHinXsys: 2D oscillation beam example-one body version           *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases, also the first case for            *
* understanding SPH method for solid simulation.                              *
* In this case, the constraint of the beam is implemented with                *
* internal constrained subregion.                                             *
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
Real PL = 0.2; 	//beam length
Real PH = 0.02; //for thick plate; =0.01 for thin plate
Real SL = 0.06; //depth of the insert
//particle spacing, at least three particles
Real resolution_ref = PH / 10.0;
Real BW = resolution_ref * 4; 	//boundary width
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
	Vec2d(PL + 3.0 * BW, PL / 2.0));

//for material properties of the beam
Real rho0_s = 1.0e3; 			//reference density
Real Youngs_modulus = 2.0e6;	//reference Youngs modulus
Real poisson = 0.3975; 			//Poisson ratio

//for initial condition on velocity
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl)*sinh(kl) - sin(kl)*cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * Pi);

/**
* @brief define geometry and initial conditions of SPH bodies
*/
/**
* @brief create a beam base shape
*/
std::vector<Vecd> CreatBeamBaseShape()
{
	//geometry
	std::vector<Vecd> beam_base_shape;
	beam_base_shape.push_back(Vecd(-SL - BW, -PH / 2 - BW));
	beam_base_shape.push_back(Vecd(-SL - BW, PH / 2 + BW));
	beam_base_shape.push_back(Vecd(0.0, PH / 2 + BW));
	beam_base_shape.push_back(Vecd(0.0, -PH / 2 - BW));
	beam_base_shape.push_back(Vecd(-SL - BW, -PH / 2 - BW));
	
	return beam_base_shape;
}
/**
* @brief create a beam shape
*/
std::vector<Vecd> CreatBeamShape()
{
	std::vector<Vecd> beam_shape;
	beam_shape.push_back(Vecd(-SL, -PH / 2));
	beam_shape.push_back(Vecd(-SL, PH / 2));
	beam_shape.push_back(Vecd(PL, PH / 2));
	beam_shape.push_back(Vecd(PL, -PH / 2));
	beam_shape.push_back(Vecd(-SL, -PH / 2));

	return beam_shape;
}
/**
* @brief define the beam body
*/
class Beam : public SolidBody
{
public:
	Beam(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> beam_base_shape = CreatBeamBaseShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		std::vector<Vecd> beam_shape = CreatBeamShape();
		body_shape_->addAPolygon(beam_shape, ShapeBooleanOps::add);
	}
};
/**
 * @brief Define beam material.
 */
class BeamMaterial : public LinearElasticSolid
{
public:
	BeamMaterial()	: LinearElasticSolid()
	{
		rho0_ = rho0_s;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};

/**
 * application dependent initial condition 
 */
class BeamInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	BeamInitialCondition(SolidBody *beam)
		: solid_dynamics::ElasticDynamicsInitialCondition(beam) {};
protected:
	void Update(size_t index_i, Real dt) override {
		/** initial velocity profile */
		Real x = pos_n_[index_i][0] / PL;
		if (x > 0.0) {
			vel_n_[index_i][1] 
				= vf * material_->ReferenceSoundSpeed()*(M*(cos(kl*x) - cosh(kl*x)) - N * (sin(kl*x) - sinh(kl*x))) / Q;
		}
	};
};
/**
* @brief define the beam base which will be constrained.
* NOTE: this class can only be instanced after body particles
* have been generated
*/
class BeamBase : public BodyPartByParticle
{
public:
	BeamBase(SolidBody *solid_body, std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		/* Geometry definition */
		std::vector<Vecd> beam_base_shape = CreatBeamBaseShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		std::vector<Vecd> beam_shape = CreatBeamShape();
		body_part_shape_->addAPolygon(beam_shape, ShapeBooleanOps::sub);

		//tag the particles within the body part
		tagBodyPart();
	}
};

//define an observer body
class BeamObserver : public FictitiousBody
{
public:
	BeamObserver(SPHSystem &system, std::string body_name)
		: FictitiousBody(system, body_name, new ParticleAdaptation(1.15, 2.0))
	{
		body_input_points_volumes_.push_back(std::make_pair(Vecd(PL, 0.0), 0.0));
	}
};
//------------------------------------------------------------------------------
//the main program
//------------------------------------------------------------------------------

int main()
{

	//build up context -- a SPHSystem
	SPHSystem system(system_domain_bounds, resolution_ref);

	//the oscillating beam
	Beam *beam_body = new Beam(system, "BeamBody");
	//Configuration of solid materials
	BeamMaterial *beam_material = new BeamMaterial();
	//creat particles for the elastic body
	ElasticSolidParticles beam_particles(beam_body, beam_material);

	BeamObserver *beam_observer = new BeamObserver(system, "BeamObserver");
	//create observer particles
	BaseParticles observer_particles(beam_observer);

	/** topology */
	BodyRelationInner* beam_body_inner = new BodyRelationInner(beam_body);
	BodyRelationContact* beam_observer_contact = new BodyRelationContact(beam_observer, { beam_body });

	//-----------------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-----------------------------------------------------------------------------
	/** initial condition */
	BeamInitialCondition beam_initial_velocity(beam_body);
	//corrected strong configuration	
	solid_dynamics::CorrectConfiguration
		beam_corrected_configuration_in_strong_form(beam_body_inner);

	//time step size calculation
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(beam_body);

	//stress relaxation for the beam
	solid_dynamics::StressRelaxationFirstHalf
		stress_relaxation_first_half(beam_body_inner);
	solid_dynamics::StressRelaxationSecondHalf
		stress_relaxation_second_half(beam_body_inner);

	// clamping a solid body part. This is softer than a driect constraint
	solid_dynamics::ClampConstrainSolidBodyRegion
		clamp_constrain_beam_base(beam_body_inner, new BeamBase(beam_body, "BeamBase"));

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_beam_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<indexVector, Vecd>
		write_beam_tip_displacement("Position", in_output, beam_observer_contact);
	/**
	 * @brief Setup geomtry and initial conditions
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	beam_initial_velocity.exec();
	beam_corrected_configuration_in_strong_form.parallel_exec();

	//-----------------------------------------------------------------------------
	//from here the time stepping begines
	//-----------------------------------------------------------------------------
	//starting time zero
	GlobalStaticVariables::physical_time_ = 0.0;
	write_beam_states.writeToFile(0);
	write_beam_tip_displacement.writeToFile(0);

	int ite = 0;
	Real T0 = 1.0;
	Real End_Time = T0;
	//time step size for ouput file
	Real D_Time = 0.01*T0;
	Real Dt = 0.1*D_Time;			/**< Time period for data observing */
	Real dt = 0.0; 					//default acoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

	//computation loop starts 
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		//integrate time (loop) until the next output time
		while (integration_time < D_Time) {

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) {

				if (ite % 100 == 0) {
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				stress_relaxation_first_half.parallel_exec(dt);
				clamp_constrain_beam_base.parallel_exec();
				stress_relaxation_second_half.parallel_exec(dt);

				ite++;
				dt = computing_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}

		write_beam_tip_displacement.writeToFile(ite);

		tick_count t2 = tick_count::now();
		write_beam_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
