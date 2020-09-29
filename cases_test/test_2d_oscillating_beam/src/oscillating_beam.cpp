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
Real particle_spacing_ref = PH / 10.0;
Real BW = particle_spacing_ref * 4; 	//boundary width

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
	Beam(SPHSystem &system, string body_name, int refinement_level)
		: SolidBody(system, body_name, refinement_level)
	{
		/** Geometry definition. */
		std::vector<Point> beam_base_shape = CreatBeamBaseShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		std::vector<Point> beam_shape = CreatBeamShape();
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
		rho_0_ = rho0_s;
		E_0_ = Youngs_modulus;
		nu_ = poisson;

		assignDerivedMaterialParameters();
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
	BeamBase(SolidBody *solid_body, string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		/* Geometry definition */
		std::vector<Point> beam_base_shape = CreatBeamBaseShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		std::vector<Point> beam_shape = CreatBeamShape();
		body_part_shape_->addAPolygon(beam_shape, ShapeBooleanOps::sub);

		//tag the particles within the body part
		tagBodyPart();
	}
};

//define an observer body
class BeamObserver : public FictitiousBody
{
public:
	BeamObserver(SPHSystem &system, string body_name, int refinement_level)
		: FictitiousBody(system, body_name, refinement_level, 1.3)
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

	//the oscillating beam
	Beam *beam_body = new Beam(system, "BeamBody", 0);
	//Configuration of solid materials
	BeamMaterial *beam_material = new BeamMaterial();
	//creat particles for the elastic body
	ElasticSolidParticles beam_particles(beam_body, beam_material);

	BeamObserver *beam_observer = new BeamObserver(system, "BeamObserver", 1);
	//create observer particles
	BaseParticles observer_particles(beam_observer);

	/** topology */
	SPHBodyInnerRelation* beam_body_inner = new SPHBodyInnerRelation(beam_body);
	SPHBodyContactRelation* beam_observer_contact = new SPHBodyContactRelation(beam_observer, { beam_body });

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
	WriteAnObservedQuantity<Vecd, BaseParticles, &BaseParticles::pos_n_>
		write_beam_tip_displacement("Displacement", in_output, beam_observer_contact);
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
	write_beam_states.WriteToFile(GlobalStaticVariables::physical_time_);
	write_beam_tip_displacement.WriteToFile(GlobalStaticVariables::physical_time_);

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
					cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}

				stress_relaxation_first_half.parallel_exec(dt);
				constrain_beam_base.parallel_exec(dt);
				stress_relaxation_second_half.parallel_exec(dt);

				ite++;
				dt = computing_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
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
