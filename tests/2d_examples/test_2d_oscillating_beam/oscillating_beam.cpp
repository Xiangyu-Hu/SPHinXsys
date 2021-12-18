/* ---------------------------------------------------------------------------*
*            SPHinXsys: 2D oscillation beam example-one body version           *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases, also the first case for            *
* understanding SPH method for solid simulation.                              *
* In this case, the constraint of the beam is implemented with                *
* internal constrained subregion.                                             *
* ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
using namespace SPH;
//------------------------------------------------------------------------------
//global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;	//beam length
Real PH = 0.02; //for thick plate; =0.01 for thin plate
Real SL = 0.06; //depth of the insert
//reference particle spacing
Real resolution_ref = PH / 10.0;
Real BW = resolution_ref * 4; //boundary width, at least three particles
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
								 Vec2d(PL + 3.0 * BW, PL / 2.0));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;		 //reference density
Real Youngs_modulus = 2.0e6; //reference Youngs modulus
Real poisson = 0.3975;		 //Poisson ratio
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * Pi);
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
// a beam base shape
std::vector<Vecd> beam_base_shape{
	Vecd(-SL - BW, -PH / 2 - BW), Vecd(-SL - BW, PH / 2 + BW), Vecd(0.0, PH / 2 + BW),
	Vecd(0.0, -PH / 2 - BW), Vecd(-SL - BW, -PH / 2 - BW)};
// a beam shape
std::vector<Vecd> beam_shape{
	Vecd(-SL, -PH / 2), Vecd(-SL, PH / 2), Vecd(PL, PH / 2), Vecd(PL, -PH / 2), Vecd(-SL, -PH / 2)};
//----------------------------------------------------------------------
//	Define the beam body
//----------------------------------------------------------------------
class Beam : public SolidBody
{
public:
	Beam(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		multi_polygon.addAPolygon(beam_shape, ShapeBooleanOps::add);
		body_shape_.add<MultiPolygonShape>(multi_polygon);
	}
};
//----------------------------------------------------------------------
//	application dependent initial condition 
//----------------------------------------------------------------------
class BeamInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit BeamInitialCondition(SolidBody &beam)
		: solid_dynamics::ElasticDynamicsInitialCondition(beam){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		/** initial velocity profile */
		Real x = pos_n_[index_i][0] / PL;
		if (x > 0.0)
		{
			vel_n_[index_i][1] = vf * material_->ReferenceSoundSpeed() *
								 (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x))) / Q;
		}
	};
};
//----------------------------------------------------------------------
//	define the beam base which will be constrained.
//----------------------------------------------------------------------
MultiPolygon createBeamConstrainShape()
{
	MultiPolygon multi_polygon;
	multi_polygon.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
	multi_polygon.addAPolygon(beam_shape, ShapeBooleanOps::sub);
	return multi_polygon;
};
//----------------------------------------------------------------------
//	Observer particle generator
//----------------------------------------------------------------------
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(Vecd(PL, 0.0), 0.0));
	}
};
//------------------------------------------------------------------------------
//the main program
//------------------------------------------------------------------------------

int main()
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	//the oscillating beam
	Beam beam_body(system, "BeamBody");
	//create particles for the elastic body
	ElasticSolidParticles beam_particles(beam_body, makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson));

	ObserverBody beam_observer(system, "BeamObserver", makeShared<SPHAdaptation>(1.15, 2.0));
	//create observer particles
	ObserverParticles observer_particles(beam_observer, makeShared<ObserverParticleGenerator>());

	/** topology */
	BodyRelationInner beam_body_inner(beam_body);
	BodyRelationContact beam_observer_contact(beam_observer, {&beam_body});

	//-----------------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-----------------------------------------------------------------------------
	/** initial condition */
	BeamInitialCondition beam_initial_velocity(beam_body);
	//corrected strong configuration
	solid_dynamics::CorrectConfiguration
		beam_corrected_configuration(beam_body_inner);

	//time step size calculation
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(beam_body);

	//stress relaxation for the beam
	solid_dynamics::StressRelaxationFirstHalf
		stress_relaxation_first_half(beam_body_inner);
	solid_dynamics::StressRelaxationSecondHalf
		stress_relaxation_second_half(beam_body_inner);

	// clamping a solid body part. This is softer than a driect constraint
	MultiPolygonShape beam_cobstrain_shape(createBeamConstrainShape());
	BodyRegionByParticle beam_base(beam_body, "BeamBase", beam_cobstrain_shape);
	solid_dynamics::ClampConstrainSolidBodyRegion clamp_constrain_beam_base(beam_body_inner, beam_base);

	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToVtp write_beam_states(in_output, system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<indexVector, Vecd>>
		write_beam_tip_displacement("Position", in_output, beam_observer_contact);

	/**
	 * @brief Setup geometry and initial conditions
	 */
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	beam_initial_velocity.exec();
	beam_corrected_configuration.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 1.0;
	Real End_Time = T0;
	//time step size for output file
	Real D_Time = 0.01 * T0;
	Real Dt = 0.1 * D_Time; /**< Time period for data observing */
	Real dt = 0.0;			//default acoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//-----------------------------------------------------------------------------
	//from here the time stepping begines
	//-----------------------------------------------------------------------------
	write_beam_states.writeToFile(0);
	write_beam_tip_displacement.writeToFile(0);

	//computation loop starts
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		//integrate time (loop) until the next output time
		while (integration_time < D_Time)
		{

			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{
				stress_relaxation_first_half.parallel_exec(dt);
				clamp_constrain_beam_base.parallel_exec();
				stress_relaxation_second_half.parallel_exec(dt);

				ite++;
				dt = computing_time_step_size.parallel_exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;

				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: "
							  << dt << "\n";
				}
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

	write_beam_tip_displacement.newResultTest();

	return 0;
}
