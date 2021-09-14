/**
* @file 	self_contact.cpp
* @brief 	This is the case file for the test of dynamic self contact.
* @author   Xiangyu Hu
*/

#include "sphinxsys.h"
using namespace SPH;

//------------------------------------------------------------------------------
//global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;	//beam length
Real PH = 0.01; //for thick plate; 0.01 for thin plate
Real SL = 0.04; //depth of the insert
Real resolution_ref = PH / 10.0;
Real BW = resolution_ref * 4; //boundary width, at least three particles
// Domain bounds of the system.
BoundingBox system_domain_bounds(Vec2d(-SL - BW, -PL / 2.0),
								 Vec2d(PL + 3.0 * BW, PL / 2.0));
//----------------------------------------------------------------------
//	Global parameters for material properties.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;		 //reference density
Real Youngs_modulus = 1.0e5; //reference Youngs modulus
Real poisson = 0.45;		 //Poisson ratio
//----------------------------------------------------------------------
//	Global parameters for initial condition
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.15;
Real R = PL / (0.5 * Pi);
//----------------------------------------------------------------------
//	Geometries used in the case.
//----------------------------------------------------------------------
std::vector<Vecd> creatBeamBaseShape()
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
std::vector<Vecd> creatBeamShape()
{
	std::vector<Vecd> beam_shape;
	beam_shape.push_back(Vecd(-SL, -PH / 2));
	beam_shape.push_back(Vecd(-SL, PH / 2));
	beam_shape.push_back(Vecd(PL, PH / 2));
	beam_shape.push_back(Vecd(PL, -PH / 2));
	beam_shape.push_back(Vecd(-SL, -PH / 2));

	return beam_shape;
}
//----------------------------------------------------------------------
//	Bodies used in the case.
//----------------------------------------------------------------------
class Beam : public SolidBody
{
public:
	Beam(SPHSystem &system, std::string body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		std::vector<Vecd> beam_base_shape = creatBeamBaseShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		std::vector<Vecd> beam_shape = creatBeamShape();
		body_shape_->addAPolygon(beam_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	Body parts usually for impose constraints.
//----------------------------------------------------------------------
class BeamBase : public BodyPartByParticle
{
public:
	BeamBase(SolidBody *solid_body, std::string constrained_region_name)
		: BodyPartByParticle(solid_body, constrained_region_name)
	{
		/* Geometry definition */
		std::vector<Vecd> beam_base_shape = creatBeamBaseShape();
		body_part_shape_ = new ComplexShape(constrained_region_name);
		body_part_shape_->addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		std::vector<Vecd> beam_shape = creatBeamShape();
		body_part_shape_->addAPolygon(beam_shape, ShapeBooleanOps::sub);

		//tag the particles within the body part
		tagBodyPart();
	}
};
// an observer body
class BeamObserver : public FictitiousBody
{
public:
	BeamObserver(SPHSystem &system, std::string body_name)
		: FictitiousBody(system, body_name, new ParticleAdaptation(1.15, 2.0))
	{
		body_input_points_volumes_.push_back(std::make_pair(Vecd(PL, 0.0), 0.0));
	}
};
//----------------------------------------------------------------------
//	Materials used in the case.
//----------------------------------------------------------------------
class BeamMaterial : public LinearElasticSolid
{
public:
	BeamMaterial() : LinearElasticSolid()
	{
		rho0_ = rho0_s;
		youngs_modulus_ = Youngs_modulus;
		poisson_ratio_ = poisson;

		assignDerivedMaterialParameters();
	}
};
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
class BeamInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit BeamInitialCondition(SolidBody *beam)
		: solid_dynamics::ElasticDynamicsInitialCondition(beam){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		/** initial velocity profile */
		Real x = pos_n_[index_i][0] / PL;
		if (x > 0.0)
		{
			vel_n_[index_i][1] = vf * material_->ReferenceSoundSpeed() / Q *
								 (M * (cos(kl * x) - cosh(kl * x)) - N * (sin(kl * x) - sinh(kl * x)));
		}
	};
};
//------------------------------------------------------------------------------
//the main program
//------------------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up -- a SPHSystem
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	Beam *beam_body = new Beam(system, "BeamBody");
	BeamMaterial *beam_material = new BeamMaterial();
	ElasticSolidParticles beam_particles(beam_body, beam_material);
	beam_particles.addAVariableToWrite<indexScalar, Real>("ContactDensity");

	BeamObserver *beam_observer = new BeamObserver(system, "BeamObserver");
	BaseParticles observer_particles(beam_observer);
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner *beam_body_inner = new BodyRelationInner(beam_body);
	BodyRelationContact *beam_observer_contact = new BodyRelationContact(beam_observer, {beam_body});
	SolidBodyRelationSelfContact *beam_self_contact = new SolidBodyRelationSelfContact(beam_body);
	//-----------------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-----------------------------------------------------------------------------
	// initial condition
	BeamInitialCondition beam_initial_velocity(beam_body);
	//corrected strong configuration
	solid_dynamics::CorrectConfiguration
		beam_corrected_configuration_in_strong_form(beam_body_inner);
	//time step size calculation
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(beam_body);
	//stress relaxation for the beam
	solid_dynamics::KirchhoffStressRelaxationFirstHalf stress_relaxation_first_half(beam_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(beam_body_inner);
	// algorithms for solid self contact
	solid_dynamics::DynamicSelfContactForce beam_self_contact_forces(beam_self_contact);
	// clamping a solid body part. This is softer than a driect constraint
	solid_dynamics::ClampConstrainSolidBodyRegion
		clamp_constrain_beam_base(beam_body_inner, new BeamBase(beam_body, "BeamBase"));
	//-----------------------------------------------------------------------------
	//	outputs
	//-----------------------------------------------------------------------------
	In_Output in_output(system);
	BodyStatesRecordingToVtu write_beam_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<indexVector, Vecd>
		write_beam_tip_displacement("Position", in_output, beam_observer_contact);
	//-----------------------------------------------------------------------------
	//	Setup particle configuration and initial conditions
	//-----------------------------------------------------------------------------
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
	Real D_Time = 0.01 * T0;
	Real Dt = 0.1 * D_Time; /**< Time period for data observing */
	Real dt = 0.0;			//default acoustic time step sizes

	//statistics for computing time
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;

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

				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
							  << GlobalStaticVariables::physical_time_ << "	dt: "
							  << dt << "\n";
				}

				beam_self_contact_forces.parallel_exec();
				beam_body->updateCellLinkedList();
				beam_self_contact->updateConfiguration();

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
