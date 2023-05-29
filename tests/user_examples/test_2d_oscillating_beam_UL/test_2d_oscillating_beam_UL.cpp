/* ---------------------------------------------------------------------------*
*            SPHinXsys: 2D oscillation beam example-update Lagrange           *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for understanding SPH method for    *
* solid simulation based on update Lagrange method                            *
* In this case, the constraint of the beam is implemented with                *
* internal constrained subregion.                                             *
* ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
#include "all_granular.h"
using namespace SPH;
//------------------------------------------------------------------------------
//global parameters for the case
//------------------------------------------------------------------------------
Real sc = 1;
Real PL = 0.2 * sc;	//beam length
Real PH = 0.02 * sc; //for thick plate; =0.01 for thin plate
Real SL = 0.06 * sc; //depth of the insert
Real resolution_ref = PH / 10;
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
Real c0 = sqrt(Youngs_modulus / (3*(1-2*poisson)*rho0_s));

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
	Vecd(0.0, -PH / 2 - BW), Vecd(-SL - BW, -PH / 2 - BW) };
// a beam shape
std::vector<Vecd> beam_shape{
	Vecd(-SL, -PH / 2), Vecd(-SL, PH / 2), Vecd(PL, PH / 2), Vecd(PL, -PH / 2), Vecd(-SL, -PH / 2) };
//Beam observer location
StdVec<Vecd> observation_location = { Vecd(PL, 0.0) };
//----------------------------------------------------------------------
//	Define the beam body
//----------------------------------------------------------------------
class Beam : public MultiPolygonShape
{
public:
	explicit Beam(const std::string &shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(beam_base_shape, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(beam_shape, ShapeBooleanOps::add);
	}
};
//----------------------------------------------------------------------
//	application dependent initial condition 
//----------------------------------------------------------------------
class BeamInitialCondition
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit BeamInitialCondition(GranularBody &granular_column)
		: fluid_dynamics::FluidInitialCondition(granular_column) {};

protected:
	void update(size_t index_i, Real dt)
	{
		/** initial velocity profile */
		Real x = pos_[index_i][0] / PL;
		if (x > 0.0)
		{
			vel_[index_i][1] = vf * c0 *
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
//------------------------------------------------------------------------------
//the main program
//------------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
	// handle command line arguments
	system.handleCommandlineOptions(ac, av);
#endif //----------------------------------------------------------------------
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	GranularBody beam_body(system, makeShared<Beam>("BeamBody"));
	beam_body.sph_adaptation_->resetKernel<KernelCubicBSpline>();
	beam_body.defineParticlesAndMaterial<GranularMaterialParticles, GranularMaterial>(rho0_s, c0, Youngs_modulus, poisson);
	beam_body.generateParticles<ParticleGeneratorLattice>();
	beam_body.addBodyStateForRecording<Real>("VonMisesStress");

	ObserverBody beam_observer(system, "BeamObserver");
	beam_observer.sph_adaptation_->resetAdaptationRatios(1.15, 2.0);
	beam_observer.generateParticles<ObserverParticleGenerator>(observation_location);

	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	ContactRelation beam_observer_contact(beam_observer, { &beam_body });
	InnerRelation beam_body_inner(beam_body);
	//-----------------------------------------------------------------------------
	//this section define all numerical methods will be used in this case
	//-----------------------------------------------------------------------------

	/** initial condition */
	SimpleDynamics<BeamInitialCondition> beam_initial_velocity(beam_body);
	InteractionWithUpdate<fluid_dynamics::FreeSurfaceIndicationInner> surface_detection(beam_body_inner);
	InteractionDynamics<fluid_dynamics::TransportVelocityCorrectionInner> beam_transport_correction(beam_body_inner);
	InteractionWithUpdate<fluid_dynamics::DensitySummationFreeSurfaceInner> beam_density_by_summation(beam_body_inner);
	Dynamics1Level<fluid_dynamics::Integration1stHalf> granular_pressure_relaxation(beam_body_inner);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfDissipativeRiemann> granular_density_relaxation(beam_body_inner);

	Dynamics1Level<granular_dynamics::ShearStressRelaxation1stHalf> granular_shear_stress_relaxation(beam_body_inner);
	InteractionWithUpdate<granular_dynamics::ShearStressRelaxation2ndHalf> granular_shear_stress_rate_relaxation(beam_body_inner);
	InteractionDynamics<granular_dynamics::ArtificialNormalShearStressRelaxation> granular_artificial_normal_shear_stress_relaxation(beam_body_inner);
	ReduceDynamics<granular_dynamics::GranularAcousticTimeStepSize> computing_time_step_size(beam_body);
	// clamping a solid body part.
	BodyRegionByParticle beam_base(beam_body, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
	SimpleDynamics<granular_dynamics::FixConstraint<BodyPartByParticle>> constraint_beam_base(beam_base);
	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_beam_states(io_environment, system.real_bodies_);
	RegressionTestEnsembleAveraged<ObservedQuantityRecording<Vecd>>
		write_beam_tip_displacement("Position", io_environment, beam_observer_contact);
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	beam_initial_velocity.exec();
	//----------------------------------------------------------------------
	//	Setup computing time-step controls.
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 1.0;
	Real End_Time = T0;
	Real D_Time = End_Time/100;  /**< Time period for data observing */
	Real Dt = 0.1 * D_Time; 
	Real dt = 0.0;			//default acoustic time step sizes
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//-----------------------------------------------------------------------------
	//from here the time stepping begines
	//-----------------------------------------------------------------------------
	write_beam_states.writeToFile(0);
	write_beam_tip_displacement.writeToFile(0);
	//computation loop starts
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < D_Time)
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt)
			{	
				granular_pressure_relaxation.exec(dt); 
				granular_artificial_normal_shear_stress_relaxation.exec();
				granular_shear_stress_relaxation.exec(dt);
				constraint_beam_base.exec();
				granular_shear_stress_rate_relaxation.exec(dt);
				granular_density_relaxation.exec(dt);

				ite++;
				dt = computing_time_step_size.exec();
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
				if (ite % 100 == 0)
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: "
						<< dt << "\n";
				}
				beam_body.updateCellLinkedList();
				beam_body_inner.updateConfiguration();
			}
		}
		write_beam_tip_displacement.writeToFile(ite);
		TickCount t2 = TickCount::now();
		write_beam_states.writeToFile();
		TickCount t3 = TickCount::now();
		interval += t3 - t2;
	}
	TickCount t4 = TickCount::now();
	TimeInterval tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	//system.generate_regression_data_ = true;
	if (system.generate_regression_data_)
	{
		write_beam_tip_displacement.generateDataBase(Vec2d(3.0e-2, 3.0e-2), Vec2d(1.0e-2, 1.0e-2));
	}
	else
	{
		write_beam_tip_displacement.testResult();
	}
	return 0;
}