/* ---------------------------------------------------------------------------*
*            SPHinXsys: 2D oscillation beam example-update Lagrange           *
* ----------------------------------------------------------------------------*
* This is the one of the basic test cases for understanding SPH method for    *
* solid simulation based on update Lagrange method                            *
* In this case, the constraint of the beam is implemented with                *
* internal constrained subregion.                                             *
* ----------------------------------------------------------------------------*/
#include "sphinxsys.h"
#include "all_continuum.h"
using namespace SPH;
//------------------------------------------------------------------------------
//global parameters for the case
//------------------------------------------------------------------------------
Real PL = 0.2;	//beam length
Real PH = 0.02; //for thick plate; =0.01 for thin plate
Real SL = 0.06; //depth of the insert
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
//Real poisson = 0.4;		 //Poisson ratio
Real c0 = sqrt(Youngs_modulus / (3 * (1 - 2 * poisson) * rho0_s));
Real gravity_g = 0.0;
//----------------------------------------------------------------------
//	Parameters for initial condition on velocity 
//----------------------------------------------------------------------
Real kl = 1.875;
Real M = sin(kl) + sinh(kl);
Real N = cos(kl) + cosh(kl);
Real Q = 2.0 * (cos(kl) * sinh(kl) - sin(kl) * cosh(kl));
Real vf = 0.05;
Real R = PL / (0.5 * Pi);
//for dual time-step
Real U_max = vf * c0 * (M * (cos(kl) - cosh(kl)) - N * (sin(kl) - sinh(kl))) / Q;
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
	explicit Beam(const std::string& shape_name) : MultiPolygonShape(shape_name)
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
	explicit BeamInitialCondition(RealBody& beam_column)
		: fluid_dynamics::FluidInitialCondition(beam_column) {};
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
int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, resolution_ref);
#ifdef BOOST_AVAILABLE
	// handle command line arguments
	system.handleCommandlineOptions(ac, av);
#endif 
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	RealBody beam_body(system, makeShared<Beam>("BeamBody"));
	beam_body.sph_adaptation_->resetKernel<KernelCubicBSpline>();
	beam_body.defineParticlesAndMaterial<ContinuumParticles, GeneralContinuum>(rho0_s, c0, Youngs_modulus, poisson);
	beam_body.generateParticles<ParticleGeneratorLattice>();
	beam_body.addBodyStateForRecording<Real>("VonMisesStress");
	beam_body.addBodyStateForRecording<Real>("VonMisesStrain");
	beam_body.addBodyStateForRecording<Real>("Pressure");
	beam_body.addBodyStateForRecording<Real>("Density");

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
	SharedPtr<Gravity> gravity_ptr = makeShared<Gravity>(Vecd(0.0, -gravity_g));
	Dynamics1Level<continuum_dynamics::Integration1stHalf> beam_pressure_relaxation(beam_body_inner);
	Dynamics1Level<fluid_dynamics::Integration2ndHalfDissipativeRiemann> beam_density_relaxation(beam_body_inner);
	InteractionDynamics<continuum_dynamics::AngularConservativeShearAccelerationRelaxation>
		beam_shear_acceleration_angular_conservative(beam_body_inner);
	InteractionWithUpdate<CorrectedConfigurationInner> correcttion_matrix(beam_body_inner);
	Dynamics1Level<continuum_dynamics::ShearStressRelaxation> beam_shear_stress_relaxation(beam_body_inner);
	//for dula timestep
	ReduceDynamics<fluid_dynamics::AdvectionTimeStepSize> fluid_advection_time_step(beam_body, U_max, 0.2);
	ReduceDynamics<fluid_dynamics::AcousticTimeStepSize> fluid_acoustic_time_step(beam_body, 0.4);
	// clamping a solid body part.
	BodyRegionByParticle beam_base(beam_body, makeShared<MultiPolygonShape>(createBeamConstrainShape()));
	SimpleDynamics<continuum_dynamics::FixConstraint<BodyPartByParticle>> constraint_beam_base(beam_base);
	//-----------------------------------------------------------------------------
	//outputs
	//-----------------------------------------------------------------------------
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_beam_states(io_environment, system.real_bodies_);
	RegressionTestEnsembleAverage<ObservedQuantityRecording<Vecd>>
		write_beam_tip_displacement("Position", io_environment, beam_observer_contact);
	RegressionTestDynamicTimeWarping<ReducedQuantityRecording<ReduceDynamics<TotalMechanicalEnergy>>>
		write_water_mechanical_energy(io_environment, beam_body, gravity_ptr);
	//----------------------------------------------------------------------
	//	Setup computing and initial conditions.
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	beam_initial_velocity.exec();
	correcttion_matrix.exec();
	//----------------------------------------------------------------------
	//	Setup computing time-step controls.
	//----------------------------------------------------------------------
	int ite = 0;
	Real T0 = 1.0;
	Real End_Time = T0;
	Real D_Time = End_Time / 100;  /**< Time period for data observing */
	TickCount t1 = TickCount::now();
	TimeInterval interval;
	//-----------------------------------------------------------------------------
	//from here the time stepping begines
	//-----------------------------------------------------------------------------
	write_beam_states.writeToFile(0);
	write_beam_tip_displacement.writeToFile(0);
	write_water_mechanical_energy.writeToFile(0);
	//computation loop starts
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < D_Time)
		{
			Real relaxation_time = 0.0;
			Real advection_dt = fluid_advection_time_step.exec();

			while (relaxation_time < advection_dt)
			{
				Real acoustic_dt = fluid_acoustic_time_step.exec();
				beam_shear_stress_relaxation.exec(acoustic_dt);
				beam_pressure_relaxation.exec(acoustic_dt);
				constraint_beam_base.exec();
				beam_density_relaxation.exec(acoustic_dt);
				//shear acceleration with angular conservative
				beam_shear_acceleration_angular_conservative.exec(acoustic_dt);
				ite++;
				relaxation_time += acoustic_dt;
				integration_time += acoustic_dt;
				GlobalStaticVariables::physical_time_ += acoustic_dt;
				if (ite % 500 == 0)
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	advection_dt: "
						<< advection_dt << "	acoustic_dt: "
						<< acoustic_dt << "\n";
				}
			}
			beam_body.updateCellLinkedList();
			beam_body_inner.updateConfiguration();
			correcttion_matrix.exec();
		}
		write_beam_tip_displacement.writeToFile(ite);
		write_water_mechanical_energy.writeToFile(ite);
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
        write_beam_tip_displacement.generateDataBase(Vec2d(1.0e-2, 1.0e-2), Vec2d(1.0e-2, 1.0e-2));
    }
    else
    {
        write_beam_tip_displacement.testResult();
    }
	return 0;
}