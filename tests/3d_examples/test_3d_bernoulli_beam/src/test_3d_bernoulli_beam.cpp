/**
 * @file test_3d_bernoulli_beam.cpp
 * @brief Bernoulli beam for verification/validation. The test is based on the analytical solution.
 * @author Bence Rochlitz, Virtonomy GmbH
 */
#include <gtest/gtest.h>
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;

/** Geometry parameters. */
Real PL = 0.2;
Real PH = 0.01;
Real PW = 0.01;
Real SL = 0.01;
Real resolution_ref = PH / 6.0; /**< Initial particle spacing. */
Real BW = resolution_ref * 4;	/**< Boundary width. */
Vecd halfsize_beam(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_beam(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
Vecd halfsize_holder(0.5 * SL, 0.5 * PH, 0.5 * PW);
Vecd translation_holder(-0.5 * SL, 0.5 * PH, 0.5 * PW);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-SL, 0, 0), Vecd(PL, PH, PH));

/**< SimTK geometric modeling resolution. */
int resolution(20);
/** For material properties of the solid. */
Real rho0_s = 1000.0;
Real poisson = 0.3;
Real Youngs_modulus = 5e8;
Real physical_viscosity = Youngs_modulus / 100.0; //physical damping
Real gravity_g = 100.0;							  /**< Value of gravity. */
Real time_to_full_gravity = 0.0;

/** Define the beam body. */
class Beam : public SolidBody
{
public:
	Beam(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{

		body_shape_.add<TriangleMeshShapeBrick>(halfsize_beam, resolution, translation_beam);
		body_shape_.add<TriangleMeshShapeBrick>(halfsize_holder, resolution, translation_holder);
	}
};
/**
 * define time dependent gravity
 */
class TimeDependentGravity : public Gravity
{
public:
	explicit TimeDependentGravity(Vecd gravity_vector) : Gravity(gravity_vector) {}
	virtual Vecd InducedAcceleration(Vecd &position) override
	{
		Real current_time = GlobalStaticVariables::physical_time_;
		return current_time < time_to_full_gravity ? current_time * global_acceleration_ / time_to_full_gravity : global_acceleration_;
	}
};
//define an observer particle generator
class ObserverParticleGenerator : public ParticleGeneratorDirect
{
public:
	ObserverParticleGenerator() : ParticleGeneratorDirect()
	{
		positions_volumes_.push_back(std::make_pair(Vecd(PL, PH, PW), 0.0));
	}
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(system_domain_bounds, resolution_ref);

	/** Creat a Beam body, corresponding material, particles and reaction model. */
	Beam beam_body(system, "BeamBody");
	ElasticSolidParticles beam_particles(beam_body, makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson));
	/** Define Observer. */
	ObserverBody beam_observer(system, "BeamObserver");
	ObserverParticles observer_particles(beam_observer, makeShared<ObserverParticleGenerator>());

	/** topology */
	BodyRelationInner beam_body_inner(beam_body);
	BodyRelationContact beam_observer_contact(beam_observer, {&beam_body});

	//-------- common particle dynamics ----------------------------------------
	TimeDependentGravity gravity(Vec3d(0.0, -gravity_g, 0.0));
	TimeStepInitialization initialize_gravity(beam_body, gravity);

	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** Corrected configuration. */
	solid_dynamics::CorrectConfiguration
		corrected_configuration(beam_body_inner);
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize
		computing_time_step_size(beam_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf
		stress_relaxation_first_half(beam_body_inner);
	/** Setup the damping stress, if you know what you are doing. */
	//stress_relaxation_first_step.setupDampingStressFactor(1.0);
	solid_dynamics::StressRelaxationSecondHalf
		stress_relaxation_second_half(beam_body_inner);
	/** Constrain the holder. */
	TriangleMeshShapeBrick holder_shape(halfsize_holder, resolution, translation_holder);
	BodyRegionByParticle holder(beam_body, "Holder", holder_shape);
	solid_dynamics::ConstrainSolidBodyRegion constrain_holder(beam_body, holder);
	DampingWithRandomChoice<DampingPairwiseInner<indexVector, Vec3d>>
		muscle_damping(beam_body_inner, 0.1, "Velocity", physical_viscosity);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	ObservedQuantityRecording<indexVector, Vecd>
		write_displacement("Position", in_output, beam_observer_contact);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.parallel_exec();
	write_states.writeToFile(0);
	write_displacement.writeToFile(0);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 0.15;
	Real output_period = end_time / 100.0;
	Real dt = 0.0;
	/** Statistics for computing time. */
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	/**
	 * Main loop
	 */
	while (GlobalStaticVariables::physical_time_ < end_time)
	{
		Real integration_time = 0.0;
		while (integration_time < output_period)
		{
			if (ite % 100 == 0)
			{
				std::cout << "N=" << ite << " Time: "
						  << GlobalStaticVariables::physical_time_ << "	dt: "
						  << dt << "\n";
			}

			initialize_gravity.parallel_exec(); // gravity force
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			muscle_damping.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}
		write_displacement.writeToFile(ite);
		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	//From here using google test to check the numerical solution
	const StdLargeVec<Vecd> &pos_0 = beam_particles.pos_0_;
	const StdLargeVec<Vecd> &pos_n = beam_particles.pos_n_;
	Real displacement_max = 0.0;
	for (size_t i = 0; i < pos_0.size(); i++)
	{
		displacement_max = SMAX((pos_n[i] - pos_0[i]).norm(), displacement_max);
	}
	Real displacement_max_analytical = 4.8e-3;														// in mm, absolute max displacement
	EXPECT_NEAR(displacement_max, displacement_max_analytical, displacement_max_analytical * 0.25); // 25% tolerance with 6 particles/thickness
	// the tolerance goes down if the number of particles / thickness is higher
	// for CI testing use 6 particles to be faster
	std::cout << "displacement_max: " << displacement_max << std::endl;

	return 0;
}
