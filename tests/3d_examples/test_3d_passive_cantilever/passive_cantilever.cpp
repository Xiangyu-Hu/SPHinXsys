/**
 * @file passive_cantilever.cpp
 * @brief This is the first example of cantilever 
 * @author Chi Zhang and Xiangyu Hu
 * @ref 	doi.org/10.1016/j.jcp.2013.12.012
 */
#include "sphinxsys.h"
/** Name space. */
using namespace SPH;
/** Geometry parameters. */
Real PL = 6.0;
Real PH = 1.0;
Real PW = 1.0;
Real SL = 0.5;
Real resolution_ref = PH / 12.0; /**< Initial particle spacing. */
Real BW = resolution_ref * 4;	 /**< Boundary width. */
Vecd halfsize_cantilever(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_cantilever(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
Vecd halfsize_holder(0.5 * SL, 0.5 * PH, 0.5 * PW);
Vecd translation_holder(-0.5 * SL, 0.5 * PH, 0.5 * PW);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vecd(-SL - BW, -BW, -BW),
								 Vecd(PL + BW, PH + BW, PH + BW));

/**< SimTK geometric modeling resolution. */
int resolution(20);
/** For material properties of the solid. */
Real rho0_s = 1100.0;
Real poisson = 0.45;
Real Youngs_modulus = 1.7e7;
Real a = Youngs_modulus / (2.0 * (1.0 + poisson));
Real a_f = 0.0 * a;
Real a0[4] = {a, a_f, 0.0, 0.0};
Real b0[4] = {1.0, 0.0, 0.0, 0.0};
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
Real bulk_modulus = Youngs_modulus / 3.0 / (1.0 - 2.0 * poisson);

/** Define the cantilever body. */
class Cantilever : public SolidBody
{
public:
	Cantilever(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		body_shape_.add<TriangleMeshShapeBrick>(halfsize_cantilever, resolution, translation_cantilever);
		body_shape_.add<TriangleMeshShapeBrick>(halfsize_holder, resolution, translation_holder);
	}
};
/**
 * application dependent initial condition 
 */
class CantileverInitialCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	explicit CantileverInitialCondition(SolidBody &cantilever)
		: solid_dynamics::ElasticDynamicsInitialCondition(cantilever){};

protected:
	void Update(size_t index_i, Real dt) override
	{
		if (pos_n_[index_i][0] > 0.0)
		{
			vel_n_[index_i][1] = 5.0 * sqrt(3.0);
			vel_n_[index_i][2] = 5.0;
		}
	};
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
	/** Setup the system. Please the make sure the global domain bounds are correctly defined. */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** create a Cantilever body, corresponding material, particles and reaction model. */
	Cantilever cantilever_body(system, "CantileverBody");
	SharedPtr<Muscle> cantilever_material = makeShared<Muscle>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
	ElasticSolidParticles particles(cantilever_body, cantilever_material);
	/** Define Observer. */
	ObserverBody cantilever_observer(system, "CantileverObserver");
	ObserverParticles observer_particles(cantilever_observer, makeShared<ObserverParticleGenerator>());

	/** topology */
	BodyRelationInner cantilever_body_inner(cantilever_body);
	BodyRelationContact cantilever_observer_contact(cantilever_observer, {&cantilever_body});

	/** 
	 * This section define all numerical methods will be used in this case.
	 */
	/** Initialization. */
	CantileverInitialCondition initialization(cantilever_body);
	/** Corrected configuration. */
	solid_dynamics::CorrectConfiguration
		corrected_configuration(cantilever_body_inner);
	/** Time step size calculation. */
	solid_dynamics::AcousticTimeStepSize
		computing_time_step_size(cantilever_body);
	/** active and passive stress relaxation. */
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(cantilever_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(cantilever_body_inner);
	/** Constrain the holder. */
	TriangleMeshShapeBrick holder_shape(halfsize_holder, resolution, translation_holder);
	BodyRegionByParticle holder(cantilever_body, "Holder", holder_shape);
	solid_dynamics::ConstrainSolidBodyRegion constrain_holder(cantilever_body, holder);
	/** Output */
	In_Output in_output(system);
	BodyStatesRecordingToVtp write_states(in_output, system.real_bodies_);
	RegressionTestDynamicTimeWarping<ObservedQuantityRecording<Vecd>>
		write_displacement("Position", in_output, cantilever_observer_contact);
	/**
	 * From here the time stepping begins.
	 * Set the starting time.
	 */
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	/** apply initial condition */
	initialization.parallel_exec();
	corrected_configuration.parallel_exec();
	write_states.writeToFile(0);
	write_displacement.writeToFile(0);
	/** Setup physical parameters. */
	int ite = 0;
	Real end_time = 3.0;
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
			stress_relaxation_first_half.parallel_exec(dt);
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

	write_displacement.newResultTest();

	return 0;
}
