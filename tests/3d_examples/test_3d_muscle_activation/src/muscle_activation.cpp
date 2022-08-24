/**
 * @file muscle_activation.cpp
 * @brief This is the first example of electro activation of myocardium
 * @author Chi Zhang and Xiangyu Hu
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real PL = 1.0;					 /**< Length of the myocardium body. */
Real PH = 1.0;					 /**< Thickness of the myocardium body. */
Real PW = 1.0;					 /**< Width of the myocardium body. */
Real resolution_ref = PH / 25.0; /**< Initial particle spacing. */
Real SL = 4.0 * resolution_ref;	 /**< Extension for holder. */
Vecd halfsize_myocardium(0.5 * (PL + SL), 0.5 * PH, 0.5 * PW);
Vecd translation_myocardium(0.5 * (PL - SL), 0.5 * PH, 0.5 * PW);
Vecd halfsize_holder(0.5 * SL, 0.5 * PH, 0.5 * PW);
Vecd translation_holder(-0.5 * SL, 0.5 * PH, 0.5 * PW);
//----------------------------------------------------------------------
//	Material parameters.
//----------------------------------------------------------------------
Real rho0_s = 1.0;
Real a0[4] = {0.059, 0.0, 0.0, 0.0};
Real b0[4] = {8.023, 0.0, 0.0, 0.0};
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
Real reference_voltage = 30.0;
Real linear_active_stress_factor = -0.5;
/** reference stress or bulk modulus to achieve weakly compressible condition */
Real bulk_modulus = 30.0 * reference_voltage * fabs(linear_active_stress_factor);
//----------------------------------------------------------------------
//	Case dependent muscle activation history.
//----------------------------------------------------------------------
class MyocardiumActivation
	: public active_muscle_dynamics::MuscleActivation
{
public:
	explicit MyocardiumActivation(SPHBody &sph_body)
		: active_muscle_dynamics::MuscleActivation(sph_body){};

	void update(size_t index_i, Real dt)
	{
		Real voltage = pos0_[index_i][0] <= 0 ? 0.0 : reference_voltage * pos0_[index_i][0] / PL;
		active_contraction_stress_[index_i] += GlobalStaticVariables::physical_time_ <= 1.0
												   ? linear_active_stress_factor * voltage * dt
												   : 0.0;
	};
};
//----------------------------------------------------------------------
//	Case dependent constraint.
//----------------------------------------------------------------------
class ConstrainHolder : public solid_dynamics::ConstrainSolidBodyRegion
{
public:
	ConstrainHolder(SolidBody &body, BodyPartByParticle &body_part, int axis_id)
		: solid_dynamics::ConstrainSolidBodyRegion(body, body_part),
		  axis_id_(axis_id){};

protected:
	int axis_id_;
	virtual Vecd getDisplacement(Vecd &pos_0, Vecd &pos_n)
	{
		Vecd pos_temp = pos_n;
		pos_temp[axis_id_] = pos_0[axis_id_];
		return pos_temp;
	};
	virtual Vecd getVelocity(Vecd &pos_0, Vecd &pos_n, Vecd &vel_n)
	{
		Vecd vel_temp = vel_n;
		vel_temp[axis_id_] = 0.0;
		return vel_temp;
	};
	virtual Vecd getAcceleration(Vecd &pos_0, Vecd &pos_n, Vecd &acc)
	{
		Vecd acc_temp = acc;
		acc_temp[axis_id_] = 0.0;
		return acc_temp;
	};
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main()
{
	//----------------------------------------------------------------------
	//	Build up an SPHSystem.
	//----------------------------------------------------------------------
	BoundingBox system_domain_bounds(Vecd(-SL, -SL, -SL), Vecd(PL + SL, PH + SL, PW + SL));
	SPHSystem system(system_domain_bounds, resolution_ref);
	//----------------------------------------------------------------------
	//	Creating bodies with corresponding materials and particles.
	//----------------------------------------------------------------------
	SolidBody myocardium_muscle_body(system, makeShared<TransformShape<GeometricShapeBox>>(
						Transform3d(translation_myocardium), halfsize_myocardium, "MyocardiumMuscleBody"));
	myocardium_muscle_body.defineParticlesAndMaterial<
		ElasticSolidParticles, ActiveMuscle<Muscle>>(rho0_s, bulk_modulus, fiber_direction, sheet_direction, a0, b0);
	myocardium_muscle_body.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner myocardium_muscle_body_inner(myocardium_muscle_body);
	//----------------------------------------------------------------------
	//	Define the numerical methods used in the simulation.
	//	Note that there may be data dependence on the sequence of constructions.
	//----------------------------------------------------------------------
	solid_dynamics::StressRelaxationFirstHalf stress_relaxation_first_half(myocardium_muscle_body_inner);
	solid_dynamics::StressRelaxationSecondHalf stress_relaxation_second_half(myocardium_muscle_body_inner);
	solid_dynamics::CorrectConfiguration corrected_configuration(myocardium_muscle_body_inner);
	solid_dynamics::AcousticTimeStepSize computing_time_step_size(myocardium_muscle_body);
	SimpleDynamics<MyocardiumActivation> myocardium_activation(myocardium_muscle_body);
	BodyRegionByParticle holder(myocardium_muscle_body, makeShared<TransformShape<GeometricShapeBox>>(translation_holder, halfsize_holder));
	ConstrainHolder constrain_holder(myocardium_muscle_body, holder, 0);
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations, observations
	//	and regression tests of the simulation.
	//----------------------------------------------------------------------
	IOEnvironment io_environment(system);
	BodyStatesRecordingToVtp write_states(io_environment, system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	GlobalStaticVariables::physical_time_ = 0.0;
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	corrected_configuration.parallel_exec();
	//----------------------------------------------------------------------
	//	Setup for time-stepping control
	//----------------------------------------------------------------------
	int ite = 0;
	Real end_time = 1.2;
	Real output_period = end_time / 60.0;
	Real dt = 0.0;
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	First output before the main loop.
	//----------------------------------------------------------------------
	write_states.writeToFile(0);
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
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
			myocardium_activation.parallel_exec(dt);
			stress_relaxation_first_half.parallel_exec(dt);
			constrain_holder.parallel_exec(dt);
			stress_relaxation_second_half.parallel_exec(dt);

			ite++;
			dt = computing_time_step_size.parallel_exec();
			integration_time += dt;
			GlobalStaticVariables::physical_time_ += dt;
		}

		tick_count t2 = tick_count::now();
		write_states.writeToFile();
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

	return 0;
}
