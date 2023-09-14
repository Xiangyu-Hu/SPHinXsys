/**
 * @file 	relaxation_evolution.cpp
 * @brief 	This is the first case by testing the relaxation with evolution method.
 * @author 	Bo Zhang
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real LL = 1.0;
Real LH = 1.0;
Real resolution_ref = LH / 50.0;
BoundingBox system_domain_bounds(Vec2d::Zero(), Vec2d(LL, LH));
//----------------------------------------------------------------------
//	Define geometries
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH);
Vec2d water_block_translation = water_block_halfsize;

class Insert : public ComplexShape
{
public:
	explicit Insert(const std::string& shape_name) : ComplexShape(shape_name)
	{
		add<TransformShape<GeometricShapeBox>>(Transform(water_block_halfsize), water_block_translation);
	}
};

class TestingInitialCondition
	: public fluid_dynamics::FluidInitialCondition
{
public:
	explicit TestingInitialCondition(SPHBody& sph_body)
		: FluidInitialCondition(sph_body), pos_(particles_->pos_),
		  p_(*particles_->getVariableByName<Real>("Pressure")) {};

	void update(size_t index_i, Real dt)
	{
		/* initial pressure distribution. */
		p_[index_i] = pos_[index_i][0];
	}

protected:
	StdLargeVec<Vecd>& pos_;
	StdLargeVec<Real>& p_;
};

int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	sph_system.setRunParticleRelaxation(true);
	sph_system.setReloadParticles(true);
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
#endif
	IOEnvironment io_environment(sph_system);

	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody body(sph_system, makeShared<Insert>("WaterBody"));
	body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(1, 1, 1);
	body.addBodyStateForRecording<Vecd>("Position");
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? body.generateParticles<ParticleGeneratorReload>(io_environment, body.getName())
		: body.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Define body relation map.
		//----------------------------------------------------------------------
		InnerRelation insert_body_inner(body);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_insert_body_to_vtp(io_environment, { &body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, { &body });

		/* Relaxation method: including based on the 0th and 1st order consistency. */
		relax_dynamics::RelaxationStepInner relaxation_0th_inner(insert_body_inner, true);
		relax_dynamics::RelaxationStepImplicitInner relaxation_0th_implicit_inner(insert_body_inner, true);
		InteractionDynamics<relax_dynamics::CalculateCorrectionMatrix> calculate_correction_matrix(insert_body_inner, true);
		relax_dynamics::RelaxationStepByCMInner relaxation_1st_inner(insert_body_inner, true);
		relax_dynamics::RelaxationStepByCMImplicitInner relaxation_1st_implicit_inner(insert_body_inner, true);

		/* Update the relaxation residual. */
		InteractionDynamics<relax_dynamics::CheckConsistencyRealization> check_consistency_realization(insert_body_inner, true);
		ReduceAverage<QuantitySummation<Real>> calculate_pressure_gradient(body, "PressureGradientErrorNorm");
		ReduceAverage<QuantitySummation<Real>> calculate_zero_order_error(body, "ZeroOrderErrorNorm");
		ReduceAverage<QuantitySummation<Real>> calculate_reproduce_gradient(body, "ReproduceGradientErrorNorm");
		InteractionDynamics<relax_dynamics::CheckReverseConsistencyRealization> check_reverse_consistency_realization(insert_body_inner, true);
		ReduceAverage<QuantitySummation<Real>> calculate_pressure_gradient_reverse(body, "PressureGradientErrorNormReverse");
		ReduceAverage<QuantitySummation<Real>> calculate_zero_order_error_reverse(body, "ZeroOrderErrorNormReverse");
		ReduceAverage<QuantitySummation<Real>> calculate_reproduce_gradient_reverse(body, "ReproduceGradientErrorNormReverse");
		SimpleDynamics<TestingInitialCondition> testing_initial_condition(body);
		//----------------------------------------------------------------------  
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_insert_body_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		sph_system.initializeSystemConfigurations();
		testing_initial_condition.exec();
		write_insert_body_to_vtp.writeToFile(0);

		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();
		int ite = 0; //iteration step for the total relaxation step.

		std::string filefullpath_all_information = io_environment.output_folder_ + "/" + "all_information.dat";
		std::ofstream out_file_all_information(filefullpath_all_information.c_str(), std::ios::app);

		GlobalStaticVariables::physical_time_ = ite;
		/* The procedure to obtain uniform particle distribution that satisfies the 0th order consistency. */
		while (ite < 5)
		{
			calculate_correction_matrix.exec();
			relaxation_1st_inner.exec();
			body.updateCellLinkedList();
			insert_body_inner.updateConfiguration();

			ite++;

			if (ite % 1 == 0)
			{
				calculate_correction_matrix.exec();
				check_consistency_realization.exec();
				check_reverse_consistency_realization.exec();
				out_file_all_information << std::fixed << std::setprecision(12) << ite << "   " <<
					calculate_pressure_gradient.exec() << "   " << calculate_pressure_gradient_reverse.exec() << "   " <<
					calculate_zero_order_error.exec() << "   " << calculate_zero_order_error_reverse.exec() << "   " <<
					calculate_reproduce_gradient.exec() << "   " << calculate_reproduce_gradient_reverse.exec() << "\n";
				std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";

				write_insert_body_to_vtp.writeToFile(ite);
			}
		}

		ite++;
		write_insert_body_to_vtp.writeToFile(ite);
		write_particle_reload_files.writeToFile(0);

		TickCount t2 = TickCount::now();
		TickCount::interval_t tt;
		tt = t2 - t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
		return 0;
	}

	//----------------------------------------------------------------------
	//	Define body relation map.
	//----------------------------------------------------------------------
	InnerRelation insert_body_inner(body);
	InteractionDynamics<relax_dynamics::CalculateCorrectionMatrix> calculate_correction_matrix(insert_body_inner, true);
	InteractionDynamics<relax_dynamics::CheckConsistencyRealization> check_consistency_realization(insert_body_inner, true);
	ReduceAverage<QuantitySummation<Real>> calculate_pressure_gradient(body, "PressureGradientErrorNorm");
	ReduceAverage<QuantitySummation<Real>> calculate_zero_order_error(body, "ZeroOrderErrorNorm");
	ReduceAverage<QuantitySummation<Real>> calculate_reproduce_gradient(body, "ReproduceGradientErrorNorm");
	InteractionDynamics<relax_dynamics::CheckReverseConsistencyRealization> check_reverse_consistency_realization(insert_body_inner, true);
	ReduceAverage<QuantitySummation<Real>> calculate_pressure_gradient_reverse(body, "PressureGradientErrorNormReverse");
	ReduceAverage<QuantitySummation<Real>> calculate_zero_order_error_reverse(body, "ZeroOrderErrorNormReverse");
	ReduceAverage<QuantitySummation<Real>> calculate_reproduce_gradient_reverse(body, "ReproduceGradientErrorNormReverse");
	SimpleDynamics<TestingInitialCondition> testing_initial_condition(body);

	body.addBodyStateForRecording<Real>("Pressure");
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_insert_body_states(io_environment, sph_system.real_bodies_);
	std::string filefullpath_all_information = io_environment.output_folder_ + "/" + "all_information.dat";
	std::ofstream out_file_all_information(filefullpath_all_information.c_str(), std::ios::app);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	/** initialize cell linked lists for all bodies. */
	sph_system.initializeSystemCellLinkedLists();
	/** initialize configurations for all bodies. */
	sph_system.initializeSystemConfigurations();
	testing_initial_condition.exec();

	calculate_correction_matrix.exec();
	check_consistency_realization.exec();
	check_reverse_consistency_realization.exec();
	out_file_all_information << std::fixed << std::setprecision(12) << 1 << "   " <<
		calculate_pressure_gradient.exec() << "   " << calculate_pressure_gradient_reverse.exec() << "   " <<
		calculate_zero_order_error.exec() << "   " << calculate_zero_order_error_reverse.exec() << "   " <<
		calculate_reproduce_gradient.exec() << "   " << calculate_reproduce_gradient_reverse.exec() << "\n";
	write_insert_body_states.writeToFile();
}
