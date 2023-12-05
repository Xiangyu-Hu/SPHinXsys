/**
 * @file 	relaxation_evolution_periodic.cpp
 * @brief   This is the first case with periodical boundary condition
 *          by testing the relaxation and consistency.
 * @author  Bo Zhang, Xiangyu Hu.
 */
#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real LL = 1.0;
Real LH = 1.0;
Real resolution_ref = LH / 40.0;
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
		: FluidInitialCondition(sph_body), pos_(particles_->pos_)
	{
         particles_->registerVariable(scalar_, "Scalar");
		 particles_->registerVariable(vector_, "Vector");
         particles_->registerVariable(matrix_, "Matrix");
	};

	void update(size_t index_i, Real dt)
	{
		/* initial pressure distribution. */
		scalar_[index_i] = sin(2 * Pi * pos_[index_i][0]);
	}

protected:
	StdLargeVec<Vecd>& pos_;
	StdLargeVec<Real> scalar_;
	StdLargeVec<Vecd> vector_;
    StdLargeVec<Matd> matrix_;
};

int main(int ac, char* av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	sph_system.setRunParticleRelaxation(true); 
	sph_system.setReloadParticles(false);
#ifdef BOOST_AVAILABLE
	sph_system.handleCommandlineOptions(ac, av); // handle command line arguments
#endif
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	FluidBody body(sph_system, makeShared<Insert>("WaterBody"));
	body.defineAdaptationRatios(1.1, 1.0);
	body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	body.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(1,1,1);
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
		InnerRelation body_inner(body);
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_body_to_vtp(io_environment, { &body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, { &body });
		/* Relaxation method: including based on the 0th and 1st order consistency. */
		InteractionWithUpdate<KernelCorrectionMatrixInner> correction_matrix(body_inner);
		body.addBodyStateForRecording<Matd>("KernelCorrectionMatrix");
		//relax_dynamics::RelaxationStepInnerImplicit relaxation_inner(body_inner, false);
		relax_dynamics::RelaxationStepInnerImplicit<CorrectionMatrixRelaxation> relaxation_inner(body_inner, false);
		PeriodicConditionUsingCellLinkedList periodic_condition_x(body, body.getBodyShapeBounds(), xAxis);
		PeriodicConditionUsingCellLinkedList periodic_condition_y(body, body.getBodyShapeBounds(), yAxis);
		/* Update the relaxation residual. */
		SimpleDynamics<TestingInitialCondition> testing_initial_condition(body);
		InteractionDynamics<relax_dynamics::CheckConsistencyRealization> check_skgc_realization(body_inner,false);
		InteractionDynamics<relax_dynamics::CheckReverseConsistencyRealization> check_rkgc_realization(body_inner, false);
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_ac(body, "ACTERMNORM");
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_as(body, "ASTERMNORM");
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_c(body, "CTERMNORM");
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_ab(body, "ABTERMNORM");
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_ar(body, "ARTERMNORM");
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_b(body, "BTERMNORM");
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_nkgc(body, "NKGCNORM");
		/* Update the kinetic energy for stopping the relaxation. */
		SimpleDynamics<relax_dynamics::UpdateParticleKineticEnergy> update_body_kinetic_energy(body_inner);
		ReduceDynamics<Average<QuantitySummation<Real>>> calculate_body_average_kinetic_energy(body, "ParticleKineticEnergy");
		ReduceDynamics<QuantityMaximum<Real>> calculate_body_maximum_kinetic_energy(body, "ParticleKineticEnergy");
		std::string filefullpath_error_analysis = io_environment.output_folder_ + "/" + "error_analysis.dat";
		std::ofstream out_file_error_analysis(filefullpath_error_analysis.c_str(), std::ios::app);
		//----------------------------------------------------------------------  
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_insert_body_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		periodic_condition_x.update_cell_linked_list_.exec();
		periodic_condition_y.update_cell_linked_list_.exec();
		sph_system.initializeSystemConfigurations();
		testing_initial_condition.exec();
		write_body_to_vtp.writeToFile(0);
		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();
		int ite = 0; //iteration step for the total relaxation step.

		Real body_average_kinetic_energy = 100.0;
		Real body_maximum_kinetic_energy = 100.0;
		Real ac_term = 0;
		Real as_term = 0;
		Real c_term = 0;
		Real ab_term = 0;
		Real ar_term = 0;
		Real b_term = 0;
		Real nkgc_norm = 0;
		
		/* Initial error analysis */
		correction_matrix.exec();
		check_skgc_realization.exec();
		check_rkgc_realization.exec();
		ac_term = calculate_body_average_ac.exec();
		as_term = calculate_body_average_as.exec();
		c_term = calculate_body_average_c.exec();
		ab_term = calculate_body_average_ab.exec();
		ar_term = calculate_body_average_ar.exec();
		b_term = calculate_body_average_b.exec();
		nkgc_norm = calculate_body_average_nkgc.exec();
		out_file_error_analysis << std::fixed << std::setprecision(12) << ite << "   " <<
			body_average_kinetic_energy << " " << nkgc_norm << " " <<
			ac_term << " " << as_term << " " << c_term << " " <<
			ab_term << " " << ar_term << " " << b_term << "\n";

		GlobalStaticVariables::physical_time_ = ite;
		while (body_average_kinetic_energy > 3e-5)
		{
			periodic_condition_x.bounding_.exec();
			periodic_condition_y.bounding_.exec();
			body.updateCellLinkedList();
			periodic_condition_x.update_cell_linked_list_.exec();
			periodic_condition_y.update_cell_linked_list_.exec();
			body_inner.updateConfiguration();

			correction_matrix.exec();
			relaxation_inner.exec();

			periodic_condition_x.bounding_.exec();
			periodic_condition_y.bounding_.exec();
			body.updateCellLinkedList();
			periodic_condition_x.update_cell_linked_list_.exec();
			periodic_condition_y.update_cell_linked_list_.exec();
			body_inner.updateConfiguration();
			ite++;

			if (ite % 100 == 0)
			{
				testing_initial_condition.exec();
				correction_matrix.exec();
				check_skgc_realization.exec();
				check_rkgc_realization.exec();
				update_body_kinetic_energy.exec();
				body_average_kinetic_energy = calculate_body_average_kinetic_energy.exec();
				body_maximum_kinetic_energy = calculate_body_maximum_kinetic_energy.exec();
				ac_term = calculate_body_average_ac.exec();
				as_term = calculate_body_average_as.exec();
				c_term = calculate_body_average_c.exec();
				ab_term = calculate_body_average_ab.exec();
				ar_term = calculate_body_average_ar.exec();
				b_term = calculate_body_average_b.exec();
				nkgc_norm = calculate_body_average_nkgc.exec();
				std::cout << std::fixed << std::setprecision(9) << "The 0th relaxation steps for the body N = " << ite << "\n";
				std::cout << "Body: "
					<< " Average: " << body_average_kinetic_energy
					<< " Maximum: " << body_maximum_kinetic_energy << std::endl;
				out_file_error_analysis << std::fixed << std::setprecision(12) << ite << "   " <<
					body_average_kinetic_energy << " " << nkgc_norm << " " <<
					ac_term << " " << as_term << " " << c_term << " " <<
					ab_term << " " << ar_term << " " << b_term << "\n";
				write_body_to_vtp.writeToFile(ite);
			}
		}
		ite++;

		write_body_to_vtp.writeToFile(ite);
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
	InnerRelation body_inner(body);
	PeriodicConditionUsingCellLinkedList periodic_condition_x(body, body.getBodyShapeBounds(), xAxis);
	PeriodicConditionUsingCellLinkedList periodic_condition_y(body, body.getBodyShapeBounds(), yAxis);
	SimpleDynamics<TestingInitialCondition> testing_initial_condition(body);
	InteractionWithUpdate<KernelCorrectionMatrixInner> correction_matrix(body_inner);
	InteractionDynamics<relax_dynamics::CheckConsistencyRealization> check_consistency_realization(body_inner, false);
	InteractionDynamics<relax_dynamics::CheckL2NormError> check_l2_error(body_inner);
	ReduceDynamics<Average<QuantitySummation<Real>>> calculate_nkgc_error(body, "NKGCError");
	ReduceDynamics<Average<QuantitySummation<Real>>> calculate_skgc_error(body, "SKGCError");
	ReduceDynamics<Average<QuantitySummation<Real>>> calculate_ckgc_error(body, "CKGCError");
	//----------------------------------------------------------------------
	//	Define the methods for I/O operations and observations of the simulation.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_body_states(io_environment, sph_system.real_bodies_);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary.
	//----------------------------------------------------------------------
	/** initialize cell linked lists for all bodies. */
	sph_system.initializeSystemCellLinkedLists();
	/** periodic condition applied after the mesh cell linked list build 
	    up but before the configuration build up. */ 
	periodic_condition_x.update_cell_linked_list_.exec();
	periodic_condition_y.update_cell_linked_list_.exec();
	/** initialize configurations for all bodies. */
	sph_system.initializeSystemConfigurations();
	correction_matrix.exec();
	testing_initial_condition.exec();
	periodic_condition_x.bounding_.exec();
	periodic_condition_y.bounding_.exec();
	body.updateCellLinkedList();
	periodic_condition_x.update_cell_linked_list_.exec();
	periodic_condition_y.update_cell_linked_list_.exec();
	body_inner.updateConfiguration();
	//----------------------------------------------------------------------
	//	Setup for error initialization
	//----------------------------------------------------------------------
	Real nkgc_error = 0.0;
	Real skgc_error = 0.0;
	Real ckgc_error = 0.0;
	//----------------------------------------------------------------------
   //	Main loop starts here.
   //----------------------------------------------------------------------
	std::string filefullpath_error_information = io_environment.output_folder_ + "/" + "error_information.dat";
	std::ofstream out_file_nonopt_temperature(filefullpath_error_information.c_str(), std::ios::app);
	check_l2_error.exec();
	nkgc_error = calculate_nkgc_error.exec();
	skgc_error = calculate_skgc_error.exec();
	ckgc_error = calculate_ckgc_error.exec();

	std::cout << "nkgc error is " << nkgc_error << std::endl;
	std::cout << "skgc error is " << skgc_error << std::endl;
	std::cout << "ckgc error is " << ckgc_error << std::endl;
	write_body_states.writeToFile(0);

	out_file_nonopt_temperature << std::fixed << std::setprecision(12) << "nkgc error is " << "   " << nkgc_error << "\n";
	out_file_nonopt_temperature << std::fixed << std::setprecision(12) << "skgc error is " << "   " << skgc_error << "\n";
	out_file_nonopt_temperature << std::fixed << std::setprecision(12) << "ckgc error is " << "   " << ckgc_error << "\n";
	out_file_nonopt_temperature.close();

	return 0;
}
