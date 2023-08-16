/**
 * @file 	relaxation_evolution.cpp
 * @brief 	This is the first case with periodical boundary condition 
            by testing the relaxation with evolution method.
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



int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up the environment of a SPHSystem with global controls.
	//----------------------------------------------------------------------
	SPHSystem sph_system(system_domain_bounds, resolution_ref);
	sph_system.setRunParticleRelaxation(true);
	IOEnvironment io_environment(sph_system);
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	SolidBody body(sph_system, makeShared<Insert>("WaterBody"));
	body.defineBodyLevelSetShape()->writeLevelSet(io_environment);
	body.defineParticlesAndMaterial();
	body.addBodyStateForRecording<Vecd>("Position");
	(!sph_system.RunParticleRelaxation() && sph_system.ReloadParticles())
		? body.generateParticles<ParticleGeneratorReload>(io_environment, body.getName())
		: body.generateParticles<ParticleGeneratorLattice>();
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	InnerRelation insert_body_inner(body);
	//----------------------------------------------------------------------
	//	Run particle relaxation for body-fitted distribution if chosen.
	//----------------------------------------------------------------------
	if (sph_system.RunParticleRelaxation())
	{
		//----------------------------------------------------------------------
		//	Methods used for particle relaxation.
		//----------------------------------------------------------------------
		/** Random reset the insert body particle position. */
		SimpleDynamics<RandomizeParticlePosition> random_insert_body_particles(body);
		/** Write the body state to Vtp file. */
		BodyStatesRecordingToVtp write_insert_body_to_vtp(io_environment, { &body });
		/** Write the particle reload files. */
		ReloadParticleIO write_particle_reload_files(io_environment, { &body });

		/* Relaxation method: including 0th and 1st order consistency. */
		relax_dynamics::RelaxationStepInner relaxation_0th_inner(insert_body_inner, false);
		relax_dynamics::RelaxationStepImplicitInner relaxation_0th_implicit_inner(insert_body_inner, false);
		InteractionDynamics<relax_dynamics::CalculateParticleStress> calculate_particle_stress(insert_body_inner, false);
		relax_dynamics::RelaxationStepByStressInner relaxation_1st_inner(insert_body_inner, false);
		relax_dynamics::RelaxationStepByStressImplicitInner relaxation_1st_implicit_inner(insert_body_inner, false);

		PeriodicConditionUsingCellLinkedList periodic_condition_x(body, body.getBodyShapeBounds(), xAxis);
		PeriodicConditionUsingCellLinkedList periodic_condition_y(body, body.getBodyShapeBounds(), yAxis);

		/* Update relaxation residual. */
		InteractionDynamics<relax_dynamics::CheckCorrectedZeroOrderConsistency> check_corrected_zero_order_consistency(insert_body_inner);
		ReduceAverage<QuantitySummation<Real>> calculate_particle_average_zero_error(body, "corrected_zero_order_error");
		ReduceDynamics<QuantityMaximum<Real>> calculate_particle_maximum_zero_error(body, "corrected_zero_order_error");
		InteractionDynamics<relax_dynamics::CheckCorrectedFirstOrderConsistency> check_corrected_first_order_consistency(insert_body_inner);
		ReduceAverage<QuantitySummation<Real>> calculate_particle_average_first_error(body, "corrected_first_order_error");
		ReduceDynamics<QuantityMaximum<Real>> calculate_particle_maximum_first_error(body, "corrected_first_order_error");
		//----------------------------------------------------------------------  
		//	Particle relaxation starts here.
		//----------------------------------------------------------------------
		random_insert_body_particles.exec(0.25);
		sph_system.initializeSystemCellLinkedLists();
		periodic_condition_x.update_cell_linked_list_.exec();
		periodic_condition_y.update_cell_linked_list_.exec();
		sph_system.initializeSystemConfigurations();
		write_insert_body_to_vtp.writeToFile(0);

		//----------------------------------------------------------------------
		//	Relax particles of the insert body.
		//----------------------------------------------------------------------
		TickCount t1 = TickCount::now();

		int ite = 0; //iteration step for the total relaxation step.

		Real last_zero_maximum_residual = 1;
		Real last_zero_average_residual = 1;
		Real last_first_maximum_residual = 1;
		Real last_first_average_residual = 1;

		Real current_zero_maximum_residual = 1; //maximum zero order consistency residual.
		Real current_zero_average_residual = 1; //average zero order consistency residual.
		Real current_first_maximum_residual = 1; //maximum first order consistency residual.
		Real current_first_average_residual = 1; //average first order consistency residual.

		GlobalStaticVariables::physical_time_ = ite;

		/* The procedure to obtain uniform particle distribution that satisfies the 0th order consistency. */
		while (current_zero_maximum_residual > 0.0001)
		{
			periodic_condition_x.bounding_.exec();
			periodic_condition_y.bounding_.exec();
			body.updateCellLinkedList();
			periodic_condition_x.update_cell_linked_list_.exec();
			periodic_condition_y.update_cell_linked_list_.exec();
			insert_body_inner.updateConfiguration();

			relaxation_0th_implicit_inner.exec(0.1);
			calculate_particle_stress.exec();
			relaxation_1st_implicit_inner.exec(0.1);

			ite++;

			if (ite % 100 == 0)
			{
				check_corrected_zero_order_consistency.exec();
				current_zero_average_residual = calculate_particle_average_zero_error.exec();
				current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();

				std::cout << std::fixed << std::setprecision(9) << "0th relaxation steps for the body N = " << ite << "\n";
				std::cout << "The 0th consistency error: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;
				write_insert_body_to_vtp.writeToFile(ite);
			}
		}

		check_corrected_zero_order_consistency.exec();
		current_zero_average_residual = calculate_particle_average_zero_error.exec();
		current_zero_maximum_residual = calculate_particle_maximum_zero_error.exec();
		std::cout << "The 0th consistency error: maximum = " << current_zero_maximum_residual << "; average = " << current_zero_average_residual << std::endl;

		ite++;
		write_insert_body_to_vtp.writeToFile(ite);
		write_particle_reload_files.writeToFile(0);
		TickCount t2 = TickCount::now();
		TickCount::interval_t tt;
		tt = t2 - t1;
		std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
		return 0;
	}
}
