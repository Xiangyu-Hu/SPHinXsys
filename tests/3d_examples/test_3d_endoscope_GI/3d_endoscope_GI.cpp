/**
* @file 	test_3d_endoscope_GI.cpp
* @brief 	This is a test to simulate an endoscope stirred in the GI system.
* @details	The encoscope will be leaded through esophagus and stomach.
* @author 	Massoud Rezavand, Virtonomy GmbH
*/

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
// std::string full_path_to_geometry = "./input/stomach_esophagus.stl";
// std::string full_path_to_file_endoscope = "./input/endoscope.stl";
std::string full_path_to_geometry = "./input/onlyEsophagus.stl";
std::string full_path_to_file_endoscope = "./input/halfEndoscope.stl";
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Vec3d domain_lower_bound(-20.0, 0.0, 200.);
Vec3d domain_upper_bound(20.0, 80.0, 650.0);
Real dp_0 = 1.5;
Real thickness = 1.0;	
Real level_set_refinement_ratio = dp_0 / (0.2 * thickness);
Real gravity_g = 0.1;
Real initial_ball_speed = 0.8;
Vec3d initial_velocity = initial_ball_speed*Vec3d(0.0, -0.1, -1.0);
//----------------------------------------------------------------------
//	Domain bounds of the system.
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
//----------------------------------------------------------------------
//	For material properties of the solid.
//----------------------------------------------------------------------
Real rho0_s = 1.0e3;
Real Youngs_modulus = 5.0e4;
Real poisson = 0.45;
//----------------------------------------------------------------------
//	Define the body.
//----------------------------------------------------------------------
class GIShellModel : public ThinStructure
{
public:
	GIShellModel(SPHSystem &system, const std::string body_name)
		: ThinStructure(system, body_name, makeShared<SPHAdaptation>(1.15, 1.0, 0.75, level_set_refinement_ratio))
	{
		/** Geometry definition. */
		Vecd translation(0.0, 0.0, 0.0);
		TriangleMeshShapeSTL triangle_mesh_geometry_shape(full_path_to_geometry, translation, 1.0);
		body_shape_.add<LevelSetShape>(this, triangle_mesh_geometry_shape, true, false);
	}
};

class endoscope : public SolidBody
{
public:
	endoscope(SPHSystem &system, const std::string &body_name)
		: SolidBody(system, body_name)
	{
		/** Geometry definition. */
		Vecd translation(0.0, 0.0, 0.0);
		TriangleMeshShapeSTL triangle_mesh_shape_stl(full_path_to_file_endoscope, translation, 1.0);
		body_shape_.add<LevelSetShape>(this, triangle_mesh_shape_stl, true);
	}
};
/** Case dependent boundary condition. */
class InitialVelCondition
	: public solid_dynamics::ElasticDynamicsInitialCondition
{
public:
	InitialVelCondition(SolidBody &body)
		: solid_dynamics::ElasticDynamicsInitialCondition(body) {};
protected:
	void Update(size_t index_i, Real dt) override 
	{
		vel_n_[index_i] = initial_velocity;
	};
};
//--------------------------------------------------------------------------
//	Main program starts here.
//--------------------------------------------------------------------------
int main(int ac, char *av[])
{
	//----------------------------------------------------------------------
	//	Build up a SPHSystem.
	//----------------------------------------------------------------------
	SPHSystem system(system_domain_bounds, dp_0);
	//----------------------------------------------------------------------
	//	Tag for run particle relaxation for the initial body fitted distribution.
	//----------------------------------------------------------------------
	system.run_particle_relaxation_ = false;
	system.reload_particles_ = true;
	//----------------------------------------------------------------------
	//	handle command line arguments.
	//----------------------------------------------------------------------
	#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
	#endif
	//----------------------------------------------------------------------
	//	reauired methdods for the simulation
	//----------------------------------------------------------------------
	/** output environment. */
	In_Output 	in_output(system);
	/** Define external force.*/
	Gravity gravity(Vecd(0.0, -0.1*gravity_g, -gravity_g));
	//----------------------------------------------------------------------
	//	Creating body, materials and particles.
	//----------------------------------------------------------------------
	endoscope endoscope_model(system, "halfEndoscopeModel");
	ElasticSolidParticles endoscope_particles(endoscope_model, 
											  makeShared<NeoHookeanSolid>(rho0_s, Youngs_modulus, poisson));
											  
	std::cout <<"Endoscope reloaded !" << std::endl;

	GIShellModel gi_model(system, "onlyEsophagusShellModel");
	ShellParticles gi_model_particles(gi_model,
									  makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
									  makeShared<ShellParticleGeneratorLattice>(thickness),
									  thickness);
	std::cout <<"EsophagusStomach reloaded !" << std::endl;
	//----------------------------------------------------------------------
	//	Define body relation map.
	//	The contact map gives the topological connections between the bodies.
	//	Basically the the range of bodies to build neighbor particle lists.
	//----------------------------------------------------------------------
	BodyRelationInner endoscope_inner(endoscope_model);
	BodyRelationInner gi_inner(gi_model);
	SolidBodyRelationContact endoscope_gi_contact(endoscope_model, {&gi_model});
	SolidBodyRelationContact gi_endoscope_contact(gi_model, {&endoscope_model});
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition  random_imported_model_particles(gi_model);
	/** A  Physics relaxation step. */
	relax_dynamics::ShellRelaxationStepInner relaxation_step_inner(gi_inner, thickness, level_set_refinement_ratio);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_imported_model_particles.parallel_exec(0.25);
	relaxation_step_inner.mid_surface_bounding_.parallel_exec();
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_p = 0;
	while (ite_p < 1000)
	{
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the esophagus_stomach body N = " << ite_p << "\n";
			// write_imported_model_to_vtp.writeToFile(ite_p);
		}
		relaxation_step_inner.parallel_exec();
		ite_p += 1;
	}
	//----------------------------------------------------------------------
	//	Methods used for particle relaxation.
	//----------------------------------------------------------------------
	RandomizePartilePosition random_endoscope_particles(endoscope_model);
	/** A  Physics relaxation step. */
	relax_dynamics::RelaxationStepInner endoscope_relaxation_step_inner(endoscope_inner, true);
	MeshRecordingToPlt write_endoscope_cell_linked_list(in_output, endoscope_model, endoscope_model.cell_linked_list_);
	//----------------------------------------------------------------------
	//	Particle relaxation starts here.
	//----------------------------------------------------------------------
	random_endoscope_particles.parallel_exec(0.25);
	endoscope_relaxation_step_inner.surface_bounding_.parallel_exec();
	//----------------------------------------------------------------------
	//	Particle relaxation time stepping start here.
	//----------------------------------------------------------------------
	int ite_2 = 0;
	while (ite_2 < 1000)
	{
		endoscope_relaxation_step_inner.parallel_exec();
		ite_2 += 1;
		if (ite_2 % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the endoscope model N = " << ite_2 << "\n";
			// write_endoscope_to_vtp.writeToFile(ite_2);
		}
	}
	//----------------------------------------------------------------------
	//	Define simple file input and outputs functions.
	//----------------------------------------------------------------------
	BodyStatesRecordingToVtp write_system_state_to_vtp(in_output, system.real_bodies_);
	write_system_state_to_vtp.writeToFile(0);
	//----------------------------------------------------------------------
	//	Define the main numerical methods used in the simultion.
	//	Note that there may be data dependence on the constructors of these methods.
	//----------------------------------------------------------------------
	TimeStepInitialization 	endoscope_initialize_timestep(endoscope_model/* , gravity */);
	solid_dynamics::CorrectConfiguration endoscope_corrected_configuration(endoscope_inner);
	solid_dynamics::AcousticTimeStepSize endoscope_get_time_step_size(endoscope_model);
	/** stress relaxation for the endoscopes. */
	solid_dynamics::StressRelaxationFirstHalf endoscope_stress_relaxation_first_half(endoscope_inner);
	solid_dynamics::StressRelaxationSecondHalf endoscope_stress_relaxation_second_half(endoscope_inner);
	/** Algorithms for solid-solid contact. */
	solid_dynamics::ShellContactDensity endoscope_gi_update_contact_density(endoscope_gi_contact);
	solid_dynamics::ContactDensitySummation gi_endoscope_update_contact_density(gi_endoscope_contact);
	solid_dynamics::ContactForce endoscope_compute_solid_contact_forces(endoscope_gi_contact);
	/** initial condition */
	InitialVelCondition initial_velocity(endoscope_model);
	//----------------------------------------------------------------------
	//	Prepare the simulation with cell linked list, configuration
	//	and case specified initial condition if necessary. 
	//----------------------------------------------------------------------
	system.initializeSystemCellLinkedLists();
	system.initializeSystemConfigurations();
	gi_model_particles.initializeNormalDirectionFromBodyShape();
	endoscope_corrected_configuration.parallel_exec();
	initial_velocity.exec();

	/** Main loop. */
	int ite 		= 0;
	Real T0 		= 400.0;
	Real End_Time 	= T0;
	Real D_Time 	= 0.01*T0;
	Real Dt 		= 0.1*D_Time;			
	Real dt 		= 0.0; 	
	//----------------------------------------------------------------------
	//	Statistics for CPU time
	//----------------------------------------------------------------------
	tick_count t1 = tick_count::now();
	tick_count::interval_t interval;
	//----------------------------------------------------------------------
	//	Main loop starts here.
	//----------------------------------------------------------------------
	while (GlobalStaticVariables::physical_time_ < End_Time)
	{
		Real integration_time = 0.0;
		while (integration_time < D_Time) 
		{
			Real relaxation_time = 0.0;
			while (relaxation_time < Dt) 
			{
				endoscope_initialize_timestep.parallel_exec();
				if (ite % 100 == 0) 
				{
					std::cout << "N=" << ite << " Time: "
						<< GlobalStaticVariables::physical_time_ << "	dt: " << dt << "\n";
				}
				endoscope_gi_update_contact_density.parallel_exec();
				gi_endoscope_update_contact_density.parallel_exec();
				endoscope_compute_solid_contact_forces.parallel_exec();
				endoscope_stress_relaxation_first_half.parallel_exec(dt);
				endoscope_stress_relaxation_second_half.parallel_exec(dt);

				endoscope_model.updateCellLinkedList();
				gi_model.updateCellLinkedList();
				endoscope_gi_contact.updateConfiguration();
				gi_endoscope_contact.updateConfiguration();

				ite++;
				Real dt_free = endoscope_get_time_step_size.parallel_exec();
				dt = dt_free;
				relaxation_time += dt;
				integration_time += dt;
				GlobalStaticVariables::physical_time_ += dt;
			}
		}
		tick_count t2 = tick_count::now();
		write_system_state_to_vtp.writeToFile(ite);
		tick_count t3 = tick_count::now();
		interval += t3 - t2;
	}
	tick_count t4 = tick_count::now();

	tick_count::interval_t tt;
	tt = t4 - t1 - interval;
	std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;
	return 0;
}


