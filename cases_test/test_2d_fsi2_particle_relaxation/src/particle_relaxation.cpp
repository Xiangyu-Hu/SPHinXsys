/**
 * @file 	particle_relaxation.cpp
 * @brief 	Particle relaxation for fsi.
 * @detaisl This is the basic test case for relaxing the particle distribution
 * 			using the Level-set and particle method.
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.1
 */
 /**
   * @brief 	SPHinXsys Library.
   */
#include "sphinxsys.h"
   /**
  * @brief Namespace cite here.
  */
using namespace SPH;

/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 11.0; 								/**< Channel length. */
Real DH = 4.1; 									/**< Channel height. */
Real particle_spacing_ref = 0.1; 				/**< Initial refernece particle spacing. */
Real DLsponge = particle_spacing_ref * 20.0;	/**< Sponge region to impose inflow condition. */
Real BW = particle_spacing_ref * 4.0; 			/**< Extending width for BCs. */

Vec2d insert_circle_center(2.0, 2.0);			/**< Center of cylinder. */
Real insert_circle_radius = 0.5;				/**< Radiu of cylinder. */

Real bh = 0.4*insert_circle_radius;				/**< Height of beam. */
Real bl = 7.0*insert_circle_radius;				/**< Lenght of beam. */
/**
 * @brief 	Define the beam geometry by using the corner of box.
 */
Real hbh = bh / 2.0;
Vec2d BLB(insert_circle_center[0], insert_circle_center[1] - hbh);
Vec2d BLT(insert_circle_center[0], insert_circle_center[1] + hbh);
Vec2d BRB(insert_circle_center[0] + insert_circle_radius + bl,
	insert_circle_center[1] - hbh);
Vec2d BRT(insert_circle_center[0] + insert_circle_radius + bl,
	insert_circle_center[1] + hbh);
/**
* @brief create a water block shape
*/
std::vector<Point> CreatWaterBlockShape()
{
	//geometry
	std::vector<Point> water_block_shape;
	water_block_shape.push_back(Point(-DLsponge, 0.0));
	water_block_shape.push_back(Point(-DLsponge, DH));
	water_block_shape.push_back(Point(DL, DH));
	water_block_shape.push_back(Point(DL, 0.0));
	water_block_shape.push_back(Point(-DLsponge, 0.0));

	return water_block_shape;
}
/**
* @brief create a beam shape
*/
std::vector<Point> CreatBeamShape()
{
	std::vector<Point> beam_shape;
	beam_shape.push_back(BLB);
	beam_shape.push_back(BLT);
	beam_shape.push_back(BRT);
	beam_shape.push_back(BRB);
	beam_shape.push_back(BLB);

	return beam_shape;
}
/**
* @brief create outer wall shape
*/
std::vector<Point> CreatOuterWallShape()
{
	std::vector<Point> outer_wall_shape;
	outer_wall_shape.push_back(Point(-DLsponge - BW, -BW));
	outer_wall_shape.push_back(Point(-DLsponge - BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, DH + BW));
	outer_wall_shape.push_back(Point(DL + BW, -BW));
	outer_wall_shape.push_back(Point(-DLsponge - BW, -BW));

	return outer_wall_shape;
}
/**
* @brief create inner wall shape
*/
std::vector<Point> CreatInnerWallShape()
{
	std::vector<Point> inner_wall_shape;
	inner_wall_shape.push_back(Point(-DLsponge - 2.0*BW, 0.0));
	inner_wall_shape.push_back(Point(-DLsponge - 2.0*BW, DH));
	inner_wall_shape.push_back(Point(DL + 2.0*BW, DH));
	inner_wall_shape.push_back(Point(DL + 2.0*BW, 0.0));
	inner_wall_shape.push_back(Point(-DLsponge - 2.0*BW, 0.0));

	return inner_wall_shape;
}
/**
 * @brief 	Water body defintion.
 */
class RelaxWater : public RelaxBody
{
public:
	RelaxWater(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: RelaxBody(system, body_name, refinement_level, op)
	{
		/** Define the water outer geometry by using the corner of a box.*/
		std::vector<Point> water_block_shape = CreatWaterBlockShape();
		body_region_.add_geometry(new Geometry(water_block_shape), RegionBooleanOps::add);

		/** Define a circle geometry*/
		Geometry *circle_geometry = new Geometry(insert_circle_center, insert_circle_radius, 100);
		body_region_.add_geometry(circle_geometry, RegionBooleanOps::sub);
		std::vector<Point> beam_shape = CreatBeamShape();
		body_region_.add_geometry(new Geometry(beam_shape), RegionBooleanOps::sub);

		/** Finish the region modeling. */
		body_region_.done_modeling();
	}
};
/**
 * @brief 	Wall body defintion.
 */
class Wall : public RelaxBody
{
public:
	Wall(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: RelaxBody(system, body_name, refinement_level, op)
	{
		/** Define a box geomerty by using corner coordinate. */
		std::vector<Point> outer_wall_shape = CreatOuterWallShape();
		Geometry *outer_wall_geometry = new Geometry(outer_wall_shape);
		body_region_.add_geometry(outer_wall_geometry, RegionBooleanOps::add);
		/** Define a box geomerty by using corner coordinate. */
		std::vector<Point> inner_wall_shape = CreatInnerWallShape();
		Geometry *inner_wall_geometry = new Geometry(inner_wall_shape);
		body_region_.add_geometry(inner_wall_geometry, RegionBooleanOps::sub);
		/** Finish the region modeling. */
		body_region_.done_modeling();
	}
};
/**
 * @brief 	Elastic body defintion.
 */
class RelaxSolid : public RelaxBody
{
public:
	RelaxSolid(SPHSystem &system, string body_name,
		int refinement_level, ParticlesGeneratorOps op)
		: RelaxBody(system, body_name, refinement_level, op)
	{	/** Create a circle geomerty. */
		Geometry *circle_geometry = new Geometry(insert_circle_center, insert_circle_radius, 100);
		body_region_.add_geometry(circle_geometry, RegionBooleanOps::add);
		/** Create a box geometry by defining cornore coordinates. */
		std::vector<Point> beam_shape = CreatBeamShape();
		body_region_.add_geometry(new Geometry(beam_shape), RegionBooleanOps::add);
		/** Finish the region modeling. */
		body_region_.done_modeling();
		/** Initilize the background mesh, e.g., Level set data. */
		InitializeBackgroundMesh();
	}
};
/**
 * @brief 	Main program starts here.
 */
int main()
{
	/**
	 * @brief 	Build up context -- a SPHSystem.
	 */
	SPHSystem system(Vec2d(-DLsponge - BW, -BW),
		Vec2d(DL + BW, DH + BW), particle_spacing_ref);
	/**
	 * @brief 	Particles and body creation for water.
	 */
	RelaxWater *relax_water = new RelaxWater(system, "WaterBody",
		0, ParticlesGeneratorOps::lattice);
	RelaxBodyParticles 			relax_water_particles(relax_water);
	/**
	 * @brief 	Particles and body creation for wall.
	 */
	Wall *wall_boundary = new Wall(system, "Wall",
		0, ParticlesGeneratorOps::lattice);
	RelaxBodyParticles 			wall_particles(wall_boundary);
	/**
	 * @brief 	Particles and body creation for Elastic structure.
	 */
	RelaxSolid *relax_solid = new RelaxSolid(system, "InsertedBody",
		1, ParticlesGeneratorOps::lattice);
	RelaxBodyParticles 			relax_solid_particles(relax_solid);
	/**
	 * @brief 	Body contact map.
	 * @details The contact map gives the data conntections between the bodies.
	 * 			Basically the the range of bidies to build neighbor particle lists.
	 */
	SPHBodyTopology body_topology = { { relax_water,{ wall_boundary, relax_solid } },
									  { wall_boundary,{} },{ relax_solid,{ relax_water } } };
	system.SetBodyTopology(&body_topology);
	/** Setting up the simulation. */
	system.SetupSPHSimulation();

	/**
	 * @brief 	Methods used for updating data structure.
	 */
	 /** Update the cell linked list system. */
	ParticleDynamicsCellLinkedList		update_water_block_cell_linked_list(relax_water);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsConfiguration 		update_water_block_configuration(relax_water);
	/** Update the cell linked list of bodies when neccessary. */
	ParticleDynamicsCellLinkedList 		update_relax_solid_cell_linked_list(relax_solid);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsInnerConfiguration 	update_relax_solid_configuration(relax_solid);
	/**
	 * @brief 	Algorithms for periodic BCs.
	 */
	 /** Periodic bounding in x direction. */
	PeriodicBoundingInAxisDirection 	periodic_bounding(relax_water, 0);
	/** Periodic BCs in x direction. */
	PeriodicConditionInAxisDirection 	periodic_condition(relax_water, 0);
	/** Random reset the relax solid particle position. */
	RandomizePartilePosition  random_relax_solid_particles(relax_solid);
	/** Random reset the water particle position. */
	RandomizePartilePosition random_relax_water_particles(relax_water);

	/**
	 * @brief 	Algorithms for particle relaxation.
	 */
	 /** Constraint particles to surface. */
	relax_dynamics::ConstriantSurfaceParticles 		
		constrain_particles_to_gemotery(relax_solid, new RelaxBodySurface(relax_solid));
	/** Compute the time step for physics relaxation. */
	relax_dynamics::GetTimeStepSize 	get_solid_relax_timestep(relax_solid);
	/** Physics relax algorith without contact interactions. */
	relax_dynamics::PhysicsRelaxationInner 			relax_process_for_solid(relax_solid);
	/** Compute the time step for physics relaxation. */
	relax_dynamics::GetTimeStepSize 	get_water_relax_timestep(relax_water);
	/** Physics relax algorith without contact interactions. */
	relax_dynamics::PhysicsRelaxationComplex 		fluid_physics_relaxation(relax_water, { wall_boundary, relax_solid });
	/**
	 * @brief Output.
	 */
	 /** Initialize the output system. */
	In_Output in_output(system);
	/** Write the body state to Vtu file. */
	WriteBodyStatesToVtu 		write_states_to_vtu(in_output, system.real_bodies_);
	/** Write the body state to plt file. */
	WriteBodyStatesToPlt 		write_states_to_plt(in_output, system.real_bodies_);
	/** Write the particle reload files. */
	WriteReloadParticle 		write_particle_reload_files(in_output, system.real_bodies_);
	/** Write mesh data. */
	WriteRelaxBodyMeshToPlt 	write_relax_solid_background_mesh(in_output, relax_solid);
	write_relax_solid_background_mesh.WriteToFile(0.0);
	/**
	 * @brief 	Physics relaxation starts here.
	 */
	 /** Relax the elastic structure. */
	random_relax_solid_particles.parallel_exec();
	random_relax_water_particles.parallel_exec();
	constrain_particles_to_gemotery.parallel_exec();
	write_states_to_plt.WriteToFile(0.0);

	int ite_p = 0;
	Real dt_p = 0.0;
	while (ite_p < 10000)
	{
		dt_p = get_solid_relax_timestep.parallel_exec();

		relax_process_for_solid.parallel_exec(dt_p);
		constrain_particles_to_gemotery.parallel_exec();
		ite_p += 1;

		update_relax_solid_cell_linked_list.parallel_exec();
		update_relax_solid_configuration.parallel_exec();
	}
	std::cout << "The physics relaxation process of inserted body finish !" << std::endl;
	/** relax the water particles. */
	periodic_condition.parallel_exec();
	update_water_block_configuration.parallel_exec();

	int ite_f = 0;
	Real dt_f = 0.0;
	while (ite_f < 10000) {
		dt_f = get_water_relax_timestep.parallel_exec();

		fluid_physics_relaxation.parallel_exec(dt_f);
		ite_f += 1;

		periodic_bounding.parallel_exec();
		update_water_block_cell_linked_list.parallel_exec();
		periodic_condition.parallel_exec();
		update_water_block_configuration.parallel_exec();
	}
	std::cout << "The physics relaxation of fluid process finish !" << std::endl;
	/** Output results. */
	write_states_to_plt.WriteToFile(1.0);
	write_states_to_vtu.WriteToFile(1.0);
	write_particle_reload_files.WriteToFile(0);

	return 0;
}
