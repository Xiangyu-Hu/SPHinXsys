/**
 * @file 	test_2d_shell_particle_relaxation.cpp
 * @brief 	This is a test for generating 2D shell particles.
 * @author 	Dong Wu and Xiangyu Hu
 */

#include "sphinxsys.h"
using namespace SPH;
/**
 * @brief Basic geometry parameters.
 */
Real radius = 24.5;                                    /** Radius of the inner boundary of the cylinder. */
Real thickness = 1.0;                                  /** Thickness of the cylinder. */
Real radius_mid_surface = radius + thickness / 2.0;    /** Radius of the mid surface. */
Real resolution_ref = 0.5; 			                   /**< Global reference resolution. */
Real level_set_refinement_ratio = resolution_ref / (0.1 * thickness);
Vec2d insert_circle_center(0.0, 0.0);		    	   /** Location of the circle center. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-radius - thickness, -radius - thickness),
	Vec2d(radius + thickness, radius + thickness));
/** For material properties of the solid. */
Real rho0_s = 3.67346939; 			                   /** Normalized density. */
Real Youngs_modulus = 4.32e8;	                       /** Normalized Youngs Modulus. */
Real poisson = 0.0; 			                       /** Poisson ratio. */
/**
* @brief define geometry of SPH bodies
*/
 /** Shell body definition */
class Cylinder : public ThinStructure
{
public:
	Cylinder(SPHSystem& system, const std::string body_name)
		: ThinStructure(system, body_name, makeShared<SPHAdaptation>(1.15, 1.0, 0.75, level_set_refinement_ratio))
	{
		/** Geomtry definition. */
		MultiPolygon multi_polygon;
		multi_polygon.addACircle(insert_circle_center, radius + thickness, 100, ShapeBooleanOps::add);
		multi_polygon.addACircle(insert_circle_center, radius, 100, ShapeBooleanOps::sub);
		MultiPolygonShape multi_polygon_shape(multi_polygon);
		body_shape_.add<LevelSetShape>(this, multi_polygon_shape);
	}
};
//--------------------------------------------------------------------------
//	Main program starts here.
//--------------------------------------------------------------------------
int main(int ac, char* av[])
{
	/** Build up a SPHSystem. */
	SPHSystem system(system_domain_bounds, resolution_ref);
	/** Tag for run particle relaxation for the initial body fitted distribution. */
	system.run_particle_relaxation_ = true;
	/** Handle command line arguments. */
	#ifdef BOOST_AVAILABLE
	system.handleCommandlineOptions(ac, av);
	#endif
	/** output environment. */
	In_Output 	in_output(system);

	/** Creating body, materials and particles. */
	Cylinder cylinder_body(system, "CylinderBody");
	ShellParticles cylinder_body_particles(cylinder_body,
										   makeShared<LinearElasticSolid>(rho0_s, Youngs_modulus, poisson),
										   makeShared<ShellParticleGeneratorLattice>(thickness), thickness);
	cylinder_body_particles.addAVariableToWrite<Vecd>("InitialNormalDirection");
	cylinder_body_particles.addAVariableToWrite<Vecd>("NormalDirection");
	/**
	 * @brief define simple data file input and outputs functions.
	 */
	BodyStatesRecordingToVtp write_real_body_states(in_output, { cylinder_body });
	MeshRecordingToPlt 	write_mesh_cell_linked_list(in_output, cylinder_body, cylinder_body.cell_linked_list_);

	/** Set body contact map
	 *  The contact map gives the data conntections between the bodies
	 *  basically the the range of bodies to build neighbor particle lists
	 */
	BodyRelationInner cylinder_body_inner(cylinder_body);

	/** Random reset the particle position. */
	RandomizePartilePosition  random_cylinder_body_particles(cylinder_body);
	/** A  Physics relaxation step. */
	relax_dynamics::ShellRelaxationStepInner relaxation_step_cylinder_body_inner
	                                                     (cylinder_body_inner, thickness, level_set_refinement_ratio);
	/**
	* @brief 	Particle relaxation starts here.
	*/
	random_cylinder_body_particles.parallel_exec(0.25);
	relaxation_step_cylinder_body_inner.mid_surface_bounding_.parallel_exec();
	write_real_body_states.writeToFile(0.0);
	cylinder_body.updateCellLinkedList();
	write_mesh_cell_linked_list.writeToFile(0.0);

	/** relax particles of the insert body. */
	int ite_p = 0;
	while (ite_p < 1000)
	{
		if (ite_p % 100 == 0)
		{
			std::cout << std::fixed << std::setprecision(9) << "Relaxation steps for the inserted body N = " << ite_p << "\n";
			write_real_body_states.writeToFile(ite_p);
		}
		relaxation_step_cylinder_body_inner.parallel_exec();
		ite_p += 1;
	}
	relaxation_step_cylinder_body_inner.mid_surface_bounding_.calculateNormalDirection();
	write_real_body_states.writeToFile(ite_p);
	std::cout << "The physics relaxation process of the cylinder finish !" << std::endl;

	return 0;
}