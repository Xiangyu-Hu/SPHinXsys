/**
 * @file 	heart_reader.cpp
 * @brief 	This is the subroutine of generating particle by reading stl files of leftventricle and mitral valve 
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 */
/** header file and namespace. */
#include "sphinxsys.h"
/** case file to setup the test case */
#include "sphere.h"

using namespace SPH;
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(domain_lower_bound, domain_upper_bound, dp_0,4);
	/** Output sytem */
	In_Output in_output(system);
	/** Creat a Heart body, corresponding material, particles and reaction model. */
	MyPolygonBody *polygon_body = new MyPolygonBody(system, "Polygon", 0, ParticlesGeneratorOps::direct);
	NetworkTree* network_tree 	= new NetworkTree(polygon_body, Point(-1.0, 0.0, 0.0), Point(-0.964, 0.0, 0.266), 0.1, 0.01);
	BodyMaterial* body_material = new BodyMaterial();
	ElasticSolidParticles 	body_particles(polygon_body, body_material);
	/** Write particle data. */
	WriteBodyStatesToPlt 		write_states(in_output, system.real_bodies_);
	write_states.WriteToFile(0.0);

	return 0;
}