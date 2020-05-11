/**
 * @file    solid_body.cpp
 * @brief 	This is the class for bodies used for solid BCs or Elastic structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */

#include "solid_body.h"
#include "sph_system.h"
#include "base_material.h"
#include "mesh_cell_linked_list.h"

namespace SPH {
	//===============================================================//
	SolidBody::SolidBody(SPHSystem &system, string body_name, 
		int refinement_level, ParticlesGeneratorOps op)
		: RealBody(system, body_name, refinement_level, 1.0, op)
	{
	
	}
}