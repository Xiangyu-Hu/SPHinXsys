/**
* @file 	collision_dynamics.h
* @brief 	Here, we define the algorithm classes for solid dynamics.
* @details 	We consider here a weakly compresible fluids. .
* @author	Luhui Han, Chi ZHang and Xiangyu Hu
* @version	0.1
*/
#pragma once
#include "all_particle_dynamics.h"
#include "base_material.h"
#include "base_kernel.h"

namespace SPH
{
	namespace solid_dynamics
	{
		typedef ParticleDynamicsContact<SolidBody, SolidParticles, Solid, SolidBody, SolidParticles, Solid> CollisionDynamicsContact;
		typedef ParticleDynamicsReduce<Real, ReduceMax, SolidBody, SolidParticles, Solid> CollisionDynamicsMaximum;
		typedef ParticleDynamicsSimple<SolidBody, SolidParticles, Solid> CollisionDynamicsSimple;

	}
}
