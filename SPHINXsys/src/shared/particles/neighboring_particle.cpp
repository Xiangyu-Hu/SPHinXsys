/**
 * @file neighboring_particles.cpp
 * @brief Definition of funcitons declared in neighboring_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */
#include "neighboring_particle.h"
#include "base_kernel.h"

namespace SPH
{
	//===============================================================//
	NeighboringParticle::NeighboringParticle(Kernel &kernel, Vecd &r_ij, size_t j_index) 
		: r_ij_(r_ij.norm()), W_ij_(kernel.W(r_ij)), dW_ij_(kernel.dW(r_ij)),
		e_ij_(- normalize(r_ij)), j_(j_index)
	{
		/* Nothing done here! */
	}
	//===============================================================//
	void NeighboringParticle::Reset(Kernel &kernel,	Vecd &r_ij, size_t j_index)
	{
		j_ 		= j_index;
		e_ij_ 	= - normalize(r_ij); /**< e_ij_ pointing from i to j */
		r_ij_ 	= r_ij.norm();
		W_ij_ 	= kernel.W(r_ij);
		dW_ij_ 	= kernel.dW(r_ij);
	}
	//===============================================================//
	ReferenceNeighboringParticle::ReferenceNeighboringParticle(Kernel &kernel, Vecd &r_ij, size_t j_index)
		: r_ij_(r_ij.norm()), e_ij_(normalize(r_ij)), W_ij_(kernel.W(r_ij)), dW_ij_(kernel.dW(r_ij)),
		j_(j_index)
	{
		/* Nothing done here! */
	}
	//===============================================================//
	void ReferenceNeighboringParticle::Reset(Kernel &kernel, Vecd &r_ij, size_t j_index)
	{
		j_ 		= j_index;
		e_ij_ 	= normalize(r_ij); /**< e_ij_ pointing from i to j. */
		r_ij_ 	= r_ij.norm();
		W_ij_ 	= kernel.W(r_ij);
		dW_ij_ 	= kernel.dW(r_ij);
	}
}
