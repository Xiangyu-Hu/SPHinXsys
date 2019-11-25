/**
 * @file neighboring_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 * @version	0.1
 */

#include "neighboring_particle.h"
#include "base_kernel.h"
//=================================================================================================//
namespace SPH
{
//=================================================================================================//
	NeighboringParticle::NeighboringParticle(Kernel &kernel, Vecd &r_ij, size_t j_index) 
		: r_ij_(r_ij.norm()), W_ij_(kernel.W(r_ij)), dW_ij_(kernel.dW(r_ij)),
		e_ij_(normalize(r_ij)), j_(j_index)
	{
		/* Nothing done here! */
	}
//=================================================================================================//
	void NeighboringParticle::Reset(Kernel &kernel,	Vecd &r_ij, size_t j_index)
	{
		j_ 		= j_index;
		e_ij_ 	= normalize(r_ij);
		r_ij_ 	= r_ij.norm();
		W_ij_ 	= kernel.W(r_ij);
		dW_ij_ 	= kernel.dW(r_ij);
	}
//=================================================================================================//
	ReferenceNeighboringParticle::ReferenceNeighboringParticle(Kernel &kernel, Vecd &r_ij, size_t j_index)
		: r_ij_(r_ij.norm()), e_ij_(normalize(r_ij)), W_ij_(kernel.W(r_ij)), dW_ij_(kernel.dW(r_ij)),
		j_(j_index)
	{
		/* Nothing done here! */
	}
//=================================================================================================//
	void ReferenceNeighboringParticle::Reset(Kernel &kernel, Vecd &r_ij, size_t j_index)
	{
		j_ 		= j_index;
		e_ij_ 	= normalize(r_ij); /**< e_ij_ pointing from j to i. */
		r_ij_ 	= r_ij.norm();
		W_ij_ 	= kernel.W(r_ij);
		dW_ij_ 	= kernel.dW(r_ij);
	}
//=================================================================================================//
	ReferenceNeighboringParticleDiffusion::ReferenceNeighboringParticleDiffusion(size_t j_index, Matd &diff_ij)
		: j_(j_index), diffusion_ij_(diff_ij)
	{
		/* Nothing done here! */
	}
//=================================================================================================//
	void ReferenceNeighboringParticleDiffusion::Reset(size_t j_index, Matd &diff_ij)
	{ 
		j_ 		= j_index;
		diffusion_ij_ = diff_ij;
	}
//=================================================================================================//
}
//=================================================================================================//
