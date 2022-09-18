/**
 * @file neighboring_particles.cpp
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "neighborhood.h"

#include "complex_body.h"
#include "base_particles.h"
#include "base_particle_dynamics.h"

namespace SPH
{
	//=================================================================================================//
	void Neighborhood::removeANeighbor(size_t neighbor_n)
	{
		current_size_--;
		j_[neighbor_n] = j_[current_size_];
		W_ij_[neighbor_n] = W_ij_[current_size_];
		dW_ij_[neighbor_n] = dW_ij_[current_size_];
		r_ij_[neighbor_n] = r_ij_[current_size_];
		e_ij_[neighbor_n] = e_ij_[current_size_];
	}
	//=================================================================================================//
	void NeighborBuilder::createNeighbor(Neighborhood &neighborhood,
										  Real &distance, Vecd &displacement, size_t j_index) const
	{
		neighborhood.j_.push_back(j_index);
		neighborhood.W_ij_.push_back(kernel_->W(distance, displacement));
		neighborhood.dW_ij_.push_back(kernel_->dW(distance, displacement));
		neighborhood.r_ij_.push_back(distance);
		neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
		neighborhood.allocated_size_++;
	}
	//=================================================================================================//
	void NeighborBuilder::initializeNeighbor(Neighborhood &neighborhood,
											  Real &distance, Vecd &displacement, size_t j_index) const
	{
		size_t current_size = neighborhood.current_size_;
		neighborhood.j_[current_size] = j_index;
		neighborhood.W_ij_[current_size] = kernel_->W(distance, displacement);
		neighborhood.dW_ij_[current_size] = kernel_->dW(distance, displacement);
		neighborhood.r_ij_[current_size] = distance;
		neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
	}
	//=================================================================================================//
	void NeighborBuilder::createNeighbor(Neighborhood &neighborhood, Real &distance,
										  Vecd &displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min) const
	{
		neighborhood.j_.push_back(j_index);
		Real weight = distance < kernel_->CutOffRadius(i_h_ratio) ? kernel_->W(i_h_ratio, distance, displacement) : 0.0;
		neighborhood.W_ij_.push_back(weight);
		neighborhood.dW_ij_.push_back(kernel_->dW(h_ratio_min, distance, displacement));
		neighborhood.r_ij_.push_back(distance);
		neighborhood.e_ij_.push_back(displacement / (distance + TinyReal));
		neighborhood.allocated_size_++;
	}
	//=================================================================================================//
	void NeighborBuilder::
		initializeNeighbor(Neighborhood &neighborhood, Real &distance,
						   Vecd &displacement, size_t j_index, Real i_h_ratio, Real h_ratio_min) const
	{
		size_t current_size = neighborhood.current_size_;
		neighborhood.j_[current_size] = j_index;
		neighborhood.W_ij_[current_size] = distance < kernel_->CutOffRadius(i_h_ratio)
											   ? kernel_->W(i_h_ratio, distance, displacement)
											   : 0.0;
		neighborhood.dW_ij_[current_size] = kernel_->dW(h_ratio_min, distance, displacement);
		neighborhood.r_ij_[current_size] = distance;
		neighborhood.e_ij_[current_size] = displacement / (distance + TinyReal);
	}
	//=================================================================================================//
	NeighborBuilderInner::NeighborBuilderInner(SPHBody &body) : NeighborBuilder()
	{
		kernel_ = body.sph_adaptation_->getKernel();
	}
	//=================================================================================================//
	void NeighborBuilderInner::operator()(Neighborhood &neighborhood,
										   Vecd &displacement, size_t i_index, size_t j_index) const
	{
		Real distance = displacement.norm();
		if (distance < kernel_->CutOffRadius() && i_index != j_index)
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_
				? createNeighbor(neighborhood, distance, displacement, j_index)
				: initializeNeighbor(neighborhood, distance, displacement, j_index);
			neighborhood.current_size_++;
		}
	};
	//=================================================================================================//
	NeighborBuilderInnerVariableSmoothingLength::
		NeighborBuilderInnerVariableSmoothingLength(SPHBody &body)
		: NeighborBuilder(),
		  h_ratio_(*body.getBaseParticles().getVariableByName<Real>("SmoothingLengthRatio"))
	{
		kernel_ = body.sph_adaptation_->getKernel();
	}
	//=================================================================================================//
	void NeighborBuilderInnerVariableSmoothingLength::
	operator()(Neighborhood &neighborhood, Vecd &displacement, size_t i_index, size_t j_index) const
	{
		Real i_h_ratio = h_ratio_[i_index];
		Real h_ratio_min = SMIN(i_h_ratio, h_ratio_[j_index]);
		Real cutoff_radius = kernel_->CutOffRadius(h_ratio_min);
		Real distance = displacement.norm();
		if (distance < cutoff_radius && i_index != j_index)
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_
				? createNeighbor(neighborhood, distance, displacement, j_index, i_h_ratio, h_ratio_min)
				: initializeNeighbor(neighborhood, distance, displacement, j_index, i_h_ratio, h_ratio_min);
			neighborhood.current_size_++;
		}
	};
	//=================================================================================================//
	NeighborBuilderSelfContact::
		NeighborBuilderSelfContact(SPHBody &body)
		: NeighborBuilder(),
		  pos0_(*body.getBaseParticles().getVariableByName<Vecd>("InitialPosition"))
	{
		kernel_ = body.sph_adaptation_->getKernel();
	}
	//=================================================================================================//
	void NeighborBuilderSelfContact::operator()(Neighborhood &neighborhood,
												 Vecd &displacement, size_t i_index, size_t j_index) const
	{
		Real distance0 = (pos0_[i_index] - pos0_[j_index]).norm();
		Real distance = displacement.norm();
		if (distance < kernel_->CutOffRadius() && distance0 > kernel_->CutOffRadius())
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_
				? createNeighbor(neighborhood, distance, displacement, j_index)
				: initializeNeighbor(neighborhood, distance, displacement, j_index);
			neighborhood.current_size_++;
		}
	};
	//=================================================================================================//
	NeighborBuilderContact::
		NeighborBuilderContact(SPHBody &body, SPHBody &contact_body) : NeighborBuilder()
	{
		Kernel *source_kernel = body.sph_adaptation_->getKernel();
		Kernel *target_kernel = contact_body.sph_adaptation_->getKernel();
		kernel_ = source_kernel->SmoothingLength() > target_kernel->SmoothingLength() ? source_kernel : target_kernel;
	}
	//=================================================================================================//
	void NeighborBuilderContact::operator()(Neighborhood &neighborhood,
											 Vecd &displacement, size_t i_index, size_t j_index) const
	{
		Real distance = displacement.norm();
		if (distance < kernel_->CutOffRadius())
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_
				? createNeighbor(neighborhood, distance, displacement, j_index)
				: initializeNeighbor(neighborhood, distance, displacement, j_index);
			neighborhood.current_size_++;
		}
	};
	//=================================================================================================//
	NeighborBuilderSolidContact::NeighborBuilderSolidContact(SPHBody &body, SPHBody &contact_body)
		: NeighborBuilderContact(body, contact_body)
	{
		Real source_smoothing_length = body.sph_adaptation_->ReferenceSmoothingLength();
		Real target_smoothing_length = contact_body.sph_adaptation_->ReferenceSmoothingLength();
		kernel_ = kernel_keeper_.createPtr<KernelWendlandC2>(0.5 * (source_smoothing_length + target_smoothing_length));
	}
	//=================================================================================================//
	NeighborBuilderContactBodyPart::
		NeighborBuilderContactBodyPart(SPHBody &body, BodyPart &contact_body_part) : NeighborBuilder()
	{
		contact_body_part.getSPHBody().getBaseParticles().registerVariable(part_indicator_, "BodyPartByParticleIndicator");
		Kernel *source_kernel = body.sph_adaptation_->getKernel();
		Kernel *target_kernel = contact_body_part.getSPHBody().sph_adaptation_->getKernel();
		kernel_ = source_kernel->SmoothingLength() > target_kernel->SmoothingLength() ? source_kernel : target_kernel;

		BodyPartByParticle &contact_body_part_by_particle = DynamicCast<BodyPartByParticle>(this, contact_body_part);
		IndexVector part_particles = contact_body_part_by_particle.body_part_particles_;

		for (size_t i = 0; i != part_particles.size(); ++i)
		{
			part_indicator_[part_particles[i]] = 1;
		}
	}
	//=================================================================================================//
	void NeighborBuilderContactBodyPart::operator()(Neighborhood &neighborhood,
													 Vecd &displacement, size_t i_index, size_t j_index) const
	{
		Real distance = displacement.norm();
		if (distance < kernel_->CutOffRadius() && part_indicator_[j_index] == 1)
		{
			neighborhood.current_size_ >= neighborhood.allocated_size_
				? createNeighbor(neighborhood, distance, displacement, j_index)
				: initializeNeighbor(neighborhood, distance, displacement, j_index);
			neighborhood.current_size_++;
		}
	}
	//=================================================================================================//
}
//=================================================================================================//
