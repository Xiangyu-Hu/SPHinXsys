/**
 * @file sph_adaptation.cpp
 * @brief Definition of functions declared in adaptation.h
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "adaptation.h"

#include "sph_system.h"
#include "all_kernels.h"
#include "base_body.h"
#include "base_particles.h"

namespace SPH
{
	//=================================================================================================//
	SPHAdaptation::SPHAdaptation(Real h_spacing_ratio, Real system_resolution_ratio, Real small_shift_factor)
		: h_spacing_ratio_(h_spacing_ratio),
		  system_resolution_ratio_(system_resolution_ratio),
		  local_refinement_level_(0),
		  local_coarse_level_(local_refinement_level_ / 2),
		  spacing_ref_(0), h_ref_(0),
		  spacing_min_(0), spacing_ratio_min_(1.0),
		  spacing_ratio_max_(1.0), h_ratio_min_(1.0), h_ratio_max_(1.0),
		  number_density_min_(1.0), number_density_max_(1.0),
		  kernel_(kernel_ptr_keeper_.createPtr<KernelWendlandC2>()),
		  sph_body_(nullptr), system_domain_bounds_(),
		  base_particles_(nullptr),
		  small_shift_factor_(small_shift_factor)
		  {};
	//=================================================================================================//
	void SPHAdaptation::initialize(SPHBody *sph_body)
	{
		sph_body_ = sph_body;
		system_domain_bounds_ = sph_body_->getSPHSystem().system_domain_bounds_;
		spacing_ref_ = sph_body_->getSPHSystem().resolution_ref_ / system_resolution_ratio_;
		h_ref_ = h_spacing_ratio_ * spacing_ref_;
		kernel_->initialize(h_ref_);
		spacing_min_ = RefinedSpacing(spacing_ref_, local_refinement_level_);
		spacing_ratio_min_ = powerN(0.5, local_refinement_level_);
		spacing_ratio_max_ = powerN(2.0, local_coarse_level_);
		h_ratio_min_ = powerN(0.5, local_coarse_level_);
		h_ratio_max_ = powerN(2.0, local_refinement_level_);
		number_density_min_ = computeReferenceNumberDensity(Vecd(0), h_ratio_max_);
		number_density_max_ = computeReferenceNumberDensity(Vecd(0), h_ratio_min_);
	}
	//=================================================================================================//
	Real SPHAdaptation::
		RefinedSpacing(Real coarse_particle_spacing, int refinement_level)
	{
		return coarse_particle_spacing / powerN(2.0, refinement_level);
	}
	//=================================================================================================//
	Real SPHAdaptation::computeReferenceNumberDensity(Vec2d zero, Real h_ratio)
	{
		Real sigma(0);
		Real cutoff_radius = kernel_->CutOffRadius(h_ratio);
		Real particle_spacing = ReferenceSpacing() / h_ratio;
		int search_depth = int(cutoff_radius / particle_spacing) + 1;
		for (int j = -search_depth; j <= search_depth; ++j)
			for (int i = -search_depth; i <= search_depth; ++i)
			{
				Vec2d particle_location(Real(i) * particle_spacing, Real(j) * particle_spacing);
				Real distance = particle_location.norm();
				if (distance < cutoff_radius)
					sigma += kernel_->W(h_ratio, distance, particle_location);
			}
		return sigma;
	}
	//=================================================================================================//
	Real SPHAdaptation::computeReferenceNumberDensity(Vec3d zero, Real h_ratio)
	{
		Real sigma(0);
		Real cutoff_radius = kernel_->CutOffRadius(h_ratio);
		Real particle_spacing = ReferenceSpacing() / h_ratio;
		int search_depth = int(cutoff_radius / particle_spacing) + 1;
		for (int k = -search_depth; k <= search_depth; ++k)
			for (int j = -search_depth; j <= search_depth; ++j)
				for (int i = -search_depth; i <= search_depth; ++i)
				{
					Vec3d particle_location(Real(i) * particle_spacing,
											Real(j) * particle_spacing, Real(k) * particle_spacing);
					Real distance = particle_location.norm();
					if (distance < cutoff_radius)
						sigma += kernel_->W(h_ratio, distance, particle_location);
				}
		return sigma;
	}
	//=================================================================================================//
	Real SPHAdaptation::ReferenceNumberDensity()
	{
		return computeReferenceNumberDensity(Vecd(0), 1.0);
	}
	//=================================================================================================//
	Real SPHAdaptation::probeNumberDensity(Vecd zero, Real h_ratio)
	{
		Real alpha = (h_ratio_max_ - h_ratio) / (h_ratio_max_ - h_ratio_min_ + TinyReal);

		return alpha * number_density_max_ + (1.0 - alpha) * number_density_min_;
	}
	//=================================================================================================//
	void SPHAdaptation::assignBaseParticles(BaseParticles *base_particles)
	{
		base_particles_ = base_particles;
	}
	//=================================================================================================//
	UniquePtr<BaseCellLinkedList> SPHAdaptation::createCellLinkedList()
	{
		return makeUnique<CellLinkedList>(system_domain_bounds_, kernel_->CutOffRadius(), *sph_body_, *this);
	}
	//=================================================================================================//
	UniquePtr<BaseLevelSet> SPHAdaptation::createLevelSet(Shape &shape)
	{
		return makeUnique<LevelSet>(shape.findBounds(), ReferenceSpacing(), shape, *this);
	}
	//=================================================================================================//
	ParticleWithLocalRefinement::
		ParticleWithLocalRefinement(Real h_spacing_ratio,
									Real system_resolution_ratio, int local_refinement_level)
		: SPHAdaptation(h_spacing_ratio, system_resolution_ratio)
	{
		local_refinement_level_ = local_refinement_level;
		local_coarse_level_ = local_refinement_level_ / 2;
	}
	//=================================================================================================//
	size_t ParticleWithLocalRefinement::getCellLinkedListTotalLevel()
	{
		return size_t(local_coarse_level_ + local_refinement_level_);
	}
	//=================================================================================================//
	size_t ParticleWithLocalRefinement::getLevelSetTotalLevel()
	{
		return getCellLinkedListTotalLevel() + 1;
	}
	//=================================================================================================//
	void ParticleWithLocalRefinement::assignBaseParticles(BaseParticles *base_particles)
	{
		SPHAdaptation::assignBaseParticles(base_particles);
		base_particles->registerAVariable<indexScalar, Real>(h_ratio_, "SmoothingLengthRatio", 1.0);
	}
	//=================================================================================================//
	UniquePtr<BaseCellLinkedList> ParticleWithLocalRefinement::createCellLinkedList()
	{
		return makeUnique<MultilevelCellLinkedList>(system_domain_bounds_, kernel_->CutOffRadius(),
													getCellLinkedListTotalLevel(),
													MaximumSpacingRatio(), *sph_body_, *this);
	}
	//=================================================================================================//
	UniquePtr<BaseLevelSet> ParticleWithLocalRefinement::createLevelSet(Shape &shape)
	{
		return makeUnique<MultilevelLevelSet>(shape.findBounds(),
											  ReferenceSpacing(), getLevelSetTotalLevel(),
											  MaximumSpacingRatio(), shape, *this);
	}
	//=================================================================================================//
	ParticleSpacingByBodyShape::
		ParticleSpacingByBodyShape(Real smoothing_length_ratio,
								   Real system_resolution_ratio, int local_refinement_level)
		: ParticleWithLocalRefinement(smoothing_length_ratio,
									  system_resolution_ratio, local_refinement_level){};
	//=================================================================================================//
	Real ParticleSpacingByBodyShape::getLocalSpacing(Shape &shape, Vecd &position)
	{
		Real phi = fabs(shape.findSignedDistance(position));
		Real ratio_ref = phi / (2.0 * spacing_ref_ * spacing_ratio_max_);
		Real target_ratio = spacing_ratio_max_;
		if (ratio_ref < kernel_->KernelSize())
		{
			Real weight = kernel_->W_1D(ratio_ref);
			target_ratio = weight * spacing_ratio_min_ + (1.0 - weight) * spacing_ratio_max_;
		}
		return target_ratio * spacing_ref_;
	}
	//=================================================================================================//
}
