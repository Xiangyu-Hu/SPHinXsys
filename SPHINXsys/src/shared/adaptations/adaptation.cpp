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
#include "base_particle_dynamics.h"
#include "mesh_with_data_packages.hpp"
#include "vector_functions.h"

namespace SPH
{
	//=================================================================================================//
	SPHAdaptation::SPHAdaptation(SPHBody *sph_body, Real h_spacing_ratio, Real system_resolution_ratio)
		: sph_body_(sph_body),
		  h_spacing_ratio_(h_spacing_ratio),
		  system_resolution_ratio_(system_resolution_ratio),
		  local_refinement_level_(0),
		  spacing_ref_(sph_body_->getSPHSystem().resolution_ref_ / system_resolution_ratio_),
		  h_ref_(h_spacing_ratio_ * spacing_ref_),
		  kernel_(kernel_ptr_keeper_.createPtr<KernelWendlandC2>(h_ref_)),
		  spacing_min_(this->RefinedSpacing(spacing_ref_, local_refinement_level_)),
		  spacing_ratio_min_(powerN(0.5, local_refinement_level_)),
		  h_ratio_max_(powerN(2.0, local_refinement_level_)),
		  number_density_max_(this->computeReferenceNumberDensity(Vecd(0), h_ratio_max_)),
		  system_domain_bounds_(sph_body->getSPHSystem().system_domain_bounds_){};
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
	void SPHAdaptation::resetAdaptationRatios(Real h_spacing_ratio, Real system_resolution_ratio)
	{
		h_spacing_ratio_ = h_spacing_ratio;
		system_resolution_ratio_ = system_resolution_ratio;
		spacing_ref_ = sph_body_->getSPHSystem().resolution_ref_ / system_resolution_ratio_;
		h_ref_ = h_spacing_ratio_ * spacing_ref_;
		kernel_ = kernel_ptr_keeper_.createPtr<KernelWendlandC2>(h_ref_);
		spacing_min_ = RefinedSpacing(spacing_ref_, local_refinement_level_);
	}
	//=================================================================================================//
	UniquePtr<BaseCellLinkedList> SPHAdaptation::createCellLinkedList()
	{
		return makeUnique<CellLinkedList>(system_domain_bounds_, kernel_->CutOffRadius(), *sph_body_, *this);
	}
	//=================================================================================================//
	UniquePtr<BaseLevelSet> SPHAdaptation::createLevelSet(Shape &shape, Real refinement_ratio)
	{
		// estimate the required mesh levels
		size_t total_levels = (int)log10(MinimumDimension(shape.findBounds()) / ReferenceSpacing()) + 2;
		Real coarsest_spacing = ReferenceSpacing() * powerN(2.0, total_levels - 1);
		MultilevelLevelSet coarser_level_sets(shape.findBounds(), coarsest_spacing / refinement_ratio,
											  total_levels - 1, shape, *this);
		// return the finest level set only
		return makeUnique<RefinedLevelSet>(shape.findBounds(), *coarser_level_sets.getMeshLevels().back(), shape, *this);
	}
	//=================================================================================================//
	ParticleWithLocalRefinement::
		ParticleWithLocalRefinement(SPHBody *sph_body, Real h_spacing_ratio,
									Real system_resolution_ratio, int local_refinement_level)
		: SPHAdaptation(sph_body, h_spacing_ratio, system_resolution_ratio)
	{
		local_refinement_level_ = local_refinement_level;
		spacing_min_ = RefinedSpacing(spacing_ref_, local_refinement_level_);
		spacing_ratio_min_ = powerN(0.5, local_refinement_level_);
		h_ratio_max_ = powerN(2.0, local_refinement_level_);
		number_density_max_ = computeReferenceNumberDensity(Vecd(0), h_ratio_max_);
	}
	//=================================================================================================//
	size_t ParticleWithLocalRefinement::getCellLinkedListTotalLevel()
	{
		return size_t(local_refinement_level_);
	}
	//=================================================================================================//
	size_t ParticleWithLocalRefinement::getLevelSetTotalLevel()
	{
		return getCellLinkedListTotalLevel() + 1;
	}
	//=================================================================================================//
	StdLargeVec<Real> &ParticleWithLocalRefinement::
		registerSmoothingLengthRatio(BaseParticles *base_particles)
	{
		base_particles->registerAVariable(h_ratio_, "SmoothingLengthRatio", 1.0);
		return h_ratio_;
	}
	//=================================================================================================//
	UniquePtr<BaseCellLinkedList> ParticleWithLocalRefinement::createCellLinkedList()
	{
		return makeUnique<MultilevelCellLinkedList>(system_domain_bounds_, kernel_->CutOffRadius(),
													getCellLinkedListTotalLevel(), *sph_body_, *this);
	}
	//=================================================================================================//
	UniquePtr<BaseLevelSet> ParticleWithLocalRefinement::createLevelSet(Shape &shape, Real refinement_ratio)
	{
		return makeUnique<MultilevelLevelSet>(shape.findBounds(), ReferenceSpacing() / refinement_ratio,
											  getLevelSetTotalLevel(), shape, *this);
	}
	//=================================================================================================//
	ParticleSpacingByBodyShape::
		ParticleSpacingByBodyShape(SPHBody *sph_body, Real smoothing_length_ratio,
								   Real system_resolution_ratio, int local_refinement_level)
		: ParticleWithLocalRefinement(sph_body, smoothing_length_ratio,
									  system_resolution_ratio, local_refinement_level){};
	//=================================================================================================//
	Real ParticleSpacingByBodyShape::getLocalSpacing(Shape &shape, const Vecd &position)
	{
		Real phi = fabs(shape.findSignedDistance(position));
		Real ratio_ref = phi / (2.0 * spacing_ref_);
		Real target_ratio = 1.0;
		if (ratio_ref < kernel_->KernelSize())
		{
			Real weight = kernel_->W_1D(ratio_ref);
			target_ratio = weight * spacing_ratio_min_ + (1.0 - weight);
		}
		return target_ratio * spacing_ref_;
	}
	//=================================================================================================//
}
