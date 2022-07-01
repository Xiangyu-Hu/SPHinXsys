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
	SPHAdaptation::SPHAdaptation(SPHBody *sph_body, Real h_spacing_ratio, Real system_refinement_ratio)
		: h_spacing_ratio_(h_spacing_ratio),
		  system_refinement_ratio_(system_refinement_ratio),
		  local_refinement_level_(0),
		  spacing_ref_(sph_body->getSPHSystem().resolution_ref_ / system_refinement_ratio_),
		  h_ref_(h_spacing_ratio_ * spacing_ref_),
		  kernel_ptr_(makeUnique<KernelWendlandC2>(h_ref_)),
		  spacing_min_(this->RefinedSpacing(spacing_ref_, local_refinement_level_)),
		  h_ratio_max_(powerN(2.0, local_refinement_level_)),
		  number_density_max_(this->computeReferenceNumberDensity(Vecd(0), h_ratio_max_)){};
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
		Real cutoff_radius = kernel_ptr_->CutOffRadius(h_ratio);
		Real particle_spacing = ReferenceSpacing() / h_ratio;
		int search_depth = int(cutoff_radius / particle_spacing) + 1;
		for (int j = -search_depth; j <= search_depth; ++j)
			for (int i = -search_depth; i <= search_depth; ++i)
			{
				Vec2d particle_location(Real(i) * particle_spacing, Real(j) * particle_spacing);
				Real distance = particle_location.norm();
				if (distance < cutoff_radius)
					sigma += kernel_ptr_->W(h_ratio, distance, particle_location);
			}
		return sigma;
	}
	//=================================================================================================//
	Real SPHAdaptation::computeReferenceNumberDensity(Vec3d zero, Real h_ratio)
	{
		Real sigma(0);
		Real cutoff_radius = kernel_ptr_->CutOffRadius(h_ratio);
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
						sigma += kernel_ptr_->W(h_ratio, distance, particle_location);
				}
		return sigma;
	}
	//=================================================================================================//
	Real SPHAdaptation::ReferenceNumberDensity()
	{
		return computeReferenceNumberDensity(Vecd(0), 1.0);
	}
	//=================================================================================================//
	void SPHAdaptation::resetAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio)
	{
		h_spacing_ratio_ = h_spacing_ratio;
		spacing_ref_ = spacing_ref_ * system_refinement_ratio_ / new_system_refinement_ratio;
		system_refinement_ratio_ = new_system_refinement_ratio;
		h_ref_ = h_spacing_ratio_ * spacing_ref_;
		kernel_ptr_.reset(new KernelWendlandC2(h_ref_));
		spacing_min_ = RefinedSpacing(spacing_ref_, local_refinement_level_);
	}
	//=================================================================================================//
	UniquePtr<BaseCellLinkedList> SPHAdaptation::
		createCellLinkedList(const BoundingBox &domain_bounds, SPHBody &sph_body)
	{
		return makeUnique<CellLinkedList>(domain_bounds, kernel_ptr_->CutOffRadius(), sph_body, *this);
	}
	//=================================================================================================//
	UniquePtr<BaseLevelSet> SPHAdaptation::createLevelSet(Shape &shape, Real refinement_ratio)
	{
		// estimate the required mesh levels
		size_t total_levels = (int)log10(MinimumDimension(shape.getBounds()) / ReferenceSpacing()) + 2;
		Real coarsest_spacing = ReferenceSpacing() * powerN(2.0, total_levels - 1);
		MultilevelLevelSet coarser_level_sets(shape.getBounds(), coarsest_spacing / refinement_ratio,
											  total_levels - 1, shape, *this);
		// return the finest level set only
		return makeUnique<RefinedLevelSet>(shape.getBounds(), *coarser_level_sets.getMeshLevels().back(), shape, *this);
	}
	//=================================================================================================//
	ParticleWithLocalRefinement::
		ParticleWithLocalRefinement(SPHBody *sph_body, Real h_spacing_ratio,
									Real system_refinement_ratio, int local_refinement_level)
		: SPHAdaptation(sph_body, h_spacing_ratio, system_refinement_ratio)
	{
		local_refinement_level_ = local_refinement_level;
		spacing_min_ = RefinedSpacing(spacing_ref_, local_refinement_level_);
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
	UniquePtr<BaseCellLinkedList> ParticleWithLocalRefinement::
		createCellLinkedList(const BoundingBox &domain_bounds, SPHBody &sph_body)
	{
		return makeUnique<MultilevelCellLinkedList>(domain_bounds, kernel_ptr_->CutOffRadius(),
													getCellLinkedListTotalLevel(), sph_body, *this);
	}
	//=================================================================================================//
	UniquePtr<BaseLevelSet> ParticleWithLocalRefinement::createLevelSet(Shape &shape, Real refinement_ratio)
	{
		return makeUnique<MultilevelLevelSet>(shape.getBounds(), ReferenceSpacing() / refinement_ratio,
											  getLevelSetTotalLevel(), shape, *this);
	}
	//=================================================================================================//
	ParticleSpacingByBodyShape::
		ParticleSpacingByBodyShape(SPHBody *sph_body, Real smoothing_length_ratio,
								   Real system_refinement_ratio, int local_refinement_level)
		: ParticleWithLocalRefinement(sph_body, smoothing_length_ratio,
									  system_refinement_ratio, local_refinement_level){};
	//=================================================================================================//
	Real ParticleSpacingByBodyShape::getLocalSpacing(Shape &shape, const Vecd &position)
	{
		Real phi = fabs(shape.findSignedDistance(position));
		Real ratio_ref = phi / (2.0 * spacing_ref_);
		Real target_spacing = spacing_ref_;
		if (ratio_ref < kernel_ptr_->KernelSize())
		{
			Real weight = kernel_ptr_->W_1D(ratio_ref) / kernel_ptr_->W_1D(0.0);
			target_spacing = weight * spacing_min_ + (1.0 - weight) * spacing_ref_;
		}
		return target_spacing;
	}
	//=================================================================================================//
}
