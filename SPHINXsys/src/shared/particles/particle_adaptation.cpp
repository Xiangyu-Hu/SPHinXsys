/**
 * @file particle_adaptation.cpp
 * @brief Definition of functions declared in particle_adaptation.h
 * @author	Xiangyu Hu and Chi Zhang
 */

#include "particle_adaptation.h"

#include "sph_system.h"
#include "all_kernels.h"
#include "base_body.h"
#include "base_particles.h"
#include "cell_linked_list.h"
#include "level_set.h"

namespace SPH
{
	//=================================================================================================//
	ParticleAdaptation::ParticleAdaptation(Real h_spacing_ratio, Real system_resolution_ratio, Real small_shift_factor)
		: h_spacing_ratio_(h_spacing_ratio),
		  system_resolution_ratio_(system_resolution_ratio),
		  local_refinement_level_(0),
		  local_coarse_level_(local_refinement_level_ / 2),
		  spacing_ref_(0), h_ref_(0),
		  spacing_min_(0), spacing_ratio_min_(1.0),
		  spacing_ratio_max_(1.0), h_ratio_min_(1.0), h_ratio_max_(1.0),
		  number_density_min_(1.0), number_density_max_(1.0),
		  kernel_(new KernelWendlandC2()),
		  sph_body_(nullptr), system_domain_bounds_(),
		  base_particles_(nullptr),
		  small_shift_factor_(small_shift_factor)
		  {};
	//=================================================================================================//
	void ParticleAdaptation::initialize(SPHBody *sph_body)
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
	void ParticleAdaptation::replaceKernel(Kernel *another_kernel)
	{
		delete kernel_;
		kernel_ = another_kernel;
	}
	//=================================================================================================//
	Real ParticleAdaptation::
		RefinedSpacing(Real coarse_particle_spacing, int refinement_level)
	{
		return coarse_particle_spacing / powerN(2.0, refinement_level);
	}
	//=================================================================================================//
	Real ParticleAdaptation::computeReferenceNumberDensity(Vec2d zero, Real h_ratio)
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
	Real ParticleAdaptation::computeReferenceNumberDensity(Vec3d zero, Real h_ratio)
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
	Real ParticleAdaptation::ReferenceNumberDensity()
	{
		return computeReferenceNumberDensity(Vecd(0), 1.0);
	}
	//=================================================================================================//
	Real ParticleAdaptation::probeNumberDensity(Vecd zero, Real h_ratio)
	{
		Real alpha = (h_ratio_max_ - h_ratio) / (h_ratio_max_ - h_ratio_min_ + TinyReal);

		return alpha * number_density_max_ + (1.0 - alpha) * number_density_min_;
	}
	//=================================================================================================//
	void ParticleAdaptation::assignBaseParticles(BaseParticles *base_particles)
	{
		base_particles_ = base_particles;
	}
	//=================================================================================================//
	BaseCellLinkedList *ParticleAdaptation::createCellLinkedList()
	{
		return new CellLinkedList(system_domain_bounds_, kernel_->CutOffRadius(), *sph_body_, *this);
	}
	//=================================================================================================//
	BaseLevelSet *ParticleAdaptation::createLevelSet(ComplexShape &complex_shape)
	{
		return new LevelSet(complex_shape.findBounds(), ReferenceSpacing(), complex_shape, *this, small_shift_factor_);
	}
	//=================================================================================================//
	ParticleWithLocalRefinement::
		ParticleWithLocalRefinement(Real h_spacing_ratio,
									Real system_resolution_ratio, int local_refinement_level)
		: ParticleAdaptation(h_spacing_ratio, system_resolution_ratio)
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
		ParticleAdaptation::assignBaseParticles(base_particles);
		base_particles->registerAVariable<indexScalar, Real>(h_ratio_, "SmoothingLengthRatio", 1.0);
	}
	//=================================================================================================//
	BaseCellLinkedList *ParticleWithLocalRefinement::createCellLinkedList()
	{
		return new MultilevelCellLinkedList(system_domain_bounds_, kernel_->CutOffRadius(),
											getCellLinkedListTotalLevel(),
											MaximumSpacingRatio(), *sph_body_, *this);
	}
	//=================================================================================================//
	BaseLevelSet *ParticleWithLocalRefinement::createLevelSet(ComplexShape &complex_shape)
	{
		return new MultilevelLevelSet(complex_shape.findBounds(),
									  ReferenceSpacing(), getLevelSetTotalLevel(),
									  MaximumSpacingRatio(), complex_shape, *this);
	}
	//=================================================================================================//
	ParticleSpacingByBodyShape::
		ParticleSpacingByBodyShape(Real smoothing_length_ratio,
								   Real system_resolution_ratio, int local_refinement_level)
		: ParticleWithLocalRefinement(smoothing_length_ratio,
									  system_resolution_ratio, local_refinement_level){};
	//=================================================================================================//
	Real ParticleSpacingByBodyShape::getLocalSpacing(ComplexShape &complex_shape, Vecd &position)
	{
		Real phi = fabs(complex_shape.findSignedDistance(position));
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
