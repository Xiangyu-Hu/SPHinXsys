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
#include "mesh_cell_linked_list.h"
#include "level_set.h"

namespace SPH
{
	//=================================================================================================//
	ParticleAdaptation::
		ParticleAdaptation(Real h_spacing_ratio, int global_refinement_level) :
		h_spacing_ratio_(h_spacing_ratio),
		global_refinement_level_(global_refinement_level),
		system_resolution_ratio_(1.0),
		local_refinement_level_(0),
		local_coarse_level_(local_refinement_level_ / 2),
		spacing_ref_(0), Vol_ref_(0), h_ref_(0),
		spacing_min_(0), spacing_ratio_min_(1.0),
		spacing_ratio_max_(1.0), h_ratio_min_(1.0), h_ratio_max_(1.0),
		number_density_min_(1.0), number_density_max_(1.0),
		kernel_(new KernelWendlandC2()), 
		sph_body_(NULL), system_domain_bounds_(),
		base_particles_(NULL) {};
	//=================================================================================================//
	void ParticleAdaptation::initialize(SPHBody* sph_body)
	{
		sph_body_ = sph_body;
		system_domain_bounds_ = sph_body_->getSPHSystem().system_domain_bounds_;
		Real body_resololution = system_resolution_ratio_ * sph_body_->getSPHSystem().resolution_ref_;
		spacing_ref_ = RefinedSpacing(body_resololution, global_refinement_level_);
		Vol_ref_ = powerN(spacing_ref_, Dimensions);
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
	void ParticleAdaptation::replaceKernel(Kernel* another_kernel)
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
		int search_range = int(cutoff_radius / particle_spacing) + 1;
		for (int j = -search_range; j <= search_range; ++j)
			for (int i = -search_range; i <= search_range; ++i)
			{
				Vec2d particle_location(Real(i) * particle_spacing, Real(j) * particle_spacing);
				Real distance = particle_location.norm();
				if (distance < cutoff_radius) sigma += kernel_->W(h_ratio, distance, particle_location);
			}
		return sigma;
	}
	//=================================================================================================//
	Real ParticleAdaptation::computeReferenceNumberDensity(Vec3d zero, Real h_ratio)
	{
		Real sigma(0);
		Real cutoff_radius = kernel_->CutOffRadius(h_ratio);
		Real particle_spacing = ReferenceSpacing() / h_ratio;
		int search_range = int(cutoff_radius / particle_spacing) + 1;
		for (int k = -search_range; k <= search_range; ++k)
			for (int j = -search_range; j <= search_range; ++j)
				for (int i = -search_range; i <= search_range; ++i)
				{
					Vec3d particle_location(Real(i) * particle_spacing,
						Real(j) * particle_spacing, Real(k) * particle_spacing);
					Real distance = particle_location.norm();
					if (distance < cutoff_radius) sigma += kernel_->W(h_ratio, distance, particle_location);
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
	void ParticleAdaptation::assignBaseParticles(BaseParticles* base_particles)
	{
		base_particles_ = base_particles;
	}
	//=================================================================================================//
	BaseMeshCellLinkedList* ParticleAdaptation::createMeshCellLinkedList()
	{
		return new MeshCellLinkedList(*sph_body_, *this, system_domain_bounds_, kernel_->CutOffRadius());
	}
	//=================================================================================================//
	BaseLevelSet* ParticleAdaptation::createLevelSet(ComplexShape& complex_shape)
	{
		return new LevelSet(complex_shape, *this, complex_shape.findBounds(), ReferenceSpacing());
	}
	//=================================================================================================//
	ParticleWithLocalRefinement::ParticleWithLocalRefinement(Real h_spacing_ratio,
		int global_refinement_level, int local_refinement_level) :
		ParticleAdaptation(h_spacing_ratio, global_refinement_level)
	{
		local_refinement_level_ = local_refinement_level;
		local_coarse_level_ = local_refinement_level_ / 2;
	}
	//=================================================================================================//
	size_t ParticleWithLocalRefinement::getMeshCellLinkedListTotalLevel()
	{
		return size_t(local_coarse_level_ + local_refinement_level_);
	}
	//=================================================================================================//
	size_t ParticleWithLocalRefinement::getLevelSetTotalLevel()
	{
		return getMeshCellLinkedListTotalLevel() + 1;
	}
	//=================================================================================================//
	void ParticleWithLocalRefinement::assignBaseParticles(BaseParticles* base_particles)
	{
		ParticleAdaptation::assignBaseParticles(base_particles);
		base_particles->registerAVariable<indexScalar, Real>(h_ratio_, "SmoothingLengthRatio", 1.0);
	}
	//=================================================================================================//
	BaseMeshCellLinkedList* ParticleWithLocalRefinement::createMeshCellLinkedList()
	{
		return new MultilevelMeshCellLinkedList(*sph_body_, *this, system_domain_bounds_, 
			kernel_->CutOffRadius(), getMeshCellLinkedListTotalLevel(), MaximumSpacingRatio());
	}
	//=================================================================================================//
	BaseLevelSet* ParticleWithLocalRefinement::createLevelSet(ComplexShape& complex_shape)
	{
		return new MultilevelLevelSet(complex_shape, *this, complex_shape.findBounds(), 
			ReferenceSpacing(), getLevelSetTotalLevel(), MaximumSpacingRatio());
	}
	//=================================================================================================//
	ParticleSpacingByBodyShape::
		ParticleSpacingByBodyShape(Real smoothing_length_ratio,
			int global_refinement_level, int local_refinement_level) :
		ParticleWithLocalRefinement(smoothing_length_ratio,
			global_refinement_level, local_refinement_level) {};
	//=================================================================================================//
	Real ParticleSpacingByBodyShape::getLocalSpacing(ComplexShape& complex_shape, Vecd& position)
	{
		Real phi = abs(complex_shape.findSignedDistance(position));
		Real ratio_ref = phi / (2.0 * spacing_ref_ * spacing_ratio_max_);
		Real target_ratio = spacing_ratio_max_;
		if (ratio_ref < kernel_->KernelSize()) {
			Real weight = kernel_->W_1D(ratio_ref);
			target_ratio = weight * spacing_ratio_min_ + (1.0 - weight) * spacing_ratio_max_;
		}
		return target_ratio * spacing_ref_;
	}
	//=================================================================================================//
	ShellParticleAdaptation::ShellParticleAdaptation(Real global_avg_thickness)
	: ParticleAdaptation(),
	global_avg_thickness_(global_avg_thickness)
	{
	}
	//=================================================================================================//
	BaseLevelSet* ShellParticleAdaptation::createLevelSet(ComplexShape& complex_shape)
	{
		return new LevelSet(complex_shape, *this, complex_shape.findBounds(), RefinedReferenceSpacing());
	}
	//=================================================================================================//
}
