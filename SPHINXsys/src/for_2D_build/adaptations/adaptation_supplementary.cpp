#include "adaptation.h"

namespace SPH
{
    //=================================================================================================//
	Real SPHAdaptation::computeReferenceNumberDensity(Vec2d zero)
	{
		Real sigma(0);
		Real cutoff_radius = kernel_ptr_->CutOffRadius();
		Real particle_spacing = ReferenceSpacing();
		int search_depth = int(cutoff_radius / particle_spacing) + 1;
		for (int j = -search_depth; j <= search_depth; ++j)
			for (int i = -search_depth; i <= search_depth; ++i)
			{
				Vec2d particle_location(Real(i) * particle_spacing, Real(j) * particle_spacing);
				Real distance = particle_location.norm();
				if (distance < cutoff_radius)
					sigma += kernel_ptr_->W(distance, particle_location);
			}
		return sigma;
	}
    //=================================================================================================//
}
