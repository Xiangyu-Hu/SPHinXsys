#include "adaptation.h"

namespace SPH
{
    //=================================================================================================//
	Real SPHAdaptation::computeReferenceNumberDensity(Real h_ratio)
	{
		Real sigma(0);
		Real cutoff_radius = kernel_ptr_->CutOffRadius(h_ratio);
		Real particle_spacing = ReferenceSpacing() / h_ratio;
		int search_depth = int(cutoff_radius / particle_spacing) + 1;
		for (int k = -search_depth; k <= search_depth; ++k)
			for (int j = -search_depth; j <= search_depth; ++j)
				for (int i = -search_depth; i <= search_depth; ++i)
				{
					Vecd particle_location(Real(i) * particle_spacing,
											Real(j) * particle_spacing, Real(k) * particle_spacing);
					Real distance = particle_location.norm();
					if (distance < cutoff_radius)
						sigma += kernel_ptr_->W(h_ratio, distance, particle_location);
				}
		return sigma;
	}
    //=================================================================================================//
}