#include "all_solid_dynamics.h"
#include "solid_body.h"
#include "solid_particles.h"
#include "neighborhood.h"
#include "base_kernel.h"
#include "base_data_package.h"
#include "elastic_solid.h"
#include "external_force.h"
#include "cell_linked_list.h"
#include "fluid_particles.h"
#include "weakly_compressible_fluid.h"

#include "Simbody.h"
#include "SimTKmath.h"
#include "SimTKcommon.h"

namespace SPH
{
	namespace solid_dynamics
	{
		//=========================================================================================//
		void UpdateElasticNormalDirection::update(size_t index_i, Real dt)
		{
			Matd &F = F_[index_i];
			// deformation tensor in 2D
			Real F00 = F(0, 0);
			Real F01 = F(0, 1);
			Real F10 = F(1, 0);
			Real F11 = F(1, 1);
			// polar decomposition
			Mat2d R{
				{F00 + F11, F01 - F10},
				{F10 - F01, F00 + F11},
			};

			n_[index_i] = R * n0_[index_i] / sqrt(pow(F00 + F11, 2) + pow(F01 - F10, 2));
		}
		//=================================================================================================//
	}
}
