#ifndef GRANULAR_DYNAMICS_COMPLEX_H
#define GRANULAR_DYNAMICS_COMPLEX_H

#include "granular_dynamics_inner.h"
#include "granular_dynamics_inner.hpp"
#include "fluid_dynamics_complex.h"
#include "granular_material.h"
#include "granular_particles.h"
#include "granular_body.h"


namespace SPH
{
	namespace granular_dynamics
	{

		typedef DataDelegateContact<GranularMaterialParticles, SolidParticles, DataDelegateEmptyBase>
			GranularWallData;
		typedef DataDelegateContact<GranularMaterialParticles, BaseParticles, DataDelegateEmptyBase>
			GranularContactData;
		typedef DataDelegateContact<GranularMaterialParticles, SolidParticles> FSIContactData;

		/**
		 * @class BaseShearStressRelaxation1stHalfType
		 */
		template <class BaseShearStressRelaxation1stHalfType>
		class BaseShearStressRelaxation1stHalfWithWall : public fluid_dynamics::InteractionWithWall<BaseShearStressRelaxation1stHalfType>
		{
		public:
			template <typename... Args>
			BaseShearStressRelaxation1stHalfWithWall(Args &&...args)
				: fluid_dynamics::InteractionWithWall<BaseShearStressRelaxation1stHalfType>(std::forward<Args>(args)...) {};
			virtual ~BaseShearStressRelaxation1stHalfWithWall() {};
			void interaction(size_t index_i, Real dt = 0.0);

		protected:
			//virtual Vecd computeNonConservativeAcceleration(size_t index_i) override;
		};

		using ShearStressRelaxation1stHalfWithWall = BaseShearStressRelaxation1stHalfWithWall<ShearStressRelaxation1stHalf>;
		

		/**
		 * @class BaseShearStressRelaxation2ndHalfType
		 */
		template <class BaseShearStressRelaxation2ndHalfType>
		class BaseShearStressRelaxation2ndHalfWithWall : public fluid_dynamics::InteractionWithWall<BaseShearStressRelaxation2ndHalfType>
		{
		public:
			template <typename... Args>
			BaseShearStressRelaxation2ndHalfWithWall(Args &&...args)
				: fluid_dynamics::InteractionWithWall<BaseShearStressRelaxation2ndHalfType>(std::forward<Args>(args)...) {};
			virtual ~BaseShearStressRelaxation2ndHalfWithWall() {};
			void interaction(size_t index_i, Real dt = 0.0);
		};

		using ShearStressRelaxation2ndHalfWithWall = BaseShearStressRelaxation2ndHalfWithWall<ShearStressRelaxation2ndHalf>;
	}
}

#endif //GRANULAR_DYNAMICS_COMPLEX_H
