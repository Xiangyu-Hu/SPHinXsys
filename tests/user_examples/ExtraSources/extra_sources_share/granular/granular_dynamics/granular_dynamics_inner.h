#ifndef GRANULAR_DYNAMICS_INNER_H
#define GRANULAR_DYNAMICS_INNER_H
#include "fluid_dynamics_inner.h"
#include "fluid_dynamics_inner.hpp"
#include "granular_body.h"
#include "granular_material.h"
#include "granular_particles.h"
#include "constraint_dynamics.h"

namespace SPH
{
	namespace granular_dynamics
	{
		typedef DataDelegateSimple<GranularMaterialParticles> GranularDataSimple;
		typedef DataDelegateInner<GranularMaterialParticles> GranularDataInner;

		class GranularInitialCondition : public LocalDynamics, public GranularDataSimple
		{
		public:
			explicit GranularInitialCondition(SPHBody& sph_body);
			virtual ~GranularInitialCondition() {};

		protected:
			StdLargeVec<Vecd>& pos_, & vel_;
		};


		class GranularAcousticTimeStepSize : public fluid_dynamics::AcousticTimeStepSize
		{
		public:
			explicit GranularAcousticTimeStepSize(SPHBody& sph_body, Real acousticCFL = 0.5);
			virtual ~GranularAcousticTimeStepSize() {};
			Real reduce(size_t index_i, Real dt = 0.0);
			virtual Real outputResult(Real reduced_value) override;
		};

		/**
		 * @class BaseIntegration
		 * @brief Pure abstract base class for all fluid relaxation schemes
		 */
		class BaseRelaxation : public LocalDynamics, public GranularDataInner
		{
		public:
			explicit BaseRelaxation(BaseInnerRelation& inner_relation);
			virtual ~BaseRelaxation() {};

		protected:
			GranularMaterial& granular_material_;
			StdLargeVec<Real>& rho_, & p_, & drho_dt_;
			StdLargeVec<Vecd>& pos_, & vel_, & acc_, & acc_prior_;
		};

		/**
		* @class BaseArtificialStressRelaxation
		*/
		class BaseArtificialStressRelaxation : public BaseRelaxation
		{
		public:
			explicit BaseArtificialStressRelaxation(BaseInnerRelation& inner_relation, Real epsilon = 0.3);
			virtual ~BaseArtificialStressRelaxation() {};
			Matd repulsiveForce(Matd stress_tensor_i, Real rho_i);
		protected:
			Real smoothing_length_, reference_spacing_, epsilon_;
		};

		/**
		* @class ArtificialNormalStressRelaxation
		*/
		class ArtificialNormalShearStressRelaxation : public BaseArtificialStressRelaxation
		{
		public:
			explicit ArtificialNormalShearStressRelaxation(BaseInnerRelation& inner_relation, Real exponent = 4);
			virtual ~ArtificialNormalShearStressRelaxation() {};
			void interaction(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Matd>& shear_stress_;
			StdLargeVec<Vecd>& acc_shear_;
			Real exponent_;
		};


		/**
		* @class ShearStressRelaxation1stHalf
		*/
		class ShearStressRelaxation1stHalf : public BaseRelaxation
		{
		public:
			explicit ShearStressRelaxation1stHalf(BaseInnerRelation& inner_relation);
			virtual ~ShearStressRelaxation1stHalf() {};
			void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Matd>& shear_stress_, & shear_stress_rate_;
			StdLargeVec<Vecd>& acc_shear_;
		};
		/**
		* @class ShearStressRelaxation2ndHalf
		*/
		class ShearStressRelaxation2ndHalf : public BaseRelaxation
		{
		public:

			explicit ShearStressRelaxation2ndHalf(BaseInnerRelation& inner_relation);
			virtual ~ShearStressRelaxation2ndHalf() {};
			//void initialization(size_t index_i, Real dt = 0.0);
			void interaction(size_t index_i, Real dt = 0.0);
			void update(size_t index_i, Real dt = 0.0);
		protected:
			StdLargeVec<Matd>& shear_stress_, & shear_stress_rate_, & velocity_gradient_, & strain_tensor_, & strain_tensor_rate_;
			StdLargeVec<Real>& von_mises_stress_;
		};

		/**
		 * @class BaseMotionConstraint
		 */
		template <class DynamicsIdentifier>
		class BaseMotionConstraint : public BaseLocalDynamics<DynamicsIdentifier>, public GranularDataSimple
		{
		public:
			explicit BaseMotionConstraint(DynamicsIdentifier& identifier)
				: BaseLocalDynamics<DynamicsIdentifier>(identifier), GranularDataSimple(identifier.getSPHBody()),
				pos_(particles_->pos_), pos0_(particles_->pos0_),
				n_(particles_->n_), n0_(particles_->n0_),
				vel_(particles_->vel_), acc_(particles_->acc_) {};

			virtual ~BaseMotionConstraint() {};

		protected:
			StdLargeVec<Vecd>& pos_, & pos0_;
			StdLargeVec<Vecd>& n_, & n0_;
			StdLargeVec<Vecd>& vel_, & acc_;
		};
		/**@class FixConstraint
		 * @brief Constraint with zero velocity.
		 */
		template <class DynamicsIdentifier>
		class FixConstraint : public BaseMotionConstraint<DynamicsIdentifier>
		{
		public:
			explicit FixConstraint(DynamicsIdentifier& identifier)
				: BaseMotionConstraint<DynamicsIdentifier>(identifier) {};
			virtual ~FixConstraint() {};

			void update(size_t index_i, Real dt = 0.0) { this->vel_[index_i] = Vecd::Zero(); };
		};
		using FixBodyConstraint = FixConstraint<SPHBody>;
		using FixBodyPartConstraint = FixConstraint<BodyPartByParticle>;
	}
}
#endif //GRANULAR_DYNAMICS_INNER_H