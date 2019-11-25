/**
 * @file 	eletro_mechanics.h
 * @brief 	In is file, we declear the electro mechanics and the corresponing function. 
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.1
 * 			From here, I will denote version a beta, e.g. 0.2.1, other than 0.1 as
 * 			we will introduce cardiac electrophysiology and cardaic mechanics herein.
 * 			Chi Zhang
 */
#pragma once

#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "base_kernel.h"

namespace SPH
{
	class ElectroPhysiology;
    namespace electro_mechanics
	{
        typedef ParticleDynamicsSimple<SolidBody, MuscleParticles, Muscle> ElectroMechanicsSimple;
        typedef ParticleDynamicsReduce<Real, ReduceMin, SolidBody, MuscleParticles, Muscle> ElectroMechanicsMinimum;
        typedef ParticleDynamicsInner1Level<SolidBody, MuscleParticles, Muscle> ElectroMechanicsInnver1Level;
        /**
		 * @class  OffsetInitialParticlePosition
		 * @brief  set initial condition for a cardiac muscle body
		*/
		class OffsetInitialParticlePosition : public ElectroMechanicsSimple
		{
			Vecd offset_;
		protected:
			//default for set all particle at rest
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			OffsetInitialParticlePosition(SolidBody *body, Vecd offset)
				: ElectroMechanicsSimple(body), offset_(offset){};
			virtual ~OffsetInitialParticlePosition() {};
		};
        /**
		 * @class ElectroMechanicsInitialCondition
		 * @brief  set initial condition for electro mechanics
		*/
		class ElectroMechanicsInitialCondition : public ElectroMechanicsSimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ElectroMechanicsInitialCondition(SolidBody *body)
				: ElectroMechanicsSimple(body) {};
			virtual ~ElectroMechanicsInitialCondition() {};
		};
        /**
		 * @class computeActiveContractionStress
		 * @brief compute the active contraction stree T_a for cardiomyocite
		*/
		class computeActiveContractionStress : public ElectroMechanicsSimple
		{
		protected:
            ElectroPhysiology* reaction_model_;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			computeActiveContractionStress(SolidBody *body);
			virtual ~computeActiveContractionStress() {};
		};
		/**
		 * @class computeLinearActiveContractionStress
		 * @brief This is linear active contraction stree, Ta = alpha * voltage_n_
		*/
		class computeLinearActiveContractionStress : public ElectroMechanicsSimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			computeLinearActiveContractionStress(SolidBody *body);
			virtual ~computeLinearActiveContractionStress() {};
		};
		/**
		* @class ActivePassiveStressRelaxationFirstStep
		* @brief computing Active-Passive stress relaxation process by verlet time stepping
		* This is the first step
		*/
		class ActivePassiveStressRelaxationFirstStep : public ElectroMechanicsInnver1Level
		{
		protected:
			Real numerical_viscosity_;

			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ActivePassiveStressRelaxationFirstStep(SolidBody *body) : ElectroMechanicsInnver1Level(body) 
			{
				numerical_viscosity_ = 0.5 * material_->rho_0_ * material_->c_0_ * 
					body_->kernel_->GetSmoothingLength(); 
			};
			virtual ~ActivePassiveStressRelaxationFirstStep() {};
		};

		/**
		* @class ActivePassiveStressRelaxationSecondStep
		* @brief computing Active-Passive stress relaxation process by verlet time stepping
		* This is the second step
		*/
		class ActivePassiveStressRelaxationSecondStep : public ElectroMechanicsInnver1Level
		{
		protected:
			virtual void Initialization(size_t index_particle_i, Real dt = 0.0) override;
			virtual void InnerInteraction(size_t index_particle_i, Real dt = 0.0) override;
			virtual void Update(size_t index_particle_i, Real dt = 0.0) override;
		public:
			ActivePassiveStressRelaxationSecondStep(SolidBody *body) : ElectroMechanicsInnver1Level(body) {};
			virtual ~ActivePassiveStressRelaxationSecondStep() {};
		};
		/**@class SpringConstrainMuscleRegion
		 * @brief Constriant a solid body part with a spring force with original position.
		 */
		class SpringConstrainMuscleRegion 
			: public ConstraintByParticle<SolidBody, MuscleParticles, SolidBodyPart>
		{
		protected:
			Vecd stiffness_;
			virtual Vecd GetAcceleration(Vecd &disp, Real mass);
			virtual void ConstraintAParticle(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			SpringConstrainMuscleRegion(SolidBody *body, SolidBodyPart *body_part)
				: ConstraintByParticle<SolidBody, MuscleParticles, SolidBodyPart>(body, body_part) {};
			virtual ~SpringConstrainMuscleRegion() {};
			void setUpSpringStiffness(Vecd stiffness){stiffness_ = stiffness;}
		};
    }
}