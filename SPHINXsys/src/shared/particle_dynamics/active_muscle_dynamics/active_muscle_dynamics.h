/**
 * @file 	active_muscle_dynamics.h
 * @brief 	In is file, we declear muscle dynamics which is driven by an external injection of energy. 
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.3.1
 *			Here, we need identify the physical differences between electrophysilogy and active muscle.
 *			The former is on the diffusion and electro-chemical reaction happens in tissue.
 *			The latter is for muscle dynamics which is driven by an external injection of energy.
 *			As the naming of class, function and variables in this code need be based on physics.
 *			We will identify physical differences by properly choosing names.
 *			Xiangyu Hu
 */
#pragma once

#include "all_particle_dynamics.h"
#include "elastic_solid.h"
#include "base_kernel.h"
#include "solid_body.h"
#include "solid_particles.h"

namespace SPH
{
    namespace active_muscle_dynamics
	{
        typedef ParticleDynamicsSimple<SolidBody, ActiveMuscleParticles, ActiveMuscle> ActiveMuscleSimple;
 
		/**
		 * @class ElectroMechanicsInitialCondition
		 * @brief  set initial condition for electro mechanics
		*/
		class ElectroMechanicsInitialCondition : public ActiveMuscleSimple
		{
		protected:
			virtual void Update(size_t index_particle_i, Real dt = 0.0) = 0;
		public:
			ElectroMechanicsInitialCondition(SolidBody *body)
				: ActiveMuscleSimple(body) {};
			virtual ~ElectroMechanicsInitialCondition() {};
		};

		/**@class SpringConstrainMuscleRegion
		 * @brief Constriant a solid body part with a spring force with original position.
		 */
		class SpringConstrainMuscleRegion 
			: public PartDynamicsByParticle<SolidBody, ActiveMuscleParticles, BodyPartByParticle>
		{
		protected:
			Vecd stiffness_;
			virtual Vecd GetAcceleration(Vecd &disp, Real mass);
			virtual void Update(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			SpringConstrainMuscleRegion(SolidBody *body, BodyPartByParticle*body_part)
				: PartDynamicsByParticle<SolidBody, ActiveMuscleParticles, BodyPartByParticle>(body, body_part) {};
			virtual ~SpringConstrainMuscleRegion() {};
			void setUpSpringStiffness(Vecd stiffness){stiffness_ = stiffness;}
		};
		
		/**@class ImposingStress
		 * @brief impose activation stress on a solid body part
		 */
		class ImposingStress
			: public PartDynamicsByParticle<SolidBody, ActiveMuscleParticles, SolidBodyPartForSimbody>
		{
		protected:
			/**
			 * @brief the constrian will be specified by the application
			 */
			virtual Matd getStress(Vecd &pos) = 0;
			virtual void Update(size_t index_particle_i,
				Real dt = 0.0) override;
		public:
			ImposingStress(SolidBody *body, SolidBodyPartForSimbody *body_part)
				: PartDynamicsByParticle<SolidBody, ActiveMuscleParticles, SolidBodyPartForSimbody>(body, body_part) {};
			virtual ~ImposingStress() {};
		};
    }
}