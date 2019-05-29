/**
 * @file    weakly_compressible_fluid_body.h
 * @brief 	This is the class for bodies used for weakly compressible fluid.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_body.h"
#include <fstream>

using namespace std;
namespace SPH {
	/**
	 * @brief Friend Class.
	 */
	class SPHSystem;
	class WeaklyCompressibleFluid;
	class Oldroyd_B_Fluid;
	class WeaklyCompressibleFluidParticles;
	class Oldroyd_B_FluidParticles;
	/**
	 * @class WeaklyCompressibleFluidBody
	 * @brief Declaration of weakly compressible fluid body which is used for fluid and derived from RealBody.
	 */
	class WeaklyCompressibleFluidBody : public RealBody
	{
	public:
		/**
		 * @brief Defaut constructor of ObserverBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] material Material property of this body.
		 * @param[in] weakly_compressible_fluid_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		explicit WeaklyCompressibleFluidBody(SPHSystem &system, string body_name,
			WeaklyCompressibleFluid* material,
			WeaklyCompressibleFluidParticles &weakly_compressible_fluid_particles,
			int refinement_level, ParticlesGeneratorOps op);
		/**
		 * @brief Default destructor.
		 */
		virtual ~WeaklyCompressibleFluidBody() {};

		WeaklyCompressibleFluid* material_;		/**< Material of this body. */
		WeaklyCompressibleFluidParticles &weakly_compressible_fluid_particles_; 	/**< Particles in this body. */

		Real signal_speed_max_;					/**< Maximum signal speed, total kinetic energy. */
		StdVec<StdVec<Vecu>> lists_of_constrained_cells_; /**< List of cells in which particles will be constrianed. */

		/**
		 * @brief Build inner configuration.
		 */
		virtual void BuildInnerConfiguration() override;
		/**
		 * @brief Build contact configuration.
		 */
		virtual void BuildContactConfiguration() override;
		/**
		 * @brief Initial condition for relax body.
		 */
		virtual void InitialCondition() = 0;
		/**
		 * @brief Set a particle at rest for easy initial condition.
		 */
		virtual void SetAllParticleAtRest() override;
		/**
		 * @brief Override the virtual function to output global basic parameters in this body.
		 * @param[in,out] out_file Ofstream of out put data.
		 */
		virtual void GlobalBasicParameters(ofstream &out_file) override;
	};

	/**
	  * @class FluidBodyPart
	  * @brief An auxillariy class for WeaklyCompressiblefluidBody to
	  * indicate a part of the body fixed with location.
	  */
	class FluidBodyPart : public EulerianBodyPart
	{
	protected:
		WeaklyCompressibleFluidBody *fluid_body_;
		Region fluid_body_part_region_;

		virtual void TagBodyPartCells() override;
	public:
		FluidBodyPart(WeaklyCompressibleFluidBody *fluid_body, string fluid_body_part_name);
		virtual~FluidBodyPart() {};
	};

	/**
	 * @class Oldroyd_B_FluidBody
	 * @brief Declaration of Oldroyd B fluid body which is used for Non-newtonian fluid 
	 * and derived from WeaklyCompressibleFluidBody.
	 */
	class Oldroyd_B_FluidBody : public WeaklyCompressibleFluidBody
	{

	public:
		/**
		 * @brief Defaut constructor of ObserverBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] oldroyd_b_material Material property of this body.
		 * @param[in] weakly_compressible_fluid_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		explicit Oldroyd_B_FluidBody(SPHSystem &system, string body_name,
			Oldroyd_B_Fluid* oldroyd_b_material,
			Oldroyd_B_FluidParticles &oldroyd_b_fluid_particles,
			int refinement_level, ParticlesGeneratorOps op);
		/**
		 * @brief Default destructor.
		 */
		virtual ~Oldroyd_B_FluidBody() {};

		Oldroyd_B_Fluid* oldroyd_b_material_;					/**< Material property of this body. */
		Oldroyd_B_FluidParticles &oldroyd_b_fluid_particles_;	/**< Particle of this body. */
		/**
		 * @brief Initial condition for relax body.
		 */
		virtual void InitialCondition() = 0;
		/**
		 * @brief Override the virtual function to output global basic parameters in this body.
		 * @param[in,out] out_file Ofstream of out put data.
		 */
		virtual void GlobalBasicParameters(ofstream &out_file) override;

	};
}