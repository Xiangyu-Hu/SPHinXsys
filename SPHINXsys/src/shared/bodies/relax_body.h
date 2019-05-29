/**
 * @file 	relax_body.h
 * @brief 	This is the class for bodies used for particle relaxation scheme.
 * @author	Chi ZHang and Xiangyu Hu
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
	class MeshCellLinkedList;
	class RelaxBodyParticles;
	class MeshBackground;
	/**
	 * @class RelaxBody
	 * @brief Declaration of relaxbody which is used for particle relaxation and derived from RealBody
	 */
	class RelaxBody : public RealBody
	{
	protected:
		//initialize the background level set, displacement to the body surface
		void InitializeBackgroundMesh();

	public:
		/**
		 * @brief Defaut constructor of ObserverBody.
		 * @param[in] sph_system SPHSystem.
		 * @param[in] body_name Name of Body.
		 * @param[in] relax_particles Particles object.
		 * @param[in] refinement_level Refinement level of this body.
		 * @param[in] op Partciel generator manner.
		 */
		RelaxBody(SPHSystem &sph_system, string body_name,
			RelaxBodyParticles &relax_particles, int refinement_level, ParticlesGeneratorOps op);
		/**
		 * @brief Default destructor.
		 */
		virtual ~RelaxBody() {};

		RelaxBodyParticles &relax_particles_;	/**< Particles in this body. */
		MeshBackground *mesh_background_;		/**< Background mesh.*/

		IndexVector lists_of_singularity_particles_;
		
		/**
		 * @brief Build inner configuration.
		*/
		virtual void BuildInnerConfiguration() override;
		/**
		 * @brief Build contact configuration.
		 */
		virtual void BuildContactConfiguration() override;
		/**
		 * @brief initial condition for relax body.
		 */
		virtual void InitialCondition() override;
		/**
		 * @brief Override the virtual function to output global basic parameters in RelaxBody.
		 * @param[in,out] out_file Ofstream of out put data.
		 */
		virtual void GlobalBasicParameters(ofstream &out_file) override;
	};

	/**
	 * @class SolidBodyPart
	 * @brief A auxillariy class for RelaxBody to
	 * indicate the surface particles
	 */
	class RelaxBodySurface : public LagrangianBodyPart
	{
	protected:
		RelaxBody *relax_body_;

		virtual void TagBodyPartParticles() override;
	public:

		RelaxBodySurface(RelaxBody *relax_body);
		virtual~RelaxBodySurface() {};
	};
}