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
	 * @brief Preclaimed class.
	 */
	class SPHSystem;
	class RelaxBodyParticles;
	class MeshBackground;
	class LevelSet;
	/**
	 * @class RelaxBody
	 * @brief Declaration of relaxbody which is used for particle relaxation and derived from RealBody
	 */
	class RelaxBody : public RealBody
	{
	protected:
		//initialize the background level set, displacement to the body surface
		void InitializeBackgroundMesh(Real mesh_size_ratio = 0.5);
	public:
		RelaxBody(SPHSystem &sph_system, string body_name,
			int refinement_level, ParticlesGeneratorOps op);
		RelaxBody(SPHSystem &sph_system, string body_name,
			int refinement_level,Real smoothing_factor, ParticlesGeneratorOps op);
		virtual ~RelaxBody() {};

		/** Build inner configuration. */
		virtual void BuildInnerConfiguration() override;
		/** Build contact configuration. */
		virtual void BuildContactConfiguration() override;
		/** The pointer to derived class object. */
		virtual RelaxBody* PointToThisObject() override { return this; };
	};

	/**
	 * @class RelaxBodySurface
	 * @brief A auxillariy class for RelaxBody to
	 * indicate the surface particles
	 */
	class RelaxBodySurface : public BodyPartByParticle
	{
	protected:
		RelaxBody *relax_body_;

		virtual void TagBodyPartParticles() override;
	public:

		RelaxBodySurface(RelaxBody *relax_body);
		virtual~RelaxBodySurface() {};
		void updateSurfaceParticles();
	};
	/**
	 * @class RelaxBodySingularity
	 * @brief A auxillariy class for RelaxBody to
	 * indicate the surface particles
	 */
	class RelaxBodySingularity : public BodyPartByParticle
	{
	protected:
		RelaxBody *relax_body_;
		Real threshold_;
		virtual void TagBodyPartParticles() override;
	public:

		RelaxBodySingularity(RelaxBody *relax_body);
		virtual~RelaxBodySingularity() {};
		void updateSingularityParticles(Real threshold);
	};
}