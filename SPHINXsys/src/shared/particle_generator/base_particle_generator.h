/**
 * @file 	base_particle_generator.h
 * @brief 	This is the base class of particle generator, which generates particles
 * 			with given positions and volumes. The direct generator simply generate
 * 			particle with given position and volume. The lattice generator generate
 * 			at lattice position by check whether the poision is contained by a SPH body.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_data_package.h"

namespace SPH {
	/**
	 * @brief Friend Class.
	 */
	class MeshCellLinkedList;
	class SPHBody;
	/**
	 * @class ParticleGenerator.
	 * @brief Base class for particle gneration.
	 */
	class ParticleGenerator
	{
	public:
		/**
	 	* @brief Default constructor.
	 	* @param[in] 
	 	*/
		ParticleGenerator(SPHBody &sph_body);
		/**
	 	* @brief Default destructor.
	 	*/
		virtual ~ParticleGenerator() {};
		/**
	 	* @brief Create particle for a body.
	 	*/
		virtual void CreateParticles(SPHBody &sph_body) = 0;
	};
	/**
	 * @class ParticleGeneratorDirect
	 * @brief generate particle directly from poistion-and-volume data.
	 */
	class ParticleGeneratorDirect
		: public ParticleGenerator
	{

	public:
		/**
	 	* @brief Default constructor.
	 	* @param[in] 
	 	*/
		ParticleGeneratorDirect(SPHBody &sph_body);
		/**
	 	* @brief Default destructor.
	 	*/
		virtual ~ParticleGeneratorDirect() {};
		/**
	 	* @brief Create particle for a body.
	 	*/
		virtual void CreateParticles(SPHBody &sph_body) override;
	};
	/**
	 * @class ReadRelaxedParticlsFromXmlFile
	 * @brief generate particle directly from poistion-and-volume data.
	 */
	class ReadRelaxedParticlsFromXmlFile : public ParticleGenerator
	{
	public:
		/**
	 	* @brief Default constructor.
	 	* @param[in] 
	 	*/
		ReadRelaxedParticlsFromXmlFile(SPHBody &sph_body);
		/**
	 	* @brief Default destructor.
	 	*/
		virtual ~ReadRelaxedParticlsFromXmlFile() {};
		/**
	 	* @brief Create particle for a body.
	 	*/
		virtual void CreateParticles(SPHBody &sph_body) override;
	};
	/**
	 * @class ReadRestartParticlsFromXmlFile
	 * @brief read restart particles from xml files.
	 */
	class ReadRestartParticlsFromXmlFile : public ParticleGenerator
	{
	public:
		/**
	 	* @brief Default constructor.
	 	* @param[in] 
	 	*/
		ReadRestartParticlsFromXmlFile(SPHBody &sph_body);
		/**
	 	* @brief Default destructor.
	 	*/
		virtual ~ReadRestartParticlsFromXmlFile() {};
		/**
	 	* @brief Create particle for a body.
	 	*/
		virtual void CreateParticles(SPHBody &sph_body) override;
	};
}