/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
 * @file 	in_output.h
 * @brief 	Classes for input and output functions.
 * @author	Chi Zhang and Xiangyu Hu
 * @version	0.1
 */

#pragma once
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_data_package.h"
#include "sph_data_conainers.h"
#include "all_physical_dynamics.h"
 
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"

#include <fstream>
/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

namespace SPH {
	/**
	 * @class In_Output
	 * @brief The base class which defines folders for output, 
	 * restart and particle reload folders.
	 */
	class In_Output
	{
	public:
		In_Output(SPHSystem &sph_system);
		virtual ~In_Output() {};

		SPHSystem &sph_system_;
		std::string output_folder_;
		std::string restart_folder_;
		std::string reload_folder_;
		std::string restart_step_;
	};

	/**
	 * @class BodyStatesIO
	 * @brief base class for write and read body states.
	 */
	class BodyStatesIO
	{
	protected:
		In_Output &in_output_;
		SPHBody *body_;
		SPHBodyVector bodies_;
	public:
		BodyStatesIO(In_Output &in_output, SPHBody *body)
			: in_output_(in_output), body_(body) {};
		BodyStatesIO(In_Output& in_output, SPHBodyVector bodies)
			: in_output_(in_output), bodies_(bodies), body_(bodies[0]) {};
		virtual ~BodyStatesIO() {};
	};

	/**
	 * @class WriteBodyStates
	 * @brief base class for write body states.
	 */
	class WriteBodyStates : public BodyStatesIO
	{
	public:
		WriteBodyStates(In_Output &in_output, SPHBody *body)
			: BodyStatesIO(in_output, body) {};
		WriteBodyStates(In_Output &in_output, SPHBodyVector bodies)
			: BodyStatesIO(in_output, bodies) {};
		virtual ~WriteBodyStates() {};

		virtual void WriteToFile(Real time) = 0;
	};

	/**
	 * @class ReadBodyStates
	 * @brief base class for read body states.
	 */
	class ReadBodyStates : public BodyStatesIO
	{
	public:
		ReadBodyStates(In_Output &in_output, SPHBody *body)
			: BodyStatesIO(in_output, body) {};
		ReadBodyStates(In_Output &in_output, SPHBodyVector bodies)
			: BodyStatesIO(in_output, bodies) {};
		virtual ~ReadBodyStates() {};

		virtual void ReadFromFile(size_t iteration_step) = 0;
	};
	/**
	 * @class SimBodyStatesIO
	 * @brief base class for write and read SimBody states.
	 */
	template<class MobilizedBodyType>
	class SimBodyStatesIO
	{
	protected:
		In_Output &in_output_;
		SimTK::RungeKuttaMersonIntegrator& integ_;
		MobilizedBodyType& mobody_;
	public:
		SimBodyStatesIO(In_Output& in_output, SimTK::RungeKuttaMersonIntegrator& integ, MobilizedBodyType& mobody)
			: in_output_(in_output),integ_(integ), mobody_(mobody) {};
		virtual ~SimBodyStatesIO() {};
	};

	/**
	 * @class WriteSimBodyStates
	 * @brief base class for write SimBody states.
	 */
	template<class MobilizedBodyType>
	class WriteSimBodyStates : public SimBodyStatesIO<MobilizedBodyType>
	{
	public:
		WriteSimBodyStates(In_Output& in_output, SimTK::RungeKuttaMersonIntegrator& integ, MobilizedBodyType& mobody)
			: SimBodyStatesIO<MobilizedBodyType>(in_output,integ, mobody) {};
		virtual ~WriteSimBodyStates() {};

		virtual void WriteToFile(Real time) = 0;
	};

	/**
	 * @class ReadSimBodyStates
	 * @brief base class for read SimBody states.
	 */
	template<class MobilizedBodyType>
	class ReadSimBodyStates : public SimBodyStatesIO<MobilizedBodyType>
	{
	public:
		ReadSimBodyStates(In_Output& in_output, MobilizedBodyType* mobody)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, mobody) {};
		ReadSimBodyStates(In_Output& in_output, StdVec<MobilizedBodyType*> mobodies)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, mobodies) {};
		virtual ~ReadSimBodyStates() {};

		virtual void ReadFromFile(size_t iteration_step) = 0;
	};

	/**
	 * @class WriteBodyStatesToVtu
	 * @brief  Write files for bodies
	 * the output file is VTK XML format can visualized by ParaView
	 * the data type vtkUnstructedGrid
	 */
	class WriteBodyStatesToVtu : public WriteBodyStates
	{
	public:
		WriteBodyStatesToVtu(In_Output& in_output, SPHBodyVector bodies)
			: WriteBodyStates(in_output, bodies) {};
		virtual ~WriteBodyStatesToVtu() {};

		virtual void WriteToFile(Real time) override;
	};
	
	/**
	 * @class WriteBodyStatesToPlt
	 * @brief  Write files for bodies
	 * the output file is dat format can visualized by TecPlot
	 */
	class WriteBodyStatesToPlt : public WriteBodyStates
	{
	public:
		WriteBodyStatesToPlt(In_Output& in_output, SPHBodyVector bodies)
			: WriteBodyStates(in_output, bodies) {};;
		virtual ~WriteBodyStatesToPlt() {};

		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteToVtuIfVelocityOutOfBound
	 * @brief  output body sates if particle velocity is
	 * out of a bound
	 */
	class WriteToVtuIfVelocityOutOfBound
		: public WriteBodyStatesToVtu
	{
	protected:
		StdVec<VelocityBoundCheck *> check_bodies_;
	public:
		WriteToVtuIfVelocityOutOfBound(In_Output& in_output,
			SPHBodyVector bodies, Real velocity_bound);
		virtual ~WriteToVtuIfVelocityOutOfBound() {};

		bool out_of_bound_;
		virtual void WriteToFile(Real time) override;
	};

	/**
	 * @class WriteBodyMeshToPlt
	 * @brief  write the background mesh data for relax body
	 */
	class WriteBodyMeshToPlt : public WriteBodyStates
	{
	protected:
		std::string filefullpath_;
	public:
		WriteBodyMeshToPlt(In_Output& in_output, SPHBody *body);
		virtual ~WriteBodyMeshToPlt() {};

		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class WriteAnObservedQuantity
	 * @brief write files for observed quantity
	 */
	template <class DataType, class TargetParticlesType, class TargetDataType,
		StdLargeVec<TargetDataType> TargetParticlesType:: * TrgtDataMemPtr, DataType TargetDataType:: * TrgtMemPtr>
	class WriteAnObservedQuantity : public WriteBodyStates,
		public observer_dynamics::ObservingAQuantity<DataType, TargetParticlesType, TargetDataType, TrgtDataMemPtr, TrgtMemPtr>
	{
	protected:
		SPHBody* observer_;
		std::string filefullpath_;

		void writeFileHead(std::ofstream& out_file, Real& observed_quantity, string quantity_name, size_t i) {
			out_file << "  " << quantity_name << "[" << i << "]" << " ";
		};
		void writeDataToFile(std::ofstream& out_file, Real& observed_quantity) {
			out_file << "  " << observed_quantity << " ";
		};

		void writeFileHead(std::ofstream& out_file, Vecd& observed_quantity, string quantity_name, size_t i) {
			for (int j = 0; j < observed_quantity.size(); ++j)
				out_file << "  " << quantity_name <<"[" << i << "][" << j << "]" << " ";
		};
		void writeDataToFile(std::ofstream& out_file, Vecd& observed_quantity) {
			for (int j = 0; j < observed_quantity.size(); ++j)
				out_file << "  " << observed_quantity[j] << " ";
		};

	public:
		WriteAnObservedQuantity(string quantity_name, In_Output& in_output, SPHBodyContactRelation* body_contact_relation)
			: WriteBodyStates(in_output, body_contact_relation->body_), observer_(body_contact_relation->body_),
			observer_dynamics::ObservingAQuantity<DataType, TargetParticlesType, TargetDataType, TrgtDataMemPtr, TrgtMemPtr>(body_contact_relation)
		{
			filefullpath_ = in_output_.output_folder_ + "/" + observer_->GetBodyName()
				+ "_" + quantity_name + "_" + in_output_.restart_step_ + ".dat";
			std::ofstream out_file(filefullpath_.c_str(), ios::app);
			out_file << "run_time" << "   ";
			for (size_t i = 0; i != observer_->number_of_particles_; ++i)
			{
				writeFileHead(out_file, this->observed_quantities_[i], quantity_name, i);
			}
			out_file << "\n";
			out_file.close();
		};
		virtual ~WriteAnObservedQuantity() {};

		virtual void WriteToFile(Real time = 0.0) override 
		{
			this->parallel_exec();
			std::ofstream out_file(filefullpath_.c_str(), ios::app);
			out_file << time << "   ";
			for (size_t i = 0; i != observer_->number_of_particles_; ++i)
			{
				writeDataToFile(out_file, this->observed_quantities_[i]);
			}
			out_file << "\n";
			out_file.close();
		};
	};

	/**
 * @class WriteObservedDiffusionReactionQuantity
 * @brief write the observed diffusion and reaction quantity to files.
 */
	template <class DiffusionReactionParticlesType>
	class WriteObservedDiffusionReactionQuantity
		: public WriteBodyStates,
		public observer_dynamics::ObservingADiffusionReactionQuantity<DiffusionReactionParticlesType>
	{
	protected:
		SPHBody* observer_;
		std::string filefullpath_;
	public:
		/** Constructor and Destructor. */
		WriteObservedDiffusionReactionQuantity(string species_name, In_Output& in_output, SPHBodyContactRelation* body_contact_relation)
			: WriteBodyStates(in_output, body_contact_relation->body_), observer_(body_contact_relation->body_),
			observer_dynamics::ObservingADiffusionReactionQuantity<DiffusionReactionParticlesType>(species_name, body_contact_relation)
		{
			filefullpath_ = in_output_.output_folder_ + "/" + observer_->GetBodyName()
				+ "_" + species_name + "_" + in_output_.restart_step_ + ".dat";
			std::ofstream out_file(filefullpath_.c_str(), ios::app);
			out_file << "run_time" << "   ";
			for (size_t i = 0; i != observer_->number_of_particles_; ++i)
			{
				out_file << "  " << species_name << "[" << i << "]" << " ";
			}
			out_file << "\n";
			out_file.close();
		};

		virtual ~WriteObservedDiffusionReactionQuantity() {};
		/**
		 * @brief Output data to files.
		 * @param[in] time Physical time.
		 */
		virtual void WriteToFile(Real time) override 
		{
			this->parallel_exec();
			std::ofstream out_file(filefullpath_.c_str(), ios::app);
			out_file << time << "   ";
			for (size_t i = 0; i != observer_->number_of_particles_; ++i)
			{
				out_file << "  " << this->observed_quantities_[i] << " ";
			}
			out_file << "\n";
			out_file.close();
		};
	};

	/**
	 * @class WriteTotalMechanicalEnergy
	 * @brief write files for the total mechanical energy of a weakly compressible fluid body
	 */
	class WriteTotalMechanicalEnergy 
		: public WriteBodyStates, public fluid_dynamics::TotalMechanicalEnergy
	{
	protected:
		std::string filefullpath_;
	public:
		WriteTotalMechanicalEnergy(In_Output& in_output, FluidBody* water_block, Gravity* gravity);
		virtual ~WriteTotalMechanicalEnergy() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class WriteMaximumSpeed
	 * @brief write files for the maximum speed within the body
	 */
	class WriteMaximumSpeed : public WriteBodyStates, public MaximumSpeed
	{
	protected:
		std::string filefullpath_;
	public:
		WriteMaximumSpeed(In_Output& in_output, SPHBody* sph_body);
		virtual ~WriteMaximumSpeed() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class WriteTotalViscousForceOnSolid
	 * @brief write total viscous force acting a solid body
	 */
	class WriteTotalViscousForceOnSolid
		: public WriteBodyStates, public solid_dynamics::TotalViscousForceOnSolid
	{
	protected:
		size_t dimension_;
		std::string filefullpath_;
	public:
		WriteTotalViscousForceOnSolid(In_Output& in_output, SolidBody *solid_body);
		virtual ~WriteTotalViscousForceOnSolid() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};

	/** 
	 * @class WriteTotalForceOnSolid
	 * @brief Write total force acting a solid body.
	 */
	class WriteTotalForceOnSolid
		: public WriteBodyStates, public solid_dynamics::TotalForceOnSolid
	{
	protected:
		size_t dimension_;
		std::string filefullpath_;
	public:
		WriteTotalForceOnSolid(In_Output& in_output, SolidBody *solid_body);
		virtual ~WriteTotalForceOnSolid() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class WriteUpperFrontInXDirection
	 * @brief write files for water front in free surface flow
	 */
	class WriteUpperFrontInXDirection
		: public WriteBodyStates, public UpperFrontInXDirection
	{
	protected:
		std::string filefullpath_;
	public:
		WriteUpperFrontInXDirection(In_Output& in_output, SPHBody* body);
		virtual ~WriteUpperFrontInXDirection() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class ReloadParticleIO
	 * @brief For write  and read particle reload.
	 */
	class ReloadParticleIO
	{
	protected:
		StdVec<std::string> file_paths_;

	public:
		ReloadParticleIO(In_Output& in_output, SPHBodyVector bodies);
		virtual ~ReloadParticleIO() {};
	};

	/**
	  * @class WriteReloadParticle
	  * @brief Write the reload particles file in XML format.
	  */
	class WriteReloadParticle : public ReloadParticleIO, public WriteBodyStates
	{
	public:
		WriteReloadParticle(In_Output& in_output, SPHBodyVector bodies);
		WriteReloadParticle(In_Output& in_output, SPHBodyVector bodies, StdVec<string> given_body_names);
		virtual ~WriteReloadParticle() {};

		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	  * @class ReadReloadParticle
	  * @brief Write the reload particles file in XML format.
	  */
	class ReadReloadParticle : public ReloadParticleIO, public ReadBodyStates
	{
	public:
		ReadReloadParticle(In_Output& in_output, SPHBodyVector bodies, StdVec<std::string> reload_body_names);
		virtual ~ReadReloadParticle() { };

		virtual void ReadFromFile(size_t iteration_step = 0) override;
	};

	/**
	 * @class RestartIO
	 * @brief Write the restart file in XML format.
	 */
	class RestartIO
	{
	protected:
		std::string overall_file_path_;
		StdVec<std::string> file_paths_;

	public:
		RestartIO(In_Output& in_output, SPHBodyVector bodies);
		virtual ~RestartIO() {};
	};

	/**
	 * @class WriteRestart
	 * @brief Write the restart file in XML format.
	 */
	class WriteRestart : public RestartIO, public WriteBodyStates
	{
	public:
		WriteRestart(In_Output& in_output, SPHBodyVector bodies)
			: RestartIO(in_output, bodies), WriteBodyStates(in_output, bodies) {};
		virtual ~WriteRestart() {};

		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class WriteRestart
	 * @brief Write the restart file in XML format.
	 */
	class ReadRestart : public RestartIO, public ReadBodyStates
	{
	protected:
		Real ReadRestartTime(size_t restart_step);
	public:
		ReadRestart(In_Output& in_output, SPHBodyVector bodies)
			: RestartIO(in_output, bodies), ReadBodyStates(in_output, bodies) {};
		virtual ~ReadRestart() {};
		virtual Real ReadRestartFiles(size_t restart_step) {
			ReadFromFile(restart_step);
			return ReadRestartTime(restart_step);
		};
		virtual void ReadFromFile(size_t iteration_step = 0) override;
	};

	/**
	 * @class WriteSimBodyPinData
	* @brief Write total force acting a solid body.
	*/
	class WriteSimBodyPinData : public WriteSimBodyStates<SimTK::MobilizedBody::Pin>
	{
	protected:
		std::string filefullpath_;
	public:
		WriteSimBodyPinData(In_Output& in_output, SimTK::RungeKuttaMersonIntegrator& integ, SimTK::MobilizedBody::Pin& pinbody);
		virtual ~WriteSimBodyPinData() {};
		virtual void WriteToFile(Real time = 0.0) override;
	};
		/**
	 * @class MaterialPropertyIO
	 * @brief base class for write and read material property.
	 */
	class MaterialPropertyIO
	{
	protected:
		In_Output &in_output_;
		BaseMaterial *material_;
		MaterialVector materials_;
	public:
		MaterialPropertyIO(In_Output &in_output, BaseMaterial *material)
			: in_output_(in_output), material_(material) {};
		MaterialPropertyIO(In_Output &in_output, MaterialVector materials)
			: in_output_(in_output), materials_(materials) {};
		virtual ~MaterialPropertyIO() {};
	};

	/**
	 * @class WriteMaterialProperty
	 * @brief base class for write material property.
	 */
	class WriteMaterialProperty : public MaterialPropertyIO
	{
	public:
		WriteMaterialProperty(In_Output &in_output, BaseMaterial *material)
			: MaterialPropertyIO(in_output, material) {};
		WriteMaterialProperty(In_Output &in_output, MaterialVector materials)
			: MaterialPropertyIO(in_output, materials) {};
		virtual ~WriteMaterialProperty() {};

		virtual void WriteToFile(Real time) = 0;
	};

	/**
	 * @class ReadBodyStates
	 * @brief base class for read material property.
	 */
	class ReadMaterialProperty : public MaterialPropertyIO
	{
	public:
		ReadMaterialProperty(In_Output &in_output, BaseMaterial *material)
			: MaterialPropertyIO(in_output, material) {};
		ReadMaterialProperty(In_Output &in_output, MaterialVector materials)
			: MaterialPropertyIO(in_output, materials) {};
		virtual ~ReadMaterialProperty() {};

		virtual void ReadFromFile(size_t iteration_step) = 0;
	};

	/**
	 * @class ReloadMaterialPropertyIO
	 * @brief For write  and read material property.
	 */
	class ReloadMaterialPropertyIO
	{
	protected:
		std::string file_path_;
	public:
		ReloadMaterialPropertyIO(In_Output& in_output, BaseMaterial *material);
		virtual ~ReloadMaterialPropertyIO() {};
	};

	/**
	  * @class WriteReloadMaterialProperty
	  * @brief Write the material property to file in XML format.
	  */
	class WriteReloadMaterialProperty : public ReloadMaterialPropertyIO, public WriteMaterialProperty
	{
	public:
		WriteReloadMaterialProperty(In_Output& in_output, BaseMaterial *material)
			: ReloadMaterialPropertyIO(in_output, material), WriteMaterialProperty(in_output, material) {};
		WriteReloadMaterialProperty(In_Output& in_output, BaseMaterial* material, string given_material_name)
			:WriteReloadMaterialProperty(in_output, material) {
			file_path_ = in_output.reload_folder_ + "/Material_" + given_material_name + "_rld.xml";
		};

		virtual ~WriteReloadMaterialProperty() {};

		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	  * @class ReadReloadMaterialProperty
	  * @brief Read the material property to file in XML format.
	  */
	class ReadReloadMaterialProperty : public ReloadMaterialPropertyIO, public ReadMaterialProperty
	{
	public:
		ReadReloadMaterialProperty(In_Output& in_output, BaseMaterial *material);
		virtual ~ReadReloadMaterialProperty() {};

		virtual void ReadFromFile(size_t iteration_step = 0) override;
	};
}
