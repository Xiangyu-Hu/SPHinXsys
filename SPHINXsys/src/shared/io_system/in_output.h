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
 */


#ifndef IN_OUTPUT_H
#define IN_OUTPUT_H


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
	 * @brief preclaimed classes.
	 */
	class BaseLevelSet;

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
		std::string input_folder_;			/**< Folder for saving input files. */
		std::string output_folder_;			/**< Folder for saving output files. */
		std::string restart_folder_;		/**< Folder for saving restart files. */
		std::string reload_folder_;			/**< Folder for saving particle reload files. */

		std::string restart_step_;
	};

	/**
	 * @class PltEngine
	 * @brief The base class which defines Tecplot file related operation.
	 */
	class PltEngine
	{
	public:
		PltEngine() {};
		virtual ~PltEngine() {};

		void writeAQuantityHeader(std::ofstream& out_file, const Real& quantity, std::string quantity_name);
		void writeAQuantityHeader(std::ofstream& out_file, const Vecd& quantity, std::string quantity_name);
		void writeAQuantity(std::ofstream& out_file, const Real& quantity);
		void writeAQuantity(std::ofstream& out_file, const Vecd& quantity);
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
			: in_output_(in_output), body_(body), bodies_({body}) {};
		BodyStatesIO(In_Output& in_output, SPHBodyVector bodies)
			: in_output_(in_output), body_(bodies[0]), bodies_(bodies) {};
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
			: WriteBodyStates(in_output, bodies) {};
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
	 * @class WriteMeshToPlt
	 * @brief  write the background mesh data for relax body
	 */
	class WriteMeshToPlt : public WriteBodyStates
	{
	protected:
		std::string filefullpath_;
		Mesh* mesh_;
	public:
		WriteMeshToPlt(In_Output& in_output, SPHBody* body, Mesh* mesh);
		virtual ~WriteMeshToPlt() {};

		virtual void WriteToFile(Real time = 0.0) override;
	};

	/**
	 * @class WriteAnObservedQuantity
	 * @brief write files for observed quantity
	 */
	template <int DataTypeIndex, typename VariableType>
	class WriteAnObservedQuantity : public WriteBodyStates,
		public observer_dynamics::InterpolatingAQuantity<DataTypeIndex, VariableType>
	{
	protected:
		SPHBody* observer_;
		PltEngine plt_engine_;
		BaseParticles* base_particles_;
		std::string filefullpath_;

	public:
		WriteAnObservedQuantity(std::string quantity_name, In_Output& in_output,
			BaseContactBodyRelation* body_contact_relation) : 
			WriteBodyStates(in_output, body_contact_relation->sph_body_), 
			observer_dynamics::InterpolatingAQuantity<DataTypeIndex, VariableType>(body_contact_relation, quantity_name),
			observer_(body_contact_relation->sph_body_), plt_engine_(), base_particles_(observer_->base_particles_)
		{
			filefullpath_ = in_output_.output_folder_ + "/" + observer_->getBodyName()
				+ "_" + quantity_name + "_" + in_output_.restart_step_ + ".dat";
			std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
			out_file << "run_time" << "   ";
			for (size_t i = 0; i != base_particles_->total_real_particles_; ++i)
			{
				std::string quantity_name_i = quantity_name + "[" + std::to_string(i) + "]";
				plt_engine_.writeAQuantityHeader(out_file, this->interpolated_quantities_[i], quantity_name_i);
			}
			out_file << "\n";
			out_file.close();
		};
		virtual ~WriteAnObservedQuantity() {};

		virtual void WriteToFile(Real time = 0.0) override 
		{
			this->parallel_exec();
			std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
			out_file << time << "   ";
			for (size_t i = 0; i != base_particles_->total_real_particles_; ++i)
			{
				plt_engine_.writeAQuantity(out_file, this->interpolated_quantities_[i]);
			}
			out_file << "\n";
			out_file.close();
		};
	};

	/**
	 * @class WriteBodyReducedQuantity
	 * @brief write reduced quantity of a body
	 */
	template<class ReduceMethodType>
	class WriteBodyReducedQuantity 
	{
	protected:
		In_Output& in_output_; 
		PltEngine plt_engine_;
		ReduceMethodType reduce_method_;
		std::string body_name_;
		std::string quantity_name_;
		std::string filefullpath_;

	public:
		template<typename... ConstructorArgs>
		WriteBodyReducedQuantity(In_Output& in_output, ConstructorArgs... constructor_args) :
			in_output_(in_output), plt_engine_(), reduce_method_(constructor_args...),
			body_name_(reduce_method_.getSPHBody()->getBodyName()),
			quantity_name_(reduce_method_.QuantityName())
		{
			filefullpath_ = in_output_.output_folder_ + "/" + body_name_
				+ "_" + quantity_name_ + "_" + in_output_.restart_step_ + ".dat";
			std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
			out_file << "\"run_time\"" << "   ";
			plt_engine_.writeAQuantityHeader(out_file, reduce_method_.InitialReference(), quantity_name_);
			out_file << "\n";
			out_file.close();
		};
		virtual ~WriteBodyReducedQuantity() {};

		virtual void WriteToFile(Real time = 0.0)
		{
			std::ofstream out_file(filefullpath_.c_str(), std::ios::app);
			out_file << time << "   ";
			plt_engine_.writeAQuantity(out_file, reduce_method_.parallel_exec());
			out_file << "\n";
			out_file.close();
		};
	};

	/**
	  * @class ReloadParticleIO
	  * @brief Write the reload particles file in XML format.
	  */
	class ReloadParticleIO : public BodyStatesIO
	{
	protected:
		StdVec<std::string> file_paths_;
	public:
		ReloadParticleIO(In_Output& in_output, SPHBodyVector bodies);
		ReloadParticleIO(In_Output& in_output, SPHBodyVector bodies, StdVec<std::string> given_body_names);
		virtual ~ReloadParticleIO() {};

		virtual void WriteToFile(Real time = 0.0);
		virtual void ReadFromFile(size_t iteration_step = 0);
	};

	/**
	 * @class RestartIO
	 * @brief Write the restart file in XML format.
	 */
	class RestartIO : public BodyStatesIO
	{
	protected:
		std::string overall_file_path_;
		StdVec<std::string> file_paths_;

		Real readRestartTime(size_t restart_step);
	public:
		RestartIO(In_Output& in_output, SPHBodyVector bodies);
		virtual ~RestartIO() {};

		virtual void WriteToFile(Real time = 0.0);
		virtual void ReadFromFile(size_t iteration_step = 0);
		virtual Real readRestartFiles(size_t restart_step) {
			ReadFromFile(restart_step);
			return readRestartTime(restart_step);
		};
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
	 * @class ReloadMaterialParameterIO
	 * @brief For write  and read material property.
	 */
	class ReloadMaterialParameterIO
	{
	protected:
		In_Output& in_output_;		
		BaseMaterial *material_;
		std::string file_path_;
	public:
		ReloadMaterialParameterIO(In_Output& in_output, BaseMaterial* material);
		ReloadMaterialParameterIO(In_Output& in_output, BaseMaterial *material, std::string given_parameters_name);
		virtual ~ReloadMaterialParameterIO() {};

		virtual void WriteToFile(Real time = 0.0);
		virtual void ReadFromFile(size_t iteration_step = 0);
	};
}
#endif //IN_OUTPUT_H