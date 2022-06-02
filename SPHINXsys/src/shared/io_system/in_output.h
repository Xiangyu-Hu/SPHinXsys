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
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 * HU1527/12-1 and HU1527/12-4.												*
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
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
 * @author	Chi Zhang, Shuoguo Zhang, Zhenxi Zhao and Xiangyu Hu
 */

#pragma once
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "all_physical_dynamics.h"
#include "parameterization.h"
#include "xml_engine.h"

#include <sstream>
#include <iomanip>
#include <fstream>
/** Macro for APPLE compilers*/
#ifdef __APPLE__
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;
#else
#include <experimental/filesystem>
namespace fs = std::experimental::filesystem;
#endif

using VtuStringData = std::map<std::string, std::string>;

namespace SPH
{

	/**
	 * @brief preclaimed classes.
	 */
	class BaseLevelSet;

	/**
	 * @class InOutput
	 * @brief The base class which defines folders for output,
	 * restart and particle reload folders.
	 */
	class InOutput
	{
	private:
		UniquePtrKeeper<ParameterizationIO> parameterization_io_ptr_keeper_;

	public:
		explicit InOutput(SPHSystem &sph_system, bool delete_output = true);
		virtual ~InOutput(){};

		SPHSystem &sph_system_;
		std::string input_folder_;
		std::string output_folder_;
		std::string restart_folder_;
		std::string reload_folder_;
		std::string restart_step_;

		ParameterizationIO &defineParameterizationIO();
	};

	/**
	 * @class PltEngine
	 * @brief The base class which defines Tecplot file related operation.
	 */
	class PltEngine
	{
	public:
		PltEngine(){};
		virtual ~PltEngine(){};

		void writeAQuantityHeader(std::ofstream &out_file, const Real &quantity, const std::string &quantity_name);
		void writeAQuantityHeader(std::ofstream &out_file, const Vecd &quantity, const std::string &quantity_name);
		void writeAQuantity(std::ofstream &out_file, const Real &quantity);
		void writeAQuantity(std::ofstream &out_file, const Vecd &quantity);
	};

	/**
	 * @class BodyStatesIO
	 * @brief base class for write and read body states.
	 */
	class BodyStatesIO
	{
	protected:
		InOutput &in_output_;
		SPHBody *body_;
		SPHBodyVector bodies_;

		std::string convertPhysicalTimeToString(Real physical_time);

	public:
		BodyStatesIO(InOutput &in_output, SPHBody &body)
			: in_output_(in_output), body_(&body), bodies_({&body}){};
		BodyStatesIO(InOutput &in_output, SPHBodyVector bodies)
			: in_output_(in_output), body_(bodies[0]), bodies_(bodies){};
		virtual ~BodyStatesIO(){};
	};

	/**
	 * @class BodyStatesRecording
	 * @brief base class for write body states.
	 */
	class BodyStatesRecording : public BodyStatesIO
	{
	public:
		BodyStatesRecording(InOutput &in_output, SPHBody &body)
			: BodyStatesIO(in_output, body){};
		BodyStatesRecording(InOutput &in_output, SPHBodyVector bodies)
			: BodyStatesIO(in_output, bodies){};
		virtual ~BodyStatesRecording(){};

		/** write with filename indicated by physical time */
		void writeToFile()
		{
			writeWithFileName(convertPhysicalTimeToString(GlobalStaticVariables::physical_time_));
		};

		/** write with filename indicated by iteration step */
		virtual void writeToFile(size_t iteration_step);

	protected:
		virtual void writeWithFileName(const std::string &sequence) = 0;
	};

	/**
	 * @class SimBodyStatesIO
	 * @brief base class for write and read SimBody states.
	 */
	template <class MobilizedBodyType>
	class SimBodyStatesIO
	{
	protected:
		InOutput &in_output_;
		SimTK::RungeKuttaMersonIntegrator &integ_;
		MobilizedBodyType &mobody_;

	public:
		SimBodyStatesIO(InOutput &in_output, SimTK::RungeKuttaMersonIntegrator &integ, MobilizedBodyType &mobody)
			: in_output_(in_output), integ_(integ), mobody_(mobody){};
		virtual ~SimBodyStatesIO(){};
	};

	/**
	 * @class WriteSimBodyStates
	 * @brief base class for write SimBody states.
	 */
	template <class MobilizedBodyType>
	class WriteSimBodyStates : public SimBodyStatesIO<MobilizedBodyType>
	{
	public:
		WriteSimBodyStates(InOutput &in_output, SimTK::RungeKuttaMersonIntegrator &integ, MobilizedBodyType &mobody)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, integ, mobody){};
		virtual ~WriteSimBodyStates(){};

		virtual void writeToFile(size_t iteration_step) = 0;
	};

	/**
	 * @class ReadSimBodyStates
	 * @brief base class for read SimBody states.
	 */
	template <class MobilizedBodyType>
	class ReadSimBodyStates : public SimBodyStatesIO<MobilizedBodyType>
	{
	public:
		ReadSimBodyStates(InOutput &in_output, MobilizedBodyType *mobody)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, mobody){};
		ReadSimBodyStates(InOutput &in_output, StdVec<MobilizedBodyType *> mobodies)
			: SimBodyStatesIO<MobilizedBodyType>(in_output, mobodies){};
		virtual ~ReadSimBodyStates(){};

		virtual void readFromFile(size_t iteration_step) = 0;
	};

	/**
	 * @class BodyStatesRecordingToVtp
	 * @brief  Write files for bodies
	 * the output file is VTK XML format can visualized by ParaView the data type vtkPolyData
	 */
	class BodyStatesRecordingToVtp : public BodyStatesRecording
	{
	public:
		BodyStatesRecordingToVtp(InOutput &in_output, SPHBody &body)
			: BodyStatesRecording(in_output, body){};
		BodyStatesRecordingToVtp(InOutput &in_output, SPHBodyVector bodies)
			: BodyStatesRecording(in_output, bodies){};
		virtual ~BodyStatesRecordingToVtp(){};

	protected:
		virtual void writeWithFileName(const std::string &sequence) override;
	};

	/**
	 * @class BodyStatesRecordingToVtpString
	 * @brief  Write strings for bodies
	 * the output is map of strings with VTK XML format can visualized by ParaView
	 * the data type vtkUnstructedGrid
	 */
	class BodyStatesRecordingToVtpString : public BodyStatesRecording
	{
	public:
		BodyStatesRecordingToVtpString(InOutput &in_output, SPHBodyVector bodies)
			: BodyStatesRecording(in_output, bodies){};
		virtual ~BodyStatesRecordingToVtpString() = default;

		const VtuStringData &GetVtuData() const;
		void clear()
		{
			_vtuData.clear();
		}

	protected:
		virtual void writeWithFileName(const std::string &sequence) override;
		virtual void writeVtu(std::ostream &stream, SPHBody *body) const;

	private:
		VtuStringData _vtuData;
	};

	/**
	 * @class BodyStatesRecordingToPlt
	 * @brief  Write files for bodies
	 * the output file is dat format can visualized by TecPlot
	 */
	class BodyStatesRecordingToPlt : public BodyStatesRecording
	{
	public:
		BodyStatesRecordingToPlt(InOutput &in_output, SPHBody &body)
			: BodyStatesRecording(in_output, body){};
		BodyStatesRecordingToPlt(InOutput &in_output, SPHBodyVector bodies)
			: BodyStatesRecording(in_output, bodies){};
		virtual ~BodyStatesRecordingToPlt(){};

	protected:
		virtual void writeWithFileName(const std::string &sequence) override;
	};

	/**
	 * @class WriteToVtpIfVelocityOutOfBound
	 * @brief  output body sates if particle velocity is
	 * out of a bound
	 */
	class WriteToVtpIfVelocityOutOfBound
		: public BodyStatesRecordingToVtp
	{
	protected:
		bool out_of_bound_;
		StdVec<VelocityBoundCheck> check_bodies_;
		virtual void writeWithFileName(const std::string &sequence) override;

	public:
		WriteToVtpIfVelocityOutOfBound(InOutput &in_output,
									   SPHBodyVector bodies, Real velocity_bound);
		virtual ~WriteToVtpIfVelocityOutOfBound(){};
	};

	/**
	 * @class MeshRecordingToPlt
	 * @brief  write the background mesh data for relax body
	 */
	class MeshRecordingToPlt : public BodyStatesRecording
	{
	protected:
		std::string filefullpath_;
		BaseMeshField *mesh_field_;
		virtual void writeWithFileName(const std::string &sequence) override;

	public:
		MeshRecordingToPlt(InOutput &in_output, SPHBody &body, BaseMeshField *mesh_field);
		virtual ~MeshRecordingToPlt(){};
	};

	/**
	 * @class ObservedQuantityRecording
	 * @brief write files for observed quantity
	 */
	template <typename VariableType>
	class ObservedQuantityRecording : public BodyStatesRecording,
									  public observer_dynamics::ObservingAQuantity<VariableType>
	{
	protected:
		SPHBody *observer_;
		PltEngine plt_engine_;
		BaseParticles *base_particles_;
		std::string body_name_;
		const std::string quantity_name_;
		std::string filefullpath_output_;

	public:
		VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

	public:
		ObservedQuantityRecording(const std::string &quantity_name, InOutput &in_output,
								  BaseBodyRelationContact &contact_relation)
			: BodyStatesRecording(in_output, *contact_relation.sph_body_),
			  observer_dynamics::ObservingAQuantity<VariableType>(contact_relation, quantity_name),
			  observer_(contact_relation.sph_body_), plt_engine_(),
			  base_particles_(observer_->base_particles_), body_name_(contact_relation.sph_body_->getBodyName()),
			  quantity_name_(quantity_name)
		{
			/** Output for .dat file. */
			filefullpath_output_ = in_output_.output_folder_ + "/" + body_name_ + "_" + quantity_name + "_" + in_output_.restart_step_ + ".dat";
			std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
			out_file << "run_time"
					 << "   ";
			for (size_t i = 0; i != base_particles_->total_real_particles_; ++i)
			{
				std::string quantity_name_i = quantity_name + "[" + std::to_string(i) + "]";
				plt_engine_.writeAQuantityHeader(out_file, (*this->interpolated_quantities_)[i], quantity_name_i);
			}
			out_file << "\n";
			out_file.close();
		};
		virtual ~ObservedQuantityRecording(){};

		virtual void writeWithFileName(const std::string &sequence) override
		{
			this->parallel_exec();
			std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
			out_file << GlobalStaticVariables::physical_time_ << "   ";
			for (size_t i = 0; i != base_particles_->total_real_particles_; ++i)
			{
				plt_engine_.writeAQuantity(out_file, (*this->interpolated_quantities_)[i]);
			}
			out_file << "\n";
			out_file.close();
		};
	};

	/**
	 * @class BodyReducedQuantityRecording
	 * @brief write reduced quantity of a body
	 */
	template <class ReduceMethodType>
	class BodyReducedQuantityRecording
	{
	protected:
		InOutput &in_output_;
		PltEngine plt_engine_;
		ReduceMethodType reduce_method_;
		std::string body_name_;
		const std::string quantity_name_;
		std::string filefullpath_output_;

	public:
		/*< deduce variable type from reduce method. */
		using VariableType = decltype(reduce_method_.InitialReference());
		VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

	public:
		template <typename... ConstructorArgs>
		BodyReducedQuantityRecording(InOutput &in_output, ConstructorArgs &&...args)
			: in_output_(in_output), plt_engine_(), reduce_method_(std::forward<ConstructorArgs>(args)...),
			  body_name_(reduce_method_.getSPHBody()->getBodyName()),
			  quantity_name_(reduce_method_.QuantityName())
		{
			/** output for .dat file. */
			filefullpath_output_ = in_output_.output_folder_ + "/" + body_name_ + "_" + quantity_name_ + "_" + in_output_.restart_step_ + ".dat";
			std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
			out_file << "\"run_time\""
					 << "   ";
			plt_engine_.writeAQuantityHeader(out_file, reduce_method_.InitialReference(), quantity_name_);
			out_file << "\n";
			out_file.close();
		};
		virtual ~BodyReducedQuantityRecording(){};

		virtual void writeToFile(size_t iteration_step = 0)
		{
			std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
			out_file << GlobalStaticVariables::physical_time_ << "   ";
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
		ReloadParticleIO(InOutput &in_output, SPHBodyVector bodies);
		ReloadParticleIO(InOutput &in_output, SPHBodyVector bodies, const StdVec<std::string> &given_body_names);
		virtual ~ReloadParticleIO(){};

		virtual void writeToFile(size_t iteration_step = 0);
		virtual void readFromFile(size_t iteration_step = 0);
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
		RestartIO(InOutput &in_output, SPHBodyVector bodies);
		virtual ~RestartIO(){};

		virtual void writeToFile(size_t iteration_step = 0);
		virtual void readFromFile(size_t iteration_step = 0);
		virtual Real readRestartFiles(size_t restart_step)
		{
			readFromFile(restart_step);
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
		WriteSimBodyPinData(InOutput &in_output, SimTK::RungeKuttaMersonIntegrator &integ, SimTK::MobilizedBody::Pin &pinbody);
		virtual ~WriteSimBodyPinData(){};
		virtual void writeToFile(size_t iteration_step = 0) override;
	};

	/**
	 * @class ReloadMaterialParameterIO
	 * @brief For write  and read material property.
	 */
	class ReloadMaterialParameterIO
	{
	protected:
		InOutput &in_output_;
		BaseMaterial *base_material_;
		std::string file_path_;

	public:
		ReloadMaterialParameterIO(InOutput &in_output, BaseMaterial *base_material);
		ReloadMaterialParameterIO(InOutput &in_output, BaseMaterial *base_material, const std::string &given_parameters_name);
		virtual ~ReloadMaterialParameterIO(){};

		virtual void writeToFile(size_t iteration_step = 0);
		virtual void readFromFile(size_t iteration_step = 0);
	};
}
