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

#pragma once
#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "all_physical_dynamics.h"
#include "xml_engine.h" 
 
#include "SimTKcommon.h"
#include "SimTKmath.h"
#include "Simbody.h"

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
		std::string input_folder_;
		std::string output_folder_;
		std::string restart_folder_;
		std::string reload_folder_;
			
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

		std::string convertPhysicalTimeToString(Real physical_time);
	public:
		BodyStatesIO(In_Output &in_output, SPHBody *body)
			: in_output_(in_output), body_(body), bodies_({body}) {};
		BodyStatesIO(In_Output& in_output, SPHBodyVector bodies)
			: in_output_(in_output), body_(bodies[0]), bodies_(bodies) {};
		virtual ~BodyStatesIO() {};
	};

	/**
	 * @class BodyStatesRecording
	 * @brief base class for write body states.
	 */
	class BodyStatesRecording : public BodyStatesIO
	{
	public:
		BodyStatesRecording(In_Output &in_output, SPHBody *body)
			: BodyStatesIO(in_output, body) {};
		BodyStatesRecording(In_Output &in_output, SPHBodyVector bodies)
			: BodyStatesIO(in_output, bodies) {};
		virtual ~BodyStatesRecording() {};

		/** write with filename indicated by physical time */
		void writeToFile() 
		{
			writeWithFileName(convertPhysicalTimeToString(GlobalStaticVariables::physical_time_));
		};

		/** write with filename indicated by iteration step */
		virtual void writeToFile(size_t iteration_step)
		{
			writeWithFileName(std::to_string(iteration_step));
		};

	protected:
		virtual void writeWithFileName(const std::string& sequence) = 0;
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

		virtual void writeToFile(size_t iteration_step) = 0;
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

		virtual void readFromFile(size_t iteration_step) = 0;
	};

	/**
	 * @class BodyStatesRecordingToVtu
	 * @brief  Write files for bodies
	 * the output file is VTK XML format can visualized by ParaView
	 * the data type vtkUnstructedGrid
	 */
	class BodyStatesRecordingToVtu : public BodyStatesRecording
	{
	public:
		BodyStatesRecordingToVtu(In_Output& in_output, SPHBodyVector bodies)
			: BodyStatesRecording(in_output, bodies) {};
		virtual ~BodyStatesRecordingToVtu() {};

	protected:
		virtual void writeWithFileName(const std::string& sequence) override;
		virtual void writeVtu(std::ostream& stream, SPHBody* body) const;
	};

	/**
	 * @class BodyStatesRecordingToVtuString
	 * @brief  Write strings for bodies
	 * the output is map of strings with VTK XML format can visualized by ParaView
	 * the data type vtkUnstructedGrid
	 */
	class BodyStatesRecordingToVtuString : public BodyStatesRecordingToVtu
	{
	public:
		BodyStatesRecordingToVtuString(In_Output& in_output, SPHBodyVector bodies);
		virtual ~BodyStatesRecordingToVtuString() = default;

		using VtuStringData = std::map<std::string, std::string>;

		const VtuStringData& GetVtuData() const;
	protected:
		virtual void writeWithFileName(const std::string& sequence) override;
	private:
		VtuStringData _vtuData;
	};

	/**
	 * @class BodyStatesRecordingToVtuStringRunTime
	 * @brief  Writes simulation results as strings for bodies in rum time
	 * the output is a std::pair of strings with VTK XML format that can be visualized by ParaView
	 * the data type vtkUnstructedGrid
	 */
	class BodyStatesRecordingToVtuStringRunTime : public BodyStatesRecordingToVtu
	{
	public:
		BodyStatesRecordingToVtuStringRunTime(In_Output& in_output, SPHBodyVector bodies);
		virtual ~BodyStatesRecordingToVtuStringRunTime() = default;

		using VtuStringDataRunTime = std::pair<std::string, std::string>;

		const VtuStringDataRunTime& GetVtuDataRunTime() const;
	protected:
		void writeWithFileName(const std::string& sequence) override;
	private:
		VtuStringDataRunTime _vtuDataRunTime;
	};

	/**
	 * @class SurfaceOnlyBodyStatesRecordingToVtu
	 * @brief  Write files for surface particles of bodies
	 * the output file is VTK XML format can visualized by ParaView
	 * the data type vtkUnstructedGrid
	 */
	class SurfaceOnlyBodyStatesRecordingToVtu : public BodyStatesRecording
	{
	public:
		SurfaceOnlyBodyStatesRecordingToVtu(In_Output& in_output, SPHBodyVector bodies);

	protected:
		virtual void writeWithFileName(const std::string& sequence) override;
		StdVec<ShapeSurface> surface_body_layer_vector_;
	};
	
	/**
	 * @class BodyStatesRecordingToPlt
	 * @brief  Write files for bodies
	 * the output file is dat format can visualized by TecPlot
	 */
	class BodyStatesRecordingToPlt : public BodyStatesRecording
	{
	public:
		BodyStatesRecordingToPlt(In_Output& in_output, SPHBodyVector bodies)
			: BodyStatesRecording(in_output, bodies) {};
		virtual ~BodyStatesRecordingToPlt() {};

	protected:
		virtual void writeWithFileName(const std::string& sequence) override;
	};

	/**
	 * @class WriteToVtuIfVelocityOutOfBound
	 * @brief  output body sates if particle velocity is
	 * out of a bound
	 */
	class WriteToVtuIfVelocityOutOfBound
		: public BodyStatesRecordingToVtu
	{
	protected:
		bool out_of_bound_;
		StdVec<VelocityBoundCheck *> check_bodies_;
		virtual void writeWithFileName(const std::string& sequence) override;
	public:
		WriteToVtuIfVelocityOutOfBound(In_Output& in_output,
			SPHBodyVector bodies, Real velocity_bound);
		virtual ~WriteToVtuIfVelocityOutOfBound() {};
	};

	/**
	 * @class MeshRecordingToPlt
	 * @brief  write the background mesh data for relax body
	 */
	class MeshRecordingToPlt : public BodyStatesRecording
	{
	protected:
		std::string filefullpath_;
		BaseMeshField* mesh_field_;
		virtual void writeWithFileName(const std::string& sequence) override;
	public:
		MeshRecordingToPlt(In_Output& in_output, SPHBody* body, BaseMeshField* mesh_field);
		virtual ~MeshRecordingToPlt() {};
	};

	/**
	 * @class ObservedQuantityRecording
	 * @brief write files for observed quantity
	 */
	template <int DataTypeIndex, typename VariableType>
	class ObservedQuantityRecording : public BodyStatesRecording,
		public observer_dynamics::InterpolatingAQuantity<DataTypeIndex, VariableType>
	{
	protected:
		SPHBody* observer_;
		PltEngine plt_engine_;
		BaseParticles* base_particles_;
		std::string body_name_;
		std::string quantity_name_;
		XmlEngine observe_xml_engine_;
		std::string filefullpath_input_;
		std::string filefullpath_output_;

		DataVec<VariableType> current_result_; /* the container of the current result. */
		StdVec<string> element_tag_;           /* the container of the current tag. */

		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element,
			string element_name, size_t particle_n, const Real& quantity)
		{
			SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
			std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
			xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity);
		};
		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element,
			string element_name, size_t particle_n, const Vecd& quantity)
		{
			SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
			std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
			xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity);
		};
		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element,
			string element_name, size_t particle_n, const Matd& quantity)
		{
			SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
			std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
			xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity);
		};
		
		void ReadDataFromXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element,
			size_t particle_n, DataVec<Real>& container)
		{
			size_t index_i_ = 0;
			SimTK::Xml::element_iterator ele_ite = element.element_begin();
			for (; ele_ite != element.element_end(); ++ele_ite)
			{
				std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
				xmlengine.getRequiredAttributeValue<Real>(ele_ite, attribute_name_, container[index_i_][particle_n]);
				index_i_++;
			}
		};
		void ReadDataFromXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element, 
			size_t particle_n, DataVec<Vecd>& container)
		{
			size_t index_i_ = 0;
			SimTK::Xml::element_iterator ele_ite = element.element_begin();
			for (; ele_ite != element.element_end(); ++ele_ite)
			{
				std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
				xmlengine.getRequiredAttributeValue<Vecd>(ele_ite, attribute_name_, container[index_i_][particle_n]);
				index_i_++;
			}
		};
		void ReadDataFromXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element, 
			size_t particle_n, DataVec<Matd>& container)
		{
			size_t index_i_ = 0;
			SimTK::Xml::element_iterator ele_ite = element.element_begin();
			for (; ele_ite != element.element_end(); ++ele_ite)
			{
				std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
				xmlengine.getRequiredAttributeMatrixValue(ele_ite, attribute_name_, container[index_i_][particle_n]);
				index_i_++;
			}
		};

		void ReadTagFromXmlMemory(SimTK::Xml::Element element, StdVec<string>& element_tag)
		{
			size_t index_i_ = 0;
			SimTK::Xml::element_iterator ele_ite = element.element_begin();
			for (; ele_ite != element.element_end(); ++ele_ite)
			{
				element_tag[index_i_] = ele_ite->getElementTag();
				index_i_++;
			}
		};
	
	public:
		ObservedQuantityRecording(std::string quantity_name, In_Output& in_output,
			BaseBodyRelationContact* body_contact_relation) :
			BodyStatesRecording(in_output, body_contact_relation->sph_body_),
			observer_dynamics::InterpolatingAQuantity<DataTypeIndex, VariableType>(body_contact_relation, quantity_name),
			observer_(body_contact_relation->sph_body_), plt_engine_(), base_particles_(observer_->base_particles_),
			body_name_(body_contact_relation->sph_body_->getBodyName()), quantity_name_(quantity_name),
			observe_xml_engine_("xml_observe", quantity_name_)
		{
			/** Output for .dat file. */
			filefullpath_output_ = in_output_.output_folder_ + "/" + body_name_
				+ "_" + quantity_name + "_" + in_output_.restart_step_ + ".dat";
			std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
			out_file << "run_time" << "   ";
			for (size_t i = 0; i != base_particles_->total_real_particles_; ++i)
			{
				std::string quantity_name_i = quantity_name + "[" + std::to_string(i) + "]";
				plt_engine_.writeAQuantityHeader(out_file, this->interpolated_quantities_[i], quantity_name_i);
			}
			out_file << "\n";
			out_file.close();

			/** Output for .xml file. */
			filefullpath_input_ = in_output_.input_folder_ + "/" + body_name_
				+ "_" + quantity_name + "_" + in_output_.restart_step_ + ".xml";
		};
		virtual ~ObservedQuantityRecording() {};

		VariableType type_indicator_; /*< this is an indicator to identify the variable type. */
	
		void WriteXmlToXmlFile() { observe_xml_engine_.writeToXmlFile(filefullpath_input_); }

		virtual void writeWithFileName(const std::string& sequence) override
		{
			this->parallel_exec();
			std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
			out_file << GlobalStaticVariables::physical_time_ << "   ";
			for (size_t i = 0; i != base_particles_->total_real_particles_; ++i)
			{
				plt_engine_.writeAQuantity(out_file, this->interpolated_quantities_[i]);
			}
			out_file << "\n";
			out_file.close();
		};

		void WriteToXml(size_t iteration = 0)
		{
			this->parallel_exec();
			std::string element_name_ = "Snapshot_" + std::to_string(iteration);
			SimTK::Xml::Element element_ = observe_xml_engine_.root_element_;
			observe_xml_engine_.addElementToXmlDoc(element_name_);
			for (size_t i = 0; i != base_particles_->total_real_particles_; ++i)
			{
				WriteDataToXmlMemory(observe_xml_engine_, element_, element_name_, i, this->interpolated_quantities_[i]);
			};
		};

		void ReadFromXml()
		{
			observe_xml_engine_.loadXmlFile(filefullpath_input_);
			size_t number_of_particle_ = base_particles_->total_real_particles_;
			size_t number_of_snapshot_ = std::distance(observe_xml_engine_.root_element_.element_begin(),
				observe_xml_engine_.root_element_.element_end());
			DataVec<VariableType> current_result_temp_(number_of_snapshot_, StdVec<VariableType>(number_of_particle_));
			StdVec<string> element_tag_temp_(number_of_snapshot_);
			current_result_ = current_result_temp_;
			element_tag_ = element_tag_temp_;
			SimTK::Xml::Element element_ = observe_xml_engine_.root_element_;
			for (size_t j = 0; j != number_of_particle_; ++j)
			{
				ReadDataFromXmlMemory(observe_xml_engine_, element_, j, current_result_);
				ReadTagFromXmlMemory(element_, element_tag_);
			}
		};
	};

	/**
	 * @class BodyReducedQuantityRecording
	 * @brief write reduced quantity of a body
	 */
	template<class ReduceMethodType>
	class BodyReducedQuantityRecording
	{
	protected:
		In_Output& in_output_;
		PltEngine plt_engine_;
		ReduceMethodType reduce_method_;
		std::string body_name_;
		std::string quantity_name_;
		XmlEngine observe_xml_engine_;
		std::string filefullpath_input_;
		std::string filefullpath_output_;

		/*< deduce variable type from reduce method. */
		using VariableType = decltype(reduce_method_.InitialReference());
		DataVec<VariableType> current_result_; /* the container of the current result. */
		StdVec<string> element_tag_;           /* the container of the current tag. */

		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element,
			string element_name, size_t particle_n, const Real& quantity)
		{
			SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
			std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
			xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity);
		};
		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element,
			string element_name, size_t particle_n, const Vecd& quantity)
		{
			SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
			std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
			xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity);
		};
		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element,
			string element_name, size_t particle_n, const Matd& quantity)
		{
			SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
			std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
			xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity);
		};

		void ReadDataFromXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element_name, 
			size_t particle_n, DataVec<Real>& container)
		{
			size_t index_i_ = 0;
			SimTK::Xml::element_iterator ele_ite = element_name.element_begin();
			for (; ele_ite != element_name.element_end(); ++ele_ite)
			{
				std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
				xmlengine.getRequiredAttributeValue<Real>(ele_ite, attribute_name_, container[index_i_][particle_n]);
				index_i_++;
			}
		};
		void ReadDataFromXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element_name, 
			size_t particle_n, DataVec<Vecd>& container)
		{
			size_t index_i_ = 0;
			SimTK::Xml::element_iterator ele_ite = element_name.element_begin();
			for (; ele_ite != element_name.element_end(); ++ele_ite)
			{
				std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
				xmlengine.getRequiredAttributeValue<Vecd>(ele_ite, attribute_name_, container[index_i_][particle_n]);
				index_i_++;
			}
		};
		void ReadDataFromXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element_name, 
			size_t particle_n, DataVec<Matd>& container)
		{
			size_t index_i_ = 0;
			SimTK::Xml::element_iterator ele_ite = element_name.element_begin();
			for (; ele_ite != element_name.element_end(); ++ele_ite)
			{
				std::string attribute_name_ = quantity_name_ + "_" + std::to_string(particle_n);
				xmlengine.getRequiredAttributeMatrixValue(ele_ite, attribute_name_, container[index_i_][particle_n]);
				index_i_++;
			}
		};

		void ReadTagFromXmlMemory(SimTK::Xml::Element element, StdVec<string>& element_tag)
		{
			size_t index_i_ = 0;
			SimTK::Xml::element_iterator ele_ite = element.element_begin();
			for (; ele_ite != element.element_end(); ++ele_ite)
			{
				element_tag[index_i_] = ele_ite->getElementTag();
				index_i_++;
			}
		};

	public:
		template<typename... ConstructorArgs>
		BodyReducedQuantityRecording(In_Output& in_output, ConstructorArgs... constructor_args) :
			in_output_(in_output), plt_engine_(), reduce_method_(constructor_args...),
			body_name_(reduce_method_.getSPHBody()->getBodyName()),
			quantity_name_(reduce_method_.QuantityName()),
			observe_xml_engine_("xml_reduce", quantity_name_)
		{
			/** output for .dat file. */
			filefullpath_output_ = in_output_.output_folder_ + "/" + body_name_
				+ "_" + quantity_name_ + "_" + in_output_.restart_step_ + ".dat";
			std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
			out_file << "\"run_time\"" << "   ";
			plt_engine_.writeAQuantityHeader(out_file, reduce_method_.InitialReference(), quantity_name_);
			out_file << "\n";
			out_file.close();

			/** output for .xml file. */
			filefullpath_input_ = in_output_.input_folder_ + "/" + body_name_
				+ "_" + quantity_name_ + "_" + in_output_.restart_step_ + ".xml";
		};
		virtual ~BodyReducedQuantityRecording() {};

		VariableType type_indicator_; /*< this is an indicator to identify the variable type. */

		void writeXmlToXmlFile() { observe_xml_engine_.writeToXmlFile(filefullpath_input_); }

		virtual void writeToFile(size_t iteration_step = 0)
		{
			std::ofstream out_file(filefullpath_output_.c_str(), std::ios::app);
			out_file << GlobalStaticVariables::physical_time_ << "   ";
			plt_engine_.writeAQuantity(out_file, reduce_method_.parallel_exec());
			out_file << "\n";
			out_file.close();
		};

		void WriteToXml(size_t iteration = 0)
		{
			std::string element_name_ = "Snapshot_" + std::to_string(iteration);
			SimTK::Xml::Element element_ = observe_xml_engine_.root_element_;
			observe_xml_engine_.addElementToXmlDoc(element_name_);
			WriteDataToXmlMemory(observe_xml_engine_, element_, element_name_, 0,  reduce_method_.parallel_exec());
		};

		void ReadFromXml()
		{
			observe_xml_engine_.loadXmlFile(filefullpath_input_);
			size_t number_of_particle_ = 1;
			size_t number_of_snapshot_ = std::distance(observe_xml_engine_.root_element_.element_begin(),
				observe_xml_engine_.root_element_.element_end());
			DataVec<VariableType> current_result_temp_(number_of_snapshot_, StdVec<VariableType>(number_of_particle_));
			StdVec<string> element_tag_temp_(number_of_snapshot_);
			current_result_ = current_result_temp_;
			element_tag_ = element_tag_temp_;
			SimTK::Xml::Element element_ = observe_xml_engine_.root_element_;
			for (size_t j = 0; j != number_of_particle_; ++j)
			{
				ReadDataFromXmlMemory(observe_xml_engine_, element_, j, current_result_);
				ReadTagFromXmlMemory(element_, element_tag_);
			}
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
		RestartIO(In_Output& in_output, SPHBodyVector bodies);
		virtual ~RestartIO() {};

		virtual void writeToFile(size_t iteration_step = 0);
		virtual void readFromFile(size_t iteration_step = 0);
		virtual Real readRestartFiles(size_t restart_step) {
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
		WriteSimBodyPinData(In_Output& in_output, SimTK::RungeKuttaMersonIntegrator& integ, SimTK::MobilizedBody::Pin& pinbody);
		virtual ~WriteSimBodyPinData() {};
		virtual void writeToFile(size_t iteration_step = 0) override;
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

		virtual void writeToFile(size_t iteration_step = 0);
		virtual void readFromFile(size_t iteration_step = 0);
	};
}