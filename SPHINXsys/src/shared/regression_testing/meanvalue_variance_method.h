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
 * @file regression_testing.h
 * @brief Classes for the comparison between validated and tested results.
 * @author Bo Zhang, Xiangyu Hu
 */

#pragma once
#include "in_output.h"
#include "xml_engine.h"
#include "all_physical_dynamics.h"

namespace SPH
{
	/**
	* @class meanvalue_variance_method.h
	* @brief the regression testing is based on the converged meanvalue and variance.
	*/

	template<class ObserveMethodType>
	class RegressionTesting : public ObserveMethodType
	{
		/* identify the the variable type from the parent class. */
		using VariableType = decltype(ObserveMethodType::type_indicator_);
	
	protected:
		std::string result_filefullpath_;                   /* the path for result. */
		std::string meanvalue_filefullpath_;                /* the path for meanvalue. */
		std::string variance_filefullpath_;                 /* the path for variance. */
		std::string runtimes_filefullpath_;                 /* the path for runtimes. */
		std::string converged;                              /* the tag for result converged. */
		
		XmlEngine result_xml_engine_in_;                    /* xml engine for result input. */
		XmlEngine meanvalue_xml_engine_in_;                 /* xml engine for meanvalue input. */
		XmlEngine variance_xml_engine_in_;                  /* xml engine for variance input. */
		XmlEngine result_xml_engine_out_;                   /* xml engine for result output. */
		XmlEngine meanvalue_xml_engine_out_;                /* xml engine for meanvalue output. */
		XmlEngine variance_xml_engine_out_;                 /* xml engine for variance output. */

		VariableType threshold_mean_, threshold_variance_;  /* the threshold container for mean and variance. */
		StdVec<DataVec<VariableType>> result_;              /* the container of results in all runs. */
		DataVec<VariableType> meanvalue_, meanvalue_new_;   /* the container of meanvalue. */
		DataVec<VariableType> variance_, variance_new_;     /* the container of variance. */

		size_t i_, j_;                                      /* the size of each layer of the result vector. */
		size_t number_of_snapshot_old_;                     /* the snapshot size of existed result. */
		size_t difference_;                                 /* the difference of snapshot between old and new result. */
		size_t number_of_run_;                              /* the times of run. */
		size_t label_for_repeat_;                           /* the int label for stable convergence. */

		/* initialize the threshold of meanvalue and variance. */
		void InitializeThreshold(Real& threshold_mean, Real& threshold_variance);
		void InitializeThreshold(Vecd& threshold_mean, Vecd& threshold_variance);
		void InitializeThreshold(Matd& threshold_mean, Matd& threshold_variance);

		/* the method for calculating the meanvalue. */
		void GetNewMeanValue(DataVec<Real>& current_result, DataVec<Real>& meanvalue, DataVec<Real>& meanvalue_new);
		void GetNewMeanValue(DataVec<Vecd>& current_result, DataVec<Vecd>& meanvalue, DataVec<Vecd>& meanvalue_new);
		void GetNewMeanValue(DataVec<Matd>& current_result, DataVec<Matd>& meanvalue, DataVec<Matd>& meanvalue_new);

		/* the method for calculating the variance. */
		void GetNewVariance(StdVec<DataVec<Real>>& result, DataVec<Real>& meanvalue_new, DataVec<Real>& variance, DataVec<Real>& variance_new);
		void GetNewVariance(StdVec<DataVec<Vecd>>& result, DataVec<Vecd>& meanvalue_new, DataVec<Vecd>& variance, DataVec<Vecd>& variance_new);
		void GetNewVariance(StdVec<DataVec<Matd>>& result, DataVec<Matd>& meanvalue_new, DataVec<Matd>& variance, DataVec<Matd>& variance_new);

		/* the method for writing data to xml memory. */
		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element, const DataVec<Real>& quantity);
		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element, const DataVec<Vecd>& quantity);
		void WriteDataToXmlMemory(XmlEngine xmlengine, SimTK::Xml::Element element, const DataVec<Matd>& quantity);

		/* the method for comparing the meanvalue and variance.  */
		size_t CompareParameter(string par_name, DataVec<Real>& parameter, DataVec<Real>& parameter_new, Real& threshold);
		size_t CompareParameter(string par_name, DataVec<Vecd>& parameter, DataVec<Vecd>& parameter_new, Vecd& threshold);
		size_t CompareParameter(string par_name, DataVec<Matd>& parameter, DataVec<Matd>& parameter_new, Matd& threshold);

		/** the method for testing the new result. */
		size_t TestingNewResult(size_t diff, DataVec<Real>& current_result, DataVec<Real>& meanvalue, DataVec<Real>& variance);
		size_t TestingNewResult(size_t diff, DataVec<Vecd>& current_result, DataVec<Vecd>& meanvalue, DataVec<Vecd>& variance);
		size_t TestingNewResult(size_t diff, DataVec<Matd>& current_result, DataVec<Matd>& meanvalue, DataVec<Matd>& variance);

        /** setup (load .xml file) and correct the number of old and new result. */
		void SettingupAndCorrection();

		/** read the result, meanvalue, variance from the .xml files. */
		void ReadResultFromXml();
		void ReadMeanAndVarianceToXml();

		/** update the meanvalue, variance to a new value. */
		void UpdateMeanValueAndVariance();

		/** write the result, meanvalue, variance to the .xml files. */
		void WriteResultToXml();
		void WriteMeanAndVarianceToXml();

		/** compare the meanvalue, variance if converged. */
		bool CompareMeanValueAndVariance();

		/** testing the new result if the converged within the range. */
		void ResultTesting();

	public:
		template<typename... ConstructorArgs>
		explicit RegressionTesting(ConstructorArgs... constructor_args) :
			ObserveMethodType(constructor_args...),
			result_xml_engine_in_("result_xml_engine_in", "result"),
			meanvalue_xml_engine_in_("meanvalue_xml_engine_in", "meanvalue"),
			variance_xml_engine_in_("variance_xml_engine_in", "variance"),
			result_xml_engine_out_("result_xml_engine_out", "result"),
			meanvalue_xml_engine_out_("meanvalue_xml_engine_out", "meanvalue"),
			variance_xml_engine_out_("variance_xml_engine_out", "variance")
		{
			result_filefullpath_ = this->in_output_.input_folder_ + "/" + this->body_name_
				+ "_" + this->quantity_name_ + "_" + this->in_output_.restart_step_ + "_result.xml";
			meanvalue_filefullpath_ = this->in_output_.input_folder_ + "/" +this-> body_name_
				+ "_" + this->quantity_name_ + "_" + this->in_output_.restart_step_ + "._meanvalue.xml";
			variance_filefullpath_ = this->in_output_.input_folder_ + "/" + this->body_name_
				+ "_" + this->quantity_name_ + "_" + this->in_output_.restart_step_ + "._variance.xml";
			runtimes_filefullpath_ = this->in_output_.input_folder_ + "/" + this->body_name_
				+ "_" + this->quantity_name_ + "_" + this->in_output_.restart_step_ + "_runtimes.dat";

			if (!fs::exists(runtimes_filefullpath_))
			{
				number_of_run_ = 1;
				converged = "false";
				label_for_repeat_ = 0;
			}
			else
			{
				std::ifstream in_file(runtimes_filefullpath_.c_str());
				in_file >> converged;
				in_file >> number_of_run_;
				in_file >> label_for_repeat_;
				in_file.close();
			};
		};
		virtual ~RegressionTesting();

		/* the interface to write date into xml memory. */
		void writeToFile(Real iteration = 0) 
		{ 
			ObserveMethodType::writeToFile();
			this->WriteToXml(iteration);
		}

		/* the interface to write XML memory into XML file. */
		void WriteXmlToXmlFile() { this->observe_xml_engine_.writeToXmlFile(this->filefullpath_input_); };

		/* read local data from XML file. */
		void ReadXmlFromXmlFile() { this->ReadFromXml(); };

		/** the interface for generating the priori converged result. */
		void generateDataBase(VariableType threshold_mean, VariableType threshold_variance)
		{
			WriteXmlToXmlFile();
			ReadXmlFromXmlFile();
			InitializeThreshold(threshold_mean, threshold_variance);
			if (converged == "false")
			{
				SettingupAndCorrection();
				ReadResultFromXml();
				ReadMeanAndVarianceToXml();
				UpdateMeanValueAndVariance();
				WriteResultToXml();
				WriteMeanAndVarianceToXml();
				CompareMeanValueAndVariance();
			}
		};

		/* the interface for testing the new result. */
		void newResultTesting()
		{
			WriteXmlToXmlFile();
			ReadXmlFromXmlFile();
			SettingupAndCorrection();
			ReadMeanAndVarianceToXml();
			ResultTesting();
		};

		/* calculate the meanvalue of a time series result. */
		void GetMeanvalueInTime(size_t iteration = 0);
	};
};