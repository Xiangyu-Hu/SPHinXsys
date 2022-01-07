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
 * @file regression_test_base.h
 * @brief Base classes for the comparison between validated and tested results.
 * @author Bo Zhang and Xiangyu Hu
 */

#pragma once
#include "in_output.h"
#include "xml_engine.h"
#include "all_physical_dynamics.h"

namespace SPH
{
	/**
	 * @class RegressionTestBase
	 * @brief The base of regression test for various method (time-averaged, ensemble-averaged, dynamic time warping)
	 */
	 template<class ObserveMethodType>
	 class RegressionTestBase : public ObserveMethodType
	 {
		 /*identify the variable type from the parent class. */
		 using VariableType = decltype(ObserveMethodType::type_indicator_);

	 protected:
		 std::string input_folder_path_;             /* the path for input folder. */
		 std::string result_filefullpath_;           /* the path for all results. */
		 std::string runtimes_filefullpath_;         /* the path for run times. */
		 std::string converged;                      /* the tag for result converged. */

		 XmlMemoryIO xmlmemory_io_;                  /* xml memory in_output operator */
		 XmlEngine result_xml_engine_in_;            /* xml engine for result input. */
		 XmlEngine result_xml_engine_out_;           /* xml engine for result output. */

		 DoubleVec<VariableType> current_result_ji_; /* the container of current result with i&j transfered. */
		 DoubleVec<VariableType> result_in_;         /* the container of each result with i&j transfered. */
		 TripleVec<VariableType> result_;            /* the container of results in all runs.*/
		 
		 int i_, j_;                                 /* the size of each layer of current result. */
		 int difference_;                            /* the difference of snapshot between old and new result. */
		 int number_of_run_;                         /* the times of run. */
		 int label_for_repeat_;                      /* the label for stable convergence. */
		 int number_of_snapshot_old_;                /* the snapshot size of existed result. */
		 
	 public:
		 template<typename... ConstructorArgs>
		 explicit RegressionTestBase(ConstructorArgs &&...args) :
			 ObserveMethodType(std::forward<ConstructorArgs>(args)...), xmlmemory_io_(),
			 result_xml_engine_in_("result_xml_engine_in", "result"),
			 result_xml_engine_out_("result_xml_engine_out", "result")
		 {
			 input_folder_path_ = this->in_output_.input_folder_;
			 result_filefullpath_ = input_folder_path_ + "/" + this->body_name_
				 + "_" + this->quantity_name_ + "_result.xml";
			 runtimes_filefullpath_ = input_folder_path_ + "/" + this->body_name_
				 + "_" + this->quantity_name_ + "_runtimes.dat";

			 if (!fs::exists(runtimes_filefullpath_))
			 {
				 converged = "false";
				 number_of_run_ = 1;
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
		 virtual ~RegressionTestBase();

		 void transferTheIndex();  /** transfer the index of a double vector. */
		 void readResultFromXml(); /** read the result from the .xml file. */
		 void writeResultToXml();  /** write the result to the .xml file. */
		 void readResultFromXml(int index_of_run_); /* read the result from the .xml file with the index. */
		 void writeResultToXml(int index_of_run_); /* write the result to the .xml file with the index. */

		 /** the interface to write observed quantity into xml memory. */
		 void writeToFile(size_t iteration = 0) override
		 {
			 ObserveMethodType::writeToFile(iteration);
			 this->writeToXml(iteration);
		 };

		 /** the interface to write xml memory into Xml file. */
		 void writeXmlToXmlFile()
		 {
			 this->observe_xml_engine_.writeToXmlFile(this->filefullpath_input_);
		 };

		 /** read current result from Xml file from Xml memory. */
		 void readXmlFromXmlFile()
		 {
			 this->readFromXml();
		 };
	 };
};