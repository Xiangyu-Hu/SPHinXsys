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
 * @file time_averaged_method.h
 * @brief Classes for the comparison between validated and tested results
          with time-averaged meanvalue and variance method.
 * @author Bo Zhang and Xiangyu Hu
 */

#pragma once
#include "in_output.h"
#include "xml_engine.h"
#include "all_physical_dynamics.h"

namespace SPH
{
	/**
	 * @class RegressionTestTimeAveraged
	 * @brief The regression test is based on the time-averaged meanvalue and variance.
	 */
	template<class ObserveMethodType>
	class RegressionTestTimeAveraged : public RegressionTestBase<ObserveMethodType>
	{
		/* identify the variable type from the parent class. */
		using VariableType = decltype(ObserveMethodType::type_indicator_);

	protected:
		int snapshot_for_converged_ ;                 /* the label for steady converged starting point. */
		std::string mean_variance_filefullpath_;      /* the path for mean and variance. */
		std::string filefullpath_filter_output_;      /* the path for filtered output. */
		XmlEngine mean_variance_xml_engine_in_;       /* xml engine for mean and variance input. */
		XmlEngine mean_variance_xml_engine_out_;      /* xml engine for mean and variance output. */
		
		VariableType threshold_mean_, threshold_variance_;  /* the threshold container for mean and variance. */
		StdVec<VariableType> meanvalue_, meanvalue_new_;    /* the container of (new) meanvalue. */
		StdVec<VariableType> variance_, variance_new_;      /* the container of (new) variance. */
		StdVec<VariableType> local_meanvalue_;              /* the container for meanvalue of current result. */
		
		/** the method used for filtering local extreme values. */
		void filterLocalResult(DoubleVec<Real> &current_result);
		void filterLocalResult(DoubleVec<Vecd> &current_result);
		void filterLocalResult(DoubleVec<Matd> &current_result);

		/** the method used for test steady starting point. */
		void searchSteadyStart(DoubleVec<Real> &current_result);
		void searchSteadyStart(DoubleVec<Vecd> &current_result);
		void searchSteadyStart(DoubleVec<Matd> &current_result);

		/** the method used for calculating the new variance. */
		void calculateNewVariance(DoubleVec<Real> &current_result, StdVec<Real> &local_meanvalue, StdVec<Real> &variance, StdVec<Real> &variance_new);
		void calculateNewVariance(DoubleVec<Vecd> &current_result, StdVec<Vecd> &local_meanvalue, StdVec<Vecd> &variance, StdVec<Vecd> &variance_new);
		void calculateNewVariance(DoubleVec<Matd> &current_result, StdVec<Matd> &local_meanvalue, StdVec<Matd> &variance, StdVec<Matd> &variance_new);

		/** the method used for comparing the meanvalue and variance. */
		int compareParameter(string par_name, StdVec<Real> &parameter, StdVec<Real> &parameter_new, Real &threshold);
		int compareParameter(string par_name, StdVec<Vecd> &parameter, StdVec<Vecd> &parameter_new, Vecd &threshold);
		int compareParameter(string par_name, StdVec<Matd> &parameter, StdVec<Matd> &parameter_new, Matd &threshold);

		/** the method used for testing the new result with meanvalue and variance. */
		int testNewResult(DoubleVec<Real> &current_result, StdVec<Real> &meanvalue, StdVec<Real> &local_meanvalue, StdVec<Real> &variance);
		int testNewResult(DoubleVec<Vecd> &current_result, StdVec<Vecd> &meanvalue, StdVec<Vecd> &local_meanvalue, StdVec<Vecd> &variance);
		int testNewResult(DoubleVec<Matd> &current_result, StdVec<Matd> &meanvalue, StdVec<Matd> &local_meanvalue, StdVec<Matd> &variance);
		
	public:
		template<typename... ConstructorArgs>
		explicit RegressionTestTimeAveraged(ConstructorArgs &&...args) :
			RegressionTestBase<ObserveMethodType>(std::forward<ConstructorArgs>(args)...),
			mean_variance_xml_engine_in_("mean_variance_xml_engine_in_", "meanvariance"),
			mean_variance_xml_engine_out_("mean_variance_xml_engine_out_", "meanvariance")
		{
			mean_variance_filefullpath_ = this->input_folder_path_ + "/" + this->body_name_
				+ "_" + this->quantity_name_ + "_time_averaged_mean_variance.xml";
		};
		virtual ~RegressionTestTimeAveraged() {};

        /** initialize the threshold of meanvalue and variance. */
		void initializeThreshold(VariableType &threshold_mean, VariableType &threshold_variance);
		void settingupTheTest();  /** setup the test and define basic variables. */
		void readMeanVarianceFromXml();  /** read the mean and variance from the .xml file. */
		void searchForStartPoint(); /** search for the starting point of the steady result. */
		void filterExtremeValues();  /** filter out the extreme values. */
		void updateMeanVariance();  /** update the meanvalue and variance from new result. */
		void writeMeanVarianceToXml();  /** write the meanvalue and variance to the .xml file. */
		bool compareMeanVariance();  /** compare the meanvalue and variance between old and new ones. */
		void resultTest();  /** test the new result if it is converged within the range. */

		/* the interface for generating the priori converged result with time-averaged meanvalue and variance. */
		void generateDataBase(VariableType threshold_mean, VariableType threshold_variance, string filter = "false")
		{
			this->writeXmlToXmlFile();
			this->readXmlFromXmlFile();
			initializeThreshold(threshold_mean, threshold_variance);
			if (this->converged == "false")
			{
				settingupTheTest();
				if (filter == "true")
					filterExtremeValues();
				searchForStartPoint();
				this->transferTheIndex();
				readMeanVarianceFromXml();
				updateMeanVariance();
				this->writeResultToXml(this->number_of_run_ - 1);
				writeMeanVarianceToXml();
				compareMeanVariance();
			}
		};

		/** the interface for testing new result. */
		void newResultTest(string filter = "false")
		{
			this->writeXmlToXmlFile();
			this->readXmlFromXmlFile();
			settingupTheTest();
			if (filter == "true")
				filterExtremeValues();
			searchForStartPoint();
			readMeanVarianceFromXml();
			resultTest();
		}
	};
}