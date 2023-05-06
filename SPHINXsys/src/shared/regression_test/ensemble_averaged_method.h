/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	ensemble_averaged_method.h
 * @brief 	Classes for the comparison between validated and tested results
          	with ensemble-averaged mean value and variance method.
 * @author	Bo Zhang , Chi ZHang and Xiangyu Hu
 */

#pragma once
#include "regression_test_base.hpp"

namespace SPH
{
	/**
	 * @class RegressionTestEnsembleAveraged
	 * @brief the regression test is based on the ensemble-averaged mean value and variance.
	 */
	template<class ObserveMethodType>
	class RegressionTestEnsembleAveraged : public RegressionTestTimeAveraged<ObserveMethodType>
	{
		/* identify the variable type from the parent class. */
		using VariableType = decltype(ObserveMethodType::type_indicator_);

	protected:
		DoubleVec<VariableType> meanvalue_, meanvalue_new_;  /* the container of (new) mean value. [different from time-averaged]*/
		DoubleVec<VariableType> variance_, variance_new_;    /* the container of (new) variance. [different from time-averaged]*/
		
		/** the method used for calculating the new variance. */
		void calculateNewVariance(TripleVec<Real> &result, DoubleVec<Real> &meanvalue_new, DoubleVec<Real> &variance, DoubleVec<Real> &variance_new);
		void calculateNewVariance(TripleVec<Vecd> &result, DoubleVec<Vecd> &meanvalue_new, DoubleVec<Vecd> &variance, DoubleVec<Vecd> &variance_new);
		void calculateNewVariance(TripleVec<Matd> &result, DoubleVec<Matd> &meanvalue_new, DoubleVec<Matd> &variance, DoubleVec<Matd> &variance_new);

		/** the method used for comparing the meanvalue and variance. */
		int compareParameter(std::string par_name, DoubleVec<Real> &parameter, DoubleVec<Real> &parameter_new, Real &threshold);
		int compareParameter(std::string par_name, DoubleVec<Vecd> &parameter, DoubleVec<Vecd> &parameter_new, Vecd &threshold);
		int compareParameter(std::string par_name, DoubleVec<Matd> &parameter, DoubleVec<Matd> &parameter_new, Matd &threshold);
		
		/** the method used for testing the new result with meanvalue and variance. */
		int testNewResult(int diff, DoubleVec<Real> &current_result, DoubleVec<Real> &meanvalue, DoubleVec<Real> &variance);
		int testNewResult(int diff, DoubleVec<Vecd> &current_result, DoubleVec<Vecd> &meanvalue, DoubleVec<Vecd> &variance);
		int testNewResult(int diff, DoubleVec<Matd> &current_result, DoubleVec<Matd> &meanvalue, DoubleVec<Matd> &variance);

	public:
		template<typename... ConstructorArgs>
		explicit RegressionTestEnsembleAveraged(ConstructorArgs &&...args) :
			RegressionTestTimeAveraged<ObserveMethodType>(std::forward<ConstructorArgs>(args)...)
		{
			this->mean_variance_filefullpath_ = this->input_folder_path_ + "/" + this->dynamics_identifier_name_ 
				+ "_" + this->quantity_name_ + "_ensemble_averaged_mean_variance.xml";
		};
		virtual ~RegressionTestEnsembleAveraged() {};

		void setupAndCorrection(); /** setup and correct the number of old and new result. */
		void readMeanVarianceFromXml(); /** read the meanvalue and variance from the .xml file. */
		void updateMeanVariance(); /** update the meanvalue and variance from new result. */
		void writeMeanVarianceToXml(); /** write the meanvalue and variance to the .xml file. */
		bool compareMeanVariance();  /** compare the meanvalue and variance between old and new ones. */
		void resultTest(); /** test the new result if it is converged within the range. */

		/* the interface for generating the priori converged result with M&V. */
		void generateDataBase(VariableType threshold_mean, VariableType threshold_variance, std::string filter = "false")
		{
			this->writeXmlToXmlFile();
			this->readXmlFromXmlFile();
			this->initializeThreshold(threshold_mean, threshold_variance);
			if (this->converged == "false")
			{
				setupAndCorrection();
				this->readResultFromXml();
				if (filter == "true")
					this->filterExtremeValues();
				readMeanVarianceFromXml();
				updateMeanVariance();
				this->writeResultToXml();
				writeMeanVarianceToXml();
				compareMeanVariance();
			};
			/*else
				std::cout << "The results have been converged." << endl;*/
		};

		/** the interface for testing new result. */
		void newResultTest(std::string filter = "false")
		{
			this->writeXmlToXmlFile();
			this->readXmlFromXmlFile();
			setupAndCorrection();
			if (filter == "true")
				this->filterExtremeValues();
			readMeanVarianceFromXml();
			resultTest();
		};
	};
};