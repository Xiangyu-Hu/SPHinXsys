/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	ensemble_averaged_method.h
 * @brief 	Classes for the comparison between validated and tested results
                with ensemble-averaged mean value and variance method.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#ifndef ENSEMBLE_AVERAGE_H
#define ENSEMBLE_AVERAGE_H

#include "regression_test_base.hpp"

namespace SPH
{
/**
 * @class RegressionTestEnsembleAverage
 * @brief the regression test is based on the ensemble-averaged mean value and variance.
 */
template <class ObserveMethodType>
class RegressionTestEnsembleAverage : public RegressionTestTimeAverage<ObserveMethodType>
{
    /* identify the variable type from the parent class. */
    using VariableType = decltype(ObserveMethodType::type_indicator_);

  protected:
    BiVector<VariableType> meanvalue_, meanvalue_new_; /* the container of (new) mean value. [different from time-averaged]*/
    BiVector<VariableType> variance_, variance_new_;   /* the container of (new) variance. [different from time-averaged]*/

    /** the method used for calculating the new variance. */
    void calculateNewVariance(TriVector<VariableType> &result);
    /** the method used for comparing the meanvalue and variance. */
    int compareParameter(std::string par_name, BiVector<VariableType> &parameter,
                         BiVector<VariableType> &parameter_new, VariableType &threshold);
    /** the method used for testing the new result with meanvalue and variance. */
    int testNewResult(int diff, BiVector<VariableType> &current_result,
                      BiVector<VariableType> &meanvalue, BiVector<VariableType> &variance);

  public:
    template <typename... Args>
    explicit RegressionTestEnsembleAverage(Args &&...args)
        : RegressionTestTimeAverage<ObserveMethodType>(std::forward<Args>(args)...)
    {
        this->mean_variance_filefullpath_ = this->input_folder_path_ + "/" + this->dynamics_identifier_name_ +
                                            "_" + this->quantity_name_ + "_ensemble_averaged_mean_variance.xml";
    };
    virtual ~RegressionTestEnsembleAverage() {};

    void setupAndCorrection();      /** setup and correct the number of old and new result. */
    void readMeanVarianceFromXml(); /** read the meanvalue and variance from the .xml file. */
    void updateMeanVariance();      /** update the meanvalue and variance from new result. */
    void writeMeanVarianceToXml();  /** write the meanvalue and variance to the .xml file. */
    bool compareMeanVariance();     /** compare the meanvalue and variance between old and new ones. */
    void resultTest();              /** test the new result if it is converged within the range. */

    /* the interface for generating the priori converged result with M&V. */
    void generateDataBase(VariableType threshold_mean, VariableType threshold_variance, const std::string &filter = "false")
    {
        this->initializeThreshold(threshold_mean, threshold_variance);
        if (this->converged_ == "false")
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
    void testResult(const std::string &filter = "false")
    {
        setupAndCorrection();
        if (filter == "true")
            this->filterExtremeValues();
        readMeanVarianceFromXml();
        resultTest();
    };
};
}; // namespace SPH
#endif // ENSEMBLE_AVERAGE_H
