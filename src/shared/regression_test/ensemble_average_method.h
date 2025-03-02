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

#include "time_average_method.h"

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
    void calculateNewVariance(TriVector<Real> &result, BiVector<Real> &meanvalue_new, BiVector<Real> &variance, BiVector<Real> &variance_new);
    void calculateNewVariance(TriVector<Vecd> &result, BiVector<Vecd> &meanvalue_new, BiVector<Vecd> &variance, BiVector<Vecd> &variance_new);
    void calculateNewVariance(TriVector<Matd> &result, BiVector<Matd> &meanvalue_new, BiVector<Matd> &variance, BiVector<Matd> &variance_new);

    /** the method used for comparing the meanvalue and variance. */
    int compareParameter(std::string par_name, BiVector<Real> &parameter, BiVector<Real> &parameter_new, Real &threshold);
    int compareParameter(std::string par_name, BiVector<Vecd> &parameter, BiVector<Vecd> &parameter_new, Vecd &threshold);
    int compareParameter(std::string par_name, BiVector<Matd> &parameter, BiVector<Matd> &parameter_new, Matd &threshold);

    /** the method used for testing the new result with meanvalue and variance. */
    int testNewResult(int diff, BiVector<Real> &current_result, BiVector<Real> &meanvalue, BiVector<Real> &variance);
    int testNewResult(int diff, BiVector<Vecd> &current_result, BiVector<Vecd> &meanvalue, BiVector<Vecd> &variance);
    int testNewResult(int diff, BiVector<Matd> &current_result, BiVector<Matd> &meanvalue, BiVector<Matd> &variance);

    template <typename T>
    void readDataFromXmlMemory(XmlParser &xml_parser, tinyxml2::XMLElement *element,
                               int observe_k, BiVector<T> &result_container, const std::string &quantity_name)
    {
        int snapshot_index = 0;
        for (auto child = element->FirstChildElement(); child; child = child->NextSiblingElement())
        {
            std::string attribute_name = quantity_name + "_" + std::to_string(observe_k);
            xml_parser.queryAttributeValue(child, attribute_name, result_container[snapshot_index][observe_k]);
            snapshot_index++;
        }
    };

    template <typename T>
    void writeDataToXmlMemory(XmlParser &xml_parser, tinyxml2::XMLElement *element, const BiVector<T> &quantity,
                              int snapshot_, int observe_, const std::string &quantity_name, StdVec<std::string> &element_tag)
    {
        for (int snapshot_index = 0; snapshot_index != snapshot_; ++snapshot_index)
        {
            std::string child_element_name = element_tag[snapshot_index];
            auto child_element = xml_parser.addNewElement(element, child_element_name);
            for (int observe_k = 0; observe_k != observe_; ++observe_k)
            {
                std::string attribute_name = quantity_name + "_" + std::to_string(observe_k);
                xml_parser.setAttributeToElement(child_element, attribute_name, quantity[snapshot_index][observe_k]);
            }
        }
    };

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
    void testResult(const std::string &filter = "false")
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
}; // namespace SPH
#endif // ENSEMBLE_AVERAGE_H
