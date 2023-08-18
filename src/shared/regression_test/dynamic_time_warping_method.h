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
 * @file 	dynamic_time_warping_method.h
 * @brief 	Classes for the comparison between validated and tested results
                        with dynamic time warping method.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "time_average_method.hpp"

namespace SPH
{
/**
 * @class RegressionTestDynamicTimeWarping
 * @brief the regression test is based on the dynamic time warping.
 */
template <class ObserveMethodType>
class RegressionTestDynamicTimeWarping : public RegressionTestTimeAverage<ObserveMethodType>
{
    /*identify the variable type from the parent class. */
    using VariableType = decltype(ObserveMethodType::type_indicator_);

  protected:
    std::string dtw_distance_filefullpath_; /* the path for DTW distance. */
    XmlEngine dtw_distance_xml_engine_in_;  /* xml engine for dtw distance input. */
    XmlEngine dtw_distance_xml_engine_out_; /* xml engine for dtw distance output. */

    StdVec<Real> dtw_distance_, dtw_distance_new_; /* the container of DTW distance between each pairs. */

    /** the method used for calculating the p_norm. (calculateDTWDistance) */
    Real calculatePNorm(Real variable_a, Real variable_b)
    {
        return std::abs(variable_a - variable_b);
    };
    template <typename Variable>
    Real calculatePNorm(const Variable &variable_a, const Variable variable_b)
    {
        return (variable_a - variable_b).norm();
    };

    /** the local constrained method used for calculating the dtw distance between two lines. */
    StdVec<Real> calculateDTWDistance(BiVector<VariableType> dataset_a_, BiVector<VariableType> dataset_b_);

  public:
    template <typename... ConstructorArgs>
    explicit RegressionTestDynamicTimeWarping(ConstructorArgs &&...args) : RegressionTestTimeAverage<ObserveMethodType>(std::forward<ConstructorArgs>(args)...),
                                                                           dtw_distance_xml_engine_in_("dtw_distance_xml_engine_in", "dtw_distance"),
                                                                           dtw_distance_xml_engine_out_("dtw_distance_xml_engine_out", "dtw_distance")
    {
        dtw_distance_filefullpath_ = this->input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + "_dtwdistance.xml";
    };
    virtual ~RegressionTestDynamicTimeWarping(){};

    void setupTheTest();                           /** setup the test and defined basic variables. */
    void readDTWDistanceFromXml();                 /** read the old DTW distance from the .xml file. */
    void updateDTWDistance();                      /** update the maximum DTWDistance with the new result. */
    void writeDTWDistanceToXml();                  /* write the updated DTWDistance to .xml file.*/
    bool compareDTWDistance(Real threshold_value); /* compare the DTWDistance if converged. */
    void resultTest();                             /** test the new result if it is converged within the range. */

    /** the interface for generating the priori converged result with DTW */
    void generateDataBase(Real threshold_value, std::string filter = "false")
    {
        this->writeXmlToXmlFile();
        this->readXmlFromXmlFile();
        this->transposeTheIndex();
        if (this->converged == "false")
        {
            setupTheTest();
            if (filter == "true")
                this->filterExtremeValues();
            readDTWDistanceFromXml();
            /* loop all existed result to get maximum dtw distance. */
            for (int n = 0; n != (this->number_of_run_ - 1); ++n)
            {
                this->readResultFromXml(n);
                updateDTWDistance();
            }
            this->writeResultToXml(this->number_of_run_ - 1);
            writeDTWDistanceToXml();
            compareDTWDistance(threshold_value);
        }
        else
            std::cout << "The results have been converged." << std::endl;
    };

    /** the interface for generating the priori converged result with DTW. */
    void testResult(std::string filter = "false")
    {
        this->writeXmlToXmlFile();
        this->readXmlFromXmlFile();
        this->transposeTheIndex();
        setupTheTest();
        if (filter == "true")
            this->filterExtremeValues();
        readDTWDistanceFromXml();
        for (int n = 0; n != this->number_of_run_; ++n)
        {
            this->result_filefullpath_ = this->input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + "_Run_" + std::to_string(n) + "_result.xml";
            if (!fs::exists(this->result_filefullpath_))
            {
                std::cout << "This result has not been preserved and will not be compared." << std::endl;
                continue;
            }
            this->readResultFromXml(n);
            resultTest();
        }
        std::cout << "The result of " << this->quantity_name_
                  << " is correct based on the dynamic time warping regression test!" << std::endl;
    };
};
} // namespace SPH