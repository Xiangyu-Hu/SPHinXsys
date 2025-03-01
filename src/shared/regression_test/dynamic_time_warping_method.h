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
 * @author	Bo Zhang, Chi Zhang and Xiangyu Hu
 */

#ifndef DYNAMIC_TIME_WARPING_H
#define DYNAMIC_TIME_WARPING_H

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

    /** the local constrained method used for calculating the dtw distance between two lines. */
    StdVec<Real> calculateDTWDistance(BiVector<VariableType> dataset_a_, BiVector<VariableType> dataset_b_);

  public:
    template <typename... Args>
    explicit RegressionTestDynamicTimeWarping(Args &&...args);
    virtual ~RegressionTestDynamicTimeWarping() {};

    void setupTest();                           /** setup the test and defined basic variables. */
    void readDTWDistanceFromXml();                 /** read the old DTW distance from the .xml file. */
    void updateDTWDistance();                      /** update the maximum DTWDistance with the new result. */
    void writeDTWDistanceToXml();                  /* write the updated DTWDistance to .xml file.*/
    bool compareDTWDistance(Real threshold_value); /* compare the DTWDistance if converged. */
    void resultTest();                             /** test the new result if it is converged within the range. */

    /** the interface for generating the priori converged result with DTW */
    void generateDataBase(Real threshold_value, const std::string &filter = "false");

    /** the interface for generating the priori converged result with DTW. */
    void testResult(const std::string &filter = "false");
};
} // namespace SPH
#endif // DYNAMIC_TIME_WARPING_H
