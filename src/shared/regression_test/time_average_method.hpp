/**
 * @file 	time_averaged_method.hpp
 * @brief 	Classes for the comparison between validated and tested results
 *         	with time-averaged meanvalue and variance method.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "time_average_method.h"

namespace SPH
{
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::filterLocalResult(BiVector<Real> &current_result)
{
    int scale = round(this->snapshot_ / 200);
    std::cout << "The filter scale is " << scale * 2 << "." << std::endl;
    for (int snapshot_index = 0; snapshot_index != this->snapshot_; ++snapshot_index)
    {
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        {
            Real filter_meanvalue = 0;
            Real filter_variance = 0;
            for (int index = SMAX(snapshot_index - scale, 0); index != SMIN(snapshot_index + scale, this->snapshot_); ++index)
            {
                filter_meanvalue += current_result[index][observation_index];
            }
            filter_meanvalue = (filter_meanvalue - current_result[snapshot_index][observation_index]) /
                               (SMIN(snapshot_index + scale, this->snapshot_) - SMAX(snapshot_index - scale, 0));
            for (int index = SMAX(snapshot_index - scale, 0); index != SMIN(snapshot_index + scale, this->snapshot_); ++index)
            {
                filter_variance += pow(current_result[index][observation_index] - filter_meanvalue, 2);
            }
            Real current_variance = pow(current_result[snapshot_index][observation_index] - filter_meanvalue, 2);
            filter_variance = (filter_variance - current_variance) / (SMIN(snapshot_index + scale, this->snapshot_) - SMAX(snapshot_index - scale, 0));
            if (current_variance > 4 * filter_variance)
            {
                current_result[snapshot_index][observation_index] = filter_meanvalue;
                std::cout << "The current value of " << this->quantity_name_ << "[" << snapshot_index << "][" << observation_index << "] is "
                          << current_result[snapshot_index][observation_index]
                          << ", but the neighbor averaged value is " << filter_meanvalue << ", and the rate is " << current_variance / filter_variance << std::endl;
            }
        }
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::filterLocalResult(BiVector<Vecd> &current_result)
{
    int scale = round(this->snapshot_ / 200);
    std::cout << "The filter scale is " << scale * 2 << "." << std::endl;
    for (int snapshot_index = 0; snapshot_index != this->snapshot_; ++snapshot_index)
    {
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        {
            for (int i = 0; i != current_result[0][0].size(); ++i)
            {
                Real filter_meanvalue = 0;
                Real filter_variance = 0;
                for (int index = SMAX(snapshot_index - scale, 0); index != SMIN(snapshot_index + scale, this->snapshot_); ++index)
                {
                    filter_meanvalue += current_result[index][observation_index][i];
                }
                filter_meanvalue = (filter_meanvalue - current_result[snapshot_index][observation_index][i]) /
                                   (SMIN(snapshot_index + scale, this->snapshot_) - SMAX(snapshot_index - scale, 0));
                for (int index = SMAX(snapshot_index - scale, 0); index != SMIN(snapshot_index + scale, this->snapshot_); ++index)
                {
                    filter_variance += pow(current_result[index][observation_index][i] - filter_meanvalue, 2);
                }
                Real current_variance = pow(current_result[snapshot_index][observation_index][i] - filter_meanvalue, 2);
                filter_variance = (filter_variance - current_variance) / (SMIN(snapshot_index + scale, this->snapshot_) - SMAX(snapshot_index - scale, 0));
                if (current_variance > 4 * filter_variance)
                {
                    current_result[snapshot_index][observation_index][i] = filter_meanvalue;
                    std::cout << "The current value of " << this->quantity_name_ << "[" << snapshot_index << "][" << observation_index << "][" << i << "] is "
                              << current_result[snapshot_index][observation_index][i]
                              << ", but the neighbor averaged value is " << filter_meanvalue << ", and the rate is " << current_variance / filter_variance << std::endl;
                }
            }
        }
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::filterLocalResult(BiVector<Matd> &current_result)
{
    int scale = round(this->snapshot_ / 200);
    std::cout << "The filter scale is " << scale * 2 << "." << std::endl;
    for (int snapshot_index = 0; snapshot_index != this->snapshot_; ++snapshot_index)
    {
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        {
            for (int i = 0; i != current_result[0][0].size(); ++i)
            {
                for (int j = 0; j != current_result[0][0].size(); ++j)
                {
                    Real filter_meanvalue = 0;
                    Real filter_variance = 0;
                    for (int index = SMAX(snapshot_index - scale, 0); index != SMIN(snapshot_index + scale, this->snapshot_); ++index)
                    {
                        filter_meanvalue += current_result[index][observation_index](i, j);
                    }
                    filter_meanvalue = (filter_meanvalue - current_result[snapshot_index][observation_index](i, j)) /
                                       (SMIN(snapshot_index + scale, this->snapshot_) - SMAX(snapshot_index - scale, 0));
                    for (int index = SMAX(snapshot_index - scale, 0); index != SMIN(snapshot_index + scale, this->snapshot_); ++index)
                    {
                        filter_variance += pow(current_result[index][observation_index](i, j) - filter_meanvalue, 2);
                    }
                    Real current_variance = pow(current_result[snapshot_index][observation_index](i, j) - filter_meanvalue, 2);
                    filter_variance = (filter_variance - current_variance) / (SMIN(snapshot_index + scale, this->snapshot_) - SMAX(snapshot_index - scale, 0));
                    if (current_variance > 4 * filter_variance)
                    {
                        current_result[snapshot_index][observation_index](i, j) = filter_meanvalue;
                        std::cout << "The current value of " << this->quantity_name_ << "[" << snapshot_index << "]["
                                  << observation_index << "][" << i << "][" << j << "] is " << current_result[snapshot_index][observation_index](i, j)
                                  << ", but the neighbor averaged value is " << filter_meanvalue << ", and the rate is " << current_variance / filter_variance << std::endl;
                    }
                }
            }
        }
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::searchSteadyStart(BiVector<Real> &current_result)
{
    /* the search is only for one value. */
    int scale = round(this->snapshot_ / 20);
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        for (int snapshot_index = this->snapshot_ - 1; snapshot_index != 3 * scale; --snapshot_index)
        {
            Real value_one = 0, value_two = 0;
            for (int index = snapshot_index; index != snapshot_index - scale; --index)
            {
                value_one += current_result[index][observation_index] / scale;
                value_two += current_result[index - 2 * scale][observation_index] / scale;
            }

            if (ABS(value_one - value_two) / ABS((value_one + value_two) / 2) > 0.1)
            {
                snapshot_for_converged_ = SMAX(snapshot_for_converged_, snapshot_index - scale);
                break;
            }
        }
    std::cout << "The scale is " << scale << "." << std::endl;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::searchSteadyStart(BiVector<Vecd> &current_result)
{
    /* the search is for each value within parameters. */
    int scale = round(this->snapshot_ / 20);
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        for (int snapshot_index = this->snapshot_ - 1; snapshot_index != 3 * scale; --snapshot_index)
        {
            Real value_one = 0, value_two = 0;
            for (int index = snapshot_index; index != snapshot_index - scale; --index)
            {
                value_one += current_result[index][observation_index][0] / scale;
                value_two += current_result[index - 2 * scale][observation_index][0] / scale;
            }

            if (ABS(value_one - value_two) / ABS((value_one + value_two) / 2) > 0.1)
            {
                snapshot_for_converged_ = SMAX(snapshot_for_converged_, snapshot_index - scale);
                break; /** This break just jump out of dimension iteration.  */
            }
        }
    std::cout << "The scale is " << scale << "." << std::endl;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::searchSteadyStart(BiVector<Matd> &current_result)
{
    int scale = round(this->snapshot_ / 20);
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        for (int snapshot_index = this->snapshot_ - 1; snapshot_index != 3 * scale; --snapshot_index)
        {
            Real value_one = 0, value_two = 0;
            for (int index = snapshot_index; index != snapshot_index - scale; --index)
            {
                value_one += current_result[index][observation_index](0, 0);
                value_two += current_result[index - 2 * scale][observation_index](0, 0);
            }

            if (ABS(value_one - value_two) / ABS((value_one + value_two) / 2) > 0.1)
            {
                snapshot_for_converged_ = SMAX(snapshot_for_converged_, snapshot_index - scale);
                break; /** This break just jump out of dimension iteration.  */
            }
        }
    std::cout << "The scale is " << scale << "." << std::endl;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::calculateNewVariance(BiVector<Real> &current_result,
                                                                        StdVec<Real> &local_meanvalue, StdVec<Real> &variance, StdVec<Real> &variance_new)
{
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
            variance_new[observation_index] += pow((current_result[observation_index][snapshot_index] - local_meanvalue[observation_index]), 2);
        variance_new[observation_index] = SMAX(
            (Real)(variance_new[observation_index] / (this->snapshot_ - snapshot_for_converged_)),
            (Real)variance[observation_index],
            (Real)pow(local_meanvalue[observation_index] * Real(0.01), 2));
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::calculateNewVariance(BiVector<Vecd> &current_result,
                                                                        StdVec<Vecd> &local_meanvalue, StdVec<Vecd> &variance, StdVec<Vecd> &variance_new)
{
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        for (int i = 0; i != current_result[0][0].size(); ++i)
        {
            for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
                variance_new[observation_index][i] += pow((current_result[observation_index][snapshot_index][i] - local_meanvalue[observation_index][i]), 2);
            variance_new[observation_index][i] = SMAX(
                (Real)(variance_new[observation_index][i] / (this->snapshot_ - snapshot_for_converged_)),
                (Real)variance[observation_index][i],
                (Real)pow(local_meanvalue[observation_index][i] * Real(0.01), 2));
        }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::calculateNewVariance(BiVector<Matd> &current_result,
                                                                        StdVec<Matd> &local_meanvalue, StdVec<Matd> &variance, StdVec<Matd> &variance_new)
{
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        for (int i = 0; i != current_result[0][0].size(); ++i)
            for (int j = 0; j != current_result[0][0].size(); ++j)
            {
                for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
                    variance_new[observation_index](i, j) += pow((current_result[observation_index][snapshot_index](i, j) - local_meanvalue[observation_index](i, j)), 2);
                variance_new[observation_index](i, j) = SMAX(
                    (Real)(variance_new[observation_index](i, j) / (this->snapshot_ - snapshot_for_converged_)),
                    (Real)variance[observation_index](i, j),
                    (Real)pow(local_meanvalue[observation_index](i, j) * Real(0.01), 2));
            }
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestTimeAverage<ObserveMethodType>::compareParameter(std::string par_name,
                                                                   StdVec<Real> &parameter, StdVec<Real> &parameter_new, Real &threshold)
{
    int count = 0;
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        if ((par_name == "meanvalue") && (ABS(parameter[observation_index]) < 0.005) && (ABS(parameter_new[observation_index]) < 0.005))
        {
            std::cout << "The old meanvalue is " << parameter[observation_index] << ", and the new meanvalue is " << parameter_new[observation_index]
                      << ". So this variable will be ignored due to its tiny effect." << std::endl;
            continue;
        }
        Real relative_value_ = ABS((parameter[observation_index] - parameter_new[observation_index]) / (parameter_new[observation_index] + TinyReal));
        if (relative_value_ > threshold)
        {
            std::cout << par_name << ": " << this->quantity_name_ << "[" << observation_index << "]"
                      << " is not converged, and difference is " << relative_value_ << std::endl;
            count++;
        }
    }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestTimeAverage<ObserveMethodType>::compareParameter(std::string par_name,
                                                                   StdVec<Vecd> &parameter, StdVec<Vecd> &parameter_new, Vecd &threshold)
{
    int count = 0;
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        for (int i = 0; i != parameter[0].size(); ++i)
        {
            if ((par_name == "meanvalue") && (ABS(parameter[observation_index][i]) < 0.001) && (ABS(parameter_new[observation_index][i]) < 0.001))
            {
                std::cout << "The old meanvalue is " << parameter[observation_index][i] << ", and the new meanvalue is " << parameter_new[observation_index][i]
                          << ". So this variable will be ignored due to its tiny effect." << std::endl;
                continue;
            }
            Real relative_value_ = ABS((parameter[observation_index][i] - parameter_new[observation_index][i]) / (parameter_new[observation_index][i] + TinyReal));
            if (relative_value_ > threshold[i])
            {
                std::cout << par_name << ": " << this->quantity_name_ << "[" << observation_index << "][" << i << "]"
                          << " is not converged, and difference is " << relative_value_ << std::endl;
                count++;
            }
        }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestTimeAverage<ObserveMethodType>::compareParameter(std::string par_name,
                                                                   StdVec<Matd> &parameter, StdVec<Matd> &parameter_new, Matd &threshold)
{
    int count = 0;
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        for (int i = 0; i != parameter[0].size(); ++i)
            for (int j = 0; j != parameter[0].size(); ++j)
            {
                if ((par_name == "meanvalue") && (ABS(parameter[observation_index](i, j)) < 0.001) && (ABS(parameter_new[observation_index](i, j)) < 0.001))
                {
                    std::cout << "The old meanvalue is " << parameter[observation_index](i, j) << ", and the new meanvalue is "
                              << parameter_new[observation_index](i, j) << ". So this variable will be ignored due to its tiny effect." << std::endl;
                    continue;
                }
                Real relative_value_ = ABS((parameter[observation_index](i, j) - parameter_new[observation_index](i, j)) /
                                           (parameter_new[observation_index](i, j) + TinyReal));
                if (relative_value_ > threshold(i, j))
                {
                    std::cout << par_name << ": " << this->quantity_name_ << "[" << observation_index << "][" << i << "][" << j << "]"
                              << " is not converged, and difference is " << relative_value_ << std::endl;
                    count++;
                }
            }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestTimeAverage<ObserveMethodType>::testNewResult(BiVector<Real> &current_result,
                                                                StdVec<Real> &meanvalue, StdVec<Real> &local_meanvalue, StdVec<Real> &variance)
{
    int count = 0;
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
        {
            variance_new_[observation_index] += pow((current_result[snapshot_index][observation_index] - local_meanvalue[observation_index]), 2);
        }
        variance_new_[observation_index] = variance_new_[observation_index] / (this->snapshot_ - snapshot_for_converged_);
        if ((ABS(meanvalue[observation_index]) < 0.005) && (ABS(local_meanvalue[observation_index]) < 0.005))
        {
            std::cout << "The old meanvalue is " << meanvalue[observation_index] << ", and the current meanvalue is " << local_meanvalue[observation_index]
                      << ". So this variable will not be tested due to its tiny effect." << std::endl;
            continue;
        }
        Real relative_value_ = ABS((meanvalue[observation_index] - local_meanvalue[observation_index]) / (meanvalue[observation_index] + TinyReal));
        if (relative_value_ > 0.1 || variance_new_[observation_index] > (1.01 * variance[observation_index]))
        {
            std::cout << this->quantity_name_ << "[" << observation_index << "] is beyond the exception !" << std::endl;
            std::cout << "The meanvalue is " << meanvalue[observation_index] << ", and the current meanvalue is " << local_meanvalue[observation_index] << std::endl;
            std::cout << "The variance is " << variance[observation_index] << ", and the current variance is " << variance_new_[observation_index] << std::endl;
            count++;
        }
    }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestTimeAverage<ObserveMethodType>::testNewResult(BiVector<Vecd> &current_result,
                                                                StdVec<Vecd> &meanvalue, StdVec<Vecd> &local_meanvalue, StdVec<Vecd> &variance)
{
    int count = 0;
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        for (int i = 0; i != meanvalue_[0].size(); ++i)
        {
            for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
            {
                variance_new_[observation_index][i] += pow((current_result[snapshot_index][observation_index][i] -
                                                            local_meanvalue[observation_index][i]),
                                                           2);
            }
            variance_new_[observation_index][i] = variance_new_[observation_index][i] / (this->snapshot_ - snapshot_for_converged_);
            if ((ABS(meanvalue[observation_index][i]) < 0.005) && (ABS(local_meanvalue[observation_index][i]) < 0.005))
            {
                std::cout << "The old meanvalue is " << meanvalue[observation_index][i] << ", and the current meanvalue is "
                          << local_meanvalue[observation_index][i] << ". So this variable will not be tested due to its tiny effect." << std::endl;
                continue;
            }
            Real relative_value_ = ABS((meanvalue[observation_index][i] - local_meanvalue[observation_index][i]) / (meanvalue[observation_index][i] + TinyReal));
            if (relative_value_ > 0.1 || (variance_new_[observation_index][i] > 1.01 * variance[observation_index][i]))
            {
                std::cout << this->quantity_name_ << "[" << observation_index << "][" << i << "] is beyond the exception !" << std::endl;
                std::cout << "The meanvalue is " << meanvalue[observation_index][i] << ", and the current meanvalue is " << local_meanvalue[observation_index][i] << std::endl;
                std::cout << "The variance is " << variance[observation_index][i] << ", and the new variance is " << variance_new_[observation_index][i] << std::endl;
                count++;
            }
        }
    }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestTimeAverage<ObserveMethodType>::testNewResult(BiVector<Matd> &current_result,
                                                                StdVec<Matd> &meanvalue, StdVec<Matd> &local_meanvalue, StdVec<Matd> &variance)
{
    int count = 0;
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        for (int i = 0; i != meanvalue[0].size(); ++i)
        {
            for (int j = 0; j != meanvalue[0].size(); ++j)
            {
                for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
                {
                    variance_new_[observation_index](i, j) += pow((current_result[snapshot_index][observation_index](i, j) -
                                                                   local_meanvalue[observation_index](i, j)),
                                                                  2);
                }
                variance_new_[observation_index](i, j) = variance_new_[observation_index](i, j) / (this->snapshot_ - snapshot_for_converged_);
                if ((ABS(meanvalue[observation_index](i, j)) < 0.005) && (ABS(local_meanvalue[observation_index](i, j)) < 0.005))
                {
                    std::cout << "The old meanvalue is " << meanvalue[observation_index](i, j) << ", and the new meanvalue is "
                              << local_meanvalue[observation_index](i, j) << ". So this variable will not be tested due to its tiny effect. " << std::endl;
                    continue;
                }
                Real relative_value_ = ABS((meanvalue_[observation_index](i, j) - local_meanvalue[observation_index](i, j)) / (meanvalue[observation_index](i, j) + TinyReal));
                if (relative_value_ > 0.1 || variance_new_[observation_index](i, j) > 1.01 * variance[observation_index](i, j))
                {
                    std::cout << this->quantity_name_ << "[" << observation_index << "][" << i << "][" << j << "] is beyond the exception !" << std::endl;
                    std::cout << "The meanvalue is " << meanvalue[observation_index](i, j) << ", and the new meanvalue is " << local_meanvalue[observation_index](i, j) << std::endl;
                    std::cout << "The variance is " << variance[observation_index](i, j) << ", and the new variance is " << variance_new_[observation_index](i, j) << std::endl;
                    count++;
                }
            }
        }
    }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::initializeThreshold(VariableType &threshold_mean, VariableType &threshold_variance)
{
    threshold_mean_ = threshold_mean;
    threshold_variance_ = threshold_variance;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::setupTheTest()
{
    this->snapshot_ = this->current_result_.size();
    this->observation_ = this->current_result_[0].size();
    StdVec<VariableType> temp(this->observation_);
    meanvalue_ = temp;
    variance_ = temp;
    local_meanvalue_ = temp;
    meanvalue_new_ = meanvalue_;
    variance_new_ = variance_;

    if ((this->number_of_run_ > 1) && (!fs::exists(mean_variance_filefullpath_)))
    {
        std::cout << "\n Error: the input file:" << mean_variance_filefullpath_ << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::readMeanVarianceFromXml()
{
    if (this->number_of_run_ > 1)
    {
        mean_variance_xml_engine_in_.loadXmlFile(mean_variance_filefullpath_);
        SimTK::Xml::Element meanvalue_element_ = mean_variance_xml_engine_in_.getChildElement("MeanValue_Element");
        SimTK::Xml::element_iterator ele_ite_mean_ = meanvalue_element_.element_begin();
        for (; ele_ite_mean_ != meanvalue_element_.element_end(); ++ele_ite_mean_)
            for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
            {
                std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(observation_index);
                mean_variance_xml_engine_in_.getRequiredAttributeValue(ele_ite_mean_, attribute_name_, meanvalue_[observation_index]);
            }

        SimTK::Xml::Element variance_element_ = mean_variance_xml_engine_in_.getChildElement("Variance_Element");
        SimTK::Xml::element_iterator ele_ite_variance_ = variance_element_.element_begin();
        for (; ele_ite_variance_ != variance_element_.element_end(); ++ele_ite_variance_)
            for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
            {
                std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(observation_index);
                mean_variance_xml_engine_in_.getRequiredAttributeValue(ele_ite_variance_, attribute_name_, variance_[observation_index]);
            }
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::searchForStartPoint()
{
    snapshot_for_converged_ = 0;
    searchSteadyStart(this->current_result_);
    std::cout << "The snapshot for converged is " << snapshot_for_converged_ << std::endl;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::filterExtremeValues()
{
    filterLocalResult(this->current_result_);
    filefullpath_filter_output_ = this->input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + ".dat";
    std::ofstream out_file(filefullpath_filter_output_.c_str(), std::ios::app);
    out_file << "run_time"
             << "   ";
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        std::string quantity_name_i = this->quantity_name_ + "[" + std::to_string(observation_index) + "]";
        this->plt_engine_.writeAQuantityHeader(out_file, this->current_result_[0][0], quantity_name_i);
    }
    out_file << "\n";
    out_file.close();

    for (int snapshot_index = 0; snapshot_index != this->snapshot_; ++snapshot_index)
    {
        std::ofstream out_file(filefullpath_filter_output_.c_str(), std::ios::app);
        out_file << this->element_tag_[snapshot_index] << "   ";
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        {
            this->plt_engine_.writeAQuantity(out_file, this->current_result_[snapshot_index][observation_index]);
        }
        out_file << "\n";
        out_file.close();
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::updateMeanVariance()
{
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
        {
            local_meanvalue_[observation_index] += this->current_result_[snapshot_index][observation_index];
        }
        local_meanvalue_[observation_index] = local_meanvalue_[observation_index] / (this->snapshot_ - snapshot_for_converged_);
        meanvalue_new_[observation_index] = (local_meanvalue_[observation_index] + meanvalue_[observation_index] * (this->number_of_run_ - 1)) / this->number_of_run_;
    }
    calculateNewVariance(this->current_result_trans_, local_meanvalue_, variance_, variance_new_);
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::writeMeanVarianceToXml()
{
    mean_variance_xml_engine_out_.addElementToXmlDoc("MeanValue_Element");
    SimTK::Xml::Element meanvalue_element_ = mean_variance_xml_engine_out_.getChildElement("MeanValue_Element");
    mean_variance_xml_engine_out_.addChildToElement(meanvalue_element_, "Snapshot_MeanValue");
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        SimTK::Xml::element_iterator ele_ite_mean = meanvalue_element_.element_begin();
        std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(observation_index);
        mean_variance_xml_engine_out_.setAttributeToElement(ele_ite_mean, attribute_name_, meanvalue_new_[observation_index]);
    }
    mean_variance_xml_engine_out_.addElementToXmlDoc("Variance_Element");
    SimTK::Xml::Element variance_element_ = mean_variance_xml_engine_out_.getChildElement("Variance_Element");
    mean_variance_xml_engine_out_.addChildToElement(variance_element_, "Snapshot_Variance");
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        SimTK::Xml::element_iterator ele_ite_variance = variance_element_.element_begin();
        std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(observation_index);
        mean_variance_xml_engine_out_.setAttributeToElement(ele_ite_variance, attribute_name_, variance_new_[observation_index]);
    }
    mean_variance_xml_engine_out_.writeToXmlFile(mean_variance_filefullpath_);
}
//=================================================================================================//
template <class ObserveMethodType>
bool RegressionTestTimeAverage<ObserveMethodType>::compareMeanVariance()
{
    int count_not_converged_m = 0;
    int count_not_converged_v = 0;
    count_not_converged_m = this->compareParameter("meanvalue", this->meanvalue_, this->meanvalue_new_, this->threshold_mean_);
    count_not_converged_v = this->compareParameter("variance", this->variance_, this->variance_new_, this->threshold_variance_);
    if (count_not_converged_m == 0)
    {
        std::cout << "The meanvalue of " << this->quantity_name_ << " are converged now." << std::endl;
        if (count_not_converged_v == 0)
        {
            if (this->label_for_repeat_ == 4)
            {
                this->converged = "true";
                this->label_for_repeat_++;
                std::cout << "The meanvalue and variance of " << this->quantity_name_ << " are converged enough times, and run will stop now." << std::endl;
                return true;
            }
            else
            {
                this->converged = "false";
                this->label_for_repeat_++;
                std::cout << "The variance of " << this->quantity_name_ << " are also converged, and this is the " << this->label_for_repeat_
                          << " times. They should be converged more times to be stable." << std::endl;
                return false;
            }
        }
        else
        {
            this->converged = "false";
            this->label_for_repeat_ = 0;
            std::cout << "The variance of " << this->quantity_name_ << " are not converged " << count_not_converged_v << " times." << std::endl;
            return false;
        };
    }
    else
    {
        this->converged = "false";
        this->label_for_repeat_ = 0;
        std::cout << "The meanvalue of " << this->quantity_name_ << " are not converged " << count_not_converged_m << " times." << std::endl;
        return false;
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::resultTest()
{
    int test_wrong = 0;

    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
            local_meanvalue_[observation_index] += this->current_result_[snapshot_index][observation_index];
        local_meanvalue_[observation_index] = local_meanvalue_[observation_index] / (this->snapshot_ - snapshot_for_converged_);
    }

    test_wrong = testNewResult(this->current_result_, meanvalue_, local_meanvalue_, variance_);
    if (test_wrong == 0)
        std::cout << "The result of " << this->quantity_name_ << " is correct based on the time-averaged regression test!" << std::endl;
    else
    {
        std::cout << "There are " << test_wrong << " particles are not within the expected range." << std::endl;
        std::cout << "Please try again. If it still post this conclusion, the result is not correct!" << std::endl;
        exit(1);
    }
};
//=================================================================================================//
}; // namespace SPH