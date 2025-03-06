/**
 * @file 	ensemble_averaged_method.hpp
 * @brief 	Classes for the comparison between validated and tested results
                with ensemble-averaged mean value and variance method.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "ensemble_average_method.h"

namespace SPH
{
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestEnsembleAverage<ObserveMethodType>::calculateNewVariance(TriVector<VariableType> &result)
{
    for (int snapshot_index = 0; snapshot_index != SMIN(this->snapshot_, this->number_of_snapshot_old_); ++snapshot_index)
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
            for (int run_index = 0; run_index != this->number_of_run_; ++run_index)
            {
                variance_new_[snapshot_index][observation_index] = transform_component(
                    variance_new_[snapshot_index][observation_index],
                    [&](Real variance_new, Real variance, Real result, Real meanvalue_new)
                    { return SMAX(
                          variance_new,
                          variance,
                          (Real)pow(result - meanvalue_new, 2),
                          (Real)pow(meanvalue_new * Real(0.01), 2)); },
                    variance_[snapshot_index][observation_index], result[run_index][snapshot_index][observation_index],
                    meanvalue_new_[snapshot_index][observation_index]);
            }
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestEnsembleAverage<ObserveMethodType>::compareParameter(
    std::string par_name, BiVector<VariableType> &parameter, BiVector<VariableType> &parameter_new, VariableType &threshold)
{
    int count = 0;
    for (int snapshot_index = 0; snapshot_index != SMIN(this->snapshot_, this->number_of_snapshot_old_); ++snapshot_index)
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        {
            for_each_component(
                parameter[snapshot_index][observation_index], [&](Real para, Real para_new, Real thr)
                {
                    Real relative_value = ABS(para - para_new) / (para + TinyReal);
                    if (relative_value > thr)
                    {
                        std::cout << par_name << ": " << this->quantity_name_ << "[" << observation_index << "]"
                                  << this->element_tag_[snapshot_index] << " is not converged, and difference is "
                                  << relative_value << std::endl;
                        count++;
                    } },
                parameter_new[snapshot_index][observation_index], threshold);
        }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestEnsembleAverage<ObserveMethodType>::testNewResult(
    int diff, BiVector<VariableType> &current_result,
    BiVector<VariableType> &meanvalue, BiVector<VariableType> &variance)
{
    int count = 0;
    std::cout << "The current length difference is " << diff << "." << std::endl;
    for (int snapshot_index = 0; snapshot_index != SMIN(this->snapshot_, this->number_of_snapshot_old_); ++snapshot_index)
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        {
            for_each_component(
                current_result[snapshot_index][observation_index], [&](Real cur_result, Real mean, Real var)
                {
                    Real relative_value = (pow(cur_result - mean, 2) - var) / var ;
                    if (relative_value > 0.01)
                    {
                        std::cout << this->quantity_name_ << "[" << observation_index << "]" 
                        << this->element_tag_[snapshot_index] << " is beyond the exception, and difference is " 
                        << relative_value << std::endl;
                        count++;
                    } },
                meanvalue[snapshot_index + diff][observation_index],
                variance[snapshot_index + diff][observation_index]);
        }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestEnsembleAverage<ObserveMethodType>::setupAndCorrection()
{
    if (this->snapshot_ == 0)
    {
        std::cout << "\n Error: the current results for ensemble-average regression test is empty!" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }

    if (this->number_of_run_ > 1)
    {
        if (this->converged_ == "false") /*< To identify the database generation or new result testing. */
        {
            if (!fs::exists(this->result_filefullpath_))
            {
                std::cout << "\n Error: the input file:" << this->result_filefullpath_ << " is not exists" << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }
            else
                this->result_xml_engine_in_.loadXmlFile(this->result_filefullpath_);
        }

        if (!fs::exists(this->mean_variance_filefullpath_))
        {
            std::cout << "\n Error: the input file:" << this->mean_variance_filefullpath_ << " is not exists" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }
        else
        {
            this->mean_variance_xml_engine_in_.loadXmlFile(this->mean_variance_filefullpath_);
            SimTK::Xml::Element mean_element_ = this->mean_variance_xml_engine_in_.getChildElement("Mean_Element");
            this->number_of_snapshot_old_ = std::distance(mean_element_.element_begin(), mean_element_.element_end());

            BiVector<VariableType> temp(SMAX(this->snapshot_, this->number_of_snapshot_old_), StdVec<VariableType>(this->observation_));
            meanvalue_ = temp;
            variance_ = temp;

            /** Unify the length of current result and previous result. */
            if (this->number_of_snapshot_old_ < this->snapshot_)
            {
                this->snapshot_difference_ = this->snapshot_ - this->number_of_snapshot_old_;
                for (int delete_ = 0; delete_ != this->snapshot_difference_; ++delete_)
                    this->current_result_.pop_back();
            }
            else if (this->number_of_snapshot_old_ > this->snapshot_)
                this->snapshot_difference_ = this->number_of_snapshot_old_ - this->snapshot_;
            else
                this->snapshot_difference_ = 0;
        }
    }
    else if (this->number_of_run_ == 1)
    {
        this->number_of_snapshot_old_ = this->snapshot_;
        BiVector<VariableType> temp(this->snapshot_, StdVec<VariableType>(this->observation_));
        this->result_.push_back(this->current_result_);
        meanvalue_ = temp;
        variance_ = temp;
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestEnsembleAverage<ObserveMethodType>::readMeanVarianceFromXml()
{
    if (this->number_of_run_ > 1)
    {
        SimTK::Xml::Element mean_element_ = this->mean_variance_xml_engine_in_.getChildElement("Mean_Element");
        SimTK::Xml::Element variance_element_ = this->mean_variance_xml_engine_in_.getChildElement("Variance_Element");
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        {
            this->xmlmemory_io_.readDataFromXmlMemory(this->mean_variance_xml_engine_in_,
                                                      mean_element_, observation_index, this->meanvalue_, this->quantity_name_);
            this->xmlmemory_io_.readDataFromXmlMemory(this->mean_variance_xml_engine_in_,
                                                      variance_element_, observation_index, this->variance_, this->quantity_name_);
        }
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestEnsembleAverage<ObserveMethodType>::updateMeanVariance()
{
    /** Unify the length of result and meanvalue. */
    if (this->number_of_run_ > 1)
    {
        for (int delete_ = 0; delete_ != this->snapshot_difference_; ++delete_)
        {
            meanvalue_.pop_back();
            variance_.pop_back();
        }
    }
    meanvalue_new_ = meanvalue_;
    variance_new_ = variance_;

    /** update the meanvalue of the result. */
    for (int snapshot_index = 0; snapshot_index != SMIN(this->snapshot_, this->number_of_snapshot_old_); ++snapshot_index)
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
            meanvalue_new_[snapshot_index][observation_index] = (meanvalue_[snapshot_index][observation_index] * (this->number_of_run_ - 1) +
                                                                 this->current_result_[snapshot_index][observation_index]) /
                                                                this->number_of_run_;
    /** Update the variance of the result. */
    calculateNewVariance(this->result_);
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestEnsembleAverage<ObserveMethodType>::writeMeanVarianceToXml()
{
    this->mean_variance_xml_engine_out_.addElementToXmlDoc("Mean_Element");
    SimTK::Xml::Element mean_element_ = this->mean_variance_xml_engine_out_.getChildElement("Mean_Element");
    this->xmlmemory_io_.writeDataToXmlMemory(this->mean_variance_xml_engine_out_, mean_element_, this->meanvalue_new_,
                                             SMIN(this->snapshot_, this->number_of_snapshot_old_), this->observation_, this->quantity_name_, this->element_tag_);
    this->mean_variance_xml_engine_out_.addElementToXmlDoc("Variance_Element");
    SimTK::Xml::Element variance_element_ = this->mean_variance_xml_engine_out_.getChildElement("Variance_Element");
    this->xmlmemory_io_.writeDataToXmlMemory(this->mean_variance_xml_engine_out_, variance_element_, this->variance_new_,
                                             SMIN(this->snapshot_, this->number_of_snapshot_old_), this->observation_, this->quantity_name_, this->element_tag_);
    this->mean_variance_xml_engine_out_.writeToXmlFile(this->mean_variance_filefullpath_);
}
//=================================================================================================//
template <class ObserveMethodType>
bool RegressionTestEnsembleAverage<ObserveMethodType>::compareMeanVariance()
{
    int count_not_converged_m = 0;
    int count_not_converged_v = 0;
    count_not_converged_m = compareParameter("meanvalue", meanvalue_, meanvalue_new_, this->threshold_mean_);
    count_not_converged_v = compareParameter("variance", variance_, variance_new_, this->threshold_variance_);
    if (count_not_converged_m == 0)
    {
        std::cout << "The meanvalue of " << this->quantity_name_ << " are converged now." << std::endl;
        if (count_not_converged_v == 0)
        {
            if (this->label_for_repeat_ == 4)
            {
                this->converged_ = "true";
                this->label_for_repeat_++;
                std::cout << "The meanvalue and variance of " << this->quantity_name_ << " are converged enough times, and run will stop now." << std::endl;
                return true;
            }
            else
            {
                this->converged_ = "false";
                this->label_for_repeat_++;
                std::cout << "The variance of " << this->quantity_name_ << " are also converged, and this is the " << this->label_for_repeat_
                          << " times. They should be converged more times to be stable." << std::endl;
                return false;
            }
        }
        else
        {
            this->converged_ = "false";
            this->label_for_repeat_ = 0;
            std::cout << "The variance of " << this->quantity_name_ << " are not converged " << count_not_converged_v << " times." << std::endl;
            return false;
        };
    }
    else
    {
        this->converged_ = "false";
        this->label_for_repeat_ = 0;
        std::cout << "The meanvalue of " << this->quantity_name_ << " are not converged " << count_not_converged_m << " times." << std::endl;
        return false;
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestEnsembleAverage<ObserveMethodType>::resultTest()
{
    /* compare the current result to the converged mean value and variance. */
    int test_wrong = 0;
    if (this->snapshot_ < this->number_of_snapshot_old_)
        test_wrong = testNewResult(this->snapshot_difference_, this->current_result_, meanvalue_, variance_);
    else
    {
        /** Unify the length of meanvalue, variance, old result, new result. */
        for (int delete_ = 0; delete_ != this->snapshot_difference_; ++delete_)
        {
            meanvalue_.pop_back();
            variance_.pop_back();
        }
        test_wrong = testNewResult(0, this->current_result_, meanvalue_, variance_);
    }
    /* draw the conclusion. */
    if (test_wrong == 0)
        std::cout << "The result of " << this->quantity_name_ << " are correct based on the ensemble averaged regression test!" << std::endl;
    else
    {
        std::cout << "There are " << test_wrong << " snapshots are not within the expected range." << std::endl;
        std::cout << "Please try again. If it still post this conclusion, the result is not correct!" << std::endl;
        exit(1);
    }
}
//=================================================================================================//
} // namespace SPH