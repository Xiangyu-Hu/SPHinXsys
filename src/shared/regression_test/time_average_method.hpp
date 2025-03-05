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
void RegressionTestTimeAverage<ObserveMethodType>::filterLocalResult(BiVector<VariableType> &current_result)
{
    int scale = round(this->snapshot_ / 200);
    std::cout << "The filter scale is " << scale * 2 << "." << std::endl;
    for (int snapshot_index = 0; snapshot_index != this->snapshot_; ++snapshot_index)
    {
        for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        {
            VariableType filter_meanvalue = ZeroData<VariableType>::value;
            VariableType filter_variance = ZeroData<VariableType>::value;
            for (int index = SMAX(snapshot_index - scale, 0); index != SMIN(snapshot_index + scale, this->snapshot_); ++index)
            {
                filter_meanvalue += current_result[index][observation_index];
            }
            filter_meanvalue = (filter_meanvalue - current_result[snapshot_index][observation_index]) /
                               (SMIN(snapshot_index + scale, this->snapshot_) - SMAX(snapshot_index - scale, 0));
            for (int index = SMAX(snapshot_index - scale, 0); index != SMIN(snapshot_index + scale, this->snapshot_); ++index)
            {
                VariableType filter_deviation = current_result[index][observation_index] - filter_meanvalue;
                filter_variance += componentSquare(filter_deviation);
            }

            VariableType current_deviation = current_result[snapshot_index][observation_index] - filter_meanvalue;
            VariableType current_variance = componentSquare(current_deviation);
            filter_variance = (filter_variance - current_variance) /
                              (SMIN(snapshot_index + scale, this->snapshot_) - SMAX(snapshot_index - scale, 0));

            current_result[snapshot_index][observation_index] = transformComponent(
                current_result[snapshot_index][observation_index],
                [&](Real cur_result, Real cur_variance, Real flt_variance, Real flt_meanvalue)
                { 
                    if(cur_variance > 4 * flt_variance)
                    {
                        std::cout << "A component of the current " << this->quantity_name_
                                  << "[" << snapshot_index << "][" << observation_index << "] is "
                                  << cur_result << ", but the neighbor averaged value is " << flt_meanvalue
                                  << ", and the rate is " << cur_variance / flt_variance << std::endl;
                        return flt_meanvalue;
                    }
                    return cur_result; },
                current_variance, filter_variance, filter_meanvalue);
        }
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestTimeAverage<ObserveMethodType>::searchSteadyStart()
{
    /* the search is for each value within parameters. */
    int scale = round(this->snapshot_ / 20);
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
        for (int snapshot_index = this->snapshot_ - 1; snapshot_index != 3 * scale; --snapshot_index)
        {
            Real value_one = 0, value_two = 0;
            for (int index = snapshot_index; index != snapshot_index - scale; --index)
            {
                value_one += FirstComponent(this->current_result_[index][observation_index]) / scale;
                value_two += FirstComponent(this->current_result_[index - 2 * scale][observation_index]) / scale;
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
void RegressionTestTimeAverage<ObserveMethodType>::calculateNewVariance(BiVector<VariableType> &current_result_trans)
{
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
        {
            VariableType deviation = current_result_trans[observation_index][snapshot_index] - local_meanvalue_[observation_index];
            variance_new_[observation_index] += componentSquare(deviation);
        }

        variance_new_[observation_index] = transformComponent(
            variance_new_[observation_index], [&](Real variance_new, Real variance, Real local_meanvalue)
            { return SMAX(
                  (Real)(variance_new / (this->snapshot_ - snapshot_for_converged_)),
                  (Real)variance,
                  (Real)pow(local_meanvalue * Real(0.01), 2)); },
            variance_[observation_index], local_meanvalue_[observation_index]);
    }
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestTimeAverage<ObserveMethodType>::compareParameter(
    std::string par_name, StdVec<VariableType> &parameter, StdVec<VariableType> &parameter_new, VariableType &threshold)
{
    int count = 0;
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        transformComponent(
            parameter[observation_index], [&](Real para, Real para_new, Real thr)
            { 
                if ((par_name == "meanvalue") && (ABS(para) < 0.005) && (ABS(para_new) < 0.005))
                {
                    std::cout << "A component of the old meanvalue is " << para << ", and the new meanvalue is " << para_new
                              << ". So this component will be ignored due to its tiny effect." << std::endl;
                }
                Real rel_value = ABS((para - para_new) / (para_new + TinyReal));
                if (rel_value > thr)
                {
                    std::cout << par_name << ": " << this->quantity_name_ << "[" << observation_index << "]"
                              << " has a component is not converged, and difference is " << rel_value << std::endl;
                    count++;
                }
                return rel_value; },
            parameter_new[observation_index], threshold);
    }
    return count;
}
//=================================================================================================//
template <class ObserveMethodType>
int RegressionTestTimeAverage<ObserveMethodType>::testNewResult(
    BiVector<VariableType> &current_result,
    StdVec<VariableType> &meanvalue, StdVec<VariableType> &local_meanvalue, StdVec<VariableType> &variance)
{
    int count = 0;
    for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
    {
        for (int snapshot_index = snapshot_for_converged_; snapshot_index != this->snapshot_; ++snapshot_index)
        {
            VariableType deviation = current_result[snapshot_index][observation_index] - local_meanvalue[observation_index];
            variance_new_[observation_index] += componentSquare(deviation);
        }
        variance_new_[observation_index] = variance_new_[observation_index] / (this->snapshot_ - snapshot_for_converged_);

        transformComponent(
            meanvalue[observation_index], [&](Real mean, Real local_mean, Real var, Real var_new)
            { 
                if ((ABS(mean) < 0.005) && (ABS(local_mean) < 0.005))
                {
                    std::cout << "A component of the old meanvalue is " << mean << ", and the new meanvalue is " << local_mean
                              << ". So this component will be ignored due to its tiny effect." << std::endl;
                }
                Real rel_value = ABS((mean - local_mean) / (mean + TinyReal));
                if (rel_value > 0.1 || var_new > (1.01 * var))
                {
                    std::cout << this->quantity_name_ << "[" << observation_index << "]"
                              << " has a component beyond the exception !" << std::endl;
                    std::cout << "The meanvalue component is " << mean << ", and the current meanvalue is " << local_mean << std::endl;
                    std::cout << "The variance component is " << var << ", and the current variance is " << var_new << std::endl;
                    count++;
                }
                return rel_value; },
            local_meanvalue[observation_index], variance[observation_index], variance_new_[observation_index]);
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
    searchSteadyStart();
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
    calculateNewVariance(this->current_result_trans_);
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