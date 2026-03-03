/**
 * @file 	dynamic_time_warping_method.hpp
 * @brief 	Classes for the comparison between validated and tested results
                with dynamic time warping method.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "dynamic_time_warping_method.h"

//=================================================================================================//
namespace SPH
{
//=================================================================================================//
template <class ObserveMethodType>
StdVec<Real> RegressionTestDynamicTimeWarping<ObserveMethodType>::
    calculateDTWDistance(BiVector<VariableType> dataset_a_, BiVector<VariableType> dataset_b_)
{
    /* define the container to hold the dtw distance.*/
    StdVec<Real> dtw_distance;
    for (int k = 0; k != this->observation_; ++k)
    {
        int window_size_ = 5;
        int a_length = dataset_a_[k].size();
        int b_length = dataset_b_[k].size();

        if (b_length > 1.2 * a_length || b_length < 0.8 * a_length)
        {
            std::cout << "\n Error: please check the time step change, because the data length changed a lot !" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        /** create a 2D vector with [a_length, b_length] to contain the local DTW distance value. */
        BiVector<Real> local_distance(a_length, StdVec<Real>(b_length, 0));
        local_distance[0][0] = calculatePNorm(dataset_a_[k][0], dataset_b_[k][0]);
        for (int i = 1; i < a_length; ++i)
            local_distance[i][0] = local_distance[i - 1][0] +
                                       calculatePNorm(dataset_a_[k][i], dataset_b_[k][0]);
        for (int j = 1; j < b_length; ++j)
            local_distance[0][j] = local_distance[0][j - 1] +
                                       calculatePNorm(dataset_a_[k][0], dataset_b_[k][j]);

        /** add locality constraint */
        window_size_ = SMAX(window_size_, ABS(a_length - b_length));
        for (int i = 1; i != a_length; ++i)
            for (int j = SMAX(1, i - window_size_); j != SMIN(b_length, i + window_size_); ++j)
                local_distance[i][j] =
                    calculatePNorm(dataset_a_[k][i], dataset_b_[k][j]) +
                    SMIN(local_distance[i - 1][j], local_distance[i][j - 1], local_distance[i - 1][j - 1]);
        dtw_distance.push_back(local_distance[a_length - 1][b_length - 1]);
    }
    return dtw_distance;
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::setupTheTest()
{
    StdVec<Real> dtw_distance_temp_(this->observation_, 0);
    dtw_distance_ = dtw_distance_temp_;
    dtw_distance_new_ = dtw_distance_;

    if ((this->number_of_run_ > 1) && (!fs::exists(dtw_distance_filefullpath_)))
    {
        std::cout << "\n Error: the input file:" << dtw_distance_filefullpath_ << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::readDTWDistanceFromXml()
{
    if (this->number_of_run_ > 1)
    {
        dtw_distance_xml_engine_in_.loadXmlFile(dtw_distance_filefullpath_);
        SimTK::Xml::Element element_name_dtw_distance = dtw_distance_xml_engine_in_.root_element_;
        SimTK::Xml::element_iterator ele_ite = element_name_dtw_distance.element_begin();
        for (; ele_ite != element_name_dtw_distance.element_end(); ++ele_ite)
            for (int k = 0; k != this->observation_; ++k)
            {
                std::string attribute_name = this->quantity_name_ + "_" + std::to_string(k);
                dtw_distance_xml_engine_in_.getRequiredAttributeValue(ele_ite, attribute_name, dtw_distance_[k]);
            }
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::updateDTWDistance()
{
    if (this->number_of_run_ > 1)
    {
        StdVec<Real> dtw_distance_local_(this->observation_, 0);
        dtw_distance_local_ = calculateDTWDistance(this->current_result_trans_, this->result_in_);
        for (int k = 0; k != this->observation_; ++k)
        {
            dtw_distance_new_[k] = SMAX(dtw_distance_local_[k], dtw_distance_[k], dtw_distance_new_[k]);
        }
        this->result_in_.clear();
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::writeDTWDistanceToXml()
{
    SimTK::Xml::Element DTWElement = dtw_distance_xml_engine_out_.root_element_;
    dtw_distance_xml_engine_out_.addChildToElement(DTWElement, "DTWDistance");
    SimTK::Xml::element_iterator ele_ite = DTWElement.element_begin();
    for (int k = 0; k != this->observation_; ++k)
    {
        std::string attribute_name = this->quantity_name_ + "_" + std::to_string(k);
        dtw_distance_xml_engine_out_.setAttributeToElement(ele_ite, attribute_name, dtw_distance_new_[k]);
    }
    dtw_distance_xml_engine_out_.writeToXmlFile(dtw_distance_filefullpath_);
};
//=================================================================================================//
template <class ObserveMethodType>
bool RegressionTestDynamicTimeWarping<ObserveMethodType>::compareDTWDistance(Real threshold_value)
{
    if (this->number_of_run_ > 1)
    {
        int count_not_converged_ = 0;
        for (int k = 0; k != this->observation_; ++k)
        {
            if (std::abs(dtw_distance_[k] - dtw_distance_new_[k]) > threshold_value)
            {
                count_not_converged_++;
                std::cout << "The DTW distance of " << this->quantity_name_ << " [" << k << "] is not converged." << std::endl;
                std::cout << "The old DTW distance is " << dtw_distance_[k] << ", and the new DTW distance is "
                          << dtw_distance_new_[k] << "." << std::endl;
            }
        };

        if (count_not_converged_ == 0)
        {
            if (this->label_for_repeat_ == 4)
            {
                this->converged_ = "true";
                std::cout << "The DTW distance of " << this->quantity_name_
                          << " are converged enough times, and run will stop now." << std::endl;
                return true;
            }
            else
            {
                this->converged_ = "false";
                this->label_for_repeat_++;
                std::cout << "The DTW distance of " << this->quantity_name_
                          << " are converged, and this is the " << this->label_for_repeat_
                          << " times. They should be converged more times to be stable." << std::endl;
            }
        }
        else if (count_not_converged_ != 0)
        {
            this->converged_ = "false";
            this->label_for_repeat_ = 0;
            return false;
        };
    }
    return true;
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::resultTest()
{
    int test_wrong = 0;
    StdVec<Real> dtw_distance_current_;
    dtw_distance_current_ = calculateDTWDistance(this->result_in_, this->current_result_trans_);
    for (int k = 0; k != this->observation_; ++k)
    {
        if (dtw_distance_current_[k] > 1.01 * dtw_distance_[k])
        {
            std::cout << "The maximum distance of " << this->quantity_name_ << "[" << k << "] is " << dtw_distance_[k]
                      << ", and the current distance is " << dtw_distance_current_[k] << "." << std::endl;
            test_wrong++;
        }
    };
    if (test_wrong == 0)
    {
        std::cout << "The DTW distance of " << this->quantity_name_
                  << " between current result and this previous local result is within exception!" << std::endl;
    }
    else
    {
        std::cout << "The DTW distance of " << this->quantity_name_
                  << " between current result and this previous local result is beyond exception!" << std::endl;
        std::cout << "Please try again. If it still post this sentence, the result is not correct!" << std::endl;
        exit(1);
    }
};
//=================================================================================================/
} // namespace SPH