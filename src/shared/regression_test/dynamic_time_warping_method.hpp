/**
 * @file 	dynamic_time_warping_method.hpp
 * @brief 	Classes for the comparison between validated and tested results
                with dynamic time warping method.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "dynamic_time_warping_method.h"

namespace SPH
{
//=================================================================================================//
template <class ObserveMethodType>
template <typename... Args>
RegressionTestDynamicTimeWarping<ObserveMethodType>::RegressionTestDynamicTimeWarping(Args &&...args)
    : RegressionTestTimeAverage<ObserveMethodType>(std::forward<Args>(args)...),
      dtw_distance_xml_parser_("dtw_distance_xml_engine_out", "dtw_distance")
{
    dtw_distance_filefullpath_ = this->input_folder_path_ + "/" + this->dynamics_identifier_name_ +
                                 "_" + this->quantity_name_ + "_dtwdistance.xml";
}
//=================================================================================================//
template <class ObserveMethodType>
StdVec<Real> RegressionTestDynamicTimeWarping<ObserveMethodType>::calculateDTWDistance(
    BiVector<VariableType> dataset_a_, BiVector<VariableType> dataset_b_)
{
    /* define the container to hold the dtw distance.*/
    StdVec<Real> dtw_distance;
    for (int observe_k = 0; observe_k != this->observe_; ++observe_k)
    {
        int window_size = 5;
        int a_length = dataset_a_[observe_k].size();
        int b_length = dataset_b_[observe_k].size();

        if (b_length > 1.1 * a_length || b_length < 0.9 * a_length)
        {
            std::cout << "\n Error: please check the time step change, because the data length changed a lot !" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            exit(1);
        }

        /** create a 2D vector with [a_length, b_length] to contain the local DTW distance value. */
        BiVector<Real> local_distance(a_length, StdVec<Real>(b_length, 0));
        local_distance[0][0] = getFirstNorm(dataset_a_[observe_k][0] - dataset_b_[observe_k][0]);
        for (int i = 1; i < a_length; ++i)
            local_distance[i][0] = local_distance[i - 1][0] +
                                   getFirstNorm(dataset_a_[observe_k][i] - dataset_b_[observe_k][0]);
        for (int j = 1; j < b_length; ++j)
            local_distance[0][j] = local_distance[0][j - 1] +
                                   getFirstNorm(dataset_a_[observe_k][0] - dataset_b_[observe_k][j]);

        /** add locality constraint */
        window_size = SMAX(window_size, ABS(a_length - b_length));
        for (int i = 1; i != a_length; ++i)
            for (int j = SMAX(1, i - window_size); j != SMIN(b_length, i + window_size); ++j)
                local_distance[i][j] = getFirstNorm(dataset_a_[observe_k][i] - dataset_b_[observe_k][j]) +
                                       SMIN(local_distance[i - 1][j], local_distance[i][j - 1], local_distance[i - 1][j - 1]);
        dtw_distance.push_back(local_distance[a_length - 1][b_length - 1]);
    }
    return dtw_distance;
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::setupTest()
{
    this->snapshot_ = this->current_result_.size();
    this->observe_ = this->current_result_[0].size();

    StdVec<Real> dtw_distance_temp_(this->observe_, 0);
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
        dtw_distance_xml_parser_.loadXmlFile(dtw_distance_filefullpath_);
        auto first_element = dtw_distance_xml_parser_.first_element_;
        auto child_element = dtw_distance_xml_parser_.findElement(first_element, "DTWDistance");
        for (int observe_k = 0; observe_k != this->observe_; ++observe_k)
        {
            std::string attribute_name = this->quantity_name_ + "_" + std::to_string(observe_k);
            dtw_distance_xml_parser_.queryAttributeValue(child_element, attribute_name, dtw_distance_[observe_k]);
        }
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::updateDTWDistance()
{
    if (this->number_of_run_ > 1)
    {
        StdVec<Real> local_distance(this->observe_, 0);
        local_distance = calculateDTWDistance(this->current_result_trans_, this->result_in_);
        for (int observe_k = 0; observe_k != this->observe_; ++observe_k)
        {
            dtw_distance_new_[observe_k] =
                SMAX(local_distance[observe_k], dtw_distance_[observe_k], dtw_distance_new_[observe_k]);
        }
        this->result_in_.clear();
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::writeDTWDistanceToXml()
{
    auto first_element = dtw_distance_xml_parser_.first_element_;
    auto child_element = dtw_distance_xml_parser_.addNewElement(first_element, "DTWDistance");
    for (int observe_k = 0; observe_k != this->observe_; ++observe_k)
    {
        std::string attribute_name = this->quantity_name_ + "_" + std::to_string(observe_k);
        dtw_distance_xml_parser_.setAttributeToElement(child_element, attribute_name, dtw_distance_new_[observe_k]);
    }
    dtw_distance_xml_parser_.writeToXmlFile(dtw_distance_filefullpath_);
};
//=================================================================================================//
template <class ObserveMethodType>
bool RegressionTestDynamicTimeWarping<ObserveMethodType>::compareDTWDistance(Real threshold_value)
{
    if (this->number_of_run_ > 1)
    {
        int count_not_converged = 0;
        for (int observe_k = 0; observe_k != this->observe_; ++observe_k)
        {
            if (std::abs(dtw_distance_[observe_k] - dtw_distance_new_[observe_k]) > threshold_value)
            {
                count_not_converged++;
                std::cout << "The DTW distance of " << this->quantity_name_
                          << " [" << observe_k << "] is not converged." << std::endl;
                std::cout << "The old DTW distance is " << dtw_distance_[observe_k]
                          << ", and the new DTW distance is " << dtw_distance_new_[observe_k] << "." << std::endl;
            }
        };

        if (count_not_converged == 0)
        {
            if (this->label_for_repeat_ == 4)
            {
                this->converged = "true";
                std::cout << "The DTW distance of " << this->quantity_name_
                          << " are converged enough times, and rum will stop now." << std::endl;
                return true;
            }
            else
            {
                this->converged = "false";
                this->label_for_repeat_++;
                std::cout << "The DTW distance of " << this->quantity_name_
                          << " are converged, and this is the " << this->label_for_repeat_
                          << " times. They should be converged more times to be stable." << std::endl;
            }
        }
        else if (count_not_converged != 0)
        {
            this->converged = "false";
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
    StdVec<Real> dtw_distance_current;
    dtw_distance_current = calculateDTWDistance(this->result_in_, this->current_result_trans_);
    for (int observe_k = 0; observe_k != this->observe_; ++observe_k)
    {
        if (dtw_distance_current[observe_k] > 1.01 * dtw_distance_[observe_k])
        {
            std::cout << "The maximum distance of " << this->quantity_name_ << "[" << observe_k << "] is "
                      << dtw_distance_[observe_k] << ", and the current distance is "
                      << dtw_distance_current[observe_k] << "." << std::endl;
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
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::generateDataBase(
    Real threshold_value, const std::string &filter)
{
    this->writeXmlToXmlFile();
    this->readXmlFromXmlFile();
    this->transposeTheIndex();
    if (this->converged == "false")
    {
        setupTest();
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
        compareDTWDistance(threshold_value); // wether the distance is convergence.
    }
    else
        std::cout << "The results have been converged." << std::endl;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestDynamicTimeWarping<ObserveMethodType>::testResult(const std::string &filter)
{
    this->writeXmlToXmlFile();
    this->readXmlFromXmlFile();
    this->transposeTheIndex();
    setupTest();
    if (filter == "true")
        this->filterExtremeValues();
    readDTWDistanceFromXml();
    for (int n = 0; n != this->number_of_run_; ++n)
    {
        this->result_filefullpath_ = this->input_folder_path_ + "/" + this->dynamics_identifier_name_ +
                                     "_" + this->quantity_name_ + "_Run_" + std::to_string(n) + "_result.xml";
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
}
//=================================================================================================/
} // namespace SPH