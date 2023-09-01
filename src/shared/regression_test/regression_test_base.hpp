/**
 * @file regression_test_base.hpp
 * @brief Base classes for comparisons between validated and tested results.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "regression_test_base.h"

//=================================================================================================//
namespace SPH
{
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::writeToXml(ObservedQuantityRecording<VariableType> *observe_method, size_t iteration)
{
    this->exec();
    std::string element_name_ = "Snapshot_" + std::to_string(iteration);
    SimTK::Xml::Element &element_ = observe_xml_engine_.root_element_;
    observe_xml_engine_.addElementToXmlDoc(element_name_);
    for (size_t i = 0; i != this->base_particles_.total_real_particles_; ++i)
    {
        xmlmemory_io_.writeDataToXmlMemory(observe_xml_engine_, element_,
                                           element_name_, i, (*this->interpolated_quantities_)[i], this->quantity_name_);
    };
};
//=================================================================================================//
template <class ObserveMethodType>
template <typename ReduceType>
void RegressionTestBase<ObserveMethodType>::writeToXml(ReducedQuantityRecording<ReduceType> *reduce_method, size_t iteration)
{
    std::string element_name_ = "Snapshot_" + std::to_string(iteration);
    SimTK::Xml::Element &element_ = observe_xml_engine_.root_element_;
    observe_xml_engine_.addElementToXmlDoc(element_name_);
    xmlmemory_io_.writeDataToXmlMemory(observe_xml_engine_, element_,
                                       element_name_, 0, this->reduce_method_.exec(), this->quantity_name_);
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::readFromXml(ObservedQuantityRecording<VariableType> *observe_method)
{
    observe_xml_engine_.loadXmlFile(in_output_filefullpath_);
    size_t number_of_particle_ = this->base_particles_.total_real_particles_;
    size_t number_of_snapshot_ = std::distance(observe_xml_engine_.root_element_.element_begin(),
                                               observe_xml_engine_.root_element_.element_end());
    BiVector<VariableType> current_result_temp_(number_of_snapshot_, StdVec<VariableType>(number_of_particle_));
    StdVec<std::string> element_tag_temp_(number_of_snapshot_);
    current_result_ = current_result_temp_;
    element_tag_ = element_tag_temp_;
    SimTK::Xml::Element &element_ = observe_xml_engine_.root_element_;
    for (size_t j = 0; j != number_of_particle_; ++j)
    {
        xmlmemory_io_.readDataFromXmlMemory(observe_xml_engine_, element_, j, current_result_, this->quantity_name_);
        xmlmemory_io_.readTagFromXmlMemory(element_, element_tag_);
    }
};
//=================================================================================================//
template <class ObserveMethodType>
template <typename ReduceType>
void RegressionTestBase<ObserveMethodType>::readFromXml(ReducedQuantityRecording<ReduceType> *reduce_method)
{
    observe_xml_engine_.loadXmlFile(in_output_filefullpath_);
    size_t number_of_particle_ = 1;
    size_t number_of_snapshot_ = std::distance(observe_xml_engine_.root_element_.element_begin(),
                                               observe_xml_engine_.root_element_.element_end());
    BiVector<VariableType> current_result_temp_(number_of_snapshot_, StdVec<VariableType>(number_of_particle_));
    StdVec<std::string> element_tag_temp_(number_of_snapshot_);
    current_result_ = current_result_temp_;
    element_tag_ = element_tag_temp_;
    SimTK::Xml::Element &element_ = observe_xml_engine_.root_element_;
    for (size_t j = 0; j != number_of_particle_; ++j)
    {
        xmlmemory_io_.readDataFromXmlMemory(observe_xml_engine_, element_, j, current_result_, this->quantity_name_);
        xmlmemory_io_.readTagFromXmlMemory(element_, element_tag_);
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::transposeTheIndex()
{
    int number_of_snapshot = this->current_result_.size();
    int number_of_observation = this->current_result_[0].size();
    BiVector<VariableType> temp(number_of_observation, StdVec<VariableType>(number_of_snapshot));
    current_result_trans_ = temp;
    for (int snapshot_index = 0; snapshot_index != number_of_snapshot; ++snapshot_index)
        for (int observation_index = 0; observation_index != number_of_observation; ++observation_index)
            current_result_trans_[observation_index][snapshot_index] = this->current_result_[snapshot_index][observation_index];
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::readResultFromXml()
{
    if (number_of_run_ > 1) /*only read the result from the 2nd run, because the 1st run doesn't have previous results. */
    {
        /*Here result_in is a temporary container that reloads each previous result.*/
        BiVector<VariableType> result_in_(SMAX(snapshot_, number_of_snapshot_old_), StdVec<VariableType>(observation_));
        for (int run_index_ = 0; run_index_ != number_of_run_ - 1; ++run_index_)
        {
            std::string node_name_ = "Round_" + std::to_string(run_index_);
            SimTK::Xml::Element father_element_ = result_xml_engine_in_.getChildElement(node_name_);
            for (int observation_index_ = 0; observation_index_ != observation_; ++observation_index_)
                xmlmemory_io_.readDataFromXmlMemory(result_xml_engine_in_, father_element_, observation_index_, result_in_, this->quantity_name_);
            BiVector<VariableType> result_temp_ = result_in_;
            for (int delete_ = 0; delete_ != difference_; ++delete_)
                result_temp_.pop_back(); /* trim the new reading result to unify the length of all results. (number of snapshots) */
            result_.push_back(result_temp_);
        }
        result_.push_back(this->current_result_); /* Finally, push back the current result into the result vector. */
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::writeResultToXml()
{
    for (int run_index_ = 0; run_index_ != number_of_run_; ++run_index_)
    {
        std::string node_name_ = "Round_" + std::to_string(run_index_);
        result_xml_engine_out_.addElementToXmlDoc(node_name_);
        SimTK::Xml::Element father_element_ =
            result_xml_engine_out_.getChildElement(node_name_);
        xmlmemory_io_.writeDataToXmlMemory(result_xml_engine_out_, father_element_, result_[run_index_],
                                           SMIN(snapshot_, number_of_snapshot_old_), observation_, this->quantity_name_, this->element_tag_);
    }
    result_xml_engine_out_.writeToXmlFile(result_filefullpath_);
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::readResultFromXml(int index_of_run_)
{
    if (number_of_run_ > 1) /*only read the result from the 2nd run, because the 1st run doesn't have previous results. */
    {
        result_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ +
                               "_Run_" + std::to_string(index_of_run_) + "_result.xml";

        /* To identify the database generation or new result test. */
        if (converged == "false")
        {
            if (!fs::exists(result_filefullpath_))
            {
                std::cout << "\n Error: the input file:" << result_filefullpath_ << " is not exists" << std::endl;
                std::cout << __FILE__ << ':' << __LINE__ << std::endl;
                exit(1);
            }
        }

        /*Each result has two elements, one records the length of this result, and the other one is itself.*/
        result_xml_engine_in_.loadXmlFile(result_filefullpath_);
        SimTK::Xml::Element snapshot_element_ = result_xml_engine_in_.getChildElement("Snapshot_Element");
        SimTK::Xml::element_iterator ele_ite = snapshot_element_.element_begin();
        result_xml_engine_in_.getRequiredAttributeValue(ele_ite, "number_of_snapshot_for_local_result_", snapshot_);

        BiVector<VariableType> result_temp_(observation_, StdVec<VariableType>(snapshot_));
        result_in_ = result_temp_;
        SimTK::Xml::Element result_element_ = result_xml_engine_in_.getChildElement("Result_Element");
        for (int snapshot_index = 0; snapshot_index != snapshot_; ++snapshot_index)
        {
            int observation_index = 0;
            SimTK::Xml::element_iterator ele_ite = result_element_.element_begin();
            for (; ele_ite != result_element_.element_end(); ++ele_ite)
            {
                std::string attribute_name_ = "snapshot_" + std::to_string(snapshot_index);
                result_xml_engine_in_.getRequiredAttributeValue(ele_ite, attribute_name_, result_in_[observation_index][snapshot_index]);
                observation_index++;
            }
        };
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::writeResultToXml(int index_of_run_)
{
    /** write result to .xml (with different data structure to Base), here is
        observation * snapshot, which can be used for TA and DTW methods. */
    int total_snapshot_ = current_result_trans_[0].size();
    result_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ +
                           "_Run_" + std::to_string(index_of_run_) + "_result.xml";
    result_xml_engine_out_.addElementToXmlDoc("Snapshot_Element");
    SimTK::Xml::Element snapshot_element_ = result_xml_engine_out_.getChildElement("Snapshot_Element");
    result_xml_engine_out_.addChildToElement(snapshot_element_, "Snapshot");
    SimTK::Xml::element_iterator ele_ite = snapshot_element_.element_begin();
    result_xml_engine_out_.setAttributeToElement(ele_ite, "number_of_snapshot_for_local_result_", total_snapshot_);

    result_xml_engine_out_.addElementToXmlDoc("Result_Element");
    SimTK::Xml::Element result_element_ = result_xml_engine_out_.getChildElement("Result_Element");
    for (int observation_index = 0; observation_index != observation_; ++observation_index)
    {
        std::string element_name_ = "Particle_" + std::to_string(observation_index);
        result_xml_engine_out_.addChildToElement(result_element_, element_name_);
        for (int snapshot_index = 0; snapshot_index != total_snapshot_; ++snapshot_index)
        {
            SimTK::Xml::element_iterator ele_ite = result_element_.element_begin(element_name_);
            std::string attribute_name_ = "snapshot_" + std::to_string(snapshot_index);
            result_xml_engine_out_.setAttributeToElement(ele_ite, attribute_name_, current_result_trans_[observation_index][snapshot_index]);
        }
    }
    result_xml_engine_out_.writeToXmlFile(result_filefullpath_);
};
//=================================================================================================//
template <class ObserveMethodType>
RegressionTestBase<ObserveMethodType>::~RegressionTestBase()
{
    if (converged == "false")
    {
        number_of_run_ += 1;
    }
    if (fs::exists(runtimes_filefullpath_))
    {
        fs::remove(runtimes_filefullpath_);
    }
    std::ofstream out_file(runtimes_filefullpath_.c_str());
    out_file << converged;
    out_file << "\n";
    out_file << number_of_run_;
    out_file << "\n";
    out_file << label_for_repeat_;
    out_file.close();
};
//=================================================================================================//
} // namespace SPH