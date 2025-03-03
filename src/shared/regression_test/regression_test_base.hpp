/**
 * @file regression_test_base.hpp
 * @brief Base classes for comparisons between validated and tested results.
 * @author	Bo Zhang , Chi Zhang and Xiangyu Hu
 */

#pragma once

#include "regression_test_base.h"

namespace SPH
{
//=================================================================================================//
template <class ObserveMethodType>
template <typename... Args>
RegressionTestBase<ObserveMethodType>::RegressionTestBase(Args &&...args)
    : ObserveMethodType(std::forward<Args>(args)...),
      generate_regression_data_(this->sph_system_.GenerateRegressionData()),
      observe_xml_engine_("xml_observe_reduce", this->quantity_name_),
      result_xml_engine_in_("result_xml_engine_in", "result"),
      result_xml_engine_out_("result_xml_engine_out", "result")
{
    input_folder_path_ = this->io_environment_.input_folder_;
    in_output_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + ".xml";
    result_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + "_result.xml";
    runtimes_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + "_runtimes.dat";

    if (!fs::exists(runtimes_filefullpath_))
    {
        if (generate_regression_data_)
        {
            converged = "false";
            number_of_run_ = 1;
            label_for_repeat_ = 0;
        }
        else
        {
            std::cout << "Error: runtimes file is missing for regression test! " << std::endl;
            exit(1);
        }
    }
    else
    {
        std::ifstream in_file(runtimes_filefullpath_.c_str());
        in_file >> converged;
        in_file >> number_of_run_;
        in_file >> label_for_repeat_;
        in_file.close();
    };
}
//=================================================================================================//
template <class ObserveMethodType>
template <typename T>
void RegressionTestBase<ObserveMethodType>::writeDataToXmlMemory(
    XmlEngine &xml_engine, SimTK::Xml::Element &element, const BiVector<T> &quantity,
    int snapshot_, int observe_, const std::string &quantity_name, StdVec<std::string> &element_tag)
{
    for (int snapshot_index = 0; snapshot_index != snapshot_; ++snapshot_index)
    {
        std::string element_name = element_tag[snapshot_index];
        auto ele_ite = xml_engine.addChildToElement(element, element_name);
        for (int observe_k = 0; observe_k != observe_; ++observe_k)
        {
            std::string attribute_name_ = quantity_name + "_" + std::to_string(observe_k);
            xml_engine.setAttributeToElement(ele_ite, attribute_name_, quantity[snapshot_index][observe_k]);
        }
    }
}
//=================================================================================================//
template <class ObserveMethodType>
template <typename T>
void RegressionTestBase<ObserveMethodType>::writeDataToXmlMemory(
    XmlEngine &xml_engine, SimTK::Xml::element_iterator ele_ite,
    int observe_k, const T &quantity, const std::string &quantity_name)
{
    std::string attribute_name_ = quantity_name + "_" + std::to_string(observe_k);
    xml_engine.setAttributeToElement(ele_ite, attribute_name_, quantity);
}
//=================================================================================================//
template <class ObserveMethodType>
template <typename T>
void RegressionTestBase<ObserveMethodType>::readDataFromXmlMemory(
    XmlEngine &xml_engine, SimTK::Xml::Element &element,
    int observe_k, BiVector<T> &result_container, const std::string &quantity_name)
{
    int snapshot_index = 0;
    for (auto ele_ite = element.element_begin(); ele_ite != element.element_end(); ++ele_ite)
    {
        std::string attribute_name_ = quantity_name + "_" + std::to_string(observe_k);
        xml_engine.getRequiredAttributeValue(ele_ite, attribute_name_, result_container[snapshot_index][observe_k]);
        snapshot_index++;
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::readTagFromXmlMemory(
    SimTK::Xml::Element &element, StdVec<std::string> &element_tag)
{
    size_t snapshot_index = 0;
    for (auto ele_ite = element.element_begin(); ele_ite != element.element_end(); ++ele_ite)
    {
        element_tag[snapshot_index] = ele_ite->getElementTag();
        snapshot_index++;
    }
}
//=================================================================================================//
template <class ObserveMethodType>
template <typename... Parameters>
void RegressionTestBase<ObserveMethodType>::
    writeToXml(ObservedQuantityRecording<Parameters...> *observe_method, size_t iteration)
{
    this->exec();
    std::string element_name_ = "Snapshot_" + std::to_string(iteration);
    auto ele_ite = observe_xml_engine_.addElementToXmlDoc(element_name_);
    VariableType *interpolated_quantities = this->dv_interpolated_quantities_->Data();
    for (size_t i = 0; i != this->base_particles_.TotalRealParticles(); ++i)
    {
        writeDataToXmlMemory(observe_xml_engine_, ele_ite,
                             i, interpolated_quantities[i], this->quantity_name_);
    };
};
//=================================================================================================//
template <class ObserveMethodType>
template <typename... Parameters>
void RegressionTestBase<ObserveMethodType>::
    writeToXml(ReducedQuantityRecording<Parameters...> *reduce_method, size_t iteration)
{
    std::string element_name_ = "Snapshot_" + std::to_string(iteration);
    auto ele_ite = observe_xml_engine_.addElementToXmlDoc(element_name_);
    writeDataToXmlMemory(observe_xml_engine_, ele_ite,
                         0, this->reduce_method_.exec(), this->quantity_name_);
};
//=================================================================================================//
template <class ObserveMethodType>
template <typename... Parameters>
void RegressionTestBase<ObserveMethodType>::
    readFromXml(ObservedQuantityRecording<Parameters...> *observe_method)
{
    observe_xml_engine_.loadXmlFile(in_output_filefullpath_);
    size_t number_of_particle_ = this->base_particles_.TotalRealParticles();
    size_t number_of_snapshot_ = std::distance(observe_xml_engine_.root_element_.element_begin(),
                                               observe_xml_engine_.root_element_.element_end());
    BiVector<VariableType> current_result_temp_(number_of_snapshot_, StdVec<VariableType>(number_of_particle_));
    StdVec<std::string> element_tag_temp_(number_of_snapshot_);
    current_result_ = current_result_temp_;
    element_tag_ = element_tag_temp_;
    SimTK::Xml::Element &element_ = observe_xml_engine_.root_element_;
    for (size_t j = 0; j != number_of_particle_; ++j)
    {
        readDataFromXmlMemory(observe_xml_engine_, element_, j, current_result_, this->quantity_name_);
        readTagFromXmlMemory(element_, element_tag_);
    }
};
//=================================================================================================//
template <class ObserveMethodType>
template <typename... Parameters>
void RegressionTestBase<ObserveMethodType>::
    readFromXml(ReducedQuantityRecording<Parameters...> *reduce_method)
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
        readDataFromXmlMemory(observe_xml_engine_, element_, j, current_result_, this->quantity_name_);
        readTagFromXmlMemory(element_, element_tag_);
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
        for (int observe_k = 0; observe_k != number_of_observation; ++observe_k)
            current_result_trans_[observe_k][snapshot_index] = this->current_result_[snapshot_index][observe_k];
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::readResultFromXml()
{
    if (number_of_run_ > 1) /*only read the result from the 2nd run, because the 1st run doesn't have previous results. */
    {
        /*Here result_in is a temporary container that reloads each previous result.*/
        BiVector<VariableType> result_in_(SMAX(snapshot_, number_of_snapshot_old_), StdVec<VariableType>(observe_));
        for (int run_index_ = 0; run_index_ != number_of_run_ - 1; ++run_index_)
        {
            std::string node_name_ = "Round_" + std::to_string(run_index_);
            SimTK::Xml::Element father_element_ = result_xml_engine_in_.getChildElement(node_name_);
            for (int observation_index_ = 0; observation_index_ != observe_; ++observation_index_)
                readDataFromXmlMemory(result_xml_engine_in_, father_element_, observation_index_, result_in_, this->quantity_name_);
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
        writeDataToXmlMemory(result_xml_engine_out_, father_element_, result_[run_index_],
                             SMIN(snapshot_, number_of_snapshot_old_), observe_, this->quantity_name_, this->element_tag_);
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

        BiVector<VariableType> result_temp_(observe_, StdVec<VariableType>(snapshot_));
        result_in_ = result_temp_;
        SimTK::Xml::Element result_element_ = result_xml_engine_in_.getChildElement("Result_Element");
        for (int snapshot_index = 0; snapshot_index != snapshot_; ++snapshot_index)
        {
            int observe_k = 0;
            for (auto ele_ite = result_element_.element_begin(); ele_ite != result_element_.element_end(); ++ele_ite)
            {
                std::string attribute_name_ = "snapshot_" + std::to_string(snapshot_index);
                result_xml_engine_in_.getRequiredAttributeValue(ele_ite, attribute_name_, result_in_[observe_k][snapshot_index]);
                observe_k++;
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
    auto snapshot_element_ = result_xml_engine_out_.addElementToXmlDoc("Snapshot_Element");
    auto ele_ite = result_xml_engine_out_.addChildToElement(*snapshot_element_, "Snapshot");
    result_xml_engine_out_.setAttributeToElement(ele_ite, "number_of_snapshot_for_local_result_", total_snapshot_);

    auto result_element_ = result_xml_engine_out_.addElementToXmlDoc("Result_Element");
    for (int observe_k = 0; observe_k != observe_; ++observe_k)
    {
        std::string element_name_ = "Particle_" + std::to_string(observe_k);
        auto ele_ite = result_xml_engine_out_.addChildToElement(*result_element_, element_name_);
        for (int snapshot_index = 0; snapshot_index != total_snapshot_; ++snapshot_index)
        {
            std::string attribute_name_ = "snapshot_" + std::to_string(snapshot_index);
            result_xml_engine_out_.setAttributeToElement(ele_ite, attribute_name_, current_result_trans_[observe_k][snapshot_index]);
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
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::ensureGenerateRegressionData()
{
    if (!generate_regression_data_)
    {
        std::cout << "Error: run test when command option is generate regression data! " << std::endl;
        exit(1);
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::ensureTestResult()
{
    if (generate_regression_data_)
    {
        std::cout << "Error: generate regression data when command option is run test! " << std::endl;
        exit(1);
    }
}
//=================================================================================================//
} // namespace SPH