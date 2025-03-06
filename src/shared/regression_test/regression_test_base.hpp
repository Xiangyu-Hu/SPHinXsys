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
      result_xml_engine_in_("result_xml_engine_in", "result"),
      result_xml_engine_out_("result_xml_engine_out", "result"),
      snapshot_(0), observation_(this->NumberOfObservedQuantity()),
      number_of_snapshot_old_(0), snapshot_difference_(0),
      converged_("false"), number_of_run_(1), label_for_repeat_(0)

{
    input_folder_path_ = this->io_environment_.input_folder_;
    in_output_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + ".xml";
    result_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + "_result.xml";
    runtimes_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ + "_runtimes.dat";

    if (!fs::exists(runtimes_filefullpath_) && !generate_regression_data_)
    {
        std::cout << "Error: " << runtimes_filefullpath_ << " is missing for regression test! " << std::endl;
        exit(1);
    }
    else
    {
        std::ifstream in_file(runtimes_filefullpath_.c_str());
        in_file >> converged_;
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
    int snapshot, int observation, const std::string &quantity_name, StdVec<std::string> &element_tag)
{
    for (int l = 0; l != snapshot; ++l)
    {
        std::string element_name = element_tag[l];
        xml_engine.addChildToElement(element, element_name);
        for (int k = 0; k != observation; ++k)
        {
            SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
            std::string attribute_name = quantity_name + "_" + std::to_string(k);
            xml_engine.setAttributeToElement(ele_ite, attribute_name, quantity[l][k]);
        }
    }
}
//=================================================================================================//
template <class ObserveMethodType>
template <typename T>
void RegressionTestBase<ObserveMethodType>::writeDataToXmlMemory(
    XmlEngine &xml_engine, SimTK::Xml::Element &element,
    std::string element_name, int k, const T &quantity, const std::string &quantity_name)
{
    SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name);
    std::string attribute_name = quantity_name + "_" + std::to_string(k);
    xml_engine.setAttributeToElement(ele_ite, attribute_name, quantity);
}
//=================================================================================================//
template <class ObserveMethodType>
template <typename T>
void RegressionTestBase<ObserveMethodType>::readDataFromXmlMemory(
    XmlEngine &xml_engine, SimTK::Xml::Element &element,
    int k, BiVector<T> &result_container, const std::string &quantity_name)
{
    int l = 0;
    SimTK::Xml::element_iterator ele_ite = element.element_begin();
    for (; ele_ite != element.element_end(); ++ele_ite)
    {
        std::string attribute_name = quantity_name + "_" + std::to_string(k);
        xml_engine.getRequiredAttributeValue(ele_ite, attribute_name, result_container[l][k]);
        l++;
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::readTagFromXmlMemory(
    SimTK::Xml::Element &element, StdVec<std::string> &element_tag)
{
    size_t l = 0;
    SimTK::Xml::element_iterator ele_ite = element.element_begin();
    for (; ele_ite != element.element_end(); ++ele_ite)
    {
        element_tag[l] = ele_ite->getElementTag();
        l++;
    }
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::rememberObservations(size_t iteration)
{
    VariableType *observed_quantities = this->getObservedQuantity();
    StdVec<VariableType> observed;
    for (int i = 0; i != observation_; ++i)
    {
        observed.push_back(observed_quantities[i]);
    };
    current_result_.push_back(observed);
    element_tag_.push_back("Snapshot_" + std::to_string(iteration));
    snapshot_++;
}
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::transposeTheIndex()
{
    current_result_trans_ = BiVector<VariableType>(observation_, StdVec<VariableType>(snapshot_));
    for (int l = 0; l != snapshot_; ++l)
        for (int k = 0; k != observation_; ++k)
            current_result_trans_[k][l] = this->current_result_[l][k];
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::readResultFromXml()
{
    if (number_of_run_ > 1) /*only read the result from the 2nd run, because the 1st run doesn't have previous results. */
    {
        /*Here result_temp is a temporary container that reloads each previous result.*/
        BiVector<VariableType> result_temp(SMAX(snapshot_, number_of_snapshot_old_), StdVec<VariableType>(observation_));
        for (int run_index = 0; run_index != number_of_run_ - 1; ++run_index)
        {
            std::string node_name = "Round_" + std::to_string(run_index);
            SimTK::Xml::Element father_element = result_xml_engine_in_.getChildElement(node_name);
            for (int k = 0; k != observation_; ++k)
                readDataFromXmlMemory(result_xml_engine_in_, father_element, k, result_temp, this->quantity_name_);
            for (int k = 0; k != snapshot_difference_; ++k)
                result_temp.pop_back(); /* trim the new reading result to unify the length of all results. (number of snapshots) */
            result_.push_back(result_temp);
        }
        result_.push_back(this->current_result_); /* Finally, push back the current result into the result vector. */
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::writeResultToXml()
{
    for (int run_index = 0; run_index != number_of_run_; ++run_index)
    {
        std::string node_name = "Round_" + std::to_string(run_index);
        result_xml_engine_out_.addElementToXmlDoc(node_name);
        SimTK::Xml::Element father_element =
            result_xml_engine_out_.getChildElement(node_name);
        writeDataToXmlMemory(
            result_xml_engine_out_, father_element, result_[run_index],
            SMIN(snapshot_, number_of_snapshot_old_), observation_, this->quantity_name_, this->element_tag_);
    }
    result_xml_engine_out_.writeToXmlFile(result_filefullpath_);
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::readResultFromXml(int run_index)
{
    if (number_of_run_ > 1) /*only read the result from the 2nd run, because the 1st run doesn't have previous results. */
    {
        result_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ +
                               "_Run_" + std::to_string(run_index) + "_result.xml";

        /* To identify the database generation or new result test. */
        if (converged_ == "false")
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
        SimTK::Xml::Element snapshot_element = result_xml_engine_in_.getChildElement("Snapshot_Element");
        SimTK::Xml::element_iterator ele_ite = snapshot_element.element_begin();
        result_xml_engine_in_.getRequiredAttributeValue(ele_ite, "number_of_snapshot_for_local_result_", snapshot_);

        result_in_ = BiVector<VariableType>(observation_, StdVec<VariableType>(snapshot_));
        SimTK::Xml::Element result_element = result_xml_engine_in_.getChildElement("Result_Element");
        for (int l = 0; l != snapshot_; ++l)
        {
            int k = 0;
            SimTK::Xml::element_iterator ele_ite = result_element.element_begin();
            for (; ele_ite != result_element.element_end(); ++ele_ite)
            {
                std::string attribute_name = "snapshot_" + std::to_string(l);
                result_xml_engine_in_.getRequiredAttributeValue(ele_ite, attribute_name, result_in_[k][l]);
                k++;
            }
        };
    }
};
//=================================================================================================//
template <class ObserveMethodType>
void RegressionTestBase<ObserveMethodType>::writeResultToXml(int run_index)
{
    /** write result to .xml (with different data structure to Base), here is
        observation * snapshot, which can be used for TA and DTW methods. */
    result_filefullpath_ = input_folder_path_ + "/" + this->dynamics_identifier_name_ + "_" + this->quantity_name_ +
                           "_Run_" + std::to_string(run_index) + "_result.xml";
    result_xml_engine_out_.addElementToXmlDoc("Snapshot_Element");
    SimTK::Xml::Element snapshot_element = result_xml_engine_out_.getChildElement("Snapshot_Element");
    result_xml_engine_out_.addChildToElement(snapshot_element, "Snapshot");
    SimTK::Xml::element_iterator ele_ite = snapshot_element.element_begin();
    result_xml_engine_out_.setAttributeToElement(ele_ite, "number_of_snapshot_for_local_result_", snapshot_);

    result_xml_engine_out_.addElementToXmlDoc("Result_Element");
    SimTK::Xml::Element result_element = result_xml_engine_out_.getChildElement("Result_Element");
    for (int k = 0; k != observation_; ++k)
    {
        std::string element_name = "Particle_" + std::to_string(k);
        result_xml_engine_out_.addChildToElement(result_element, element_name);
        for (int l = 0; l != snapshot_; ++l)
        {
            SimTK::Xml::element_iterator ele_ite = result_element.element_begin(element_name);
            std::string attribute_name = "snapshot_" + std::to_string(l);
            result_xml_engine_out_.setAttributeToElement(ele_ite, attribute_name, current_result_trans_[k][l]);
        }
    }
    result_xml_engine_out_.writeToXmlFile(result_filefullpath_);
};
//=================================================================================================//
template <class ObserveMethodType>
RegressionTestBase<ObserveMethodType>::~RegressionTestBase()
{
    if (converged_ == "false")
    {
        number_of_run_ += 1;
    }
    if (fs::exists(runtimes_filefullpath_))
    {
        fs::remove(runtimes_filefullpath_);
    }
    std::ofstream out_file(runtimes_filefullpath_.c_str());
    out_file << converged_;
    out_file << "\n";
    out_file << number_of_run_;
    out_file << "\n";
    out_file << label_for_repeat_;
    out_file.close();
};
//=================================================================================================//
} // namespace SPH