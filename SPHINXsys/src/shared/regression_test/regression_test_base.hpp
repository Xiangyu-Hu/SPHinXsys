/**
 * @file 	regression_test_base.cpp
 * @author	Bo Zhang and Xiangyu Hu
 */

#include "regression_test_base.h"

 //=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestBase<ObserveMethodType>::transferTheIndex()
	{
		int number_of_index_i = this->current_result_.size();
		int number_of_index_j = this->current_result_[0].size();
		DoubleVec<VariableType> temp(number_of_index_j, StdVec<VariableType>(number_of_index_i));
		current_result_ji_ = temp;
		for (int index_i = 0; index_i != number_of_index_i; ++index_i)
			for (int index_j = 0; index_j != number_of_index_j; ++index_j)
				current_result_ji_[index_j][index_i] = this->current_result_[index_i][index_j];
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestBase<ObserveMethodType>::readResultFromXml()
	{
		if (number_of_run_ > 1)
		{
			DoubleVec<VariableType> result_in_(SMAX(i_, number_of_snapshot_old_), StdVec<VariableType>(j_));
			for (int run_n_ = 0; run_n_ != number_of_run_ - 1; ++run_n_)
			{
				std::string node_name_ = "Round_" + std::to_string(run_n_);
				SimTK::Xml::Element father_element_ = result_xml_engine_in_.getChildElement(node_name_);
				for (int particle_n_ = 0; particle_n_ != j_; ++particle_n_)
					xmlmemory_io_.readDataFromXmlMemory(result_xml_engine_in_, father_element_, particle_n_, result_in_, this->quantity_name_);
				DoubleVec<VariableType> result_temp_ = result_in_;
				for (int delete_ = 0; delete_ != difference_; ++delete_)
					result_temp_.pop_back(); /* Unify the length of all results. (number of snapshots)*/
				result_.push_back(result_temp_);
			}
			result_.push_back(this->current_result_);
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestBase<ObserveMethodType>::writeResultToXml()
	{
		for (int run_n_ = 0; run_n_ != number_of_run_; ++run_n_)
		{
			std::string node_name_ = "Round_" + std::to_string(run_n_);
			result_xml_engine_out_.addElementToXmlDoc(node_name_);
			SimTK::Xml::Element father_element_ =
				result_xml_engine_out_.getChildElement(node_name_);
			xmlmemory_io_.writeDataToXmlMemory(result_xml_engine_out_, father_element_, result_[run_n_],
				SMIN(i_, number_of_snapshot_old_), j_, this->quantity_name_, this->element_tag_);
		}
		result_xml_engine_out_.writeToXmlFile(result_filefullpath_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestBase<ObserveMethodType>::readResultFromXml(int index_of_run_)
	{
		if (number_of_run_ > 1)
		{
			result_filefullpath_ = input_folder_path_ + "/" + this->body_name_ + "_" + this->quantity_name_ +
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

			result_xml_engine_in_.loadXmlFile(result_filefullpath_);
			SimTK::Xml::Element snapshot_element_ = result_xml_engine_in_.getChildElement("Snapshot_Element");
			SimTK::Xml::element_iterator ele_ite = snapshot_element_.element_begin();
			result_xml_engine_in_.getRequiredAttributeValue(ele_ite, "number_of_snapshot_for_local_result_", i_);

			DoubleVec<VariableType> result_temp_(j_, StdVec<VariableType>(i_));
			result_in_ = result_temp_;
			SimTK::Xml::Element result_element_ = result_xml_engine_in_.getChildElement("Result_Element");
			for (int snapshot_n_ = 0; snapshot_n_ != i_; ++snapshot_n_)
			{
				int index_i_ = 0;
				SimTK::Xml::element_iterator ele_ite = result_element_.element_begin();
				for (; ele_ite != result_element_.element_end(); ++ele_ite)
				{
					std::string attribute_name_ = "snapshot_" + std::to_string(snapshot_n_);
					result_xml_engine_in_.getRequiredAttributeValue(ele_ite, attribute_name_, result_in_[index_i_][snapshot_n_]);
					index_i_++;
				}
			};
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestBase<ObserveMethodType>::writeResultToXml(int index_of_run_)
	{
		/** write result to .xml (with different data structure to Base). */
		int total_snapshot_ = current_result_ji_[0].size();
		result_filefullpath_ = input_folder_path_ + "/" + this->body_name_ + "_" + this->quantity_name_ +
			"_Run_" + std::to_string(index_of_run_) + "_result.xml";
		result_xml_engine_out_.addElementToXmlDoc("Snapshot_Element");
		SimTK::Xml::Element snapshot_element_ = result_xml_engine_out_.getChildElement("Snapshot_Element");
		result_xml_engine_out_.addChildToElement(snapshot_element_, "Snapshot");
		SimTK::Xml::element_iterator ele_ite = snapshot_element_.element_begin();
		result_xml_engine_out_.setAttributeToElement(ele_ite, "number_of_snapshot_for_local_result_", total_snapshot_);

		result_xml_engine_out_.addElementToXmlDoc("Result_Element");
		SimTK::Xml::Element result_element_ = result_xml_engine_out_.getChildElement("Result_Element");
		for (int particle_n_ = 0; particle_n_ != j_; ++particle_n_)
		{
			std::string element_name_ = "Particle_" + std::to_string(particle_n_);
			result_xml_engine_out_.addChildToElement(result_element_, element_name_);
			for (int snapshot_n_ = 0; snapshot_n_ != total_snapshot_; ++snapshot_n_)
			{
				SimTK::Xml::element_iterator ele_ite = result_element_.element_begin(element_name_);
				std::string attribute_name_ = "snapshot_" + std::to_string(snapshot_n_);
				result_xml_engine_out_.setAttributeToElement(ele_ite, attribute_name_, current_result_ji_[particle_n_][snapshot_n_]);
			}
		}
		result_xml_engine_out_.writeToXmlFile(result_filefullpath_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
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
}