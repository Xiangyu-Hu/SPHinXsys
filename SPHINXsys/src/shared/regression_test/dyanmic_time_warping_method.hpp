/**
 * @file 	dynamic_time_warping_method.hpp
 * @author	Bo Zhang and Xiangyu Hu
 */

#pragma once

#include "dynamic_time_warping_method.h"

 //=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	template<class ObserveMethodType>
	Real RegressionTestDynamicTimeWarping<ObserveMethodType>::calculatePNorm(Real variable_a, Real variable_b)
	{
		return std::abs(variable_a - variable_b);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	Real RegressionTestDynamicTimeWarping<ObserveMethodType>::calculatePNorm(Vecd variable_a, Vecd variable_b)
	{
		Real distance = 0;
		for (int dimension_index = 0; dimension_index < variable_a.size(); ++dimension_index)
			distance = std::pow(std::abs(variable_a[dimension_index] - variable_b[dimension_index]), 2);
		return std::pow(distance, 0.5);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	Real RegressionTestDynamicTimeWarping<ObserveMethodType>::calculatePNorm(Matd variable_a, Matd variable_b)
	{
		/** the current method is TBD. */
		Real distance = 0;
		for (int dimension_index_i = 0; dimension_index_i != variable_a.size(); ++dimension_index_i)
			for (int dimension_index_j = 0; dimension_index_j != variable_a.size(); ++dimension_index_j)
				distance = std::pow(std::abs(variable_a[dimension_index_i][dimension_index_j] - variable_b[dimension_index_i][dimension_index_j]), 2);
		return std::pow(distance, 0.5);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	StdVec<Real> RegressionTestDynamicTimeWarping<ObserveMethodType>::calculateDTWDistance
	(DoubleVec<VariableType> dataset_a_, DoubleVec<VariableType> dataset_b_)
	{
		/* define the container to hold the dtw distance.*/
		StdVec<Real> dtw_distance;
		for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
		{
			int window_size_ = 5;
			int a_length = dataset_a_[observation_index].size();
			int b_length = dataset_b_[observation_index].size();
			/** create a 2D vector with [a_length, b_length] to contain the local DTW distance value. */
			DoubleVec<Real> local_dtw_distance(a_length, StdVec<Real>(b_length, 0));
			local_dtw_distance[0][0] = calculatePNorm(dataset_a_[observation_index][0], dataset_b_[observation_index][0]);
			for (int index_i = 1; index_i < a_length; ++index_i)
				local_dtw_distance[index_i][0] = local_dtw_distance[index_i - 1][0] + calculatePNorm(dataset_a_[observation_index][index_i], dataset_b_[observation_index][0]);
			for (int index_j = 1; index_j < b_length; ++index_j)
				local_dtw_distance[0][index_j] = local_dtw_distance[0][index_j - 1] + calculatePNorm(dataset_a_[observation_index][0], dataset_b_[observation_index][index_j]);

			/** add locality constraint */
			window_size_ = SMAX(window_size_, ABS(a_length - b_length));
			for (int index_i = 1; index_i != a_length; ++index_i)
				for (int index_j = SMAX(1, index_i - window_size_); index_j != SMIN(b_length, index_i + window_size_); ++index_j)
					local_dtw_distance[index_i][index_j] = calculatePNorm(dataset_a_[observation_index][index_i], dataset_b_[observation_index][index_j]) +
						SMIN(local_dtw_distance[index_i - 1][index_j], local_dtw_distance[index_i][index_j - 1], local_dtw_distance[index_i - 1][index_j - 1]);
			dtw_distance.push_back(local_dtw_distance[a_length - 1][b_length - 1]);
		}
		return dtw_distance;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestDynamicTimeWarping<ObserveMethodType>::setupTheTest()
	{
		this->snapshot_ = this->current_result_.size();
		this->observation_ = this->current_result_[0].size();

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
	template<class ObserveMethodType>
	void RegressionTestDynamicTimeWarping<ObserveMethodType>::readDTWDistanceFromXml()
	{
		if (this->number_of_run_ > 1)
		{
			dtw_distance_xml_engine_in_.loadXmlFile(dtw_distance_filefullpath_);
			SimTK::Xml::Element element_name_dtw_distance_ = dtw_distance_xml_engine_in_.root_element_;
			SimTK::Xml::element_iterator ele_ite = element_name_dtw_distance_.element_begin();
			for (; ele_ite != element_name_dtw_distance_.element_end(); ++ele_ite)
				for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
				{
					std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(observation_index);
					dtw_distance_xml_engine_in_.getRequiredAttributeValue<Real>(ele_ite, attribute_name_, dtw_distance_[observation_index]);
				}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestDynamicTimeWarping<ObserveMethodType>::updateDTWDistance()
	{
		if (this->number_of_run_ > 1)
		{
			StdVec<Real> dtw_distance_local_(this->observation_, 0);
			dtw_distance_local_ = calculateDTWDistance(this->current_result_trans_, this->result_in_);
			for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
			{
				dtw_distance_new_[observation_index] = SMAX(dtw_distance_local_[observation_index], dtw_distance_[observation_index], dtw_distance_new_[observation_index]);
			}	
			this->result_in_.clear();
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestDynamicTimeWarping<ObserveMethodType>::writeDTWDistanceToXml()
	{
		SimTK::Xml::Element DTWElement = dtw_distance_xml_engine_out_.root_element_;
		dtw_distance_xml_engine_out_.addChildToElement(DTWElement, "DTWDistance");
		SimTK::Xml::element_iterator ele_ite = DTWElement.element_begin();
		for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
		{
			std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(observation_index);
			dtw_distance_xml_engine_out_.setAttributeToElement(ele_ite, attribute_name_, dtw_distance_new_[observation_index]);
		}
		dtw_distance_xml_engine_out_.writeToXmlFile(dtw_distance_filefullpath_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	bool RegressionTestDynamicTimeWarping<ObserveMethodType>::compareDTWDistance(Real threshold_value)
	{
		if (this->number_of_run_ > 1)
		{
			int count_not_converged_ = 0;
			for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
			{
				if (std::abs(dtw_distance_[observation_index] - dtw_distance_new_[observation_index]) > threshold_value)
				{
					count_not_converged_++;
					std::cout << "The DTW distance of " << this->quantity_name_ << " [" << observation_index << "] is not converged." << endl;
					std::cout << "The old DTW distance is " << dtw_distance_[observation_index] << ", and the new DTW distance is " << dtw_distance_new_[observation_index] << "." << endl;
				}
			};

			if (count_not_converged_ == 0)
			{
				if (this->label_for_repeat_ == 4)
				{
					this->converged = "true";
					std::cout << "The DTW distance of " << this->quantity_name_ << " are converged enough times, and rum will stop now." << endl;
					return true;
				}
				else
				{
					this->converged = "false";
					this->label_for_repeat_++;
					std::cout << "The DTW distance of " << this->quantity_name_ << " are converged, and this is the " << this->label_for_repeat_
						<< " times. They should be converged more times to be stable." << endl;
				}
			}
			else if (count_not_converged_ != 0)
			{
				this->converged = "false";
				this->label_for_repeat_ = 0;
				return false;
			};
		}
		return true;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestDynamicTimeWarping<ObserveMethodType>::resultTest()
	{
		int test_wrong = 0;
		StdVec<Real> dtw_distance_current_;
		dtw_distance_current_ = calculateDTWDistance(this->result_in_, this->current_result_trans_);
		for (int observation_index = 0; observation_index != this->observation_; ++observation_index)
		{
			if (dtw_distance_current_[observation_index] > 1.01 * dtw_distance_[observation_index])
			{
				std::cout << "The maximum distance of " << this->quantity_name_ << "[" << observation_index << "] is " << dtw_distance_[observation_index] 
					<< ", and the current distance is " << dtw_distance_current_[observation_index] << "." << endl;
				test_wrong++;
			}
		};
		if (test_wrong == 0)
		{
			std::cout << "The DTW distance of " << this->quantity_name_ << " between current result and this previous local result is within exception!" << endl;
		}
		else
		{
			std::cout << "The DTW distance of " << this->quantity_name_ << " between current result and this previous local result is beyond exception!" << endl;
			std::cout << "Please try again. If it still post this sentence, the result is not correct!" << endl;
			exit(1);
		}
	};
	//=================================================================================================/
}