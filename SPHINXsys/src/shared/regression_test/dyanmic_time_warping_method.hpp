/**
 * @file 	dynamic_time_warping_method.hpp
 * @author	Bo Zhang and Xiangyu Hu
 */

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
		for (int i = 0; i < variable_a.size(); ++i)
			distance = std::pow(std::abs(variable_a[i] - variable_b[i]), 2);
		return std::pow(distance, 0.5);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	Real RegressionTestDynamicTimeWarping<ObserveMethodType>::calculatePNorm(Matd variable_a, Matd variable_b)
	{
		/** the current method is TBD. */
		Real distance = 0;
		for (int i = 0; i != variable_a.size(); ++i)
			for (int j = 0; j != variable_a.size(); ++j)
				distance = std::pow(std::abs(variable_a[i][j] - variable_b[i][j]), 2);
		return std::pow(distance, 0.5);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	StdVec<Real> RegressionTestDynamicTimeWarping<ObserveMethodType>::calculateDTWDistance
	(DoubleVec<VariableType> dataset_a_, DoubleVec<VariableType> dataset_b_)
	{
		/* define the container to hold the dtw distance.*/
		StdVec<Real> dtw_distance;
		for (int n = 0; n != this->j_; ++n)
		{
			int window_size_ = 5;
			int a_length = dataset_a_[n].size();
			int b_length = dataset_b_[n].size();
			/** create a 2D vector with [a_length, b_length] to contain the local DTW distance value. */
			DoubleVec<Real> local_dtw_distance(a_length, StdVec<Real>(b_length, 0));
			local_dtw_distance[0][0] = calculatePNorm(dataset_a_[n][0], dataset_b_[n][0]);
			for (int i = 1; i < a_length; ++i)
				local_dtw_distance[i][0] = local_dtw_distance[i - 1][0] + calculatePNorm(dataset_a_[n][i], dataset_b_[n][0]);
			for (int j = 1; j < b_length; ++j)
				local_dtw_distance[0][j] = local_dtw_distance[0][j - 1] + calculatePNorm(dataset_a_[n][0], dataset_b_[n][j]);

			/** add locality constraint */
			window_size_ = SMAX(window_size_, ABS(a_length - b_length));
			for (int i = 1; i != a_length; ++i)
				for (int j = SMAX(1, i - window_size_); j != SMIN(b_length, i + window_size_); ++j)
					local_dtw_distance[i][j] = calculatePNorm(dataset_a_[n][i], dataset_b_[n][j]) +
						SMIN(local_dtw_distance[i - 1][j], local_dtw_distance[i][j - 1], local_dtw_distance[i - 1][j - 1]);
			dtw_distance.push_back(local_dtw_distance[a_length - 1][b_length - 1]);
		}
		return dtw_distance;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestDynamicTimeWarping<ObserveMethodType>::settingupTheTest()
	{
		this->i_ = this->current_result_.size();
		this->j_ = this->current_result_[0].size();

		StdVec<Real> dtw_distance_temp_(this->j_, 0);
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
				for (int particle_n_ = 0; particle_n_ != this->j_; ++particle_n_)
				{
					std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
					dtw_distance_xml_engine_in_.getRequiredAttributeValue<Real>(ele_ite, attribute_name_, dtw_distance_[particle_n_]);
				}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestDynamicTimeWarping<ObserveMethodType>::updateDTWDistance()
	{
		if (this->number_of_run_ > 1)
		{
			StdVec<Real> dtw_distance_local_(this->j_, 0);
			dtw_distance_local_ = calculateDTWDistance(this->current_result_ji_, this->result_in_);
			for (int j = 0; j != this->j_; ++j)
			{
				dtw_distance_new_[j] = SMAX(dtw_distance_local_[j], dtw_distance_[j], dtw_distance_new_[j]);
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
		for (int particle_n_ = 0; particle_n_ != this->j_; ++particle_n_)
		{
			std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
			dtw_distance_xml_engine_out_.setAttributeToElement(ele_ite, attribute_name_, dtw_distance_new_[particle_n_]);
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
			for (int j = 0; j != this->j_; ++j)
			{
				if (std::abs(dtw_distance_[j] - dtw_distance_new_[j]) > threshold_value)
				{
					count_not_converged_++;
					std::cout << "The DTW distance of " << this->quantity_name_ << " [" << j << "] is not converged." << endl;
					std::cout << "The old DTW distance is " << dtw_distance_[j] << ", and the new DTW distance is " << dtw_distance_new_[j] << "." << endl;
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
		dtw_distance_current_ = calculateDTWDistance(this->result_in_, this->current_result_ji_);
		for (int i = 0; i != this->j_; ++i)
		{
			if (dtw_distance_current_[i] > 1.01 * dtw_distance_[i])
			{
				std::cout << "The maximum distance of " << this->quantity_name_ << "[" << i << "] is " << dtw_distance_[i] 
					<< ", and the current distance is " << dtw_distance_current_[i] << "." << endl;
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