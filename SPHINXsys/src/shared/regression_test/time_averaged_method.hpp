/**
 * @file 	time_averaged_method.cpp
 * @author	Bo Zhang and Xiangyu Hu
 */

#include "time_averaged_method.h"

 //=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::filterLocalResult(DoubleVec<Real> &current_result)
	{
		int scale = round(this->i_ / 200);
		std::cout << "The filter scale is " << scale * 2 << "." << endl;
		for (int i = 0; i != this->i_; ++i)
		{
			for (int j = 0; j != this->j_; ++j)
			{
				Real filter_meanvalue = 0;
				Real filter_variance = 0;
				for (int index = SMAX(i - scale, 0); index != SMIN(i + scale, this->i_); ++index)
				{
					filter_meanvalue += current_result[index][j];
				}
				filter_meanvalue = (filter_meanvalue - current_result[i][j]) / (SMIN(i + scale, this->i_) - SMAX(i - scale, 0));
				for (int index = SMAX(i - scale, 0); index != SMIN(i + scale, this->i_); ++index)
				{
					filter_variance += std::pow(current_result[index][j] - filter_meanvalue, 2);
				}
				Real current_variance = std::pow(current_result[i][j] - filter_meanvalue, 2);
				filter_variance = (filter_variance - current_variance) / (SMIN(i + scale, this->i_) - SMAX(i - scale, 0));
				if (current_variance > 4 * filter_variance)
				{
					current_result[i][j] = filter_meanvalue;
					std::cout << "The current value of " << this->quantity_name_ << "[" << i << "][" << j << "] is " << current_result[i][j]
						<< ", but the neighbor averaged value is " << filter_meanvalue << ", and the rate is " << current_variance / filter_variance << endl;
				}
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::filterLocalResult(DoubleVec<Vecd> &current_result)
	{
		int scale = round(this->i_ / 200);
		std::cout << "The filter scale is " << scale * 2 << "." << endl;
		for (int i = 0; i != this->i_; ++i)
		{
			for (int j = 0; j != this->j_; ++j)
			{
				for (int index_i = 0; index_i != current_result[0][0].size(); ++index_i)
				{
					Real filter_meanvalue = 0;
					Real filter_variance = 0;
					for (int index = SMAX(i - scale, 0); index != SMIN(i + scale, this->i_); ++index)
					{
						filter_meanvalue += current_result[index][j][index_i];
					}
					filter_meanvalue = (filter_meanvalue - current_result[i][j][index_i]) / (SMIN(i + scale, this->i_) - SMAX(i - scale, 0));
					for (int index = SMAX(i - scale, 0); index != SMIN(i + scale, this->i_); ++index)
					{
						filter_variance += std::pow(current_result[index][j][index_i] - filter_meanvalue, 2);
					}
					Real current_variance = std::pow(current_result[i][j][index_i] - filter_meanvalue, 2);
					filter_variance = (filter_variance - current_variance) / (SMIN(i + scale, this->i_) - SMAX(i - scale, 0));
					if (current_variance > 4 * filter_variance)
					{
						current_result[i][j][index_i] = filter_meanvalue;
						std::cout << "The current value of " << this->quantity_name_ << "[" << i << "][" << j << "][" << index_i << "] is " << current_result[i][j][index_i]
							<< ", but the neighbor averaged value is " << filter_meanvalue << ", and the rate is " << current_variance / filter_variance << endl;
					}
				}
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::filterLocalResult(DoubleVec<Matd> &current_result)
	{
		int scale = round(this->i_ / 200);
		std::cout << "The filter scale is " << scale * 2 << "." << endl;
		for (int i = 0; i != this->i_; ++i)
		{
			for (int j = 0; j != this->j_; ++j)
			{
				for (int index_i = 0; index_i != current_result[0][0].size(); ++index_i)
				{
					for (int index_j = 0; index_j != current_result[0][0].size(); ++index_j)
					{
						Real filter_meanvalue = 0;
						Real filter_variance = 0;
						for (int index = SMAX(i - scale, 0); index != SMIN(i + scale, this->i_); ++index)
						{
							filter_meanvalue += current_result[index][j][index_i][index_j];
						}
						filter_meanvalue = (filter_meanvalue - current_result[i][j][index_i][index_j]) / (SMIN(i + scale, this->i_) - SMAX(i - scale, 0));
						for (int index = SMAX(i - scale, 0); index != SMIN(i + scale, this->i_); ++index)
						{
							filter_variance += std::pow(current_result[index][j][index_i][index_j] - filter_meanvalue, 2);
						}
						Real current_variance = std::pow(current_result[i][j][index_i][index_j] - filter_meanvalue, 2);
						filter_variance = (filter_variance - current_variance) / (SMIN(i + scale, this->i_) - SMAX(i - scale, 0));
						if (current_variance > 4 * filter_variance)
						{
							current_result[i][j][index_i][index_j] = filter_meanvalue;
							std::cout << "The current value of " << this->quantity_name_ << "[" << i << "][" << j << "][" << index_i << "][" << index_j << "] is " << current_result[i][j][index_i][index_j]
								<< ", but the neighbor averaged value is " << filter_meanvalue << ", and the rate is " << current_variance / filter_variance << endl;
						}
					}
				}
			}	
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::searchSteadyStart(DoubleVec<Real> &current_result)
	{
		int scale = round(this->i_ / 20);
		for (int j = 0; j != this->j_; ++j)
			for (int i = this->i_ - 1; i != 3 * scale; --i)
			{
				Real value_one = 0, value_two = 0;
				for (int index = i; index != i - scale; --index)
				{
					value_one += current_result[index][j];
					value_two += current_result[index - 2 * scale][j];
				}
				value_one = value_one / scale;
				value_two = value_two / scale;
				if (ABS(value_one - value_two) / ABS((value_one + value_two) / 2) > 0.1)
				{
					snapshot_for_converged_ = SMAX(snapshot_for_converged_, i - scale);
					break;
				}
			}
		std::cout << "The scale is " << scale << "." << endl;
	};
	//=================================================================================================// 
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::searchSteadyStart(DoubleVec<Vecd> &current_result)
	{
		int scale = round(this->i_ / 20);
		for (int j = 0; j != this->j_; ++j)
			for (int i = this->i_ - 1; i != 3 * scale; --i)
			{
				Real value_one = 0, value_two = 0;
				for (int index = i; index != i - scale; --index)
				{
					value_one += current_result[index][j][0];
					value_two += current_result[index - 2 * scale][j][0];
				}
				value_one = value_one / scale;
				value_two = value_two / scale;
				if (ABS(value_one - value_two) / ABS((value_one + value_two) / 2) > 0.1)
				{
					snapshot_for_converged_ = SMAX(snapshot_for_converged_, i - scale);
					break;
				}
			}
		std::cout << "The scale is " << scale << "." << endl;
	};
	//=================================================================================================// 
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::searchSteadyStart(DoubleVec<Matd> &current_result)
	{
		int scale = round(this->i_ / 20);
		for (int j = 0; j != this->j_; ++j)
			for (int i = this->i_ - 1; i != 3 * scale; --i)
			{
				Real value_one = 0, value_two = 0;
				for (int index = i; index != i - 20; --index)
				{
					value_one += current_result[index][j][0][0];
					value_two += current_result[index - 2 * scale][j][0][0];
				}
				value_one = value_one / scale;
				value_two = value_two / scale;
				if (ABS(value_one - value_two) / ABS((value_one + value_two) / 2) > 0.1)
				{
					snapshot_for_converged_ = SMAX(snapshot_for_converged_, i - scale);
					break;
				} 
			}
		std::cout << "The scale is " << scale << "." << endl;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::calculateNewVariance(DoubleVec<Real> &current_result, 
		StdVec<Real> &local_meanvalue, StdVec<Real> &variance, StdVec<Real> &variance_new)
	{
		for (int j = 0; j != this->j_; ++j)
		{
			for (int i = snapshot_for_converged_; i != this->i_; ++i)
				variance_new[j] += std::pow((current_result[j][i] - local_meanvalue[j]), 2);
			variance_new[j] = SMAX((variance_new[j] / (this->i_ - snapshot_for_converged_)), variance[j], std::pow(local_meanvalue[j] * 1.0e-2, 2));
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::calculateNewVariance(DoubleVec<Vecd> &current_result, 
		StdVec<Vecd> &local_meanvalue, StdVec<Vecd> &variance, StdVec<Vecd> &variance_new)
	{
		for (int j = 0; j != this->j_; ++j)
			for (int index_i = 0; index_i != current_result[0][0].size(); ++index_i)
			{
				for (int i = snapshot_for_converged_; i != this->i_; ++i)
					variance_new[j][index_i] += std::pow((current_result[j][i][index_i] - local_meanvalue[j][index_i]), 2);
				variance_new[j][index_i] = SMAX((variance_new[j][index_i] / (this->i_ - snapshot_for_converged_)), variance[j][index_i], std::pow(local_meanvalue[j][index_i] * 1.0e-2, 2));
			}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::calculateNewVariance(DoubleVec<Matd> &current_result, 
		StdVec<Matd> &local_meanvalue, StdVec<Matd> &variance, StdVec<Matd> &variance_new)
	{
		for (int j = 0; j != this->j_; ++j)
			for (int index_i = 0; index_i != current_result[0][0].size(); ++index_i)
				for (int index_j = 0; index_j != current_result[0][0].size(); ++index_j)
				{
					for (int i = snapshot_for_converged_; i != this->i_; ++i)
						variance_new[j][index_i][index_j] += std::pow((current_result[j][i][index_i][index_j] - local_meanvalue[j][index_i][index_j]), 2);
					variance_new[j][index_i][index_j] = SMAX((variance_new[j][index_i][index_j] / (this->i_ - snapshot_for_converged_)), variance[j][index_i][index_j], std::pow(local_meanvalue[j][index_i][index_j] * 1.0e-2, 2));
				}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestTimeAveraged<ObserveMethodType>::compareParameter(string par_name,
		StdVec<Real> &parameter, StdVec<Real> &parameter_new, Real &threshold)
	{
		int count = 0;
		for (int j = 0; j != this->j_; ++j)
		{
			if ((par_name == "meanvalue") && (ABS(parameter[j]) < 0.005) && (ABS(parameter_new[j]) < 0.005))
			{
				std::cout << "The old meanvalue is " << parameter[j] << ", and the new meanvalue is " << parameter_new[j]
					<< ". So this variable will be ignored due to its tiny effect." << endl;
				continue;
			}
			Real relative_value_ = ABS((parameter[j] - parameter_new[j]) / (parameter_new[j] + TinyReal));
			if (relative_value_ > threshold)
			{
				std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "]"
					<< " is not converged, and difference is " << relative_value_ << endl;
				count++;
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestTimeAveraged<ObserveMethodType>::compareParameter(string par_name,
		StdVec<Vecd> &parameter, StdVec<Vecd> &parameter_new, Vecd &threshold)
	{
		int count = 0;
		for (int j = 0; j != this->j_; ++j)
			for (int index_i = 0; index_i != parameter[0].size(); ++index_i)
			{
				if ((par_name == "meanvalue") && (ABS(parameter[j][index_i]) < 0.001) && (ABS(parameter_new[j][index_i]) < 0.001))
				{
					std::cout << "The old meanvalue is " << parameter[j][index_i] << ", and the new meanvalue is " << parameter_new[j][index_i]
						<< ". So this variable will be ignored due to its tiny effect." << endl;
					continue;
				}
				Real relative_value_ = ABS((parameter[j][index_i] - parameter_new[j][index_i]) / (parameter_new[j][index_i] + TinyReal));
				if (relative_value_ > threshold[index_i])
				{
					std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "][" << index_i << "]"
						<< " is not converged, and difference is " << relative_value_ << endl;
					count++;
				}
			}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestTimeAveraged<ObserveMethodType>::compareParameter(string par_name,
		StdVec<Matd> &parameter, StdVec<Matd> &parameter_new, Matd &threshold)
	{
		int count = 0;
		for (int j = 0; j != this->j_; ++j)
			for (int index_i = 0; index_i != parameter[0].size(); ++index_i)
				for (int index_j = 0; index_j != parameter[0].size(); ++index_i)
				{
					if ((par_name == "meanvalue") && (ABS(parameter[j][index_i][index_j]) < 0.001) && (ABS(parameter_new[j][index_i][index_j]) < 0.001))
					{
						std::cout << "The old meanvalue is " << parameter[j][index_i][index_j] << ", and the new meanvalue is " << parameter_new[j][index_i][index_j] 
							<< ". So this variable will be ignored due to its tiny effect." << endl;
						continue;
					}
					Real relative_value_ = ABS((parameter[j][index_i][index_j] - parameter_new[j][index_i][index_j]) / (parameter_new[j][index_i][index_j] + TinyReal));
					if (relative_value_ > threshold[index_i][index_j])
					{
						std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "][" << index_i << "][" << index_j << "]"
							<< " is not converged, and difference is " << relative_value_ << endl;
						count++;
					}
				}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestTimeAveraged<ObserveMethodType>::testNewResult(DoubleVec<Real> &current_result, 
		StdVec<Real> &meanvalue, StdVec<Real> &local_meanvalue, StdVec<Real> &variance)
	{
		int count = 0;
		for (int j = 0; j != this->j_; ++j)
		{
			for (int i = snapshot_for_converged_; i != this->i_; ++i)
			{
				variance_new_[j] += std::pow((current_result[i][j] - local_meanvalue[j]), 2);
			}
			variance_new_[j] = variance_new_[j] / (this->i_ - snapshot_for_converged_);
			if ((ABS(meanvalue[j]) < 0.005) && (ABS(local_meanvalue[j]) < 0.005))
			{
				std::cout << "The old meanvalue is " << meanvalue[j] << ", and the current meanvalue is " << local_meanvalue[j]
					<< ". So this variable will not be tested due to its tiny effect." << endl;
				continue;
			}
			Real relative_value_ = ABS((meanvalue[j] - local_meanvalue[j]) / meanvalue[j]);
			if (relative_value_ > 0.1 || variance_new_[j] > (1.01 * variance[j]))
			{
				std::cout << this->quantity_name_ << "[" << j << "] is beyond the exception !" << endl;
				std::cout << "The meanvalue is " << meanvalue[j] << ", and the current meanvalue is " << local_meanvalue[j] << endl;
				std::cout << "The variance is " << variance[j] << ", and the current variance is " << variance_new_[j] << endl;
				count++;
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestTimeAveraged<ObserveMethodType>::testNewResult(DoubleVec<Vecd> &current_result, 
		StdVec<Vecd> &meanvalue, StdVec<Vecd> &local_meanvalue, StdVec<Vecd> &variance)
	{
		int count = 0;
		for (int j = 0; j != this->j_; ++j)
		{
			for (int index_i = 0; index_i != meanvalue_[0].size(); ++index_i)
			{
				for (int i = snapshot_for_converged_; i != this->i_; ++i)
				{
					variance_new_[j][index_i] += std::pow((current_result[i][j][index_i] - local_meanvalue[j][index_i]), 2);
				}
				variance_new_[j][index_i] = variance_new_[j][index_i] / (this->i_ - snapshot_for_converged_);
				if ((ABS(meanvalue[j][index_i]) < 0.005) && (ABS(local_meanvalue[j][index_i]) < 0.005))
				{
					std::cout << "The old meanvalue is " << meanvalue[j][index_i] << ", and the current meanvalue is " << local_meanvalue[j][index_i]
						<< ". So this variable will not be tested due to its tiny effect." << endl;
					continue;
				}
				Real relative_value_ = ABS((meanvalue[j][index_i] - local_meanvalue[j][index_i]) / meanvalue[j][index_i]);
				if (relative_value_ > 0.1 || (variance_new_[j][index_i] > 1.01 * variance[j][index_i]))
				{
					std::cout << this->quantity_name_ << "[" << j << "][" << index_i << "] is beyond the exception !" << endl;
					std::cout << "The meanvalue is " << meanvalue[j][index_i] << ", and the current meanvalue is " << local_meanvalue[j][index_i] << endl;
					std::cout << "The variance is " << variance[j][index_i] << ", and the new variance is " << variance_new_[j][index_i] << endl;
					count++;
				}
			}
		}
		return count;
	};
	//=================================================================================================// 
	template<class ObserveMethodType>
	int RegressionTestTimeAveraged<ObserveMethodType>::testNewResult(DoubleVec<Matd> &current_result, 
		StdVec<Matd> &meanvalue, StdVec<Matd> &local_meanvalue, StdVec<Matd> &variance)
	{
		int count = 0;
		for (int j = 0; j != this->j_; ++j)
		{
			for (int index_i = 0; index_i != meanvalue[0].size(); ++index_i)
			{
				for (int index_j = 0; index_j != meanvalue[0].size(); ++index_j)
				{
					for (int i = snapshot_for_converged_; i != this->i_; ++i)
					{
						variance_new_[j][index_i][index_j] += std::pow((current_result[i][j][index_i][index_j] - local_meanvalue[j][index_i][index_j]), 2);
					}
					variance_new_[j][index_i][index_j] = variance_new_[j][index_i][index_j] / (this->i_ - snapshot_for_converged_);
					if ((ABS(meanvalue[j][index_i][index_j]) < 0.005) && (ABS(local_meanvalue[j][index_i][index_j]) < 0.005))
					{
						std::cout << "The old meanvalue is " << meanvalue[j][index_i][index_j] << ", and the new meanvalue is " << local_meanvalue[j][index_i][index_j]
							<< ". So this variable will not be tested due to its tiny effect. " << endl;
						continue;
					}
					Real relative_value_ = ABS((meanvalue_[j][index_i][index_j] - local_meanvalue[j][index_i][index_j]) / meanvalue[j][index_i][index_j]);
					if (relative_value_ > 0.1 || variance_new_[j][index_i][index_j] > 1.01 * variance[j][index_i][index_j])
					{
						std::cout << this->quantity_name_ << "[" << j << "][" << index_i << "][" << index_j << "] is beyond the exception !" << endl;
						std::cout << "The meanvalue is " << meanvalue[j][index_i][index_j] << ", and the new meanvalue is " << local_meanvalue[j][index_i][index_j] << endl;
						std::cout << "The variance is " << variance[j][index_i][index_j] << ", and the new variance is " << variance_new_[j][index_i][index_j] << endl;
						count++;
					}
				}
			}
		}
		return count;
	};
	//=================================================================================================//	
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::initializeThreshold(VariableType &threshold_mean, VariableType &threshold_variance)
	{
		threshold_mean_ = threshold_mean;
		threshold_variance_ = threshold_variance;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::settingupTheTest()
	{
		this->i_ = this->current_result_.size();
		this->j_ = this->current_result_[0].size();
		StdVec<VariableType> temp(this->j_);
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
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::readMeanVarianceFromXml()
	{
		if (this->number_of_run_ > 1)
		{
			mean_variance_xml_engine_in_.loadXmlFile(mean_variance_filefullpath_);
			SimTK::Xml::Element meanvalue_element_ = mean_variance_xml_engine_in_.getChildElement("MeanValue_Element");
			SimTK::Xml::element_iterator ele_ite_mean_ = meanvalue_element_.element_begin();
			for (; ele_ite_mean_ != meanvalue_element_.element_end(); ++ele_ite_mean_)
				for (int particle_n_ = 0; particle_n_ != this->j_; ++particle_n_)
				{
					std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
					mean_variance_xml_engine_in_.getRequiredAttributeValue(ele_ite_mean_, attribute_name_, meanvalue_[particle_n_]);
				}

			SimTK::Xml::Element variance_element_ = mean_variance_xml_engine_in_.getChildElement("Variance_Element");
			SimTK::Xml::element_iterator ele_ite_variance_ = variance_element_.element_begin();
			for (; ele_ite_variance_ != variance_element_.element_end(); ++ele_ite_variance_)
				for (int particle_n_ = 0; particle_n_ != this->j_; ++particle_n_)
				{
					std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
					mean_variance_xml_engine_in_.getRequiredAttributeValue(ele_ite_variance_, attribute_name_, variance_[particle_n_]);
				}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::searchForStartPoint()
	{
		snapshot_for_converged_ = 0;
		searchSteadyStart(this->current_result_);
		std::cout << "The snapshot for converged is " << snapshot_for_converged_ << endl;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::filterExtremeValues()
	{
		filterLocalResult(this->current_result_);
		filefullpath_filter_output_ = this->input_folder_path_ + "/" + this->body_name_
			+ "_" + this->quantity_name_ + ".dat";
		std::ofstream out_file(filefullpath_filter_output_.c_str(), std::ios::app);
		out_file << "run_time" << "   ";
		for (int i = 0; i != this->j_; ++i)
		{
			std::string quantity_name_i = this->quantity_name_ + "[" + std::to_string(i) + "]";
			this->plt_engine_.writeAQuantityHeader(out_file, this->current_result_[0][0], quantity_name_i);
		}
		out_file << "\n";
		out_file.close();

		for (int i = 0; i != this->i_; ++i)
		{
			std::ofstream out_file(filefullpath_filter_output_.c_str(), std::ios::app);
			out_file << this->element_tag_[i] << "   ";
			for (int j = 0; j != this->j_; ++j)
			{
				this->plt_engine_.writeAQuantity(out_file, this->current_result_[i][j]);
			}
			out_file << "\n";
			out_file.close();
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::updateMeanVariance()
	{
		for (int j = 0; j != this->j_; ++j)
		{
			for (int i = snapshot_for_converged_; i != this->i_; ++i)
			{
				local_meanvalue_[j] += this->current_result_[i][j];
			}
			local_meanvalue_[j] = local_meanvalue_[j] / (this->i_ - snapshot_for_converged_);
			meanvalue_new_[j] = (local_meanvalue_[j] + meanvalue_[j] * (this->number_of_run_ - 1)) / this->number_of_run_;
		}
		calculateNewVariance(this->current_result_ji_, local_meanvalue_, variance_, variance_new_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::writeMeanVarianceToXml()
	{
		mean_variance_xml_engine_out_.addElementToXmlDoc("MeanValue_Element");
		SimTK::Xml::Element meanvalue_element_ = mean_variance_xml_engine_out_.getChildElement("MeanValue_Element");
		mean_variance_xml_engine_out_.addChildToElement(meanvalue_element_, "Snapshot_MeanValue");
		for (int particle_n_ = 0; particle_n_ != this->j_; ++particle_n_)
		{
			SimTK::Xml::element_iterator ele_ite_mean = meanvalue_element_.element_begin();
			std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
			mean_variance_xml_engine_out_.setAttributeToElement(ele_ite_mean, attribute_name_, meanvalue_new_[particle_n_]);
		}
		mean_variance_xml_engine_out_.addElementToXmlDoc("Variance_Element");
		SimTK::Xml::Element variance_element_ = mean_variance_xml_engine_out_.getChildElement("Variance_Element");
		mean_variance_xml_engine_out_.addChildToElement(variance_element_, "Snapshot_Variance");
		for (int particle_n_ = 0; particle_n_ != this->j_; ++particle_n_)
		{
			SimTK::Xml::element_iterator ele_ite_variance = variance_element_.element_begin();
			std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
			mean_variance_xml_engine_out_.setAttributeToElement(ele_ite_variance, attribute_name_, variance_new_[particle_n_]);
		}
		mean_variance_xml_engine_out_.writeToXmlFile(mean_variance_filefullpath_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	bool RegressionTestTimeAveraged<ObserveMethodType>::compareMeanVariance()
	{
		int count_not_converged_m = 0;
		int count_not_converged_v = 0;
		count_not_converged_m = this->compareParameter("meanvalue", this->meanvalue_, this->meanvalue_new_, this->threshold_mean_);
		count_not_converged_v = this->compareParameter("variance", this->variance_, this->variance_new_, this->threshold_variance_);
		if (count_not_converged_m == 0)
		{
			std::cout << "The meanvalue of " << this->quantity_name_ << " are converged now." << endl;
			if (count_not_converged_v == 0)
			{
				if (this->label_for_repeat_ == 4)
				{
					this->converged = "true";
					this->label_for_repeat_++;
					std::cout << "The meanvalue and variance of " << this->quantity_name_ << " are converged enough times, and run will stop now." << endl;
					return true;
				}
				else
				{
					this->converged = "false";
					this->label_for_repeat_++;
					std::cout << "The variance of " << this->quantity_name_ << " are also converged, and this is the " << this->label_for_repeat_
						<< " times. They should be converged more times to be stable." << endl;
					return false;
				}
			}
			else if (count_not_converged_v != 0)
			{
				this->converged = "false";
				this->label_for_repeat_ = 0;
				std::cout << "The variance of " << this->quantity_name_ << " are not converged " << count_not_converged_v << " times." << endl;
				return false;
			};
		}
		else if (count_not_converged_m != 0)
		{
			this->converged = "false";
			this->label_for_repeat_ = 0;
			std::cout << "The meanvalue of " << this->quantity_name_ << " are not converged " << count_not_converged_m << " times." << endl;
			return false;
		}
	};
	//=================================================================================================//	
	template<class ObserveMethodType>
	void RegressionTestTimeAveraged<ObserveMethodType>::resultTest()
	{
		int test_wrong = 0;
		
		for (int j = 0; j != this->j_; ++j)
		{
			for (int i = snapshot_for_converged_; i != this->i_; ++i)
				local_meanvalue_[j] += this->current_result_[i][j];
			local_meanvalue_[j] = local_meanvalue_[j] / (this->i_-snapshot_for_converged_);
		}

		test_wrong = testNewResult(this->current_result_, meanvalue_, local_meanvalue_, variance_);
		if (test_wrong == 0)
			std::cout << "The result of " << this->quantity_name_ << "is correct based on the time-averaged regression test!" << endl;
		else
		{
			std::cout << "There are " << test_wrong << " particles are not within the expected range." << endl;
			std::cout << "Please try again. If it still post this conclusion, the result is not correct!" << endl;
			exit(1);
		}
	};
	//=================================================================================================//
};