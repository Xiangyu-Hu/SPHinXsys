/**
 * @file 	regression_testing.cpp
 * @author	Bo Zhang, Xiangyu Hu
 */

#include "meanvalue_variance_method.h"

 //=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::InitializeThreshold(
		Real& threshold_mean, Real& threshold_variance)
	{
		threshold_mean_ = threshold_mean;
		threshold_variance_ = threshold_variance;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::InitializeThreshold(
		Vecd& threshold_mean, Vecd& threshold_variance)
	{
		for (size_t index_i = 0; index_i != threshold_mean_.size(); ++index_i)
		{
			threshold_mean_[index_i] = threshold_mean[index_i];
			threshold_variance_[index_i] = threshold_variance[index_i];
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::InitializeThreshold(
		Matd& threshold_mean, Matd& threshold_variance)
	{
		for (size_t index_i = 0; index_i != threshold_mean_.size(); ++index_i)
		{
			for (size_t index_j = 0; index_j != threshold_mean_.size(); ++index_j)
			{
				threshold_mean_[index_i][index_j] = threshold_mean[index_i][index_j];
				threshold_variance_[index_i][index_j] = threshold_variance[index_i][index_j];
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::GetNewMeanValue(DataVec<Real>& current_result, 
		DataVec<Real>& meanvalue, DataVec<Real>& meanvalue_new)
	{
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				meanvalue_new[i][j] = (meanvalue[i][j] * (number_of_run_ - 1) 
					+ current_result[i][j]) / number_of_run_;
			}
		}	
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::GetNewMeanValue(DataVec<Vecd>& current_result, 
		DataVec<Vecd>& meanvalue, DataVec<Vecd>& meanvalue_new)
	{
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t index_i = 0; index_i != meanvalue[0][0].size(); ++index_i)
				{
					meanvalue_new[i][j][index_i] = (meanvalue[i][j][index_i] * (number_of_run_ - 1)
						+ current_result[i][j][index_i]) / number_of_run_;
				}
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::GetNewMeanValue(DataVec<Matd>& current_result, 
		DataVec<Matd>& meanvalue, DataVec<Matd>& meanvalue_new)
	{	
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t index_i = 0; index_i != meanvalue[0][0].size(); ++index_i)
				{
					for (size_t index_j = 0; index_j != meanvalue[0][0].size(); ++index_j)
					{
						meanvalue_new[i][j][index_i][index_j] = (meanvalue[i][j][index_i][index_j] * 
							(number_of_run_ - 1) + current_result[i][j][index_i][index_j]) / number_of_run_;
					}
				}
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::GetNewVariance(StdVec<DataVec<Real>>& result, 
		DataVec<Real>& meanvalue_new, DataVec<Real>& variance, DataVec<Real>& variance_new)
	{
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t n = 0; n != number_of_run_; ++n)
				{
					variance_new[i][j] = SMAX(variance[i][j], variance_new[i][j],
						std::pow((result[n][i][j] - meanvalue_new[i][j]), 2));
					variance_new[i][j] = variance_new[i][j] < std::pow(meanvalue_new[i][j] * 1.0e-6, 2) ?
						std::pow(meanvalue_new[i][j] * 1.0e-6, 2) : variance_new[i][j];
				}
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::GetNewVariance(StdVec<DataVec<Vecd>>& result,
		DataVec<Vecd>& meanvalue_new, DataVec<Vecd>& variance, DataVec<Vecd>& variance_new)
	{
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t n = 0; n != number_of_run_; ++n)
				{
					for (size_t index_i = 0; index_i != variance[0][0].size(); ++index_i)
					{
						variance_new[i][j][index_i] = SMAX(variance[i][j][index_i], variance_new[i][j][index_i],
							std::pow((result[n][i][j][index_i] - meanvalue_new[i][j][index_i]), 2));
						variance_new[i][j][index_i] = variance_new[i][j][index_i] < 
							std::pow(meanvalue_new[i][j][index_i] * 1.0e-6, 2) ?
							std::pow(meanvalue_new[i][j][index_i] * 1.0e-6, 2) : variance_new[i][j][index_i];
					}
				}
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::GetNewVariance(StdVec<DataVec<Matd>>& result,
		DataVec<Matd>& meanvalue_new, DataVec<Matd>& variance, DataVec<Matd>& variance_new)
	{
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t n = 0; n != number_of_run_; ++n)
				{
					for (size_t index_i = 0; index_i != variance[0][0].size(); ++index_i)
					{
						for (size_t index_j = 0; index_j != variance[0][0].size(); ++index_j)
						{
							variance_new[i][j][index_i][index_j] = SMAX(variance[i][j][index_i][index_j],
								variance_new[i][j][index_i][index_j], std::pow((result[n][i][j][index_i][index_j]
									- meanvalue_new[i][j][index_i][index_j]), 2));
							variance_new[i][j][index_i][index_j] = variance_new[i][j][index_i][index_j] < 
								std::pow(meanvalue_new[i][j][index_i][index_j] * 1.0e-6, 2) ? 
								std::pow(meanvalue_new[i][j][index_i][index_j] * 1.0e-6, 2) : 
								variance_new[i][j][index_i][index_j];
						}
					}
				}
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::WriteDataToXmlMemory(XmlEngine xmlengine, 
		SimTK::Xml::Element element, const DataVec<Real>& quantity)
	{
		for (size_t snapshot_n_ = 0; snapshot_n_ != SMIN(i_, number_of_snapshot_old_); ++snapshot_n_)
		{
			std::string element_name_ = this->element_tag_[snapshot_n_];
			xmlengine.addChildToElement(element, element_name_);
			for (size_t particle_n_ = 0; particle_n_ != j_; ++particle_n_)
			{
				SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name_);
				std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
				xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity[snapshot_n_][particle_n_]);
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::WriteDataToXmlMemory(XmlEngine xmlengine,
		SimTK::Xml::Element element, const DataVec<Vecd>& quantity)
	{
		for (size_t snapshot_n_ = 0; snapshot_n_ != SMIN(i_, number_of_snapshot_old_); ++snapshot_n_)
		{
			std::string element_name_ = this->element_tag_[snapshot_n_];
			xmlengine.addChildToElement(element, element_name_);
			for (size_t particle_n_ = 0; particle_n_ != j_; ++particle_n_)
			{
				SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name_);
				std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
				xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity[snapshot_n_][particle_n_]);
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::WriteDataToXmlMemory(XmlEngine xmlengine,
		SimTK::Xml::Element element, const DataVec<Matd>& quantity)
	{
		for (size_t snapshot_n_ = 0; snapshot_n_ != SMIN(i_, number_of_snapshot_old_); ++snapshot_n_)
		{
			std::string element_name_ = this->element_tag_[snapshot_n_];
			xmlengine.addChildToElement(element, element_name_);
			for (size_t particle_n_ = 0; particle_n_ != j_; ++particle_n_)
			{
				SimTK::Xml::element_iterator ele_ite = element.element_begin(element_name_);
				std::string attribute_name_ = this->quantity_name_ + "_" + std::to_string(particle_n_);
				xmlengine.setAttributeToElement(ele_ite, attribute_name_, quantity[snapshot_n_][particle_n_]);
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	size_t RegressionTesting<ObserveMethodType>::CompareParameter(string par_name, 
		DataVec<Real>& parameter, DataVec<Real>& parameter_new, Real& threshold)
	{
		size_t count = 0;
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				if ((ABS(parameter[i][j] - parameter_new[i][j]) / (parameter_new[i][j] + TinyReal)) > threshold)
				{
					std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "] in "
						<< this->element_tag_[i] << "is not converged, and difference is "
						<< (ABS(parameter[i][j] - parameter_new[i][j]) / (parameter_new[i][j] + TinyReal)) << endl;
					count++;
				}
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	size_t RegressionTesting<ObserveMethodType>::CompareParameter(string par_name, 
		DataVec<Vecd>& parameter, DataVec<Vecd>& parameter_new, Vecd& threshold)
	{
		size_t count = 0;
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t index_i = 0; index_i != parameter[0][0].size(); ++index_i)
				{
					if ((ABS(parameter[i][j][index_i] - parameter_new[i][j][index_i]) 
						/ (parameter_new[i][j][index_i] + TinyReal)) > threshold[index_i])
					{
						std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "][" <<
							index_i << "] in " << this->element_tag_[i] << "is not converged, and difference is "
							<< (ABS(parameter[i][j][index_i] - parameter_new[i][j][index_i]) /
							(parameter_new[i][j][index_i] + TinyReal)) << endl;
						count++;
					}
				}
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	size_t RegressionTesting<ObserveMethodType>::CompareParameter(string par_name, 
		DataVec<Matd>& parameter, DataVec<Matd>& parameter_new, Matd& threshold)
	{
		size_t count = 0;
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t index_i = 0; index_i != parameter[0][0].size(); ++index_i)
				{
					for (size_t index_j = 0; index_j != parameter[0][0].size(); ++index_j)
					{
						if ((ABS(parameter[i][j][index_i][index_j] - parameter_new[i][j][index_i][index_j]) /
							(parameter_new[i][j][index_i][index_j] + TinyReal)) > threshold[index_i][index_j])
						{
							std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "][" <<
								index_i << "][" << index_j << " ] in " << this->element_tag_[i] <<
								"is not converged, and difference is " << (ABS(parameter[i][j][index_i][index_j]
									- parameter_new[i][j][index_i][index_j]) /
								(parameter_new[i][j][index_i][index_j] + TinyReal)) << endl;
							count++;
						}
					}
				}
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	size_t RegressionTesting<ObserveMethodType>::TestingNewResult(size_t diff,
		DataVec<Real>& current_result, DataVec<Real>& meanvalue, DataVec<Real>& variance)
	{
		size_t count = 0;
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				if (std::pow(current_result[i][j] - meanvalue[i + diff][j], 2) > variance[i + diff][j])
				{
					std::cout << this->quantity_name_ << "[" << j << "] in " << this->element_tag_[i] <<
						" is beyond the expection, and difference is " <<
						(ABS((std::pow(current_result[i][j] - meanvalue[i + diff][j], 2) - variance[i + diff][j]))
						/ variance[i + diff][j]) << endl;
					count++;
				}
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	size_t RegressionTesting<ObserveMethodType>::TestingNewResult(size_t diff, 
		DataVec<Vecd>& current_result, DataVec<Vecd>& meanvalue, DataVec<Vecd>& variance)
	{
		size_t count = 0;
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t index_i = 0; index_i != meanvalue[0][0].size(); ++index_i)
				{
					if (std::pow(current_result[i][j][index_i] - meanvalue[i + diff][j][index_i], 2) 
						> variance[i + diff][j][index_i])
					{
						std::cout << this->quantity_name_ << "[" << j << "][" << index_i << "] in "
							<< this->element_tag_[i] << " is beyond the expection, and difference is "
							<< (ABS((std::pow(current_result[i][j][index_i] - meanvalue[i + diff][j][index_i], 2)
								- variance[i + diff][j][index_i])) / variance[i + diff][j][index_i]) << endl;
						count++;
					}
				}
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	size_t RegressionTesting<ObserveMethodType>::TestingNewResult(size_t diff, 
		DataVec<Matd>& current_result, DataVec<Matd>& meanvalue, DataVec<Matd>& variance)
	{
		size_t count = 0;
		for (size_t i = 0; i != SMIN(i_, number_of_snapshot_old_); ++i)
		{
			for (size_t j = 0; j != j_; ++j)
			{
				for (size_t index_i = 0; index_i != meanvalue[0][0].size(); ++index_i)
				{
					for (size_t index_j = 0; index_j != meanvalue[0][0].size(); ++index_j)
					{
						if (std::pow(current_result[i][j][index_i][index_j] - 
							meanvalue[i + diff][j][index_i][index_j], 2) >
							variance[i + diff][j][index_i][index_j])
						{
							std::cout << this->quantity_name_ << "[" << j << "][" << index_i << "] in "
								<< this->element_tag_[i] << " is beyond the expection, and difference is "
								<< (ABS((std::pow(current_result[i][j][index_i][index_j] - 
									meanvalue[i + diff][j][index_i][index_j], 2)
									- variance[i + diff][j][index_i][index_j])) / 
								variance[i + diff][j][index_i][index_j]) << endl;
							count++;
						}
					}
					
				}
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::SettingupAndCorrection()
	{
		/* obtain the size of result. */
		i_ = this->current_result_.size();
		j_ = this->current_result_[0].size();

		if (number_of_run_ > 1)
		{
			if (!fs::exists(result_filefullpath_))
			{
				std::cout << "\n Error: the input file:" << result_filefullpath_ << " is not exists" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
			else if (!fs::exists(meanvalue_filefullpath_))
			{
				std::cout << "\n Error: the input file:" << meanvalue_filefullpath_ << " is not exists" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
			else if (!fs::exists(variance_filefullpath_))
			{
				std::cout << "\n Error: the input file:" << variance_filefullpath_ << " is not exists" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
			else
			{
				/** load the result from .xml file. */
				result_xml_engine_in_.loadXmlFile(result_filefullpath_);
				/** load the meanvalue from .xml file. */
				meanvalue_xml_engine_in_.loadXmlFile(meanvalue_filefullpath_);
				/** load the variance from .xml file. */
				variance_xml_engine_in_.loadXmlFile(variance_filefullpath_);

				number_of_snapshot_old_ = std::distance(meanvalue_xml_engine_in_.root_element_.element_begin(),
					meanvalue_xml_engine_in_.root_element_.element_end());

				DataVec<VariableType> meanvalue_temp_(SMAX(i_, number_of_snapshot_old_), StdVec<VariableType>(j_));
				DataVec<VariableType> variance_temp_(SMAX(i_, number_of_snapshot_old_), StdVec<VariableType>(j_));
				meanvalue_ = meanvalue_temp_;
				variance_ = variance_temp_;

				if (number_of_snapshot_old_ < i_)
				{
					difference_ = i_ - number_of_snapshot_old_;
					for (size_t delete_ = 0; delete_ != difference_; ++delete_)
					{
						this->current_result_.pop_back();
					}
				}
				else if (number_of_snapshot_old_ > i_)
				{
					difference_ = number_of_snapshot_old_ - i_;
				}
				else
				{
					difference_ = 0;
				}
			}
		}
		else if (number_of_run_ == 1)
		{
			number_of_snapshot_old_ = i_;
			DataVec<VariableType> meanvalue_temp_(i_, StdVec<VariableType>(j_));
			DataVec<VariableType> variance_temp_(i_, StdVec<VariableType>(j_));
			result_.push_back(this->current_result_);
			meanvalue_ = meanvalue_temp_;
			variance_ = variance_temp_;
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::ReadResultFromXml()
	{
		if (number_of_run_ > 1)
		{
			/** read result from .xml */
			DataVec<VariableType> result_in_(SMAX(i_, number_of_snapshot_old_), StdVec<VariableType>(j_));
			for (size_t run_n_ = 0; run_n_ != number_of_run_ - 1; ++run_n_)
			{
				std::string node_name_ = "Round_" + std::to_string(run_n_);
				SimTK::Xml::Element father_element_ = result_xml_engine_in_.getChildElement(node_name_);
				for (size_t particle_n_ = 0; particle_n_ != j_; ++particle_n_)
				{
					this->ReadDataFromXmlMemory(result_xml_engine_in_, father_element_, particle_n_, result_in_);
				}
				DataVec<VariableType> result_temp_ = result_in_;
				for (size_t delete_ = 0; delete_ != difference_; ++delete_)
				{
					result_temp_.pop_back();
				}
				result_.push_back(result_temp_);
			}
			result_.push_back(this->current_result_);
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::ReadMeanAndVarianceToXml()
	{
		if (number_of_run_ > 1)
		{
			/** read mean value from .xml file. */
			SimTK::Xml::Element element_name_meanvalue_ = meanvalue_xml_engine_in_.root_element_;
			SimTK::Xml::Element element_name_variance_ = variance_xml_engine_in_.root_element_;
			for (size_t particle_n_ = 0; particle_n_ != j_; ++particle_n_)
			{
				this->ReadDataFromXmlMemory(meanvalue_xml_engine_in_, element_name_meanvalue_, particle_n_, meanvalue_);
				this->ReadDataFromXmlMemory(variance_xml_engine_in_, element_name_variance_, particle_n_, variance_);
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::UpdateMeanValueAndVariance()
	{
		if (number_of_run_ > 1)
		{
			for (size_t delete_ = 0; delete_ != difference_; ++delete_)
			{
				meanvalue_.pop_back();
				variance_.pop_back();
			}
		}
		meanvalue_new_ = meanvalue_;
		variance_new_ = variance_;
		/** update meanvalue of result. */
		GetNewMeanValue(this->current_result_, meanvalue_, meanvalue_new_);
		/** update variance of result. */
		GetNewVariance(result_, meanvalue_new_, variance_, variance_new_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	bool RegressionTesting<ObserveMethodType>::CompareMeanValueAndVariance()
	{
		size_t count_not_converged_m = 0;
		size_t count_not_converged_v = 0;
		/** determine if the average value has converged. */
		count_not_converged_m = CompareParameter("meanvalue", meanvalue_, meanvalue_new_, threshold_mean_);
		count_not_converged_v = CompareParameter("variance", variance_, variance_new_, threshold_variance_);
		if (count_not_converged_m == 0)
		{
			std::cout << "The meanvalue of " << this->quantity_name_ << " are converged now." << endl;
			if (count_not_converged_v == 0)
			{
				if (label_for_repeat_ == 1)
				{
					converged = "true";
					std::cout << "The meanvalue and variance of " << this->quantity_name_ <<
						" are converged enough times, and rum will stop now." << endl;
					return true;
				}
				else
				{
					converged = "false";
					label_for_repeat_++;
					std::cout << "The variance of " << this->quantity_name_ << " are also converged,"
						"but they should be converged again to be stable." << endl;
					return false;
				}
			}
			else if (count_not_converged_v != 0)
			{
				converged = "false";
				label_for_repeat_ = 0;
				std::cout << "The variance of " << this->quantity_name_ << " is not converged "<<
					count_not_converged_v << " times." << endl;
				return false;
			}
		}
		else if (count_not_converged_m != 0)
		{
			converged = "false";
			label_for_repeat_ = 0;
			std::cout << "The meanvalue of " << this->quantity_name_ << " is not converged " << 
				count_not_converged_m << " times." << endl;
			return false;

		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::WriteResultToXml()
	{
		/** Write result to .xml file */
		for (size_t run_n_ = 0; run_n_ != number_of_run_; ++run_n_)
		{
			std::string node_name_ = "Round_" + std::to_string(run_n_);
			result_xml_engine_out_.addElementToXmlDoc(node_name_);
			SimTK::Xml::Element father_element_ =
				result_xml_engine_out_.getChildElement(node_name_);
			WriteDataToXmlMemory(result_xml_engine_out_, father_element_, result_[run_n_]);
		}
		/** write the result in XmlEngine to xml file. */
		result_xml_engine_out_.writeToXmlFile(result_filefullpath_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::WriteMeanAndVarianceToXml()
	{
		/** Write meanvalue and variance to .xml file */
		SimTK::Xml::Element element_meanvalue_ = meanvalue_xml_engine_out_.root_element_;
		SimTK::Xml::Element element_variance_ = variance_xml_engine_out_.root_element_;
		WriteDataToXmlMemory(meanvalue_xml_engine_out_, element_meanvalue_, meanvalue_new_);
		WriteDataToXmlMemory(variance_xml_engine_out_, element_variance_, variance_new_);
		/** write the meanvalue and variance in XmlEngine to xml file. */
		meanvalue_xml_engine_out_.writeToXmlFile(meanvalue_filefullpath_);
		variance_xml_engine_out_.writeToXmlFile(variance_filefullpath_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	RegressionTesting<ObserveMethodType>::~RegressionTesting()
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
	}
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::ResultTesting()
	{
		/* compare the current result to the converged mean value and variance. */
		size_t test_wrong = 0;
		if (i_ < number_of_snapshot_old_)
		{
			test_wrong = TestingNewResult(difference_, this->current_result_, meanvalue_, variance_);
			if (test_wrong == 0)
				std::cout << "The result of " << this->quantity_name_ << 
				" is correct based on the regression testing!" << endl;
			else
			{
				std::cout << "There are " << test_wrong << 
					" snapshots are not within the expected range." << endl;
				std::cout << "Please try again. If it still post this sectence, "
					"the result is not correct!" << endl;
			}
		}
		else if (i_ >= number_of_snapshot_old_)
		{
			for (size_t delete_ = 0; delete_ != difference_; ++delete_)
			{
				meanvalue_.pop_back();
				variance_.pop_back();
			}
			test_wrong = TestingNewResult(0, this->current_result_, meanvalue_, variance_);
			if (test_wrong == 0)
				std::cout << "The result of " << this->quantity_name_ <<
				" is correct based on the regression testing!" << endl;
			else
			{
				std::cout << "There are " << test_wrong << 
					" snapshots are not within the expected range." << endl;
				std::cout << "Please try again. If it still post this sectence, "
					"the result is not correct!" << endl;
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTesting<ObserveMethodType>::GetMeanvalueInTime(size_t iteration)
	{
		VariableType average = { 0 };
		for (size_t i = iteration; i != i_; ++i)
		{
			average += this->current_result_[i][0];
		}
		average = average / (i_ - iteration);

		meanvalue_xml_engine_out_.addElementToXmlDoc("meanvalue");
		SimTK::Xml::element_iterator element_iterator = meanvalue_xml_engine_out_.root_element_.element_begin();
		meanvalue_xml_engine_out_.setAttributeToElement(element_iterator, "AverageViscousForce", average);
		meanvalue_xml_engine_out_.writeToXmlFile(meanvalue_filefullpath_);
	};
	//=================================================================================================//
}