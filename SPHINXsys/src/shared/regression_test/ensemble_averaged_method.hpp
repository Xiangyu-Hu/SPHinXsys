/**
 * @file 	ensemble_averaged_method.cpp
 * @author	Bo Zhang and Xiangyu Hu
 */

#include "ensemble_averaged_method.h"

 //=================================================================================================//
namespace SPH
{
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestEnsembleAveraged<ObserveMethodType>::calculateNewVariance(TripleVec<Real> &result,
		DoubleVec<Real> &meanvalue_new, DoubleVec<Real> &variance, DoubleVec<Real> &variance_new)
	{
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i)
			for (int j = 0; j != this->j_; ++j)
				for (int n = 0; n != this->number_of_run_; ++n)
					variance_new[i][j] = SMAX(variance[i][j], variance_new[i][j],
						std::pow((result[n][i][j] - meanvalue_new[i][j]), 2), std::pow(meanvalue_new[i][j] * 1.0e-2, 2));
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestEnsembleAveraged<ObserveMethodType>::calculateNewVariance(TripleVec<Vecd> &result, 
		DoubleVec<Vecd> &meanvalue_new, DoubleVec<Vecd> &variance, DoubleVec<Vecd> &variance_new)
	{
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i) 
			for (int j = 0; j != this->j_; ++j) 
				for (int n = 0; n != this->number_of_run_; ++n) 
					for (int index_i = 0; index_i != variance[0][0].size(); ++index_i) 
						variance_new[i][j][index_i] = SMAX(variance[i][j][index_i], variance_new[i][j][index_i],
							std::pow((result[n][i][j][index_i] - meanvalue_new[i][j][index_i]), 2),
							std::pow(meanvalue_new[i][j][index_i] * 1.0e-2, 2));
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestEnsembleAveraged<ObserveMethodType>::calculateNewVariance(TripleVec<Matd> &result,
		DoubleVec<Matd> &meanvalue_new, DoubleVec<Matd> &variance, DoubleVec<Matd> &variance_new)
	{
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i)
			for (int j = 0; j != this->j_; ++j)
				for (int n = 0; n != this->number_of_run_; ++n)
					for (size_t index_i = 0; index_i != variance[0][0].size(); ++index_i)
						for (size_t index_j = 0; index_j != variance[0][0].size(); ++index_j)
							variance_new[i][j][index_i][index_j] = SMAX(variance[i][j][index_i][index_j], variance_new[i][j][index_i][index_j],
								std::pow((result[n][i][j][index_i][index_j] - meanvalue_new[i][j][index_i][index_j]), 2),
								std::pow(meanvalue_new[i][j][index_i][index_j] * 1.0e-2, 2));
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestEnsembleAveraged<ObserveMethodType>::compareParameter(string par_name, 
		DoubleVec<Real> &parameter, DoubleVec<Real> &parameter_new, Real &threshold)
	{
		int count = 0;
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i) 
			for (int j = 0; j != this->j_; ++j)
			{
				Real relative_value_ = ABS((parameter[i][j] - parameter_new[i][j]) / (parameter_new[i][j] + TinyReal));
				if (relative_value_ > threshold)
				{
					std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "] in " << this->element_tag_[i]
						<< " is not converged, and difference is " << relative_value_ << endl;
					count++;
				}
			}	
		return count;
	};
	////=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestEnsembleAveraged<ObserveMethodType>::compareParameter(string par_name,
		DoubleVec<Vecd> &parameter, DoubleVec<Vecd> &parameter_new, Vecd &threshold)
	{
		int count = 0;
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i) 
			for (int j = 0; j != this->j_; ++j) 
				for (int index_i = 0; index_i != parameter[0][0].size(); ++index_i)
				{
					Real relative_value_ = ABS((parameter[i][j][index_i] - parameter_new[i][j][index_i]) / (parameter_new[i][j][index_i] + TinyReal));
					if (relative_value_ > threshold[index_i])
					{
						std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "][" << index_i << "] in " << this->element_tag_[i]
							<< " is not converged, and difference is " << relative_value_ << endl;
						count++;
					}
				}
		return count;
	};
	////=================================================================================================// 
	template<class ObserveMethodType>
	int RegressionTestEnsembleAveraged<ObserveMethodType>::compareParameter(string par_name,
		DoubleVec<Matd> &parameter, DoubleVec<Matd> &parameter_new, Matd &threshold)
	{
		int count = 0;
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i)
			for (int j = 0; j != this->j_; ++j)
				for (int index_i = 0; index_i != parameter[0][0].size(); ++index_i)
					for (int index_j = 0; index_j != parameter[0][0].size(); ++index_j)
					{
						Real relative_value_ = ABS(parameter[i][j][index_i][index_j] - parameter_new[i][j][index_i][index_j])
							/ (parameter_new[i][j][index_i][index_j] + TinyReal);
						if (relative_value_ > threshold[index_i][index_j])
						{
							std::cout << par_name << ": " << this->quantity_name_ << "[" << j << "][" << index_i << "][" << index_j << " ] in "
								<< this->element_tag_[i] << " is not converged, and difference is " << relative_value_ << endl;
							count++;
						}
					}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestEnsembleAveraged<ObserveMethodType>::testNewResult(int diff, DoubleVec<Real> &current_result,
		DoubleVec<Real> &meanvalue, DoubleVec<Real> &variance)
	{
		int count = 0;
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i)
		{
			for (int j = 0; j != this->j_; ++j)
			{
				Real relative_value_ = (std::pow(current_result[i][j] - meanvalue[i + diff][j], 2) - variance[i + diff][j]) / variance[i + diff][j];
				if (relative_value_ > 0.01)
				{
					std::cout << this->quantity_name_ << "[" << j << "] in " << this->element_tag_[i] << " is beyond the exception, and difference is "
						<< relative_value_ << endl;
					count++;
				}
			}
		}
		return count;
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	int RegressionTestEnsembleAveraged<ObserveMethodType>::testNewResult(int diff, DoubleVec<Vecd> &current_result,
		DoubleVec<Vecd> &meanvalue, DoubleVec<Vecd> &variance)
	{
		int count = 0;
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i)
		{
			for (int j = 0; j != this->j_; ++j)
			{
				for (int index_i = 0; index_i != meanvalue[0][0].size(); ++index_i)
				{
					Real relative_value_ = (std::pow(current_result[i][j][index_i] - meanvalue[i + diff][j][index_i], 2) - variance[i + diff][j][index_i]) / variance[i + diff][j][index_i];
					if (relative_value_ > 0.01)
					{
						std::cout << this->quantity_name_ << "[" << j << "][" << index_i << "] in " << this->element_tag_[i] << " is beyond the exception, and difference is "
							<< relative_value_ << endl;
						count++;
					}
				}
			}	
		}
		return count;
	};
	//=================================================================================================// 
	template<class ObserveMethodType>
	int RegressionTestEnsembleAveraged<ObserveMethodType>::testNewResult(int diff, DoubleVec<Matd> &current_result,
		DoubleVec<Matd> &meanvalue, DoubleVec<Matd> &variance)
	{
		int count = 0;
		std::cout << "The current length difference is " << diff << "." << endl;
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i)
		{
			for (int j = 0; j != this->j_; ++j)
			{
				for (int index_i = 0; index_i != meanvalue[0][0].size(); ++index_i)
				{
					for (int index_j = 0; index_j != meanvalue[0][0].size(); ++index_j)
					{
						Real relative_value_ = (std::pow(current_result[i][j][index_i][index_j] - meanvalue[i + diff][j][index_i][index_j], 2) - variance[i + diff][j][index_i][index_j]) / variance[i + diff][j][index_i][index_j];
						if (relative_value_ > 0.01)
						{
							std::cout << this->quantity_name_ << "[" << j << "][" << index_i << "] in " << this->element_tag_[i] << " is beyond the exception, and difference is "
								<< relative_value_ << endl;
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
	void RegressionTestEnsembleAveraged<ObserveMethodType>::settingupAndCorrection()
	{
		this->i_ = this->current_result_.size();
		this->j_ = this->current_result_[0].size();

		if (this->number_of_run_ > 1)
		{
			if (this->converged == "false" ) /*< To identify the database generation or new result testing. */
			{
				if (!fs::exists(this->result_filefullpath_))
				{
					std::cout << "\n Error: the input file:" << this->result_filefullpath_ << " is not exists" << std::endl;
					std::cout << __FILE__ << ':' << __LINE__ << std::endl;
					exit(1);
				}
				else
					this->result_xml_engine_in_.loadXmlFile(this->result_filefullpath_);
			}

			if (!fs::exists(this->mean_variance_filefullpath_))
			{
				std::cout << "\n Error: the input file:" << this->mean_variance_filefullpath_ << " is not exists" << std::endl;
				std::cout << __FILE__ << ':' << __LINE__ << std::endl;
				exit(1);
			}
			else
			{
				this->mean_variance_xml_engine_in_.loadXmlFile(this->mean_variance_filefullpath_);
				SimTK::Xml::Element mean_element_ = this->mean_variance_xml_engine_in_.getChildElement("Mean_Element");
				this->number_of_snapshot_old_ = std::distance(mean_element_.element_begin(), mean_element_.element_end());

				DoubleVec<VariableType> temp(SMAX(this->i_, this->number_of_snapshot_old_), StdVec<VariableType>(this->j_));
				meanvalue_ = temp;
				variance_ = temp;

				/** Unify the length of current result and previous result. */
				if (this->number_of_snapshot_old_ < this->i_)
				{
					this->difference_ = this->i_ - this->number_of_snapshot_old_;
					for (int delete_ = 0; delete_ != this->difference_; ++delete_)
						this->current_result_.pop_back();
				}
				else if (this->number_of_snapshot_old_ > this->i_)
					this->difference_ = this->number_of_snapshot_old_ - this->i_;
				else
					this->difference_ = 0;
			}
		}
		else if (this->number_of_run_ == 1)
		{
			this->number_of_snapshot_old_ = this->i_;
			DoubleVec<VariableType> temp(this->i_, StdVec<VariableType>(this->j_));
			this->result_.push_back(this->current_result_);
			meanvalue_ = temp;
			variance_ = temp;
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestEnsembleAveraged<ObserveMethodType>::readMeanVarianceFromXml()
	{
		if (this->number_of_run_ > 1)
		{
			SimTK::Xml::Element mean_element_ = this->mean_variance_xml_engine_in_.getChildElement("Mean_Element");
			SimTK::Xml::Element variance_element_ = this->mean_variance_xml_engine_in_.getChildElement("Variance_Element");
			for (int particle_n_ = 0; particle_n_ != this->j_; ++particle_n_)
			{
				this->xmlmemory_io_.readDataFromXmlMemory(this->mean_variance_xml_engine_in_, 
					mean_element_, particle_n_, this->meanvalue_, this->quantity_name_);
				this->xmlmemory_io_.readDataFromXmlMemory(this->mean_variance_xml_engine_in_,
					variance_element_, particle_n_, this->variance_, this->quantity_name_);
			}
		}
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestEnsembleAveraged<ObserveMethodType>::updateMeanVariance()
	{
		/** Unify the length of result and meanvalue. */
		if (this->number_of_run_ > 1)
		{
			for (int delete_ = 0; delete_ != this->difference_; ++delete_)
			{
				meanvalue_.pop_back(); 
				variance_.pop_back();
			}
		}
		meanvalue_new_ = meanvalue_;
		variance_new_ = variance_;

		/** update the meanvalue of the result. */
		for (int i = 0; i != SMIN(this->i_, this->number_of_snapshot_old_); ++i)
			for (int j = 0; j != this->j_; ++j)
				meanvalue_new_[i][j] = (meanvalue_[i][j] * (this->number_of_run_ - 1) + this->current_result_[i][j]) / this->number_of_run_;
		/** Update the variance of the result. */
		calculateNewVariance(this->result_, meanvalue_new_, variance_, variance_new_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	void RegressionTestEnsembleAveraged<ObserveMethodType>::writeMeanVarianceToXml()
	{
		this->mean_variance_xml_engine_out_.addElementToXmlDoc("Mean_Element");
		SimTK::Xml::Element mean_element_ = this->mean_variance_xml_engine_out_.getChildElement("Mean_Element");
		this->xmlmemory_io_.writeDataToXmlMemory(this->mean_variance_xml_engine_out_, mean_element_, this->meanvalue_new_,
			SMIN(this->i_, this->number_of_snapshot_old_), this->j_, this->quantity_name_, this->element_tag_);
		this->mean_variance_xml_engine_out_.addElementToXmlDoc("Variance_Element");
		SimTK::Xml::Element variance_element_ = this->mean_variance_xml_engine_out_.getChildElement("Variance_Element");
		this->xmlmemory_io_.writeDataToXmlMemory(this->mean_variance_xml_engine_out_, variance_element_, this->variance_new_,
			SMIN(this->i_, this->number_of_snapshot_old_), this->j_, this->quantity_name_, this->element_tag_);
		this->mean_variance_xml_engine_out_.writeToXmlFile(this->mean_variance_filefullpath_);
	};
	//=================================================================================================//
	template<class ObserveMethodType>
	bool RegressionTestEnsembleAveraged<ObserveMethodType>::compareMeanVariance()
	{
		int count_not_converged_m = 0;
		int count_not_converged_v = 0;
		count_not_converged_m = compareParameter("meanvalue", meanvalue_, meanvalue_new_, this->threshold_mean_);
		count_not_converged_v = compareParameter("variance", variance_, variance_new_, this->threshold_variance_);
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
	void RegressionTestEnsembleAveraged<ObserveMethodType>::resultTest()
	{
		/* compare the current result to the converged mean value and variance. */
		int test_wrong = 0;
		if (this->i_ < this->number_of_snapshot_old_)
			test_wrong = testNewResult(this->difference_, this->current_result_, meanvalue_, variance_);
		else
		{
			/** Unify the length of meanvalue, variance, old result, new result. */
			for (int delete_ = 0; delete_ != this->difference_; ++delete_)
			{
				meanvalue_.pop_back();
				variance_.pop_back();
			}
			test_wrong = testNewResult(0, this->current_result_, meanvalue_, variance_);
		}
		/* draw the conclusion. */
		if (test_wrong == 0)
			std::cout << "The result of " << this->quantity_name_ << " are correct based on the ensemble averaged regression test!" << endl;
		else
		{
			std::cout << "There are " << test_wrong << " snapshots are not within the expected range." << endl;
			std::cout << "Please try again. If it still post this conclusion, the result is not correct!" << endl;
			exit(1);
		}
	};
	//=================================================================================================//
}