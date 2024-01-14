#include "k-epsilon_turbulent_model_inner.hpp"
#include "k-epsilon_turbulent_model_inner.h"
namespace SPH
{
	//=================================================================================================//
	namespace fluid_dynamics
	{
		//=================================================================================================//
		BaseTurbuClosureCoeffInner::BaseTurbuClosureCoeffInner()
			: Karman(0.4187), C_mu(0.09), TurbulentIntensity(5.0e-2), sigma_k(1.0),
			C_l(1.44), C_2(1.92), sigma_E(1.3), turbu_const_E(9.793){}
		//=================================================================================================//
		BaseTurtbulentModelInner::BaseTurtbulentModelInner(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
			particle_spacing_min_(inner_relation.real_body_->sph_adaptation_->MinimumSpacing()),
			 rho_(particles_->rho_), vel_(particles_->vel_), 
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()), dimension_(Vecd(0).size()),
			smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()) {}
		//=================================================================================================//
		K_TurtbulentModelInner::K_TurtbulentModelInner(BaseInnerRelation& inner_relation, const StdVec<Real>& initial_values)
			: BaseTurtbulentModelInner(inner_relation)
		{
			//** Obtain Initial values for transport equations *
			turbu_k_initial_ = initial_values[0];
			turbu_ep_initial_ = initial_values[1];
			turbu_mu_initial_ = initial_values[2];

			particles_->registerVariable(dk_dt_, "ChangeRateOfTKE");
			particles_->registerSortableVariable<Real>("ChangeRateOfTKE");

			//particles_->registerVariable(turbu_k_, "TurbulenceKineticEnergy", 0.000180001);
			particles_->registerVariable(turbu_k_, "TurbulenceKineticEnergy", turbu_k_initial_);

			particles_->registerSortableVariable<Real>("TurbulenceKineticEnergy");
			particles_->addVariableToWrite<Real>("TurbulenceKineticEnergy");

			//particles_->registerVariable(turbu_mu_, "TurbulentViscosity", 1.0e-9);
			particles_->registerVariable(turbu_mu_, "TurbulentViscosity", turbu_mu_initial_);
			particles_->registerSortableVariable<Real>("TurbulentViscosity");
			particles_->addVariableToWrite<Real>("TurbulentViscosity");

			//particles_->registerVariable(turbu_epsilon_, "TurbulentDissipation", 3.326679e-5);
			particles_->registerVariable(turbu_epsilon_, "TurbulentDissipation", turbu_ep_initial_);

			particles_->registerSortableVariable<Real>("TurbulentDissipation");
			particles_->addVariableToWrite<Real>("TurbulentDissipation");

			particles_->registerVariable(k_production_, "K_Production");
			particles_->registerSortableVariable<Real>("K_Production");
			particles_->addVariableToWrite<Real>("K_Production");

			//particles_->registerVariable(B_, "CorrectionMatrix");
			//particles_->registerSortableVariable<Matd>("CorrectionMatrix");
			//particles_->addVariableToWrite<Matd>("CorrectionMatrix");

			particles_->registerVariable(is_near_wall_P1_, "IsNearWallP1");
			particles_->registerSortableVariable<int>("IsNearWallP1");
			particles_->addVariableToWrite<int>("IsNearWallP1");

			//** for test */
			particles_->registerVariable(k_diffusion_, "K_Diffusion");
			particles_->registerSortableVariable<Real>("K_Diffusion");
			particles_->addVariableToWrite<Real>("K_Diffusion");

			particles_->addVariableToWrite<Real>("ChangeRateOfTKE");

			particles_->registerVariable(velocity_gradient_, "VelocityGradient");
			particles_->registerSortableVariable<Matd>("VelocityGradient");
			particles_->addVariableToWrite<Matd>("VelocityGradient");

			particles_->registerVariable(vel_x_, "Velocity_X");
			particles_->registerSortableVariable<Real>("Velocity_X");


		}
		//=================================================================================================//
		BaseGetTimeAverageData::BaseGetTimeAverageData(BaseInnerRelation& inner_relation,int num_observer_points)
			: BaseTurtbulentModelInner(inner_relation), pos_(particles_->pos_), num_cell(num_observer_points),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation")), plt_engine_()
		{
			num_data = 5;
			file_name_.push_back("vel_x_sto_");
			file_name_.push_back("turbu_k_sto_");
			file_name_.push_back("turbu_epsilon_sto_");
			file_name_.push_back("turbu_mu_sto_");
			file_name_.push_back("vel_sto_");

			num_in_cell_.resize(num_cell);
			data_time_aver_sto_.resize(num_cell); //Rows

			data_sto_.resize(num_cell); //Rows
			for (size_t i = 0; i != num_cell; ++i)
			{
				data_sto_[i].resize(num_data); //Cols
			}

			for (size_t j = 0; j != num_data; ++j)
			{
				file_path_output_ = "../bin/output/" + file_name_[j] + ".dat";
				std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
				out_file << "run_time" << "   ";
				for (size_t i = 0; i != num_cell; ++i)
				{
					std::string quantity_name_i = file_name_[j] + "[" + std::to_string(i) + "]";
					plt_engine_.writeAQuantityHeader(out_file, data_sto_[i][j], quantity_name_i);
				}
				out_file << "\n";
				out_file.close();
			}
		}
		//=================================================================================================//
		void BaseGetTimeAverageData::output_time_history_data(Real cutoff_time)
		{
			/** Output for .dat file. */
			for (size_t j = 0; j != num_data; ++j)
			{
				file_path_output_ = "../bin/output/" + file_name_[j] + ".dat";
				std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
				out_file << GlobalStaticVariables::physical_time_ << "   ";
				for (size_t i = 0; i != num_cell; ++i)
				{
					//if (num_in_cell_[i] == 0 && GlobalStaticVariables::physical_time_ > cutoff_time)
					//{
					//	std::cout << "There is a empaty monitoring cell, cell number=" << i << std::endl;
					//	system("pause");
					//}
					num_in_cell_[i] == 0 ? plt_engine_.writeAQuantity(out_file, 0.0) :
						plt_engine_.writeAQuantity(out_file, data_sto_[i][j] / num_in_cell_[i]);
				}
				out_file << "\n";
				out_file.close();
			}
			//** Clear data *
			for (int i = 0; i < num_cell; i++)
			{
				num_in_cell_[i] = 0;
				for (size_t j = 0; j != num_data; ++j)
				{
					data_sto_[i][j] = 0.0;
				}
			}
		}
		//=================================================================================================//
		void BaseGetTimeAverageData::get_time_average_data(Real cutoff_time)
		{
			for (size_t j = 0; j != num_data; ++j)
			{
				data_loaded_.clear();
				int num_line_data = 0;
				//** Load data *
				file_path_input_ = "../bin/output/" + file_name_[j] + ".dat";
				std::ifstream in_file(file_path_input_.c_str());
				bool skipFirstLine = true;
				std::string line;
				while (std::getline(in_file, line))
				{
					if (skipFirstLine)
					{
						skipFirstLine = false;
						continue;
					}
					num_line_data++;
					std::vector<Real> data_point;
					std::istringstream iss(line);
					Real value;
					while (iss >> value)
					{
						data_point.push_back(value);
					}
					data_loaded_.push_back(data_point);
				}

				in_file.close();
				//** Deal with data *
				for (size_t k = 0; k != num_cell; ++k)
				{
					Real sum = 0.0;
					int count = 0;
					for (size_t i = 0; i != num_line_data; ++i)
					{
						if (data_loaded_[i][0] > cutoff_time)
						{
							count++;
							Real delta_t = data_loaded_[i][0] - data_loaded_[i - 1][0];
							sum += data_loaded_[i][k + 1] * delta_t;
							//sum += data_loaded_[i][k + 1]; //**the first col is time*
						}
					}
					//data_time_aver_sto_[k] = sum / count;
					data_time_aver_sto_[k] = sum / (data_loaded_[num_line_data-1][0] - cutoff_time);
				}
				//** Output data *
				file_path_output_ = "../bin/output/TimeAverageData.dat";
				std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
				out_file << file_name_[j] << "\n";
				for (size_t k = 0; k != num_cell; ++k)
				{
					plt_engine_.writeAQuantity(out_file, data_time_aver_sto_[k]);
				}
				out_file << "\n";
				out_file.close();
			}
			std::cout << "The cutoff_time is " << cutoff_time << std::endl;
		}
		//=================================================================================================//
		GetTimeAverageCrossSectionData::GetTimeAverageCrossSectionData(BaseInnerRelation& inner_relation, int num_observer_points, const StdVec<Real>& bound_x, Real offset_dist_y)
			: BaseGetTimeAverageData(inner_relation, num_observer_points)
		{
			x_min_ = bound_x[0];
			x_max_ = bound_x[1];
			offset_dist_y_ = offset_dist_y;
			//** Get the center coordinate of the monitoring cell *
			for (int i = 0; i < num_cell; i++)
			{
				Real upper_bound = ((i + 1) * particle_spacing_min_ + offset_dist_y_);
				Real lower_bound = (i * particle_spacing_min_ + offset_dist_y_);
				monitor_cellcenter_y.push_back((lower_bound + upper_bound) / 2.0);
			}
			file_path_output_ = "../bin/output/monitor_cell_center_y.dat";
			std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
			for (size_t i = 0; i != num_cell; ++i)
			{
				plt_engine_.writeAQuantity(out_file, monitor_cellcenter_y[i]);
				out_file << "\n";
			}
			out_file << "\n";
			out_file.close();
		}
		//=================================================================================================//
		void GetTimeAverageCrossSectionData::update(size_t index_i, Real dt)
		{
			//** Get data *
			if (pos_[index_i][0] > x_min_ && pos_[index_i][0] <= x_max_)
			{
				for (int i = 0; i < num_cell; i++)
				{
					if (pos_[index_i][1] > (i* particle_spacing_min_+ offset_dist_y_) &&
						pos_[index_i][1] <= ((i+1)* particle_spacing_min_+ offset_dist_y_))
					{
						num_in_cell_[i] += 1;
						data_sto_[i][0] += vel_[index_i][0];
						data_sto_[i][1] += turbu_k_[index_i];
						data_sto_[i][2] += turbu_epsilon_[index_i];
						data_sto_[i][3] += turbu_mu_[index_i];
						data_sto_[i][4] += vel_[index_i].norm();
					}
				}
			}
		}
		//=================================================================================================//
		GetTimeAverageCenterLineData::GetTimeAverageCenterLineData(BaseInnerRelation& inner_relation, 
			int num_observer_points, Real observe_x_ratio, const StdVec<Real>& bound_y, const StdVec<Real>& bound_x_f, const StdVec<Real>& bound_x_b)
			: BaseGetTimeAverageData(inner_relation, num_observer_points), observe_x_ratio_(observe_x_ratio),
			bound_x_f_(bound_x_f),bound_x_b_(bound_x_b), bound_y_(bound_y)
		{
			observe_x_spacing_ = particle_spacing_min_ * observe_x_ratio_;
		}
		//=================================================================================================//
		void GetTimeAverageCenterLineData::update(size_t index_i, Real dt)
		{
			//** Get data *
			if (pos_[index_i][1] > bound_y_[0] && pos_[index_i][1] <= bound_y_[1])
			{
				for (int i = 0; i < num_cell; i++) 
				{ 
					if (i < bound_x_f_.size() - 1) //* Front of cylinder
					{
						if (pos_[index_i][0] > bound_x_f_[i] && pos_[index_i][0] <= bound_x_f_[i + 1])
						{
							num_in_cell_[i] += 1;
							data_sto_[i][0] += vel_[index_i][0];
							data_sto_[i][1] += turbu_k_[index_i];
							data_sto_[i][2] += turbu_epsilon_[index_i];
							data_sto_[i][3] += turbu_mu_[index_i];
							data_sto_[i][4] += vel_[index_i].norm();
						}
					}
					else if (i >= bound_x_f_.size() - 1) //* behind of cylinder
					{
						int j = i - (bound_x_f_.size() - 1);
						if (pos_[index_i][0] > bound_x_b_[j] && pos_[index_i][0] <= bound_x_b_[j + 1])
						{
							num_in_cell_[i] += 1;
							data_sto_[i][0] += vel_[index_i][0];
							data_sto_[i][1] += turbu_k_[index_i];
							data_sto_[i][2] += turbu_epsilon_[index_i];
							data_sto_[i][3] += turbu_mu_[index_i];
							data_sto_[i][4] += vel_[index_i].norm();
						}
					}
				}

			}
		}
		//=================================================================================================//
		void GetTimeAverageCenterLineData::output_monitor_x_coordinate()
		{
			StdVec<Real> monitor_cellcenter_x;
			for (int i = 0; i < bound_x_f_.size() - 1; i++)
			{
				monitor_cellcenter_x.push_back((bound_x_f_[i] + bound_x_f_[i + 1]) / 2.0);
			}
			for (int i = 0; i < bound_x_b_.size() - 1; i++)
			{
				monitor_cellcenter_x.push_back((bound_x_b_[i] + bound_x_b_[i + 1]) / 2.0);
			}

			file_path_output_ = "../bin/output/monitor_cell_center_x.dat";
			std::ofstream out_file(file_path_output_.c_str(), std::ios::app);
			for (size_t i = 0; i != num_cell; ++i)
			{
				plt_engine_.writeAQuantity(out_file, monitor_cellcenter_x[i]);
				out_file << "\n";
			}
			out_file << "\n";
			out_file.close();
		}

		//=================================================================================================//
		//GetVelocityGradientInner::GetVelocityGradientInner(BaseInnerRelation& inner_relation)
		//	: LocalDynamics(inner_relation.getSPHBody()), FluidDataInner(inner_relation),
		//	vel_(particles_->vel_), pos_(particles_->pos_),
		//	velocity_gradient_(*particles_->getVariableByName<Matd>("VelocityGradient")),
		//	is_near_wall_P1_(*particles_->getVariableByName<int>("IsNearWallP1"))
		//{
		//	//for test
		//	particles_->registerVariable(velocity_gradient_wall, "Velocity_Gradient_Wall");
		//	particles_->registerSortableVariable<Matd>("Velocity_Gradient_Wall");
		//	particles_->addVariableToWrite<Matd>("Velocity_Gradient_Wall");
		//}
		//=================================================================================================//
		void K_TurtbulentModelInner::update(size_t index_i, Real dt)
		{
			turbu_k_[index_i] += dk_dt_[index_i] * dt;
		}
		//=================================================================================================//
		E_TurtbulentModelInner::E_TurtbulentModelInner(BaseInnerRelation& inner_relation)
			: BaseTurtbulentModelInner(inner_relation), 
			k_production_(*particles_->getVariableByName<Real>("K_Production")),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation"))
		{
			particles_->registerVariable(dE_dt_, "ChangeRateOfTDR");
			particles_->registerSortableVariable<Real>("ChangeRateOfTDR");
			particles_->addVariableToWrite<Real>("ChangeRateOfTDR");

			particles_->registerVariable(ep_production, "Ep_Production");
			particles_->registerSortableVariable<Real>("Ep_Production");
			particles_->addVariableToWrite<Real>("Ep_Production");
			particles_->registerVariable(ep_dissipation_, "Ep_Dissipation_");
			particles_->registerSortableVariable<Real>("Ep_Dissipation_");
			particles_->addVariableToWrite<Real>("Ep_Dissipation_");
			particles_->registerVariable(ep_diffusion_, "Ep_Diffusion_");
			particles_->registerSortableVariable<Real>("Ep_Diffusion_");
			particles_->addVariableToWrite<Real>("Ep_Diffusion_");

		}
		//=================================================================================================//
		void E_TurtbulentModelInner::update(size_t index_i, Real dt)
		{
			turbu_epsilon_[index_i] += dE_dt_[index_i] * dt;
		}
		//=================================================================================================//
		//TKEnergyAccInner::
		//	TKEnergyAccInner(BaseInnerRelation& inner_relation)
		//	: BaseTurtbulentModelInner(inner_relation), acc_prior_(particles_->acc_prior_),
		//	indicator_(particles_->indicator_), pos_(particles_->pos_),
		//	turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy"))
		//{
		//	particles_->registerVariable(tke_acc_inner_, "TkeAccInner");
		//	particles_->addVariableToWrite<Vecd>("TkeAccInner");
		//	particles_->registerVariable(tke_acc_wall_, "TkeAccWall");
		//	particles_->addVariableToWrite<Vecd>("TkeAccWall");

		//	particles_->registerVariable(test_k_grad_rslt_, "TkeGradResult");
		//	particles_->addVariableToWrite<Vecd>("TkeGradResult");
		//}
		//=================================================================================================//
		//TurbuViscousAccInner::TurbuViscousAccInner(BaseInnerRelation& inner_relation)
		//	: BaseViscousAccelerationInner(inner_relation),
		//	turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
		//	wall_Y_plus_(*particles_->getVariableByName<Real>("WallYplus")),
		//	velo_friction_(*particles_->getVariableByName<Vecd>("FrictionVelocity")),
		//	y_p_(*particles_->getVariableByName<Real>("Y_P"))
		//{
		//	particles_->registerVariable(visc_acc_inner_, "ViscousAccInner");
		//	particles_->addVariableToWrite<Vecd>("ViscousAccInner");
		//	particles_->registerVariable(visc_acc_wall_, "ViscousAccWall");
		//	particles_->addVariableToWrite<Vecd>("ViscousAccWall");
		//}
		//=================================================================================================//
		TurbulentEddyViscosity::
			TurbulentEddyViscosity(SPHBody& sph_body)
			: LocalDynamics(sph_body), FluidDataSimple(sph_body),
			rho_(particles_->rho_), wall_Y_star_(*particles_->getVariableByName<Real>("WallYstar")),
			wall_Y_plus_(*particles_->getVariableByName<Real>("WallYplus")),
			mu_(DynamicCast<Fluid>(this, particles_->getBaseMaterial()).ReferenceViscosity()),
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation")){}
		//=================================================================================================//
		void TurbulentEddyViscosity::update(size_t index_i, Real dt)
		{
			turbu_mu_[index_i] = rho_[index_i] * C_mu * turbu_k_[index_i] * turbu_k_[index_i] / (turbu_epsilon_[index_i]);
		}
		//=================================================================================================//
		TurbulentAdvectionTimeStepSize::TurbulentAdvectionTimeStepSize(SPHBody& sph_body, Real U_max, Real advectionCFL)
			: LocalDynamicsReduce<Real, ReduceMax>(sph_body, U_max* U_max),FluidDataSimple(sph_body), 
			vel_(particles_->vel_), advectionCFL_(advectionCFL),
			smoothing_length_min_(sph_body.sph_adaptation_->MinimumSmoothingLength()),
			fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
			turbu_mu_(*particles_->getVariableByName<Real>("TurbulentViscosity"))
		{
			Real viscous_speed = fluid_.ReferenceViscosity() / fluid_.ReferenceDensity() / smoothing_length_min_;
			reference_ = SMAX(viscous_speed * viscous_speed, reference_);
		}
		//=================================================================================================//
		Real TurbulentAdvectionTimeStepSize::reduce(size_t index_i, Real dt)
		{
			Real turbu_viscous_speed = (fluid_.ReferenceViscosity() + turbu_mu_[index_i]) 
				/ fluid_.ReferenceDensity() / smoothing_length_min_;
			Real turbu_viscous_speed_squre = turbu_viscous_speed * turbu_viscous_speed;
			Real vel_n_squre = vel_[index_i].squaredNorm();
			Real vel_bigger = SMAX(turbu_viscous_speed_squre, vel_n_squre);

			return vel_bigger;
		}
		//=================================================================================================//
		Real TurbulentAdvectionTimeStepSize::outputResult(Real reduced_value)
		{
			Real speed_max = sqrt(reduced_value);
			return advectionCFL_ * smoothing_length_min_ / (speed_max + TinyReal);
		}
		//=================================================================================================//
		InflowTurbulentCondition::InflowTurbulentCondition(BodyPartByCell& body_part
			, Real CharacteristicLength, Real relaxation_rate): 
			BaseFlowBoundaryCondition(body_part),
			relaxation_rate_(relaxation_rate), 
			CharacteristicLength_(CharacteristicLength), 
			turbu_k_(*particles_->getVariableByName<Real>("TurbulenceKineticEnergy")),
			turbu_epsilon_(*particles_->getVariableByName<Real>("TurbulentDissipation"))
		{
			TurbulentLength_ = 0.07 * CharacteristicLength_ / pow(C_mu, 0.75);
		}
		//=================================================================================================//
		void InflowTurbulentCondition::update(size_t index_i, Real dt)
		{
			Real target_in_turbu_k = getTurbulentInflowK(pos_[index_i], vel_[index_i], turbu_k_[index_i]);
			turbu_k_[index_i] += relaxation_rate_ * (target_in_turbu_k - turbu_k_[index_i]);
			Real target_in_turbu_E = getTurbulentInflowE(pos_[index_i], turbu_k_[index_i], turbu_epsilon_[index_i]);
			turbu_epsilon_[index_i] += relaxation_rate_ * (target_in_turbu_E - turbu_epsilon_[index_i]);
		}
		//=================================================================================================//
		Real InflowTurbulentCondition:: getTurbulentInflowK(Vecd& position, Vecd& velocity, Real& turbu_k)
		{
			Real u = velocity[0];
			Real temp_in_turbu_k = 1.5 * pow((TurbulentIntensity * u), 2);
			Real turbu_k_original = turbu_k;
			if (position[0] < 0.0)
			{
				turbu_k_original = temp_in_turbu_k;
				//std::cout << "temp_in_turbu_k="<< temp_in_turbu_k << std::endl;
			}
			return turbu_k_original;
		}
		//=================================================================================================//
		Real InflowTurbulentCondition::getTurbulentInflowE(Vecd& position, Real& turbu_k, Real& turbu_E)
		{
			//Real temp_in_turbu_E = C_mu * pow(turbu_k, 1.5) / (0.1*getTurbulentLength());
			Real temp_in_turbu_E = pow(turbu_k, 1.5) / TurbulentLength_;
			Real turbu_E_original = turbu_E;
			if (position[0] < 0.0)
			{
				turbu_E_original = temp_in_turbu_E;
			}
			return turbu_E_original;
		}
		//=================================================================================================//
		ClearYPositionForTest::
			ClearYPositionForTest(SPHBody& sph_body)
			: LocalDynamics(sph_body), FluidDataSimple(sph_body),
			pos_(particles_->pos_), vel_(particles_->vel_) {}
		//=================================================================================================//
		void ClearYPositionForTest::update(size_t index_i, Real dt)
		{
			vel_[index_i][1] = 0.0;
		}
	}
	//=================================================================================================//
}
//=================================================================================================//