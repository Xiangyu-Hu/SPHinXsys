#include "granular_dynamics_inner.h"
#include "granular_dynamics_inner.hpp"
#include <numeric>
namespace SPH
{
	//=================================================================================================//
	namespace granular_dynamics
	{
		GranularInitialCondition::GranularInitialCondition(SPHBody& sph_body)
			: LocalDynamics(sph_body), GranularDataSimple(sph_body),
			pos_(particles_->pos_), vel_(particles_->vel_) {}
		//=================================================================================================//
		GranularAcousticTimeStepSize::GranularAcousticTimeStepSize(SPHBody& sph_body, Real acousticCFL)
			: fluid_dynamics::AcousticTimeStepSize(sph_body, acousticCFL) {}
		//=================================================================================================//
		Real GranularAcousticTimeStepSize::reduce(size_t index_i, Real dt)
		{
			return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
		}
		//=================================================================================================//
		Real GranularAcousticTimeStepSize::outputResult(Real reduced_value)
		{
			return acousticCFL_ * smoothing_length_min_ / (fluid_.ReferenceSoundSpeed() + TinyReal);
		}
		//=================================================================================================//
		BaseRelaxation::BaseRelaxation(BaseInnerRelation& inner_relation)
			: LocalDynamics(inner_relation.getSPHBody()), GranularDataInner(inner_relation),
			granular_material_(particles_->granular_material_), rho_(particles_->rho_),
			p_(particles_->p_), drho_dt_(particles_->drho_dt_),
			pos_(particles_->pos_), vel_(particles_->vel_),
			acc_(particles_->acc_), acc_prior_(particles_->acc_prior_) {}
		//=================================================================================================//
		//===============================ArtificialStressRelaxation======================================//
		//=================================================================================================//
		BaseArtificialStressRelaxation ::
			BaseArtificialStressRelaxation(BaseInnerRelation& inner_relation, Real epsilon)
			: BaseRelaxation(inner_relation), smoothing_length_(sph_body_.sph_adaptation_->ReferenceSmoothingLength()),
			reference_spacing_(sph_body_.sph_adaptation_->ReferenceSpacing()), epsilon_(epsilon) { }
		Matd BaseArtificialStressRelaxation::repulsiveForce(Matd stress_tensor_i, Real rho_i)
		{
			Real sigma_xx = stress_tensor_i(0, 0);
			Real sigma_xy = stress_tensor_i(0, 1);
			Real sigma_yy = stress_tensor_i(1, 1);
			Real tiny_real(0);
			sigma_xx - sigma_yy > 0 ? tiny_real = TinyReal : tiny_real = -TinyReal;
			Real tan_sita_2 = 2 * sigma_xy / (sigma_xx - sigma_yy + tiny_real);
			Real sita_2 = atan(tan_sita_2);
			Real sita = sita_2 / 2.0;
			Real sigma_xx_dot = cos(sita) * cos(sita) * sigma_xx + 2 * cos(sita) * sin(sita) * sigma_xy + sin(sita) * sin(sita) * sigma_yy;
			Real sigma_yy_dot = sin(sita) * sin(sita) * sigma_xx + 2 * cos(sita) * sin(sita) * sigma_xy + cos(sita) * cos(sita) * sigma_yy;
			Real R_xx_dot = 0;
			Real R_yy_dot = 0;
			if (sigma_xx_dot > 0)
			{
				R_xx_dot = -epsilon_ * sigma_xx_dot / (rho_i * rho_i);
			}
			if (sigma_yy_dot > 0)
			{
				R_yy_dot = -epsilon_ * sigma_yy_dot / (rho_i * rho_i);
			}
			Matd R = Matd::Zero();
			R(0, 0) = R_xx_dot * cos(sita) * cos(sita) + R_yy_dot * sin(sita) * sin(sita);
			R(1, 1) = R_xx_dot * sin(sita) * sin(sita) + R_yy_dot * cos(sita) * cos(sita);
			R(0, 1) = (R_xx_dot - R_yy_dot) * cos(sita) * sin(sita);
			R(1, 0) = R(0, 1);
			return R;
		}
		ArtificialNormalShearStressRelaxation ::
			ArtificialNormalShearStressRelaxation(BaseInnerRelation& inner_relation, Real exponent)
			: BaseArtificialStressRelaxation(inner_relation), shear_stress_(particles_->shear_stress_), acc_shear_(particles_->acc_shear_), exponent_(exponent) {}
		void ArtificialNormalShearStressRelaxation::interaction(size_t index_i, Real dt)
		{
			Vecd acceleration = Vecd::Zero();
			Real rho_i = rho_[index_i];
			Matd stress_i = shear_stress_[index_i] - p_[index_i] * Matd::Identity();
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			Real W_ini = sph_body_.sph_adaptation_->getKernel()->W_2D(reference_spacing_ / smoothing_length_);
			Matd R_i = repulsiveForce(stress_i, rho_i);
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Real r_ij = inner_neighborhood.r_ij_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

				Real W_ij = sph_body_.sph_adaptation_->getKernel()->W_2D(r_ij / smoothing_length_);
				Real f_ij = W_ij / W_ini;
				Matd stress_j = shear_stress_[index_j] - p_[index_j] * Matd::Identity();
				Matd R_j = repulsiveForce(stress_j, rho_[index_j]);
				Matd repulsive_force = pow(f_ij, exponent_) * (R_i + R_j);

				acceleration += rho_[index_j] * repulsive_force * nablaW_ijV_j;
			}
			acc_shear_[index_i] = acceleration;
		}

		//=================================================================================================//
		//===============================ShearStressRelaxation1stHalf======================================//
		//=================================================================================================/
		ShearStressRelaxation1stHalf ::
			ShearStressRelaxation1stHalf(BaseInnerRelation& inner_relation)
			: BaseRelaxation(inner_relation),
			shear_stress_(particles_->shear_stress_), shear_stress_rate_(particles_->shear_stress_rate_),
			acc_shear_(particles_->acc_shear_) {}
		//=================================================================================================//
		void ShearStressRelaxation1stHalf::initialization(size_t index_i, Real dt)
		{
			shear_stress_[index_i] += shear_stress_rate_[index_i] * dt * 0.5;
		}
		//=================================================================================================//
		void ShearStressRelaxation1stHalf::interaction(size_t index_i, Real dt)
		{
			Real rho_i = rho_[index_i];
			Matd shear_stress_i = shear_stress_[index_i];
			Vecd acceleration = Vecd::Zero();
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
				acceleration += rho_[index_j] * (shear_stress_i / (rho_i * rho_i) + shear_stress_[index_j] / (rho_[index_j] * rho_[index_j])) * nablaW_ijV_j;
			}
			acc_shear_[index_i] += acceleration;  //for with artificial stress
		}

		void ShearStressRelaxation1stHalf::update(size_t index_i, Real dt)
		{
			vel_[index_i] += acc_shear_[index_i] * dt;
		}
		//=================================================================================================//
		//===============================ShearStressRelaxation2ndHalf======================================//
		//=================================================================================================//
		ShearStressRelaxation2ndHalf ::
			ShearStressRelaxation2ndHalf(BaseInnerRelation& inner_relation)
			: BaseRelaxation(inner_relation),
			shear_stress_(particles_->shear_stress_), shear_stress_rate_(particles_->shear_stress_rate_),
			velocity_gradient_(particles_->velocity_gradient_), strain_tensor_(particles_->strain_tensor_),
			strain_tensor_rate_(particles_->strain_tensor_rate_), von_mises_stress_(particles_->von_mises_stress_) {}
		//=================================================================================================//
		void ShearStressRelaxation2ndHalf::interaction(size_t index_i, Real dt)
		{
			Matd velocity_gradient = Matd::Zero();
			Neighborhood& inner_neighborhood = inner_configuration_[index_i];
			for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
			{
				size_t index_j = inner_neighborhood.j_[n];
				Vecd nablaW_ijV_j = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];

				Matd velocity_gradient_ij = -(vel_[index_i] - vel_[index_j]) * nablaW_ijV_j.transpose();
				velocity_gradient += velocity_gradient_ij;
			}
			velocity_gradient_[index_i] = velocity_gradient;
		}

		void ShearStressRelaxation2ndHalf::update(size_t index_i, Real dt)
		{
			shear_stress_rate_[index_i] = granular_material_.ConstitutiveRelationShearStress(velocity_gradient_[index_i], shear_stress_[index_i]);
			shear_stress_[index_i] += shear_stress_rate_[index_i] * dt * 0.5;
			Matd stress_tensor_i = shear_stress_[index_i] + p_[index_i] * Matd::Identity();
			von_mises_stress_[index_i] = particles_->getVonMisesStress(stress_tensor_i);
		}
	}
}

