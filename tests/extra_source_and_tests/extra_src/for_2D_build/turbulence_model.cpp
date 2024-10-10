#ifndef TURBULENCEMODEL_CPP
#define TURBULENCEMODEL_CPP
#include "turbulence_model.h"
#include "sphinxsys.h"
namespace SPH
{  
    namespace fluid_dynamics
    { 
    //=================================================================================================//
        BaseTurbulence::BaseTurbulence(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator)
            : BaseIntegration<DataDelegateInner>(inner_relation),
          mom_(this->particles_->template registerStateVariable<Vecd>("Momentum")),
          dmom_dt_(this->particles_->template registerStateVariable<Vecd>("MomentumChangeRate")),
          dmass_dt_(this->particles_->template registerStateVariable<Real>("MassChangeRate")),
          K_prod_p_(this->particles_->template registerStateVariable<Real>("TKEProductionInWallAdjCell")),
          K_prod_(this->particles_->template registerStateVariable<Real>("TKEProduction")),
          Eps_p_(this->particles_->template registerStateVariable<Real>("DissipationRateInWallAdjCell")),
          K_adv_(this->particles_->template registerStateVariable<Real>("TKEAdvection")),
          K_lap_(this->particles_->template registerStateVariable<Real>("TKELaplacian")),
          Eps_adv_(this->particles_->template registerStateVariable<Real>("DissipationAdvection")),
          Eps_lap_(this->particles_->template registerStateVariable<Real>("DissipationLaplacian")),
          Eps_prod_(this->particles_->template registerStateVariable<Real>("DissipationProd")),
          Eps_destruction_(this->particles_->template registerStateVariable<Real>("DissipationDestruction")),
          Tau_wall_(this->particles_->template registerStateVariable<Real>("WallShearStress")),
          C_mu_(0.09), sigma_k_(1.0), sigma_eps_(1.3), C1_eps_(1.44), C2_eps_(1.92),
          K_(this->particles_->template registerStateVariable<Real>("TKE")),
          Eps_(this->particles_->template registerStateVariable<Real>("Dissipation")),
          mu_t_(this->particles_->template registerStateVariable<Real>("TurblunetViscosity")),
          ghost_creator_(ghost_creator)
          {}
        //=================================================================================================//
            WallAdjacentCells::WallAdjacentCells(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
              : BaseTurbulence(inner_relation, ghost_creator), wall_normal_(this->particles_->template registerStateVariable<Vecd>("WallNormal")),
                wall_adjacent_cell_flag_(this->particles_->template registerStateVariable<Real>("FlagForWallAdjacentCells")),
                yp_(this->particles_->template registerStateVariable<Real>("WallNormalDistance")),
                corner_cell_flag_(this->particles_->template registerStateVariable<Real>("CornerCellFlag")), 
                boundary_type_(this->particles_->template registerStateVariable<Real>("BoundaryType")), ymax_(0.0), 
                bounds_(inner_relation.getSPHBody())
                {
                    walladjacentcellyp();
                }
        //=================================================================================================//
        void WallAdjacentCells::walladjacentcellyp()
        {
            wall_adjacent_index_.resize(particles_->ParticlesBound());
            wall_ghost_index_.resize(particles_->ParticlesBound());
            wall_eij_.resize(particles_->ParticlesBound());
            
            for (size_t boundary_type = 0; boundary_type < ghost_creator_.each_boundary_type_with_all_ghosts_index_.size(); ++boundary_type)
            {
                if (!ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].empty())
                {
                    for (size_t ghost_number = 0; ghost_number != ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type].size(); ++ghost_number)
                    {
                        size_t ghost_index = ghost_creator_.each_boundary_type_with_all_ghosts_index_[boundary_type][ghost_number];
                        size_t index_real = ghost_creator_.each_boundary_type_contact_real_index_[boundary_type][ghost_number];
                        wall_eij_[index_real] = ghost_creator_.each_boundary_type_with_all_ghosts_eij_[boundary_type][ghost_number];

                        if (boundary_type == 3)
                        {
                            wall_adjacent_index_[index_real] = index_real;
                            wall_ghost_index_[index_real] = ghost_index;
                            wall_adjacent_cell_flag_[index_real] = 1.0;
                            if ((pos_[index_real] - pos_[ghost_index]).dot(wall_eij_[index_real]) > ymax_)
                            {
                                ymax_ = (pos_[index_real] - pos_[ghost_index]).dot(wall_eij_[index_real]);
                            }
                        }
                        if (wall_adjacent_cell_flag_[index_real] == 1.0 && boundary_type != 3)
                        {
                            corner_cell_flag_[index_real] = 1.0;
                            boundary_type_[ghost_index] = boundary_type;
                        }
                    }
                }
            }
        }
        //=================================================================================================//
        void WallAdjacentCells::update(size_t index_i, Real dt)
        {
            Vecd lower_wall = {pos_[index_i][0], 0.0};
            Vecd upper_wall = {pos_[index_i][0], 2.0};
            Vecd lower_wall_normal = {0.0, 1.0};
            Vecd upper_wall_normal = {0.0, -1.0};
            /*
            BoundingBox bounds = bounds_.getSPHSystemBounds();
            Real channelheight = bounds.second_[1];
            Real halfwidth = 0.5 * channelheight;
            */ 

            bool lower_wall_condition = ((pos_[index_i] - lower_wall).dot(lower_wall_normal) <= 1.0 * ymax_);
            bool upper_wall_condition = ((pos_[index_i] - upper_wall).dot(upper_wall_normal) <= 1.0 * ymax_);

            if (lower_wall_condition)
            {
                yp_[index_i] = (pos_[index_i] - lower_wall).dot(lower_wall_normal);
                wall_normal_[index_i] = lower_wall_normal;
            }
            else if (upper_wall_condition)
            {
                yp_[index_i] = (pos_[index_i] - upper_wall).dot(upper_wall_normal);   
                wall_normal_[index_i] = upper_wall_normal;
            }
        }
        //=================================================================================================//
        TurbuleceVariablesGradient::TurbuleceVariablesGradient(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
            : BaseTurbulence(inner_relation, ghost_creator), K_grad_(this->particles_->template registerStateVariable<Vecd>("TKEGradient")),
              Eps_grad_(this->particles_->template registerStateVariable<Vecd>("DissipationGradient"))
        {}
        //=================================================================================================//
        void TurbuleceVariablesGradient::update(size_t index_i, Real dt)
        {
            K_grad_[index_i] = Vecd::Zero(), Eps_grad_[index_i] = Vecd::Zero();
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Real dW_ij = inner_neighborhood.dW_ij_[n];
                Real r_ij = inner_neighborhood.r_ij_[n];
                Vecd &e_ij = inner_neighborhood.e_ij_[n];

                K_grad_[index_i] += dW_ij * Vol_[index_j] * (K_[index_i] - K_[index_j]) * e_ij;
                Eps_grad_[index_i] += dW_ij * Vol_[index_j] * (Eps_[index_i] - Eps_[index_j]) * e_ij;
                if (index_i == 41468)
                {
                    Vecd tkegrad = K_grad_[index_i];
                    Real ki = K_[index_i];
                    Real kj = K_[index_j];
                    Vecd epsgrad = Eps_grad_[index_i];
                    Real epsi = Eps_[index_i];
                    Real epsj = Eps_[index_j];
                    Real c = 0.0;
                }
                
            }
        }
        //=================================================================================================//

        KEpsilonStd1stHalf::KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator) 
            : StdWallFunctionFVM(inner_relation, ghost_creator),
              dK_dt_(this->particles_->template registerStateVariable<Real>("TKEChangeRate")),
              wall_adjacent_cell_flag_(this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells")),
              strain_rate_(this->particles_->template registerStateVariable<Real>("StrainRate")),
              corner_cell_flag_(this->particles_->template getVariableDataByName<Real>("CornerCellFlag")),
              boundary_type_(this->particles_->template getVariableDataByName<Real>("BoundaryType")),
              dudx_(this->particles_->template registerStateVariable<Real>("dudx")),
              dudy_(this->particles_->template registerStateVariable<Real>("dudy")),
              dvdx_(this->particles_->template registerStateVariable<Real>("dvdx")),
              dvdy_(this->particles_->template registerStateVariable<Real>("dvdy")),
              K_grad_(this->particles_->template getVariableDataByName<Vecd>("TKEGradient")),
              Eps_grad_(this->particles_->template getVariableDataByName<Vecd>("DissipationGradient"))
              {}
        //=================================================================================================//
        void KEpsilonStd1stHalf::interaction(size_t index_i, Real dt)
        {
            
            Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            Matd K_prod = Matd::Zero();
            Matd vel_matrix = Matd::Zero();
            Matd strain_tensor = Matd::Zero(), strain_rate_modulus = Matd::Zero();
            K_prod_[index_i] = 0.0, K_adv_[index_i] = 0.0, K_lap_[index_i] = 0.0, strain_rate_[index_i] = 0.0;
            dudx_[index_i] = 0.0, dudy_[index_i] = 0.0, dvdx_[index_i] = 0.0, dvdy_[index_i] = 0.0;
            vel_gradient_mat_[index_i] = Matd::Zero();
            Real mu_t_upperlimit = 1e4 * fluid_.ReferenceViscosity();
            Real mu_t_lowerlimit = 1e-3 * fluid_.ReferenceViscosity();
            Real mu_t = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
            mu_t_[index_i] = std::max(std::min(mu_t_upperlimit, mu_t), mu_t_lowerlimit);

            if (wall_adjacent_cell_flag_[index_i] == 1.0)
            {
                nearwallquantities(index_i);
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);
                     
                    if ((vel_[index_i]).dot(e_ij) > 0.0)
                    {    
                        Vecd distance = pos_[index_i] - pos_[index_j];
                        K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i]) * ((K_[index_j] - K_[index_i]) + (K_grad_[index_j]).dot(distance)) * ((vel_[index_i]).dot(e_ij));

                        if (index_i == 41468)
                        {
                            Vecd veli = vel_[index_i];
                            Real Kadv = K_adv_[index_i];
                            Vecd tkegrad = K_grad_[index_j];
                            Real ki = K_[index_i];
                            Real kj = K_[index_j];
                            Vecd epsgrad = Eps_grad_[index_i];
                            Real epsi = Eps_[index_i];
                            Real epsj = Eps_[index_j];
                            Real c = 0.0;
                        }
                    }
                    else
                    {
                        Vecd distance = pos_[index_j] - pos_[index_i];
                        K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i]) * ((K_[index_i] - K_[index_j]) + (K_grad_[index_i]).dot(distance)) * ((vel_[index_i]).dot(e_ij));

                        if (index_i == 41468)
                        {
                            Vecd veli = vel_[index_i];
                            Real Kadv = K_adv_[index_i];
                            Vecd tkegrad = K_grad_[index_i];
                            Vecd d = distance;
                            Real ki = K_[index_i];
                            Real kj = K_[index_j];
                            Vecd epsgrad = Eps_grad_[index_i];
                            Real epsi = Eps_[index_i];
                            Real epsj = Eps_[index_j];
                            Real c = 0.0;
                        }
                    }

                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij));
                    /*
                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (K_[index_i] - K_[index_j])) * (((vel_[index_i]) - vel_[index_j]).dot(e_ij));      // For better results
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_[index_i] / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij)); //For better results
                    */ 
                }
                strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
                strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
                strain_rate_[index_i] = std::sqrt(strain_rate_modulus.sum());
                
                K_prod_[index_i] = K_prod_p_[index_i];
                Eps_[index_i] = Eps_p_[index_i];
                
                dudx_[index_i] = vel_gradient_mat_[index_i](0, 0);
                dudy_[index_i] = vel_gradient_mat_[index_i](0, 1);
                dvdx_[index_i] = vel_gradient_mat_[index_i](1, 0);
                dvdy_[index_i] = vel_gradient_mat_[index_i](1, 1);

                dK_dt_[index_i] = K_adv_[index_i] + K_prod_[index_i] - rho_[index_i] * Eps_[index_i] + K_lap_[index_i];
            }
            else
            {
                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);
                    
                    if ((vel_[index_i]).dot(e_ij) > 0.0)
                    {
                        Vecd distance = pos_[index_i] - pos_[index_j];
                        K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i]) * ((K_[index_j] - K_[index_i]) + (K_grad_[index_j]).dot(distance)) * ((vel_[index_i]).dot(e_ij));
                    }
                    else
                    {
                        Vecd distance = pos_[index_j] - pos_[index_i];
                        K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i]) * ((K_[index_i] - K_[index_j]) + (K_grad_[index_i]).dot(distance)) * ((vel_[index_i]).dot(e_ij));
                    }

                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij));
                    /*
                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (K_[index_i] - K_[index_j])) * (((vel_[index_i]) - vel_[index_j]).dot(e_ij));      // For better results
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_[index_i] / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij)); // For better results*/
                    
                    vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                    vel_gradient_mat_[index_i] += dW_ij * Vol_[index_j] * vel_matrix;

                }
                strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
                strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
                //strain_rate_modulus = strain_tensor - strain_tensor.trace() / 2.0 * Matd::Identity();
                //Matd S = strain_rate_modulus.array() * strain_rate_modulus.array();
                strain_rate_[index_i] = std::sqrt(strain_rate_modulus.sum());
                
                K_prod = (mu_t_[index_i] * strain_rate_modulus);
                K_prod_[index_i] = K_prod.sum();

                dudx_[index_i] = vel_gradient_mat_[index_i](0, 0);
                dudy_[index_i] = vel_gradient_mat_[index_i](0, 1);
                dvdx_[index_i] = vel_gradient_mat_[index_i](1, 0);
                dvdy_[index_i] = vel_gradient_mat_[index_i](1, 1);

                dK_dt_[index_i] = K_adv_[index_i] + K_prod_[index_i] - rho_[index_i] * Eps_[index_i] + K_lap_[index_i];
            }
        }
        //=================================================================================================//
        void KEpsilonStd1stHalf::update(size_t index_i, Real dt)
        {
            K_[index_i] += (dK_dt_[index_i] / rho_[index_i]) * dt;
            if (K_[index_i] < 0.0)
            {
                K_[index_i] = 1e-7;
            }
        }
        //=================================================================================================//
        KEpsilonStd2ndHalf::KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator) 
            : BaseTurbulence(inner_relation, ghost_creator),
            dEps_dt_(this->particles_->template registerStateVariable<Real>("DissipationChangeRate")),
            wall_adjacent_cell_flag_(this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells")),
            K_grad_(this->particles_->template getVariableDataByName<Vecd>("TKEGradient")),
            Eps_grad_(this->particles_->template getVariableDataByName<Vecd>("DissipationGradient"))
        {}
        //=================================================================================================//
        
        void KEpsilonStd2ndHalf::interaction(size_t index_i, Real dt)
        {
            Real Eps_changerate = 0.0;
            Eps_adv_[index_i] = 0.0, Eps_lap_[index_i] = 0.0, Eps_prod_[index_i] = 0.0, Eps_destruction_[index_i] = 0.0;
            Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            if (wall_adjacent_cell_flag_[index_i] != 1)
            {
                Real mu_t_upperlimit = 1e4 * fluid_.ReferenceViscosity();
                Real mu_t_lowerlimit = 1e-3 * fluid_.ReferenceViscosity();
                Real mu_t = rho_[index_i] * C_mu_ * ((K_[index_i] * K_[index_i]) / (Eps_[index_i]));
                mu_t_[index_i] = std::max(std::min(mu_t_upperlimit, mu_t), mu_t_lowerlimit);

                for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
                {
                    size_t index_j = inner_neighborhood.j_[n];
                    Real dW_ij = inner_neighborhood.dW_ij_[n];
                    Real r_ij = inner_neighborhood.r_ij_[n];
                    Vecd &e_ij = inner_neighborhood.e_ij_[n];
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);

                    if ((vel_[index_i]).dot(e_ij) > 0.0)
                    {
                        Vecd distance = pos_[index_i] - pos_[index_j];
                        Eps_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i]) * ((Eps_[index_j] - Eps_[index_i]) + (Eps_grad_[index_j]).dot(distance)) * ((vel_[index_i]).dot(e_ij));
                    }
                    else
                    {
                        Vecd distance = pos_[index_j] - pos_[index_i];
                        Eps_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i]) * ((Eps_[index_i] - Eps_[index_j]) + (Eps_grad_[index_i]).dot(distance)) * ((vel_[index_i]).dot(e_ij));
                    }

                    Eps_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * (fluid_.ReferenceViscosity() + mu_t_avg / sigma_eps_) * ((Eps_[index_i] - Eps_[index_j]) / (r_ij));
                    /*
                    Eps_adv_[index_i] += -(dW_ij * Vol_[index_j] * rho_[index_i] * (Eps_[index_i] - Eps_[index_j])) * (((vel_[index_i]) - vel_[index_j]).dot(e_ij));
                    Eps_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * (fluid_.ReferenceViscosity() + mu_t_[index_i] / sigma_eps_) * ((Eps_[index_i] - Eps_[index_j]) / (r_ij)); //For better results*/ 

                    Eps_changerate = Eps_adv_[index_i] + Eps_lap_[index_i];
                }
                Eps_prod_[index_i] = C1_eps_ * (Eps_[index_i] / (K_[index_i])) * K_prod_[index_i];
                Eps_destruction_[index_i] = -C2_eps_ * rho_[index_i] * (Eps_[index_i] * Eps_[index_i]) / (K_[index_i]);
                dEps_dt_[index_i] = Eps_changerate + Eps_prod_[index_i] + Eps_destruction_[index_i];
            }
        }
        //=================================================================================================//
        void KEpsilonStd2ndHalf::update(size_t index_i, Real dt)
        {
            if (wall_adjacent_cell_flag_[index_i] != 1)
            {
                Eps_[index_i] += (dEps_dt_[index_i] / rho_[index_i]) * dt;
            }
            if (Eps_[index_i] < 0.0)
            {
                Eps_[index_i] = 1e-7;
            }
        }
        //=================================================================================================//
         StdWallFunctionFVM::StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator)
            : BaseTurbulence(inner_relation, ghost_creator), y_star_(this->particles_->template registerStateVariable<Real>("Ystar")),
              yp_(this->particles_->template getVariableDataByName<Real>("WallNormalDistance")),
              wall_normal_(this->particles_->template getVariableDataByName<Vecd>("WallNormal")),
              vel_gradient_mat_(this->particles_->template registerStateVariable<Matd>("VelocityGradient")), 
              von_kar_(0.4187), E_(9.793)
              {}
        //=================================================================================================//
        void StdWallFunctionFVM::nearwallquantities(size_t index_i)
        {
            y_star_[index_i] = (rho_[index_i] * std::pow(C_mu_, 0.25) * std::pow(K_[index_i], 0.5) * yp_[index_i]) / (fluid_.ReferenceViscosity());
            Real u_star;
            Vecd veltangential = (vel_[index_i] - wall_normal_[index_i].dot(vel_[index_i]) * (wall_normal_[index_i]));

            if (y_star_[index_i] >= 11.225)
            {
                u_star = (1.0 / von_kar_) * std::log(E_ * y_star_[index_i]);
                mu_t_[index_i] = fluid_.ReferenceViscosity() * ((y_star_[index_i]) / (1 / von_kar_ * std::log(E_ * y_star_[index_i])) - 1.0);
                
                Tau_wall_[index_i] = (veltangential.norm() * std::pow(C_mu_, 0.25) * std::pow(K_[index_i], 0.5) * rho_[index_i]) / (u_star);
                vel_gradient_mat_[index_i] = Matd::Zero();
                vel_gradient_mat_[index_i](0, 1) = Tau_wall_[index_i] / (rho_[index_i] * pow(C_mu_, 0.25) * pow(K_[index_i], 0.5) * von_kar_ * yp_[index_i]);
                //vel_gradient_mat_[index_i](0, 1) = veltangential.norm() / (yp_[index_i] * log(E_ * ystar_[index_i]));

                K_prod_p_[index_i] = std::pow(Tau_wall_[index_i], 2.0) / (von_kar_ * rho_[index_i] * std::pow(C_mu_, 0.25) * std::pow(K_[index_i], 0.5) * yp_[index_i]);
                Eps_p_[index_i] = (std::pow(C_mu_, 3.0 / 4.0) * std::pow(K_[index_i], 1.5)) / (von_kar_ * yp_[index_i]);
               
            }
            else if (y_star_[index_i] < 11.225)
            {
                u_star = y_star_[index_i];
                Tau_wall_[index_i] = fluid_.ReferenceViscosity() * veltangential.norm() / yp_[index_i];
                vel_gradient_mat_[index_i] = Matd::Zero();
                vel_gradient_mat_[index_i](0, 1) = Tau_wall_[index_i] / fluid_.ReferenceViscosity();
                K_prod_p_[index_i] = 0.0;
                Eps_p_[index_i] = (K_[index_i] * 2.0 * fluid_.ReferenceViscosity()) / (rho_[index_i] * yp_[index_i] * yp_[index_i]);
            }  
        }
        //=================================================================================================// 
    }// namespace fluid_dynamics

}// namespace SPH
#endif // TURBULENCEMODEL_CPP