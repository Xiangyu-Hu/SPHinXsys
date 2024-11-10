#ifndef TURBULENCEMODEL_HPP
#define TURBULENCEMODEL_HPP
#include "turbulence_model.h"


namespace SPH
{  
    namespace fluid_dynamics
    { 
        //=================================================================================================//
        template <class RiemannSolverType>
        KEpsilonStd1stHalf<RiemannSolverType>::KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator, Real limiter_parameter) 
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
              riemann_solver_(this->fluid_, this->fluid_),
              K_grad_(this->particles_->template getVariableDataByName<Vecd>("TKEGradient")), 
              Eps_grad_(this->particles_->template getVariableDataByName<Vecd>("DissipationGradient"))
              {}
        //=================================================================================================//
        template <class RiemannSolverType>
        void KEpsilonStd1stHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
        {
            ExendedFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], K_[index_i], Eps_[index_i]);
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
                    ExendedFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], K_[index_j], Eps_[index_j]);
                    ExtendedFluidStarState interface_state = riemann_solver_.getExtendedInterfaceState(state_i, state_j, e_ij);
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);
                    Real K_avg = 0.5 * (K_[index_i] + K_[index_j]);
                    
                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * interface_state.rho_) * (interface_state.K_) * (interface_state.vel_.dot(e_ij));
                    Real x = 1.0;
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij));
                    
                    
                    /*
                    K_adv_[index_i] += -2.0 * (dW_ij * Vol_[index_j] * interface_state.rho_ * K_avg * interface_state.vel_).dot(e_ij); // Riemann solver
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij)); */ 
                    
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
                    ExendedFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], K_[index_j], Eps_[index_j]);
                    ExtendedFluidStarState interface_state = riemann_solver_.getExtendedInterfaceState(state_i, state_j, e_ij);
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);
                    Real K_avg = 0.5 * (K_[index_i] + K_[index_j]);
                    
                    K_adv_[index_i] += -(dW_ij * Vol_[index_j] * interface_state.rho_) * (interface_state.K_) * (interface_state.vel_.dot(e_ij));
                    Real C = 1.0;
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij));
                    
                    
                    /*
                    K_adv_[index_i] += -2.0 * (dW_ij * Vol_[index_j] * interface_state.rho_ * K_avg * interface_state.vel_).dot(e_ij);                        // For better results
                    K_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * ((fluid_.ReferenceViscosity() + mu_t_avg / sigma_k_) * (K_[index_i] - K_[index_j]) / (r_ij)); // For better results*/ 
                    
                    vel_matrix = (vel_[index_i] - vel_[index_j]) * e_ij.transpose();
                    vel_gradient_mat_[index_i] += dW_ij * Vol_[index_j] * vel_matrix;

                }
                strain_tensor = 0.5 * (vel_gradient_mat_[index_i] + vel_gradient_mat_[index_i].transpose());
                strain_rate_modulus = 2.0 * strain_tensor.array() * strain_tensor.array();
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
        template <class RiemannSolverType>
        void KEpsilonStd1stHalf<RiemannSolverType> ::update(size_t index_i, Real dt)
        {
            K_[index_i] += (dK_dt_[index_i] / rho_[index_i]) * dt;
            if (K_[index_i] < 0.0)
            {
                K_[index_i] = 1e-6;
            }
        }
        //=================================================================================================//
        template <class RiemannSolverType>
        KEpsilonStd2ndHalf<RiemannSolverType>::KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator, Real limiter_parameter) 
            : BaseTurbulence(inner_relation, ghost_creator),
            dEps_dt_(this->particles_->template registerStateVariable<Real>("DissipationChangeRate")),
            wall_adjacent_cell_flag_(this->particles_->template getVariableDataByName<Real>("FlagForWallAdjacentCells")),
            riemann_solver_(this->fluid_, this->fluid_),
            K_grad_(this->particles_->template getVariableDataByName<Vecd>("TKEGradient")),
            Eps_grad_(this->particles_->template getVariableDataByName<Vecd>("DissipationGradient"))
        {}
        //=================================================================================================//
        template <class RiemannSolverType>
        void KEpsilonStd2ndHalf<RiemannSolverType>::interaction(size_t index_i, Real dt)
        {
            ExendedFluidState state_i(rho_[index_i], vel_[index_i], p_[index_i], K_[index_i], Eps_[index_i]);
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
                    ExendedFluidState state_j(rho_[index_j], vel_[index_j], p_[index_j], K_[index_i], Eps_[index_i]);
                    ExtendedFluidStarState interface_state = riemann_solver_.getExtendedInterfaceState(state_i, state_j, e_ij);
                    Real mu_t_avg = (2.0 * mu_t_[index_i] * mu_t_[index_j]) / (mu_t_[index_i] + mu_t_[index_j]);
                    Real Eps_avg = 0.5 * (Eps_[index_i] + Eps_[index_j]);
                    
                    Eps_adv_[index_i] += -(dW_ij * Vol_[index_j] * interface_state.rho_) * (interface_state.Eps_) * (interface_state.vel_.dot(e_ij));
                    Real V = 1.2; 
                    Eps_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * (fluid_.ReferenceViscosity() + mu_t_avg / sigma_eps_) * ((Eps_[index_i] - Eps_[index_j]) / (r_ij));
                    
                    /* Eps_adv_[index_i] += -2.0 * (dW_ij * Vol_[index_j] * interface_state.rho_ * Eps_avg * interface_state.vel_).dot(e_ij);
                    Eps_lap_[index_i] += 2.0 * dW_ij * Vol_[index_j] * (fluid_.ReferenceViscosity() + mu_t_avg / sigma_eps_) * ((Eps_[index_i] - Eps_[index_j]) / (r_ij)); // For better results
                    */
                   

                    Eps_changerate = Eps_adv_[index_i] + Eps_lap_[index_i];
                }
                Eps_prod_[index_i] = C1_eps_ * (Eps_[index_i] / (K_[index_i])) * K_prod_[index_i];
                Eps_destruction_[index_i] = -C2_eps_ * rho_[index_i] * (Eps_[index_i] * Eps_[index_i]) / (K_[index_i]);
                dEps_dt_[index_i] = Eps_changerate + Eps_prod_[index_i] + Eps_destruction_[index_i];
            }
        }
        //=================================================================================================//
        template <class RiemannSolverType>
        void KEpsilonStd2ndHalf<RiemannSolverType>::update(size_t index_i, Real dt)
        {
            if (wall_adjacent_cell_flag_[index_i] != 1)
            {
                Eps_[index_i] += (dEps_dt_[index_i] / rho_[index_i]) * dt;
            }
            if (Eps_[index_i] < 0.0)
            {
                Eps_[index_i] = 1e-6;
            }
        }



    }// namespace fluid_dynamics

}// namespace SPH
#endif // TURBULENCEMODEL_HPP