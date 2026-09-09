
#include "udf_P_refinement.hpp"
namespace SPH
{

namespace fluid_dynamics
{

namespace udf
{

    P_refinement_GetVelocityGradient<Inner<>>::
        P_refinement_GetVelocityGradient(BaseInnerRelation& inner_relation)
        : P_refinement_GetVelocityGradient<DataDelegateInner>(inner_relation),
        turbu_B_(particles_->getVariableDataByName<Matd>("TurbulentLinearGradientCorrectionMatrix")),
        B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")) {}

    void P_refinement_GetVelocityGradient<Inner<>>::interaction(size_t index_i, Real dt)
    {
        velocity_gradient_only_P_[index_i] = Matd::Zero();
        k_gradient_only_P_[index_i] = Vecd::Zero();
        omega_gradient_only_P_[index_i] = Vecd::Zero();
        if (is_near_wall_P1_[index_i] == 1)
        {
            Vecd vel_i = vel_[index_i];
            Real k_i = turbu_k_[index_i];
            Real w_i = turbu_omega_[index_i];
            const Neighborhood& inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                Vecd nablaW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
                velocity_gradient_only_P_[index_i] += -(vel_i - vel_[index_j]) * nablaW_ijV_j.transpose();
                k_gradient_only_P_[index_i] += -1.0 * (k_i - turbu_k_[index_j]) * nablaW_ijV_j;
                omega_gradient_only_P_[index_i] += -1.0 * (w_i - turbu_omega_[index_j]) * nablaW_ijV_j;
            }
        }
    }

    void P_refinement_GetVelocityGradient<Inner<>>::update(size_t index_i, Real dt)
    {
        if (is_near_wall_P1_[index_i] == 1)
        {
            velocity_gradient_only_P_[index_i] *= turbu_B_[index_i];
        }
    }

    P_refinement_GetVelocityGradient<Contact<Wall>>::P_refinement_GetVelocityGradient(BaseContactRelation& contact_relation)
        : InteractionWithWall<P_refinement_GetVelocityGradient>(contact_relation) {}

    void P_refinement_GetVelocityGradient<Contact<Wall>>::interaction(size_t index_i, Real dt)
    {
        if (is_near_wall_P1_[index_i] == 1)
        {
            Matd vel_grad = Matd::Zero();
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Vecd* vel_ave_k = wall_vel_ave_[k];
                Real* Vol_k = wall_Vol_[k];
                Neighborhood& contact_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    Vecd nablaW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
                    vel_grad += -1.0 * 2.0 * (vel_[index_i] - vel_ave_k[index_j]) * nablaW_ijV_j.transpose();
                }
            }
            velocity_gradient_only_P_[index_i] += vel_grad;
        }
    }

    template <int Ny, int TypeTDMA>
    P_refinement<Ny, TypeTDMA>::
        P_refinement(SPHBody& sph_body, Real constant_y_p, Real constant_y_p_node)
        : LocalDynamics(sph_body),
        num_sub_node_(ny),
        friction_velocity_from_sublayer_(particles_->registerStateVariableData<Real>("FrictionVelocityFromSublayer")),
        target_flow_rate_in_sublayer_(particles_->registerStateVariableData<Real>("TargetFlowRateInSublayer")),
        vel_ps_magnitude_(particles_->registerStateVariableData<Real>("VelPS")),
        dudn_for_local_flow_rate_(particles_->registerStateVariableData<Real>("dudnForLocalFlowRate")),
        utau_node_(particles_->registerStateVariableData<Real>("utauNode")),
        node_value_vel_(particles_->registerStateVariableData<Vec6d>("NodeValue")),
        node_value_k_(particles_->registerStateVariableData<Vec6d>("NodeValueTKE")),
        dUdn_P_sublayer_magnitude_(particles_->registerStateVariableData<Real>("dUdnFromSublayerMagnitude")),
        dUdn_P_sublayer_(particles_->registerStateVariableData<Matd>("dUdnFromSublayer")),
        vel_nodeO_(particles_->registerStateVariableData<Real>("VelNodeO")),
        vel_nodeUM_(particles_->registerStateVariableData<Real>("VelNodeUM")),
        global_flow_rate_over_P_(particles_->registerStateVariableData<Real>("global_flow_rate_over_P_")),
        half_flow_rate_over_P_(particles_->registerStateVariableData<Real>("half_flow_rate_over_P_")),
        is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")),
        vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
        turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
        turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
        rho_(particles_->getVariableDataByName<Real>("Density")),
        viscosity_(sph_body_->getMaterialProperty<Viscosity>()),
        mu_(viscosity_.ReferenceViscosity()),
        velocity_gradient_only_P_(particles_->getVariableDataByName<Matd>("VelocityGradientInnerOnlyP")),
        turbu_mu_(particles_->getVariableDataByName<Real>("TurbulentViscosity")),
        distance_to_dummy_interface_(particles_->getVariableDataByName<Real>("DistanceToDummyInterface")),
        wall_shear_stress_(particles_->getVariableDataByName<Real>("WallShearStress")),
        e_nearest_normal_(particles_->getVariableDataByName<Vecd>("WallNearestNormalUnitVector")),
        fluid_particle_spacing_(sph_body.getSPHAdaptation().ReferenceSpacing()),
        physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")),
        velocity_gradient_(particles_->getVariableDataByName<Matd>("TurbulentVelocityGradient")),
        k_gradient_only_P_(particles_->getVariableDataByName<Vecd>("TurbulentKineticEnergyGradientOnlyP")),
        omega_gradient_only_P_(particles_->getVariableDataByName<Vecd>("TurbulentSpecificDissipationGradientOnlyP"))
    {
        check_num_node_consistency();
        construct_node_distribution(constant_y_p, constant_y_p_node);
        output_node_distribution();
    }

    template <int Ny, int TypeTDMA>
    void P_refinement<Ny, TypeTDMA>::update(size_t index_i, Real dt)
    {
        Real vel_nodeO_i_prior = 0.0;
        if (is_near_wall_P1_[index_i] == 1)
        {
            vel_nodeO_i_prior = vel_nodeO_[index_i];
        }

        friction_velocity_from_sublayer_[index_i] = 0.0;
        target_flow_rate_in_sublayer_[index_i] = 0.0;
        global_flow_rate_over_P_[index_i] = 0.0;
        half_flow_rate_over_P_[index_i] = 0.0;
        vel_ps_magnitude_[index_i] = 0.0;
        dudn_for_local_flow_rate_[index_i] = 0.0;
        utau_node_[index_i] = 0.0;
        node_value_vel_[index_i] = Vec6d::Zero();
        node_value_k_[index_i] = Vec6d::Zero();
        dUdn_P_sublayer_magnitude_[index_i] = 0.0;
        dUdn_P_sublayer_[index_i] = Matd::Zero();
        vel_nodeO_[index_i] = 0.0;
        vel_nodeUM_[index_i] = 0.0;
        double U_nodeO = 0.0;
        double U_nodeUM = 0.0;
        if (is_near_wall_P1_[index_i] == 1)
        {
            Real dudn_outer = 0.0;
            Real dkdn_outer = 0.0;
            Real dwdn_outer = 0.0;
            Real flow_rate_local = 0.0;
            Real u_outer = 0.0;
            Real k_outer = 0.0;
            Real omega_outer = 0.0;
            Real nut_outer = 0.0;
            Real friction_vel_magnitude_outer = 0.0;
            Vecd normal = e_nearest_normal_[index_i];
            Real nu = mu_ / rho_[index_i];
            u_outer = obtainTangentialComponent(vel_[index_i], normal);
            k_outer = turbu_k_[index_i];
            omega_outer = turbu_omega_[index_i];
            nut_outer = turbu_mu_[index_i] / rho_[index_i];
            friction_vel_magnitude_outer = std::sqrt(wall_shear_stress_[index_i] / rho_[index_i]);
            Real distance_to_wall = sublayer_height_contant_;
            Vecd velocity_gradient_only_P_normal = velocity_gradient_only_P_[index_i] * normal;
            Real dudn_from_SPH = obtainTangentialComponent(velocity_gradient_only_P_normal, normal);
            dudn_outer = dudn_from_SPH;
            dkdn_outer = k_gradient_only_P_[index_i].dot(normal);
            dwdn_outer = omega_gradient_only_P_[index_i].dot(normal);
            Real averaged_vel_over_P = u_outer;
            SublayerResult sublayer_result{};
            vel_nodeO_i_prior = u_outer;
            flow_rate_local = get_loacal_flow_rate(averaged_vel_over_P * fluid_particle_spacing_, dudn_outer, vel_nodeO_i_prior, fluid_particle_spacing_);
            U_nodeO = 0.0;
            U_nodeUM = 0.0;
            sublayer_result = solve_1D_sublayer_Dirichlet(nu, u_outer, k_outer, omega_outer, std::abs(dudn_outer),
                nut_outer, distance_to_wall, friction_vel_magnitude_outer, std::abs(flow_rate_local), dkdn_outer, dwdn_outer, U_nodeO, U_nodeUM);
            friction_velocity_from_sublayer_[index_i] = sublayer_result.sublayer_utau;
            vel_nodeO_[index_i] = U_nodeO;
            vel_nodeUM_[index_i] = U_nodeUM;
            if (is_truncated_output_)
            {

                node_value_vel_[index_i][0] = sublayer_result.sublayer_vel[0];
                node_value_k_[index_i][0] = sublayer_result.sublayer_k[0];

                node_value_vel_[index_i][1] = sublayer_result.sublayer_vel[1];
                node_value_k_[index_i][1] = sublayer_result.sublayer_k[1];

                for (int i = 2; i < node_dim_output_limit_; ++i) {
                    node_value_vel_[index_i][i] = sublayer_result.sublayer_vel[i + num_skip_output_];
                    node_value_k_[index_i][i] = sublayer_result.sublayer_k[i + num_skip_output_];
                }
            }
            else
            {
                for (int i = 0; i < num_sub_node_; ++i) {
                    node_value_vel_[index_i][i] = sublayer_result.sublayer_vel[i];
                    node_value_k_[index_i][i] = sublayer_result.sublayer_k[i];
                }
            }

            target_flow_rate_in_sublayer_[index_i] = flow_rate_local;
            global_flow_rate_over_P_[index_i] = averaged_vel_over_P * fluid_particle_spacing_;
            half_flow_rate_over_P_[index_i] = (U_nodeO + (U_nodeO + dudn_outer * 0.5 * fluid_particle_spacing_)) * (0.5 * fluid_particle_spacing_) / 2.0;
            vel_ps_magnitude_[index_i] = U_nodeO + dudn_outer * 0.5 * fluid_particle_spacing_;
            dudn_for_local_flow_rate_[index_i] = dudn_outer;
            utau_node_[index_i] = friction_vel_magnitude_outer;

            Vecd vel_tangential = vel_[index_i] - vel_[index_i].dot(normal) * normal;
            Real tangential_velocity_P_magnitude = vel_tangential.norm();
            Vecd tangential = vel_tangential / (tangential_velocity_P_magnitude + TinyReal);
            Real dUdn_P_sublayer_magnitude = 0.0;
            int index_nodeU = num_sub_node_ - 1;
            Real tangential_velocity_node_U = sublayer_result.sublayer_vel[index_nodeU];
            Real dist_nodeU_P = sublayer_height_contant_ - sublayer_y_[index_nodeU];
            dUdn_P_sublayer_magnitude = std::abs(tangential_velocity_P_magnitude - tangential_velocity_node_U) / (dist_nodeU_P + TinyReal);
            Matd dUdn_P_sublayer = dUdn_P_sublayer_magnitude * (tangential * normal.transpose());
            dUdn_P_sublayer_magnitude_[index_i] = dUdn_P_sublayer_magnitude;
            dUdn_P_sublayer_[index_i] = dUdn_P_sublayer;
        }
    }

    template <int Ny, int TypeTDMA>
    typename P_refinement<Ny, TypeTDMA>::SublayerResult P_refinement<Ny, TypeTDMA>::solve_1D_sublayer_Dirichlet(double kinematic_viscosity, double u_p_outer, double k_p_outer,
        double w_p_outer, double vel_grad_p_outer, double nut_p_outer, double h_sublayer, double utau_outer,
        double Q_target, double k_grad_p_outer, double w_grad_p_outer, double& vel_nodeO, double& vel_nodeUM)
    {

        double utau_init = utau_outer;
        double nu = kinematic_viscosity;

        double u_init = u_p_outer;
        double k_init = k_p_outer;
        double turbu_omega_init = w_p_outer;

        double convergence_criteria_outer = 1.0e-3;
        double tiny = 1.0e-6;
        double relax_u = 0.9;
        double relax_k = 0.9;
        double relax_w = 0.6;
        double relax_utau = 0.4;
        double yplus_min = 0.01;

        double flow_rate_target = Q_target;
        double utau = utau_init;

        int type_tdma = type_tdma_;

        double y[ny];
        double dist_node_i_to_wall[ny]{};
        for (int i = 0; i < ny; ++i)
        {
            y[i] = sublayer_y_[i];
            dist_node_i_to_wall[i] = y[i];
        }

        double hy = sublayer_node_uniform_distance_;
        double y_p = sublayer_y_p_constant_;
        double utau_min = yplus_min * nu / y_p;
        double yplus = utau * y_p / nu;
        double u_p = utau * yplus;
        double turbu_omega_p = 6.0 * nu / (std_kw_beta_i_ * y_p * y_p);
        double u_init_value[ny];
        double k_init_value[ny];
        double turbu_omega_init_value[ny];
        std::fill_n(u_init_value, ny, u_init);
        std::fill_n(k_init_value, ny, k_init);
        std::fill_n(turbu_omega_init_value, ny, turbu_omega_init);
        double u_nodeUM = u_p_outer;
        double k_nodeUM = k_p_outer;
        double w_nodeUM = w_p_outer;
        double nut_nodeUM = nut_p_outer;
        double u_nodeO = u_nodeUM;
        double phi_current[3 * ny];
        double phi_solved[3 * ny];
        for (int j = 0; j < ny; ++j)
        {
            phi_current[j] = u_init_value[j];
            phi_current[j + ny] = k_init_value[j];
            phi_current[j + 2 * ny] = turbu_omega_init_value[j];
        }
        std::copy_n(phi_current, 3 * ny, phi_solved);
        double differ = 1.0;
        int num_iter_out = 0;
        int n_start = 0;
        int last = 0;
        double flow_rate_current = 0.0;
        while (differ > convergence_criteria_outer)
        {
            double* u_star = phi_solved;
            double* k_star = phi_solved + ny;
            double* turbu_omega_star = phi_solved + 2 * ny;
            if (num_iter_out != 0)
            {
                yplus = utau * y_p / nu;
                u_p = utau * yplus;
                turbu_omega_p = 6.0 * nu / (std_kw_beta_i_ * y_p * y_p);
                u_nodeUM = u_p_outer;
                k_nodeUM = k_p_outer;
                w_nodeUM = w_p_outer;
                nut_nodeUM = nut_p_outer;
                u_nodeO = u_nodeUM;
            }
            double dudy_discretized_central[ny]{};
            double dkdy[ny]{};
            double dwdy[ny]{};
            for (int i = 1; i < ny - 1; ++i)
            {
                dudy_discretized_central[i] = (u_star[i + 1] - u_star[i - 1]) / (2.0 * hy);
                dkdy[i] = (k_star[i + 1] - k_star[i - 1]) / (2.0 * hy);
                dwdy[i] = (turbu_omega_star[i + 1] - turbu_omega_star[i - 1]) / (2.0 * hy);
            }
            dudy_discretized_central[0] = (u_star[1] + u_star[0]) / (2.0 * hy);
            dkdy[0] = (k_star[1] - k_star[0]) / (2.0 * hy);
            dwdy[0] = (turbu_omega_star[1] - turbu_omega_star[0]) / (2.0 * hy);
            dudy_discretized_central[ny - 1] = (u_nodeUM - u_star[ny - 2]) / (2.0 * hy);
            dkdy[ny - 1] = (k_nodeUM - k_star[ny - 2]) / (2.0 * hy);
            dwdy[ny - 1] = (w_nodeUM - turbu_omega_star[ny - 2]) / (2.0 * hy);
            double diffusion_coefficient_k[ny]{};
            double diffusion_coefficient_turbu_omega[ny]{};
            for (int i = 0; i < ny; ++i) {
                diffusion_coefficient_k[i] = nu + std_kw_sigma_star_ * k_star[i] / (turbu_omega_star[i] + tiny);
                diffusion_coefficient_turbu_omega[i] = nu + std_kw_sigma_ * k_star[i] / (turbu_omega_star[i] + tiny);
            }
            double turbu_omega_tilde[ny]{};
            double nut_star[ny]{};
            for (int i = 0; i < ny; ++i)
            {
                turbu_omega_tilde[i] = std::max(
                    turbu_omega_star[i],
                    std_kw_C_lim_ * dudy_discretized_central[i] / std_kw_beta_star_5_
                );
                nut_star[i] = k_star[i] / (turbu_omega_tilde[i] + tiny);
            }
            double dudy_star[ny]{};
            for (int i = 0; i < ny; ++i)
            {
                dudy_star[i] = dudy_discretized_central[i];
            }
            double grad_prod[ny]{};
            double part_cross_diffusion[ny]{};
            for (int i = 0; i < ny; ++i) {
                grad_prod[i] = dkdy[i] * dwdy[i];
                double sigma_d = (grad_prod[i] > 0.0) ? std_kw_sigma_do_ : 0.0;
                part_cross_diffusion[i] = sigma_d / (turbu_omega_star[i] + tiny) * grad_prod[i];
            }
            double a_u[ny]{};
            double b_u[ny]{};
            double c_u[ny]{};
            double d_u[ny]{};
            int type_u_discretization = 2;
            double pressure_gradient_constant_part[ny]{};
            double additional_coefficient_iminus1[ny]{};
            double additional_coefficient_iplus1[ny]{};
            for (int i = 0; i < ny; ++i) {
                pressure_gradient_constant_part[i] = (utau * utau) / dist_node_i_to_wall[i];
                additional_coefficient_iminus1[i] = -hy * (nu + nut_star[i]) / dist_node_i_to_wall[i] / (2.0);
                additional_coefficient_iplus1[i] = hy * (nu + nut_star[i]) / dist_node_i_to_wall[i] / (2.0);
            }
            for (int i = 1; i < ny - 1; ++i)
            {
                double nu_eff_i = nu + nut_star[i];
                double nu_eff_i_plus = nu + nut_star[i + 1];
                double nu_eff_i_minus = nu + nut_star[i - 1];
                double nu_eff_i_plus_half = 2.0 * nu_eff_i_plus * nu_eff_i / std::max((nu_eff_i_plus + nu_eff_i), tiny);
                double nu_eff_i_minus_half = 2.0 * nu_eff_i_minus * nu_eff_i / std::max((nu_eff_i_minus + nu_eff_i), tiny);
                if (type_u_discretization == 2)
                {
                    a_u[i] = -nu_eff_i_minus_half + additional_coefficient_iminus1[i];
                    b_u[i] = std::max((nu_eff_i_plus_half + nu_eff_i_minus_half), tiny);
                    c_u[i] = -nu_eff_i_plus_half + additional_coefficient_iplus1[i];
                    d_u[i] = pressure_gradient_constant_part[i] * hy * hy;
                }
                else
                {
                    std::cout << "type_u_discretization: Type not define! Stop here." << std::endl;
                    std::cin.get();
                }
            }
            a_u[0] = 0.0;
            b_u[0] = 1.0;
            c_u[0] = 0.0;
            d_u[0] = u_p;
            last = ny - 1;
            double nu_eff_last = nu + nut_star[last];
            double nu_eff_last_plus = nu + nut_nodeUM;
            double nu_eff_last_minus = nu + nut_star[last - 1];
            double nu_eff_last_plus_half = 2.0 * nu_eff_last_plus * nu_eff_last / std::max((nu_eff_last_plus + nu_eff_last), tiny);
            double nu_eff_last_minus_half = 2.0 * nu_eff_last_minus * nu_eff_last / std::max((nu_eff_last_minus + nu_eff_last), tiny);
            if (type_u_discretization == 2)
            {
                a_u[last] = -nu_eff_last_minus_half + additional_coefficient_iminus1[last];
                b_u[last] = std::max((nu_eff_last_plus_half + nu_eff_last_minus_half), tiny);
                c_u[last] = 0.0;
                d_u[last] = pressure_gradient_constant_part[last] * hy * hy + (nu_eff_last_plus_half - additional_coefficient_iplus1[last]) * u_nodeUM;
            }
            else
            {
                std::cout << "type_u_discretization: Type not define! Stop here." << std::endl;
                std::cin.get();
            }
            double U_new[ny]{};

            if (type_tdma == 0)
            {
                tdma(ny, a_u, b_u, c_u, d_u, U_new);
            }
            else if (type_tdma == 5)
            {
                tdma5(a_u, b_u, c_u, d_u, U_new);
            }
            else if (type_tdma == 10)
            {
                tdma10(a_u, b_u, c_u, d_u, U_new);
            }
            else
            {
                std::cout << "TDMA: Type not define! Stop here." << std::endl;
                std::cin.get();
            }
            double a_k[ny]{};
            double b_k[ny]{};
            double c_k[ny]{};
            double d_k[ny]{};
            for (int i = 1; i < ny - 1; ++i) {
                double Dk_i = diffusion_coefficient_k[i];
                double Dk_i_plus = diffusion_coefficient_k[i + 1];
                double Dk_i_minus = diffusion_coefficient_k[i - 1];
                double Dk_i_plus_half = 2.0 * Dk_i_plus * Dk_i / std::max((Dk_i_plus + Dk_i), tiny);
                double Dk_i_minus_half = 2.0 * Dk_i_minus * Dk_i / std::max((Dk_i_minus + Dk_i), tiny);
                a_k[i] = -1.0 * Dk_i_minus_half;
                b_k[i] = Dk_i_plus_half + Dk_i_minus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[i];
                c_k[i] = -1.0 * Dk_i_plus_half;
                d_k[i] = hy * hy * nut_star[i] * dudy_star[i] * dudy_star[i];
            }
            double Dk_first = diffusion_coefficient_k[0];
            double Dk_first_plus = diffusion_coefficient_k[0 + 1];
            double Dk_first_plus_half = 2.0 * Dk_first_plus * Dk_first / std::max((Dk_first_plus + Dk_first), tiny);
            a_k[0] = 0.0;
            b_k[0] = Dk_first_plus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[0];
            c_k[0] = -1.0 * Dk_first_plus_half;
            d_k[0] = hy * hy * nut_star[0] * dudy_star[0] * dudy_star[0];
            last = ny - 1;
            double Dk_last = diffusion_coefficient_k[last];
            double Dk_last_plus = nu + std_kw_sigma_star_ * k_nodeUM / (w_nodeUM + tiny);
            double Dk_last_minus = diffusion_coefficient_k[last - 1];
            double Dk_last_plus_half = 2.0 * Dk_last_plus * Dk_last / std::max((Dk_last_plus + Dk_last), tiny);
            double Dk_last_minus_half = 2.0 * Dk_last_minus * Dk_last / std::max((Dk_last_minus + Dk_last), tiny);
            a_k[last] = -1.0 * Dk_last_minus_half;
            b_k[last] = Dk_last_plus_half + Dk_last_minus_half + hy * hy * std_kw_beta_star_ * turbu_omega_star[last];
            c_k[last] = 0.0;
            d_k[last] = hy * hy * nut_star[last] * dudy_star[last] * dudy_star[last] + Dk_last_plus_half * k_nodeUM;
            double K_new[ny]{};
            if (type_tdma == 0)
            {
                tdma(ny, a_k, b_k, c_k, d_k, K_new);
            }
            else if (type_tdma == 5)
            {
                tdma5(a_k, b_k, c_k, d_k, K_new);
            }
            else if (type_tdma == 10)
            {
                tdma10(a_k, b_k, c_k, d_k, K_new);
            }
            else
            {
                std::cout << "TDMA: Type not define! Stop here." << std::endl;
                std::cin.get();
            }
            double k_min = 1e-10;
            for (int i = 0; i < ny; ++i) {
                K_new[i] = std::max(K_new[i], k_min);
            }
            double a_w[ny]{};
            double b_w[ny]{};
            double c_w[ny]{};
            double d_w[ny]{};
            for (int i = 1; i < ny - 1; ++i) {
                double Dw_i = diffusion_coefficient_turbu_omega[i];
                double Dw_i_plus = diffusion_coefficient_turbu_omega[i + 1];
                double Dw_i_minus = diffusion_coefficient_turbu_omega[i - 1];
                double Dw_i_plus_half = 2.0 * Dw_i_plus * Dw_i / std::max((Dw_i_plus + Dw_i), tiny);
                double Dw_i_minus_half = 2.0 * Dw_i_minus * Dw_i / std::max((Dw_i_minus + Dw_i), tiny);
                a_w[i] = -1.0 * Dw_i_minus_half;
                b_w[i] = Dw_i_plus_half + Dw_i_minus_half + hy * hy * std_kw_beta_ * turbu_omega_star[i];
                c_w[i] = -1.0 * Dw_i_plus_half;
                double part_production = std_kw_alpha_ * turbu_omega_star[i] / (k_star[i] + tiny) * nut_star[i] * dudy_star[i] * dudy_star[i];
                d_w[i] = hy * hy * (part_production + part_cross_diffusion[i]);
            }
            a_w[0] = 0.0;
            b_w[0] = 1.0;
            c_w[0] = 0.0;
            d_w[0] = turbu_omega_p;
            last = ny - 1;
            double Dw_last = diffusion_coefficient_turbu_omega[last];
            double Dw_last_plus = nu + std_kw_sigma_ * k_nodeUM / (w_nodeUM + tiny);
            double Dw_last_minus = diffusion_coefficient_turbu_omega[last - 1];
            double Dw_last_plus_half = 2.0 * Dw_last_plus * Dw_last / std::max((Dw_last_plus + Dw_last), tiny);
            double Dw_last_minus_half = 2.0 * Dw_last_minus * Dw_last / std::max((Dw_last_minus + Dw_last), tiny);
            a_w[last] = -1.0 * Dw_last_minus_half;
            b_w[last] = Dw_last_plus_half + Dw_last_minus_half + hy * hy * std_kw_beta_ * turbu_omega_star[last];
            c_w[last] = 0.0;
            double part_production_last = std_kw_alpha_ * turbu_omega_star[last] / (k_star[last] + tiny) * nut_star[last] * dudy_star[last] * dudy_star[last];
            d_w[last] = hy * hy * (part_production_last + part_cross_diffusion[last]) + Dw_last_plus_half * w_nodeUM;
            double Turbu_omega_new[ny]{};
            if (type_tdma == 0)
            {
                tdma(ny, a_w, b_w, c_w, d_w, Turbu_omega_new);
            }
            else if (type_tdma == 5)
            {
                tdma5(a_w, b_w, c_w, d_w, Turbu_omega_new);
            }
            else if (type_tdma == 10)
            {
                tdma10(a_w, b_w, c_w, d_w, Turbu_omega_new);
            }
            else
            {
                std::cout << "TDMA: Type not define! Stop here." << std::endl;
                std::cin.get();
            }
            double omega_min = 1e-10;
            for (int i = 0; i < ny; ++i) {
                Turbu_omega_new[i] = std::max(Turbu_omega_new[i], omega_min);
            }
            n_start = 0;
            for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_u) * u_star[i] + relax_u * U_new[i];
            n_start += ny;
            for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_k) * k_star[i] + relax_k * K_new[i];
            n_start += ny;
            for (int i = 0; i < ny; ++i) phi_solved[n_start + i] = (1.0 - relax_w) * turbu_omega_star[i] + relax_w * Turbu_omega_new[i];
            flow_rate_current = std::accumulate(U_new, U_new + ny, 0.0) * hy;
            flow_rate_current += u_nodeO * hy * 0.5;
            double ratio = flow_rate_target / (flow_rate_current + 1e-12);
            ratio = 0.5 * (ratio + std::sqrt(ratio * ratio + 1e-12));
            ratio = std::min(ratio, 10.0);
            double utau_new = utau * std::sqrt(ratio);
            utau = (1.0 - relax_utau) * utau + relax_utau * utau_new;
            if (!std::isfinite(utau))
            {
                utau = utau_init;
            }
            utau = std::max(utau, utau_min);
            if (!std::isfinite(utau))
            {
                utau = utau_init;
            }
            utau = std::max(utau, utau_min);
            differ = 0.0;
            for (int i = 0; i < 3 * ny; ++i) {
                double diff = phi_current[i] - phi_solved[i];
                differ += diff * diff;
            }
            differ = std::sqrt(differ);
            for (int i = 0; i < 3 * ny; ++i) {
                phi_current[i] = phi_solved[i];
            }
            num_iter_out += 1;
            int num_iter_out_limit = 100000;
            if (num_iter_out > num_iter_out_limit)
            {
                std::cout << "num_iter_out = " << num_iter_out << std::endl;
                std::cout << "Hard to achieve convergence in sublayer solver!" << std::endl;
                std::cout << "differ: " << differ << std::endl;
                std::cout << "flow_rate_current = " << flow_rate_current
                    << ", target = " << flow_rate_target << std::endl;
                std::cout << "updated utau = " << utau << std::endl;
                std::cout << "------------" << std::endl;
                std::cin.get();
                if (num_iter_out > 1.5 * num_iter_out_limit)
                {
                    std::cout << "Too many iterations, stop here!" << std::endl;
                    std::cin.get();
                }
            }
        }

        SublayerResult res;
        res.sublayer_utau = utau;

        for (int i = 0; i < ny; ++i) {
            res.sublayer_vel[i] = phi_solved[i];
            res.sublayer_k[i] = phi_solved[i + ny];
            res.sublayer_omega[i] = phi_solved[i + 2 * ny];
        }

        vel_nodeO = u_nodeO;
        vel_nodeUM = u_nodeUM;
        return res;
    }

    template <int Ny, int TypeTDMA>
    void P_refinement<Ny, TypeTDMA>::tdma(int N, const double* a, const double* b, const double* c, const double* d, double* x)
    {
        std::vector<double> cp(N, 0.0);
        std::vector<double> dp(N, 0.0);

        if (std::abs(b[0]) < 1e-14)
            throw std::runtime_error("TDMA: b[0] too small!");

        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];

        for (int i = 1; i < N; ++i)
        {
            double denom = b[i] - a[i] * cp[i - 1];
            if (std::abs(denom) < 1e-14)
                throw std::runtime_error("TDMA: denom too small");

            cp[i] = (i == N - 1) ? 0.0 : c[i] / denom;
            dp[i] = (d[i] - a[i] * dp[i - 1]) / denom;
        }

        x[N - 1] = dp[N - 1];
        for (int i = N - 2; i >= 0; --i)
        {
            x[i] = dp[i] - cp[i] * x[i + 1];
        }
    }

    template <int Ny, int TypeTDMA>
    void P_refinement<Ny, TypeTDMA>::tdma5(const double a[5], const double b[5], const double c[5], const double d[5], double x[5])
    {
        if (num_sub_node_ != 5)
        {
            std::cout << "TDMA5: Node number mismatch! Stop here." << std::endl;
            std::cin.get();
        }

        double cp[5]{ 0.0 };
        double dp[5]{ 0.0 };

        if (std::abs(b[0]) < 1e-14) throw std::runtime_error("TDMA5: b[0] too small!");
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];

        {
            double denom = b[1] - a[1] * cp[0];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA5: denom too small at row 1");
            cp[1] = c[1] / denom;
            dp[1] = (d[1] - a[1] * dp[0]) / denom;
        }

        {
            double denom = b[2] - a[2] * cp[1];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA5: denom too small at row 2");
            cp[2] = c[2] / denom;
            dp[2] = (d[2] - a[2] * dp[1]) / denom;
        }

        {
            double denom = b[3] - a[3] * cp[2];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA5: denom too small at row 3");
            cp[3] = c[3] / denom;
            dp[3] = (d[3] - a[3] * dp[2]) / denom;
        }

        {
            double denom = b[4] - a[4] * cp[3];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA5: denom too small at row 4");
            cp[4] = 0.0;
            dp[4] = (d[4] - a[4] * dp[3]) / denom;
        }

        x[4] = dp[4];
        x[3] = dp[3] - cp[3] * x[4];
        x[2] = dp[2] - cp[2] * x[3];
        x[1] = dp[1] - cp[1] * x[2];
        x[0] = dp[0] - cp[0] * x[1];
    }

    template <int Ny, int TypeTDMA>
    void P_refinement<Ny, TypeTDMA>::tdma10(const double a[10], const double b[10], const double c[10], const double d[10], double x[10])
    {
        if (num_sub_node_ != 10)
        {
            std::cout << "TDMA10: Node number mismatch! Stop here." << std::endl;
            std::cin.get();
        }

        double cp[10]{ 0.0 };
        double dp[10]{ 0.0 };

        if (std::abs(b[0]) < 1e-14) throw std::runtime_error("TDMA10: b[0] too small!");
        cp[0] = c[0] / b[0];
        dp[0] = d[0] / b[0];

        {
            double denom = b[1] - a[1] * cp[0];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 1");
            cp[1] = c[1] / denom;
            dp[1] = (d[1] - a[1] * dp[0]) / denom;
        }

        {
            double denom = b[2] - a[2] * cp[1];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 2");
            cp[2] = c[2] / denom;
            dp[2] = (d[2] - a[2] * dp[1]) / denom;
        }

        {
            double denom = b[3] - a[3] * cp[2];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 3");
            cp[3] = c[3] / denom;
            dp[3] = (d[3] - a[3] * dp[2]) / denom;
        }

        {
            double denom = b[4] - a[4] * cp[3];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 4");
            cp[4] = c[4] / denom;
            dp[4] = (d[4] - a[4] * dp[3]) / denom;
        }

        {
            double denom = b[5] - a[5] * cp[4];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 5");
            cp[5] = c[5] / denom;
            dp[5] = (d[5] - a[5] * dp[4]) / denom;
        }

        {
            double denom = b[6] - a[6] * cp[5];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 6");
            cp[6] = c[6] / denom;
            dp[6] = (d[6] - a[6] * dp[5]) / denom;
        }

        {
            double denom = b[7] - a[7] * cp[6];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 7");
            cp[7] = c[7] / denom;
            dp[7] = (d[7] - a[7] * dp[6]) / denom;
        }

        {
            double denom = b[8] - a[8] * cp[7];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 8");
            cp[8] = c[8] / denom;
            dp[8] = (d[8] - a[8] * dp[7]) / denom;
        }

        {
            double denom = b[9] - a[9] * cp[8];
            if (std::abs(denom) < 1e-14) throw std::runtime_error("TDMA10: denom too small at row 9");
            cp[9] = 0.0;
            dp[9] = (d[9] - a[9] * dp[8]) / denom;
        }

        x[9] = dp[9];
        x[8] = dp[8] - cp[8] * x[9];
        x[7] = dp[7] - cp[7] * x[8];
        x[6] = dp[6] - cp[6] * x[7];
        x[5] = dp[5] - cp[5] * x[6];
        x[4] = dp[4] - cp[4] * x[5];
        x[3] = dp[3] - cp[3] * x[4];
        x[2] = dp[2] - cp[2] * x[3];
        x[1] = dp[1] - cp[1] * x[2];
        x[0] = dp[0] - cp[0] * x[1];
    }

    template class P_refinement<5, 5>;
    template class P_refinement<10, 10>;
    template class P_refinement<5, 0>;
    template class P_refinement<10, 0>;
}

}

}
