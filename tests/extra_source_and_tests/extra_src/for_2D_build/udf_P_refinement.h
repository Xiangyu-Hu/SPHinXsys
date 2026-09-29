#ifndef UDF_P_REFINEMENT_H
#define UDF_P_REFINEMENT_H

#include "udf_k-omega_turbulent_model.h"
#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
namespace udf
{

    template <typename... InteractionTypes>
    class P_refinement_GetVelocityGradient;

    template <class DataDelegationType>
    class P_refinement_GetVelocityGradient<DataDelegationType>
        : public LocalDynamics, public DataDelegationType
    {
    public:
        template <class BaseRelationType>
        explicit P_refinement_GetVelocityGradient(BaseRelationType& base_relation);
        virtual ~P_refinement_GetVelocityGradient() {};

    protected:
        Matd* velocity_gradient_only_P_;
        Vecd* k_gradient_only_P_;
        Vecd* omega_gradient_only_P_;

        Real* Vol_;
        Vecd* vel_;
        int* is_near_wall_P1_;
        int* is_near_wall_P2_;
        Real* turbu_k_;
        Real* turbu_omega_;
    };

    template <>
    class P_refinement_GetVelocityGradient<Inner<>> : public P_refinement_GetVelocityGradient<DataDelegateInner>
    {
    public:
        explicit P_refinement_GetVelocityGradient(BaseInnerRelation& inner_relation);
        virtual ~P_refinement_GetVelocityGradient() {};
        void interaction(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);

    protected:
        Matd* turbu_B_;
        Matd* B_;
    };
    using P_refinement_GetVelocityGradientInner = P_refinement_GetVelocityGradient<Inner<>>;

    template <>
    class P_refinement_GetVelocityGradient<Contact<Wall>> : public InteractionWithWall<P_refinement_GetVelocityGradient>
    {
    public:
        explicit P_refinement_GetVelocityGradient(BaseContactRelation& contact_relation);
        virtual ~P_refinement_GetVelocityGradient() {};
        void interaction(size_t index_i, Real dt = 0.0);

    protected:

    };

    using P_refinement_GetVelocityGradientComplex = ComplexInteraction<P_refinement_GetVelocityGradient<Inner<>, Contact<Wall>>>;

    template <int Ny = 5, int TypeTDMA = 0>
    class P_refinement : public LocalDynamics, public kOmega_BaseTurbuClosureCoeff, public WallFunctionCoefficient
    {
    public:
        explicit P_refinement(SPHBody& sph_body, Real constant_y_p, Real constant_y_p_node = 0.0);
        virtual ~P_refinement(){};

        void update(size_t index_i, Real dt = 0.0);

        void test_sublayer_model_half_channel_height();
        void test_sublayer_model_specific_channel_height();

        static constexpr int ny = Ny;
        static constexpr int type_tdma_ = TypeTDMA;

        using Vec_ny_d = Eigen::Matrix<Real, ny, 1>;
        inline void check_num_node_consistency()
        {
            node_dim_output_limit_ = std::remove_pointer_t<decltype(node_value_vel_)>::RowsAtCompileTime;
            if (ny > node_dim_output_limit_)
            {
                std::cout << "**** The number of sublayer nodes is "
                    << num_sub_node_ << ". ****" << std::endl;
                std::cout << "**** Node number is larger than the output limit. "
                    << "Only the two wall-nearest nodes and the last "
                    << node_dim_output_limit_ - 2
                    << "nodes toward U will be output. ****"
                    << std::endl;
                is_truncated_output_ = true;
                num_skip_output_ = num_sub_node_ - node_dim_output_limit_;
            }
            else
            {
                is_truncated_output_ = false;
                num_skip_output_ = 0;
            }
        }
        inline void construct_node_distribution(Real constant_y_p, Real constant_y_p_node)
        {
            sublayer_height_contant_ = constant_y_p;

            if (constant_y_p_node < TinyReal)
            {
                sublayer_y_p_constant_ = sublayer_height_contant_ / (Real(ny) + 0.5) / 2.0;
            }
            else
            {
                sublayer_y_p_constant_ = constant_y_p_node;
            }
            sublayer_node_uniform_distance_ = (sublayer_height_contant_ - sublayer_y_p_constant_) / Real(ny);
            for (int i = 0; i < ny; ++i)
            {
                sublayer_y_[i] = sublayer_y_p_constant_ + i * sublayer_node_uniform_distance_;
                std::cout << "**** sublayer_y_[" << i << "]= " << sublayer_y_[i] << ". ****" << std::endl;
            }
            std::cout << "**** sublayer_height_contant_= " << sublayer_height_contant_ << ". ****" << std::endl;
        }
        inline void output_node_distribution()
        {
            std::string filename = "../bin/output/Node_sublayer_y.dat";
            std::ofstream outfile(filename);
            if (!outfile.is_open())
            {
                std::cerr << "Error: Unable to open file " << filename << " for writing." << std::endl;
                return;
            }
            if (is_truncated_output_)
            {

                outfile << sublayer_y_[0] << "\n";

                outfile << sublayer_y_[1] << "\n";

                for (int i = 2; i < node_dim_output_limit_; ++i) {
                    outfile << sublayer_y_[i + num_skip_output_] << "\n";
                }
            }
            else
            {
                for (int i = 0; i < ny; ++i)
                {
                    outfile << sublayer_y_[i] << "\n";
                }
            }
            outfile.close();
        }

        struct SublayerResult {
            double sublayer_utau;
            Vec_ny_d sublayer_vel;
            Vec_ny_d sublayer_k;
            Vec_ny_d sublayer_omega;
        };

        SublayerResult solve_1D_sublayer_Dirichlet(double kinematic_viscosity, double u_p_outer, double k_p_outer,
            double w_p_outer, double vel_grad_p_outer, double nut_p_outer, double h_sublayer,
            double utau_outer, double Q_target, double k_grad_p_outer, double w_grad_p_outer, double& vel_nodeO, double& vel_nodeUM);
        void tdma(int N, const double* a, const double* b, const double* c, const double* d, double* x);
        void tdma5(const double a[5], const double b[5], const double c[5], const double d[5], double x[5]);
        void tdma10(const double a[10], const double b[10], const double c[10], const double d[10], double x[10]);

        inline Real get_loacal_flow_rate(Real average_flow_rate_over_particle_P, Real SPH_vel_grad_P, Real nodeO_U, Real dp)
        {
            Real nodeOS_U = nodeO_U + SPH_vel_grad_P * 0.5 * dp;
            Real flow_rate_half = (nodeO_U + nodeOS_U) * dp / 4.0;
            Real flow_rate_whole = average_flow_rate_over_particle_P;
            Real flow_rate_local = flow_rate_whole - flow_rate_half;
            return flow_rate_local;
        }

        inline Real obtainTangentialComponent(const Vecd& vec, const Vecd& normal)
        {
            Real dot_un = vec.dot(normal);
            Real norm_sqr = vec.dot(vec) - dot_un * dot_un;
            return std::sqrt(std::max(0.0, norm_sqr));
        }

    protected:
        int num_sub_node_, node_dim_output_limit_, num_skip_output_;
        Real sublayer_height_contant_;
        Vec_ny_d sublayer_y_;
        Real sublayer_y_p_constant_;
        Real sublayer_node_uniform_distance_;
        bool is_truncated_output_;
        bool output_detailed_info_sublayer_ = false;
        Real* friction_velocity_from_sublayer_;
        Real* target_flow_rate_in_sublayer_;
        Real* vel_ps_magnitude_;
        Real* dudn_for_local_flow_rate_;
        Real* utau_node_;
        Vec6d* node_value_vel_;
        Vec6d* node_value_k_;
        Real* dUdn_P_sublayer_magnitude_;
        Matd* dUdn_P_sublayer_;
        Real* vel_nodeO_;
        Real* vel_nodeUM_;
        Real* global_flow_rate_over_P_;
        Real* half_flow_rate_over_P_;
        int* is_near_wall_P1_;
        Vecd* vel_;
        Real* turbu_k_;
        Real* turbu_omega_;
        Real* rho_;
        Viscosity& viscosity_;
        Real mu_;
        Matd* velocity_gradient_only_P_;
        Real* turbu_mu_;
        Real* distance_to_dummy_interface_;
        Real* wall_shear_stress_;
        Vecd* e_nearest_normal_;
        Real fluid_particle_spacing_;
        Real* physical_time_;
        Matd* velocity_gradient_;
        Vecd* k_gradient_only_P_;
        Vecd* omega_gradient_only_P_;
    };

}
}
}
#endif