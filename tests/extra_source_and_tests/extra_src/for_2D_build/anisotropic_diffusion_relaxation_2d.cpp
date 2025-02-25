/**
 * @file anisotropic_diffusion_relaxation.hpp
 * @brief This the methods related on the corrections of SPH anisotropic smoothing kernel and anisotropic diffusion algorithm.
 * @details Two template are included, one is the improved kernel matrix, and the other one is anisotropic diffusion relaxation.
 * @author Xiaojing Tang and Xiangyu Hu
 */

#ifndef ANISOTROPIC_KERNEL
#define ANISOTROPIC_KERNEL
 
#include "anisotropic_diffusion_relaxation.hpp"
namespace SPH
{
     
    void AnisotropicKernelCorrectionMatrix<Inner<>>::kernel_correction_function(Vec2d r_ji, Vec2d grad_, size_t index_i)
    {  
        kernel_correction1_[index_i] += r_ji[0] * r_ji[0] * grad_;
        kernel_correction2_[index_i] += r_ji[1] * r_ji[1] * grad_;
        kernel_correction3_[index_i] += r_ji[0] * r_ji[1] * grad_;

    };

    void AnisotropicKernelCorrectionMatrix<Contact<>>::kernel_correction_function_contact(Vec2d r_ji, Vec2d grad_, size_t index_i)
    {  
        kernel_correction1_[index_i] += r_ji[0] * r_ji[0] * grad_;
        kernel_correction2_[index_i] += r_ji[1] * r_ji[1] * grad_;
        kernel_correction3_[index_i] += r_ji[0] * r_ji[1] * grad_; 
 

    };
 
    void AnisotropicDiffusionRelaxation<Inner<>>::relaxation(const Neighborhood &inner_neighborhood, Vec2d dimension, size_t index_i )
    {  
        Vec3d total_right_rate = Vec3d::Zero();
        Mat3d total_left_rate = Mat3d::Zero();
 
       for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)  
       {
           size_t index_j = inner_neighborhood.j_[n]; 
           Vecd r_ji = pos_[index_j] - pos_[index_i];
           Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
           Real modified_func_  = r_ji.dot(B_[index_i].transpose() * gradW_ijV_j) / pow(r_ji.norm(), 4.0);
           Real right_ = 2.0 * (phi_[index_j] - phi_[index_i] - r_ji.dot(species_correction_[index_i]));
        
          
           Vec3d disp_quad_ = Vec3d(r_ji[0] * r_ji[0], r_ji[1] * r_ji[1], r_ji[0] * r_ji[1]);
           Vec3d left_ = Vec3d::Zero();
           left_[0] = (r_ji[0] * r_ji[0]- r_ji.dot(kernel_correction1_[index_i]));
           left_[1] = (r_ji[1] * r_ji[1]- r_ji.dot(kernel_correction2_[index_i]));
           left_[2] = (r_ji[0] * r_ji[1]- r_ji.dot(kernel_correction3_[index_i]));

           total_left_rate += disp_quad_ * modified_func_ * left_.transpose();   
           total_right_rate += disp_quad_ * modified_func_  * right_;

       }
 
       total_right_2d[index_i] = total_right_rate;
       total_left_2d[index_i] = total_left_rate;

    };

    
    void AnisotropicDiffusionRelaxation<Contact<>>::relaxation_contact(Vec2d dimension, size_t index_i,size_t k )
    {

        Real modified_func_contact = 1.0;
        Mat3d total_left_2d_rate_contact = Mat3d::Zero();
        Vec3d total_right_2drate_contact = Vec3d::Zero();
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
             size_t index_j = contact_neighborhood.j_[n]; 
             Vecd r_ji = contact_pos_[k][index_j] - pos_[index_i];
             Vecd gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];

            Vec3d disp_quad_ = Vec3d(r_ji[0] * r_ji[0], r_ji[1] * r_ji[1], r_ji[0] * r_ji[1]);
            Real right_ = 2.0 * (0.0 - r_ji.dot(species_correction_[index_i]));  
            modified_func_contact = r_ji.dot(B_[index_i].transpose() * gradW_ij) / pow(r_ji.norm(), 4.0);
     
            Vec3d left_ = Vec3d::Zero();
            left_[0] = (r_ji[0] * r_ji[0]- r_ji.dot(kernel_correction1_[index_i]));
            left_[1] = (r_ji[1] * r_ji[1]- r_ji.dot(kernel_correction2_[index_i]));
            left_[2] = (r_ji[0] * r_ji[1]- r_ji.dot(kernel_correction3_[index_i]));

            total_left_2d_rate_contact += disp_quad_ * modified_func_contact * left_.transpose();
            total_right_2drate_contact += disp_quad_ * modified_func_contact * right_;
        
        }
        total_left_2d[index_i] += total_left_2d_rate_contact;
        total_right_2d[index_i] += total_right_2drate_contact;
    }
 
};
  
  
 
#endif // ANISOTROPIC_KERNEL
