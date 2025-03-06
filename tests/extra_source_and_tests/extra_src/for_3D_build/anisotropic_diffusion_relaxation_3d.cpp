/**
 * @file anisotropic_diffusion_relaxation.hpp
 * @brief This the methods related on the corrections of SPH anisotropic smoothing kernel and anisotropic diffusion algorithm.
 * @details Two template are included, one is the improved kernel matrix, and the other one is anisotropic diffusion relaxation.
 * @author Xiaojing Tang and Xiangyu Hu
 */

#ifndef ANISOTROPIC_KERNEL_CORRECTION
#define ANISOTROPIC_KERNEL_CORRECTION
 
#include "anisotropic_diffusion_relaxation.hpp"
namespace SPH
{
     
    void AnisotropicKernelCorrectionMatrix<Inner<>>::kernel_correction_function(Vec3d r_ji, Vec3d grad_, size_t index_i)
    {  
        kernel_correction1_[index_i] += r_ji[0] * r_ji[0] * grad_;
        kernel_correction2_[index_i] += r_ji[1] * r_ji[1] * grad_;
        kernel_correction3_[index_i] += r_ji[2] * r_ji[2] * grad_;

        kernel_correction4_[index_i] += r_ji[0] * r_ji[1] * grad_;
        kernel_correction5_[index_i] += r_ji[1] * r_ji[2] * grad_;
        kernel_correction6_[index_i] += r_ji[2] * r_ji[0] * grad_;

    }; 

    void AnisotropicKernelCorrectionMatrix<Contact<>>::kernel_correction_function_contact(Vec3d r_ji, Vec3d grad_, size_t index_i)
    {  
        
        kernel_correction1_[index_i] += r_ji[0] * r_ji[0] * grad_;
        kernel_correction2_[index_i] += r_ji[1] * r_ji[1] * grad_;
        kernel_correction3_[index_i] += r_ji[2] * r_ji[2] * grad_;

        kernel_correction4_[index_i] += r_ji[0] * r_ji[1] * grad_;
        kernel_correction5_[index_i] += r_ji[1] * r_ji[2] * grad_;
        kernel_correction6_[index_i] += r_ji[2] * r_ji[0] * grad_;
 

    };
    

    void AnisotropicDiffusionRelaxation<Inner<>>::relaxation(const Neighborhood &inner_neighborhood,Vec3d dimension, size_t index_i )
    {  
      
        Vec6d total_right_rate = Vec6d::Zero();
        Mat6d total_left_rate = Mat6d::Zero();
 
       for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)  
       {
           size_t index_j = inner_neighborhood.j_[n]; 
           Vec3d r_ji = pos_[index_j] - pos_[index_i];
           Vec3d gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
           Real modified_func_  = r_ji.dot(B_[index_i].transpose() * gradW_ijV_j) / pow(r_ji.norm(), 4.0);
          
           Real right_ = 2.0 * (species_[index_j] - species_[index_i] - r_ji.dot(species_correction_[index_i]));
        
          
           Vec6d disp_quad_ =  Vec6d(r_ji[0] * r_ji[0], r_ji[1] * r_ji[1], r_ji[2] * r_ji[2], r_ji[0] * r_ji[1], r_ji[1] * r_ji[2], r_ji[2] * r_ji[0]);
           Vec6d left_ = Vec6d::Zero();
           left_[0] = (r_ji[0] * r_ji[0]- r_ji.dot(kernel_correction1_[index_i]));
           left_[1] = (r_ji[1] * r_ji[1]- r_ji.dot(kernel_correction2_[index_i]));
           left_[2] = (r_ji[2] * r_ji[2]- r_ji.dot(kernel_correction3_[index_i]));
           left_[3] = (r_ji[0] * r_ji[1]- r_ji.dot(kernel_correction4_[index_i]));
           left_[4] = (r_ji[1] * r_ji[2]- r_ji.dot(kernel_correction5_[index_i]));
           left_[5] = (r_ji[2] * r_ji[0]- r_ji.dot(kernel_correction6_[index_i]));
        

           total_left_rate += disp_quad_ * modified_func_ * left_.transpose();   
           total_right_rate += disp_quad_ * modified_func_  * right_;

       }
 
       total_right_3d[index_i] = total_right_rate;
       total_left_3d[index_i] = total_left_rate;

    };

    void AnisotropicDiffusionRelaxation<Contact<>>::relaxation_contact(Vec3d dimension, size_t index_i )
    {
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {   

            Real modified_func_contact = 1.0;
            Mat6d total_left_3d_rate_contact = Mat6d::Zero();
            Vec6d total_right_3d_rate_contact = Vec6d::Zero();
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n]; 
                Vec3d r_ji = contact_pos_[k][index_j] - pos_[index_i];
                Vec3d gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];

                Vec6d disp_quad_ =  Vec6d(r_ji[0] * r_ji[0], r_ji[1] * r_ji[1], r_ji[2] * r_ji[2], r_ji[0] * r_ji[1], r_ji[1] * r_ji[2], r_ji[2] * r_ji[0]);
                Real right_ = 2.0 * (0.0 - r_ji.dot(species_correction_[index_i]));  
                modified_func_contact = r_ji.dot(B_[index_i].transpose() * gradW_ij) / pow(r_ji.norm(), 4.0);
        
                Vec6d left_ = Vec6d::Zero();
                left_[0] = (r_ji[0] * r_ji[0]- r_ji.dot(kernel_correction1_[index_i]));
                left_[1] = (r_ji[1] * r_ji[1]- r_ji.dot(kernel_correction2_[index_i]));
                left_[2] = (r_ji[2] * r_ji[2]- r_ji.dot(kernel_correction3_[index_i]));
                left_[3] = (r_ji[0] * r_ji[1]- r_ji.dot(kernel_correction4_[index_i]));
                left_[4] = (r_ji[1] * r_ji[2]- r_ji.dot(kernel_correction5_[index_i]));
                left_[5] = (r_ji[2] * r_ji[0]- r_ji.dot(kernel_correction6_[index_i]));
             

                total_left_3d_rate_contact += disp_quad_ * modified_func_contact * left_.transpose();
                total_right_3d_rate_contact += disp_quad_ * modified_func_contact * right_;
            
            }
            total_left_3d[index_i] += total_left_3d_rate_contact;
            total_right_3d[index_i] += total_right_3d_rate_contact;
       }       
   
        Laplacian_3d[index_i] = diffusion_coeff_ * total_left_3d[index_i].inverse() * total_right_3d[index_i];

        Laplacian_x[index_i] = Laplacian_3d[index_i][0];
        Laplacian_y[index_i] = Laplacian_3d[index_i][1]; 
        Laplacian_z[index_i] = Laplacian_3d[index_i][2]; 
		diffusion_dt_[index_i] = Laplacian_3d[index_i][0] + Laplacian_3d[index_i][1] + Laplacian_3d[index_i][2];
      
    }
  
 
};
  
  
 
#endif // ANISOTROPIC_KERNEL_CORRECTION
