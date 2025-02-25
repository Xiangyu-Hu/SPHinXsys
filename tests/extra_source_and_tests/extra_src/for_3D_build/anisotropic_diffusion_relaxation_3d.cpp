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
 
  
 
};
  
  
 
#endif // ANISOTROPIC_KERNEL_CORRECTION
