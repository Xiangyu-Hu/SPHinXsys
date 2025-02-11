/**
 * @file anisotropic_diffusion_relaxation.hpp
 * @brief This the methods related on the corrections of SPH anisotropic smoothing kernel and anisotropic diffusion algorithm.
 * @details Two template are included, one is the improved kernel matrix, and the other one is anisotropic diffusion relaxation.
 * @author Xiaojing Tang and Xiangyu Hu
 */

#ifndef ANISOTROPIC
#define ANISOTROPIC

#include "base_general_dynamics.h"
#include "sphinxsys.h"  
namespace SPH
{


//TODO:: The class name. The datatype should be applicable to 6d.

//material definition
class AnisotropicDiffusionSolid : public Solid
{
  public:
    AnisotropicDiffusionSolid(Real rho0, Real coeff)
        : Solid(rho0), diffusion_coeff(coeff)
    {
        material_type_name_ = "AnisotropicDiffusionSolid";
    };
    virtual ~AnisotropicDiffusionSolid(){};

    Real diffusion_coeff;
    Real DiffusivityCoefficient() { return diffusion_coeff; };
};
 
//----------------------------------------------------------------------
//	calculate correction matrix B to keep the accuracy in anisotropic case
//----------------------------------------------------------------------
template <typename... InteractionTypes>
class AnisotropicLinearGradientCorrectionMatrix;

template <class DataDelegationType>
class AnisotropicLinearGradientCorrectionMatrix<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit AnisotropicLinearGradientCorrectionMatrix(BaseRelationType &base_relation)
     : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
       Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
       B_(particles_->registerStateVariable<Matd>("LinearGradientCorrectionMatrix", IdentityMatrix<Matd>::value)),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position"))  {}

    virtual ~AnisotropicLinearGradientCorrectionMatrix(){};

  protected:
    Real *Vol_;
    Matd *B_;
    Vecd *pos_;

};

template <>
class AnisotropicLinearGradientCorrectionMatrix<Inner<>>
    : public AnisotropicLinearGradientCorrectionMatrix<DataDelegateInner>
{

  public:
    explicit AnisotropicLinearGradientCorrectionMatrix(BaseInnerRelation &inner_relation)
        : AnisotropicLinearGradientCorrectionMatrix<DataDelegateInner>(inner_relation)
        {};
    template <typename BodyRelationType, typename FirstArg>
    explicit AnisotropicLinearGradientCorrectionMatrix(InteractArgs<BodyRelationType, FirstArg> parameters)
        : AnisotropicLinearGradientCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~AnisotropicLinearGradientCorrectionMatrix(){};


    void interaction(size_t index_i, Real dt = 0.0)
    {
        Matd local_configuration = ZeroData<Matd>::value;
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            Vecd r_ji = pos_[index_i] - pos_[index_j];
            local_configuration -= r_ji * gradW_ij.transpose();
        }
        B_[index_i] = local_configuration;
    };
    void update(size_t index_i, Real dt = 0.0)
    {
        Matd B_T = B_[index_i].transpose();  
        Matd inverse = (B_T * B_[index_i] + SqrtEps * Matd::Identity()).inverse() * B_T;
        B_[index_i] = inverse ;
    };
    

};
using AnisotropicLinearGradientCorrectionMatrixInner = AnisotropicLinearGradientCorrectionMatrix<Inner<>>;

template <>
class AnisotropicLinearGradientCorrectionMatrix<Contact<>>
    : public AnisotropicLinearGradientCorrectionMatrix<DataDelegateContact>
{
  public:
    explicit AnisotropicLinearGradientCorrectionMatrix(BaseContactRelation &contact_relation)
    : AnisotropicLinearGradientCorrectionMatrix<DataDelegateContact>(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_pos_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Position"));
      
    }
}virtual ~AnisotropicLinearGradientCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0)
    {
    Matd local_configuration = ZeroData<Matd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real *Vol_k = contact_Vol_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Vecd gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
            Vecd r_ji = pos_[index_i] - contact_pos_[k][index_j];
            local_configuration -= r_ji * gradW_ij.transpose();
        }
    }
    B_[index_i] += local_configuration;
};

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_pos_;
    
};

using AnisotropicLinearGradientCorrectionMatrixComplex = ComplexInteraction<AnisotropicLinearGradientCorrectionMatrix<Inner<>, Contact<>>>;



template <typename... InteractionTypes>
class AnisotropicKernelCorrectionMatrix;

template <class DataDelegationType>
class AnisotropicKernelCorrectionMatrix<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit AnisotropicKernelCorrectionMatrix(BaseRelationType &base_relation)
      : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
     kernel_correction1_(this->particles_->template registerStateVariable<Vec2d>("FirstOrderCorrectionVector1")), 
      kernel_correction2_(this->particles_->template registerStateVariable<Vec2d>("FirstOrderCorrectionVector2")),
      kernel_correction3_(this->particles_->template registerStateVariable<Vec2d>("FirstOrderCorrectionVector3")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),  
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position"))  {}
          
    virtual ~AnisotropicKernelCorrectionMatrix(){};

  protected:
   
    Vec2d *kernel_correction1_,*kernel_correction2_, *kernel_correction3_;
    Real *Vol_;
    Mat2d *B_;
    Vecd *pos_;

}; 

template <>
class AnisotropicKernelCorrectionMatrix<Inner<>>
    : public AnisotropicKernelCorrectionMatrix<DataDelegateInner>
{
  public:
    explicit AnisotropicKernelCorrectionMatrix(BaseInnerRelation &inner_relation)
        : AnisotropicKernelCorrectionMatrix<DataDelegateInner>(inner_relation){};
    template <typename BodyRelationType, typename FirstArg>
    explicit AnisotropicKernelCorrectionMatrix(InteractArgs<BodyRelationType, FirstArg> parameters)
        : AnisotropicKernelCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~AnisotropicKernelCorrectionMatrix(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ik
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
         
            // Vec2d r_ji = -inner_neighborhood.r_ij_vector_[n];
            Vecd r_ji = pos_[index_j] - pos_[index_i];
            kernel_correction1_[index_i] += r_ji[0] * r_ji[0] * (B_[index_i].transpose() * gradW_ij);
            kernel_correction2_[index_i] += r_ji[1] * r_ji[1] * (B_[index_i].transpose() * gradW_ij);
            kernel_correction3_[index_i] += r_ji[0] * r_ji[1] * (B_[index_i].transpose() * gradW_ij);
         
        }
       
    };
    void update(size_t index_i, Real dt = 0.0){};
};
  
template <>
class AnisotropicKernelCorrectionMatrix<Contact<>>
    : public AnisotropicKernelCorrectionMatrix<DataDelegateContact>
{
  public:
    explicit AnisotropicKernelCorrectionMatrix(BaseContactRelation &contact_relation)
    : AnisotropicKernelCorrectionMatrix<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
            contact_pos_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Position"));
        }
    };

    virtual ~AnisotropicKernelCorrectionMatrix(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
               // Vec2d r_ji = -contact_neighborhood.r_ij_vector_[n];
                 Vecd r_ji = contact_pos_[k][index_j] - pos_[index_i];

                kernel_correction1_[index_i] += r_ji[0] * r_ji[0] * (B_[index_i].transpose() * gradW_ij);
                kernel_correction2_[index_i] += r_ji[1] * r_ji[1] * (B_[index_i].transpose() * gradW_ij);
                kernel_correction3_[index_i] += r_ji[0] * r_ji[1] * (B_[index_i].transpose() * gradW_ij);
            }
        }

    };

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_pos_;
};

using AnisotropicKernelCorrectionMatrixComplex = ComplexInteraction<AnisotropicKernelCorrectionMatrix<Inner<>, Contact<>>>;


/*

Class Anisotropic diffusion relaxation

*/

template <typename... InteractionTypes>
class AnisotropicDiffusionRelaxation;

template <class DataDelegationType>
class AnisotropicDiffusionRelaxation<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit AnisotropicDiffusionRelaxation(BaseRelationType &base_relation)
      : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation), 
      B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      phi_(particles_->registerStateVariable<Real>("Phi")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      kernel_correction1_(this->particles_->template getVariableDataByName<Vec2d>("FirstOrderCorrectionVector1")),
      kernel_correction2_(this->particles_->template getVariableDataByName<Vec2d>("FirstOrderCorrectionVector2")),
      kernel_correction3_(this->particles_->template getVariableDataByName<Vec2d>("FirstOrderCorrectionVector3")),
       anisotropic_diffusion_solid_(DynamicCast<AnisotropicDiffusionSolid>(this, base_relation.getSPHBody().getBaseMaterial())), 
       pos_(this->particles_->template getVariableDataByName<Vecd>("Position")) 
      { 
        
        total_left_ = particles_->registerStateVariable<Mat3d>( "TotalLeft", [&](size_t i) -> Mat3d { return Eps * Mat3d::Identity(); });
        species_correction_ = particles_->registerStateVariable<Vec2d>("SpeciesCorrection ", [&](size_t i) -> Vec2d { return Eps * Vec2d::Identity(); });
        
        total_right_ = particles_->registerStateVariable<Vec3d>( "TotalRight", [&](size_t i) -> Vec3d { return Eps * Vec3d::Identity(); });
        Laplacian_=  particles_->registerStateVariable<Vec3d>( "Laplacian", [&](size_t i) -> Vec3d { return Vec3d::Zero(); });

        Laplacian_x =   particles_->registerStateVariable<Real>("Laplacian_x", [&](size_t i) -> Real { return Real(0.0); });
        Laplacian_y = particles_->registerStateVariable<Real>("Laplacian_y", [&](size_t i) -> Real { return Real(0.0); });
        Laplacian_xy = particles_->registerStateVariable<Real>( "Laplacian_xy", [&](size_t i) -> Real { return Real(0.0); });
        diffusion_dt_=particles_->registerStateVariable<Real>( "diffusion_dt", [&](size_t i) -> Real { return Real(0.0); });
		
        diffusion_coeff_ = anisotropic_diffusion_solid_.DiffusivityCoefficient();

      
    }
          
    virtual ~AnisotropicDiffusionRelaxation(){};

  protected: 
    Mat2d *B_;
    Real  *phi_,*Vol_;
    Vec2d *kernel_correction1_, *kernel_correction2_, *kernel_correction3_;
    
    Mat3d *total_left_;
    Vec2d *species_correction_;
    Vec3d *total_right_;
    Vec3d *Laplacian_;

    Real *Laplacian_x, *Laplacian_y, *Laplacian_xy, *diffusion_dt_;
   
    Real diffusion_coeff_;
 

    AnisotropicDiffusionSolid &anisotropic_diffusion_solid_;
      Vecd *pos_;
};

template <>
class AnisotropicDiffusionRelaxation<Inner<>>
    : public AnisotropicDiffusionRelaxation<DataDelegateInner>
{
  public:
    explicit AnisotropicDiffusionRelaxation(BaseInnerRelation &inner_relation)
        : AnisotropicDiffusionRelaxation<DataDelegateInner>(inner_relation){};

    template <typename BodyRelationType, typename FirstArg>
    explicit AnisotropicDiffusionRelaxation(InteractArgs<BodyRelationType, FirstArg> parameters)
        : AnisotropicDiffusionRelaxation(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~AnisotropicDiffusionRelaxation(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        Vec2d species_correction_rate = Vec2d::Zero();
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) 
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
           species_correction_rate += (phi_[index_j] - phi_[index_i]) * (B_[index_i].transpose() * gradW_ij);      
        }

         species_correction_[index_i] = species_correction_rate;  
         Vec3d total_right_rate = Vec3d::Zero();
         Mat3d total_left_rate = Mat3d::Zero();
         Real modified_func_  = 1.0;

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)  
        {
            size_t index_j = inner_neighborhood.j_[n];
           // Vec2d r_ji = -inner_neighborhood.r_ij_vector_[n];
            Vecd r_ji = pos_[index_j] - pos_[index_i];
            Vec2d gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            Vec3d disp_quad_ = Vec3d(r_ji[0] * r_ji[0], r_ji[1] * r_ji[1], r_ji[0] * r_ji[1]);
            modified_func_  = r_ji.dot(B_[index_i].transpose() * gradW_ijV_j) / pow(r_ji.norm(), 4.0);
	
            Real right_ = 2.0 * (phi_[index_j] - phi_[index_i] - r_ji.dot(species_correction_[index_i]));
            total_right_rate += disp_quad_ * modified_func_  * right_;
             
            Vec3d left_ = Vec3d::Zero();
            left_[0] = (r_ji[0] * r_ji[0]- r_ji.dot(kernel_correction1_[index_i]));
            left_[1] = (r_ji[1] * r_ji[1]- r_ji.dot(kernel_correction2_[index_i]));
            left_[2] = (r_ji[0] * r_ji[1]- r_ji.dot(kernel_correction3_[index_i]));
            total_left_rate += disp_quad_ * modified_func_ * left_.transpose();   

        }
        
        total_right_[index_i] = total_right_rate;
        total_left_[index_i] = total_left_rate;
    };
 
    void update(size_t index_i, Real dt = 0.0)
    {
        phi_[index_i] += dt * diffusion_dt_[index_i];
    };

};
 

template <>
class AnisotropicDiffusionRelaxation<Contact<>>
    : public AnisotropicDiffusionRelaxation<DataDelegateContact>
{
  public:
    explicit AnisotropicDiffusionRelaxation(BaseContactRelation &contact_relation)
    : AnisotropicDiffusionRelaxation<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
            contact_pos_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Position"));
      
        }
    };

    virtual ~AnisotropicDiffusionRelaxation(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Mat3d total_left_rate_contact = Mat3d::Zero();
        Vec3d total_right_rate_contact = Vec3d::Zero();
        Real modified_func_contact = 1.0;
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                 size_t index_j = contact_neighborhood.j_[n];
                 //Vec2d r_ji = -contact_neighborhood.r_ij_vector_[n];
                 Vecd r_ji = contact_pos_[k][index_j] - pos_[index_i];

                 Vecd gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];

                Vec3d disp_quad_ = Vec3d(r_ji[0] * r_ji[0], r_ji[1] * r_ji[1], r_ji[0] * r_ji[1]);
                Real right_ = 2.0 * (0.0 - r_ji.dot(species_correction_[index_i])); ///here when it is periodic boundary condition, should notice the 0.0
                modified_func_contact = r_ji.dot(B_[index_i].transpose() * gradW_ij) / pow(r_ji.norm(), 4.0);
		 
                Vec3d left_ = Vec3d::Zero();
                left_[0] = (r_ji[0] * r_ji[0]- r_ji.dot(kernel_correction1_[index_i]));
                left_[1] = (r_ji[1] * r_ji[1]- r_ji.dot(kernel_correction2_[index_i]));
                left_[2] = (r_ji[0] * r_ji[1]- r_ji.dot(kernel_correction3_[index_i]));

                total_left_rate_contact += disp_quad_ * modified_func_contact * left_.transpose();
                total_right_rate_contact += disp_quad_ * modified_func_contact * right_;
            
            }
            total_left_[index_i] += total_left_rate_contact;
            total_right_[index_i] += total_right_rate_contact;
        }
        Laplacian_[index_i] = diffusion_coeff_ * total_left_[index_i].inverse() * total_right_[index_i];

        Laplacian_x[index_i] = Laplacian_[index_i][0];
        Laplacian_y[index_i] = Laplacian_[index_i][1];
        Laplacian_xy[index_i] = Laplacian_[index_i][2];
		diffusion_dt_[index_i] = Laplacian_[index_i][0] + Laplacian_[index_i][1];

    };

  protected:
    StdVec<Real *> contact_Vol_; 
    StdVec<Vecd *> contact_pos_;
};

 
using AnisotropicDiffusionRelaxationComplex = ComplexInteraction<AnisotropicDiffusionRelaxation<Inner<>, Contact<>>>;

 

}
#endif // ANISOTROPIC
