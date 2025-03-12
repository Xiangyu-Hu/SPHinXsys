/**
 * @file anisotropic_diffusion_relaxation.hpp
 * @brief This the methods related on the corrections of SPH anisotropic smoothing kernel and anisotropic diffusion algorithm.
 * @details Two template are included, one is the improved kernel matrix, and the other one is anisotropic diffusion relaxation.
 * @author Xiaojing Tang and Xiangyu Hu
 */

#ifndef ANISOTROPIC
#define ANISOTROPIC

#include "base_data_type.h"
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
    explicit AnisotropicLinearGradientCorrectionMatrix(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : AnisotropicLinearGradientCorrectionMatrix(parameters.identifier_, std::get<0>(parameters.others_)){};
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
      kernel_correction1_(this->particles_->template registerStateVariable<Vecd>("FirstOrderCorrectionVector1")), 
      kernel_correction2_(this->particles_->template registerStateVariable<Vecd>("FirstOrderCorrectionVector2")),
      kernel_correction3_(this->particles_->template registerStateVariable<Vecd>("FirstOrderCorrectionVector3")),
      kernel_correction4_(this->particles_->template registerStateVariable<Vecd>("FirstOrderCorrectionVector4")),
      kernel_correction5_(this->particles_->template registerStateVariable<Vecd>("FirstOrderCorrectionVector5")),
      kernel_correction6_(this->particles_->template registerStateVariable<Vecd>("FirstOrderCorrectionVector6")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),  
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position"))  {}
          
    virtual ~AnisotropicKernelCorrectionMatrix(){};

  protected:
   
    Vecd *kernel_correction1_,*kernel_correction2_, *kernel_correction3_;// for 2d case, and first three variables in 3d case
    Vecd *kernel_correction4_,*kernel_correction5_, *kernel_correction6_; // another three variables in 3d case 
    Real *Vol_;
    Matd *B_;
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
    explicit AnisotropicKernelCorrectionMatrix(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : AnisotropicKernelCorrectionMatrix(parameters.identifier_, std::get<0>(parameters.others_)){};
    
    virtual ~AnisotropicKernelCorrectionMatrix(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ik
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
         
            Vecd r_ji = pos_[index_j] - pos_[index_i];   
            Vecd grad_modified = B_[index_i].transpose() * gradW_ij;
            kernel_correction_function(r_ji, grad_modified, index_i);
                
        }
       
    };
    void update(size_t index_i, Real dt = 0.0){};

    protected:
    void kernel_correction_function(Vecd r_ji, Vecd grad_, size_t index_i);
     
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
                Vecd r_ji = contact_pos_[k][index_j] - pos_[index_i];
               
                         
                Vecd grad_modified = B_[index_i].transpose() * gradW_ij;
                kernel_correction_function_contact(r_ji, grad_modified, index_i);
                
              
            }
        }

    };

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_pos_;

  
    void kernel_correction_function_contact(Vecd r_ji, Vecd grad_, size_t index_i);
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
    explicit AnisotropicDiffusionRelaxation(BaseRelationType &base_relation )
      : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation), 
      B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      species_(particles_->registerStateVariable<Real>("Species")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      kernel_correction1_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector1")),
      kernel_correction2_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector2")),
      kernel_correction3_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector3")),
      kernel_correction4_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector4")),
      kernel_correction5_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector5")),
      kernel_correction6_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector6")),
      anisotropic_diffusion_solid_(DynamicCast<AnisotropicDiffusionSolid>(this, base_relation.getSPHBody().getBaseMaterial())), 
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")) 
      { 
        
        species_correction_ = particles_->registerStateVariable<Vecd>("SpeciesCorrection ", [&](size_t i) -> Vecd { return Eps * Vecd::Identity(); });
        total_left_2d = particles_->registerStateVariable<Mat3d>( "TotalLeft2D", [&](size_t i) -> Mat3d { return Eps * Mat3d::Identity(); });
        total_right_2d = particles_->registerStateVariable<Vec3d>( "TotalRight2D", [&](size_t i) -> Vec3d { return Eps * Vec3d::Identity(); });
        Laplacian_2d =  particles_->registerStateVariable<Vec3d>( "Laplacian2D", [&](size_t i) -> Vec3d { return Vec3d::Zero(); });
        

        total_left_3d = particles_->registerStateVariable<Mat6d>( "TotalLeft3D", [&](size_t i) -> Mat6d { return Eps * Mat6d::Identity(); });
        total_right_3d = particles_->registerStateVariable<Vec6d>( "TotalRight3D", [&](size_t i) -> Vec6d { return Eps * Vec6d::Identity(); });
        Laplacian_3d =  particles_->registerStateVariable<Vec6d>( "Laplacian3D", [&](size_t i) -> Vec6d { return Vec6d::Zero(); });

        Laplacian_x =   particles_->registerStateVariable<Real>("Laplacian_x", [&](size_t i) -> Real { return Real(0.0); });
        Laplacian_y = particles_->registerStateVariable<Real>("Laplacian_y", [&](size_t i) -> Real { return Real(0.0); });
        Laplacian_z = particles_->registerStateVariable<Real>( "Laplacian_z", [&](size_t i) -> Real { return Real(0.0); });
        diffusion_dt_=particles_->registerStateVariable<Real>( "diffusion_dt", [&](size_t i) -> Real { return Real(0.0); });
		
        diffusion_coeff_ = anisotropic_diffusion_solid_.DiffusivityCoefficient();

      
    }
          
    virtual ~AnisotropicDiffusionRelaxation(){};

  protected: 
    Matd *B_;
    Real  *species_,*Vol_;
    Vecd *kernel_correction1_, *kernel_correction2_, *kernel_correction3_;
    Vecd *kernel_correction4_, *kernel_correction5_, *kernel_correction6_;
    Vecd *species_correction_;

    Mat3d *total_left_2d;
    Vec3d *total_right_2d;
    Vec3d *Laplacian_2d;

    Mat6d *total_left_3d;
    Vec6d *total_right_3d;
    Vec6d *Laplacian_3d;


    Real *Laplacian_x, *Laplacian_y, *Laplacian_z, *diffusion_dt_;
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
    explicit AnisotropicDiffusionRelaxation(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : AnisotropicDiffusionRelaxation(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~AnisotropicDiffusionRelaxation(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        Vecd species_correction_rate = Vecd::Zero();
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) 
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
           species_correction_rate += (species_[index_j] - species_[index_i]) * (B_[index_i].transpose() * gradW_ij);      
        }

         species_correction_[index_i] = species_correction_rate;  
 
         relaxation(inner_neighborhood, Vecd::Identity(), index_i );
 
    };
 
    void update(size_t index_i, Real dt = 0.0)
    {
        species_[index_i] += dt * diffusion_dt_[index_i];
    };
   
    protected:
    void relaxation(const Neighborhood &inner_neighborhood, Vecd dimension,size_t index_i);
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
        
            relaxation_contact(Vecd::Identity(), index_i);      

    };

  protected:
    StdVec<Real *> contact_Vol_; 
    StdVec<Vecd *> contact_pos_;

    void relaxation_contact(Vecd dimension, size_t index_i );
};

 
using AnisotropicDiffusionRelaxationComplex = ComplexInteraction<AnisotropicDiffusionRelaxation<Inner<>, Contact<>>>;

 

}
#endif // ANISOTROPIC
