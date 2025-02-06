/**
 * @file anisotropic_relaxation.hpp
 * @brief This the methods related on the corrections of SPH smoothing kernel.
 * @details The corrections aim to increase the numerical consistency
 * or accuracy for kernel approximations.
 * @author Xiaojing Tang and Xiangyu Hu
 */

#ifndef ANISOTROPIC
#define ANISOTROPIC

#include "base_general_dynamics.h"
namespace SPH
{
template <typename... InteractionTypes>
class NonisotropicKernelCorrectionMatrixAC;

template <class DataDelegationType>
class NonisotropicKernelCorrectionMatrixAC<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit NonisotropicKernelCorrectionMatrixAC(BaseRelationType &base_relation)
      : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      phi_(particles_->registerStateVariable<Real>("Phi")),
      A1_(this->particles_->template registerStateVariable<Vec2d>(
          "FirstOrderCorrectionVectorA1")),
      A2_(this->particles_->template registerStateVariable<Vec2d>(
          "FirstOrderCorrectionVectorA2")),
      A3_(this->particles_->template registerStateVariable<Vec2d>(
          "FirstOrderCorrectionVectorA3")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      B_(this->particles_->template getVariableDataByName<Matd>(
          "LinearGradientCorrectionMatrix"))
    {  
        particles_->addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA1");
        particles_->addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA2");
        particles_->addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA3");
    }
          
    virtual ~NonisotropicKernelCorrectionMatrixAC(){};

  protected:
    Real *phi_;
    Vec2d *A1_;
    Vec2d *A2_;
    Vec2d *A3_;
    Real *Vol_;
    Mat2d *B_;

};


template <>
class NonisotropicKernelCorrectionMatrixAC<Inner<>>
    : public NonisotropicKernelCorrectionMatrixAC<DataDelegateInner>
{
  public:
    explicit NonisotropicKernelCorrectionMatrixAC(BaseInnerRelation &inner_relation)
        : NonisotropicKernelCorrectionMatrixAC<DataDelegateInner>(inner_relation){};
    template <typename BodyRelationType, typename FirstArg>
    explicit NonisotropicKernelCorrectionMatrixAC(InteractArgs<BodyRelationType, FirstArg> parameters)
        : NonisotropicKernelCorrectionMatrixAC(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~NonisotropicKernelCorrectionMatrixAC(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ik
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
           //Vec2d gradW_ikV_k = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
       
            Vecd r_ji = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];

            A1_[index_i] += r_ji[0] * r_ji[0] * (B_[index_i].transpose() * gradW_ij);
            A2_[index_i] += r_ji[1] * r_ji[1] * (B_[index_i].transpose() * gradW_ij);
            A3_[index_i] += r_ji[0] * r_ji[1] * (B_[index_i].transpose() * gradW_ij);
         
        }
       
    };
    void update(size_t index_i, Real dt = 0.0){};
};
  
template <>
class NonisotropicKernelCorrectionMatrixAC<Contact<>>
    : public NonisotropicKernelCorrectionMatrixAC<DataDelegateContact>
{
  public:
    explicit NonisotropicKernelCorrectionMatrixAC(BaseContactRelation &contact_relation)
    : NonisotropicKernelCorrectionMatrixAC<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_mass_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Mass"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };

    virtual ~NonisotropicKernelCorrectionMatrixAC(){};

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
                Vecd r_ji = -contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
                
                A1_[index_i] += r_ji[0] * r_ji[0] * (B_[index_i].transpose() * gradW_ij);
                A2_[index_i] += r_ji[1] * r_ji[1] * (B_[index_i].transpose() * gradW_ij);
                A3_[index_i] += r_ji[0] * r_ji[1] * (B_[index_i].transpose() * gradW_ij);
            }
        }

    };

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<Real *> contact_mass_;
};

using NonisotropicKernelCorrectionMatrixACComplex = ComplexInteraction<NonisotropicKernelCorrectionMatrixAC<Inner<>, Contact<>>>;


/*

Class nonisotropic diffusion relaxation

*/

template <typename... InteractionTypes>
class NonisotropicDiffusionRelaxation;

template <class DataDelegationType>
class NonisotropicDiffusionRelaxation<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit NonisotropicDiffusionRelaxation(BaseRelationType &base_relation)
      : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
      pos_(particles_->getVariableDataByName<Vecd>("Position")),
      B_(this->particles_->template getVariableDataByName<Matd>(
          "LinearGradientCorrectionMatrix")),
      phi_(particles_->registerStateVariable<Real>("Phi")),
      A1_(this->particles_->template getVariableDataByName<Vec2d>(
          "FirstOrderCorrectionVectorA2")),
      A2_(this->particles_->template getVariableDataByName<Vec2d>(
          "FirstOrderCorrectionVectorA2")),
      A3_(this->particles_->template getVariableDataByName<Vec2d>(
          "FirstOrderCorrectionVectorA3")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")) 
    {  
        
       SC_ = particles_->registerStateVariable<Mat3d>( "FirstOrderCorrectionMatrixSC", [&](size_t i) -> Mat3d { return Eps * Mat3d::Identity(); });
        
        
       E_= particles_->registerStateVariable<Vec2d>("FirstOrderCorrectionVectorE", [&](size_t i) -> Vec2d { return Eps * Vec2d::Identity(); });
    
       G_ =particles_->registerStateVariable<Vec3d>( "FirstOrderCorrectionVectorG", [&](size_t i) -> Vec3d { return Eps * Vec3d::Identity(); });
     Laplacian_=  particles_->registerStateVariable<Vec3d>( "Laplacian", [&](size_t i) -> Vec3d { return Vec3d::Zero(); });

     Laplacian_x=   particles_->registerStateVariable<Real>("Laplacian_x", [&](size_t i) -> Real { return Real(0.0); });
       Laplacian_y= particles_->registerStateVariable<Real>("Laplacian_y", [&](size_t i) -> Real { return Real(0.0); });
      Laplacian_xy= particles_->registerStateVariable<Real>( "Laplacian_xy", [&](size_t i) -> Real { return Real(0.0); });
		diffusion_dt_=particles_->registerStateVariable<Real>( "diffusion_dt", [&](size_t i) -> Real { return Real(0.0); });
		

        particles_->addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA1");
        particles_->addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA2");
        particles_->addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA3");
    }
          
    virtual ~NonisotropicDiffusionRelaxation(){};

  protected:
    Vec2d *pos_;
    Mat2d *B_;
    Real  *phi_;

    Vec2d *A1_, *A2_, *A3_;
    Real  *Vol_;

    Mat3d *SC_;
    Vec2d *E_;
    Vec3d *G_;
    Vec3d *Laplacian_;

    Real *Laplacian_x, *Laplacian_y, *Laplacian_xy, *diffusion_dt_;

    Real diffusion_coeff_;

};




template <>
class NonisotropicDiffusionRelaxation<Inner<>>
    : public NonisotropicDiffusionRelaxation<DataDelegateInner>
{
  public:
    explicit NonisotropicDiffusionRelaxation(BaseInnerRelation &inner_relation)
        : NonisotropicDiffusionRelaxation<DataDelegateInner>(inner_relation){};
    template <typename BodyRelationType, typename FirstArg>
    explicit NonisotropicDiffusionRelaxation(InteractArgs<BodyRelationType, FirstArg> parameters)
        : NonisotropicDiffusionRelaxation(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~NonisotropicDiffusionRelaxation(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        Vec2d E_rate = Vec2d::Zero();

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ik
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
           //Vec2d gradW_ikV_k = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
    
           E_rate += (phi_[index_j] - phi_[index_i]) * (B_[index_i].transpose() * gradW_ij);      
        }
         E_[index_i] = E_rate;
         
         Vec3d G_rate = Vec3d::Zero();
         Mat3d SC_rate = Mat3d::Zero();
         Real H_rate = 1.0;

        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ik
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd r_ji = -inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];

            Vec2d gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            Vec3d S_ = Vec3d(r_ji[0] * r_ji[0], r_ji[1] * r_ji[1], r_ji[0] * r_ji[1]);
            H_rate = r_ji.dot(B_[index_i].transpose() * gradW_ijV_j) / pow(r_ji.norm(), 4.0);
			 
            Real FF_ = 2.0 * (phi_[index_j] - phi_[index_i] - r_ji.dot(E_[index_i]));
            G_rate += S_ *H_rate * FF_;
            
             //TO DO
            Vec3d C_ = Vec3d::Zero();
            C_[0] = (r_ji[0] * r_ji[0]- r_ji.dot(A1_[index_i]));
            C_[1] = (r_ji[1] * r_ji[1]- r_ji.dot(A2_[index_i]));
            C_[2] = (r_ji[0] * r_ji[1]- r_ji.dot(A3_[index_i]));
            SC_rate += S_ *H_rate* C_.transpose();   

        }
        
        G_[index_i] = G_rate;
        SC_[index_i] = SC_rate;
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        phi_[index_i] += dt * (Laplacian_[index_i][0] + Laplacian_[index_i][1]);
    };

};




template <>
class NonisotropicDiffusionRelaxation<Contact<>>
    : public NonisotropicDiffusionRelaxation<DataDelegateContact>
{
  public:
    explicit NonisotropicDiffusionRelaxation(BaseContactRelation &contact_relation)
    : NonisotropicDiffusionRelaxation<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_mass_.push_back(contact_particles_[k]->getVariableDataByName<Real>("Mass"));
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };

    virtual ~NonisotropicDiffusionRelaxation(){};

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Mat3d SC_rate_contact = Mat3d::Zero();
        Vec3d G_rate_contact = Vec3d::Zero();
        Real H_rate_contact = 1.0;
        for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Real *Vol_k = contact_Vol_[k];
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                 size_t index_j = contact_neighborhood.j_[n];
                 Vecd r_ji = -contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
                 Vecd gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];

                Vec3d S_ = Vec3d(r_ji[0] * r_ji[0], r_ji[1] * r_ji[1], r_ji[0] * r_ji[1]);
                Real FF_ = 2.0 * (0.0 - r_ji.dot(E_[index_i])); ///here when it is periodic boundary condition, should notice the 0.0
                H_rate_contact = r_ji.dot(B_[index_i].transpose() * gradW_ij) / pow(r_ji.norm(), 4.0);
		 
                Vec3d C_ = Vec3d::Zero();
                C_[0] = (r_ji[0] * r_ji[0]- r_ji.dot(A1_[index_i]));
                C_[1] = (r_ji[1] * r_ji[1]- r_ji.dot(A2_[index_i]));
                C_[2] = (r_ji[0] * r_ji[1]- r_ji.dot(A3_[index_i]));

                SC_rate_contact += S_ * H_rate_contact * C_.transpose();
                G_rate_contact += S_ * H_rate_contact * FF_;
            
            }
            SC_[index_i] += SC_rate_contact;
            G_[index_i] += G_rate_contact;
        }
         Laplacian_[index_i] = diffusion_coeff_ * SC_[index_i].inverse() * G_[index_i];

        Laplacian_x[index_i] = Laplacian_[index_i][0];
        Laplacian_y[index_i] = Laplacian_[index_i][1];
        Laplacian_xy[index_i] = Laplacian_[index_i][2];
		diffusion_dt_[index_i] = Laplacian_[index_i][0] + Laplacian_[index_i][1];

    };

  protected:
    StdVec<Real *> contact_Vol_;
    StdVec<Real *> contact_mass_;
};

 
using NonisotropicDiffusionRelaxationComplex = ComplexInteraction<NonisotropicDiffusionRelaxation<Inner<>, Contact<>>>;




}
#endif // ANISOTROPIC
