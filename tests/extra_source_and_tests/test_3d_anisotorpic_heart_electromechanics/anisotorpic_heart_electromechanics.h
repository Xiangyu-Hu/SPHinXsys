/**
 * @file anisotorpic_heart_electromechanics.h
 * @brief This is for anisotropic voltage diffusion relaxation.
 * @author Xiaojing Tang  
 */

 #ifndef ANISOTROPICVOLTAGE
 #define ANISOTROPICVOLTAGE
 
 #include "anisotropic_diffusion_relaxation.h"
 namespace SPH
 {
 
/*
Class Anisotropic diffusion relaxation
*/
template <typename... InteractionTypes>
class AnisotropicVoltageDiffusionRelaxation;

template <class DataDelegationType>
class AnisotropicVoltageDiffusionRelaxation<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:

    template <class BaseRelationType, class AnisotropicDiffusionSolidType = AnisotropicDiffusionSolid >
    explicit AnisotropicVoltageDiffusionRelaxation(BaseRelationType &base_relation)
      : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation), 
      B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
      species_(particles_->registerStateVariable<Real>("Voltage")),
      Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
      kernel_correction1_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector1")),
      kernel_correction2_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector2")),
      kernel_correction3_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector3")),
      kernel_correction4_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector4")),
      kernel_correction5_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector5")),
      kernel_correction6_(this->particles_->template getVariableDataByName<Vecd>("FirstOrderCorrectionVector6")),
      pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
      diffusion_coeff_tensor_(this->particles_->template getVariableDataByName<Matd>("DiffusionTensor")) 
      { 
        
        species_correction_ = particles_->registerStateVariable<Vecd>("SpeciesCorrection ", [&](size_t i) -> Vecd { return Eps * Vecd::Identity(); });
        total_left_3d = particles_->registerStateVariable<Mat6d>( "TotalLeft3D", [&](size_t i) -> Mat6d { return Eps * Mat6d::Identity(); });
        total_right_3d = particles_->registerStateVariable<Vec6d>( "TotalRight3D", [&](size_t i) -> Vec6d { return Eps * Vec6d::Identity(); });
        Laplacian_3d =  particles_->registerStateVariable<Vec6d>( "Laplacian3D", [&](size_t i) -> Vec6d { return Vec6d::Zero(); });
        diffusion_dt_=particles_->registerStateVariable<Real>( "diffusion_dt", [&](size_t i) -> Real { return Real(0.0); });

    }
          
    virtual ~AnisotropicVoltageDiffusionRelaxation(){};

  protected: 
    Matd *B_;
    Real  *species_,*Vol_;
    Vecd *kernel_correction1_, *kernel_correction2_, *kernel_correction3_;
    Vecd *kernel_correction4_, *kernel_correction5_, *kernel_correction6_;
    Vecd *species_correction_;

    Mat6d *total_left_3d;
    Vec6d *total_right_3d;
    Vec6d *Laplacian_3d;


    Real  *diffusion_dt_;

    Vecd *pos_;
    Matd  *diffusion_coeff_tensor_;

    
};

template <>
class AnisotropicVoltageDiffusionRelaxation<Inner<>>
    : public AnisotropicVoltageDiffusionRelaxation<DataDelegateInner>
{
  public:
    explicit AnisotropicVoltageDiffusionRelaxation(BaseInnerRelation &inner_relation)
        : AnisotropicVoltageDiffusionRelaxation<DataDelegateInner>(inner_relation){};

    template <typename BodyRelationType, typename FirstArg>
    explicit AnisotropicVoltageDiffusionRelaxation(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : AnisotropicVoltageDiffusionRelaxation(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~AnisotropicVoltageDiffusionRelaxation(){};

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
            Vec6d total_right_rate = Vec6d::Zero();
            Mat6d total_left_rate = Mat6d::Zero();
 
           for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)  
           {
              size_t index_j = inner_neighborhood.j_[n]; 
              Vec3d r_ji = - inner_neighborhood.r_ij_[n] *inner_neighborhood.e_ij_[n];
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
              total_right_rate += disp_quad_ * modified_func_ * right_;
            }
    
          total_right_3d[index_i] = total_right_rate;
          total_left_3d[index_i] = total_left_rate;
 
    };
 
    void update(size_t index_i, Real dt = 0.0)
    {
        species_[index_i] += dt * diffusion_dt_[index_i];
    };
   
    
};
 

template <>
class AnisotropicVoltageDiffusionRelaxation<Contact<>>
    : public AnisotropicVoltageDiffusionRelaxation<DataDelegateContact>
{
  public:
    explicit AnisotropicVoltageDiffusionRelaxation(BaseContactRelation &contact_relation)
    : AnisotropicVoltageDiffusionRelaxation<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k != contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
            contact_pos_.push_back(contact_particles_[k]->getVariableDataByName<Vecd>("Position"));  
        }
    };

    virtual ~AnisotropicVoltageDiffusionRelaxation(){};

    void interaction(size_t index_i, Real dt = 0.0)
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
              Vec3d r_ji = -contact_neighborhood.r_ij_[n] *contact_neighborhood.e_ij_[n];
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
          Laplacian_3d[index_i] =  total_left_3d[index_i].inverse() * total_right_3d[index_i];
    
          Mat3d Laplacian_transform = Mat3d { 
            { Laplacian_3d[index_i][0], 0.5 * Laplacian_3d[index_i][3],  0.5 * Laplacian_3d[index_i][5] },  
             {0.5 * Laplacian_3d[index_i][3],   Laplacian_3d[index_i][1], 0.5 * Laplacian_3d[index_i][4]},
              {0.5 * Laplacian_3d[index_i][5],  0.5 * Laplacian_3d[index_i][4],  Laplacian_3d[index_i][2]}
           }; 
         
         Mat3d Laplacian_transform_aniso = diffusion_coeff_tensor_[index_i]* Laplacian_transform *  diffusion_coeff_tensor_[index_i].transpose() ;
 
      
        diffusion_dt_[index_i] = Laplacian_transform_aniso(0,0) + Laplacian_transform_aniso(1,1)
                                 + Laplacian_transform_aniso(2,2);    
      

    };

  protected:
    StdVec<Real *> contact_Vol_; 
    StdVec<Vecd *> contact_pos_;
 
};

using AnisotropicVoltageDiffusionRelaxationComplex = ComplexInteraction<AnisotropicVoltageDiffusionRelaxation<Inner<>, Contact<>>>;

 
}
#endif // ANISOTROPICVOLTAGE
