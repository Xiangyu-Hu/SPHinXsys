#ifndef K_OMEGA_TURBULENT_MODEL_H
#define K_OMEGA_TURBULENT_MODEL_H

#include "udf_common_turbulence_model.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
namespace udf
{
class kOmega_BaseTurbuClosureCoeff
{
  public:
    explicit kOmega_BaseTurbuClosureCoeff();
    virtual ~kOmega_BaseTurbuClosureCoeff(){};

  protected:
    Real std_kw_beta_star_;
    Real std_kw_sigma_star_;
    Real std_kw_alpha_;
    Real std_kw_sigma_;
    Real std_kw_beta_;
    Real std_kw_f_beta_;
    Real std_kw_beta_0_;
    Real std_kw_sigma_do_;
    Real std_kw_sigma_d_;
    Real std_kw_C_lim_;
    Real std_kw_beta_i_;
    Real std_kw_beta_star_25_;
    Real std_kw_beta_star_5_;

    Real turbulent_intensity_for_k_inlet_;
    Real turbulent_length_ratio_for_omega_inlet_;
    Real C_mu_75_for_omega_inlet_;
};

template <typename... InteractionTypes>
class kOmega_GetVelocityGradient;

template <class DataDelegationType>
class kOmega_GetVelocityGradient<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit kOmega_GetVelocityGradient(BaseRelationType &base_relation);
    virtual ~kOmega_GetVelocityGradient() {};

  protected:
    Real *Vol_;
    Vecd *vel_, *pos_;
    int *is_near_wall_P1_;
    int *is_near_wall_P2_;

    Matd *velocity_gradient_;
    Matd *velocity_gradient_wall;
};
template <>
class kOmega_GetVelocityGradient<Inner<>> : public kOmega_GetVelocityGradient<DataDelegateInner>
{
  public:
    explicit kOmega_GetVelocityGradient(BaseInnerRelation &inner_relation);
    virtual ~kOmega_GetVelocityGradient() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *velocity_gradient_;
    Matd *B_;
    Matd *turbu_B_;
    Matd *turbu_strain_rate_;
    Real *turbu_strain_rate_magnitude_;
    Real *turbu_strain_rate_traceless_magnitude_;
};
using kOmega_GetVelocityGradientInner = kOmega_GetVelocityGradient<Inner<>>;

template <>
class kOmega_GetVelocityGradient<Contact<Wall>> : public InteractionWithWall<kOmega_GetVelocityGradient>
{
  public:
    explicit kOmega_GetVelocityGradient(BaseContactRelation &contact_relation);
    virtual ~kOmega_GetVelocityGradient() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Matd *velocity_gradient_;
};
using kOmega_GetVelocityGradientComplex = ComplexInteraction<kOmega_GetVelocityGradient<Inner<>, Contact<Wall>>>;

template <typename... T>
class kOmega_BaseTurbulentModel;

template <class DataDelegationType>
class kOmega_BaseTurbulentModel<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType, public kOmega_BaseTurbuClosureCoeff
{
  public:
    template <class BaseRelationType>
    explicit kOmega_BaseTurbulentModel(BaseRelationType &base_relation);
    virtual ~kOmega_BaseTurbulentModel(){};

  protected:
    Viscosity &viscosity_;
    Real mu_, smoothing_length_, particle_spacing_min_;
    Real *rho_;
    Real *Vol_;
    Vecd *vel_;
    int dimension_;
};

class kOmega_kTransportEquationInner : public kOmega_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kOmega_kTransportEquationInner(BaseInnerRelation &inner_relation, const StdVec<Real> &initial_values, int is_extr_visc_dissipa, int is_blended = 0);
    virtual ~kOmega_kTransportEquationInner(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *dk_dt_;
    Real *dk_dt_without_dissipation_;
    Real *k_production_;
    Real *turbu_k_;
    Real *turbu_omega_;
    Real *turbu_mu_;
    int *is_extra_viscous_dissipation_;
    int *is_blended_;
    int *turbu_indicator_;
    Real *k_diffusion_;

    Matd *turbu_strain_rate_;
    int *is_near_wall_P1_;
    Matd *velocity_gradient_;
};

class kOmega_TKE_Diffusion : public kOmega_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kOmega_TKE_Diffusion(BaseInnerRelation &inner_relation);
    virtual ~kOmega_TKE_Diffusion(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *turbu_k_;
    Real *turbu_omega_;
    Real *k_diffusion_;
};

class kOmega_omegaTransportEquationInner : public kOmega_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kOmega_omegaTransportEquationInner(BaseInnerRelation &inner_relation);
    virtual ~kOmega_omegaTransportEquationInner(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *domega_dt_;
    Real *domega_dt_without_dissipation_;
    Real *omega_production_;
    Real *omega_dissipation_;
    Real *omega_diffusion_;
    Real *gradient_dot_k_omega_;
    Real *omega_cross_diffusion_;

    Real *turbu_mu_;
    Real *turbu_k_;
    Real *turbu_omega_;
    Real *k_production_;
    int *is_near_wall_P1_;
};

class kOmega_TSDR_Diffusion_and_Gradient_Dot_Inner : public kOmega_BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit kOmega_TSDR_Diffusion_and_Gradient_Dot_Inner(BaseInnerRelation &inner_relation);
    virtual ~kOmega_TSDR_Diffusion_and_Gradient_Dot_Inner(){};

    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *gradient_dot_k_omega_;
    Real *omega_diffusion_;
    Real *turbu_omega_;
    Real *turbu_k_;
    Matd *B_;
};

class kOmegaTurbulentEddyViscosity : public LocalDynamics,
                                     public kOmega_BaseTurbuClosureCoeff
{
  public:
    explicit kOmegaTurbulentEddyViscosity(SPHBody &sph_body);
    virtual ~kOmegaTurbulentEddyViscosity(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *rho_;
    Real *turbu_mu_;
    Real *turbu_k_;
    Real *turbu_omega_;
    Real *wall_Y_plus_;
    Real *wall_Y_star_;
    Real *turbu_strain_rate_traceless_magnitude_;
    Viscosity &viscosity_;
    Real mu_;
};

class kOmega_WallFunctionCorrection : public LocalDynamics,
                                      public DataDelegateContact,
                                      public WallFunction,
                                      public kOmega_BaseTurbuClosureCoeff
{
  public:
    kOmega_WallFunctionCorrection(BaseInnerRelation &inner_relation,
                                  BaseContactRelation &contact_relation);
    virtual ~kOmega_WallFunctionCorrection(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *y_p_;
    Real *wall_Y_plus_;
    Real *wall_Y_star_;
    Real *velo_tan_;
    Vecd *velo_friction_;
    Real *wall_shear_stress_;

    Vecd *vel_;
    Real *rho_;
    Viscosity &viscosity_;
    Real molecular_viscosity_;
    Real *turbu_k_;
    Real *turbu_omega_;
    Real *turbu_mu_;
    int *is_near_wall_P1_;
    int *is_near_wall_P2_;
    Matd *velocity_gradient_;
    Real *k_production_;
    Real *distance_to_dummy_interface_;
    Real *distance_to_dummy_interface_up_average_;
    int *index_nearest;
    Vecd *e_nearest_tau_;
    Vecd *e_nearest_normal_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_n_;
    Real *physical_time_;
    int *is_blended_;
    Real *turbu_strain_rate_magnitude_;
    Real *laminar_fraction_for_blend_;
};

class kOmega_InflowTurbulentCondition : public BaseFlowBoundaryCondition,
                                        public kOmega_BaseTurbuClosureCoeff
{
  public:
    explicit kOmega_InflowTurbulentCondition(BodyPartByCell &body_part,
                                             Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet_omega, int type_turbu_inlet_k);
    virtual ~kOmega_InflowTurbulentCondition(){};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int type_turbu_inlet_omega_, type_turbu_inlet_k_;
    Real relaxation_rate_;
    Real CharacteristicLength_;
    Real *turbu_k_;
    Real *turbu_omega_;
    Real TurbulentLength_;

    virtual Real getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k);
    virtual Real getTurbulentInflowTemporaryEpsilon(Vecd &position, Real &turbu_k, Real turbu_E);
    virtual Real getTurbulentInflowOmega(Vecd& position, Vecd& velocity, Real& turbu_omega);

    Real polyEval(const std::vector<Real>& a, Real x);
};

template <typename TargetTKE>
class kOmega_InflowTurbulentCondition_TKE : public BaseFlowBoundaryCondition,
    public kOmega_BaseTurbuClosureCoeff
{
public:
    explicit kOmega_InflowTurbulentCondition_TKE(OrientedBoxByCell& aligned_box_part, Real relaxation_rate)
        : BaseFlowBoundaryCondition(aligned_box_part),
        relaxation_rate_(relaxation_rate),
        turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
        oriented_box_(aligned_box_part.getOrientedBox()),
        transform_(oriented_box_.getTransform()), halfsize_(oriented_box_.HalfSize()),
        target_tke(*this),
        physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")) {}

    virtual ~kOmega_InflowTurbulentCondition_TKE() {};
    OrientedBox& getOrientedBox() { return oriented_box_; };
    void update(size_t index_i, Real dt = 0.0)
    {
        if (oriented_box_.checkContain(pos_[index_i]))
        {
            Real current_k = turbu_k_[index_i];
            Vecd frame_position = transform_.shiftBaseStationToFrame(pos_[index_i]);
            Vecd frame_velocity = transform_.xformBaseVecToFrame(vel_[index_i]);
            Real relaxed_frame_tke = target_tke(frame_position, frame_velocity, current_k, *physical_time_) * relaxation_rate_ + current_k * (1.0 - relaxation_rate_);
            turbu_k_[index_i] = relaxed_frame_tke;
        }
    };

protected:
    Real relaxation_rate_;
    Real* turbu_k_;
    OrientedBox& oriented_box_;
    Transform& transform_;
    Vecd halfsize_;
    TargetTKE target_tke;
    Real* physical_time_;
};

template <typename TargetTSDR>
class kOmega_InflowTurbulentCondition_TSDR : public BaseFlowBoundaryCondition,
    public kOmega_BaseTurbuClosureCoeff
{
public:
    explicit kOmega_InflowTurbulentCondition_TSDR(OrientedBoxByCell& aligned_box_part, Real relaxation_rate)
        : BaseFlowBoundaryCondition(aligned_box_part),
        relaxation_rate_(relaxation_rate),
        turbu_k_(particles_->getVariableDataByName<Real>("TurbulenceKineticEnergy")),
        turbu_omega_(particles_->getVariableDataByName<Real>("TurbulentSpecificDissipation")),
        oriented_box_(aligned_box_part.getOrientedBox()),
        transform_(oriented_box_.getTransform()), halfsize_(oriented_box_.HalfSize()),
        target_tsdr(*this),
        physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime")) {
    }

    virtual ~kOmega_InflowTurbulentCondition_TSDR() {};
    OrientedBox& getOrientedBox() { return oriented_box_; };
    void update(size_t index_i, Real dt = 0.0)
    {
        if (oriented_box_.checkContain(pos_[index_i]))
        {
            Real current_k = turbu_k_[index_i];
            Real current_omega = turbu_omega_[index_i];
            Vecd frame_position = transform_.shiftBaseStationToFrame(pos_[index_i]);
            Vecd frame_velocity = transform_.xformBaseVecToFrame(vel_[index_i]);
            Real relaxed_frame_tsdr = target_tsdr(frame_position, frame_velocity, current_omega, *physical_time_, current_k) * relaxation_rate_ + current_omega * (1.0 - relaxation_rate_);
            turbu_omega_[index_i] = relaxed_frame_tsdr;
        }
    };

protected:
    Real relaxation_rate_;
    Real* turbu_k_;
    Real* turbu_omega_;
    OrientedBox& oriented_box_;
    Transform& transform_;
    Vecd halfsize_;
    TargetTSDR target_tsdr;
    Real* physical_time_;
};

}
}
}
#endif