#ifndef UDF_COMMON_TURBULENCE_MODEL_H
#define UDF_COMMON_TURBULENCE_MODEL_H

#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
namespace udf
{

class WallFunctionCoefficient
{
  public:
    explicit WallFunctionCoefficient();
    virtual ~WallFunctionCoefficient(){};

  protected:
    Real Karman_;
    Real turbu_const_E_, inv_turbu_E_;
    Real C_mu_wf_, C_mu_wf_25_, C_mu_wf_75_;
    Real start_time_laminar_;
    Real y_star_threshold_laminar_;
};

class WallFunction : public WallFunctionCoefficient
{
  public:
    explicit WallFunction(){};
    virtual ~WallFunction(){};

    Real get_dimensionless_velocity(Real y_star, Real time, Real u_star_previous, int is_blended);
    Real get_near_wall_velocity_gradient_magnitude(Real y_star, Real vel_fric_mag, Real denominator_log_law, Real dynamic_viscosity);
    Real get_distance_from_P_to_wall(Real y_p_constant);

    Real log_law_wall_function(Real y_star);
    Real laminar_law_wall_function(Real y_star);
    Real log_law_velocity_gradient(Real vel_fric_mag, Real denominator_log_law);
    Real laminar_law_velocity_gradient(Real vel_fric_mag, Real dynamic_viscosity);
    Real Spalding_wall_function(Real y_star, Real u_star_guess);
};

template <typename... InteractionTypes>
class kEpsilon_GetVelocityGradient;

template <class DataDelegationType>
class kEpsilon_GetVelocityGradient<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit kEpsilon_GetVelocityGradient(BaseRelationType &base_relation);
    virtual ~kEpsilon_GetVelocityGradient(){};

  protected:
    Real *Vol_;
    Vecd *vel_, *pos_;
    int *is_near_wall_P1_;
    int *is_near_wall_P2_;

    Matd *velocity_gradient_;

    Matd *velocity_gradient_wall;
};

template <>
class kEpsilon_GetVelocityGradient<Inner<>> : public kEpsilon_GetVelocityGradient<DataDelegateInner>
{
  public:
    explicit kEpsilon_GetVelocityGradient(BaseInnerRelation &inner_relation, Real weight_sub);
    virtual ~kEpsilon_GetVelocityGradient(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *velocity_gradient_;
    Matd *B_;
    Matd *turbu_B_;
    Real weight_sub_nearwall_;
};
using kEpsilon_GetVelocityGradientInner = kEpsilon_GetVelocityGradient<Inner<>>;

template <typename... InteractionTypes>
class TKEnergyForce;

template <class DataDelegationType>
class TKEnergyForce<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit TKEnergyForce(BaseRelationType &base_relation);
    virtual ~TKEnergyForce(){};

  protected:
    Vecd *force_;
    Real *mass_;
    int *indicator_;
    Vecd *pos_;
    Real *turbu_k_;
    Real *Vol_;
    Vecd *test_k_grad_rslt_;
};

template <>
class TKEnergyForce<Inner<>> : public TKEnergyForce<Base, DataDelegateInner>
{
  public:
    explicit TKEnergyForce(BaseInnerRelation &inner_relation);
    virtual ~TKEnergyForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *test_k_grad_rslt_;
    Matd *B_;
};

template <>
class TKEnergyForce<Contact<>> : public TKEnergyForce<Base, DataDelegateContact>
{
  public:
    explicit TKEnergyForce(BaseContactRelation &contact_relation);
    virtual ~TKEnergyForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *test_k_grad_rslt_;
    Matd *B_;
};

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseTKEnergyForceComplex = ComplexInteraction<TKEnergyForce<InnerInteractionType, ContactInteractionTypes...>>;

using TKEnergyForceComplex = BaseTKEnergyForceComplex<Inner<>, Contact<>>;

template <typename... InteractionTypes>
class TurbuViscousForce;

template <class DataDelegationType>
class TurbuViscousForce<DataDelegationType> : public ViscousForce<DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit TurbuViscousForce(BaseRelationType &base_relation);
    virtual ~TurbuViscousForce(){};

  protected:

    Real *turbu_k_;
    Real *turbu_mu_;
    Real *wall_Y_plus_;
    Real *wall_Y_star_;
    Vecd *velo_friction_;
    Real *y_p_;
    int *is_near_wall_P2_;
    Viscosity &viscosity_;
    Real molecular_viscosity_;
    Real c0_;
};

template <>
class TurbuViscousForce<Inner<>> : public TurbuViscousForce<DataDelegateInner>
{
  public:
    explicit TurbuViscousForce(BaseInnerRelation &inner_relation);
    virtual ~TurbuViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    int *turbu_indicator_;
    int *is_extra_viscous_dissipation_;
    Matd *B_;
    Matd *turbu_B_;
    int* is_near_wall_P1_;
    Matd* dUdn_P_sublayer_;
    Real* physical_time_;
};

using BaseTurbuViscousForceWithWall = InteractionWithWall<TurbuViscousForce>;
template <>
class TurbuViscousForce<Contact<Wall>> : public BaseTurbuViscousForceWithWall, public WallFunction
{
  public:
    explicit TurbuViscousForce(BaseContactRelation &wall_contact_relation);
    virtual ~TurbuViscousForce(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real wall_particle_spacing_;
    Matd *B_;
    Real *physical_time_;
    int *is_blended_;
    int* is_near_wall_P1_;
    Real* friction_velocity_from_sublayer_;
    Matd* turbu_B_;
    Matd* B_only_wall_;
};

using TurbulentViscousForceWithWall = ComplexInteraction<TurbuViscousForce<Inner<>, Contact<Wall>>>;

class TurbulentAdvectionTimeStepSize : public LocalDynamicsReduce<ReduceMax<Real>>
{
  public:
    explicit TurbulentAdvectionTimeStepSize(SPHBody &sph_body, Real U_max, Real advectionCFL = 0.25);
    virtual ~TurbulentAdvectionTimeStepSize(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;

  protected:
    Vecd *vel_;
    Real smoothing_length_min_;
    Real speed_ref_turbu_, advectionCFL_;
    Real *turbu_mu_;
    Fluid &fluid_;
    Viscosity &viscosity_;
};

class JudgeIsNearWall : public LocalDynamics, public DataDelegateContact
{
  public:
    JudgeIsNearWall(BaseInnerRelation &inner_relation,
                    BaseContactRelation &contact_relation, Real constant_y_p);
    virtual ~JudgeIsNearWall(){};
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *distance_to_dummy_interface_;
    Real *distance_to_dummy_interface_up_average_;
    int *is_near_wall_P1_;
    int *is_near_wall_P2_;
    int *index_nearest_;
    Vecd *e_nearest_tau_, *e_nearest_normal_;
    Real *y_p_;
    Real constant_y_p_;

    Vecd *pos_;
    Real fluid_particle_spacing_, wall_particle_spacing_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_n_;
};

class ConstrainNormalVelocityInRegionP : public LocalDynamics, public WallFunctionCoefficient
{
  public:
    explicit ConstrainNormalVelocityInRegionP(SPHBody &sph_body);
    virtual ~ConstrainNormalVelocityInRegionP(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *vel_;
    int *is_near_wall_P1_;
    Vecd *e_nearest_normal_;
    Real *wall_Y_star_;
};

template <typename... InteractionTypes>
class TurbulentLinearGradientCorrectionMatrix;

template <class DataDelegationType>
class TurbulentLinearGradientCorrectionMatrix<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit TurbulentLinearGradientCorrectionMatrix(BaseRelationType &base_relation);
    virtual ~TurbulentLinearGradientCorrectionMatrix(){};

  protected:
    Real *Vol_;
    Matd *turbu_B_;
    Matd *B_only_wall_;
};

template <>
class TurbulentLinearGradientCorrectionMatrix<Inner<>>
    : public TurbulentLinearGradientCorrectionMatrix<DataDelegateInner>
{
    Real turbu_alpha_;
    int *is_near_wall_P1_;

  public:
    explicit TurbulentLinearGradientCorrectionMatrix(BaseInnerRelation &inner_relation, Real alpha = Real(0))
        : TurbulentLinearGradientCorrectionMatrix<DataDelegateInner>(inner_relation), turbu_alpha_(alpha),
        is_near_wall_P1_(particles_->getVariableDataByName<int>("IsNearWallP1")) {};
    template <typename BodyRelationType, typename FirstArg>
    explicit TurbulentLinearGradientCorrectionMatrix(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : TurbulentLinearGradientCorrectionMatrix(parameters.body_relation_, std::get<0>(parameters.others_)){};
    virtual ~TurbulentLinearGradientCorrectionMatrix(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using TurbulentLinearGradientCorrectionMatrixInner = TurbulentLinearGradientCorrectionMatrix<Inner<>>;

template <>
class TurbulentLinearGradientCorrectionMatrix<Contact<>>
    : public TurbulentLinearGradientCorrectionMatrix<DataDelegateContact>
{
public:
    explicit TurbulentLinearGradientCorrectionMatrix(BaseContactRelation& contact_relation);
    virtual ~TurbulentLinearGradientCorrectionMatrix() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    StdVec<Real*> contact_Vol_;
    StdVec<Real*> contact_mass_;
};

using TurbulentLinearGradientCorrectionMatrixComplex = ComplexInteraction<TurbulentLinearGradientCorrectionMatrix<Inner<>, Contact<>>>;

class GetLimiterOfTransportVelocityCorrection : public LocalDynamics
{
  public:
    explicit GetLimiterOfTransportVelocityCorrection(SPHBody &sph_body, Real slope);
    virtual ~GetLimiterOfTransportVelocityCorrection(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    const Real h_ref_;
    Vecd *zero_gradient_residue_;
    Real slope_;
    Real *limiter_tvc_;
};

template <class ParticleScope>
using TVC_Limited_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, TruncatedLinear, LinearGradientCorrection, ParticleScope>;
template <class ParticleScope>
using TVC_NoLimiter_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, NoLimiter, LinearGradientCorrection, ParticleScope>;

class ModifiedTruncatedLinear : public Limiter
{
    Real slope_;

  public:
    ModifiedTruncatedLinear(Real slope = 1000.0)
        : Limiter(), slope_(slope) {};
    Real operator()(Real dimensionless_measure)
    {
        return SMIN(slope_ * dimensionless_measure, Real(1));
    };
};
template <class ParticleScope>
using TVC_ModifiedLimited_NoRKGC =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, ModifiedTruncatedLinear, NoKernelCorrection, ParticleScope>;

template <class ParticleScope>
using TVC_ModifiedLimited_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, ModifiedTruncatedLinear, LinearGradientCorrection, ParticleScope>;

template <class ParticleScope>
using TVC_ModifiedLimited_RKGC_OBFCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, ModifiedTruncatedLinear, LinearGradientCorrectionWithBulkScope, ParticleScope>;

template <class ParticleScope>
using TVC_NotLimited_RKGC_OBFCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, NoLimiter, LinearGradientCorrectionWithBulkScope, ParticleScope>;

template <class ParticleScope>
using TVC_ModifiedLimited_withoutLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, ModifiedTruncatedLinear, NoKernelCorrection, ParticleScope>;

class NonDimensionalisePressure : public LocalDynamics
{
  public:
    explicit NonDimensionalisePressure(SPHBody &sph_body);
    virtual ~NonDimensionalisePressure(){};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *rho_;
    Real *p_;
    Real *p_dimensionless_;
};

template <typename... InteractionTypes>
class TurbulentIntegration2ndHalf;

template <class RiemannSolverType>
class TurbulentIntegration2ndHalf<Contact<Wall>, RiemannSolverType>
    : public BaseIntegrationWithWall
{
  public:
    explicit TurbulentIntegration2ndHalf(BaseContactRelation &wall_contact_relation);
    virtual ~TurbulentIntegration2ndHalf(){};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    RiemannSolverType riemann_solver_;
};
using Integration2ndHalfOnlyWallAcousticRiemannAdjusted = TurbulentIntegration2ndHalf<Contact<Wall>, DissipativeRiemannSolver>;

}
}
}
#endif