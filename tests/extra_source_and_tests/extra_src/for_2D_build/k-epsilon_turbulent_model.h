/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4													*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	k-epsilon_turbulent_model.h
 * @brief
 * @details
 * @author Xiangyu Hu
 */

#ifndef K_EPSILON_TURBULENT_MODEL_H
#define K_EPSILON_TURBULENT_MODEL_H

#include "sphinxsys.h"
#include <mutex>

namespace SPH
{
namespace fluid_dynamics
{
class BaseTurbuClosureCoeff
{
  public:
    explicit BaseTurbuClosureCoeff();
    virtual ~BaseTurbuClosureCoeff() {};

  protected:
    Real Karman_;
    Real turbu_const_E_;
    Real C_mu_, C_mu_25_, C_mu_75_;
    Real turbulent_intensity_;

    Real sigma_k_;

    Real C_l_, C_2_;
    Real sigma_E_;
    Real turbulent_length_ratio_for_epsilon_inlet_;

    Real start_time_laminar_;
    Real y_star_threshold_laminar_;
};
//=================================================================================================//
class WallFunction : public BaseTurbuClosureCoeff
{
  public:
    explicit WallFunction() {};
    virtual ~WallFunction() {};

    Real get_dimensionless_velocity(Real y_star, Real time);
    Real get_near_wall_velocity_gradient_magnitude(Real y_star, Real vel_fric_mag, Real denominator_log_law, Real dynamic_viscosity);
    Real get_distance_from_P_to_wall(Real y_p_constant);

    Real log_law_wall_function(Real y_star);
    Real laminar_law_wall_function(Real y_star);
    Real log_law_velocity_gradient(Real vel_fric_mag, Real denominator_log_law);
    Real laminar_law_velocity_gradient(Real vel_fric_mag, Real dynamic_viscosity);
};
//=================================================================================================//
template <typename... InteractionTypes>
class GetVelocityGradient;

template <class DataDelegationType>
class GetVelocityGradient<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit GetVelocityGradient(BaseRelationType &base_relation);
    virtual ~GetVelocityGradient() {};

  protected:
    Real *Vol_;
    Vecd *vel_, *pos_;
    int *is_near_wall_P1_;
    int *is_near_wall_P2_;

    Matd *velocity_gradient_;

    Matd *velocity_gradient_wall;
};

template <>
class GetVelocityGradient<Inner<>> : public GetVelocityGradient<DataDelegateInner>
{
  public:
    explicit GetVelocityGradient(BaseInnerRelation &inner_relation, Real weight_sub);
    virtual ~GetVelocityGradient() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Matd *velocity_gradient_;
    Matd *B_;
    Matd *turbu_B_;
    Real weight_sub_nearwall_;
};
using GetVelocityGradientInner = GetVelocityGradient<Inner<>>;

template <>
class GetVelocityGradient<Contact<Wall>> : public InteractionWithWall<GetVelocityGradient>
{
  public:
    explicit GetVelocityGradient(BaseContactRelation &contact_relation);
    virtual ~GetVelocityGradient() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Matd *velocity_gradient_;
};

using GetVelocityGradientComplex = ComplexInteraction<GetVelocityGradient<Inner<>, Contact<Wall>>>;
//=================================================================================================//
template <typename... T>
class BaseTurbulentModel;

template <class DataDelegationType>
class BaseTurbulentModel<Base, DataDelegationType>
    : public LocalDynamics, public DataDelegationType, public BaseTurbuClosureCoeff
{
  public:
    template <class BaseRelationType>
    explicit BaseTurbulentModel(BaseRelationType &base_relation);
    virtual ~BaseTurbulentModel() {};

  protected:
    Matd *turbu_strain_rate_;
    Viscosity &viscosity_;
    Real mu_, smoothing_length_, particle_spacing_min_;
    Real *rho_, *Vol_;
    Vecd *vel_;
    int dimension_;
};
//=================================================================================================//
/**
 * @class K_TurbulentModelInner
 * @brief  K_TurbulentModelInner
 */
class K_TurbulentModelInner : public BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit K_TurbulentModelInner(BaseInnerRelation &inner_relation, const StdVec<Real> &initial_values, int is_extr_visc_dissipa, bool is_STL);
    virtual ~K_TurbulentModelInner() {};

    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *dk_dt_;
    Real *dk_dt_without_dissipation_;
    Real *k_production_;

    int *is_near_wall_P1_;
    Matd *velocity_gradient_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real *turbu_mu_;
    Matd *turbu_strain_rate_;
    int *is_extra_viscous_dissipation_;
    bool is_STL_;

    int *turbu_indicator_;
    Real *k_diffusion_, *vel_x_;
};
//=================================================================================================//
/**
 * @class E_TurbulentModelInner
 * @brief  E_TurbulentModelInner
 */
class E_TurbulentModelInner : public BaseTurbulentModel<Base, DataDelegateInner>
{
  public:
    explicit E_TurbulentModelInner(BaseInnerRelation &inner_relation, bool is_STL);
    virtual ~E_TurbulentModelInner() {};

    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *depsilon_dt_;
    Real *depsilon_dt_without_dissipation_;
    Real *ep_production;
    Real *ep_dissipation_;
    Real *ep_diffusion_;

    Real *turbu_mu_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real *k_production_;
    int *is_near_wall_P1_;
    bool is_STL_;
};
//=================================================================================================//

template <typename... InteractionTypes>
class TKEnergyForce;

template <class DataDelegationType>
class TKEnergyForce<Base, DataDelegationType>
    : public BaseTurbulentModel<Base, DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit TKEnergyForce(BaseRelationType &base_relation);
    virtual ~TKEnergyForce() {};

  protected:
    Vecd *force_;
    Real *mass_;
    int *indicator_;
    Vecd *pos_;
    Real *turbu_k_;
    Vecd *test_k_grad_rslt_;
};

template <>
class TKEnergyForce<Inner<>> : public TKEnergyForce<Base, DataDelegateInner>
{
  public:
    explicit TKEnergyForce(BaseInnerRelation &inner_relation);
    virtual ~TKEnergyForce() {};
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
    virtual ~TKEnergyForce() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *test_k_grad_rslt_;
    Matd *B_;
};

template <class InnerInteractionType, class... ContactInteractionTypes>
using BaseTKEnergyForceComplex = ComplexInteraction<TKEnergyForce<InnerInteractionType, ContactInteractionTypes...>>;

using TKEnergyForceComplex = BaseTKEnergyForceComplex<Inner<>, Contact<>>;
//=================================================================================================//
template <typename... InteractionTypes>
class TurbuViscousForce;

template <class DataDelegationType>
class TurbuViscousForce<DataDelegationType> : public ViscousForce<DataDelegationType>
{
  public:
    template <class BaseRelationType>
    explicit TurbuViscousForce(BaseRelationType &base_relation);
    virtual ~TurbuViscousForce() {};

  protected:
    Real *turbu_k_;
    Real *turbu_mu_;
    Real *wall_Y_plus_;
    Real *wall_Y_star_;
    Vecd *velo_friction_;
    Real *y_p_;
    int *is_near_wall_P2_;
    Viscosity viscosity_;
    Real molecular_viscosity_;
    Real c0_;
};

template <>
class TurbuViscousForce<Inner<>> : public TurbuViscousForce<DataDelegateInner>
{
  public:
    explicit TurbuViscousForce(BaseInnerRelation &inner_relation);
    virtual ~TurbuViscousForce() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    int *turbu_indicator_;
    int *is_extra_viscous_dissipation_;
    Matd *B_;
};

using BaseTurbuViscousForceWithWall = InteractionWithWall<TurbuViscousForce>;
template <>
class TurbuViscousForce<Contact<Wall>> : public BaseTurbuViscousForceWithWall, public WallFunction
{
  public:
    explicit TurbuViscousForce(BaseContactRelation &wall_contact_relation);
    virtual ~TurbuViscousForce() {};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real wall_particle_spacing_;
    Matd *B_;
    Real *physical_time_;
};

using TurbulentViscousForceWithWall = ComplexInteraction<TurbuViscousForce<Inner<>, Contact<Wall>>>;
//=================================================================================================//
class TurbulentEddyViscosity : public LocalDynamics, public BaseTurbuClosureCoeff
{
  public:
    explicit TurbulentEddyViscosity(SPHBody &sph_body);
    virtual ~TurbulentEddyViscosity() {};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *rho_;
    Real *turbu_mu_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real *wall_Y_plus_, *wall_Y_star_;
    Viscosity &viscosity_;
    Real mu_;
};
//=================================================================================================//
/**
 * @class TurbulentAdvectionTimeStepSize
 * @brief Computing the turbulent advection time step size
 */
class TurbulentAdvectionTimeStepSize : public LocalDynamicsReduce<ReduceMax>
{
  public:
    explicit TurbulentAdvectionTimeStepSize(SPHBody &sph_body, Real U_max, Real advectionCFL = 0.25);
    virtual ~TurbulentAdvectionTimeStepSize() {};
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
//=================================================================================================//
/**
 * @class   InflowTurbulentCondition
 * @brief   Inflow boundary condition which imposes directly to a given velocity profile.
 *          TargetVelocity gives the velocity profile along the inflow direction,
 *          i.e. x direction in local frame.
 */
class InflowTurbulentCondition : public BaseFlowBoundaryCondition, public BaseTurbuClosureCoeff
{
  public:
    explicit InflowTurbulentCondition(BodyPartByCell &body_part,
                                      Real CharacteristicLength, Real relaxation_rate, int type_turbu_inlet);
    virtual ~InflowTurbulentCondition() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    int type_turbu_inlet_;
    Real relaxation_rate_;
    Real CharacteristicLength_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
    Real TurbulentLength_;

    virtual Real getTurbulentInflowK(Vecd &position, Vecd &velocity, Real &turbu_k);
    virtual Real getTurbulentInflowE(Vecd &position, Real &turbu_k, Real &turbu_E);
};
//=================================================================================================//
class JudgeIsNearWall : public LocalDynamics, public DataDelegateContact, public BaseTurbuClosureCoeff
{
  public:
    JudgeIsNearWall(BaseInnerRelation &inner_relation,
                    BaseContactRelation &contact_relation);
    virtual ~JudgeIsNearWall() {};
    inline void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real *distance_to_dummy_interface_;
    Real *distance_to_dummy_interface_up_average_;
    int *is_near_wall_P1_;
    int *is_near_wall_P2_;
    int *index_nearest_;
    Vecd *e_nearest_tau_, *e_nearest_normal_;

    Vecd *pos_;
    int dimension_;
    Real fluid_particle_spacing_, wall_particle_spacing_;
    Vecd *distance_from_wall_;
    StdVec<Real *> contact_Vol_;
    StdVec<Vecd *> contact_n_;
};
//=================================================================================================//
class StandardWallFunctionCorrection : public LocalDynamics, public DataDelegateContact, public WallFunction
{
  public:
    StandardWallFunctionCorrection(BaseInnerRelation &inner_relation,
                                   BaseContactRelation &contact_relation, Real y_p_constant);
    virtual ~StandardWallFunctionCorrection() {};
    inline void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Real *y_p_;
    Real *wall_Y_plus_, *wall_Y_star_;
    Real *velo_tan_;
    Vecd *velo_friction_;

    Vecd *vel_;
    Real *rho_;
    Viscosity &viscosity_;
    Real molecular_viscosity_;
    Real *turbu_k_;
    Real *turbu_epsilon_;
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
};
//=================================================================================================//
class ConstrainNormalVelocityInRegionP : public LocalDynamics
{
  public:
    explicit ConstrainNormalVelocityInRegionP(SPHBody &sph_body);
    virtual ~ConstrainNormalVelocityInRegionP() {};

    void update(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *vel_;
    int *is_near_wall_P1_;
    Vecd *e_nearest_normal_;
};
//=================================================================================================//
template <typename... InteractionTypes>
class TurbulentLinearGradientCorrectionMatrix;

template <class DataDelegationType>
class TurbulentLinearGradientCorrectionMatrix<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit TurbulentLinearGradientCorrectionMatrix(BaseRelationType &base_relation);
    virtual ~TurbulentLinearGradientCorrectionMatrix() {};

  protected:
    Real *Vol_;
    Matd *turbu_B_;
    Matd *B_;
};

template <>
class TurbulentLinearGradientCorrectionMatrix<Inner<>>
    : public TurbulentLinearGradientCorrectionMatrix<DataDelegateInner>
{
    Real turbu_alpha_;

  public:
    explicit TurbulentLinearGradientCorrectionMatrix(BaseInnerRelation &inner_relation, Real alpha = Real(0))
        : TurbulentLinearGradientCorrectionMatrix<DataDelegateInner>(inner_relation), turbu_alpha_(alpha) {};
    template <typename BodyRelationType, typename FirstArg>
    explicit TurbulentLinearGradientCorrectionMatrix(DynamicsArgs<BodyRelationType, FirstArg> parameters)
        : TurbulentLinearGradientCorrectionMatrix(parameters.identifier_, std::get<0>(parameters.others_)){};
    virtual ~TurbulentLinearGradientCorrectionMatrix() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);
};
using TurbulentLinearGradientCorrectionMatrixInner = TurbulentLinearGradientCorrectionMatrix<Inner<>>;
//=================================================================================================//
template <class ParticleScope>
using TVC_Limited_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, TruncatedLinear, LinearGradientCorrection, ParticleScope>;
template <class ParticleScope>
using TVC_NoLimiter_withLinearGradientCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, NoLimiter, LinearGradientCorrection, ParticleScope>;
//=================================================================================================//
class ModifiedTruncatedLinear : public Limiter
{
    Real ref_, slope_;

  public:
    ModifiedTruncatedLinear(Real ref, Real slope = 1000.0)
        : Limiter(), ref_(ref), slope_(slope) {};
    Real operator()(Real measure)
    {
        Real measure_scale = measure * ref_;
        return SMIN(slope_ * measure_scale, Real(1));
    };
};

template <class ParticleScope>
using TVC_ModifiedLimited_RKGC_OBFCorrection =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, ModifiedTruncatedLinear, LinearGradientCorrectionWithBulkScope, ParticleScope>;
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
#endif // K_EPSILON_TURBULENT_MODEL_H