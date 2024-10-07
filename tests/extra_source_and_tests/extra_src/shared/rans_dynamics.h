#ifndef RANS_DYNAMICS_H
#define RANS_DYNAMICS_H
#include "eulerian_fluid_integration.h"
#include "all_particle_dynamics.h"
#include "riemann_solver.h"
#include "fluid_integration.hpp"

#include "base_fluid_dynamics.h"
#include "force_prior.h"
#include "viscous_dynamics.h"
namespace SPH
{
namespace fluid_dynamics
{ 
template <typename... InteractionTypes>
class TurbulentViscousForceInFVM;

    template <class DataDelegationType>
class TurbulentViscousForceInFVM<DataDelegationType>
    : public ForcePrior, public DataDelegationType
    {
      public:
        template <class BaseRelationType>
        explicit TurbulentViscousForceInFVM(BaseRelationType &base_relation);
        virtual ~TurbulentViscousForceInFVM(){};

      protected:
        Real *rho_, *mass_, *Vol_, *mu_t_, *walladjacentcellflag_;
        Vecd *vel_, *turbulent_viscous_force_;
        Real smoothing_length_;
};

    template <>
    class TurbulentViscousForceInFVM<Inner<>>
    : public TurbulentViscousForceInFVM<DataDelegateInner>
    {
      public:
        explicit TurbulentViscousForceInFVM(BaseInnerRelation &inner_relation);
        virtual ~TurbulentViscousForceInFVM(){};
        void interaction(size_t index_i, Real dt = 0.0);
    };
    using TurbulentViscousForceInner = TurbulentViscousForceInFVM<Inner<>>;
    //=================================================================================================//
    template <typename... InteractionTypes>
    class TkeGradientForceInFVM;

    template <class DataDelegationType>
    class TkeGradientForceInFVM<DataDelegationType>
        : public ForcePrior, public DataDelegationType
    {
      public:
        template <class BaseRelationType>
        explicit TkeGradientForceInFVM(BaseRelationType &base_relation);
        virtual ~TkeGradientForceInFVM(){};

      protected:
        Real *rho_, *mass_, *Vol_, *K_;
        Vecd *tke_gradient_force_;
    };

    template <>
    class TkeGradientForceInFVM<Inner<>>
        : public TkeGradientForceInFVM<DataDelegateInner>
    {
      public:
        explicit TkeGradientForceInFVM(BaseInnerRelation &inner_relation);
        virtual ~TkeGradientForceInFVM(){};
        void interaction(size_t index_i, Real dt = 0.0);
    };
    using TkeGradientForceInner = TkeGradientForceInFVM<Inner<>>;
//=================================================================================================//  
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // RANS_DYNAMICS_H
