#ifndef DENSITY_GRADIENT_H
#define DENSITY_GRADIENT_H
#include "particle_functors.h"
#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{

template <typename... InteractionTypes>
class DensityGradient;

template <class DataDelegationType>
class DensityGradient<DataDelegationType> : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit DensityGradient(BaseRelationType &base_relation);
    virtual ~DensityGradient(){};

  protected:
    Real *Vol_;
    Real *rho_;
    Vecd *rho_grad_;
};

template <class KernelCorrectionType>
class DensityGradient<Inner<KernelCorrectionType>> : public DensityGradient<DataDelegateInner>
{
  public:
    explicit DensityGradient(BaseInnerRelation &inner_relation);
    virtual ~DensityGradient(){};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

  protected:
    KernelCorrectionType kernel_correction_;
};

using DensityGradientInner = DensityGradient<Inner<LinearGradientCorrection>>;

template <>
class DensityGradient<Contact<Wall>> : public InteractionWithWall<DensityGradient>
{
  public:
    explicit DensityGradient(BaseContactRelation &wall_contact_relation);
    virtual ~DensityGradient(){};
    void interaction(size_t index_i, Real dt = 0.0);

  protected:
    Vecd *distance_from_wall_;
    Real *rho_;
    Vecd *rho_grad_;
    StdVec<Real *> wall_rho_;
};

template <class KernelCorrectionType>
using DensityGradientWithWall =
    ComplexInteraction<DensityGradient<Inner<KernelCorrectionType>, Contact<Wall>>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif // DENSITY_GRADIENT_H
