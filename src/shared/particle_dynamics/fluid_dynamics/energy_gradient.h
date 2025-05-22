#ifndef ENERGY_GRADIENT_H
#define ENERGY_GRADIENT_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{

template <typename... InteractionTypes>
class EnergyGradient;

template <class DataDelegationType>
class EnergyGradient<DataDelegationType> : public LocalDynamics, public DataDelegationType
{
public:
    template <class BaseRelationType>
    explicit EnergyGradient(BaseRelationType &base_relation);
    virtual ~EnergyGradient() {};

protected:
    Real *Vol_;
    Real *energy_;
    Vecd *energy_grad_;
};

template <class KernelCorrectionType>
class EnergyGradient<Inner<KernelCorrectionType>> : public EnergyGradient<DataDelegateInner>
{
public:
    explicit EnergyGradient(BaseInnerRelation &inner_relation);
    virtual ~EnergyGradient() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

protected:
    KernelCorrectionType kernel_correction_;
};

using EnergyGradientInner = EnergyGradient<Inner<NoKernelCorrection>>;

template <>
class EnergyGradient<Contact<Wall>> : public InteractionWithWall<EnergyGradient>
{
public:
    explicit EnergyGradient(BaseContactRelation &wall_contact_relation);
    virtual ~EnergyGradient() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    Vecd *distance_from_wall_;
    Real *energy_;
    Vecd *energy_grad_;
};

template <class KernelCorrectionType>
using EnergyGradientWithWall = ComplexInteraction<EnergyGradient<Inner<KernelCorrectionType>, Contact<Wall>>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif // ENERGY_GRADIENT_H
