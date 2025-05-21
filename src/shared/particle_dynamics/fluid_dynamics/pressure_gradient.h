// pressure_gradient.h
#ifndef PRESSURE_GRADIENT_H
#define PRESSURE_GRADIENT_H

#include "base_fluid_dynamics.h"

namespace SPH
{
namespace fluid_dynamics
{

template <typename... InteractionTypes>
class PressureGradient;

template <class DataDelegationType>
class PressureGradient<DataDelegationType> : public LocalDynamics, public DataDelegationType
{
public:
    template <class BaseRelationType>
    explicit PressureGradient(BaseRelationType &base_relation);
    virtual ~PressureGradient() {};

protected:
    Real *Vol_;
    Real *p_;
    Vecd *p_grad_;
    Real *rho_;
};

template <class KernelCorrectionType>
class PressureGradient<Inner<KernelCorrectionType>> : public PressureGradient<DataDelegateInner>
{
public:
    explicit PressureGradient(BaseInnerRelation &inner_relation);
    virtual ~PressureGradient() {};
    void interaction(size_t index_i, Real dt = 0.0);
    void update(size_t index_i, Real dt = 0.0);

protected:
    KernelCorrectionType kernel_correction_;
};

using PressureGradientInner = PressureGradient<Inner<NoKernelCorrection>>;

template <>
class PressureGradient<Contact<Wall>> : public InteractionWithWall<PressureGradient>
{
public:
    explicit PressureGradient(BaseContactRelation &wall_contact_relation);
    virtual ~PressureGradient() {};
    void interaction(size_t index_i, Real dt = 0.0);

protected:
    Vecd *distance_from_wall_;
    Real *p_;
    Real *rho_;
    Vecd *p_grad_;
};

template <class KernelCorrectionType>
using PressureGradientWithWall = ComplexInteraction<
    PressureGradient<Inner<KernelCorrectionType>, Contact<Wall>>>;

} // namespace fluid_dynamics
} // namespace SPH

#endif // PRESSURE_GRADIENT_H
