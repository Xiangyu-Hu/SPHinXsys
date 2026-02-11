#ifndef ELASTIC_SOLID_HPP
#define ELASTIC_SOLID_HPP

#include "elastic_solid.h"

namespace SPH
{
//=================================================================================================//
template <typename ExecutionPolicy>
NeoHookeanSolid::ConstituteKernel::ConstituteKernel(
    const ExecutionPolicy &ex_policy, NeoHookeanSolid &encloser)
    : rho0_(encloser.ReferenceDensity()), K0_(encloser.BulkModulus()), G0_(encloser.ShearModulus()),
      c0_(encloser.ReferenceSoundSpeed()), cs0_(encloser.ShearWaveSpeed()) {}
//=================================================================================================//
inline Matd NeoHookeanSolid::ConstituteKernel::ElasticLeftCauchy(
    const Matd &F, size_t index_i, Real dt)
{
    return F * F.transpose();
}
//=================================================================================================//
inline Real NeoHookeanSolid::ConstituteKernel::VolumetricKirchhoff(Real J)
{
    return 0.5 * K0_ * (J * J - 1);
}
//=================================================================================================//
inline Matd NeoHookeanSolid::ConstituteKernel::DeviatoricKirchhoff(const Matd &deviatoric_be)
{
    return G0_ * deviatoric_be;
}
//=================================================================================================//
template <typename ScalingType>
Matd NeoHookeanSolid::ConstituteKernel::NumericalDampingLeftCauchy(
    const Matd &deformation, const Matd &deformation_rate, const ScalingType &scaling, size_t index_i)
{
    Matd strain_rate = 0.5 * (deformation_rate * deformation.transpose() +
                              deformation * deformation_rate.transpose());
    Matd normal_rate = getDiagonal(strain_rate);
    return 0.5 * rho0_ * (cs0_ * (strain_rate - normal_rate) + c0_ * normal_rate) * scaling;
}
//=================================================================================================//
} // namespace SPH
#endif // ELASTIC_SOLID_HPP