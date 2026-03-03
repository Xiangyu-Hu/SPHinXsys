#include "fluid_time_step.h"

#include "viscosity.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
AcousticTimeStep::AcousticTimeStep(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      fluid_(DynamicCast<Fluid>(this, particles_->getBaseMaterial())),
      rho_(particles_->getVariableDataByName<Real>("Density")),
      p_(particles_->getVariableDataByName<Real>("Pressure")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      acousticCFL_(acousticCFL) {}
//=================================================================================================//
Real AcousticTimeStep::reduce(size_t index_i, Real dt)
{
    return fluid_.getSoundSpeed(p_[index_i], rho_[index_i]) + vel_[index_i].norm();
}
//=================================================================================================//
Real AcousticTimeStep::outputResult(Real reduced_value)
{
    // since the particle does not change its configuration in pressure relaxation step
    // I chose a time-step size according to Eulerian method
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
SurfaceTensionTimeStep::SurfaceTensionTimeStep(SPHBody &sph_body, Real acousticCFL)
    : AcousticTimeStep(sph_body, acousticCFL),
      rho0_(sph_body_->getBaseMaterial().ReferenceDensity()),
      surface_tension_coeff_(*(particles_->getSingularVariableByName<Real>("SurfaceTensionCoef")->Data())) {}
//=================================================================================================//
Real SurfaceTensionTimeStep::outputResult(Real reduced_value)
{
    reduced_value = SMAX(reduced_value, (Real)sqrt(2.0 * Pi * surface_tension_coeff_ / (rho0_ * h_min_)));
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
AdvectionTimeStep::
    AdvectionTimeStep(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax>(sph_body),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      vel_(particles_->getVariableDataByName<Vecd>("Velocity")),
      force_(particles_->getVariableDataByName<Vecd>("Force")),
      force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
      h_min_(sph_body.getSPHAdaptation().MinimumSmoothingLength()),
      speed_ref_(U_ref), advectionCFL_(advectionCFL) {}
//=================================================================================================//
Real AdvectionTimeStep::reduce(size_t index_i, Real dt)
{
    Real acceleration_scale = 4.0 * h_min_ *
                              (force_[index_i] + force_prior_[index_i]).norm() / mass_[index_i];
    return SMAX(vel_[index_i].squaredNorm(), acceleration_scale);
}
//=================================================================================================//
Real AdvectionTimeStep::outputResult(Real reduced_value)
{
    Real speed_max = sqrt(reduced_value);
    return advectionCFL_ * h_min_ / (SMAX(speed_max, speed_ref_) + TinyReal);
}
//=================================================================================================//
AdvectionViscousTimeStep::AdvectionViscousTimeStep(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : AdvectionTimeStep(sph_body, U_ref, advectionCFL)
{
    Fluid &fluid = DynamicCast<Fluid>(this, particles_->getBaseMaterial());
    Viscosity &viscosity = DynamicCast<Viscosity>(this, particles_->getBaseMaterial());
    Real viscous_speed = viscosity.ReferenceViscosity() / fluid.ReferenceDensity() / h_min_;
    speed_ref_ = SMAX(viscous_speed, speed_ref_);
}
//=================================================================================================//
Real AdvectionViscousTimeStep::reduce(size_t index_i, Real dt)
{
    return AdvectionTimeStep::reduce(index_i, dt);
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
