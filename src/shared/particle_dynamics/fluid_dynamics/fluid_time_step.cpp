#include "fluid_time_step.h"

#include "adaptation.h"
#include "base_body.hpp"
#include "viscosity.h"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
AcousticTimeStep::AcousticTimeStep(SPHBody &sph_body, Real acousticCFL)
    : LocalDynamicsReduce<ReduceMax<Real>>(sph_body),
      fluid_(DynamicCast<Fluid>(this, sph_body_->getMatterMaterial())),
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
AcousticTimeStepWithAcceleration::AcousticTimeStepWithAcceleration(SPHBody &sph_body, Real acousticCFL)
    : AcousticTimeStep(sph_body, acousticCFL),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      force_(particles_->getVariableDataByName<Vecd>("Force")),
      force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")) {}
//=================================================================================================//
Real AcousticTimeStepWithAcceleration::reduce(size_t index_i, Real dt)
{
    Real force_norm = (force_[index_i] + force_prior_[index_i]).norm();
    Real acceleration_scale = sqrt(4.0 * h_min_ * force_norm / mass_[index_i]);
    return SMAX(AcousticTimeStep::reduce(index_i, dt), acceleration_scale);
}
//=================================================================================================//
WallAccelerationTimeStep::WallAccelerationTimeStep(
    BaseContactRelation &wall_contact_relation, Real wallCFL)
    : LocalDynamicsReduce<ReduceMax<Real>>(wall_contact_relation.getSPHBody()),
      DataDelegateContact(wall_contact_relation),
      mass_(particles_->getVariableDataByName<Real>("Mass")),
      force_prior_(particles_->getVariableDataByName<Vecd>("ForcePrior")),
      h_min_(wall_contact_relation.getSPHBody().getSPHAdaptation().MinimumSmoothingLength()),
      wallCFL_(wallCFL)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Solid &wall_material =
            DynamicCast<Solid>(this, contact_bodies_[k]->getMatterMaterial());
        wall_acc_ave_.push_back(
            wall_material.AverageAcceleration(contact_particles_[k]));
    }
}
//=================================================================================================//
Real WallAccelerationTimeStep::reduce(size_t index_i, Real dt)
{
    Real face_acceleration = 0.0;
    const Vecd fluid_prior_acceleration = force_prior_[index_i] / mass_[index_i];
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Vecd *wall_acceleration = wall_acc_ave_[k];
        const Neighborhood &wall_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != wall_neighborhood.current_size_; ++n)
        {
            const size_t index_j = wall_neighborhood.j_[n];
            const Vecd &e_ij = wall_neighborhood.e_ij_[n];
            const Real candidate =
                SMAX(Real(0), (fluid_prior_acceleration - wall_acceleration[index_j]).dot(-e_ij));
            face_acceleration = SMAX(face_acceleration, candidate);
        }
    }
    return face_acceleration;
}
//=================================================================================================//
Real WallAccelerationTimeStep::outputResult(Real reduced_value)
{
    return wallCFL_ * sqrt(h_min_ / (reduced_value + TinyReal));
}
//=================================================================================================//
SurfaceTensionTimeStep::SurfaceTensionTimeStep(SPHBody &sph_body, Real acousticCFL)
    : AcousticTimeStep(sph_body, acousticCFL),
      rho0_(sph_body_->getMatterMaterial().ReferenceDensity()),
      surface_tension_coeff_(*(particles_->getSingleVariableByName<Real>("SurfaceTensionCoef")->Data())) {}
//=================================================================================================//
Real SurfaceTensionTimeStep::outputResult(Real reduced_value)
{
    reduced_value = SMAX(reduced_value, (Real)sqrt(2.0 * Pi * surface_tension_coeff_ / (rho0_ * h_min_)));
    return acousticCFL_ * h_min_ / (reduced_value + TinyReal);
}
//=================================================================================================//
AdvectionTimeStep::
    AdvectionTimeStep(SPHBody &sph_body, Real U_ref, Real advectionCFL)
    : LocalDynamicsReduce<ReduceMax<Real>>(sph_body),
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
    Fluid &fluid = DynamicCast<Fluid>(this, sph_body_->getMatterMaterial());
    Viscosity &viscosity = sph_body_->getMaterialProperty<Viscosity>();
    Real viscous_speed = viscosity.ReferenceViscosity() / fluid.ReferenceDensity() / h_min_;
    speed_ref_ = SMAX(viscous_speed, speed_ref_);
}
//=================================================================================================//
Real AdvectionViscousTimeStep::reduce(size_t index_i, Real dt)
{
    return AdvectionTimeStep::reduce(index_i, dt);
}
//=================================================================================================//
Real AdvectionTimeStepWithoutAcceleration::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].squaredNorm();
}
//=================================================================================================//
Real AdvectionViscousTimeStepWithoutAcceleration::reduce(size_t index_i, Real dt)
{
    return vel_[index_i].squaredNorm();
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
