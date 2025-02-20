#include "surface_tension.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
SurfaceTensionStress::
    SurfaceTensionStress(BaseContactRelation &contact_relation, StdVec<Real> contact_surface_tension)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      color_gradient_(particles_->registerStateVariable<Vecd>("ColorGradient")),
      surface_tension_stress_(particles_->registerStateVariable<Matd>("SurfaceTensionStress")),
      contact_surface_tension_(contact_surface_tension)
{
    particles_->addEvolvingVariable<Vecd>("ColorGradient");
    particles_->addVariableToWrite<Vecd>("ColorGradient");
    particles_->addEvolvingVariable<Matd>("SurfaceTensionStress");
    particles_->addVariableToWrite<Matd>("SurfaceTensionStress");
    Real rho0 = getSPHBody().base_material_->ReferenceDensity();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_surface_tension_.push_back(contact_surface_tension[k]);
        Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
    }
}
//=================================================================================================//
void SurfaceTensionStress::interaction(size_t index_i, Real dt)
{
    color_gradient_[index_i] = ZeroData<Vecd>::value;
    surface_tension_stress_[index_i] = ZeroData<Matd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Vecd weighted_color_gradient = ZeroData<Vecd>::value;
        Real contact_fraction_k = contact_fraction_[k];
        Real surface_tension_k = contact_surface_tension_[k];
        Real *Vol_k = contact_Vol_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            weighted_color_gradient -= contact_fraction_k *
                                       contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
        }
        color_gradient_[index_i] = weighted_color_gradient;
        Real norm = weighted_color_gradient.norm();
        surface_tension_stress_[index_i] += surface_tension_k / (norm + Eps) *
                                            (norm * norm * Matd::Identity() -
                                             weighted_color_gradient * weighted_color_gradient.transpose());
    }
}
//=================================================================================================//
SurfaceStressForce<Inner<>>::SurfaceStressForce(BaseInnerRelation &inner_relation)
    : SurfaceStressForce<DataDelegateInner>(inner_relation) {}
//=================================================================================================//
void SurfaceStressForce<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd summation = ZeroData<Vecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        summation += mass_[index_i] * inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
                     (surface_tension_stress_[index_i] + surface_tension_stress_[index_j]) *
                     inner_neighborhood.e_ij_[n];
    }
    surface_tension_force_[index_i] = summation / rho_[index_i];
}
//=================================================================================================//
SurfaceStressForce<Contact<>>::SurfaceStressForce(BaseContactRelation &contact_relation)
    : SurfaceStressForce<DataDelegateContact>(contact_relation)
{
    Real rho0 = getSPHBody().base_material_->ReferenceDensity();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real rho0_k = contact_bodies_[k]->base_material_->ReferenceDensity();
        contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_color_gradient_.push_back(
            contact_particles_[k]->getVariableDataByName<Vecd>("ColorGradient"));
        contact_surface_tension_stress_.push_back(
            contact_particles_[k]->getVariableDataByName<Matd>("SurfaceTensionStress"));
    }
}
//=================================================================================================//
void SurfaceStressForce<Contact<>>::interaction(size_t index_i, Real dt)
{
    Vecd summation = ZeroData<Vecd>::value;
    for (size_t k = 0; k < contact_configuration_.size(); ++k)
    {
        Real contact_fraction_k = contact_fraction_[k];
        Real *Vol_k = contact_Vol_[k];
        Vecd *contact_color_gradient_k = contact_color_gradient_[k];
        Matd *contact_surface_tension_stress_k = contact_surface_tension_stress_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Real mismatch = 1.0 - 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]).dot(e_ij) * r_ij;
            summation += mass_[index_i] * contact_neighborhood.dW_ij_[n] * Vol_k[index_j] *
                         (-0.1 * mismatch * Matd::Identity() +
                          (Real(1) - contact_fraction_k) * surface_tension_stress_[index_i] +
                          contact_surface_tension_stress_k[index_j] * contact_fraction_k) *
                         contact_neighborhood.e_ij_[n];
        }
    }
    surface_tension_force_[index_i] += summation / rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
