#include "surface_tension.hpp"

namespace SPH
{
namespace fluid_dynamics
{
//=================================================================================================//
SurfaceTensionStress::
    SurfaceTensionStress(BaseContactRelation &contact_relation, Real surface_tension_coeff)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      color_gradient_(particles_->registerStateVariable<Vecd>("ColorGradient")),
      norm_direction_(particles_->registerStateVariable<Vecd>("NormDirection")),
      surface_tension_stress_(particles_->registerStateVariable<Matd>("SurfaceTensionStress")),
      surface_tension_coeff_(*(particles_->registerSingularVariable<Real>("SurfaceTensionCoef", surface_tension_coeff)->Data()))
{
    particles_->addEvolvingVariable<Vecd>("ColorGradient");
    particles_->addVariableToWrite<Vecd>("ColorGradient");
    particles_->addEvolvingVariable<Matd>("SurfaceTensionStress");
    particles_->addVariableToWrite<Matd>("SurfaceTensionStress");
    Real rho0 = sph_body_.getBaseMaterial().ReferenceDensity();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        Real rho0_k = contact_bodies_[k]->getBaseMaterial().ReferenceDensity();
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
        Real *Vol_k = contact_Vol_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            weighted_color_gradient -= 2 * contact_fraction_k *
                                       contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
        }
        color_gradient_[index_i] = weighted_color_gradient;
        norm_direction_[index_i] = weighted_color_gradient / (weighted_color_gradient.norm() + Eps);
        surface_tension_stress_[index_i] += surface_tension_coeff_ *
                                            (Matd::Identity() - norm_direction_[index_i] * norm_direction_[index_i].transpose()) *
                                            weighted_color_gradient.norm();
    }
}
//=================================================================================================//
SurfaceStressForce<Inner<>>::SurfaceStressForce(BaseInnerRelation &inner_relation, Real hourglass_control_coeff)
    : SurfaceStressForce<DataDelegateInner>(inner_relation), hourglass_control_coeff_(hourglass_control_coeff) {}
//=================================================================================================//
void SurfaceStressForce<Inner<>>::interaction(size_t index_i, Real dt)
{
    Vecd summation = ZeroData<Vecd>::value;
    const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    Matd tangential_direction_i = Matd::Identity() - norm_direction_[index_i] * norm_direction_[index_i].transpose();
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Vecd e_ij = inner_neighborhood.e_ij_[n];
        Matd tangential_direction_j = Matd::Identity() - norm_direction_[index_j] * norm_direction_[index_j].transpose();
        Vecd color_gradient_average = 0.5 * (color_gradient_[index_i] + color_gradient_[index_j]);
        Matd mismatch = Matd::Zero() - color_gradient_average * e_ij.transpose() * r_ij * (color_gradient_average * e_ij.transpose() * r_ij) / ((color_gradient_average * e_ij.transpose() * r_ij).norm() + Eps);
        Matd hourglass_correction = hourglass_control_coeff_ * surface_tension_coeff_ * 0.5 * (tangential_direction_i + tangential_direction_j) * mismatch / (r_ij + Eps);
        summation += mass_[index_i] * inner_neighborhood.dW_ij_[n] * Vol_[index_j] *
                     (surface_tension_stress_[index_i] + surface_tension_stress_[index_j] + hourglass_correction) * inner_neighborhood.e_ij_[n];
    }
    surface_tension_force_[index_i] = summation / rho_[index_i];
}
//=================================================================================================//
SurfaceStressForce<Contact<>>::SurfaceStressForce(BaseContactRelation &contact_relation, Real hourglass_control_coeff)
    : SurfaceStressForce<DataDelegateContact>(contact_relation), hourglass_control_coeff_(hourglass_control_coeff)
{
    Real rho0 = sph_body_.getBaseMaterial().ReferenceDensity();
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        Real rho0_k = contact_bodies_[k]->getBaseMaterial().ReferenceDensity();
        contact_fraction_.push_back(rho0 / (rho0 + rho0_k));
        contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_color_gradient_.push_back(
            contact_particles_[k]->getVariableDataByName<Vecd>("ColorGradient"));
        contact_surface_tension_stress_.push_back(
            contact_particles_[k]->getVariableDataByName<Matd>("SurfaceTensionStress"));
        contact_norm_direction_.push_back(
            contact_particles_[k]->getVariableDataByName<Vecd>("NormDirection"));
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
        Vecd *contact_norm_direction_k = contact_norm_direction_[k];
        const Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Vecd e_ij = contact_neighborhood.e_ij_[n];
            Vecd color_gradient_average = 0.5 * (color_gradient_[index_i] + contact_color_gradient_k[index_j]);
            Matd mismatch = Matd::Identity() - color_gradient_average * e_ij.transpose() * r_ij * (color_gradient_average * e_ij.transpose() * r_ij) / ((color_gradient_average * e_ij.transpose() * r_ij).norm() + Eps);
            Matd hourglass_correction = -4 * contact_fraction_k * (1 - contact_fraction_k) * hourglass_control_coeff_ * 0.5 * (norm_direction_[index_i] * norm_direction_[index_i].transpose() + contact_norm_direction_k[index_j] * contact_norm_direction_k[index_j].transpose()) * mismatch * surface_tension_coeff_ / r_ij;
            summation += mass_[index_i] * (2 * (Real(1) - contact_fraction_k) * surface_tension_stress_[index_i] + 2 * contact_fraction_k *             
                        contact_surface_tension_stress_k[index_j] + hourglass_correction) 
                        * contact_neighborhood.dW_ij_[n] * contact_neighborhood.e_ij_[n] * Vol_k[index_j];
        }
    }
    surface_tension_force_[index_i] += summation / rho_[index_i];
}
//=================================================================================================//
} // namespace fluid_dynamics
} // namespace SPH
