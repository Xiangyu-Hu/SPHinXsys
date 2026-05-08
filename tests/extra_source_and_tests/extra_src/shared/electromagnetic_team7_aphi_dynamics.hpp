#ifndef ELECTROMAGNETIC_TEAM7_APHI_DYNAMICS_HPP
#define ELECTROMAGNETIC_TEAM7_APHI_DYNAMICS_HPP

#include "electromagnetic_team7_aphi_dynamics.h"
#include <cmath>

namespace SPH
{
namespace electromagnetics
{
//=================================================================================================//
namespace
{
inline Vec3d CurlNuBContribution(const Vec3d &nu_b_diff, const Vec3d &grad_w_ijV_j)
{
    return nu_b_diff.cross(grad_w_ijV_j);
}

inline Vec2d CurlNuBContribution(const Real &nu_b_diff, const Vec2d &grad_w_ijV_j)
{
    return Vec2d(nu_b_diff * grad_w_ijV_j[1], -nu_b_diff * grad_w_ijV_j[0]);
}
} // namespace
//=================================================================================================//
InitializeAphiElectromagneticVariables::
    InitializeAphiElectromagneticVariables(SPHBody &sph_body,
                                           Real default_conductivity,
                                           Real default_rho_cp,
                                           Real default_magnetic_reluctivity)
    : LocalDynamics(sph_body),
      default_conductivity_(default_conductivity),
      default_rho_cp_(default_rho_cp),
      default_magnetic_reluctivity_(default_magnetic_reluctivity),
      vector_potential_(particles_->registerStateVariableData<Vecd>("VectorPotential")),
      vector_potential_prev_(particles_->registerStateVariableData<Vecd>("VectorPotentialPrevious")),
      vector_potential_dt_(particles_->registerStateVariableData<Vecd>("VectorPotentialTimeDerivative")),
      electric_field_(particles_->registerStateVariableData<Vecd>("ElectricField")),
      current_density_(particles_->registerStateVariableData<Vecd>("CurrentDensity")),
      electric_potential_gradient_(particles_->registerStateVariableData<Vecd>("ElectricPotentialGradient")),
      source_current_density_(particles_->registerStateVariableData<Vecd>("SourceCurrentDensity")),
      vector_potential_change_rate_(particles_->registerStateVariableData<Vecd>("VectorPotentialChangeRate")),
      curl_nu_b_(particles_->registerStateVariableData<Vecd>("CurlNuB")),
      vector_potential_curl_(particles_->registerStateVariableData<AngularVecd>("VectorPotentialCurl")),
      electric_potential_(particles_->registerStateVariableData<Real>("ElectricPotential")),
      electric_potential_source_(particles_->registerStateVariableData<Real>("ElectricPotentialSource")),
      electric_potential_change_rate_(particles_->registerStateVariableData<Real>("ElectricPotentialChangeRate")),
      joule_heat_source_(particles_->registerStateVariableData<Real>("JouleHeatSource")),
      temperature_change_rate_by_joule_(particles_->registerStateVariableData<Real>("TemperatureChangeRateByJoule")),
      electrical_conductivity_(particles_->registerStateVariableData<Real>("ElectricalConductivity", default_conductivity_)),
      rho_cp_(particles_->registerStateVariableData<Real>("RhoCp", default_rho_cp_)),
      magnetic_reluctivity_(particles_->registerStateVariableData<Real>("MagneticReluctivity", default_magnetic_reluctivity_))
{
    particles_->addVariableToWrite<Vecd>("VectorPotential");
    particles_->addVariableToWrite<Vecd>("VectorPotentialTimeDerivative");
    particles_->addVariableToWrite<Vecd>("SourceCurrentDensity");
    particles_->addVariableToWrite<Vecd>("VectorPotentialChangeRate");
    particles_->addVariableToWrite<Vecd>("CurlNuB");
    particles_->addVariableToWrite<AngularVecd>("VectorPotentialCurl");
    particles_->addVariableToWrite<Real>("ElectricPotential");
    particles_->addVariableToWrite<Vecd>("ElectricField");
    particles_->addVariableToWrite<Vecd>("CurrentDensity");
    particles_->addVariableToWrite<Real>("JouleHeatSource");
    particles_->addVariableToWrite<Real>("TemperatureChangeRateByJoule");
    particles_->addVariableToWrite<Real>("ElectricalConductivity");
    particles_->addVariableToWrite<Real>("RhoCp");
    particles_->addVariableToWrite<Real>("MagneticReluctivity");
}
//=================================================================================================//
void InitializeAphiElectromagneticVariables::update(size_t index_i, Real dt)
{
    if (electrical_conductivity_[index_i] <= TinyReal)
    {
        electrical_conductivity_[index_i] = default_conductivity_;
    }
    if (rho_cp_[index_i] <= TinyReal)
    {
        rho_cp_[index_i] = default_rho_cp_;
    }
    if (magnetic_reluctivity_[index_i] <= TinyReal)
    {
        magnetic_reluctivity_[index_i] = default_magnetic_reluctivity_;
    }

    vector_potential_prev_[index_i] = vector_potential_[index_i];
    vector_potential_dt_[index_i] = ZeroData<Vecd>::value;
    electric_field_[index_i] = ZeroData<Vecd>::value;
    current_density_[index_i] = ZeroData<Vecd>::value;
    electric_potential_gradient_[index_i] = ZeroData<Vecd>::value;
    source_current_density_[index_i] = ZeroData<Vecd>::value;
    vector_potential_change_rate_[index_i] = ZeroData<Vecd>::value;
    curl_nu_b_[index_i] = ZeroData<Vecd>::value;
    vector_potential_curl_[index_i] = ZeroData<AngularVecd>::value;
    electric_potential_[index_i] = 0.0;
    electric_potential_source_[index_i] = 0.0;
    electric_potential_change_rate_[index_i] = 0.0;
    joule_heat_source_[index_i] = 0.0;
    temperature_change_rate_by_joule_[index_i] = 0.0;
}
//=================================================================================================//
UpdateVectorPotentialTimeDerivative::
    UpdateVectorPotentialTimeDerivative(SPHBody &sph_body)
    : LocalDynamics(sph_body),
      vector_potential_(particles_->getVariableDataByName<Vecd>("VectorPotential")),
      vector_potential_prev_(particles_->getVariableDataByName<Vecd>("VectorPotentialPrevious")),
      vector_potential_dt_(particles_->getVariableDataByName<Vecd>("VectorPotentialTimeDerivative")) {}
//=================================================================================================//
void UpdateVectorPotentialTimeDerivative::update(size_t index_i, Real dt)
{
    if (dt > TinyReal)
    {
        vector_potential_dt_[index_i] = (vector_potential_[index_i] - vector_potential_prev_[index_i]) / dt;
    }
    else
    {
        vector_potential_dt_[index_i] = ZeroData<Vecd>::value;
    }
    vector_potential_prev_[index_i] = vector_potential_[index_i];
}
//=================================================================================================//
ElectricPotentialSourceTermInner::
    ElectricPotentialSourceTermInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      electric_potential_source_(particles_->getVariableDataByName<Real>("ElectricPotentialSource")),
      vector_potential_dt_(particles_->getVariableDataByName<Vecd>("VectorPotentialTimeDerivative")) {}
//=================================================================================================//
void ElectricPotentialSourceTermInner::interaction(size_t index_i, Real dt)
{
    Real source_i = 0.0;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * electrical_conductivity_[index_j] /
                        (electrical_conductivity_[index_i] + electrical_conductivity_[index_j] + TinyReal);
        Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        if (!std::isfinite(gradW_ijV_j.squaredNorm()))
        {
            continue;
        }
        Vecd adot_diff = vector_potential_dt_[index_i] - vector_potential_dt_[index_j];
        source_i -= sigma_ij * adot_diff.dot(gradW_ijV_j);
    }
    electric_potential_source_[index_i] = source_i;
}
//=================================================================================================//
ElectricPotentialSourceTermContact::
    ElectricPotentialSourceTermContact(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      electric_potential_source_(particles_->getVariableDataByName<Real>("ElectricPotentialSource")),
      vector_potential_dt_(particles_->getVariableDataByName<Vecd>("VectorPotentialTimeDerivative"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_electrical_conductivity_.push_back(contact_particles->getVariableDataByName<Real>("ElectricalConductivity"));
        contact_vector_potential_dt_.push_back(contact_particles->getVariableDataByName<Vecd>("VectorPotentialTimeDerivative"));
    }
}
//=================================================================================================//
void ElectricPotentialSourceTermContact::interaction(size_t index_i, Real dt)
{
    Real source_contact = 0.0;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real *contact_vol_k = contact_vol_[k];
        Real *contact_sigma_k = contact_electrical_conductivity_[k];
        Vecd *contact_adot_k = contact_vector_potential_dt_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * contact_sigma_k[index_j] /
                            (electrical_conductivity_[index_i] + contact_sigma_k[index_j] + TinyReal);
            Vecd gradW_ijV_j = contact_neighborhood.dW_ij_[n] *
                               contact_vol_k[index_j] * contact_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }
            Vecd adot_diff = vector_potential_dt_[index_i] - contact_adot_k[index_j];
            source_contact -= sigma_ij * adot_diff.dot(gradW_ijV_j);
        }
    }
    electric_potential_source_[index_i] += source_contact;
}
//=================================================================================================//
ConstrainElectricPotential::
    ConstrainElectricPotential(BodyPartByParticle &body_part, Real reference_value)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      reference_value_(reference_value),
      electric_potential_(particles_->getVariableDataByName<Real>("ElectricPotential")) {}
//=================================================================================================//
void ConstrainElectricPotential::update(size_t index_i, Real dt)
{
    electric_potential_[index_i] = reference_value_;
}
//=================================================================================================//
ConstrainVectorPotential::
    ConstrainVectorPotential(BodyPartByParticle &body_part, const Vecd &reference_value)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      reference_value_(reference_value),
      vector_potential_(particles_->getVariableDataByName<Vecd>("VectorPotential")) {}
//=================================================================================================//
void ConstrainVectorPotential::update(size_t index_i, Real dt)
{
    vector_potential_[index_i] = reference_value_;
}
//=================================================================================================//
SetConstantElectromagneticMaterialProperties::
    SetConstantElectromagneticMaterialProperties(SPHBody &sph_body,
                                                 Real electrical_conductivity,
                                                 Real rho_cp,
                                                 Real magnetic_reluctivity)
    : LocalDynamics(sph_body),
      electrical_conductivity_value_(electrical_conductivity),
      rho_cp_value_(rho_cp),
      magnetic_reluctivity_value_(magnetic_reluctivity),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      rho_cp_(particles_->getVariableDataByName<Real>("RhoCp")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")) {}
//=================================================================================================//
void SetConstantElectromagneticMaterialProperties::update(size_t index_i, Real dt)
{
    electrical_conductivity_[index_i] = electrical_conductivity_value_;
    rho_cp_[index_i] = rho_cp_value_;
    magnetic_reluctivity_[index_i] = magnetic_reluctivity_value_;
}
//=================================================================================================//
UpdateElectricalConductivityByLinearTemperature::
    UpdateElectricalConductivityByLinearTemperature(SPHBody &sph_body,
                                                    Real sigma_ref,
                                                    Real alpha,
                                                    Real T_ref)
    : LocalDynamics(sph_body),
      sigma_ref_(sigma_ref), alpha_(alpha), T_ref_(T_ref),
      temperature_(particles_->getVariableDataByName<Real>("Temperature")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")) {}
//=================================================================================================//
void UpdateElectricalConductivityByLinearTemperature::update(size_t index_i, Real dt)
{
    Real sigma = sigma_ref_ * (1.0 - alpha_ * (temperature_[index_i] - T_ref_));
    electrical_conductivity_[index_i] = SMAX(sigma, TinyReal);
}
//=================================================================================================//
PrescribedSourceCurrentDensity::
    PrescribedSourceCurrentDensity(SPHBody &sph_body, const Vecd &source_current_density)
    : LocalDynamics(sph_body),
      source_current_density_value_(source_current_density),
      source_current_density_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensity")) {}
//=================================================================================================//
void PrescribedSourceCurrentDensity::update(size_t index_i, Real dt)
{
    source_current_density_[index_i] = source_current_density_value_;
}
//=================================================================================================//
PrescribedSourceCurrentDensityByBodyPart::
    PrescribedSourceCurrentDensityByBodyPart(BodyPartByParticle &body_part,
                                             const Vecd &source_current_density)
    : BaseLocalDynamics<BodyPartByParticle>(body_part),
      source_current_density_value_(source_current_density),
      source_current_density_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensity")) {}
//=================================================================================================//
void PrescribedSourceCurrentDensityByBodyPart::update(size_t index_i, Real dt)
{
    source_current_density_[index_i] = source_current_density_value_;
}
//=================================================================================================//
ElectricPotentialRelaxationInner::
    ElectricPotentialRelaxationInner(BaseInnerRelation &inner_relation,
                                     Real relaxation_scaling)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      smoothing_length_(inner_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength()),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      electric_potential_(particles_->getVariableDataByName<Real>("ElectricPotential")),
      electric_potential_source_(particles_->getVariableDataByName<Real>("ElectricPotentialSource")),
      electric_potential_change_rate_(particles_->getVariableDataByName<Real>("ElectricPotentialChangeRate"))
{
    (void)relaxation_scaling;
}
//=================================================================================================//
void ElectricPotentialRelaxationInner::interaction(size_t index_i, Real dt)
{
    Real potential_change_rate = electric_potential_source_[index_i];
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * electrical_conductivity_[index_j] /
                        (electrical_conductivity_[index_i] + electrical_conductivity_[index_j] + TinyReal);
        Real dW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j];
        Real r_ij = inner_neighborhood.r_ij_[n];
        Real surface_area_ij = 2.0 * dW_ijV_j / (r_ij + 0.01 * smoothing_length_);
        Real phi_diff = electric_potential_[index_j] - electric_potential_[index_i];
        potential_change_rate += sigma_ij * phi_diff * surface_area_ij;
    }
    electric_potential_change_rate_[index_i] = potential_change_rate;
}
//=================================================================================================//
ElectricPotentialRelaxationContact::
    ElectricPotentialRelaxationContact(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      smoothing_length_(contact_relation.getSPHBody().getSPHAdaptation().ReferenceSmoothingLength()),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      electric_potential_change_rate_(particles_->getVariableDataByName<Real>("ElectricPotentialChangeRate")),
      electric_potential_(particles_->getVariableDataByName<Real>("ElectricPotential"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_electrical_conductivity_.push_back(contact_particles->getVariableDataByName<Real>("ElectricalConductivity"));
        contact_electric_potential_.push_back(contact_particles->getVariableDataByName<Real>("ElectricPotential"));
    }
}
//=================================================================================================//
void ElectricPotentialRelaxationContact::interaction(size_t index_i, Real dt)
{
    Real potential_change_rate_contact = 0.0;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real *contact_vol_k = contact_vol_[k];
        Real *contact_sigma_k = contact_electrical_conductivity_[k];
        Real *contact_phi_k = contact_electric_potential_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            Real sigma_ij = 2.0 * electrical_conductivity_[index_i] * contact_sigma_k[index_j] /
                            (electrical_conductivity_[index_i] + contact_sigma_k[index_j] + TinyReal);
            Real dW_ijV_j = contact_neighborhood.dW_ij_[n] * contact_vol_k[index_j];
            Real r_ij = contact_neighborhood.r_ij_[n];
            Real surface_area_ij = 2.0 * dW_ijV_j / (r_ij + 0.01 * smoothing_length_);
            Real phi_diff = contact_phi_k[index_j] - electric_potential_[index_i];
            potential_change_rate_contact += sigma_ij * phi_diff * surface_area_ij;
        }
    }
    electric_potential_change_rate_[index_i] += potential_change_rate_contact;
}
//=================================================================================================//
UpdateElectricPotentialByRelaxationRate::
    UpdateElectricPotentialByRelaxationRate(SPHBody &sph_body,
                                            Real relaxation_scaling,
                                            Real max_abs_potential)
    : LocalDynamics(sph_body),
      relaxation_scaling_(relaxation_scaling),
      max_abs_potential_(max_abs_potential),
      electric_potential_(particles_->getVariableDataByName<Real>("ElectricPotential")),
      electric_potential_change_rate_(particles_->getVariableDataByName<Real>("ElectricPotentialChangeRate")) {}
//=================================================================================================//
void UpdateElectricPotentialByRelaxationRate::update(size_t index_i, Real dt)
{
    Real delta_phi = relaxation_scaling_ * dt * electric_potential_change_rate_[index_i];
    if (!std::isfinite(delta_phi))
    {
        delta_phi = 0.0;
    }
    electric_potential_[index_i] += delta_phi;
    if (!std::isfinite(electric_potential_[index_i]))
    {
        electric_potential_[index_i] = 0.0;
    }
    electric_potential_[index_i] =
        SMIN(max_abs_potential_, SMAX(-max_abs_potential_, electric_potential_[index_i]));
}
//=================================================================================================//
VectorPotentialCurlInner::
    VectorPotentialCurlInner(BaseInnerRelation &inner_relation,
                             const std::string &vector_potential_name,
                             const std::string &output_name,
                             Real curl_scaling)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      curl_scaling_(curl_scaling),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name)),
      vector_potential_curl_(particles_->registerStateVariableData<AngularVecd>(output_name)) {}
//=================================================================================================//
void VectorPotentialCurlInner::interaction(size_t index_i, Real dt)
{
    AngularVecd curl_i = ZeroData<AngularVecd>::value;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        if (!std::isfinite(gradW_ijV_j.squaredNorm()))
        {
            continue;
        }
        Vecd vector_potential_diff = vector_potential_[index_i] - vector_potential_[index_j];
        curl_i += getCrossProduct(vector_potential_diff, gradW_ijV_j);
    }
    vector_potential_curl_[index_i] = curl_scaling_ * curl_i;
}
//=================================================================================================//
VectorPotentialCurlContact::
    VectorPotentialCurlContact(BaseContactRelation &contact_relation,
                               const std::string &vector_potential_name,
                               const std::string &output_name,
                               Real curl_scaling)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      curl_scaling_(curl_scaling),
      vector_potential_(particles_->getVariableDataByName<Vecd>(vector_potential_name)),
      vector_potential_curl_(particles_->getVariableDataByName<AngularVecd>(output_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_vector_potential_.push_back(contact_particles->getVariableDataByName<Vecd>(vector_potential_name));
    }
}
//=================================================================================================//
void VectorPotentialCurlContact::interaction(size_t index_i, Real dt)
{
    AngularVecd curl_contact = ZeroData<AngularVecd>::value;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real *contact_vol_k = contact_vol_[k];
        Vecd *contact_vector_potential_k = contact_vector_potential_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Vecd gradW_ijV_j = contact_neighborhood.dW_ij_[n] *
                               contact_vol_k[index_j] * contact_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }
            Vecd vector_potential_diff = vector_potential_[index_i] - contact_vector_potential_k[index_j];
            curl_contact += getCrossProduct(vector_potential_diff, gradW_ijV_j);
        }
    }
    vector_potential_curl_[index_i] += curl_scaling_ * curl_contact;
}
//=================================================================================================//
CurlNuBInner::
    CurlNuBInner(BaseInnerRelation &inner_relation,
                 const std::string &magnetic_flux_density_name,
                 const std::string &output_name,
                 Real curl_scaling,
                 bool symmetric_volume_weight)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      curl_scaling_(curl_scaling),
      symmetric_volume_weight_(symmetric_volume_weight),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      magnetic_flux_density_(particles_->getVariableDataByName<AngularVecd>(magnetic_flux_density_name)),
      curl_nu_b_(particles_->registerStateVariableData<Vecd>(output_name)) {}
//=================================================================================================//
void CurlNuBInner::interaction(size_t index_i, Real dt)
{
    Vecd curl_nu_b_i = ZeroData<Vecd>::value;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Real vol_factor = symmetric_volume_weight_
                              ? 0.5 * (Vol_[index_i] + Vol_[index_j])
                              : Vol_[index_j];
        Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * vol_factor * inner_neighborhood.e_ij_[n];
        if (!std::isfinite(gradW_ijV_j.squaredNorm()))
        {
            continue;
        }
        AngularVecd nu_b_diff =
            magnetic_reluctivity_[index_i] * magnetic_flux_density_[index_i] -
            magnetic_reluctivity_[index_j] * magnetic_flux_density_[index_j];
        curl_nu_b_i += CurlNuBContribution(nu_b_diff, gradW_ijV_j);
    }
    curl_nu_b_[index_i] = curl_scaling_ * curl_nu_b_i;
}
//=================================================================================================//
CurlNuBContact::
    CurlNuBContact(BaseContactRelation &contact_relation,
                   const std::string &magnetic_flux_density_name,
                   const std::string &output_name,
                   Real curl_scaling,
                   bool symmetric_volume_weight)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      curl_scaling_(curl_scaling),
      symmetric_volume_weight_(symmetric_volume_weight),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      magnetic_reluctivity_(particles_->getVariableDataByName<Real>("MagneticReluctivity")),
      magnetic_flux_density_(particles_->getVariableDataByName<AngularVecd>(magnetic_flux_density_name)),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>(output_name))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_magnetic_reluctivity_.push_back(contact_particles->getVariableDataByName<Real>("MagneticReluctivity"));
        contact_magnetic_flux_density_.push_back(contact_particles->getVariableDataByName<AngularVecd>(magnetic_flux_density_name));
    }
}
//=================================================================================================//
void CurlNuBContact::interaction(size_t index_i, Real dt)
{
    Vecd curl_contact = ZeroData<Vecd>::value;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real *contact_vol_k = contact_vol_[k];
        Real *contact_nu_k = contact_magnetic_reluctivity_[k];
        AngularVecd *contact_b_k = contact_magnetic_flux_density_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Real vol_factor = symmetric_volume_weight_
                                  ? 0.5 * (Vol_[index_i] + contact_vol_k[index_j])
                                  : contact_vol_k[index_j];
            Vecd gradW_ijV_j =
                contact_neighborhood.dW_ij_[n] * vol_factor * contact_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }
            AngularVecd nu_b_diff =
                magnetic_reluctivity_[index_i] * magnetic_flux_density_[index_i] -
                contact_nu_k[index_j] * contact_b_k[index_j];
            curl_contact += CurlNuBContribution(nu_b_diff, gradW_ijV_j);
        }
    }
    curl_nu_b_[index_i] += curl_scaling_ * curl_contact;
}
//=================================================================================================//
ElectricPotentialGradientInner::
    ElectricPotentialGradientInner(BaseInnerRelation &inner_relation)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electric_potential_(particles_->getVariableDataByName<Real>("ElectricPotential")),
      electric_potential_gradient_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradient")) {}
//=================================================================================================//
void ElectricPotentialGradientInner::interaction(size_t index_i, Real dt)
{
    Vecd grad_phi = ZeroData<Vecd>::value;
    Neighborhood &inner_neighborhood = inner_configuration_[index_i];
    for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
    {
        size_t index_j = inner_neighborhood.j_[n];
        if (inner_neighborhood.r_ij_[n] < TinyReal)
        {
            continue;
        }
        Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
        if (!std::isfinite(gradW_ijV_j.squaredNorm()))
        {
            continue;
        }
        grad_phi -= (electric_potential_[index_i] - electric_potential_[index_j]) * gradW_ijV_j;
    }
    electric_potential_gradient_[index_i] = grad_phi;
}
//=================================================================================================//
ElectricPotentialGradientContact::
    ElectricPotentialGradientContact(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()), DataDelegateContact(contact_relation),
      electric_potential_(particles_->getVariableDataByName<Real>("ElectricPotential")),
      electric_potential_gradient_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradient"))
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        BaseParticles *contact_particles = contact_particles_[k];
        contact_vol_.push_back(contact_particles->getVariableDataByName<Real>("VolumetricMeasure"));
        contact_electric_potential_.push_back(contact_particles->getVariableDataByName<Real>("ElectricPotential"));
    }
}
//=================================================================================================//
void ElectricPotentialGradientContact::interaction(size_t index_i, Real dt)
{
    Vecd grad_phi_contact = ZeroData<Vecd>::value;
    for (size_t k = 0; k != contact_configuration_.size(); ++k)
    {
        Real *contact_vol_k = contact_vol_[k];
        Real *contact_phi_k = contact_electric_potential_[k];
        Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
        for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
        {
            size_t index_j = contact_neighborhood.j_[n];
            if (contact_neighborhood.r_ij_[n] < TinyReal)
            {
                continue;
            }
            Vecd gradW_ijV_j = contact_neighborhood.dW_ij_[n] *
                               contact_vol_k[index_j] * contact_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }
            grad_phi_contact -= (electric_potential_[index_i] - contact_phi_k[index_j]) * gradW_ijV_j;
        }
    }
    electric_potential_gradient_[index_i] += grad_phi_contact;
}
//=================================================================================================//
VectorPotentialEquationInner::
    VectorPotentialEquationInner(BaseInnerRelation &inner_relation,
                                 Real relaxation_scaling,
                                 Real max_change_rate)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      vector_potential_(particles_->getVariableDataByName<Vecd>("VectorPotential")),
      source_current_density_(particles_->getVariableDataByName<Vecd>("SourceCurrentDensity")),
      electric_potential_gradient_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradient")),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>("CurlNuB")),
      vector_potential_change_rate_(particles_->getVariableDataByName<Vecd>("VectorPotentialChangeRate")) {}
//=================================================================================================//
void VectorPotentialEquationInner::interaction(size_t index_i, Real dt)
{
    Real sigma_i = electrical_conductivity_[index_i];
    Vecd rhs = source_current_density_[index_i] - curl_nu_b_[index_i] -
               sigma_i * electric_potential_gradient_[index_i];
    Vecd a_rate = rhs / (sigma_i + TinyReal);
    Real a_rate_norm_sq = a_rate.squaredNorm();
    if (!std::isfinite(a_rate_norm_sq))
    {
        vector_potential_change_rate_[index_i] = ZeroData<Vecd>::value;
        return;
    }
    Real a_rate_norm = sqrt(a_rate_norm_sq);
    if (a_rate_norm > max_change_rate_)
    {
        a_rate *= max_change_rate_ / (a_rate_norm + TinyReal);
    }
    vector_potential_change_rate_[index_i] = a_rate;
}
//=================================================================================================//
void VectorPotentialEquationInner::update(size_t index_i, Real dt)
{
    vector_potential_[index_i] += relaxation_scaling_ * dt * vector_potential_change_rate_[index_i];
    Real vector_potential_norm_sq = vector_potential_[index_i].squaredNorm();
    if (!std::isfinite(vector_potential_norm_sq))
    {
        vector_potential_[index_i] = ZeroData<Vecd>::value;
    }
}
//=================================================================================================//
VectorPotentialEquationMagneticOnlyInner::
    VectorPotentialEquationMagneticOnlyInner(BaseInnerRelation &inner_relation,
                                             Real relaxation_scaling,
                                             Real max_change_rate)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      relaxation_scaling_(relaxation_scaling),
      max_change_rate_(max_change_rate),
      vector_potential_(particles_->getVariableDataByName<Vecd>("VectorPotential")),
      curl_nu_b_(particles_->getVariableDataByName<Vecd>("CurlNuB")),
      vector_potential_change_rate_(particles_->getVariableDataByName<Vecd>("VectorPotentialChangeRate")) {}
//=================================================================================================//
void VectorPotentialEquationMagneticOnlyInner::interaction(size_t index_i, Real dt)
{
    Vecd a_rate = -curl_nu_b_[index_i];
    Real a_rate_norm_sq = a_rate.squaredNorm();
    if (!std::isfinite(a_rate_norm_sq))
    {
        vector_potential_change_rate_[index_i] = ZeroData<Vecd>::value;
        return;
    }
    Real a_rate_norm = sqrt(a_rate_norm_sq);
    if (a_rate_norm > max_change_rate_)
    {
        a_rate *= max_change_rate_ / (a_rate_norm + TinyReal);
    }
    vector_potential_change_rate_[index_i] = a_rate;
}
//=================================================================================================//
void VectorPotentialEquationMagneticOnlyInner::update(size_t index_i, Real dt)
{
    vector_potential_[index_i] += relaxation_scaling_ * dt * vector_potential_change_rate_[index_i];
}
//=================================================================================================//
ElectricFieldCurrentAndJouleHeatInner::
    ElectricFieldCurrentAndJouleHeatInner(BaseInnerRelation &inner_relation,
                                          bool recompute_potential_gradient)
    : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
      recompute_potential_gradient_(recompute_potential_gradient),
      Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
      electrical_conductivity_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      rho_cp_(particles_->getVariableDataByName<Real>("RhoCp")),
      electric_potential_(particles_->getVariableDataByName<Real>("ElectricPotential")),
      joule_heat_source_(particles_->getVariableDataByName<Real>("JouleHeatSource")),
      temperature_change_rate_by_joule_(particles_->getVariableDataByName<Real>("TemperatureChangeRateByJoule")),
      vector_potential_dt_(particles_->getVariableDataByName<Vecd>("VectorPotentialTimeDerivative")),
      electric_field_(particles_->getVariableDataByName<Vecd>("ElectricField")),
      current_density_(particles_->getVariableDataByName<Vecd>("CurrentDensity")),
      electric_potential_gradient_(particles_->getVariableDataByName<Vecd>("ElectricPotentialGradient")) {}
//=================================================================================================//
void ElectricFieldCurrentAndJouleHeatInner::interaction(size_t index_i, Real dt)
{
    Vecd grad_phi = electric_potential_gradient_[index_i];
    if (recompute_potential_gradient_)
    {
        grad_phi = ZeroData<Vecd>::value;
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            if (!std::isfinite(gradW_ijV_j.squaredNorm()))
            {
                continue;
            }
            grad_phi -= (electric_potential_[index_i] - electric_potential_[index_j]) * gradW_ijV_j;
        }
        electric_potential_gradient_[index_i] = grad_phi;
    }
    Vecd electric_field = -vector_potential_dt_[index_i] - grad_phi;
    if (!std::isfinite(electric_field.squaredNorm()))
    {
        electric_field_[index_i] = ZeroData<Vecd>::value;
        current_density_[index_i] = ZeroData<Vecd>::value;
        joule_heat_source_[index_i] = 0.0;
        temperature_change_rate_by_joule_[index_i] = 0.0;
        return;
    }
    electric_field_[index_i] = electric_field;
    Vecd current_density = electrical_conductivity_[index_i] * electric_field;
    if (!std::isfinite(current_density.squaredNorm()))
    {
        current_density_[index_i] = ZeroData<Vecd>::value;
        joule_heat_source_[index_i] = 0.0;
        temperature_change_rate_by_joule_[index_i] = 0.0;
        return;
    }
    current_density_[index_i] = current_density;
    Real joule_heat = current_density.dot(electric_field);
    if (!std::isfinite(joule_heat))
    {
        joule_heat_source_[index_i] = 0.0;
        temperature_change_rate_by_joule_[index_i] = 0.0;
        return;
    }
    joule_heat_source_[index_i] = joule_heat;
    temperature_change_rate_by_joule_[index_i] = joule_heat / (rho_cp_[index_i] + TinyReal);
}
//=================================================================================================//
AddJouleHeatToDiffusionSpeciesRate::
    AddJouleHeatToDiffusionSpeciesRate(SPHBody &sph_body,
                                       const std::string &target_species_change_rate_name)
    : LocalDynamics(sph_body),
      temperature_change_rate_by_joule_(particles_->getVariableDataByName<Real>("TemperatureChangeRateByJoule")),
      target_species_change_rate_(particles_->getVariableDataByName<Real>(target_species_change_rate_name)) {}
//=================================================================================================//
void AddJouleHeatToDiffusionSpeciesRate::update(size_t index_i, Real dt)
{
    target_species_change_rate_[index_i] += temperature_change_rate_by_joule_[index_i];
}
//=================================================================================================//
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_TEAM7_APHI_DYNAMICS_HPP
