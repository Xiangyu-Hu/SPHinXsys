#include "johnson_cook_solid.hpp"
#include "base_particles.hpp"
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace SPH
{
//=================================================================================================//
JohnsonCookSolid::JohnsonCookSolid(Real rho0, Real youngs_modulus, Real poisson_ratio,
                                   const JohnsonCookParameters &parameters)
    : ElasticSolid(rho0), parameters_(parameters), model_(rho0, youngs_modulus, poisson_ratio, parameters), stress_(nullptr), previous_deformation_(nullptr),
      previous_time_step_(nullptr), plastic_strain_(nullptr), plastic_strain_rate_(nullptr),
      temperature_(nullptr), damage_(nullptr), equivalent_stress_(nullptr), mean_stress_(nullptr),
      plastic_work_specific_(nullptr), material_active_(nullptr)
{
    const auto &p = parameters_;
    const Real values[] = {rho0, youngs_modulus, poisson_ratio, p.yield_stress, p.hardening_coefficient,
                           p.hardening_exponent, p.strain_rate_coefficient, p.reference_strain_rate,
                           p.thermal_softening_exponent, p.reference_temperature, p.melting_temperature,
                           p.specific_heat, p.failure_D1, p.failure_D2, p.failure_D3, p.failure_D4, p.failure_D5};
    for (Real value : values)
        if (!std::isfinite(value))
            throw std::invalid_argument("Johnson-Cook parameters must be finite");
    if (!(rho0 > 0 && youngs_modulus > 0 && poisson_ratio > -1 && poisson_ratio < Real(0.5) &&
          p.yield_stress >= 0 && p.hardening_coefficient >= 0 && p.hardening_exponent > 0 &&
          p.strain_rate_coefficient >= 0 && p.reference_strain_rate > 0 &&
          p.thermal_softening_exponent > 0 && p.melting_temperature > p.reference_temperature && p.specific_heat > 0))
        throw std::invalid_argument("Invalid Johnson-Cook parameter range");
    E0_ = youngs_modulus;
    nu_ = poisson_ratio;
    G0_ = E0_ / (2 * (1 + nu_));
    K0_ = E0_ / (3 * (1 - 2 * nu_));
    setSoundSpeeds();
    contact_stiffness_ = rho0_ * c0_ * c0_;
    material_type_name_ = "JohnsonCookSolid";
}
//=================================================================================================//
JohnsonCookState JohnsonCookSolid::InitialState() const
{
    JohnsonCookState state;
    state.temperature = parameters_.reference_temperature;
    return state;
}
//=================================================================================================//
Real JohnsonCookSolid::NormalizedTemperature(Real temperature) const
{
    return model_.NormalizedTemperature(temperature);
}
//=================================================================================================//
Real JohnsonCookSolid::FractureStrain(Real ratio, Real rate, Real temperature) const
{
    return model_.FractureStrain(ratio, rate, temperature);
}
//=================================================================================================//
void JohnsonCookSolid::Integrate(const Mat3d &increment, Real J, Real dt, JohnsonCookState &state) const
{
    if (model_.Integrate(increment, J, dt, state) != JohnsonCookIntegrationStatus::success)
        throw std::invalid_argument("Invalid Johnson-Cook increment or nonfinite state");
}
//=================================================================================================//
void JohnsonCookSolid::initializeLocalParameters(BaseParticles *particles)
{
    ElasticSolid::initializeLocalParameters(particles);
    particles_ = particles;
    particles->registerStateVariable<int>("JCIntegrationStatus");
    stress_ = particles->registerStateVariableData<Mat3d>("StressCauchy");
    previous_deformation_ = particles->registerStateVariableData<Matd>("JCPreviousDeformation", Matd::Identity().eval());
    previous_time_step_ = particles->registerStateVariableData<Real>("JCPreviousTimeStep");
    plastic_strain_ = particles->registerStateVariableData<Real>("JCPlasticStrain");
    plastic_strain_rate_ = particles->registerStateVariableData<Real>("JCPlasticStrainRate");
    temperature_ = particles->registerStateVariableData<Real>("JCTemperature", parameters_.reference_temperature);
    damage_ = particles->registerStateVariableData<Real>("JCDamage");
    equivalent_stress_ = particles->registerStateVariableData<Real>("VonMisesStress");
    mean_stress_ = particles->registerStateVariableData<Real>("MeanStress");
    plastic_work_specific_ = particles->registerStateVariableData<Real>("PlasticWorkSpecific");
    material_active_ = particles->registerStateVariableData<int>("MaterialActive", 1);
    for (const char *name : {"JCPreviousTimeStep", "JCPlasticStrain", "JCPlasticStrainRate", "JCTemperature",
                             "JCDamage", "VonMisesStress", "MeanStress", "PlasticWorkSpecific"})
        particles->addEvolvingVariable<Real>(name);
    particles->addEvolvingVariable<Mat3d>("StressCauchy");
    particles->addEvolvingVariable<Matd>("JCPreviousDeformation");
    particles->addEvolvingVariable<int>("MaterialActive");
    host_kernel_ = std::make_unique<ConstituteKernel>(execution::seq, *this);
}
//=================================================================================================//
Matd JohnsonCookSolid::UpdateStress(const Matd &F, size_t i, Real dt)
{
    const Matd stress = host_kernel_->UpdateStress(F, i, dt);
    if (host_kernel_->Status(i) != 0)
        throw std::domain_error("Invalid Johnson-Cook deformation or nonfinite state");
    return stress;
}
//=================================================================================================//
Matd JohnsonCookSolid::StressPK1(Matd &, size_t) { throw std::logic_error("Use JohnsonCookIntegration1stHalf"); }
Matd JohnsonCookSolid::StressPK2(Matd &, size_t) { throw std::logic_error("Use JohnsonCookIntegration1stHalf"); }
Matd JohnsonCookSolid::StressCauchy(Matd &, size_t) { throw std::logic_error("Use JohnsonCookIntegration1stHalf"); }
Real JohnsonCookSolid::VolumetricKirchhoff(Real) { throw std::logic_error("Johnson-Cook requires material history"); }
} // namespace SPH
