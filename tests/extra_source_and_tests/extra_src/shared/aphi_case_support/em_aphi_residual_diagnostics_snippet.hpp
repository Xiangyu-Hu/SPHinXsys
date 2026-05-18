// Add to extra_src/shared/aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.h
// after VectorPotentialFrequencyCoupledBlockEquationComplex or near diagnostics utilities.
class FrequencyAEquationResidualDiagnostic : public LocalDynamics
{
  public:
    explicit FrequencyAEquationResidualDiagnostic(
        SPHBody &sph_body,
        Real angular_frequency,
        const std::string &a_real_name = "VectorPotentialReal",
        const std::string &a_imag_name = "VectorPotentialImag",
        const std::string &source_real_name = "SourceCurrentDensityReal",
        const std::string &source_imag_name = "SourceCurrentDensityImag",
        const std::string &grad_phi_real_name = "ElectricPotentialGradientReal",
        const std::string &grad_phi_imag_name = "ElectricPotentialGradientImag",
        const std::string &curl_nu_b_real_name = "CurlNuBReal",
        const std::string &curl_nu_b_imag_name = "CurlNuBImag",
        const std::string &residual_real_name = "AEquationResidualReal",
        const std::string &residual_imag_name = "AEquationResidualImag",
        const std::string &relative_residual_name = "AEquationRelativeResidual");
    virtual ~FrequencyAEquationResidualDiagnostic() {};
    void update(size_t index_i, Real dt = 0.0);

  protected:
    Real angular_frequency_;
    Real *sigma_;
    Vecd *a_real_, *a_imag_;
    Vecd *source_real_, *source_imag_;
    Vecd *grad_phi_real_, *grad_phi_imag_;
    Vecd *curl_nu_b_real_, *curl_nu_b_imag_;
    Vecd *residual_real_, *residual_imag_;
    Real *relative_residual_;
};

// Add to extra_src/shared/aphi_case_support/electromagnetic_team7_aphi_frequency_dynamics.hpp
FrequencyAEquationResidualDiagnostic::FrequencyAEquationResidualDiagnostic(
    SPHBody &sph_body,
    Real angular_frequency,
    const std::string &a_real_name,
    const std::string &a_imag_name,
    const std::string &source_real_name,
    const std::string &source_imag_name,
    const std::string &grad_phi_real_name,
    const std::string &grad_phi_imag_name,
    const std::string &curl_nu_b_real_name,
    const std::string &curl_nu_b_imag_name,
    const std::string &residual_real_name,
    const std::string &residual_imag_name,
    const std::string &relative_residual_name)
    : LocalDynamics(sph_body), angular_frequency_(angular_frequency),
      sigma_(particles_->getVariableDataByName<Real>("ElectricalConductivity")),
      a_real_(particles_->getVariableDataByName<Vecd>(a_real_name)),
      a_imag_(particles_->getVariableDataByName<Vecd>(a_imag_name)),
      source_real_(particles_->getVariableDataByName<Vecd>(source_real_name)),
      source_imag_(particles_->getVariableDataByName<Vecd>(source_imag_name)),
      grad_phi_real_(particles_->getVariableDataByName<Vecd>(grad_phi_real_name)),
      grad_phi_imag_(particles_->getVariableDataByName<Vecd>(grad_phi_imag_name)),
      curl_nu_b_real_(particles_->getVariableDataByName<Vecd>(curl_nu_b_real_name)),
      curl_nu_b_imag_(particles_->getVariableDataByName<Vecd>(curl_nu_b_imag_name)),
      residual_real_(particles_->registerStateVariableData<Vecd>(residual_real_name)),
      residual_imag_(particles_->registerStateVariableData<Vecd>(residual_imag_name)),
      relative_residual_(particles_->registerStateVariableData<Real>(relative_residual_name))
{
    particles_->addVariableToWrite<Vecd>(residual_real_name);
    particles_->addVariableToWrite<Vecd>(residual_imag_name);
    particles_->addVariableToWrite<Real>(relative_residual_name);
}

void FrequencyAEquationResidualDiagnostic::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real sigma_i = sigma_[index_i];
    Vecd omega_a_imag = sigma_i * angular_frequency_ * a_imag_[index_i];
    Vecd omega_a_real = sigma_i * angular_frequency_ * a_real_[index_i];

    residual_real_[index_i] = source_real_[index_i]
        - curl_nu_b_real_[index_i]
        - sigma_i * grad_phi_real_[index_i]
        + omega_a_imag;

    residual_imag_[index_i] = source_imag_[index_i]
        - curl_nu_b_imag_[index_i]
        - sigma_i * grad_phi_imag_[index_i]
        - omega_a_real;

    Real numerator = std::sqrt(residual_real_[index_i].squaredNorm() +
                               residual_imag_[index_i].squaredNorm());
    Real denominator = source_real_[index_i].norm() + source_imag_[index_i].norm() +
                       curl_nu_b_real_[index_i].norm() + curl_nu_b_imag_[index_i].norm() +
                       (sigma_i * grad_phi_real_[index_i]).norm() +
                       (sigma_i * grad_phi_imag_[index_i]).norm() +
                       omega_a_imag.norm() + omega_a_real.norm() + TinyReal;
    relative_residual_[index_i] = numerator / denominator;
}
