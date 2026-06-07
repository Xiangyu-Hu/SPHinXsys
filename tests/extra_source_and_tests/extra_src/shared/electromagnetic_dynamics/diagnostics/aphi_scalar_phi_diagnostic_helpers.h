#ifndef APHI_SCALAR_PHI_DIAGNOSTIC_HELPERS_H
#define APHI_SCALAR_PHI_DIAGNOSTIC_HELPERS_H

#include "base_general_dynamics.h"
#include "electromagnetic_dynamics/aphi_coupling_modes_ck.h"
#include "electromagnetic_dynamics/aphi_field_names_ck.h"

#include <cmath>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

inline AphiLhsAssemblyOptions scalarPhiLaplacePenaltyOptions(Real phi_gauge_penalty, bool use_phi_gauge_penalty = true)
{
    AphiLhsAssemblyOptions options;
    options.terms.laplace_a = false;
    options.terms.laplace_phi = true;
    options.terms.reaction = false;
    options.terms.grad_phi_coupling = false;
    options.terms.div_sigma_a_coupling = false;
    options.use_phi_gauge_penalty = use_phi_gauge_penalty;
    options.phi_gauge_penalty = phi_gauge_penalty;
    return options;
}

/** Deterministic pseudo-random scalar in [-1, 1] from particle index and seed (no runtime RNG on device). */
inline Real deterministicPseudoRandomUnit(size_t index_i, UnsignedInt seed)
{
    UnsignedInt x = static_cast<UnsignedInt>(index_i + 1);
    x ^= seed * 0x9e3779b9u;
    x ^= x >> 16;
    x *= 0x7feb352du;
    x ^= x >> 15;
    x *= 0x846ca68bu;
    x ^= x >> 16;
    const Real unit01 = static_cast<Real>(x & 0x00ffffffu) / static_cast<Real>(0x01000000u);
    return Real(2) * unit01 - Real(1);
}

/** High-frequency rough field using seed-dependent wavenumbers. */
inline Real deterministicRoughScalarField(size_t index_i, UnsignedInt seed, const Vecd &position)
{
    const Real wave_x = Real(1 + (seed % 7));
    const Real wave_y = Real(2 + ((seed / 7) % 5));
    const Real wave_z = Real(1 + ((seed / 35) % 4));
    const Real hash = deterministicPseudoRandomUnit(index_i, seed + 911u);
    return hash * std::sin(wave_x * Pi * position[0]) * std::cos(wave_y * Pi * position[1]) +
           (Real(1) - std::abs(hash)) * std::sin(wave_z * Pi * position[2]);
}

class AssignRandomScalarPhiFieldsCK : public LocalDynamics
{
  public:
    AssignRandomScalarPhiFieldsCK(SPHBody &sph_body, const AphiBlockNames &block_x, const AphiBlockNames &block_y,
                                  UnsignedInt seed_x, UnsignedInt seed_y)
        : LocalDynamics(sph_body), seed_x_(seed_x), seed_y_(seed_y),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_x_a_real_(particles_->template getVariableByName<Vecd>(block_x.a_real)),
          dv_x_a_imag_(particles_->template getVariableByName<Vecd>(block_x.a_imag)),
          dv_x_phi_real_(particles_->template getVariableByName<Real>(block_x.phi_real)),
          dv_x_phi_imag_(particles_->template getVariableByName<Real>(block_x.phi_imag)),
          dv_y_a_real_(particles_->template getVariableByName<Vecd>(block_y.a_real)),
          dv_y_a_imag_(particles_->template getVariableByName<Vecd>(block_y.a_imag)),
          dv_y_phi_real_(particles_->template getVariableByName<Real>(block_y.phi_real)),
          dv_y_phi_imag_(particles_->template getVariableByName<Real>(block_y.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : seed_x_(encloser.seed_x_), seed_y_(encloser.seed_y_),
              position_(encloser.dv_position_->DelegatedData(ex_policy)),
              x_a_real_(encloser.dv_x_a_real_->DelegatedData(ex_policy)),
              x_a_imag_(encloser.dv_x_a_imag_->DelegatedData(ex_policy)),
              x_phi_real_(encloser.dv_x_phi_real_->DelegatedData(ex_policy)),
              x_phi_imag_(encloser.dv_x_phi_imag_->DelegatedData(ex_policy)),
              y_a_real_(encloser.dv_y_a_real_->DelegatedData(ex_policy)),
              y_a_imag_(encloser.dv_y_a_imag_->DelegatedData(ex_policy)),
              y_phi_real_(encloser.dv_y_phi_real_->DelegatedData(ex_policy)),
              y_phi_imag_(encloser.dv_y_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            const Vecd zero_vec = Vecd::Zero();
            x_a_real_[index_i] = zero_vec;
            x_a_imag_[index_i] = zero_vec;
            y_a_real_[index_i] = zero_vec;
            y_a_imag_[index_i] = zero_vec;
            x_phi_real_[index_i] = deterministicRoughScalarField(index_i, seed_x_, position);
            x_phi_imag_[index_i] = deterministicRoughScalarField(index_i, seed_x_ + 104729u, position);
            y_phi_real_[index_i] = deterministicRoughScalarField(index_i, seed_y_, position);
            y_phi_imag_[index_i] = deterministicRoughScalarField(index_i, seed_y_ + 104729u, position);
        }

      protected:
        UnsignedInt seed_x_;
        UnsignedInt seed_y_;
        Vecd *position_;
        Vecd *x_a_real_;
        Vecd *x_a_imag_;
        Real *x_phi_real_;
        Real *x_phi_imag_;
        Vecd *y_a_real_;
        Vecd *y_a_imag_;
        Real *y_phi_real_;
        Real *y_phi_imag_;
    };

  protected:
    UnsignedInt seed_x_;
    UnsignedInt seed_y_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_x_a_real_;
    DiscreteVariable<Vecd> *dv_x_a_imag_;
    DiscreteVariable<Real> *dv_x_phi_real_;
    DiscreteVariable<Real> *dv_x_phi_imag_;
    DiscreteVariable<Vecd> *dv_y_a_real_;
    DiscreteVariable<Vecd> *dv_y_a_imag_;
    DiscreteVariable<Real> *dv_y_phi_real_;
    DiscreteVariable<Real> *dv_y_phi_imag_;
};

class AssignRandomScalarPhiSingleBlockCK : public LocalDynamics
{
  public:
    explicit AssignRandomScalarPhiSingleBlockCK(SPHBody &sph_body, const AphiBlockNames &block_names, UnsignedInt seed)
        : LocalDynamics(sph_body), seed_(seed),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
          dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : seed_(encloser.seed_), position_(encloser.dv_position_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd zero_vec = Vecd::Zero();
            a_real_[index_i] = zero_vec;
            a_imag_[index_i] = zero_vec;
            phi_real_[index_i] = deterministicRoughScalarField(index_i, seed_, position_[index_i]);
            phi_imag_[index_i] = deterministicRoughScalarField(index_i, seed_ + 104729u, position_[index_i]);
        }

      protected:
        UnsignedInt seed_;
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    UnsignedInt seed_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

class AssignScalarPhiBasisVectorCK : public LocalDynamics
{
  public:
    AssignScalarPhiBasisVectorCK(SPHBody &sph_body, const AphiBlockNames &block_names, size_t basis_index,
                                 bool use_phi_imag)
        : LocalDynamics(sph_body), basis_index_(basis_index), use_phi_imag_(use_phi_imag),
          dv_a_real_(particles_->template getVariableByName<Vecd>(block_names.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(block_names.a_imag)),
          dv_phi_real_(particles_->template getVariableByName<Real>(block_names.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(block_names.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : basis_index_(encloser.basis_index_), use_phi_imag_(encloser.use_phi_imag_),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd zero_vec = Vecd::Zero();
            a_real_[index_i] = zero_vec;
            a_imag_[index_i] = zero_vec;
            phi_real_[index_i] = (!use_phi_imag_ && index_i == basis_index_) ? Real(1) : Real(0);
            phi_imag_[index_i] = (use_phi_imag_ && index_i == basis_index_) ? Real(1) : Real(0);
        }

      protected:
        size_t basis_index_;
        bool use_phi_imag_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    size_t basis_index_;
    bool use_phi_imag_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_SCALAR_PHI_DIAGNOSTIC_HELPERS_H
