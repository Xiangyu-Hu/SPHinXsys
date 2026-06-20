#ifndef ELECTROMAGNETIC_OPHELIE_BIOT_SAVART_H
#define ELECTROMAGNETIC_OPHELIE_BIOT_SAVART_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/**
 * Global coil-to-glass Biot-Savart summation (not compact-support SPH neighbors).
 * For each glass particle i, sum over all coil source particles j.
 */
class ComputeOphelieCoilToGlassBiotSavartCK : public LocalDynamics
{
  public:
    ComputeOphelieCoilToGlassBiotSavartCK(SPHBody &glass_body, SPHBody &coil_body, const OphelieGlassFieldNames &glass_names,
                                          const OphelieCoilFieldNames &coil_names, const OphelieParameters &params);
    virtual ~ComputeOphelieCoilToGlassBiotSavartCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real coeff_;
        Real eps2_;
        size_t n_coil_;
        Vecd *glass_pos_;
        Vecd *a_coil_real_;
        Vecd *a_coil_imag_;
        Vecd *b_coil_real_;
        Vecd *b_coil_imag_;
        Vecd *coil_pos_;
        Vecd *j_src_real_;
        Real *coil_vol_;
    };

  protected:
    Real coeff_;
    Real eps2_;
    size_t n_coil_;
    DiscreteVariable<Vecd> *dv_glass_pos_;
    DiscreteVariable<Vecd> *dv_a_coil_real_;
    DiscreteVariable<Vecd> *dv_a_coil_imag_;
    DiscreteVariable<Vecd> *dv_b_coil_real_;
    DiscreteVariable<Vecd> *dv_b_coil_imag_;
    DiscreteVariable<Vecd> *dv_coil_pos_;
    DiscreteVariable<Vecd> *dv_j_src_real_;
    DiscreteVariable<Real> *dv_coil_vol_;
};

/** Glass self-induction: A_ind from conductive current density on glass particles. */
class ComputeOphelieGlassSelfInducedBiotSavartCK : public LocalDynamics
{
  public:
    ComputeOphelieGlassSelfInducedBiotSavartCK(SPHBody &glass_body, const OphelieGlassFieldNames &names,
                                               const OphelieParameters &params,
                                               const std::string &j_real_field = std::string(),
                                               const std::string &j_imag_field = std::string());

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real coeff_;
        Real eps2_;
        size_t n_glass_;
        Vecd *glass_pos_;
        Vecd *j_real_;
        Vecd *j_imag_;
        Real *glass_vol_;
        Vecd *a_ind_real_;
        Vecd *a_ind_imag_;
        Vecd *b_ind_real_;
        Vecd *b_ind_imag_;
    };

  protected:
    Real coeff_;
    Real eps2_;
    size_t n_glass_;
    DiscreteVariable<Vecd> *dv_glass_pos_;
    DiscreteVariable<Vecd> *dv_j_real_;
    DiscreteVariable<Vecd> *dv_j_imag_;
    DiscreteVariable<Real> *dv_glass_vol_;
    DiscreteVariable<Vecd> *dv_a_ind_real_;
    DiscreteVariable<Vecd> *dv_a_ind_imag_;
    DiscreteVariable<Vecd> *dv_b_ind_real_;
    DiscreteVariable<Vecd> *dv_b_ind_imag_;
};

class CombineOphelieCoilAndInducedVectorPotentialCK : public LocalDynamics
{
  public:
    CombineOphelieCoilAndInducedVectorPotentialCK(SPHBody &sph_body, const OphelieGlassFieldNames &names);

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Vecd *a_coil_real_;
        Vecd *a_coil_imag_;
        Vecd *b_coil_real_;
        Vecd *b_coil_imag_;
        Vecd *a_ind_real_;
        Vecd *a_ind_imag_;
        Vecd *b_ind_real_;
        Vecd *b_ind_imag_;
        Vecd *a_src_real_;
        Vecd *a_src_imag_;
        Vecd *b_src_real_;
        Vecd *b_src_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_a_coil_real_;
    DiscreteVariable<Vecd> *dv_a_coil_imag_;
    DiscreteVariable<Vecd> *dv_b_coil_real_;
    DiscreteVariable<Vecd> *dv_b_coil_imag_;
    DiscreteVariable<Vecd> *dv_a_ind_real_;
    DiscreteVariable<Vecd> *dv_a_ind_imag_;
    DiscreteVariable<Vecd> *dv_b_ind_real_;
    DiscreteVariable<Vecd> *dv_b_ind_imag_;
    DiscreteVariable<Vecd> *dv_a_src_real_;
    DiscreteVariable<Vecd> *dv_a_src_imag_;
    DiscreteVariable<Vecd> *dv_b_src_real_;
    DiscreteVariable<Vecd> *dv_b_src_imag_;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#include "electromagnetic_ophelie_biot_savart.hpp"
#endif // ELECTROMAGNETIC_OPHELIE_BIOT_SAVART_H
