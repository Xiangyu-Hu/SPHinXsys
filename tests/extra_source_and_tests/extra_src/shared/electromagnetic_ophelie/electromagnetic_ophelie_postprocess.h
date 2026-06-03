#ifndef ELECTROMAGNETIC_OPHELIE_POSTPROCESS_H
#define ELECTROMAGNETIC_OPHELIE_POSTPROCESS_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Level 0: E = -i omega A_src (no phi correction). */
class ComputeOphelieEJQFromASrcNoPhiCK : public LocalDynamics
{
  public:
    ComputeOphelieEJQFromASrcNoPhiCK(SPHBody &sph_body, const OphelieGlassFieldNames &names, const OphelieParameters &params);
    virtual ~ComputeOphelieEJQFromASrcNoPhiCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real omega_;
        Real *sigma_;
        Vecd *a_src_real_;
        Vecd *a_src_imag_;
        Vecd *e_real_;
        Vecd *e_imag_;
        Vecd *j_real_;
        Vecd *j_imag_;
        Real *joule_heat_;
    };

  protected:
    Real omega_;
    DiscreteVariable<Real> *dv_sigma_;
    DiscreteVariable<Vecd> *dv_a_src_real_;
    DiscreteVariable<Vecd> *dv_a_src_imag_;
    DiscreteVariable<Vecd> *dv_e_real_;
    DiscreteVariable<Vecd> *dv_e_imag_;
    DiscreteVariable<Vecd> *dv_j_real_;
    DiscreteVariable<Vecd> *dv_j_imag_;
    DiscreteVariable<Real> *dv_joule_heat_;
};

/** Scale E/J/A/B/Phi by field_scale and JouleHeat by power_scale so Q = 0.5 J·E stays consistent. */
class ScaleOphelieElectromagneticFieldsCK : public LocalDynamics
{
  public:
    ScaleOphelieElectromagneticFieldsCK(SPHBody &sph_body, const OphelieGlassFieldNames &names, Real field_scale,
                                        Real power_scale);
    virtual ~ScaleOphelieElectromagneticFieldsCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real field_scale_;
        Real power_scale_;
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
        Vecd *e_real_;
        Vecd *e_imag_;
        Vecd *j_real_;
        Vecd *j_imag_;
        Real *phi_imag_;
        Vecd *grad_phi_imag_;
        Real *joule_heat_;
    };

  protected:
    Real field_scale_;
    Real power_scale_;
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
    DiscreteVariable<Vecd> *dv_e_real_;
    DiscreteVariable<Vecd> *dv_e_imag_;
    DiscreteVariable<Vecd> *dv_j_real_;
    DiscreteVariable<Vecd> *dv_j_imag_;
    DiscreteVariable<Real> *dv_phi_imag_;
    DiscreteVariable<Vecd> *dv_grad_phi_imag_;
    DiscreteVariable<Real> *dv_joule_heat_;
};

/** @deprecated Use ScaleOphelieElectromagneticFieldsCK with field_scale=sqrt(scale), power_scale=scale. */
class ScaleOphelieJouleHeatCK : public LocalDynamics
{
  public:
    ScaleOphelieJouleHeatCK(SPHBody &sph_body, const OphelieGlassFieldNames &names, Real scale);
    virtual ~ScaleOphelieJouleHeatCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real scale_;
        Real *joule_heat_;
    };

  protected:
    Real scale_;
    DiscreteVariable<Real> *dv_joule_heat_;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#include "electromagnetic_ophelie_postprocess.hpp"
#endif // ELECTROMAGNETIC_OPHELIE_POSTPROCESS_H
