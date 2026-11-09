/**
 * @file gpu_paper_diffusion_robin.h
 * @brief GPU-paper overlay: Robin contact BC for DiffusionRelaxationCK (A3).
 * @note Lives under extra_src — does not modify src/ diffusion_dynamics_ck.*.
 */
#ifndef SPHINXSYS_GPU_PAPER_DIFFUSION_ROBIN_H
#define SPHINXSYS_GPU_PAPER_DIFFUSION_ROBIN_H

#include "diffusion_dynamics_ck.h"

namespace SPH
{
/**
 * @class Robin
 * @brief Contact Robin (convection) BC for DiffusionRelaxationCK.
 * transfer ~ h * (T_inf - T_i) * surface_area · (n_i - n_j), matching classic Robin form.
 *
 * Specializes the library forward `template <typename...> class Robin` (see sphinxsys_containers.h).
 */
template <class DiffusionType>
class Robin<DiffusionType>
{
  public:
    template <class DiffusionDynamics>
    Robin(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles);

    class ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
        Vecd operator()(UnsignedInt m, UnsignedInt index_i, UnsignedInt index_j, const Vecd &e_ij,
                        const Vecd &vec_r_ij);

      protected:
        Vecd *n_;
        Vecd *contact_n_;
        MultiEntryView<Real> species_;
        MultiEntryView<Real> contact_convection_;
        Real *contact_species_infinity_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_n_;
    DiscreteVariable<Vecd> *dv_contact_n_;
    DiscreteVariable<Real> *dv_species_;
    DiscreteVariable<Real> *dv_contact_convection_;
    SingleVariable<Real> *sv_contact_species_infinity_;
};
} // namespace SPH

#include "gpu_paper_diffusion_robin.hpp"

#endif
