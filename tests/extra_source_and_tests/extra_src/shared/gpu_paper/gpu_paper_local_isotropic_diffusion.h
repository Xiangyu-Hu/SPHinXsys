/**
 * @file gpu_paper_local_isotropic_diffusion.h
 * @brief GPU-paper overlay: LocalIsotropicDiffusion with CK InterParticleCoeff.
 * @note Lives under extra_src — does not modify src/. Used by A3 windows-frame.
 */
#ifndef SPHINXSYS_GPU_PAPER_LOCAL_ISOTROPIC_DIFFUSION_H
#define SPHINXSYS_GPU_PAPER_LOCAL_ISOTROPIC_DIFFUSION_H

#include "diffusion_reaction.h"

namespace SPH
{
/**
 * @class GpuPaperLocalIsotropicDiffusion
 * @brief Same physics as LocalIsotropicDiffusion, plus DiscreteVariable + InterParticleCoeff
 *        required by DiffusionRelaxationCK (paper A3).
 */
class GpuPaperLocalIsotropicDiffusion : public IsotropicDiffusion
{
  protected:
    Real diff_max_;
    Real *local_diffusivity_;
    DiscreteVariable<Real> *dv_thermal_conductivity_;

  public:
    GpuPaperLocalIsotropicDiffusion(const std::string &species_name,
                                    const std::string &gradient_species_name,
                                    Real diff_background, Real diff_max, Real cv = 1.0);
    GpuPaperLocalIsotropicDiffusion(const std::string &species_name, Real diff_background, Real diff_max,
                                    Real cv = 1.0);
    explicit GpuPaperLocalIsotropicDiffusion(ConstructArgs<std::string, Real, Real> args);
    virtual ~GpuPaperLocalIsotropicDiffusion() {};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    DiscreteVariable<Real> *dvThermalConductivity() { return dv_thermal_conductivity_; };

    virtual Real getReferenceDiffusivity() override { return diff_max_; };
    virtual Real getDiffusionCoeffWithBoundary(size_t index_i) override { return local_diffusivity_[index_i]; };
    virtual Real getInterParticleDiffusionCoeff(size_t index_i, size_t index_j, const Vecd &e_ij) override
    {
        return 0.5 * (local_diffusivity_[index_i] + local_diffusivity_[index_j]);
    };

    class InterParticleCoeff
    {
        Real *local_diffusivity_;

      public:
        template <class ExecutionPolicy, class EncloserType>
        InterParticleCoeff(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : local_diffusivity_(encloser.dvThermalConductivity()->DelegatedData(ex_policy)){};
        Real operator()(size_t index_i, size_t index_j, const Vecd &e_ij)
        {
            return 0.5 * (local_diffusivity_[index_i] + local_diffusivity_[index_j]);
        };
        Real operator()(size_t index_i, size_t index_j)
        {
            return 0.5 * (local_diffusivity_[index_i] + local_diffusivity_[index_j]);
        };
    };
};
} // namespace SPH

#endif
