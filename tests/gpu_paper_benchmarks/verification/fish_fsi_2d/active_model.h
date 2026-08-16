#ifndef ACTIVE_MODEL_H
#define ACTIVE_MODEL_H

#include "complex_solid.h"
#include "elastic_dynamics.h"

namespace SPH
{
/**
 * @class ActiveModelSolid
 * @brief Elastic solid with active strain for fish muscle simulation.
 *        Uses multiplicative decomposition F = F_e * F0 where F0 comes
 *        from the imposed active strain tensor.
 */
class ActiveModelSolid : public SaintVenantKirchhoffSolid
{
  public:
    explicit ActiveModelSolid(Real rho0, Real youngs_modulus, Real poisson_ratio);
    virtual ~ActiveModelSolid(){};

    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    virtual Matd StressPK1(Matd &deformation, size_t particle_index_i) override;

    class ConstituteKernel : public SaintVenantKirchhoffSolid::ConstituteKernel
    {
      public:
        template <typename ExecutionPolicy>
        ConstituteKernel(const ExecutionPolicy &ex_policy, ActiveModelSolid &encloser);
        inline Matd StressPK1(const Matd &F, size_t index_i);

      protected:
        Matd *active_strain_;
    };

    DiscreteVariable<Matd> *dv_active_strain_;

  protected:
    Matd *active_strain_;
};
} // namespace SPH

#include "active_model.hpp"

#endif // ACTIVE_MODEL_H
