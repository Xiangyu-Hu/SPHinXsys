#ifndef ACTIVE_MODEL_H
#define ACTIVE_MODEL_H

#include "complex_solid.h"
#include "elastic_dynamics.h"

namespace SPH
{
/**
 * @class ActiveModelSolid
 */
class ActiveModelSolid : public SaintVenantKirchhoffSolid
{
    StdLargeVec<Matd> active_strain_;

  public:
    explicit ActiveModelSolid(Real rho0, Real youngs_modulus, Real poisson_ratio);
    virtual ~ActiveModelSolid(){};

    /** initialize the local properties, fiber and sheet direction. */
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    /** second Piola-Kirchhoff stress related with green-lagrangian deformation tensor */
    virtual Matd StressPK1(Matd &deformation, size_t particle_index_i) override;
};
} // namespace SPH
#endif // ACTIVE_MODEL_H
