/**
 * @file spdb_contact_scaled.h
 * @brief Contact stiffness scaled materials for SPDB cases (extra_src only).
 */
#ifndef SPDB_CONTACT_SCALED_H
#define SPDB_CONTACT_SCALED_H

#include "spdb_ren2014_solid.h"

#include <cmath>

namespace SPH
{

class ContactScaledLinearElastic : public LinearElasticSolid
{
  public:
    ContactScaledLinearElastic(Real rho0, Real youngs_modulus, Real poisson_ratio, Real stiffness_scale)
        : LinearElasticSolid(rho0, youngs_modulus, poisson_ratio), stiffness_scale_(stiffness_scale)
    {
        setContactStiffness(rho0, ReferenceSoundSpeed() * std::sqrt(stiffness_scale_));
        material_type_name_ = "ContactScaledLinearElastic";
    }
    Real stiffnessScale() const { return stiffness_scale_; }

  private:
    Real stiffness_scale_;
};

class ContactScaledRen2014 : public Ren2014X20Cr13Solid
{
  public:
    ContactScaledRen2014(Real rho0, Real youngs_modulus, Real poisson_ratio, Real temperature_K,
                         bool inelastic_enabled, Real stiffness_scale)
        : Ren2014X20Cr13Solid(rho0, youngs_modulus, poisson_ratio, temperature_K, inelastic_enabled),
          stiffness_scale_(stiffness_scale)
    {
        setContactStiffness(rho0, ReferenceSoundSpeed() * std::sqrt(stiffness_scale_));
        material_type_name_ = "ContactScaledRen2014";
    }
    Real stiffnessScale() const { return stiffness_scale_; }

  private:
    Real stiffness_scale_;
};

} // namespace SPH

#endif // SPDB_CONTACT_SCALED_H
