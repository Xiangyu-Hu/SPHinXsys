#ifndef ELECTROMAGNETIC_OPHELIE_SOURCE_H
#define ELECTROMAGNETIC_OPHELIE_SOURCE_H

#include "base_general_dynamics.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_parameters.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

/** Azimuthal current density J = J0 e_theta on coil source particles (real phasor only). */
class InitializeOphelieCoilSourceCK : public LocalDynamics
{
  public:
    InitializeOphelieCoilSourceCK(SPHBody &sph_body, const OphelieCoilFieldNames &names, const OphelieParameters &params,
                                  const Vecd &current_center = Vecd::Zero());
    virtual ~InitializeOphelieCoilSourceCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real j0_;
        Vecd current_center_;
        bool outer_shell_only_;
        Real outer_shell_radius_fraction_;
        Real coil_max_xy_radius_;
        Vecd *position_;
        Vecd *j_src_real_;
        Vecd *j_src_imag_;
    };

  protected:
    Real j0_;
    Vecd current_center_;
    bool outer_shell_only_;
    Real outer_shell_radius_fraction_;
    Real coil_max_xy_radius_;
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_j_src_real_;
    DiscreteVariable<Vecd> *dv_j_src_imag_;
};

class AssignOphelieGlassSigmaCK : public LocalDynamics
{
  public:
    AssignOphelieGlassSigmaCK(SPHBody &sph_body, const OphelieGlassFieldNames &names, Real sigma);
    virtual ~AssignOphelieGlassSigmaCK() = default;

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);

        void update(size_t index_i, Real dt = 0.0);

      protected:
        Real sigma_;
        Real *sigma_field_;
    };

  protected:
    Real sigma_;
    DiscreteVariable<Real> *dv_sigma_;
};

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#include "electromagnetic_ophelie_source.hpp"
#endif // ELECTROMAGNETIC_OPHELIE_SOURCE_H
