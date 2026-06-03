#ifndef ELECTROMAGNETIC_OPHELIE_SOURCE_HPP
#define ELECTROMAGNETIC_OPHELIE_SOURCE_HPP

#include "electromagnetic_ophelie_source.h"

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

inline InitializeOphelieCoilSourceCK::InitializeOphelieCoilSourceCK(SPHBody &sph_body, const OphelieCoilFieldNames &names,
                                                                  const OphelieParameters &params,
                                                                  const Vecd &current_center)
    : LocalDynamics(sph_body),
      j0_(params.coil_j0_override_ > 0.0 ? params.coil_j0_override_ : params.equivalentCurrentDensity()),
      current_center_(current_center), outer_shell_only_(params.coil_j_outer_shell_only_),
      outer_shell_radius_fraction_(params.coil_j_outer_shell_radius_fraction_),
      coil_max_xy_radius_(params.coil_max_xy_radius_m_),
      dv_position_(particles_->template getVariableByName<Vecd>("Position")),
      dv_j_src_real_(particles_->template getVariableByName<Vecd>(names.j_src_real)),
      dv_j_src_imag_(particles_->template getVariableByName<Vecd>(names.j_src_imag))
{
}

template <class ExecutionPolicy, class EncloserType>
inline InitializeOphelieCoilSourceCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : j0_(encloser.j0_), current_center_(encloser.current_center_),
      outer_shell_only_(encloser.outer_shell_only_),
      outer_shell_radius_fraction_(encloser.outer_shell_radius_fraction_),
      coil_max_xy_radius_(encloser.coil_max_xy_radius_),
      position_(encloser.dv_position_->DelegatedData(ex_policy)),
      j_src_real_(encloser.dv_j_src_real_->DelegatedData(ex_policy)),
      j_src_imag_(encloser.dv_j_src_imag_->DelegatedData(ex_policy))
{
}

inline void InitializeOphelieCoilSourceCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    const Real dx = position_[index_i][0] - current_center_[0];
    const Real dy = position_[index_i][1] - current_center_[1];
    const Real r = std::sqrt(dx * dx + dy * dy) + TinyReal;
    Real j_mag = j0_;
    if (outer_shell_only_ && coil_max_xy_radius_ > TinyReal)
    {
        if (r < outer_shell_radius_fraction_ * coil_max_xy_radius_)
        {
            j_mag = 0.0;
        }
    }
    const Vecd e_theta(-dy / r, dx / r, 0.0);
    j_src_real_[index_i] = j_mag * e_theta;
    j_src_imag_[index_i] = Vecd::Zero();
}

inline AssignOphelieGlassSigmaCK::AssignOphelieGlassSigmaCK(SPHBody &sph_body, const OphelieGlassFieldNames &names,
                                                            Real sigma)
    : LocalDynamics(sph_body), sigma_(sigma), dv_sigma_(particles_->template getVariableByName<Real>(names.sigma))
{
}

template <class ExecutionPolicy, class EncloserType>
inline AssignOphelieGlassSigmaCK::UpdateKernel::UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
    : sigma_(encloser.sigma_), sigma_field_(encloser.dv_sigma_->DelegatedData(ex_policy))
{
}

inline void AssignOphelieGlassSigmaCK::UpdateKernel::update(size_t index_i, Real dt)
{
    (void)dt;
    sigma_field_[index_i] = sigma_;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_SOURCE_HPP
