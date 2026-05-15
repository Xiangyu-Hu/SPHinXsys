#ifndef ACTIVE_MODEL_HPP
#define ACTIVE_MODEL_HPP

#include "active_model.h"

namespace SPH
{
//=================================================================================================//
template <typename ExecutionPolicy>
ActiveModelSolid::ConstituteKernel::ConstituteKernel(
    const ExecutionPolicy &ex_policy, ActiveModelSolid &encloser)
    : SaintVenantKirchhoffSolid::ConstituteKernel(ex_policy, encloser),
      active_strain_(encloser.dv_active_strain_->DelegatedData(ex_policy)) {}
//=================================================================================================//
inline Matd ActiveModelSolid::ConstituteKernel::StressPK1(const Matd &F, size_t index_i)
{
    // GPU-safe analytical Cholesky for 2x2 symmetric PD matrix A = I + 2*active_strain.
    // Eigen's .llt() uses static panel-kernel variables forbidden in SYCL device code.
    Matd A = 2.0 * active_strain_[index_i] + Matd::Identity();
    Matd F0 = Matd::Zero();
    F0(0, 0) = math::sqrt(A(0, 0));
    F0(1, 0) = A(1, 0) / F0(0, 0);
    F0(1, 1) = math::sqrt(A(1, 1) - F0(1, 0) * F0(1, 0));

    Matd F0_inv = F0.inverse();
    Matd F_e = F * F0_inv;
    Matd F0_star = F0.determinant() * F0_inv.transpose();
    Matd E_e = 0.5 * (F.transpose() * F - Matd::Identity()) - active_strain_[index_i];
    return F_e * (lambda0_ * E_e.trace() * Matd::Identity() + 2.0 * G0_ * E_e) * F0_star;
}
//=================================================================================================//
} // namespace SPH
#endif // ACTIVE_MODEL_HPP
