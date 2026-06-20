#ifndef ELECTROMAGNETIC_OPHELIE_PHI_DEVICE_VECTOR_OPS_H
#define ELECTROMAGNETIC_OPHELIE_PHI_DEVICE_VECTOR_OPS_H

#include "electromagnetic_ophelie_observables.h"
#include "electromagnetic_ophelie_phi_krylov_ck.h"
#include "implementation_sycl.h"
#include "sphinxsys.h"

#include <cmath>
#include <memory>
#include <vector>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

#if SPHINXSYS_USE_SYCL

inline void deviceSubtractScaledUSM(Real *values, const Real *direction, Real scale, size_t n)
{
    if (n == 0)
    {
        return;
    }
    auto &sycl_queue = execution_instance.getQueue();
    sycl_queue
        .submit([&](sycl::handler &cgh)
                { cgh.parallel_for(execution_instance.getUniformNdRange(n),
                                   [=](sycl::nd_item<1> item)
                                   {
                                       const size_t i = item.get_global_id(0);
                                       if (i < n)
                                       {
                                           values[i] -= scale * direction[i];
                                       }
                                   }); })
        .wait_and_throw();
}

inline void deviceScaleUSM(Real *values, Real scale, size_t n)
{
    if (n == 0)
    {
        return;
    }
    auto &sycl_queue = execution_instance.getQueue();
    sycl_queue
        .submit([&](sycl::handler &cgh)
                { cgh.parallel_for(execution_instance.getUniformNdRange(n),
                                   [=](sycl::nd_item<1> item)
                                   {
                                       const size_t i = item.get_global_id(0);
                                       if (i < n)
                                       {
                                           values[i] *= scale;
                                       }
                                   }); })
        .wait_and_throw();
}

inline void deviceCopyUSM(Real *dst, const Real *src, size_t n)
{
    if (n == 0 || dst == src)
    {
        return;
    }
    auto &sycl_queue = execution_instance.getQueue();
    sycl_queue.memcpy(dst, src, n * sizeof(Real)).wait_and_throw();
}

/**
 * Krylov basis/workspace in USM shared; dot/norm via SphinxSys particle_reduce on staged particle fields;
 * subtract/scale/copy on USM (shared, validated on this backend).
 */
class OpheliePhiDeviceKrylovStorage
{
  public:
    ~OpheliePhiDeviceKrylovStorage() { release(); }

    template <class ExecutionPolicy>
    void setup(SPHBody &glass_body, const OphelieGlassFieldNames &names, size_t n, UnsignedInt restart_dimension)
    {
        (void)ExecutionPolicy{};
        release();
        glass_body_ = &glass_body;
        names_ = names;
        n_ = n;
        restart_dimension_ = restart_dimension;
        workspace_shared_ = allocateDeviceShared<Real>(n_);
        basis_shared_.resize(restart_dimension_ + 1, nullptr);
        for (UnsignedInt j = 0; j <= restart_dimension_; ++j)
        {
            basis_shared_[j] = allocateDeviceShared<Real>(n_);
        }
    }

    void writeBasisFromHost(UnsignedInt index, const Real *host_values)
    {
        copyToDevice(host_values, basis_shared_[index], n_);
    }

    void readBasisToHost(UnsignedInt index, Real *host_values)
    {
        copyFromDevice(host_values, basis_shared_[index], n_);
    }

    void writeWorkspaceFromHost(const Real *host_values)
    {
        copyToDevice(host_values, workspace_shared_, n_);
    }

    Real dotBasisWithWorkspace(UnsignedInt basis_index) const
    {
        BaseParticles &particles = glass_body_->getBaseParticles();
        hostAssignScalarField(particles, names_.krylov_basis, basis_shared_[basis_index], n_);
        hostAssignScalarField(particles, names_.krylov_workspace, workspace_shared_, n_);
        syncVariableToDevice<Real>(particles, names_.krylov_basis);
        syncVariableToDevice<Real>(particles, names_.krylov_workspace);
        return ophelieVolWeightedDotParticleFields(*glass_body_, names_.krylov_basis, names_.krylov_workspace);
    }

    Real normWorkspace() const
    {
        BaseParticles &particles = glass_body_->getBaseParticles();
        hostAssignScalarField(particles, names_.krylov_workspace, workspace_shared_, n_);
        syncVariableToDevice<Real>(particles, names_.krylov_workspace);
        return ophelieVolWeightedNormParticleField(*glass_body_, names_.krylov_workspace);
    }

    void subtractBasisFromWorkspace(UnsignedInt basis_index, Real scale)
    {
        deviceSubtractScaledUSM(workspace_shared_, basis_shared_[basis_index], scale, n_);
    }

    void normalizeWorkspaceIntoBasis(UnsignedInt target_basis_index, Real inverse_norm)
    {
        deviceScaleUSM(workspace_shared_, inverse_norm, n_);
        deviceCopyUSM(basis_shared_[target_basis_index], workspace_shared_, n_);
    }

    size_t size() const { return n_; }

  private:
    SPHBody *glass_body_ = nullptr;
    OphelieGlassFieldNames names_;
    size_t n_ = 0;
    UnsignedInt restart_dimension_ = 0;
    Real *workspace_shared_ = nullptr;
    StdVec<Real *> basis_shared_;

    void release()
    {
        if (workspace_shared_ != nullptr)
        {
            freeDeviceData(workspace_shared_);
            workspace_shared_ = nullptr;
        }
        for (Real *ptr : basis_shared_)
        {
            if (ptr != nullptr)
            {
                freeDeviceData(ptr);
            }
        }
        basis_shared_.clear();
        glass_body_ = nullptr;
        n_ = 0;
        restart_dimension_ = 0;
    }
};

struct OpheliePhiDeviceKrylovArnoldiCheck
{
    Real host_subdiagonal_norm = 0.0;
    Real device_subdiagonal_norm = 0.0;
    Real norm_rel_err = 0.0;
    bool passed = false;
};

template <class ExecutionPolicy>
inline OpheliePhiDeviceKrylovArnoldiCheck checkHostDeviceKrylovArnoldiStep(
    SPHBody &glass_body, const OphelieGlassFieldNames &names, const StdVec<const Real *> &basis_host,
    const Real *workspace_host, size_t inner_step, size_t n)
{
    OpheliePhiDeviceKrylovArnoldiCheck check;
    StdVec<Real> workspace_copy(workspace_host, workspace_host + n);
    const UnsignedInt restart_dimension = static_cast<UnsignedInt>(basis_host.size() - 1);
    BaseParticles &particles = glass_body.getBaseParticles();

    OpheliePhiDeviceKrylovStorage device_krylov;
    device_krylov.setup<ExecutionPolicy>(glass_body, names, n, restart_dimension);
    for (UnsignedInt j = 0; j <= inner_step; ++j)
    {
        device_krylov.writeBasisFromHost(j, basis_host[j]);
    }
    device_krylov.writeWorkspaceFromHost(workspace_host);

    for (UnsignedInt i = 0; i <= inner_step; ++i)
    {
        const Real projection = hostVolWeightedDot(particles, basis_host[i], workspace_copy.data(), n);
        hostSubtractScaledVector(workspace_copy.data(), basis_host[i], projection, n);
    }
    check.host_subdiagonal_norm = hostVolWeightedNorm(particles, workspace_copy.data(), n);

    for (UnsignedInt i = 0; i <= inner_step; ++i)
    {
        const Real projection = device_krylov.dotBasisWithWorkspace(i);
        device_krylov.subtractBasisFromWorkspace(i, projection);
    }
    check.device_subdiagonal_norm = device_krylov.normWorkspace();

    const Real denom = std::max(check.host_subdiagonal_norm, Real(1));
    check.norm_rel_err = std::abs(check.host_subdiagonal_norm - check.device_subdiagonal_norm) / denom;
    check.passed = std::isfinite(check.device_subdiagonal_norm) && check.norm_rel_err < Real(1e-5);
    return check;
}

class OpheliePhiDeviceVectorWorkspace
{
  public:
    void bind(SPHBody &glass_body, const OphelieGlassFieldNames &names, size_t n)
    {
        glass_body_ = &glass_body;
        names_ = names;
        size_ = n;
    }

    Real volWeightedDot(const Real *lhs_host, const Real *rhs_host) const
    {
        ophelieStageHostArrayToParticleField(glass_body_->getBaseParticles(), names_.krylov_scratch_a, lhs_host, size_);
        ophelieStageHostArrayToParticleField(glass_body_->getBaseParticles(), names_.krylov_scratch_b, rhs_host, size_);
        return ophelieVolWeightedDotParticleFields(*glass_body_, names_.krylov_scratch_a, names_.krylov_scratch_b);
    }

    Real volWeightedNorm(const Real *values_host) const
    {
        ophelieStageHostArrayToParticleField(glass_body_->getBaseParticles(), names_.krylov_scratch_a, values_host,
                                             size_);
        return ophelieVolWeightedNormParticleField(*glass_body_, names_.krylov_scratch_a);
    }

    template <class ExecutionPolicy>
    void volWeightedAxpy(Real alpha, const Real *x_host, Real *y_host) const
    {
        ophelieStageHostArrayToParticleField(glass_body_->getBaseParticles(), names_.krylov_scratch_a, x_host, size_);
        ophelieStageHostArrayToParticleField(glass_body_->getBaseParticles(), names_.krylov_scratch_b, y_host, size_);
        ophelieKrylovAxpyFields<ExecutionPolicy>(*glass_body_, alpha, names_.krylov_scratch_a, names_.krylov_scratch_b);
        ophelieReadParticleFieldToHostBuffer(glass_body_->getBaseParticles(), names_.krylov_scratch_b, y_host, size_);
    }

    size_t size() const { return size_; }

  private:
    SPHBody *glass_body_ = nullptr;
    OphelieGlassFieldNames names_;
    size_t size_ = 0;
};

#endif // SPHINXSYS_USE_SYCL

struct OpheliePhiDeviceVectorOpsCheck
{
    Real host_dot = 0.0;
    Real host_norm = 0.0;
    Real device_dot = 0.0;
    Real device_norm = 0.0;
    Real dot_rel_err = 0.0;
    Real norm_rel_err = 0.0;
    Real axpy_max_abs_err = 0.0;
    bool passed = false;
};

template <class ExecutionPolicy>
inline OpheliePhiDeviceVectorOpsCheck checkHostDeviceVolWeightedVectorOps(SPHBody &glass_body,
                                                                          const OphelieGlassFieldNames &names,
                                                                          BaseParticles &particles, const Real *lhs,
                                                                          const Real *rhs, size_t n,
                                                                          Real axpy_alpha = Real(0.25))
{
    OpheliePhiDeviceVectorOpsCheck check;
    check.host_dot = hostVolWeightedDot(particles, lhs, rhs, n);
    check.host_norm = hostVolWeightedNorm(particles, lhs, n);

#if SPHINXSYS_USE_SYCL
    OpheliePhiDeviceVectorWorkspace workspace;
    workspace.bind(glass_body, names, n);
    check.device_dot = workspace.volWeightedDot(lhs, rhs);
    check.device_norm = workspace.volWeightedNorm(lhs);

    StdVec<Real> y_host(lhs, lhs + n);
    StdVec<Real> y_ref(lhs, lhs + n);
    for (size_t i = 0; i < n; ++i)
    {
        y_ref[i] += axpy_alpha * rhs[i];
    }
    workspace.volWeightedAxpy<ExecutionPolicy>(axpy_alpha, rhs, y_host.data());
    Real axpy_max_err = 0.0;
    for (size_t i = 0; i < n; ++i)
    {
        axpy_max_err = std::max(axpy_max_err, std::abs(y_host[i] - y_ref[i]));
    }
    check.axpy_max_abs_err = axpy_max_err;

    const Real dot_denom = std::max(std::abs(check.host_dot), Real(1));
    const Real norm_denom = std::max(check.host_norm, Real(1));
    check.dot_rel_err = std::abs(check.host_dot - check.device_dot) / dot_denom;
    check.norm_rel_err = std::abs(check.host_norm - check.device_norm) / norm_denom;
    check.passed = check.dot_rel_err < Real(1e-5) && check.norm_rel_err < Real(1e-5) && check.axpy_max_abs_err < Real(1e-5);
#else
    check.device_dot = check.host_dot;
    check.device_norm = check.host_norm;
    check.passed = true;
    (void)glass_body;
    (void)names;
    (void)axpy_alpha;
#endif
    return check;
}

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_DEVICE_VECTOR_OPS_H
