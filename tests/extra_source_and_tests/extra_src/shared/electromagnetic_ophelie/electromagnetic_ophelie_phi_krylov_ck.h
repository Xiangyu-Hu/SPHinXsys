#ifndef ELECTROMAGNETIC_OPHELIE_PHI_KRYLOV_CK_H
#define ELECTROMAGNETIC_OPHELIE_PHI_KRYLOV_CK_H

#include "electromagnetic_ophelie_phi.h"
#include "electromagnetic_ophelie_device_sync.h"
#include "electromagnetic_ophelie_field_names.h"
#include "electromagnetic_ophelie_observables.h"
#include "loop_range.h"
#include "particle_iterators_sycl.h"
#include "reduce_functors.h"
#include "simple_algorithms_ck.h"
#include "sphinxsys.h"

#include <cmath>
#include <string>

namespace SPH
{
namespace electromagnetics
{
namespace ophelie
{

#if SPHINXSYS_USE_SYCL

using OphelieKrylovDevicePolicy = ParallelDevicePolicy;

/** y <- alpha * x + y on particle scalar fields (StateDynamics + CK). */
class OphelieKrylovAxpyFieldsCK : public LocalDynamics
{
  public:
    OphelieKrylovAxpyFieldsCK(SPHBody &sph_body, const std::string &y_field, const std::string &x_field, Real alpha)
        : LocalDynamics(sph_body), alpha_(alpha), dv_y_(particles_->template getVariableByName<Real>(y_field)),
          dv_x_(particles_->template getVariableByName<Real>(x_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : alpha_(encloser.alpha_), y_(encloser.dv_y_->DelegatedData(ex_policy)),
              x_(encloser.dv_x_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            y_[index_i] += alpha_ * x_[index_i];
        }

      protected:
        Real alpha_;
        Real *y_;
        Real *x_;
    };

  protected:
    Real alpha_;
    DiscreteVariable<Real> *dv_y_;
    DiscreteVariable<Real> *dv_x_;
};

template <class ExecutionPolicy>
inline void ophelieKrylovAxpyFields(SPHBody &sph_body, Real alpha, const std::string &x_field, const std::string &y_field)
{
    StateDynamics<ExecutionPolicy, OphelieKrylovAxpyFieldsCK> axpy(sph_body, y_field, x_field, alpha);
    axpy.exec();
}

/** workspace[index] -= scale * basis[index] on device (StateDynamics + CK). */
class OphelieKrylovSubtractScaledFieldsCK : public LocalDynamics
{
  public:
    OphelieKrylovSubtractScaledFieldsCK(SPHBody &sph_body, const std::string &workspace_field,
                                        const std::string &basis_field, Real scale)
        : LocalDynamics(sph_body), scale_(scale),
          dv_workspace_(particles_->template getVariableByName<Real>(workspace_field)),
          dv_basis_(particles_->template getVariableByName<Real>(basis_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : scale_(encloser.scale_), workspace_(encloser.dv_workspace_->DelegatedData(ex_policy)),
              basis_(encloser.dv_basis_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            workspace_[index_i] -= scale_ * basis_[index_i];
        }

      protected:
        Real scale_;
        Real *workspace_;
        Real *basis_;
    };

  protected:
    Real scale_;
    DiscreteVariable<Real> *dv_workspace_;
    DiscreteVariable<Real> *dv_basis_;
};

/** workspace[index] = src[index] on device. */
class OphelieKrylovCopyFieldCK : public LocalDynamics
{
  public:
    OphelieKrylovCopyFieldCK(SPHBody &sph_body, const std::string &dst_field, const std::string &src_field)
        : LocalDynamics(sph_body), dv_dst_(particles_->template getVariableByName<Real>(dst_field)),
          dv_src_(particles_->template getVariableByName<Real>(src_field))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : dst_(encloser.dv_dst_->DelegatedData(ex_policy)), src_(encloser.dv_src_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            dst_[index_i] = src_[index_i];
        }

      protected:
        Real *dst_;
        Real *src_;
    };

  protected:
    DiscreteVariable<Real> *dv_dst_;
    DiscreteVariable<Real> *dv_src_;
};

inline Real *ophelieDelegatedRealField(BaseParticles &particles, const std::string &field_name)
{
    return particles.template getVariableByName<Real>(field_name)->DelegatedData(OphelieKrylovDevicePolicy{});
}

inline void ophelieStageHostArrayToParticleField(BaseParticles &particles, const std::string &field_name,
                                               const Real *host_values, size_t n)
{
    hostAssignScalarField(particles, field_name, host_values, n);
    syncVariableToDevice<Real>(particles, field_name);
}

inline void ophelieStageSharedBufferToParticleField(BaseParticles &particles, const std::string &field_name,
                                                    const Real *shared_buffer, size_t n)
{
    Real *device_field = ophelieDelegatedRealField(particles, field_name);
    copyToDevice(shared_buffer, device_field, n);
}

inline void ophelieReadParticleFieldToHostBuffer(BaseParticles &particles, const std::string &field_name, Real *buffer,
                                                 size_t n)
{
    syncVariableToHost<Real>(particles, field_name);
    hostReadScalarField(particles, field_name, buffer, n);
}

/** Volume-weighted dot via particle_reduce + LoopRangeCK (SphinxSys SYCL native). */
inline Real ophelieVolWeightedDotParticleFields(SPHBody &sph_body, const std::string &lhs_field,
                                                const std::string &rhs_field)
{
    BaseParticles &particles = sph_body.getBaseParticles();
    syncVariableToDevice<Real>(particles, "VolumetricMeasure");
    syncVariableToDevice<Real>(particles, lhs_field);
    syncVariableToDevice<Real>(particles, rhs_field);

    Real *vol = ophelieDelegatedRealField(particles, "VolumetricMeasure");
    Real *lhs = ophelieDelegatedRealField(particles, lhs_field);
    Real *rhs = ophelieDelegatedRealField(particles, rhs_field);

    return particle_reduce<ReduceSum<Real>>(
        LoopRangeCK<OphelieKrylovDevicePolicy, SPHBody>(sph_body), ReduceReference<ReduceSum<Real>>::value,
        [=](size_t index_i) { return vol[index_i] * lhs[index_i] * rhs[index_i]; });
}

inline Real ophelieVolWeightedNormParticleField(SPHBody &sph_body, const std::string &values_field)
{
    const Real dot_value = ophelieVolWeightedDotParticleFields(sph_body, values_field, values_field);
    return std::sqrt(std::max(dot_value, Real(0)));
}

template <class ExecutionPolicy>
inline void ophelieKrylovSubtractScaledFields(SPHBody &sph_body, const std::string &workspace_field,
                                              const std::string &basis_field, Real scale)
{
    StateDynamics<ExecutionPolicy, OphelieKrylovSubtractScaledFieldsCK> subtract(sph_body, workspace_field, basis_field,
                                                                                scale);
    subtract.exec();
}

template <class ExecutionPolicy>
inline void ophelieKrylovScaleField(SPHBody &sph_body, const std::string &field_name, Real scale_factor)
{
    StateDynamics<ExecutionPolicy, OphelieScaleScalarFieldCK> scale(sph_body, field_name, scale_factor);
    scale.exec();
}

template <class ExecutionPolicy>
inline void ophelieKrylovCopyField(SPHBody &sph_body, const std::string &dst_field, const std::string &src_field)
{
    StateDynamics<ExecutionPolicy, OphelieKrylovCopyFieldCK> copy(sph_body, dst_field, src_field);
    copy.exec();
}

#endif // SPHINXSYS_USE_SYCL

} // namespace ophelie
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_OPHELIE_PHI_KRYLOV_CK_H
