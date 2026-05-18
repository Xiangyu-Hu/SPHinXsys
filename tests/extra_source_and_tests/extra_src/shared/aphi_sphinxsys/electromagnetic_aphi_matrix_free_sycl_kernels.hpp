#ifndef ELECTROMAGNETIC_APHI_MATRIX_FREE_SYCL_KERNELS_HPP
#define ELECTROMAGNETIC_APHI_MATRIX_FREE_SYCL_KERNELS_HPP

#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_residuals.h"
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_types.h"

#if SPHINXSYS_USE_SYCL
#include "aphi_sphinxsys/electromagnetic_aphi_matrix_free_sycl_queue.hpp"
#include "implementation_sycl.h"
#include <memory>
#include <sycl/sycl.hpp>
#endif

namespace SPH
{
namespace electromagnetics
{
namespace matrix_free
{
#if SPHINXSYS_USE_SYCL
namespace
{
inline Real hypot_complex(Real re, Real im)
{
    return sycl::sqrt(re * re + im * im);
}

/** (ar + i ai) / (br + i bi) */
inline void complex_div(Real ar, Real ai, Real br, Real bi, Real &out_re, Real &out_im)
{
    const Real denom = br * br + bi * bi;
    if (denom <= 0.0)
    {
        out_re = 0.0;
        out_im = 0.0;
        return;
    }
    out_re = (ar * br + ai * bi) / denom;
    out_im = (ai * br - ar * bi) / denom;
}

struct MatrixFreeSyclJacobiValueBuffers
{
    size_t size_ = 0;
    std::unique_ptr<sycl::buffer<Real, 1>> field_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> field_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> residual_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> residual_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> reaction_real_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> reaction_imag_buffer_;
    std::unique_ptr<sycl::buffer<Real, 1>> diagonal_scale_buffer_;
    const Complex *reaction_source_ = nullptr;
    size_t reaction_source_size_ = 0;
    size_t reaction_upload_count_ = 0;
    Real *dev_field_real_ = nullptr;
    Real *dev_field_imag_ = nullptr;

    void releaseFieldUsm()
    {
        if (dev_field_real_)
        {
            freeDeviceData(dev_field_real_);
            dev_field_real_ = nullptr;
        }
        if (dev_field_imag_)
        {
            freeDeviceData(dev_field_imag_);
            dev_field_imag_ = nullptr;
        }
    }

    bool deviceFieldUsmReady() const
    {
        return dev_field_real_ != nullptr && matrixFreeSyclFieldValuesUseDeviceUsm();
    }

    void syncHostFieldBuffersFromDeviceUsm()
    {
        if (!deviceFieldUsmReady())
        {
            return;
        }
        StdVec<Real> re(size_), im(size_);
        copyFromDevice(re.data(), dev_field_real_, size_);
        copyFromDevice(im.data(), dev_field_imag_, size_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> fr(*field_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> fi(*field_imag_buffer_);
        for (size_t i = 0; i != size_; ++i)
        {
            fr[i] = re[i];
            fi[i] = im[i];
        }
    }

    Real *mutableDeviceFieldRealUsmForJacobi()
    {
        return deviceFieldUsmReady() ? dev_field_real_ : nullptr;
    }
    Real *mutableDeviceFieldImagUsmForJacobi()
    {
        return deviceFieldUsmReady() ? dev_field_imag_ : nullptr;
    }

    void resize(size_t size)
    {
        if (size_ == size && field_real_buffer_ && field_imag_buffer_ && residual_real_buffer_ &&
            residual_imag_buffer_ && reaction_real_buffer_ && reaction_imag_buffer_ && diagonal_scale_buffer_)
        {
            return;
        }

        releaseFieldUsm();
        size_ = size;
        field_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        field_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        residual_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        residual_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        reaction_real_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        reaction_imag_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        diagonal_scale_buffer_ = std::make_unique<sycl::buffer<Real, 1>>(sycl::range<1>(size));
        reaction_source_ = nullptr;
        reaction_source_size_ = 0;
        if (matrixFreeSyclFieldValuesUseDeviceUsm())
        {
            dev_field_real_ = allocateDeviceOnly<Real>(size);
            dev_field_imag_ = allocateDeviceOnly<Real>(size);
        }
    }

    void loadInputs(const StdVec<Complex> &field, const StdVec<Complex> &reaction_coefficient,
                    const ScalarComplexHelmholtzResiduals &residuals)
    {
        sycl::host_accessor<Real, 1, sycl::access::mode::write> field_real(*field_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> field_imag(*field_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> residual_real(*residual_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> residual_imag(*residual_imag_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::write> diagonal_scale(*diagonal_scale_buffer_);

        for (size_t i = 0; i != field.size(); ++i)
        {
            field_real[i] = field[i].real();
            field_imag[i] = field[i].imag();
            residual_real[i] = residuals.residual_[i].real();
            residual_imag[i] = residuals.residual_[i].imag();
            diagonal_scale[i] = residuals.diagonal_scale_[i];
        }

        if (reaction_source_ != reaction_coefficient.data() || reaction_source_size_ != reaction_coefficient.size())
        {
            sycl::host_accessor<Real, 1, sycl::access::mode::write> reaction_real(*reaction_real_buffer_);
            sycl::host_accessor<Real, 1, sycl::access::mode::write> reaction_imag(*reaction_imag_buffer_);
            for (size_t i = 0; i != reaction_coefficient.size(); ++i)
            {
                reaction_real[i] = reaction_coefficient[i].real();
                reaction_imag[i] = reaction_coefficient[i].imag();
            }
            reaction_source_ = reaction_coefficient.data();
            reaction_source_size_ = reaction_coefficient.size();
            ++reaction_upload_count_;
        }

        if (deviceFieldUsmReady())
        {
            StdVec<Real> re(field.size()), im(field.size());
            for (size_t i = 0; i != field.size(); ++i)
            {
                re[i] = field_real[i];
                im[i] = field_imag[i];
            }
            copyToDevice(re.data(), dev_field_real_, field.size());
            copyToDevice(im.data(), dev_field_imag_, field.size());
        }
    }

    void readField(StdVec<Complex> &field)
    {
        syncHostFieldBuffersFromDeviceUsm();
        sycl::host_accessor<Real, 1, sycl::access::mode::read> field_real(*field_real_buffer_);
        sycl::host_accessor<Real, 1, sycl::access::mode::read> field_imag(*field_imag_buffer_);
        for (size_t i = 0; i != field.size(); ++i)
        {
            field[i] = Complex(field_real[i], field_imag[i]);
        }
    }
};

inline MatrixFreeSyclJacobiValueBuffers &matrixFreeSyclJacobiValueBuffers()
{
    static MatrixFreeSyclJacobiValueBuffers value_buffers;
    return value_buffers;
}
} // namespace

inline size_t matrixFreeJacobiReactionUploadCount()
{
    return matrixFreeSyclJacobiValueBuffers().reaction_upload_count_;
}

/**
 * @brief Device Jacobi step for scalar complex Helmholtz (matches CPU applyScalarComplexHelmholtzJacobiUpdate).
 * @details Reuses persistent SYCL buffers; host/device synchronization remains explicit at this prototype stage.
 */
inline void applyScalarComplexHelmholtzJacobiUpdateSyclWithWorkspace(
    StdVec<Complex> &field, const StdVec<Complex> &reaction_coefficient,
    const ScalarComplexHelmholtzResiduals &residuals, Real relaxation_factor, Real diagonal_regularization,
    MatrixFreeSyclJacobiValueBuffers &value_buffers)
{
    const size_t n = field.size();
    if (n == 0 || reaction_coefficient.size() != n || residuals.residual_.size() != n ||
        residuals.diagonal_scale_.size() != n)
    {
        return;
    }

    value_buffers.resize(n);
    value_buffers.loadInputs(field, reaction_coefficient, residuals);

    Real *field_re_usm = value_buffers.mutableDeviceFieldRealUsmForJacobi();
    Real *field_im_usm = value_buffers.mutableDeviceFieldImagUsmForJacobi();
    if (field_re_usm != nullptr && field_im_usm != nullptr)
    {
        matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
            auto a_rr = value_buffers.residual_real_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto a_ri = value_buffers.residual_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto a_rcr = value_buffers.reaction_real_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto a_rci = value_buffers.reaction_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto a_dr = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::read>(cgh);

            const Real relax = relaxation_factor;
            const Real reg = diagonal_regularization;
            Real *const f_re = field_re_usm;
            Real *const f_im = field_im_usm;

            cgh.parallel_for(sycl::range<1>(n), [=](sycl::id<1> id) {
                const size_t i = id[0];
                const Real dscale = a_dr[i];
                const Real ld_re = dscale + reg + a_rcr[i];
                const Real ld_im = a_rci[i];
                if (hypot_complex(ld_re, ld_im) <= reg)
                {
                    return;
                }
                Real inv_rr = 0.0, inv_ii = 0.0;
                complex_div(a_rr[i], a_ri[i], ld_re, ld_im, inv_rr, inv_ii);
                f_re[i] -= relax * inv_rr;
                f_im[i] -= relax * inv_ii;
            });
        });
    }
    else
    {
        matrixFreeSyclSubmitAndWait([&](sycl::handler &cgh) {
            auto a_fr = value_buffers.field_real_buffer_->get_access<sycl::access::mode::read_write>(cgh);
            auto a_fi = value_buffers.field_imag_buffer_->get_access<sycl::access::mode::read_write>(cgh);
            auto a_rr = value_buffers.residual_real_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto a_ri = value_buffers.residual_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto a_rcr = value_buffers.reaction_real_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto a_rci = value_buffers.reaction_imag_buffer_->get_access<sycl::access::mode::read>(cgh);
            auto a_dr = value_buffers.diagonal_scale_buffer_->get_access<sycl::access::mode::read>(cgh);

            const Real relax = relaxation_factor;
            const Real reg = diagonal_regularization;

            cgh.parallel_for(sycl::range<1>(n), [=](sycl::id<1> id) {
                const size_t i = id[0];
                const Real dscale = a_dr[i];
                const Real ld_re = dscale + reg + a_rcr[i];
                const Real ld_im = a_rci[i];
                if (hypot_complex(ld_re, ld_im) <= reg)
                {
                    return;
                }
                Real inv_rr = 0.0, inv_ii = 0.0;
                complex_div(a_rr[i], a_ri[i], ld_re, ld_im, inv_rr, inv_ii);
                a_fr[i] -= relax * inv_rr;
                a_fi[i] -= relax * inv_ii;
            });
        });
    }

    value_buffers.readField(field);
}

inline void applyScalarComplexHelmholtzJacobiUpdateSycl(StdVec<Complex> &field,
                                                       const StdVec<Complex> &reaction_coefficient,
                                                       const ScalarComplexHelmholtzResiduals &residuals,
                                                       Real relaxation_factor, Real diagonal_regularization)
{
    applyScalarComplexHelmholtzJacobiUpdateSyclWithWorkspace(field, reaction_coefficient, residuals, relaxation_factor,
                                                            diagonal_regularization,
                                                            matrixFreeSyclJacobiValueBuffers());
}
#endif // SPHINXSYS_USE_SYCL

} // namespace matrix_free
} // namespace electromagnetics
} // namespace SPH

#endif // ELECTROMAGNETIC_APHI_MATRIX_FREE_SYCL_KERNELS_HPP
