#include "base_kernel.h"

namespace SPH
{
//=================================================================================================//
Kernel::Kernel(Real h, Real kernel_size, Real truncation, const std::string &name)
    : kernel_name_(name), h_(h), inv_h_(1.0 / h), kernel_size_(kernel_size),
      truncation_(truncation), rc_ref_(truncation * h), rc_ref_sqr_(rc_ref_ * rc_ref_),
      h_factor_W_1D_(std::bind(&Kernel::factorW1D, this, _1)),
      h_factor_W_2D_(std::bind(&Kernel::factorW2D, this, _1)),
      h_factor_W_3D_(std::bind(&Kernel::factorW3D, this, _1)),
      h_factor_dW_1D_(std::bind(&Kernel::factordW1D, this, _1)),
      h_factor_dW_2D_(std::bind(&Kernel::factordW2D, this, _1)),
      h_factor_dW_3D_(std::bind(&Kernel::factordW3D, this, _1)),
      h_factor_d2W_1D_(std::bind(&Kernel::factord2W1D, this, _1)),
      h_factor_d2W_2D_(std::bind(&Kernel::factord2W2D, this, _1)),
      h_factor_d2W_3D_(std::bind(&Kernel::factord2W3D, this, _1)){};
//=================================================================================================//
void Kernel::setDerivativeParameters()
{
    factor_dW_1D_ = inv_h_ * factor_W_1D_;
    factor_dW_2D_ = inv_h_ * factor_W_2D_;
    factor_dW_3D_ = inv_h_ * factor_W_3D_;
    factor_d2W_1D_ = inv_h_ * factor_dW_1D_;
    factor_d2W_2D_ = inv_h_ * factor_dW_2D_;
    factor_d2W_3D_ = inv_h_ * factor_dW_3D_;
}
//=================================================================================================//
void Kernel::resetSmoothingLength(Real h)
{
    Real ratio = h_ / h;
    h_ /= ratio;
    inv_h_ *= ratio;
    rc_ref_ /= ratio;
    rc_ref_sqr_ /= ratio * ratio;

    factor_W_1D_ *= ratio;
    factor_W_2D_ *= ratio * ratio;
    factor_W_3D_ *= ratio * ratio * ratio;

    setDerivativeParameters();
}
//=================================================================================================//
Real Kernel::W(const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_W_1D_ * W_1D(q);
}
//=================================================================================================//
Real Kernel::W(const Real &r_ij, const Vec2d &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_W_2D_ * W_2D(q);
}
//=================================================================================================//
Real Kernel::W(const Real &r_ij, const Vec3d &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_W_3D_ * W_3D(q);
}
//=================================================================================================//
Real Kernel::dW(const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_dW_1D_ * dW_1D(q);
}
//=================================================================================================//
Real Kernel::dW(const Real &r_ij, const Vec2d &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_dW_2D_ * dW_2D(q);
}
//=================================================================================================//
Real Kernel::dW(const Real &r_ij, const Vec3d &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_dW_3D_ * dW_3D(q);
}
//=================================================================================================//
Real Kernel::d2W(const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_d2W_1D_ * d2W_1D(q);
}
//=================================================================================================//
Real Kernel::d2W(const Real &r_ij, const Vec2d &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_d2W_2D_ * d2W_2D(q);
}
//=================================================================================================//
Real Kernel::d2W(const Real &r_ij, const Vec3d &displacement) const
{
    Real q = r_ij * inv_h_;
    return factor_d2W_3D_ * d2W_3D(q);
}
//=================================================================================================//
Real Kernel::W(const Real &h_ratio, const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_W_1D_ * W_1D(q) * h_factor_W_1D_(h_ratio);
}
//=================================================================================================//
Real Kernel::W(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_W_2D_ * W_2D(q) * h_factor_W_2D_(h_ratio);
}
//=================================================================================================//
Real Kernel::W(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_W_3D_ * W_3D(q) * h_factor_W_3D_(h_ratio);
}
//=================================================================================================//
Real Kernel::W0(const Real &h_ratio, const Real &point_i) const
{
    return factor_W_1D_ * h_factor_W_1D_(h_ratio);
};
//=================================================================================================//
Real Kernel::W0(const Real &h_ratio, const Vec2d &point_i) const
{
    return factor_W_2D_ * h_factor_W_2D_(h_ratio);
};
//=================================================================================================//
Real Kernel::W0(const Real &h_ratio, const Vec3d &point_i) const
{
    return factor_W_3D_ * h_factor_W_3D_(h_ratio);
};
//=================================================================================================//
Real Kernel::dW(const Real &h_ratio, const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_dW_1D_ * dW_1D(q) * h_factor_dW_1D_(h_ratio);
}
//=================================================================================================//
Real Kernel::dW(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_dW_2D_ * dW_2D(q) * h_factor_dW_2D_(h_ratio);
}
//=================================================================================================//
Real Kernel::dW(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_dW_3D_ * dW_3D(q) * h_factor_dW_1D_(h_ratio);
}
//=================================================================================================//
Real Kernel::d2W(const Real &h_ratio, const Real &r_ij, const Real &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_d2W_1D_ * d2W_1D(q) * h_factor_d2W_1D_(h_ratio);
}
//=================================================================================================//
Real Kernel::d2W(const Real &h_ratio, const Real &r_ij, const Vec2d &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_d2W_2D_ * d2W_2D(q) * h_factor_d2W_2D_(h_ratio);
}
//=================================================================================================//
Real Kernel::d2W(const Real &h_ratio, const Real &r_ij, const Vec3d &displacement) const
{
    Real q = r_ij * inv_h_ * h_ratio;
    return factor_d2W_3D_ * d2W_3D(q) * h_factor_d2W_3D_(h_ratio);
}
//=================================================================================================//
void Kernel::reduceOnce()
{
    factor_W_3D_ = factor_W_2D_;
    factor_W_2D_ = factor_W_1D_;
    factor_W_1D_ = 0.0;
    setDerivativeParameters();

    h_factor_W_3D_ = std::bind(&Kernel::factorW2D, this, _1);
    h_factor_W_2D_ = std::bind(&Kernel::factorW1D, this, _1);
    h_factor_dW_3D_ = std::bind(&Kernel::factordW2D, this, _1);
    h_factor_dW_2D_ = std::bind(&Kernel::factordW1D, this, _1);
    h_factor_d2W_3D_ = std::bind(&Kernel::factord2W2D, this, _1);
    h_factor_d2W_2D_ = std::bind(&Kernel::factord2W1D, this, _1);
}
//=================================================================================================//
void Kernel::reduceTwice()
{
    factor_W_3D_ = factor_W_1D_;
    factor_W_2D_ = 0.0;
    factor_W_1D_ = 0.0;
    setDerivativeParameters();

    h_factor_W_3D_ = std::bind(&Kernel::factorW1D, this, _1);
    h_factor_dW_3D_ = std::bind(&Kernel::factordW1D, this, _1);
    h_factor_d2W_3D_ = std::bind(&Kernel::factord2W1D, this, _1);
}
//=================================================================================================//
} // namespace SPH
