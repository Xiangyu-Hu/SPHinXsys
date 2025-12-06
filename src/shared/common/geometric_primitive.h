/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	geometric_primitive.h
 * @brief 	TBD.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef GEOMETRIC_PRIMITIVE_H
#define GEOMETRIC_PRIMITIVE_H

#include <base_data_type.h>

namespace SPH
{
template <int N>
using VecdBound = Eigen::Matrix<Real, N, 1>;

template <int N>
using ArrayiBound = Eigen::Array<int, N, 1>;

template <template <int> typename BoundType, int N>
class BoundingBox
{
    using Vectype = BoundType<N>;

  public:
    Vectype lower_, upper_;

    BoundingBox() : lower_(Vectype::Zero()), upper_(Vectype::Zero()) {};
    BoundingBox(const Vectype &lower, const Vectype &upper)
        : lower_(lower), upper_(upper) {};
    BoundingBox(const Vectype &hlfsize) : lower_(-hlfsize), upper_(hlfsize) {};

    BoundingBox translate(const Vectype &translate) const
    {
        return BoundingBox(translate + lower_, translate + upper_);
    };

    bool checkContain(const Vectype &point) const
    {
        bool is_contain = true;
        for (int i = 0; i < N; ++i)
        {
            if (point[i] < lower_[i] || point[i] > upper_[i])
            {
                is_contain = false;
                break;
            }
        }
        return is_contain;
    };

    bool checkIntersect(const BoundingBox &another) const
    {
        bool is_intersect = false;
        for (int i = 0; i < N; ++i)
        {
            if (another.lower_[i] > upper_[i] || another.upper_[i] < lower_[i])
            {
                is_intersect = true;
                break;
            }
        }
        return is_intersect;
    };

    BoundingBox getIntersect(const BoundingBox &another) const
    {
        if (!checkIntersect(another))
        {
            std::cerr << "BoundingBox::getIntersect: no overlap!" << std::endl;
        }

        BoundingBox output(*this);
        for (int i = 0; i < N; ++i)
        {
            /** If the lower limit is inside change the lower limit. */
            if (lower_[i] < another.lower_[i] && another.lower_[i] < upper_[i])
                output.lower_[i] = another.lower_[i];
            /**  If the upper limit is inside, change the upper limit. */
            if (upper_[i] > another.upper_[i] && another.upper_[i] > lower_[i])
                output.upper_[i] = another.upper_[i];
        }
        return output;
    };

    static constexpr int DataSize() { return N; }

    Vectype BoundSize() const
    {
        return upper_ - lower_;
    };

    BoundingBox expand(const Vectype &expand_size) const
    {
        Vectype new_lower = lower_ - expand_size;
        Vectype new_upper = upper_ + expand_size;
        return BoundingBox(new_lower, new_upper);
    };

    auto MinimumDimension() const
    {
        return BoundSize().cwiseAbs().minCoeff();
    };
};
/** Operator define. */
template <template <int> typename BoundType, int N>
bool operator==(const BoundingBox<BoundType, N> &bb1, const BoundingBox<BoundType, N> &bb2)
{
    return bb1.lower_ == bb2.lower_ && bb1.upper_ == bb2.upper_ ? true : false;
};

using Rotation2d = Eigen::Rotation2D<Real>;
using Rotation3d = Eigen::AngleAxis<Real>;

template <typename RotationType, typename VecType>
class BaseTransform
{
  private:
    using MatType = decltype(RotationType().toRotationMatrix());
    MatType rotation_, inv_rotation_;
    VecType translation_;

  public:
    explicit BaseTransform(const RotationType &rotation, const VecType &translation = VecType::Zero())
        : rotation_(rotation.toRotationMatrix()), inv_rotation_(rotation_.transpose()), translation_(translation) {};
    explicit BaseTransform(const VecType &translation)
        : rotation_(MatType::Identity()), inv_rotation_(rotation_.transpose()), translation_(translation) {};
    BaseTransform() : BaseTransform(VecType::Zero()) {};

    /** Forward rotation. */
    VecType xformFrameVecToBase(const VecType &origin)
    {
        return rotation_ * origin;
    };

    /** Forward transformation. Note that the rotation operation is carried out first. */
    VecType shiftFrameStationToBase(const VecType &origin)
    {
        return translation_ + xformFrameVecToBase(origin);
    };

    /** Inverse rotation. */
    VecType xformBaseVecToFrame(const VecType &target)
    {
        return inv_rotation_ * target;
    };

    /** Inverse transformation. Note that the inverse translation operation is carried out first. */
    VecType shiftBaseStationToFrame(const VecType &target)
    {
        return xformBaseVecToFrame(target - translation_);
    };
};
} // namespace SPH

#endif // GEOMETRIC_PRIMITIVE_H
