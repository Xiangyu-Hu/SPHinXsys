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
 * @file sdf_extension.h
 * @brief Extension for signed distance function (SDF),
 * which is used to modify the primitive SDFs to generate complex geometry.
 * @author	Xiangyu Hu
 */

#ifndef SDF_EXTENSION_H
#define SDF_EXTENSION_H

#include "base_data_type.h"
#include "geometric_primitive.h"
#include "scalar_functions.h"

namespace SPH
{
template <typename InputType, typename ExtensionType>
class SDFExtension
{
    InputType input_;
    ExtensionType extension_;

  public:
    explicit SDFExtension(const InputType &input, const ExtensionType &extension);
    ExtensionType &getExtension() { return extension_; }
    InputType &getInput() { return input_; }
    Real operator()(const Vec3d &point) const { return extension_(input_, point); }
    BoundingBox3d findBounds() const { return extension_.findBounds(input_); }
};

class SDFRound
{
    Real radius_;

  public:
    explicit SDFRound(Real radius) : radius_(radius) {}
    void setParameters(Real radius) { radius_ = radius; }
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const { return input(point) - radius_; }
    template <typename InputType>
    auto findBounds(const InputType &input) const { return input.findBounds().expand(radius_); }
};

class SDFOnion
{
    Real radius_;

  public:
    explicit SDFOnion(Real radius) : radius_(radius) {}
    void setParameters(Real radius) { radius_ = radius; }
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const { return ABS(input(point)) - radius_; }
    template <typename InputType>
    auto findBounds(const InputType &input) const { return input.findBounds().expand(SMAX(radius_, Real(0))); }
};

class SDFScale
{
    Real scale_factor_;

  public:
    explicit SDFScale(Real scale_factor) : scale_factor_(scale_factor) {}
    Real getParameters() const { return scale_factor_; }
    void setParameters(Real scale_factor) { scale_factor_ = scale_factor; }
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const;
    template <typename InputType>
    auto findBounds(const InputType &input) const;
};

class SDFTransform
{
    Transform3d transform_;

  public:
    template <typename... Args>
    explicit SDFTransform(Args &&...args) : transform_(std::forward<Args>(args)...) {}
    template <typename... Args>
    void setParameters(Args &&...args);
    template <typename InputType>
    Real operator()(const InputType &input, const Vec3d &point) const;
    template <typename InputType>
    auto findBounds(const InputType &input) const;
};
} // namespace SPH

#endif // SDF_EXTENSION_H
