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
 * @file sdf_operation.h
 * @brief Operation for signed ditance function (SDF),
 * which is used to combine the primitive SDFs to generate complex geometry.
 * @author Xiangyu Hu
 */

#ifndef SDF_OPERATION_H
#define SDF_OPERATION_H

#include "base_data_type.h"
#include "geometric_primitive.h"
#include "scalar_functions.h"

namespace SPH
{
template <typename OperationType, typename InputType1, typename InputType2>
class SDFOperation
{
    OperationType operation_;
    InputType1 input1_;
    InputType2 input2_;

  public:
    explicit SDFOperation(const OperationType &operation, const InputType1 &input1, const InputType2 &input2)
        : operation_(operation), input1_(input1), input2_(input2) {}
    OperationType &getOpertion() { return operation_; }
    InputType1 &getInput1() { return input1_; }
    InputType2 &getInput2() { return input2_; }
    Real operator()(const Vec3d &point) const { return operation_(point, input1_, input2_); }
    auto findBounds() const { return operation_.findBounds(input1_, input2_); }
};

struct SDFAddition
{
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};

struct SDFSubtraction
{
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const { return input1.findBounds(); }
};

struct SDFIntersection
{
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};

class SDFSmoothAddition
{
    Real finest_grid_spacing_; // scaled with mesh size

  public:
    explicit SDFSmoothAddition(Real finest_grid_spacing) : finest_grid_spacing_(finest_grid_spacing) {}
    void setParameters(Real finest_grid_spacing) { finest_grid_spacing_ = finest_grid_spacing; }
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};

class SDFSmoothSubtraction
{
    Real finest_grid_spacing_; // scaled with mesh size

  public:
    explicit SDFSmoothSubtraction(Real finest_grid_spacing) : finest_grid_spacing_(finest_grid_spacing) {}
    void setParameters(Real finest_grid_spacing) { finest_grid_spacing_ = finest_grid_spacing; }
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};

class SDFSmoothIntersection
{
    Real finest_grid_spacing_; // scaled with mesh size

  public:
    explicit SDFSmoothIntersection(Real finest_grid_spacing) : finest_grid_spacing_(finest_grid_spacing) {}
    void setParameters(Real finest_grid_spacing) { finest_grid_spacing_ = finest_grid_spacing; }
    template <typename Input1, typename Input2>
    Real operator()(const Vec3d &point, const Input1 &input1, const Input2 &input2) const;
    template <typename Input1, typename Input2>
    auto findBounds(const Input1 &input1, const Input2 &input2) const;
};
//----------------------------------------------------------------------
// 3D geometric primitives derived from 2D primitives
//----------------------------------------------------------------------
class SDFExtrusion
{
    Real height_;

  public:
    explicit SDFExtrusion(Real height) : height_(height) {}
    void setParameters(Real height) { height_ = height; }
    template <typename Input2D>
    Real operator()(const Input2D &input, const Vec3d &point) const;
};

class SDFRotation
{
    Real angle_;

  public:
    explicit SDFRotation(Real angle) : angle_(angle) {}
    void setParameters(Real angle) { angle_ = angle; }
    template <typename Input2D>
    Real operator()(const Input2D &input, const Vec3d &point) const;
};
//----------------------------------------------------------------------
// 3D geometric primitives derived from 3D primitives
//----------------------------------------------------------------------
class SDFElongation
{
    Real elongation_factor_;

  public:
    explicit SDFElongation(Real elongation_factor) : elongation_factor_(elongation_factor) {}
    Real getParameters() const { return elongation_factor_; }
    void setParameters(Real elongation_factor) { elongation_factor_ = elongation_factor; }
    template <typename Input3D>
    Real operator()(const Input3D &input, const Vec3d &point) const;
};
} // namespace SPH

#endif // SDF_OPERATION_H
