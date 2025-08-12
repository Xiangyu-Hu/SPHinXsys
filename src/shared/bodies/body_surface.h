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
 * @file 	body_surface.h
 * @brief 	Body surfaces for different purposes.
 * @details	tbd
 * @author	Xiangyu Hu
 */

#ifndef BODY_SURFACE_H
#define BODY_SURFACE_H

#include "base_body_part.h"

namespace SPH
{
/**
 * @class BodySurface
 * @brief A  body part as the outer most layer of particles
 */
class BodySurface : public BodyPartByParticle
{
  public:
    explicit BodySurface(SPHBody &sph_body);
    virtual ~BodySurface() {};

  protected:
    Real particle_spacing_min_;
    bool tagNearSurface(size_t particle_index);
};

/**
 * @class BodySurfaceLayer
 * @brief A  body part for several layers of particle near the surface.
 */
class BodySurfaceLayer : public BodyPartByParticle
{
  public:
    explicit BodySurfaceLayer(SPHBody &sph_body, Real layer_thickness = 3.0);
    virtual ~BodySurfaceLayer() {};

  private:
    Real thickness_threshold_;
    bool tagSurfaceLayer(size_t particle_index);
};

/**
 * @class NearShapeSurface
 * @brief A body part with the cell lists near the surface of a prescribed shape.
 * @ details The body part shape can be that of the body,
 * or a sub shape of the body shape, or a shape independent of the body shape.
 * Note that only cells near the surface of the body part shape are included.
 */
class NearShapeSurface : public BodyPartByCell
{
  private:
    UniquePtrKeeper<LevelSetShape> level_set_shape_keeper_;

  public:
    NearShapeSurface(RealBody &real_body, SharedPtr<Shape> shape_ptr);
    NearShapeSurface(RealBody &real_body, LevelSetShape &level_set_shape);
    explicit NearShapeSurface(RealBody &real_body);
    NearShapeSurface(RealBody &real_body, const std::string &sub_shape_name);
    virtual ~NearShapeSurface() {};
    LevelSetShape &getLevelSetShape() { return level_set_shape_; };

  private:
    LevelSetShape &level_set_shape_;
    bool checkNearSurface(Vecd cell_position, Real threshold);
};
} // namespace SPH
#endif // BODY_SURFACE_H
