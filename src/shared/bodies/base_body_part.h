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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	base_body_part.h
 * @brief 	This is the base classes of body parts.
 * @details	There two main type of body parts. One is part by particle.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_BODY_PART_H
#define BASE_BODY_PART_H

#include "base_body.h"

#include <string>

namespace SPH
{
/**
 * @class BodyPart
 * @brief An auxillary class for SPHBody to indicate a part of the body.
 */
using namespace std::placeholders;
class BodyPart
{
  public:
    BodyPart(SPHBody &sph_body, const std::string &body_part_name)
        : sph_body_(sph_body), body_part_name_(body_part_name){};
    virtual ~BodyPart(){};

    SPHBody &getSPHBody() { return sph_body_; };
    std::string getName() { return body_part_name_; };

  protected:
    SPHBody &sph_body_;
    std::string body_part_name_;
};

/**
 * @class BodyPartByParticle
 * @brief A body part with a collection of particles.
 */
class BodyPartByParticle : public BodyPart
{
  public:
    IndexVector body_part_particles_; /**< Collection particle in this body part. */
    BaseParticles &getBaseParticles() { return base_particles_; };
    IndexVector &LoopRange() { return body_part_particles_; };
    size_t SizeOfLoopRange() { return body_part_particles_.size(); };

    BodyPartByParticle(SPHBody &sph_body, const std::string &body_part_name)
        : BodyPart(sph_body, body_part_name), base_particles_(sph_body.getBaseParticles()),
          body_part_bounds_(Vecd::Zero(), Vecd::Zero()), body_part_bounds_set_(false){};
    virtual ~BodyPartByParticle(){};

    void setBodyPartBounds(BoundingBox bbox)
    {
        body_part_bounds_ = bbox;
        body_part_bounds_set_ = true;
    };

    BoundingBox getBodyPartBounds()
    {
        if (!body_part_bounds_set_)
            std::cout << "WARNING: the body part bounds are not set for BodyPartByParticle." << std::endl;
        return body_part_bounds_;
    }

  protected:
    BaseParticles &base_particles_;
    BoundingBox body_part_bounds_;
    bool body_part_bounds_set_;

    typedef std::function<void(size_t)> TaggingParticleMethod;
    void tagParticles(TaggingParticleMethod &tagging_particle_method);
};

/**
 * @class BodyPartByCell
 * @brief A body part with a collection of cell lists.
 */
class BodyPartByCell : public BodyPart
{
  public:
    ConcurrentCellLists body_part_cells_; /**< Collection of cells to indicate the body part. */
    ConcurrentCellLists &LoopRange() { return body_part_cells_; };
    size_t SizeOfLoopRange();

    BodyPartByCell(RealBody &real_body, const std::string &body_part_name)
        : BodyPart(real_body, body_part_name), cell_linked_list_(real_body.getCellLinkedList()){};
    virtual ~BodyPartByCell(){};

  protected:
    BaseCellLinkedList &cell_linked_list_;
    typedef std::function<bool(Vecd, Real)> TaggingCellMethod;
    void tagCells(TaggingCellMethod &tagging_cell_method);
};

/**
 * @class BodyRegionByParticle
 * @brief A  body part with the collection of particles within by a prescribed shape.
 */
class BodyRegionByParticle : public BodyPartByParticle
{
  private:
    SharedPtrKeeper<Shape> shape_ptr_keeper_;

  public:
    Shape &body_part_shape_;

    BodyRegionByParticle(SPHBody &sph_body, SharedPtr<Shape> shape_ptr);
    virtual ~BodyRegionByParticle(){};

  private:
    void tagByContain(size_t particle_index);
};

/**
 * @class BodySurface
 * @brief A  body part with the collection of particles at surface of a body
 */
class BodySurface : public BodyPartByParticle
{
  public:
    explicit BodySurface(SPHBody &sph_body);
    virtual ~BodySurface(){};

  private:
    Real particle_spacing_min_;
    void tagNearSurface(size_t particle_index);
};

/**
 * @class BodySurfaceLayer
 * @brief A  body part with the collection of particles within the surface layers of a body.
 */
class BodySurfaceLayer : public BodyPartByParticle
{
  public:
    explicit BodySurfaceLayer(SPHBody &sph_body, Real layer_thickness = 3.0);
    virtual ~BodySurfaceLayer(){};

  private:
    Real thickness_threshold_;
    void tagSurfaceLayer(size_t particle_index);
};

/**
 * @class BodyRegionByCell
 * @brief A body part with the cell lists within a prescribed shape.
 */
class BodyRegionByCell : public BodyPartByCell
{
  private:
    SharedPtrKeeper<Shape> shape_ptr_keeper_;

  public:
    Shape &body_part_shape_;

    BodyRegionByCell(RealBody &real_body, SharedPtr<Shape> shape_ptr);
    virtual ~BodyRegionByCell(){};

  private:
    bool checkNotFar(Vecd cell_position, Real threshold);
};

/**
 * @class NearShapeSurface
 * @brief A body part with the cell lists near the surface of a prescribed shape.
 */
class NearShapeSurface : public BodyPartByCell
{
  private:
    UniquePtrKeeper<LevelSetShape> level_set_shape_keeper_;

  public:
    LevelSetShape &level_set_shape_;

    /** for the case that the body part shape is not that of the body */
    NearShapeSurface(RealBody &real_body, SharedPtr<Shape> shape_ptr);
    /** for the case that the body part shape is the surface of the body shape */
    explicit NearShapeSurface(RealBody &real_body);
    /** for the case that the body part shape is one part of the surface of the body shape */
    NearShapeSurface(RealBody &real_body, const std::string &shape_name);
    virtual ~NearShapeSurface(){};

  private:
    /** only cells near the surface of the body part shape are included */
    bool checkNearSurface(Vecd cell_position, Real threshold);
};

/**
 * @class AlignedBoxRegion
 * @brief A template body part with the collection of particles within by an AlignedBoxShape.
 */
template <class BodyRegionType>
class AlignedBoxRegion : public BodyRegionType
{
  public:
    AlignedBoxShape &aligned_box_;

    AlignedBoxRegion(RealBody &real_body, SharedPtr<AlignedBoxShape> aligned_box_ptr)
        : BodyRegionType(real_body, aligned_box_ptr), aligned_box_(*aligned_box_ptr.get()){};
    virtual ~AlignedBoxRegion(){};
};

using BodyAlignedBoxByParticle = AlignedBoxRegion<BodyRegionByParticle>;
using BodyAlignedBoxByCell = AlignedBoxRegion<BodyRegionByCell>;
} // namespace SPH
#endif // BASE_BODY_PART_H
