/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * --------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
 * and HU1527/12-1.															*
 *                                                                           *
 * Portions copyright (c) 2017-2020 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * --------------------------------------------------------------------------*/
/**
 * @file boundary_face.h
 * @brief Here, we define several data structures to be used in constructing complex boundary conditions.
 * @author	Huiqiang Yue and Anyong Zhang
 */

#ifndef BOUNDARY_FACE_H
#define BOUNDARY_FACE_H

#include "all_particle_dynamics.h"

namespace SPH
{
    /**
     * @brief The class represents a face at inlet or outlet boundary
     */
    class SegmentFace
    {
    public:
        /**
         * @brief Construct a new Segment Face object
         * @param[in] boundary_points    Points to describe the border of face
         * @param[in] direction          Direction pointing from face to computational domain
         * @param[in] center             Center of the face
         */
        SegmentFace(StdVec<Vecd> boundary_points, Vecd direction, Vecd center);
        SegmentFace(StdVec<Vecd> boundary_points, Vecd direction);
        SegmentFace(const Vecd& center, const Vecd& direction);

        const StdVec<Vecd>& getBoundaryPoints()
        {
            return boundary_points_;
        }

        const Vecd& direction()
        {
            return direction_;
        }

        const Vecd& center()
        {
            return center_;
        }

    private:
        StdVec<Vecd> boundary_points_;
        Vecd direction_;
        Vecd center_;
    };

    /**
     * @brief The class constructs a body region with a defined SegmentFace object
     */
    class BodyRegionWithFace
    {
    public:
        /**
         * @brief Construct a new Body Region With a SegmentFace object
         * @param[in] real_body          Fluid body to define boundary condition
         * @param[in] segment_face       SegmentFace object represents the inlet or outlet face
         * @param[in] scale              Factor to define the length of the boundary region
         */
        BodyRegionWithFace(RealBody& real_body, SegmentFace& segment_face, Real scale = 0.0);
        Real getRegionWidth() const
        {
            return region_width_;
        }

        /**
         * @brief Get the distance between the center of face and particles,
         *          positive for pointing from border to computational domain
         * @param[in] position           Positon of particle
         * @return Real                 Signed distance
         */
        Real getSignedDistance(const Vecd& position) const;

        const Vecd& getDirectionToFluid() const
        {
            return segment_face_.direction();
        }

    protected:
        SegmentFace& segment_face_;
        Real region_width_;
        BoundingBox region_bounds;

    private:
        /**
         * @brief Construct the body region and calculate its border.
         */
        void calculate_region_bounds();
    };

    /**
     * @brief A body part with the cell lists containing a boundary face.
     */
    class BodyRegionByCellsWithFace : public BodyRegionWithFace
    {
    public:
        BodyRegionByCellsWithFace(RealBody& real_body, SegmentFace& segment_face, Real scale = 0.0);
        void tagBodyDomainBoundingCells(CellLists& bound_cells);

    protected:
        RealBody& real_body_;
    };

    class PartDynamicsByCellsWithFace : public ParticleDynamics<void>
    {
    public:
        PartDynamicsByCellsWithFace(RealBody& real_body, BodyRegionByCellsWithFace& body_region);
        virtual ~PartDynamicsByCellsWithFace() {}

        virtual void exec(Real dt = 0.0) override;
        virtual void parallel_exec(Real dt = 0.0) override { exec(dt); }

    protected:
        BodyRegionByCellsWithFace& body_region_;
        CellLists bound_cells_;
        ParticleFunctor particle_functor_;
    };

    class PartSimpleDynamicsByCellsWithFace : public PartDynamicsByCellsWithFace
    {
    public:
        PartSimpleDynamicsByCellsWithFace(RealBody& real_body, BodyRegionByCellsWithFace& body_region);
        virtual ~PartSimpleDynamicsByCellsWithFace() {}

    protected:
        virtual void Update(size_t index_i, Real dt = 0.0) = 0;
    };

} // namespace SPH


#endif // !BOUNDARY_FACE_H

