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
 * @file 	outflow_boundary.h
 * @brief 	Here, we define the algorithm classes for fluid dynamics within the outlet region.
 * @details 	These classes allow the boundary conditon to be constructed in any direction.
 * @author	Huiqiang Yue
 */

#ifndef OUTFLOW_BOUNDARY_H
#define OUTFLOW_BOUNDARY_H

#include "boundary_face.h"

namespace SPH
{
    /**
     * @brief Define the boundary condition at outlet with the boundary face.
     */
    class OutflowConditionWithFace : public PartSimpleDynamicsByCellsWithFace,
        public DataDelegateSimple<FluidBody, FluidParticles, Fluid>
    {
    public:
        OutflowConditionWithFace(FluidBody& fluid_body, BodyRegionByCellsWithFace& body_region);
        virtual ~OutflowConditionWithFace() {}

    protected:
        StdLargeVec<Vecd>& pos_n_;
        virtual void Update(size_t index_i, Real dt = 0.0) final;

        virtual void UpdateOutletRegionParticles(size_t index_i, Real dt = 0) = 0;
    };

    /**
     * @brief An open boundary condition at outlet.
     * transfer fluid particles to buffer particles at outlet.
     */
    class DoNothingConditionWithFace : public OutflowConditionWithFace
    {
    public:
        explicit DoNothingConditionWithFace(FluidBody& fluid_body, BodyRegionByCellsWithFace& body_region)
            : OutflowConditionWithFace(fluid_body, body_region)
        {
        }

    protected:
        void UpdateOutletRegionParticles(size_t index_i, Real dt) override
        {
        }
    };

    /**
     * @brief A modification to the do-nothing condition
     * In this case the velocity of outlet region particles is extrapolated from the fluid.
     * @note For more information, please see the parper: 
     * An improved non-reflecting outlet boundary condition for weakly-compressible SPH
     * urls: https://arxiv.org/pdf/1907.04034.pdf
     */
    class ModifiedDoNothingConditionWithFace : public OutflowConditionWithFace
    {
    public:
        explicit ModifiedDoNothingConditionWithFace(FluidBody& fluid_body, BodyRegionByCellsWithFace& body_region);
        virtual ~ModifiedDoNothingConditionWithFace() {}

    protected:
        StdLargeVec<Vecd>& vel_n_;
        CellLinkedList* cell_linked_list_;
        NeighborRelationInner relation_inner_;

    protected:
        virtual void UpdateOutletRegionParticles(size_t index_i, Real dt = 0.0) override;
    };

} // end namespace SPH

#endif // !OUTFLOW_BOUNDARY_H
