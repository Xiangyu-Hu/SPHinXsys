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
* @file implicit_surface.h
* @brief A representation of implicit surface based on structured grid is define in this file
* @details A implicit surface means surface defined by signed distance function. Here, we define the ImplicitSurface
*          class to supply some useful utilities to costruct the implicit surface. 
* @author	Huiqiang Yue
*/

#ifndef __IMPLICIT_SURFACE_H__
#define __IMPLICIT_SURFACE_H__

#include "base_geometry.h"
#include "base_mesh.h"

#include <iostream>
#include <string>
#include <fstream>
#include <array>
#include <vector>

namespace SPH
{
    /**
     * The class represents a implicit surface defined on axis aligned structured mesh.
     * 
     * The data value stored at the center of the mesh cell represents the "level set", i.e.,
     * the distance to the surface. This is similar to the distance map which is used widely in 
     * image processing. The main goal of this class is to process the 3D image data.
     * 
     */
    class ImplicitSurface : public Shape
    {
    public:
        ImplicitSurface() : Shape("ImplicitSurface") { }
        virtual ~ImplicitSurface();

        virtual Vec3d findClosestPoint(const Vec3d& input_pnt) override;

        virtual BoundingBox findBounds() override;

        Vecd getDataOrigin() const;

        /**
         * Construct a sphere
         * @param center Center of the sphere
         * @param radius Radius of the sphere
         * @param bounding_box Axis aligned bounding box of the sphere. Note: the bounding box should more
         *                      larger than the minmum bouding box of the sphere.
         * @param grid_spacing The grid spacing of the mesh which signed distance defined on it.
         */
        void addAImplicitSphere(Vec3d center, Real radius, BoundingBox bounding_box, Real grid_spacing);
        /**
         * Load image data from .mhd image file
         * @note: At present, we only support the ".mhd" image files.
         */
        void loadFromMHDImge(const std::string &filename);

        Real getValueAtArbitraryPoint(Vec3d point);
        Vec3d getGradientAtArbitraryPoint(Vec3d point);

        Real findSignedDistance(const Vec3d &input_pnt);
        bool checkContain(const Vec3d &input_pnt, bool BOUNDARY_INCLUDED = true);
        bool checkNotFar(const Vec3d &input_pnt, Real threshold);
        bool checkNearSurface(const Vec3d &input_pnt, Real threshold);
        Vec3d findNormalDirection(const Vec3d &input_pnt);

        void writeMeshFieldToPlt(const std::string &file);
        void writeMeshFieldToVtu(const std::string &file);
    private:
        void initialize(BoundingBox bounding_box, Real grid_spacing);
        size_t convertToOneDIndex(Vec3u index);
        Real getCellCenterData(Vec3u index);   
        Vec3d getCellCenterGrad(Vec3u index);  
        /**
         * Compute the gradient at cell center using central finite difference scheme.
         */
        void computeGradient();
        /**
         * Compute the cell indices and the interpolation weights using the bilinear interpolation.
         */
        std::pair<std::array<Vec3u, 8>, std::array<Real, 8>> computeIndexAndWeights(const Vec3d &point);
    private:
        Real grid_spacing_;
        BoundingBox bounding_box_;
        
        Vecu number_of_cells_;
        Vecu number_of_grid_points_;
        Vecd mesh_lower_bound_;
        Vecd mesh_upper_bound_;
        /**
         * The data is stored at center of the mesh, so data origin is different from the mesh origin.
         * Therefore, dataOrigin_ is the center of the lower left cell.
         */
        Vecd dataOrigin_;

        std::vector<Real> data_;
        std::vector<Real> gradDataX_;
        std::vector<Real> gradDataY_;
        std::vector<Real> gradDataZ_;

        int maxNumberOfIterations_ = 5;
        Real epsilon_ = 1e-6;
    };


}





#endif // __IMPLICIT_SHAPE_H__