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
#ifndef SURFACE_SHAPE_H
#define SURFACE_SHAPE_H

#include "sphinxsys.h"
#include "vector.h"

#include <opencascade/Geom_Surface.hxx>
#include <opencascade/Standard_TypeDef.hxx>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
namespace fs = std::filesystem;

namespace SPH
{
class SurfaceShape : public Shape
{
  public:
    explicit SurfaceShape(const std::string &shape_name)
        : Shape(shape_name) {};
    virtual bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override;
    virtual Vecd findClosestPoint(const Vecd &input_pnt) override;
    Vecd getCartesianPoint(Standard_Real u, Standard_Real v);

    Handle_Geom_Surface surface_;

  protected:
    virtual BoundingBoxd findBounds() override;
};

class SurfaceShapeSTEP : public SurfaceShape
{

  public:
    // constructor for load STEP file from out side
    explicit SurfaceShapeSTEP(Standard_CString &filepathname,
                              const std::string &shape_name = "SurfaceShapeSTEP");
    virtual ~SurfaceShapeSTEP() {};
};

} // namespace SPH

#endif // SURFACE_SHAPE_H
