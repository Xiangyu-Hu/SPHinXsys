/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*                                                                        *
 * ------------------------------------------------------------------------*/


#ifndef SURFACE_SHAPE_H
#define SURFACE_SHAPE_H

#include "base_geometry.h"
#include "simbody_middle.h"
#include "opencascade/STEPCAFControl_Reader.hxx"
#include "opencascade/XSControl_WorkSession.hxx"
#include"opencascade/Transfer_TransientProcess.hxx"
#include"opencascade/TopExp_Explorer.hxx"
#include"opencascade/TopoDS.hxx"
#include"opencascade/BRep_Builder.hxx"
#include"opencascade/GeomAPI_ProjectPointOnSurf.hxx"
#include"opencascade/ProjLib_ProjectOnSurface.hxx"



#include <iostream>
#include <string>
#include <fstream>
#include <filesystem>
namespace fs = std::filesystem;

namespace SPH
{
	
	
	class SurfaceShape : public Shape
        {
              public:
                explicit SurfaceShape(const std::string &shape_name)
                    : Shape(shape_name){};
                virtual bool checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED = true) override;
                virtual Vecd findClosestPoint(const Vecd &input_pnt) override;
                Vecd findActualPoint(Standard_Real u, Standard_Real v);
              
                Handle_Geom_Surface surface_;
              protected:
                virtual BoundingBox findBounds() override;   
        };

        class SurfaceShapeSTEP : public SurfaceShape
        {
               
              public:
                // constructor for load STEP file from out side
                explicit SurfaceShapeSTEP(Standard_CString &filepathname,
                                          const std::string &shape_name = "SurfaceShapeSTEP");
                virtual ~SurfaceShapeSTEP(){};
               
        };
	
}

#endif //SURFACE_SHAPE_H
