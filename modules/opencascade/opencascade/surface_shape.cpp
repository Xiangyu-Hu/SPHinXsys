#include "surface_shape.h"

#include <opencascade/STEPCAFControl_Reader.hxx>
#include <opencascade/TopExp_Explorer.hxx>
#include <opencascade/TopoDS.hxx>
#include <opencascade/BRep_Builder.hxx>
#include <opencascade/GeomAPI_ProjectPointOnSurf.hxx>
#include <opencascade/gp_Pnt.hxx>

namespace SPH
{
	//=================================================================================================//
Vecd SurfaceShape::findClosestPoint(const Vecd &input_pnt)
{

    Vecd closest_pnt;
    gp_Pnt point1 = EigenToOcct(input_pnt);
    Extrema_ExtAlgo Algo = Extrema_ExtAlgo_Tree;
    gp_Pnt point2 = GeomAPI_ProjectPointOnSurf(point1, surface_, Algo);
    closest_pnt[0] = point2.X();
    closest_pnt[1] = point2.Y();
    closest_pnt[2] = point2.Z();

    return closest_pnt;
}

////=================================================================================================//
bool SurfaceShape::checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED) { return 0; }
BoundingBoxd SurfaceShape::findBounds() { return BoundingBoxd(); }
//=================================================================================================//

Vecd SurfaceShape::getCartesianPoint(Standard_Real u, Standard_Real v)
{
    gp_Pnt point;
    Vec3d actual_pnt;
    point = surface_->Value(u, v);
    actual_pnt[0] = point.X();
    actual_pnt[1] = point.Y();
    actual_pnt[2] = point.Z();
    return Vecd(actual_pnt[0], actual_pnt[1], actual_pnt[2]);
}

//=================================================================================================//
SurfaceShapeSTEP::
    SurfaceShapeSTEP(Standard_CString &filepathname, const std::string &shape_name)
    : SurfaceShape(shape_name)
{
    if (!std::filesystem::exists(filepathname))
    {
        std::cout << "\n Error: the input file:" << filepathname << " is not exists" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    STEPControl_Reader step_reader;
    step_reader.ReadFile(filepathname);

    for (Standard_Integer i = 1; i <= step_reader.NbRootsForTransfer(); i++)
        step_reader.TransferRoot(i);

    TopoDS_Shape step_shape;
    for (Standard_Integer i = 1; i <= step_reader.NbShapes(); i++)
        step_shape = step_reader.Shape(i);

    TopExp_Explorer explorer(step_shape, TopAbs_FACE);
    TopoDS_Face face = TopoDS::Face(explorer.Current());
    surface_ = BRep_Tool::Surface(face);
}

	//=================================================================================================//
}
