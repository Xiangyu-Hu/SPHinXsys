#include "multi_polygon_shape.h"

#include "io_environment.h"

using namespace bg;

#include "earcut.hpp"

#ifdef SPHINXSYS_USE_VTK
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkTriangle.h>
#include <vtkUnsignedIntArray.h>
#include <vtkXMLPolyDataWriter.h>
#endif

namespace SPH
{
namespace
{
struct TriangulatedPolygon
{
    NRings rings;
    std::vector<uint32_t> indices;
};

TriangulatedPolygon triangulatePolygon(const boost_poly &poly)
{
    TriangulatedPolygon triangulated;

    // 1. Add Outer Ring
    std::vector<NPoint> outer_ring;
    for (auto const &p : poly.outer())
    {
        outer_ring.push_back({bg::get<0>(p), bg::get<1>(p)});
    }
    triangulated.rings.push_back(outer_ring);

    // 2. Add Interior Rings (Holes)
    for (auto const &hole : poly.inners())
    {
        std::vector<NPoint> inner_ring;
        for (auto const &p : hole)
        {
            inner_ring.push_back({bg::get<0>(p), bg::get<1>(p)});
        }
        triangulated.rings.push_back(inner_ring);
    }

    // 3. Run triangulation (3 indices per triangle)
    triangulated.indices = mapbox::earcut<uint32_t>(triangulated.rings);

    return triangulated;
}
} // namespace

//=================================================================================================//
MultiPolygon::MultiPolygon(const std::vector<Vecd> &points)
    : MultiPolygon()
{
    addPolygon(points, GeometricOps::add);
}
//=================================================================================================//
MultiPolygon::MultiPolygon(const Vecd &center, Real radius, int resolution)
    : MultiPolygon()
{
    addCircle(center, radius, resolution, GeometricOps::add);
}
//=================================================================================================//
boost_multi_poly MultiPolygon::
    MultiPolygonByBooleanOps(boost_multi_poly multi_poly_in,
                             boost_multi_poly multi_poly_op, GeometricOps boolean_op)
{
    boost_multi_poly multi_poly_tmp_in = multi_poly_in;
    /**
     * Out multi-poly need to be emtpy
     * otherwise the operation is not valid.
     */
    boost_multi_poly multi_poly_tmp_out;

    switch (boolean_op)
    {
    case GeometricOps::add:
    {
        bg::union_(multi_poly_tmp_in, multi_poly_op, multi_poly_tmp_out);
        break;
    }

    case GeometricOps::sub:
    {
        bg::difference(multi_poly_tmp_in, multi_poly_op, multi_poly_tmp_out);
        break;
    }
    case GeometricOps::sym_diff:
    {
        bg::sym_difference(multi_poly_tmp_in, multi_poly_op, multi_poly_tmp_out);
        break;
    }
    case GeometricOps::intersect:
    {
        bg::intersection(multi_poly_tmp_in, multi_poly_op, multi_poly_tmp_out);
        break;
    }
    default:
    {
        std::cout << "\n FAILURE: the type of boolean operation is undefined!" << std::endl;
        std::cout << "\n Please check the boost library reference." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }
    }
    return multi_poly_tmp_out;
}
//=================================================================================================//
void MultiPolygon::addMultiPolygon(const MultiPolygon &multi_polygon_op, GeometricOps op)
{
    multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, multi_polygon_op.getBoostMultiPoly(), op);
}
//=================================================================================================//
void MultiPolygon::addBoostMultiPoly(boost_multi_poly &boost_multi_poly_op, GeometricOps op)
{
    multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, boost_multi_poly_op, op);
}
//=================================================================================================//
void MultiPolygon::addBox(Transform transform, const Vecd &halfsize, GeometricOps op)
{
    Vecd point0 = transform.shiftFrameStationToBase(-halfsize);
    Vecd point1 = transform.shiftFrameStationToBase(Vecd(-halfsize[0], halfsize[1]));
    Vecd point2 = transform.shiftFrameStationToBase(halfsize);
    Vecd point3 = transform.shiftFrameStationToBase(Vecd(halfsize[0], -halfsize[1]));

    std::vector<Vecd> points = {point0, point1, point2, point3, point0};
    addPolygon(points, op);
}
//=================================================================================================//
void MultiPolygon::addBox(const BoundingBox2d &bounding_box, GeometricOps op)
{
    Vecd point0 = bounding_box.lower_;
    Vecd point1 = Vecd(bounding_box.lower_[0], bounding_box.upper_[1]);
    Vecd point2 = bounding_box.upper_;
    Vecd point3 = Vecd(bounding_box.upper_[0], bounding_box.lower_[1]);

    std::vector<Vecd> points = {point0, point1, point2, point3, point0};
    addPolygon(points, op);
}
//=================================================================================================//
void MultiPolygon::addContainerBox(const BoundingBox2d &bounding_box, Real thickness, GeometricOps op)
{
    BoundingBox2d outer_box = bounding_box.expand(thickness);
    MultiPolygon container_box;
    container_box.addBox(outer_box, GeometricOps::add);
    container_box.addBox(bounding_box, GeometricOps::sub);
    addMultiPolygon(container_box, op);
}
//=================================================================================================//
void MultiPolygon::addCircle(const Vecd &center, Real radius, int resolution, GeometricOps op)
{
    Vecd buffer_center = center;
    Real buffer_radius = radius;
    int buffer_res = resolution;

    // Declare the point_circle strategy
    bg::strategy::buffer::join_round join_strategy;
    bg::strategy::buffer::end_round end_strategy;
    bg::strategy::buffer::side_straight side_strategy;
    bg::strategy::buffer::point_circle circle_strategy(buffer_res);
    bg::strategy::buffer::distance_symmetric<Real> circle_dist_strategy(buffer_radius);

    // Create the buffer of a multi point
    boost_point circle_center_pnt;

    bg::set<0>(circle_center_pnt, buffer_center[0]);
    bg::set<1>(circle_center_pnt, buffer_center[1]);

    boost_multi_poly multi_poly_circle;
    buffer(circle_center_pnt, multi_poly_circle,
           circle_dist_strategy, side_strategy,
           join_strategy, end_strategy, circle_strategy);

    if (!is_valid(multi_poly_circle))
    {
        std::cout << "\n Error: the multi polygon is not valid." << std::endl;
        std::cout << "\n The points must be in clockwise. Please check the boost library reference." << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        throw;
    }

    multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, multi_poly_circle, op);
}
//=================================================================================================//
void MultiPolygon::addPolygon(const std::vector<Vecd> &points, GeometricOps op)
{
    std::vector<boost_point> pts;
    for (const Vecd &pnt : points)
    {
        pts.push_back(boost_point(pnt[0], pnt[1]));
    }

    boost_poly poly;
    append(poly, pts);
    if (!is_valid(poly))
    {
        std::cout << "\n Try to reverse the points to clockwise." << std::endl;
        poly.clear();
        std::vector<boost_point> pts_reverse(pts.rbegin(), pts.rend());
        append(poly, pts_reverse);
        if (!is_valid(poly))
        {
            std::cout << "\n Error: the multi polygon is still not valid. Please check the boost library reference." << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            throw;
        }
    }

    boost_multi_poly multi_poly_polygon;
    convert(poly, multi_poly_polygon);

    multi_poly_ = MultiPolygonByBooleanOps(multi_poly_, multi_poly_polygon, op);
}
//=================================================================================================//
void MultiPolygon::
    addPolygonFromFile(std::string file_path_name, GeometricOps op, Vecd translation, Real scale_factor)
{
    std::fstream dataFile(file_path_name);
    Vecd temp_point;
    std::vector<Vecd> coordinates;
    Real temp1 = 0.0, temp2 = 0.0;
    if (dataFile.fail())
    {
        std::cout << "File can not open.\n"
                  << std::endl;
        ;
    }

    while (!dataFile.fail() && !dataFile.eof())
    {
        dataFile >> temp1 >> temp2;
        temp_point[0] = temp1 * scale_factor + translation[0];
        temp_point[1] = temp2 * scale_factor + translation[1];
        coordinates.push_back(temp_point);
    }
    dataFile.close();

    addPolygon(coordinates, op);
}
//=================================================================================================//
bool MultiPolygon::checkContain(const Vec2d &probe_point, bool BOUNDARY_INCLUDED /*= true*/)
{
    if (BOUNDARY_INCLUDED)
    {
        return covered_by(boost_point(probe_point[0], probe_point[1]), multi_poly_);
    }
    else
    {
        return within(boost_point(probe_point[0], probe_point[1]), multi_poly_);
    }
}
//=================================================================================================//
Vecd MultiPolygon::findClosestPoint(const Vecd &probe_point)
{
    /**
     * typedef model::segment<model::d2::point_xy<Real>> boost_seg;
     * From the documentation on segment and referring_segment, the only difference between the two is that
     * referring_segment holds a reference to the points.
     * This is what is needed in a for each that modifies the segment since the points modified should be
     * reflected in the line string. In a for each that does not modify the points, it should still take a
     * reference (most likely a const reference) since it reduces the amount of copying.
     */
    boost_point input_p(probe_point[0], probe_point[1]);
    bg::model::segment<boost_point> closest_seg;
    Real closest_dist_2seg = boost::numeric::bounds<Real>::highest();
    std::function<void(boost_seg)> findclosestsegment = [&closest_seg, &closest_dist_2seg, &input_p](boost_seg seg)
    {
        Real dist = bg::distance(input_p, seg);
        if (dist < closest_dist_2seg)
        {
            closest_dist_2seg = dist;
            // closest_seg.append(seg);
            Real x0 = bg::get<0, 0>(seg);
            Real y0 = bg::get<0, 1>(seg);
            Real x1 = bg::get<1, 0>(seg);
            Real y1 = bg::get<1, 1>(seg);
            bg::set<0, 0>(closest_seg, x0);
            bg::set<0, 1>(closest_seg, y0);
            bg::set<1, 0>(closest_seg, x1);
            bg::set<1, 1>(closest_seg, y1);
        }
    };
    bg::for_each_segment(multi_poly_, findclosestsegment);

    Vecd p_find = Vecd::Zero();

    Real x0 = bg::get<0, 0>(closest_seg);
    Real y0 = bg::get<0, 1>(closest_seg);
    Real x1 = bg::get<1, 0>(closest_seg);
    Real y1 = bg::get<1, 1>(closest_seg);
    Vecd p_0(x0, y0);
    Vecd p_1(x1, y1);
    Vecd vec_v = p_1 - p_0;
    Vecd vec_w = probe_point - p_0;

    Real c1 = vec_v.dot(vec_w);
    if (c1 <= 0)
    {
        p_find = p_0;
    }
    else
    {
        Real c2 = vec_v.dot(vec_v);
        if (c2 <= c1)
        {
            p_find = p_1;
        }
        else
        {
            p_find = p_0 + vec_v * c1 / c2;
        }
    }

    return p_find;
}
//=================================================================================================//
BoundingBoxd MultiPolygon::findBounds()
{
    Vecd lower_bound = Vecd::Zero();
    Vecd upper_bound = Vecd::Zero();
    typedef bg::model::box<model::d2::point_xy<Real>> box;
    lower_bound[0] = bg::return_envelope<box>(multi_poly_).min_corner().get<0>();
    lower_bound[1] = bg::return_envelope<box>(multi_poly_).min_corner().get<1>();
    upper_bound[0] = bg::return_envelope<box>(multi_poly_).max_corner().get<0>();
    upper_bound[1] = bg::return_envelope<box>(multi_poly_).max_corner().get<1>();
    return BoundingBoxd(lower_bound, upper_bound);
}
//=================================================================================================//
void MultiPolygonShape::writeMultiPolygonShapeToVtp()
{
    std::string filefullpath = IO::getEnvironment().OutputFolder() + "/" + getName() + ".vtp";

#ifdef SPHINXSYS_USE_VTK
    vtkNew<vtkPoints> vtk_points;
    vtkNew<vtkCellArray> vtk_cells;

    vtkIdType global_point_offset = 0;
    const auto &multi_poly = multi_polygon_.getBoostMultiPoly();

    for (const auto &poly : multi_poly)
    {
        TriangulatedPolygon triangulated = triangulatePolygon(poly);

        // 4. Add points to VTK and create Triangle Cells
        // We add all points from all rings of THIS polygon
        vtkIdType poly_start_offset = global_point_offset;
        for (const auto &ring : triangulated.rings)
        {
            for (const auto &p : ring)
            {
                vtk_points->InsertNextPoint(p[0], p[1], 0.0);
                global_point_offset++;
            }
        }

        // 5. Create the VTK triangles using the indices
        for (size_t i = 0; i < triangulated.indices.size(); i += 3)
        {
            vtkNew<vtkTriangle> vtk_tri;
            vtk_tri->GetPointIds()->SetId(0, poly_start_offset + triangulated.indices[i]);
            vtk_tri->GetPointIds()->SetId(1, poly_start_offset + triangulated.indices[i + 1]);
            vtk_tri->GetPointIds()->SetId(2, poly_start_offset + triangulated.indices[i + 2]);
            vtk_cells->InsertNextCell(vtk_tri);
        }
    }

    vtkNew<vtkPolyData> polyData;
    polyData->SetPoints(vtk_points);
    polyData->SetPolys(vtk_cells);

    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetInputData(polyData);
    writer->SetFileName(filefullpath.c_str());
    writer->SetDataModeToAscii();
    writer->Write();
#else
    if (fs::exists(filefullpath))
    {
        fs::remove(filefullpath);
    }
    std::ofstream out_file(filefullpath.c_str(), std::ios::trunc);

    std::vector<NPoint> vtk_points;
    std::vector<int> connectivity;
    std::vector<int> offsets;

    int global_point_offset = 0;
    const auto &multi_poly = multi_polygon_.getBoostMultiPoly();

    for (const auto &poly : multi_poly)
    {
        TriangulatedPolygon triangulated = triangulatePolygon(poly);

        // 4. Add points and create triangle connectivity
        // We add all points from all rings of THIS polygon
        int poly_start_offset = global_point_offset;
        for (const auto &ring : triangulated.rings)
        {
            for (const auto &p : ring)
            {
                vtk_points.push_back(p);
                global_point_offset++;
            }
        }

        // 5. Create triangle cells from the earcut indices
        for (size_t i = 0; i < triangulated.indices.size(); i += 3)
        {
            connectivity.push_back(poly_start_offset + static_cast<int>(triangulated.indices[i]));
            connectivity.push_back(poly_start_offset + static_cast<int>(triangulated.indices[i + 1]));
            connectivity.push_back(poly_start_offset + static_cast<int>(triangulated.indices[i + 2]));
            offsets.push_back(static_cast<int>(connectivity.size()));
        }
    }

    out_file << "<?xml version=\"1.0\"?>\n";
    out_file << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\">\n";
    out_file << "<PolyData>\n";
    out_file << "<Piece NumberOfPoints=\"" << vtk_points.size() << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" "
             << "NumberOfStrips=\"0\" NumberOfPolys=\"" << offsets.size() << "\">\n";

    out_file << "<Points>\n";
    out_file << "<DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n";
    for (const auto &p : vtk_points)
    {
        out_file << p[0] << " " << p[1] << " 0\n";
    }
    out_file << "</DataArray>\n";
    out_file << "</Points>\n";

    out_file << "<Polys>\n";
    out_file << "<DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n";
    for (size_t i = 0; i < connectivity.size(); ++i)
    {
        out_file << connectivity[i] << ((i + 1 < connectivity.size()) ? " " : "\n");
    }
    out_file << "</DataArray>\n";

    out_file << "<DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n";
    for (size_t i = 0; i < offsets.size(); ++i)
    {
        out_file << offsets[i] << ((i + 1 < offsets.size()) ? " " : "\n");
    }
    out_file << "</DataArray>\n";
    out_file << "</Polys>\n";

    out_file << "</Piece>\n";
    out_file << "</PolyData>\n";
    out_file << "</VTKFile>\n";

    out_file.close();
#endif
}
//=================================================================================================//
bool MultiPolygonShape::isValid()
{
    return !multi_polygon_.getBoostMultiPoly().empty();
}
//=================================================================================================//
bool MultiPolygonShape::checkContain(const Vecd &probe_point, bool BOUNDARY_INCLUDED)
{
    return multi_polygon_.checkContain(probe_point, BOUNDARY_INCLUDED);
}
//=================================================================================================//
Vecd MultiPolygonShape::findClosestPoint(const Vecd &probe_point)
{
    return multi_polygon_.findClosestPoint(probe_point);
}
//=================================================================================================//
BoundingBoxd MultiPolygonShape::findBounds()
{
    return multi_polygon_.findBounds();
}
//=================================================================================================//
} // namespace SPH