#include "base_geometry.h"

namespace SPH
{
//=================================================================================================//
BoundingBox Shape::getBounds()
{
    if (!is_bounds_found_)
    {
        bounding_box_ = findBounds();
        is_bounds_found_ = true;
    }
    return bounding_box_;
}
//=================================================================================================//
bool Shape::checkNotFar(const Vecd &probe_point, Real threshold)
{
    return checkContain(probe_point) || checkNearSurface(probe_point, threshold) ? true : false;
}
//=================================================================================================//
bool Shape::checkNearSurface(const Vecd &probe_point, Real threshold)
{
    Vecd distance = probe_point - findClosestPoint(probe_point);
    return distance.cwiseAbs().maxCoeff() < threshold ? true : false;
}
//=================================================================================================//
Real Shape::findSignedDistance(const Vecd &probe_point)
{
    Real distance_to_surface = (probe_point - findClosestPoint(probe_point)).norm();
    return checkContain(probe_point) ? -distance_to_surface : distance_to_surface;
}
//=================================================================================================//
Vecd Shape::findNormalDirection(const Vecd &probe_point)
{
    bool is_contain = checkContain(probe_point);
    Vecd displacement_to_surface = findClosestPoint(probe_point) - probe_point;
    while (displacement_to_surface.norm() < Eps)
    {
        Vecd jittered = probe_point;
        for (int l = 0; l != probe_point.size(); ++l)
            jittered[l] = probe_point[l] + (((Real)rand() / (RAND_MAX)) - 0.5) * 100.0 * Eps;
        if (checkContain(jittered) == is_contain)
            displacement_to_surface = findClosestPoint(jittered) - jittered;
    }
    Vecd direction_to_surface = displacement_to_surface.normalized();
    return is_contain ? direction_to_surface : -1.0 * direction_to_surface;
}
//=================================================================================================//
bool BinaryShapes::isValid()
{
    return shapes_and_ops_.size() == 0 ? false : true;
}
//=================================================================================================//
BoundingBox BinaryShapes::findBounds()
{
    // initial reference values
    Vecd lower_bound = Infinity * Vecd::Ones();
    Vecd upper_bound = -Infinity * Vecd::Ones();

    for (auto &shape_and_op : shapes_and_ops_)
    {
        BoundingBox shape_bounds = shape_and_op.first->getBounds();
        for (int j = 0; j != Dimensions; ++j)
        {
            lower_bound[j] = SMIN(lower_bound[j], shape_bounds.first_[j]);
            upper_bound[j] = SMAX(upper_bound[j], shape_bounds.second_[j]);
        }
    }
    return BoundingBox(lower_bound, upper_bound);
}
//=================================================================================================//
bool BinaryShapes::checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED)
{
    bool exist = false;
    bool inside = false;

    for (auto &shape_and_op : shapes_and_ops_)
    {
        Shape *geometry = shape_and_op.first;
        ShapeBooleanOps operation_string = shape_and_op.second;
        switch (operation_string)
        {
        case ShapeBooleanOps::add:
        {
            inside = geometry->checkContain(pnt);
            exist = exist || inside;
            break;
        }
        case ShapeBooleanOps::sub:
        {
            inside = geometry->checkContain(pnt);
            exist = exist && (!inside);
            break;
        }
        default:
        {
            std::cout << "\n FAILURE: the boolean operation is not applicable!" << std::endl;
            std::cout << __FILE__ << ':' << __LINE__ << std::endl;
            throw;
        }
        }
    }
    return exist;
}
//=================================================================================================//
Vecd BinaryShapes::findClosestPoint(const Vecd &probe_point)
{
    // a big positive number
    Real large_number(Infinity);
    Real dist_min = large_number;
    Vecd pnt_closest = Vecd::Zero();
    Vecd pnt_found = Vecd::Zero();

    for (auto &shape_and_op : shapes_and_ops_)
    {
        Shape *geometry = shape_and_op.first;
        pnt_found = geometry->findClosestPoint(probe_point);
        Real dist = (probe_point - pnt_found).norm();

        if (dist <= dist_min)
        {
            dist_min = dist;
            pnt_closest = pnt_found;
        }
    }

    return pnt_closest;
}
//=================================================================================================//
ShapeAndOp *BinaryShapes::getShapeAndOpByName(const std::string &shape_name)
{
    for (auto &shape_and_op : shapes_and_ops_)
    {
        Shape *shape = shape_and_op.first;
        if (shape->getName() == shape_name)
            return &shape_and_op;
    }
    std::cout << "\n FAILURE: the shape " << shape_name << " has not been created!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;

    return nullptr;
}
//=================================================================================================//
Shape *BinaryShapes::getShapeByName(const std::string &shape_name)
{
    return getShapeAndOpByName(shape_name)->first;
}
//=================================================================================================//
size_t BinaryShapes::getShapeIndexByName(const std::string &shape_name)
{
    for (size_t index = 0; index != shapes_and_ops_.size(); ++index)
    {
        if (shapes_and_ops_[index].first->getName() == shape_name)
        {
            return index;
        }
    }
    std::cout << "\n FAILURE: the shape " << shape_name << " has not been created!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;

    return MaxSize_t;
}
//=================================================================================================//
} // namespace SPH