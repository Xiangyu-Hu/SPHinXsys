#include "base_geometry.h"

#include "io_log.h"
namespace SPH
{
//=================================================================================================//
Shape::Shape(const std::string &shape_name)
    : name_(shape_name), is_bounds_found_(false), logger_(Log::init()) {}
//=================================================================================================//
BoundingBoxd Shape::getBounds()
{
    if (!is_bounds_found_)
    {
        bounding_box_ = findBounds();
        is_bounds_found_ = true;
    }

    if (bounding_box_.BoundSize().norm() < SqrtEps)
    {
        std::cout << "\n Error: the Bounding box is unreasonably small! " << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
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
            jittered[l] = probe_point[l] + rand_uniform(-0.5, 0.5) * 100.0 * Eps;
        if (checkContain(jittered) == is_contain)
            displacement_to_surface = findClosestPoint(jittered) - jittered;
    }
    Vecd direction_to_surface = displacement_to_surface.normalized();
    return is_contain ? direction_to_surface : -1.0 * direction_to_surface;
}
//=================================================================================================//
bool BinaryShapes::isValid()
{
    return sub_shapes_and_ops_.size() == 0 ? false : true;
}
//=================================================================================================//
BoundingBoxd BinaryShapes::findBounds()
{
    // initial reference values
    Vecd lower_bound = MaxReal * Vecd::Ones();
    Vecd upper_bound = MinReal * Vecd::Ones();

    for (auto &sub_shape_and_op : sub_shapes_and_ops_)
    {
        BoundingBoxd shape_bounds = sub_shape_and_op.first->getBounds();
        for (int j = 0; j != Dimensions; ++j)
        {
            lower_bound[j] = SMIN(lower_bound[j], shape_bounds.lower_[j]);
            upper_bound[j] = SMAX(upper_bound[j], shape_bounds.upper_[j]);
        }
    }
    return BoundingBoxd(lower_bound, upper_bound);
}
//=================================================================================================//
bool BinaryShapes::checkContain(const Vecd &pnt, bool BOUNDARY_INCLUDED)
{
    bool exist = false;
    bool inside = false;

    for (auto &sub_shape_and_op : sub_shapes_and_ops_)
    {
        Shape *geometry = sub_shape_and_op.first;
        ShapeBooleanOps operation_string = sub_shape_and_op.second;
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
    Real large_number(MaxReal);
    Real dist_min = large_number;
    Vecd pnt_closest = Vecd::Zero();
    Vecd pnt_found = Vecd::Zero();

    for (auto &sub_shape_and_op : sub_shapes_and_ops_)
    {
        Shape *geometry = sub_shape_and_op.first;
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
SubShapeAndOp *BinaryShapes::getSubShapeAndOpByName(const std::string &name)
{
    for (auto &sub_shape_and_op : sub_shapes_and_ops_)
    {
        Shape *shape = sub_shape_and_op.first;
        if (shape->getName() == name)
            return &sub_shape_and_op;
    }
    std::cout << "\n FAILURE: the shape " << name << " has not been created!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;

    return nullptr;
}
//=================================================================================================//
Shape *BinaryShapes::getSubShapeByName(const std::string &name)
{
    return getSubShapeAndOpByName(name)->first;
}
//=================================================================================================//
size_t BinaryShapes::getSubShapeIndexByName(const std::string &name)
{
    for (size_t index = 0; index != sub_shapes_and_ops_.size(); ++index)
    {
        if (sub_shapes_and_ops_[index].first->getName() == name)
        {
            return index;
        }
    }
    std::cout << "\n FAILURE: the shape " << name << " has not been created!" << std::endl;
    std::cout << __FILE__ << ':' << __LINE__ << std::endl;

    return MaxUnsignedInt;
}
//=================================================================================================//
} // namespace SPH