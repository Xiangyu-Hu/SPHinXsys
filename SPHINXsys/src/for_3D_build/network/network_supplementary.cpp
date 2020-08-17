/**
 * @file 	branch.cpp
 * @brief 	Here, Functions belong to BaseBranch are given.
 * @author	chi ZHang and Xiangyu Hu
 * @version	0.2.0
 */
 #include "network.h"
 //=================================================================================================//
namespace SPH
{
    //=================================================================================================//
    Vecd Node::getGradient(Point pt, Real delta)
    {
        Vec3d point = pt;
        Vec3d dx(delta, 0.0, 0.0);
        Vec3d dy(0.0, delta, 0.0);
        Vec3d dz(0.0, 0.0, delta);

        Real dist_x_m = getDistanceFromPoint(Point(point - dx));
        Real dist_x_p = getDistanceFromPoint(Point(point + dx));
        Real dist_y_m = getDistanceFromPoint(Point(point - dy));
        Real dist_y_p = getDistanceFromPoint(Point(point + dy));
        Real dist_z_m = getDistanceFromPoint(Point(point - dz));
        Real dist_z_p = getDistanceFromPoint(Point(point + dz));

        Real delta_x = (dist_x_p - dist_x_m) / (2.0 * delta);
        Real delta_y = (dist_y_p - dist_y_m) / (2.0 * delta);
        Real delta_z = (dist_z_p - dist_z_m) / (2.0 * delta);

        Vec3d grad(delta_x, delta_y, delta_z); 
        return grad;
    }
    //=================================================================================================//
}
//=================================================================================================//