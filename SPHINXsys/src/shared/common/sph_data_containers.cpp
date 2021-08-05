/**
 * @file 	sph_data_containers.cpp
 * @brief 	Global functions related to sph_data_containers.h
 * @author	Bence Rochlitz, Luhui Han, Chi ZHang and Xiangyu Hu
*/

#include "sph_data_containers.h"

namespace SPH {
//=================================================================================================//
bool checkIfPointInBoundingBox(Vec3d point, BoundingBox& bbox)
{
    return point[0] >= bbox.first[0] && point[0] <= bbox.second[0] &&
        point[1] >= bbox.first[1] && point[1] <= bbox.second[1] &&
        point[2] >= bbox.first[2] && point[2] <= bbox.second[2];
}
//=================================================================================================//
bool checkIfPointInBoundingBox(Vec2d point, BoundingBox& bbox)
{
    return point[0] >= bbox.first[0] && point[0] <= bbox.second[0] &&
        point[1] >= bbox.first[1] && point[1] <= bbox.second[1];
}
//=================================================================================================//
}