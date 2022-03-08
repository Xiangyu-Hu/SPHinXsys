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
BoundingBox getIntersectionOfBoundingBoxes(BoundingBox &bb1, BoundingBox &bb2)
{
    // check that the inputs are correct
    int dimension = bb1.first.size();
    if (dimension != bb1.second.size() || dimension != bb2.first.size() || dimension != bb2.second.size())
        std::runtime_error("getIntersectionOfBoundingBoxes: wrong input!");
    // Get the Bounding Box of the intersection of the two meshes
    BoundingBox bb(bb1);
	// #1 check that there is overlap, if not, exception
	for (int i = 0; i < dimension; ++i)
		if (bb2.first[i] > bb1.second[i] || bb2.second[i] < bb1.first[i])
			std::runtime_error("getIntersectionOfBoundingBoxes: no overlap!");
	// #2 otherwise modify the first one to get the intersection
	for (int i = 0; i < dimension; ++i)
	{	
		// if the lower limit is inside change the lower limit
		if (bb1.first[i] < bb2.first[i] && bb2.first[i] < bb1.second[i])
			bb.first[i] = bb2.first[i];
		// if the upper limit is inside, change the upper limit
		if (bb1.second[i] > bb2.second[i] && bb2.second[i] > bb1.first[i])
			bb.second[i] = bb2.second[i];
	}
    return bb;
}
//=================================================================================================//
}