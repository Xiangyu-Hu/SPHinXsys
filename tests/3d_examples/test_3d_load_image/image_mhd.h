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
/**
 * @file 	image_mesh_shape.h
 * @brief 	Image process for geometry representation.
 * @author	Yijin Mao
 */

#ifndef IMAGE_MHD_H
#define IMAGE_MHD_H

#ifndef __EMSCRIPTEN__

#include "sphinxsys_containers.h"
#include "vector_functions.h"

#include <fstream>
#include <iostream>
#include <string>

namespace SPH
{

enum Image_Data_Type
{
    MET_FLOAT,
    MET_UCHAR,
    MET_LONG
};

enum Output_Mode
{
    BINARY,
    ASCII
};

template <typename T, int nDims>
class ImageMHD
{
  public:
    ImageMHD() {};
    // constructor for input files
    explicit ImageMHD(std::string full_path_file);
    // constructor for sphere
    ImageMHD(Real radius, Array3i dxdydz, Vec3d spacings);
    ~ImageMHD();

    void set_objectType(const std::string &objectType)
    {
        objectType_ = objectType;
    };

    void set_binaryData(bool binaryData)
    {
        binaryData_ = binaryData;
    };
    void set_binaryDataByteOrderMSB(bool binaryDataByteOrderMSB)
    {
        binaryDataByteOrderMSB_ = binaryDataByteOrderMSB;
    };
    void set_compressedData(bool compressedData)
    {
        compressedData_ = compressedData;
    };
    void set_transformMatrix(Mat3d transformMatrix)
    {
        transformMatrix_ = transformMatrix;
    };
    void set_offset(Vec3d offset)
    {
        offset_ = offset;
    };
    void set_centerOfRotation(Vec3d centerOfRotation)
    {
        centerOfRotation_ = centerOfRotation;
    };
    void set_elementSpacing(Vec3d elementSpacing)
    {
        elementSpacing_ = elementSpacing;
    };
    void set_dimSize(Vec3d dimSize)
    {
        dimSize_ = dimSize;
    };
    void set_anatomicalOrientation(const std::string &anatomicalOrientation)
    {
        anatomicalOrientation_ = anatomicalOrientation;
    };
    void set_elementType(const Image_Data_Type &elementType)
    {
        elementType_ = elementType;
    };
    void set_elementDataFile(std::string elementDataFile)
    {
        elementDataFile_ = elementDataFile;
    };

    T *get_data() { return data_; };

    int get_size() { return size_; }

    Real get_min_value() { return min_value_; };
    Real get_max_value() { return max_value_; };

    Vec3d findClosestPoint(const Vec3d &probe_point);
    BoundingBoxd findBounds();
    Real findValueAtPoint(const Vec3d &probe_point);
    Vec3d findNormalAtPoint(const Vec3d &probe_point);

    void write(std::string filename, Output_Mode = BINARY);

  private:
    std::string objectType_;
    int nDims_;
    bool binaryData_;
    bool binaryDataByteOrderMSB_;
    bool compressedData_;
    Mat3d transformMatrix_;
    Vec3d offset_;
    Vec3d centerOfRotation_;
    Vec3d elementSpacing_;
    Array3i dimSize_;
    int width_;
    int height_;
    int depth_;
    int size_;
    std::string anatomicalOrientation_;
    Image_Data_Type elementType_;
    std::string elementDataFile_;
    Real min_value_;
    Real max_value_;
    T *data_;

    std::vector<int> findNeighbors(const Vec3d &probe_point, Array3i &this_cell);
    Vec3d computeGradientAtCell(int i);
    Vec3d computeNormalAtCell(int i);
    T getValueAtCell(int i);
    Vec3d convertToPhysicalSpace(Vec3d p);
    void split(const std::string &s, char delim, std::vector<std::string> &elems);
};

} // namespace SPH

#include "image_mhd.hpp"

#endif //__EMSCRIPTEN__

#endif // IMAGE_MHD_H