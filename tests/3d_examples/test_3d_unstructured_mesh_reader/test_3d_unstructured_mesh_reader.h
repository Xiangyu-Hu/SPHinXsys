/**
 * @file 	test_3d_unstructured_mesh_reader.h
 * @brief 	This is a test to show the mesh parser for ICEM and Fluent .msh files.
 * @author 	Zhentong Wang and Xiangyu Hu
 */

#ifndef TEST_3D_UNSTRUCTURED_MESH_READER_H
#define TEST_3D_UNSTRUCTURED_MESH_READER_H
#include "unstructured_mesh_3d.h"             
using namespace SPH;
using namespace std;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 4.0;                                               /**< Computation domain length. */
Real DH = 1.0;                                              /**< Computation domain height. */
Real DW = 0.1;                                              /**< Computation domain width. */
Real particle_spacing_ref = 1.0 / 240.0;                /**< Initial reference particle spacing. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec3d(0.0, 0.0, 0.0), Vec3d(DL, DH, DW));
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho = 1.0;                        /**< initial density of another. */
Real u = 1.0;                           /**< initial velocity of another in X axis. */
Real v = 0.0;                           /**< initial velocity of another in Y axis. */
Real w = 0.0;                          /**< initial velocity of another in Z axis. */
Real p = 140.2 / 1.2;                /**< initial pressure of another. */
Real heat_capacity_ratio = 1.4;              /**< heat capacity ratio. */

//----------------------------------------------------------------------
//	Set the file path to the data file.
//----------------------------------------------------------------------
std::string mesh_fullpath = "./input/3D_ICEM_MESH.msh";
//	Define geometries and body shapes
//----------------------------------------------------------------------
class AirBody : public ComplexShape
{
public:
    explicit AirBody(const std::string& shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wave(0.5 * DH, 0.5 * DL, 0.5 * DW);
        Transform translation_wave(halfsize_wave);
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wave), halfsize_wave);
    }
};

#endif // TEST_3D_UNSTRUCTURED_MESH_READER_H