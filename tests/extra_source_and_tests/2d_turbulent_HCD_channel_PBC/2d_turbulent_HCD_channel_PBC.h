/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "sphinxsys.h"
#include "k-epsilon_turbulent_model.cpp"
using namespace SPH;


/**
 * @ Parameters for validation with [1997 Joseph JCP].
 */
/*Real DH = 1.0e-3;
Real DL = 4.0 * DH;
Real rho0_f = 1000.0;
Real gravity_g = 1.0e-4;
Real nu_f = 1.0e-6;
Real mu_f = nu_f * rho0_f;
Real U_f = 1.25e-5;
Real c_f = 10.0 * U_f;*/             

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 2.0;                         /**< Channel height. */
Real DL = 6.0;                         /**< Channel length. */
Real resolution_ref = DH / 40;              /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
Real DL1 = 1.0;
Real DL2 = 1.5;
Real DL3 = 1.0;
Real DH1 = DL2 / sqrt(3.0);
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
StdVec<int> id_exclude;
Real y_p_theo = 0.0;
Real offset_dist_ref = 0.0;
//** Intial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = { 0.01 ,0.1 ,0.001 };
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------

//Real U_max = 1.0;
//Real U_f = 1.0; //*Characteristic velo is regarded as average velo here
//Real c_f = 10.0 * U_max;                                        /**< Speed of sound. */
//Real rho0_f = 1.0;                                            /**< Density. */
//Real mu_f = 0.0001;
//Real gravity_g = 0.00872269;

//** Compare with FVM3 *
//Real U_max = 1.0;
//Real U_f = 1.0; //*Characteristic velo is regarded as average velo here
//Real c_f = 10.0 * U_max;                                        /**< Speed of sound. */
//Real rho0_f = 1.0;                                            /**< Density. */
//Real mu_f = 0.0001;
//Real gravity_g = 0.01079379;

//** Compare with FVM4 *
Real U_max = 1.1;
Real U_f = 1.0; //*Characteristic velo is regarded as average velo here
Real c_f = 10.0 * U_max;                                        /**< Speed of sound. */
Real rho0_f = 1.0;                                            /**< Density. */
Real mu_f = 0.0005;
Real gravity_g = 0.01806807;

Real Re = U_f * DH * rho0_f / mu_f;

//----------------------------------------------------------------------
// Observation.
//----------------------------------------------------------------------

Real pos_x_cross_line = 1.25 ;
Real size_cell_x = 1.8 * resolution_ref;
StdVec<Real> monitoring_bound = { pos_x_cross_line - size_cell_x/2.0, pos_x_cross_line + size_cell_x / 2.0 }; //**The 2nd periodic wave*
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));

    water_block_shape.push_back(Vecd(DL - DL1, 0.0));
    water_block_shape.push_back(Vecd(DL1 + DL2 + DL3, DH1));
    water_block_shape.push_back(Vecd(DL1 + DL2, DH1));
    water_block_shape.push_back(Vecd(DL1, 0.0));

    water_block_shape.push_back(Vecd(0.0, 0.0));

    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createWaterBlockShape());
        add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
    }
};

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public MultiPolygonShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(-BW, -BW));
        outer_wall_shape.push_back(Vecd(-BW, DH + BW));
        outer_wall_shape.push_back(Vecd(DL + BW, DH + BW));
        outer_wall_shape.push_back(Vecd(DL + BW, -BW));
        outer_wall_shape.push_back(Vecd(-BW, -BW));
        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(-2.0 * BW, 0.0));
        inner_wall_shape.push_back(Vecd(-2.0 * BW, DH));
        inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, DH));
        inner_wall_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));

        inner_wall_shape.push_back(Vecd(DL - DL1, 0.0));
        inner_wall_shape.push_back(Vecd(DL1 + DL2 + DL3, DH1));
        inner_wall_shape.push_back(Vecd(DL1 + DL2, DH1));
        inner_wall_shape.push_back(Vecd(DL1, 0.0));

        inner_wall_shape.push_back(Vecd(-2.0 * BW, 0.0));

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
    }
};