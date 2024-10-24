/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "k-epsilon_turbulent_model.cpp"
#include "sphinxsys.h"
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
Real DH = 1.0; /**< Channel height. */
Real DL = 4.0; /**< Channel length. */
Real amplitude = 0.1;
Real wave_length = 1;
Real resolution_ref = DH / 20.0; /**< Initial reference particle spacing. */
Real BW = resolution_ref * 4;    /**< Reference size of the emitter. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-2.0 * BW - wave_length / 2.0, -DH), Vec2d(DL + 2.0 * BW + wave_length / 2.0, DH + 2.0 * BW));
StdVec<int> id_exclude;
Real y_p_theo = 0.0;
Real offset_dist_ref = 0.0;
//** Initial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = {0.01, 0.1, 0.001};
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------

//Real u_max = 1.5;
//Real U_f = u_max * 2.0 / 3.0; //*Characteristic velo is regarded as average velo here
Real U_f = 0.816 * 1.5; //*Characteristic velo is regarded as average velo here
Real c_f = 10.0 * U_f;  /**< Speed of sound. */

Real rho0_f = 1.0; /**< Density. */

//Real Re = 40000.0;
//Real mu_f = rho0_f * U_f * DH  / Re;
Real mu_f = 0.0001;
Real gravity_g = 0.01874594;
Real Re = U_f * DH * rho0_f / mu_f;
//Real mu_f = rho0_f * U_f * DH / Re; /**< Dynamics viscosity. */
//Real mu_f = rho0_f * u_max * DH / Re; /**< Dynamics viscosity. */
//Real gravity_g = 2.0 * mu_f * u_max / rho0_f / (DH / 2.0) / (DH / 2.0);

//----------------------------------------------------------------------
// Observation.
//----------------------------------------------------------------------

//** Trough *
int num_observer_points = std::round((DH + amplitude) / resolution_ref); //**Every particle is regarded as a cell monitor*
Real pos_x_cross_line = 1.25;
Real size_cell_x = 1.8 * resolution_ref;
StdVec<Real> monitoring_bound = {pos_x_cross_line - size_cell_x / 2.0, pos_x_cross_line + size_cell_x / 2.0}; //**The 2nd periodic wave*
Real observe_offset_y = -1.0 * amplitude;
//----------------------------------------------------------------------
//	Cases-dependent geometries
//----------------------------------------------------------------------
/** the water block . */
std::vector<Vecd> water_block_shape{
    Vecd(0.0, 0.0), Vecd(0.0, DH), Vecd(DL, DH), Vecd(DL, 0.0), Vecd(0.0, 0.0)};
/** the water in the trough part. */
std::vector<Vecd> water_shape_trough{
    Vecd(0.0, 0.0),
    Vecd(DL, 0.0),
    Vecd(DL, -amplitude),
    Vecd(0.0, -amplitude),
    Vecd(0.0, 0.0),
};
/** the top wall polygon. */
std::vector<Vecd> top_wall_shape{
    Vecd(0.0 - 2.0 * BW, DH),
    Vecd(0.0 - 2.0 * BW, DH + BW),
    Vecd(DL + 2.0 * BW, DH + BW),
    Vecd(DL + 2.0 * BW, DH),
    Vecd(0.0 - 2.0 * BW, DH),
};
/** the bottom wall polygon. */
std::vector<Vecd> bottom_wall_shape{};
class WavyShape
{
    int num_wavy_points, num_wavy_extend_left, num_wavy_extend_right;
    int num_total;
    Real wavy_interval = 0.01; //** Specify the wavy resolution here *
  public:
    explicit WavyShape()
    {
        //bottom_wall_shape.push_back(Vecd(0.0 - 2.0 * BW, 0.0));
        //bottom_wall_shape.push_back(Vecd(0.0, 0.0));
        num_wavy_extend_left = wave_length / 2.0 / wavy_interval;
        num_wavy_extend_right = wave_length / 2.0 / wavy_interval;
        num_wavy_points = (DL - 0.0) / wavy_interval;
        //num_total = num_wavy_extend_left + num_wavy_points + num_wavy_extend_right;
        for (int k = -num_wavy_extend_left; k <= num_wavy_points + num_wavy_extend_right; ++k)
        {
            Real x_coordinate = k * wavy_interval;
            Real y_coordinate = -1.0 * amplitude * std::sin(2.0 * M_PI * x_coordinate);
            bottom_wall_shape.push_back(Vecd(x_coordinate, y_coordinate));
        }
        Real x_right = (num_wavy_points + num_wavy_extend_right) * wavy_interval;
        Real x_left = -num_wavy_extend_left * wavy_interval;
        Real y_left = -1.0 * amplitude * std::sin(2.0 * M_PI * x_left);

        //bottom_wall_shape.push_back(Vecd(DL, 0.0));
        //bottom_wall_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
        bottom_wall_shape.push_back(Vecd(x_right, -BW - amplitude));
        bottom_wall_shape.push_back(Vecd(x_left, -BW - amplitude));
        bottom_wall_shape.push_back(Vecd(x_left, y_left));
    }
};

//----------------------------------------------------------------------
// @brief 	Fluid body definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(water_shape_trough, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(bottom_wall_shape, ShapeBooleanOps::sub);
    }
};

//----------------------------------------------------------------------
// @brief 	Wall boundary body definition.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        multi_polygon_.addAPolygon(top_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(bottom_wall_shape, ShapeBooleanOps::add);
    }
};