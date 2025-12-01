/**
 * @file 	 stlw.h
 * @brief 	 This is the case file for still water.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 5.0; /**< TOTAL SIMULATION TIME*/
Real DW = 3.0;                  /**< Water length. */
Real DL = 3.0;                  /**< Tank length. */
Real DH = 2.5;                  /**< Tank height. */
Real WH = 2.0;                  /**< Water block height. */
Real particle_spacing_ref = 0.1;
Real BW = particle_spacing_ref * 4.0;          /**< Extending width for BCs. */
Real Maker_width = particle_spacing_ref * 4.0; /**< Width of the wave_maker. */

BoundingBoxd system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DW + BW, DL + BW, DH + BW));

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                  /**< Reference density of fluid. */
Real gravity_g = 9.81;                 /**< Value of gravity. */
Real U_f = 2.0 * sqrt(WH * gravity_g); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                 /**< Reference sound speed. */
Real mu_f = 1.0e-3;
//----------------------------------------------------------------------
//	Water block
//----------------------------------------------------------------------
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        /** Geometry definition. */
        Vecd halfsize_water(0.5 * DW, 0.5 * DL, 0.5 * WH);
        Vecd water_pos(0.5 * DW, 0.5 * DL, 0.5 * WH);
        Transform translation_water(water_pos);
        add<GeometricShapeBox>(Transform(translation_water), halfsize_water);
    }
};
//----------------------------------------------------------------------
//	Wall geometries.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wall_outer(0.5 * DW + BW, 0.5 * DL + BW, 0.5 * DH + BW);
        Vecd wall_outer_pos(0.5 * DW, 0.5 * DL, 0.5 * DH);
        Transform translation_wall_outer(wall_outer_pos);
        add<GeometricShapeBox>(Transform(translation_wall_outer), halfsize_wall_outer);

        Vecd halfsize_wall_inner(0.5 * DW, 0.5 * DL, 0.5 * DH + BW);
        Vecd wall_inner_pos(0.5 * DW, 0.5 * DL, 0.5 * DH + BW);
        Transform translation_wall_inner(wall_inner_pos);
        subtract<GeometricShapeBox>(Transform(translation_wall_inner), halfsize_wall_inner);
    }
};
//----------------------------------------------------------------------
//	create measuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
Vecd FS_gaugeDim(0.5 * h, 0.5 * h, 0.5 * DH);
Vecd FS_gauge(DW / 3, DL / 3, 0.5 * DH);
Transform translation_FS_gauge(FS_gauge);
