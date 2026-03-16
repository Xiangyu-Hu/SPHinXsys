/**
 * @file 	 stfb.h
 * @brief 	 This is the case file for 2D still floating body.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 10.0; /**< TOTAL SIMULATION TIME*/
Real DL = 3.0;                   /**< Tank length. */
Real DH = 4.0;                   /**< Tank height. */
Real WH = 2.0;                   /**< Water block height. */
Real L = 1.0;                    /**< Base of the floating body. */
Real particle_spacing_ref = L / 20;
Real BW = particle_spacing_ref * 4.0; /**< Extending width for BCs. */
BoundingBoxd system_domain_bounds(Vec2d(-DL - BW, -DH - BW), Vec2d(DL + BW, DH + BW));
Vec2d offset = Vec2d::Zero();
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                    /**< Reference density of fluid. */
Real gravity_g = 9.81;                   /**< Value of gravity. */
Real U_f = 2.0 * sqrt(0.79 * gravity_g); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                   /**< Reference sound speed. */
Real mu_f = 1.0e-3;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d water_block_halfsize = Vec2d(0.5 * DL, 0.5 * WH);
Vec2d water_block_translation = Vec2d(0.0, -0.5 * WH);
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(0.0, 0.0);
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = Vec2d(0.0, 0.0);
//----------------------------------------------------------------------
//	Structure Properties G and Inertia
//----------------------------------------------------------------------
/* Weight of the solid structure*/
Real StructureMass = 700;
/**< Area of the solid structure*/
Real FlStA = L * L;
/**< Density of the solid structure*/
Real rho_s = StructureMass / FlStA;
/* Equilibrium position of the solid structure*/
Real H = -(rho_s / rho0_f * L - L / 2); /**< Strart placement of Flt Body*/

Real bcmx = 0;
Real bcmy = H + 0;
Vec2d G(bcmx, bcmy);
Real Ix = L * L * L * L / 3;
Real Iy = L * L * L * L / 3;
Real Iz = StructureMass / 12 * (L * L + L * L);

/** Structure observer position*/
Vec2d obs = G;
/** Structure definition*/
Vec2d structure_halfsize = Vec2d(0.5 * L, 0.5 * L);
Vec2d structure_translation = Vec2d(0.0, H);
//------------------------------------------------------------------------------
// geometric shape elements used in the case
//------------------------------------------------------------------------------

class StructureSystemForSimbody : public SolidBodyPartForSimbody
{
  public:
    StructureSystemForSimbody(SPHBody &sph_body, Shape &shape)
        : SolidBodyPartForSimbody(sph_body, shape)
    {
        // Vec2d mass_center(G[0], G[1]);
        // initial_mass_center_ = SimTKVec3(mass_center[0], mass_center[1], 0.0);
        body_part_mass_properties_ =
            mass_properties_keeper_
                .createPtr<SimTK::MassProperties>(StructureMass, SimTKVec3(0.0), SimTK::UnitInertia(Ix, Iy, Iz));
    }
};
//----------------------------------------------------------------------
//	Dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(water_block_translation), water_block_halfsize);
        subtract<GeometricShapeBox>(Transform(structure_translation), structure_halfsize);
    }
};
//----------------------------------------------------------------------
//	create measuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
Vec2d gauge_halfsize = Vec2d(0.5 * h, 0.5 * DH);
Vec2d gauge_translation = Vec2d(DL / 3, 0.5 * DH);