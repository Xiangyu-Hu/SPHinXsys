/**
 * @file 	 stfb.h
 * @brief 	 This is the case file for 3D still floaing body.
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
Real L = 1.0;                   /**< Base of the floating body. */
Real particle_spacing_ref = L / 10;
Real BW = particle_spacing_ref * 4.0;          /**< Extending width for BCs. */
Real Maker_width = particle_spacing_ref * 4.0; /**< Width of the wavemaker. */
BoundingBox system_domain_bounds(Vecd(-BW, -BW, -BW), Vecd(DW + BW, DL + BW, DH + BW));
Vecd offset = Vecd::Zero();
//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                  /**< Reference density of fluid. */
Real gravity_g = 9.81;                 /**< Value of gravity. */
Real U_f = 2.0 * sqrt(WH * gravity_g); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                 /**< Reference sound speed. */
Real mu_f = 1.0e-3;
//----------------------------------------------------------------------
//	Structure Properties G and Inertia
//----------------------------------------------------------------------
/* Weight of the solid structure*/
Real StructureMass = 700;
/**< Area of the solid structure*/
Real FlStA = L * L * L;
/**< Density of the solid structure*/
Real rho_s = StructureMass / FlStA;
/* Equilibrium position of the solid structure*/
Real H = WH - (rho_s / rho0_f * L - L / 2); /**< Strart placemnt of Flt Body*/

Real bcmx = DL / 2;
Real bcmy = DL / 2;
Real bcmz = H;
Vecd G(bcmx, bcmy, bcmz);
Real Ix = StructureMass / 12 * (L * L + L * L);
Real Iy = StructureMass / 12 * (L * L + L * L);
Real Iz = StructureMass / 12 * (L * L + L * L);
/**
 *
 * Structure observer position
 *
 * */
Vecd obs = G;

/** Geometry definition. */
Vecd halfsize_structure(0.5 * L, 0.5 * L, 0.5 * L);
Vecd structure_pos(G[0], G[1], G[2]);
Transform translation_str(structure_pos);

class FloatingStructure : public ComplexShape
{
  public:
    explicit FloatingStructure(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TransformShape<GeometricShapeBox>>(Transform(translation_str), halfsize_structure);
    }
};

class StructureSystemForSimbody : public SolidBodyPartForSimbody
{
  public:
    StructureSystemForSimbody(SPHBody &sph_body, SharedPtr<Shape> shape_ptr)
        : SolidBodyPartForSimbody(sph_body, shape_ptr)
    {
        // Vec2d mass_center(G[0], G[1]);
        // initial_mass_center_ = SimTKVec3(mass_center[0], mass_center[1], 0.0);
        body_part_mass_properties_ =
            mass_properties_ptr_keeper_
                .createPtr<SimTK::MassProperties>(StructureMass, SimTKVec3(0.0), SimTK::UnitInertia(Ix, Iy, Iz));
    }
};
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
        add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_str), halfsize_structure);
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
        add<TransformShape<GeometricShapeBox>>(Transform(translation_wall_outer), halfsize_wall_outer);

        Vecd halfsize_wall_inner(0.5 * DW, 0.5 * DL, 0.5 * DH + BW);
        Vecd wall_inner_pos(0.5 * DW, 0.5 * DL, 0.5 * DH + BW);
        Transform translation_wall_inner(wall_inner_pos);
        subtract<TransformShape<GeometricShapeBox>>(Transform(translation_wall_inner), halfsize_wall_inner);
    }
};
//----------------------------------------------------------------------
//	create mesuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
Vecd FS_gaugeDim(0.5 * h, 0.5 * h, 0.5 * DH);
Vecd FS_gauge(DW / 3, DL / 3, 0.5 * DH);
Transform translation_FS_gauge(FS_gauge);

/**
 * @class FreeSurfaceHeightZ
 * @brief Probe the free surface profile for a fluid body part by reduced operation.
 */
class FreeSurfaceHeightZ : public BaseLocalDynamicsReduce<Real, ReduceMax, BodyPartByCell>,
                           public SPH::fluid_dynamics::FluidDataSimple
{
  protected:
    StdLargeVec<Vecd> &pos_;

  public:
    FreeSurfaceHeightZ(BodyPartByCell &body_part)
        : BaseLocalDynamicsReduce<Real, ReduceMax, BodyPartByCell>(body_part, Real(MinRealNumber)),
          SPH::fluid_dynamics::FluidDataSimple(sph_body_), pos_(particles_->pos_)
    {
        quantity_name_ = "FreeSurfaceHeight";
    }
    virtual ~FreeSurfaceHeightZ(){};

    Real reduce(size_t index_i, Real dt = 0.0) { return pos_[index_i][2]; };
};