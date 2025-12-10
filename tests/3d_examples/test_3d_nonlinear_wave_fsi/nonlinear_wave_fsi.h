/**
 * @file 	nonlinear_wave_fsi.h
 * @brief 	This is the 3d case file for wave impact with tension leg moored floating structure.
 * @author   Nicolò Salis
 */
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 6.0; /**< TOTAL SIMULATION TIME*/
Real DL = 10.0;                 /**< Tank length. */
Real DH = 1.0;                  /**< Tank height. */
Real DW = 1.6;                  /**< Tank width. */
Real WH = 0.8;                  /**< Water block height. */
Real BEH = 2.0;                 /**< Flume height. */
Real Wmk_p = 0.0;               /**< Wave_maker initial position. */
Real EXS = 2.0;                 /**< Extra space behind the wave_maker*/
Real HWM = 1.0;                 /**< Wave_maker height*/
Real d = 0.8;                   /**< Still water level*/
Real dDx = 10;
Real particle_spacing_ref = d / dDx;           /**< Main water particle spacing. */
Real particle_spacing_structure = 0.1;         /**< Structure particle spacing. */
Real BW = particle_spacing_ref * 4.0;          /**< Extending width for BCs. */
Real Maker_width = particle_spacing_ref * 4.0; /**< Width of the wave_maker. */

BoundingBoxd system_domain_bounds(Vecd(-BW, -EXS - BW, -BW), Vecd(EXS + BW, DL + BW, DH + BW));

Vecd offset = Vecd::Zero();
// water block parameters
Vec2d watS_lb(0.0, 0.0); /**< Left bottom. */
Vec2d watS_lt(0.0, WH);  /**< Left top. */
/** Define the corner points of the gate geometry. */
Vec2d Wmak_lb(Wmk_p - Maker_width, 0.0); /**< Left bottom. */
Vec2d Wmak_lt(Wmk_p - Maker_width, 1.5); /**< Left top. */
Vec2d Wmak_rt(Wmk_p, 1.5);               /**< Right top. */
Vec2d Wmak_rb(Wmk_p, 0.0);               /**< Right bottom. */
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
Real Strx = 0.8;
Real Stry = 5.0;
Real Strz = 0.573;

// LEG
Real Lw = 0.19;
Real Ld = 1.3;
Real Lh = 0.12;

Real Law = Ld * Lh;
Real Lad = Lw * Lh;
Real LVol = Lw * Ld * Lh;

// PILLAR
Real Pw = 0.19;
Real Pd = 0.2;
Real Ph = 0.24;

Real Pad = Pw * Ph;
Real Paw = Pd * Ph;
Real PVol = Pw * Pd * Ph;

// TOP
Real Tw = 0.81;
Real Td = 0.94;
Real Th = 0.11;

Real Tad = Tw * Th;
Real Taw = Td * Th;
Real TVol = Tw * Td * Th;

// STIFFENER
Real Sw = 0.43;
Real Sd = 0.08;
Real Sh = 0.04;

Real Saw = Sd * Sh;
Real Sad = Sw * Sh;
Real SVol = Sw * Sd * Sh;

/**
 * Structure Volume and Density
 */

Real StructureVol = 2 * LVol + 4 * PVol + TVol + 2 * SVol;
Real StructureMass = 62.036;                          /* Weight of the solid structure (Luo share)*/
Real StructureDensity = StructureMass / StructureVol; /* Structure density */

/**
 * Inertia of Structure Parts
 */

// LEG
Real LIx = (StructureDensity * LVol) / 12 * (Ld * Ld + Lh * Lh);
Real LIy = (StructureDensity * LVol) / 12 * (Lw * Lw + Lh * Lh);
Real LIz = (StructureDensity * LVol) / 12 * (Ld * Ld + Lw * Lw);
// PILLAR
Real PIx = (StructureDensity * PVol) / 12 * (Pd * Pd + Ph * Ph);
Real PIy = (StructureDensity * PVol) / 12 * (Pw * Pw + Ph * Ph);
Real PIz = (StructureDensity * PVol) / 12 * (Pd * Pd + Pw * Pw);
// TOP
Real TIx = (StructureDensity * TVol) / 12 * (Td * Td + Th * Th);
Real TIy = (StructureDensity * TVol) / 12 * (Tw * Tw + Th * Th);
Real TIz = (StructureDensity * TVol) / 12 * (Td * Td + Tw * Tw);
// STIFFENER
Real SIx = (StructureDensity * SVol) / 12 * (Sd * Sd + Sh * Sh);
Real SIy = (StructureDensity * SVol) / 12 * (Sw * Sw + Sh * Sh);
Real SIz = (StructureDensity * SVol) / 12 * (Sd * Sd + Sw * Sw);

/*
 * Structure center of mass
 */

/* Legs */
Real dx = 0.215;
Vec3d RLcm(Strx + dx + Lw / 2, Stry + Ld / 2, Strz + Lh / 2); /* Right Leg */
Vec3d LLcm(Strx - dx - Lw / 2, Stry + Ld / 2, Strz + Lh / 2); /* Left Leg */

/* Pillars */
Real yPsP = 0.25;
Real ySsP = 0.85;
Vec3d RPsPcm(Strx + dx + Pw / 2, Stry + yPsP + Pd / 2, Strz + Lh + Ph / 2); /* Right Port side Pillar */
Vec3d LPsPcm(Strx - dx - Pw / 2, Stry + yPsP + Pd / 2, Strz + Lh + Ph / 2); /* Left Port side Pillar */
Vec3d RSsPcm(Strx + dx + Pw / 2, Stry + ySsP + Pd / 2, Strz + Lh + Ph / 2); /* Right Seaside Pillar */
Vec3d LSsPcm(Strx - dx - Pw / 2, Stry + ySsP + Pd / 2, Strz + Lh + Ph / 2); /* Left Seaside Pillar */

/* StiffEners */
Real zS = 0.14;
Real ySsS = 0.31;
Real yPsS = 0.91;
Vec3d PsScm(Strx, Stry + yPsS + Sd / 2, Strz + zS + Sh / 2); /* Port side StiffEner */
Vec3d SsScm(Strx, Stry + ySsS + Sd / 2, Strz + zS + Sh / 2); /* Seaside StiffEner */

/* Top Plate */
Real yT = 0.18;
Real zT = 0.36;
Vec3d Tcm(Strx, Stry + yT + Td / 2, Strz + zT + Th / 2); /* Top Plate */

Real bcmx = (StructureDensity * LVol * RLcm[0] +
             StructureDensity * LVol * LLcm[0] +
             StructureDensity * PVol * RSsPcm[0] +
             StructureDensity * PVol * RPsPcm[0] +
             StructureDensity * PVol * LSsPcm[0] +
             StructureDensity * PVol * LPsPcm[0] +
             StructureDensity * SVol * PsScm[0] +
             StructureDensity * SVol * SsScm[0] +
             StructureDensity * TVol * Tcm[0]) /
            (StructureDensity * StructureVol);

Real bcmy = (StructureDensity * LVol * RLcm[1] +
             StructureDensity * LVol * LLcm[1] +
             StructureDensity * PVol * RSsPcm[1] +
             StructureDensity * PVol * RPsPcm[1] +
             StructureDensity * PVol * LSsPcm[1] +
             StructureDensity * PVol * LPsPcm[1] +
             StructureDensity * SVol * PsScm[1] +
             StructureDensity * SVol * SsScm[1] +
             StructureDensity * TVol * Tcm[1]) /
            (StructureDensity * StructureVol);

Real bcmz = (StructureDensity * LVol * RLcm[2] +
             StructureDensity * LVol * LLcm[2] +
             StructureDensity * PVol * RSsPcm[2] +
             StructureDensity * PVol * RPsPcm[2] +
             StructureDensity * PVol * LSsPcm[2] +
             StructureDensity * PVol * LPsPcm[2] +
             StructureDensity * SVol * PsScm[2] +
             StructureDensity * SVol * SsScm[2] +
             StructureDensity * TVol * Tcm[2]) /
            (StructureDensity * StructureVol);

Vecd G(bcmx, bcmy, bcmz);

/**
 * Inertia of the body
 */

/* squared distance of body parts in x direction */
Real dLX = (G[1] - RLcm[1]) * (G[1] - RLcm[1]) +
           (G[2] - RLcm[2]) * (G[2] - RLcm[2]);
Real dPX = (G[1] - RPsPcm[1]) * (G[1] - RPsPcm[1]) +
           (G[2] - RPsPcm[2]) * (G[2] - RPsPcm[2]);
Real dSX = (G[1] - PsScm[1]) * (G[1] - PsScm[1]) +
           (G[2] - PsScm[2]) * (G[2] - PsScm[2]);
Real dTX = (G[1] - Tcm[1]) * (G[1] - Tcm[1]) +
           (G[2] - Tcm[2]) * (G[2] - Tcm[2]);

/* squared distance of body parts in y direction */
Real dLY = (G[0] - RLcm[0]) * (G[0] - RLcm[0]) +
           (G[2] - RLcm[2]) * (G[2] - RLcm[2]);
Real dPY = (G[0] - RPsPcm[0]) * (G[0] - RPsPcm[0]) +
           (G[2] - RPsPcm[2]) * (G[2] - RPsPcm[2]);
Real dSY = (G[0] - PsScm[0]) * (G[0] - PsScm[0]) +
           (G[2] - PsScm[2]) * (G[2] - PsScm[2]);
Real dTY = (G[0] - Tcm[0]) * (G[0] - Tcm[0]) +
           (G[2] - Tcm[2]) * (G[2] - Tcm[2]);

/* squared distance of body parts in z direction */
Real dLZ = (G[0] - RLcm[0]) * (G[0] - RLcm[0]) +
           (G[1] - RLcm[1]) * (G[1] - RLcm[1]);
Real dPZ = (G[0] - RPsPcm[0]) * (G[0] - RPsPcm[0]) +
           (G[1] - RPsPcm[1]) * (G[1] - RPsPcm[1]);
Real dSZ = (G[0] - PsScm[0]) * (G[0] - PsScm[0]) +
           (G[1] - PsScm[1]) * (G[1] - PsScm[1]);
Real dTZ = (G[0] - Tcm[0]) * (G[0] - Tcm[0]) +
           (G[1] - Tcm[1]) * (G[1] - Tcm[1]);

Real Ix = 2 * LIx + 2 * (dLX * LVol * StructureDensity) +
          4 * PIx + 4 * (dPX * PVol * StructureDensity) +
          2 * SIx + 2 * (dSX * SVol * StructureDensity) +
          TIx + (dTX * TVol * StructureDensity);

Real Iy = 2 * LIy + 2 * (dLY * LVol * StructureDensity) +
          4 * PIy + 4 * (dPY * PVol * StructureDensity) +
          2 * SIy + 2 * (dSY * SVol * StructureDensity) +
          TIy + (dTY * TVol * StructureDensity);

Real Iz = 2 * LIz + 2 * (dLZ * LVol * StructureDensity) +
          4 * PIz + 4 * (dPZ * PVol * StructureDensity) +
          2 * SIz + 2 * (dSZ * SVol * StructureDensity) +
          TIz + (dTZ * TVol * StructureDensity);

/**
 * Topology of the tethers.
 * */

// Right seaside cable topology
Vecd ground_tethering_AR(G[0] + 0.31, G[1] - 0.3, 0.0);
Vecd structure_tethering_AR(G[0] + 0.31, G[1] - 0.3, Strz);
// Left seaside cable topology
Vecd ground_tethering_AL(G[0] - 0.31, G[1] - 0.3, 0.0);
Vecd structure_tethering_AL(G[0] - 0.31, G[1] - 0.3, Strz);
// Right port side cable topology
Vecd ground_tethering_BR(G[0] + 0.31, G[1] + 0.3, 0.0);
Vecd structure_tethering_BR(G[0] + 0.31, G[1] + 0.3, Strz);
// Left port side cable topology
Vecd ground_tethering_BL(G[0] - 0.31, G[1] + 0.3, 0.0);
Vecd structure_tethering_BL(G[0] - 0.31, G[1] + 0.3, Strz);

Real cable_length = Strz;

/**
 * Structure observer position
 * */

Vecd obs(Strx, Stry + yT + Td / 2, Strz + zT + Th);

//------------------------------------------------------------------------------
// geometric shape elements used in the case
//------------------------------------------------------------------------------
/** Set the file path to the stl file. */
std::string stl_structure_path = "./input/structure_cm.stl";
Real StructureScale = 1;
Vecd translation_str = G;

class FloatingStructure : public ComplexShape
{
  public:
    explicit FloatingStructure(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(stl_structure_path, translation_str, StructureScale);
    }
};

class StructureSystemForSimbody : public SolidBodyPartForSimbody
{
  public:
    StructureSystemForSimbody(SPHBody &sph_body, SharedPtr<Shape> shape_ptr)
        : SolidBodyPartForSimbody(sph_body, shape_ptr)
    {
        // Vecd mass_center(G[0], G[1], G[2]);
        // initial_mass_center_ = SimTK::Vec3(mass_center[0], mass_center[1], mass_center[2]);
        body_part_mass_properties_ =
            mass_properties_ptr_keeper_
                .createPtr<SimTK::MassProperties>(StructureMass, SimTK::Vec3(0.0), SimTK::UnitInertia(Ix, Iy, Iz));
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
        Vecd halfsize_water(0.5 * DW, 0.5 * (DL - EXS), 0.5 * WH);
        Vecd water_pos(0.5 * DW, 0.5 * (DL - EXS), 0.5 * WH);
        Transform translation_water(water_pos);
        add<GeometricShapeBox>(Transform(translation_water), halfsize_water);
        subtract<TriangleMeshShapeSTL>(stl_structure_path, translation_str, StructureScale);
    }
};

//----------------------------------------------------------------------
//	create a wave_maker shape
//----------------------------------------------------------------------
Vecd wave_maker_shape(0.5 * DW, 0.5 * Maker_width, 0.5 * DH);
Vecd wmk_pos(0.5 * DW, -0.5 * Maker_width, 0.5 * HWM);
Transform translation_wave_maker(wmk_pos);

//----------------------------------------------------------------------
//	Wall geometries.
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_wall_outer(0.5 * DW + BW, 0.5 * DL + BW, 0.5 * DH + BW);
        Vecd wall_outer_pos(0.5 * DW, 0.5 * DL - EXS, 0.5 * DH);
        Transform translation_wall_outer(wall_outer_pos);
        add<GeometricShapeBox>(Transform(translation_wall_outer), halfsize_wall_outer);

        Vecd halfsize_wall_inner(0.5 * DW, 0.5 * DL, 0.5 * DH + BW);
        Vecd wall_inner_pos(0.5 * DW, 0.5 * DL - EXS, BW + 0.5 * DH);
        Transform translation_wall_inner(wall_inner_pos);
        subtract<GeometricShapeBox>(Transform(translation_wall_inner), halfsize_wall_inner);

        add<GeometricShapeBox>(Transform(translation_wave_maker), wave_maker_shape);
    }
};

//----------------------------------------------------------------------
//	Boundary condition for wave_maker
//----------------------------------------------------------------------
class WaveMaking : public BodyPartMotionConstraint
{
    Real h;
    Real tf;
    Real xf;
    Real fmn;
    Real fmx;
    Real a;
    int N;
    std::vector<Real> om;
    std::vector<Real> S;
    std::vector<Real> k;
    Real g;

    Vec3d getDisplacement(const Real &time)
    {
        Real dp = 0;
        for (int jj = 0; jj < (N); jj++)
        {
            dp = dp + 0.5 * S[jj] * cos(-k[jj] * xf - om[jj] * (time - tf));
        };

        Vec3d displacement{Vecd::Zero()};
        displacement[1] = dp;

        return displacement;
    }

    Vec3d getVelocity(const Real &time)
    {
        Real vl = 0;
        for (int jj = 0; jj < (N); jj++)
        {
            vl = vl + 0.5 * om[jj] * S[jj] * sin(-k[jj] * xf - om[jj] * (time - tf));
        };
        Vec3d velocity{Vecd::Zero()};
        velocity[1] = vl;
        return velocity;
    }

    Vec3d getAcceleration(const Real &time)
    {
        Real ax = 0;
        for (int jj = 0; jj < (N); jj++)
        {
            ax = ax - 0.5 * om[jj] * om[jj] * S[jj] * cos(-k[jj] * xf - om[jj] * (time - tf));
        };
        Vec3d acceleration{Vecd::Zero()};
        acceleration[1] = ax;
        return acceleration;
    }

    Real OBJ(Real wnmb, Real omsq)
    {
        return omsq - g * wnmb * tanh(wnmb * h);
    };

    void ComputeWaveChar()
    {
        std::vector<Real> f(N);

        for (int i = 0; i < (N); i++)
        {
            f[i] = fmn + i * (fmx - fmn) / N;
        };

        om = f;
        for (int i = 0; i < (N); i++)
        {
            om[i] = 2 * PI * f[i];
        };

        k = om;
        std::vector<Real> OBJ_(N);
        OBJ_ = om;
        for (int i = 0; i < (N); i++)
        {
            // Solve dispersion equation
            Real omsq = om[i] * om[i];
            Real kmin = 0;
            Real kmax = 20;

            Real tol = 1E-15;
            Real wnmb = kmin;
            while ((kmax - kmin) >= tol)
            {

                wnmb = (kmin + kmax) / 2;
                OBJ_[i] = OBJ(wnmb, omsq);

                Real OBJmin = OBJ(kmin, omsq);
                Real OBJmax = OBJ(wnmb, omsq);

                if (abs(OBJ_[i]) <= tol)
                {
                    break;
                }
                else if (OBJmin * OBJmax < 0)
                {
                    kmax = wnmb;
                }
                else
                {
                    kmin = wnmb;
                };
            };

            k[i] = wnmb;
        };

        S = f;
        for (int i = 0; i < (N); i++)
        {
            S[i] = a * (sinh(k[i] * h) * cosh(k[i] * h) + k[i] * h) / (sinh(k[i] * h) * sinh(k[i] * h));
        };
    };

  public:
    WaveMaking(BodyPartByParticle &body_part)
        : BodyPartMotionConstraint(body_part),
          h(WH), tf(5), xf(4.5), fmn(0.32), fmx(0.96), a(0.0078), N(32), g(gravity_g),
          acc_(particles_->registerStateVariableData<Vecd>("Acceleration")),
          physical_time_(sph_system_->getSystemVariableDataByName<Real>("PhysicalTime"))

    {
        ComputeWaveChar();
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        Real time = *physical_time_;
        pos_[index_i] = pos0_[index_i] + getDisplacement(time);
        vel_[index_i] = getVelocity(time);
        acc_[index_i] = getAcceleration(time);
    };

  protected:
    Vecd *acc_;
    Real *physical_time_;
};

//----------------------------------------------------------------------
//	create measuring probes
//----------------------------------------------------------------------
Real h = 1.3 * particle_spacing_ref;
Vecd WGaugeDim(0.5 * DW, 0.5 * h, 0.5 * DH);
Vecd WGauge(0.0, 4., 0.5 * DH);
Transform translation_WGauge(WGauge);
