/**
 * @file 	wfsi.h
 * @brief 	This is the 2d case header for wave impact with tension leg floating structure.
 * @author  Nicolò Salis
 */
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real total_physical_time = 20.0; /**< TOTAL SIMULATION TIME*/
Real DL = 35.0;                  /**< Tank length. */
Real DH = 2.0;                   /**< Tank height. */
Real WL = 20.0;                  /**< Water block width. */
Real WH = 0.8;                   /**< Water block height. */
Real TB = 15.0;                  /**< Beach start. */
Real DB = 25.0;                  /**< Beach end. */
Real BEH = 2.0;                  /**< Beach end height. */
Real Wmk_p = 0.0;                /**< Wavemaker initial position. */
Real EXS = 2.0;                  /**< etra space behind the wavemaker*/
Real StructureBasePlateH = 0.12;
/**
 * Initial reference particle spacing
 * It is a multiple of the structure baseplate height.
 * */
Real particle_spacing_ref = StructureBasePlateH / 4;
Real BW = particle_spacing_ref * 4.0;          /**< Extending width for BCs. */
Real Maker_width = particle_spacing_ref * 4.0; /**< Width of the wavemaker. */

BoundingBoxd system_domain_bounds(Vec2d(-EXS - BW, -BW), Vec2d(DL + BW, DH + BW));

Vec2d offset = Vec2d::Zero();

// water block parameters
Vec2d watS_lb(0.0, 0.0); /**< Left bottom. */
Vec2d watS_lt(0.0, WH);  /**< Left top. */
/** Define the corner points of the gate geometry. */
Vec2d Wmak_lb(Wmk_p - Maker_width, 0.0); /**< Left bottom. */
Vec2d Wmak_lt(Wmk_p - Maker_width, 1.5); /**< Left top. */
Vec2d Wmak_rt(Wmk_p, 1.5);               /**< Right top. */
Vec2d Wmak_rb(Wmk_p, 0.0);               /**< Right bottom. */
/* BEACH */
Vec2d b1a(TB, 0);              /**< beach start. */
Vec2d bwh((WH) * 10 + TB, WH); /**< water height at beach. */
Vec2d b1b(25, 1);              /**< beach end. */

/*
 * Buoyant Structure Shape
 */

/** baseplate. */
Real bp_x = 12.286;
Real bp_y = 0.573;
Real bp_l = 1.3;
Real bp_h = 0.12;

Vec2d bp_lb(bp_x, bp_y);               /**< Left bottom. */
Vec2d bp_lt(bp_x, bp_y + bp_h);        /**< Left top. */
Vec2d bp_rt(bp_x + bp_l, bp_y + bp_h); /**< Right top. */
Vec2d bp_rb(bp_x + bp_l, bp_y);        /**< Right bottom. */

/** Seaside pillar . */
Real ssp_x = bp_x + 0.25;
Real ssp_y = bp_y + bp_h;
Real ssp_l = 0.2;
Real ssp_h = 0.24;

Vec2d ssp_lb(ssp_x, ssp_y);                 /**< Left bottom. */
Vec2d ssp_lt(ssp_x, ssp_y + ssp_h);         /**< Left top. */
Vec2d ssp_rt(ssp_x + ssp_l, ssp_y + ssp_h); /**< Right top. */
Vec2d ssp_rb(ssp_x + ssp_l, ssp_y);         /**< Right bottom. */

/** Portside pillar . */
Real psp_x = bp_x + bp_l - 0.45;
Real psp_y = bp_y + bp_h;
Real psp_l = 0.2;
Real psp_h = 0.24;

Vec2d psp_lb(psp_x, psp_y);                 /**< Left bottom. */
Vec2d psp_lt(psp_x, psp_y + psp_h);         /**< Left top. */
Vec2d psp_rt(psp_x + psp_l, psp_y + psp_h); /**< Right top. */
Vec2d psp_rb(psp_x + psp_l, psp_y);         /**< Right bottom. */
/** Topplate. */
Real tp_x = bp_x + 0.18;
Real tp_y = bp_y + 0.36;
Real tp_l = 0.94;
Real tp_h = 0.11;

Vec2d tp_lb(tp_x, tp_y);               /**< Left bottom. */
Vec2d tp_lt(tp_x, tp_y + tp_h);        /**< Left top. */
Vec2d tp_rt(tp_x + tp_l, tp_y + tp_h); /**< Right top. */
Vec2d tp_rb(tp_x + tp_l, tp_y);        /**< Right bottom. */

/**
 * Topology of the tethers.
 * */

Real cxA = bp_x + 0.35;           /**< Center of buoy in x direction. */
Real cyA = bp_y;                  /**< Center of buoy in y direction. */
Vec2d tethering_pointA(cxA, 0.0); /**< The ground tethering point. */
Vec2d cable_endA(cxA, cyA);       /**< The structure tethering point. */

Real cxB = bp_x + bp_l - 0.35;    /**< Center of buoy in x direction. */
Real cyB = bp_y;                  /**< Center of buoy in y direction. */
Vec2d tethering_pointB(cxB, 0.0); /**< The ground tethering point. */
Vec2d cable_endB(cxB, cyB);       /**< The structure tethering point. */

Real cablength = bp_y;

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
Real rho0_f = 1000.0;                    /**< Reference density of fluid. */
Real gravity_g = 9.81;                   /**< Value of gravity. */
Real U_f = 2.0 * sqrt(0.79 * gravity_g); /**< Characteristic velocity. */
Real c_f = 10.0 * U_f;                   /**< Reference sound speed. */
Real mu_f = 1.0e-3;

//----------------------------------------------------------------------
//	Structure Properties G and Inertia
//----------------------------------------------------------------------
/* Weight of the solid structure*/
Real StructureMass = 62.036;

Real bp_area = bp_l * bp_h;
Real ssp_area = ssp_l * ssp_h;
Real psp_area = psp_l * psp_h;
Real tp_area = tp_l * tp_h;

/* Surface of the solid structure*/
Real Area = bp_area + ssp_area + psp_area + tp_area;
/**< Density of the solid structure*/
Real rho_s = StructureMass / Area;

Real rho_bp = rho_s;
Vec2d bp_cm(bp_x + bp_l / 2, bp_y + bp_h / 2);
Vec3d Ibp(
    (rho_bp * bp_area) / 12 * (bp_l * bp_h * bp_h * bp_h),
    (rho_bp * bp_area) / 12 * (bp_h * bp_l * bp_l * bp_l),
    (rho_bp * bp_area) / 12 * (bp_l * bp_l + bp_h * bp_h));

Real rho_ssp = rho_s;
Vec2d ssp_cm(ssp_x + ssp_l / 2, ssp_y + ssp_h / 2);
Vec3d Issp(
    (rho_s * ssp_area) / 12 * (ssp_l * ssp_h * ssp_h * ssp_h),
    (rho_s * ssp_area) / 12 * (ssp_h * ssp_l * ssp_l * ssp_l),
    (rho_s * ssp_area) / 12 * (ssp_l * ssp_l + ssp_h * ssp_h));

Real rho_psp = rho_s;
Vec2d psp_cm(psp_x + psp_l / 2, psp_y + psp_h / 2);
Vec3d Ipsp(
    (rho_s * psp_area) / 12 * (psp_l * psp_h * psp_h * psp_h),
    (rho_s * psp_area) / 12 * (psp_h * psp_l * psp_l * psp_l),
    (rho_s * psp_area) / 12 * (psp_l * psp_l + psp_h * psp_h));

Real rho_tp = rho_s;
Vec2d tp_cm(tp_x + tp_l / 2, tp_y + tp_h / 2);
Vec3d Itp(
    (rho_s * tp_area) / 12 * (tp_l * tp_h * tp_h * tp_h),
    (rho_s * tp_area) / 12 * (tp_h * tp_l * tp_l * tp_l),
    (rho_s * tp_area) / 12 * (tp_l * tp_l + tp_h * tp_h));

Real bcmx = (rho_bp * bp_area * bp_cm[0] +
             rho_ssp * ssp_area * ssp_cm[0] +
             rho_psp * psp_area * psp_cm[0] +
             rho_tp * tp_area * tp_cm[0]) /
            (rho_bp * bp_area +
             rho_ssp * ssp_area +
             rho_psp * psp_area +
             rho_tp * tp_area);
Real bcmy = (rho_bp * bp_area * bp_cm[1] +
             rho_ssp * ssp_area * ssp_cm[1] +
             rho_psp * psp_area * psp_cm[1] +
             rho_tp * tp_area * tp_cm[1]) /
            (rho_bp * bp_area +
             rho_ssp * ssp_area +
             rho_psp * psp_area +
             rho_tp * tp_area);

Vec2d G(bcmx, bcmy);
Real d_bp = sqrt((G[0] - bp_cm[0]) * (G[0] - bp_cm[0]) + (G[1] - bp_cm[1]) * (G[1] - bp_cm[1]));
Real d_ssp = sqrt((G[0] - ssp_cm[0]) * (G[0] - ssp_cm[0]) + (G[1] - ssp_cm[1]) * (G[1] - ssp_cm[1]));
Real d_psp = sqrt((G[0] - psp_cm[0]) * (G[0] - psp_cm[0]) + (G[1] - psp_cm[1]) * (G[1] - psp_cm[1]));
Real d_tp = sqrt((G[0] - tp_cm[0]) * (G[0] - tp_cm[0]) + (G[1] - tp_cm[1]) * (G[1] - tp_cm[1]));

Real Ix = (Ibp[0] + rho_bp * bp_area * (d_bp * d_bp) +
           Issp[0] + rho_ssp * ssp_area * (d_ssp * d_ssp) +
           Ipsp[0] + rho_psp * psp_area * (d_psp * d_psp) +
           Itp[0] + rho_tp * tp_area * (d_tp * d_tp));
Real Iy = (Ibp[1] + rho_bp * bp_area * (d_bp * d_bp) +
           Issp[1] + rho_ssp * ssp_area * (d_ssp * d_ssp) +
           Ipsp[1] + rho_psp * psp_area * (d_psp * d_psp) +
           Itp[1] + rho_tp * tp_area * (d_tp * d_tp));
Real Iz = (Ibp[2] + rho_bp * bp_area * (d_bp * d_bp) +
           Issp[2] + rho_ssp * ssp_area * (d_ssp * d_ssp) +
           Ipsp[2] + rho_psp * psp_area * (d_psp * d_psp) +
           Itp[2] + rho_tp * tp_area * (d_tp * d_tp));

/**
 * Structure observer position
 * */

Vec2d obs = bp_cm;

//------------------------------------------------------------------------------
// geometric shape elements used in the case
//------------------------------------------------------------------------------
MultiPolygon createStructureShape()
{
    /** Geometry definition. */
    std::vector<Vecd> sructure_bp;
    sructure_bp.push_back(bp_lb);
    sructure_bp.push_back(bp_lt);
    sructure_bp.push_back(bp_rt);
    sructure_bp.push_back(bp_rb);
    sructure_bp.push_back(bp_lb);

    std::vector<Vecd> sructure_ssp;
    sructure_ssp.push_back(ssp_lb);
    sructure_ssp.push_back(ssp_lt);
    sructure_ssp.push_back(ssp_rt);
    sructure_ssp.push_back(ssp_rb);
    sructure_ssp.push_back(ssp_lb);

    std::vector<Vecd> sructure_psp;
    sructure_psp.push_back(psp_lb);
    sructure_psp.push_back(psp_lt);
    sructure_psp.push_back(psp_rt);
    sructure_psp.push_back(psp_rb);
    sructure_psp.push_back(psp_lb);

    std::vector<Vecd> sructure_tp;
    sructure_tp.push_back(tp_lb);
    sructure_tp.push_back(tp_lt);
    sructure_tp.push_back(tp_rt);
    sructure_tp.push_back(tp_rb);
    sructure_tp.push_back(tp_lb);

    MultiPolygon multi_polygon_;

    multi_polygon_.addAPolygon(sructure_bp, ShapeBooleanOps::add);
    multi_polygon_.addAPolygon(sructure_ssp, ShapeBooleanOps::add);
    multi_polygon_.addAPolygon(sructure_psp, ShapeBooleanOps::add);
    multi_polygon_.addAPolygon(sructure_tp, ShapeBooleanOps::add);

    return multi_polygon_;
}

class FloatingStructure : public MultiPolygonShape
{
  public:
    explicit FloatingStructure(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> sructure_bp;
        sructure_bp.push_back(bp_lb);
        sructure_bp.push_back(bp_lt);
        sructure_bp.push_back(bp_rt);
        sructure_bp.push_back(bp_rb);
        sructure_bp.push_back(bp_lb);

        std::vector<Vecd> sructure_ssp;
        sructure_ssp.push_back(ssp_lb);
        sructure_ssp.push_back(ssp_lt);
        sructure_ssp.push_back(ssp_rt);
        sructure_ssp.push_back(ssp_rb);
        sructure_ssp.push_back(ssp_lb);

        std::vector<Vecd> sructure_psp;
        sructure_psp.push_back(psp_lb);
        sructure_psp.push_back(psp_lt);
        sructure_psp.push_back(psp_rt);
        sructure_psp.push_back(psp_rb);
        sructure_psp.push_back(psp_lb);

        std::vector<Vecd> sructure_tp;
        sructure_tp.push_back(tp_lb);
        sructure_tp.push_back(tp_lt);
        sructure_tp.push_back(tp_rt);
        sructure_tp.push_back(tp_rb);
        sructure_tp.push_back(tp_lb);

        multi_polygon_.addAPolygon(sructure_bp, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(sructure_ssp, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(sructure_psp, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(sructure_tp, ShapeBooleanOps::add);
    }
};

class StructureSystemForSimbody : public SolidBodyPartForSimbody
{
  public:
    StructureSystemForSimbody(SPHBody &sph_body, SharedPtr<Shape> shape_ptr)
        : SolidBodyPartForSimbody(sph_body, shape_ptr)
    {
        body_part_mass_properties_ =
            mass_properties_ptr_keeper_
                .createPtr<SimTK::MassProperties>(StructureMass, SimTK::Vec3(0.0), SimTK::UnitInertia(0, 0, Iz));
    }
};

//----------------------------------------------------------------------
//	Water block
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> water_block_shape;
        water_block_shape.push_back(watS_lb);
        water_block_shape.push_back(watS_lt);
        water_block_shape.push_back(bwh);
        water_block_shape.push_back(b1a);
        water_block_shape.push_back(watS_lb);

        /*Structure substract*/

        std::vector<Vecd> sructure_bp;
        sructure_bp.push_back(bp_lb);
        sructure_bp.push_back(bp_lt);
        sructure_bp.push_back(bp_rt);
        sructure_bp.push_back(bp_rb);
        sructure_bp.push_back(bp_lb);

        std::vector<Vecd> sructure_ssp;
        sructure_ssp.push_back(ssp_lb);
        sructure_ssp.push_back(ssp_lt);
        sructure_ssp.push_back(ssp_rt);
        sructure_ssp.push_back(ssp_rb);
        sructure_ssp.push_back(ssp_lb);

        std::vector<Vecd> sructure_psp;
        sructure_psp.push_back(psp_lb);
        sructure_psp.push_back(psp_lt);
        sructure_psp.push_back(psp_rt);
        sructure_psp.push_back(psp_rb);
        sructure_psp.push_back(psp_lb);

        std::vector<Vecd> sructure_tp;
        sructure_tp.push_back(tp_lb);
        sructure_tp.push_back(tp_lt);
        sructure_tp.push_back(tp_rt);
        sructure_tp.push_back(tp_rb);
        sructure_tp.push_back(tp_lb);

        std::vector<Vecd> sructure_mdp;
        sructure_mdp.push_back(ssp_rb);
        sructure_mdp.push_back(ssp_rt);
        sructure_mdp.push_back(psp_lt);
        sructure_mdp.push_back(psp_lb);
        sructure_mdp.push_back(ssp_rb);

        multi_polygon_.addAPolygon(water_block_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(sructure_bp, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(sructure_ssp, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(sructure_psp, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(sructure_tp, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(sructure_mdp, ShapeBooleanOps::sub);
    }
};

//----------------------------------------------------------------------
//	Wall cases-dependent geometries.
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        /** Geometry definition. */
        std::vector<Vecd> outer_wall_shape;
        outer_wall_shape.push_back(Vecd(-EXS, 0.0) + Vec2d(-BW, -BW));
        outer_wall_shape.push_back(Vec2d(-EXS, DH) + Vec2d(-BW, 0.0));
        outer_wall_shape.push_back(Vec2d(DL, DH) + Vec2d(+BW, 0.0));
        outer_wall_shape.push_back(Vec2d(DL, BEH) + Vec2d(+BW, -BW));
        outer_wall_shape.push_back(b1b + Vec2d(0.0, -BW));
        outer_wall_shape.push_back(b1a + Vec2d(-BW, -BW));
        outer_wall_shape.push_back(Vec2d(-EXS, 0.0) + Vec2d(-BW, -BW));

        std::vector<Vecd> inner_wall_shape;
        inner_wall_shape.push_back(Vecd(-EXS, 0.0));
        inner_wall_shape.push_back(Vec2d(-EXS, DH));
        inner_wall_shape.push_back(Vec2d(DL, DH));
        inner_wall_shape.push_back(Vec2d(DL, BEH));
        inner_wall_shape.push_back(b1b);
        inner_wall_shape.push_back(b1a);
        inner_wall_shape.push_back(Vec2d(-EXS, 0.0));

        std::vector<Vecd> Wmak_shape;
        Wmak_shape.push_back(Wmak_lb);
        Wmak_shape.push_back(Wmak_lt);
        Wmak_shape.push_back(Wmak_rt);
        Wmak_shape.push_back(Wmak_rb);
        Wmak_shape.push_back(Wmak_lb);

        multi_polygon_.addAPolygon(outer_wall_shape, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(inner_wall_shape, ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(Wmak_shape, ShapeBooleanOps::add);
    }
};

//----------------------------------------------------------------------
//	create a wavemaker shape
//----------------------------------------------------------------------
MultiPolygon createWaveMakerShape()
{
    std::vector<Vecd> Wmak_shape;
    Wmak_shape.push_back(Wmak_lb);
    Wmak_shape.push_back(Wmak_lt);
    Wmak_shape.push_back(Wmak_rt);
    Wmak_shape.push_back(Wmak_rb);
    Wmak_shape.push_back(Wmak_lb);

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(Wmak_shape, ShapeBooleanOps::add);
    return multi_polygon;
}

//----------------------------------------------------------------------
//	Boundary condition for wavemaker
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

    Vecd getDisplacement(const Real &time)
    {
        Real dp = 0;
        for (int j = 0; j < (N); j++)
        {
            dp = dp + 0.5 * S[j] * cos(-k[j] * xf - om[j] * (time - tf));
        };

        Vecd displacement{Vecd::Zero()};
        displacement[0] = dp;
        return displacement;
    }

    Vec2d getVelocity(const Real &time)
    {
        Real vl = 0;
        for (int j = 0; j < (N); j++)
        {
            vl = vl + 0.5 * om[j] * S[j] * sin(-k[j] * xf - om[j] * (time - tf));
        };
        Vec2d velocity{Vecd::Zero()};
        velocity[0] = vl;
        return velocity;
    }

    Vec2d getAcceleration(const Real &time)
    {
        Real ax = 0;
        for (int j = 0; j < (N); j++)
        {
            ax = ax - 0.5 * om[j] * om[j] * S[j] * cos(-k[j] * xf - om[j] * (time - tf));
        };
        Vec2d acceleration{Vecd::Zero()};
        acceleration[0] = ax;
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
          h(WH), tf(20.480), xf(12.0), fmn(0.32), fmx(0.96), a(0.0068), N(32), g(gravity_g),
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
MultiPolygon createWaveGauge()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(10.848 - h, 0.0));
    pnts.push_back(Vecd(10.848 - h, 1.4));
    pnts.push_back(Vecd(10.848 + h, 1.4));
    pnts.push_back(Vecd(10.848 + h, 0.0));
    pnts.push_back(Vecd(10.848 - h, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}