/**
 * @file 	owsc.h
 * @brief 	This is the case file for the test of Oscillating Wave Surge Converter (OWSC).
 * @author   Chi Zhang and Xiangyu Hu
 */
#ifndef TEST_2D_OWSC_CASE_H
#define TEST_2D_OWSC_CASE_H

#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real total_physical_time = 12.0;

Real DL = 18.42;      // tank length
Real DH = 1.0;        // tank height
Real DL_Extra = 1.0;  // for wave maker
Real Water_H = 0.691; /**< Water height. */

Real Flap_width = 0.12;
Real Flap_x = 7.92;
Real Flap_H = 0.48;

Real Base_bottom_position = 0.155;
Real Base_height = 0.1;
Real particle_spacing_ref = Flap_width / 4.0; // particle spacing
Real BW = particle_spacing_ref * 4.0;         // boundary width

BoundingBox system_domain_bounds(Vec2d(-DL_Extra - BW, -BW), Vec2d(DL + BW, DH + BW));

// the offset that the rubber flap shifted above the tank
// Real flap_off = Flap_x - 0.5 * Flap_width + DL_Extra + BW;
// Real off_set = particle_spacing_ref + floor(flap_off / particle_spacing_ref) * particle_spacing_ref - flap_off;
Vec2d offset = Vec2d::Zero();

// water block parameters
Vec2d Water_lb(0.0, 0.0);     // left bottom
Vec2d Water_lt(0.0, Water_H); // left top
Vec2d Water_rt(DL, Water_H);  // right top
Vec2d Water_rb(DL, 0.356);    // right bottom
Vec2d Water_slope_1(DL - 6.2, 0.356);
Vec2d Water_slope_2(DL - 6.2 - 3.7, 0.155);
Vec2d Water_slope_3(DL - 6.2 - 3.7 - 2.4, 0.155);
Vec2d Water_slope_4(DL - 6.2 - 3.7 - 2.4 - 1.3, 0.0);

// flap constrain region parameter
Vec2d Base_lb(Flap_x - 0.5 * Flap_width, Base_bottom_position);               // left bottom
Vec2d Base_lt(Flap_x - 0.5 * Flap_width, Base_bottom_position + Base_height); // left top
Vec2d Base_rt(Flap_x + 0.5 * Flap_width, Base_bottom_position + Base_height); // right top
Vec2d Base_rb(Flap_x + 0.5 * Flap_width, Base_bottom_position);               // right bottom

// flap geometric parameters
Vec2d Flap_lb(Flap_x - 0.5 * Flap_width, Base_bottom_position + Base_height + 0.5 * Flap_width);          // left bottom
Vec2d Flap_lt(Flap_x - 0.5 * Flap_width, Base_bottom_position + Base_height + 0.5 * Flap_width + Flap_H); // left top
Vec2d Flap_rt(Flap_x + 0.5 * Flap_width, Base_bottom_position + Base_height + 0.5 * Flap_width + Flap_H); // right top
Vec2d Flap_rb(Flap_x + 0.5 * Flap_width, Base_bottom_position + Base_height + 0.5 * Flap_width);          // right bottom

// gravity
Real gravity_g = 9.81;

// for material properties of the fluid
Real rho0_f = 1000.0;
Real U_f = 2.0 * sqrt(0.79 * gravity_g);
Real c_f = 10.0 * U_f;
Real mu_f = 1.0e-6;

// for material properties of the solid
Real flap_mass = 33.04;
Real flap_volume = 0.0579;
Real rho0_s = flap_mass / flap_volume;

//------------------------------------------------------------------------------
// geometric shape elements used in the case
//------------------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Water_lb);
    pnts.push_back(Water_lt);
    pnts.push_back(Water_rt);
    pnts.push_back(Water_rb);
    pnts.push_back(Water_slope_1);
    pnts.push_back(Water_slope_2);
    pnts.push_back(Water_slope_3);
    pnts.push_back(Water_slope_4);
    pnts.push_back(Water_lb);

    return pnts;
}

std::vector<Vecd> createFlapConstrainShape()
{
    std::vector<Vecd> pnts2;
    pnts2.push_back(Base_lb);
    pnts2.push_back(Base_lt);
    pnts2.push_back(Base_rt);
    pnts2.push_back(Base_rb);
    pnts2.push_back(Base_lb);

    return pnts2;
}

std::vector<Vecd> createFlapShape()
{
    std::vector<Vecd> pnts3;
    pnts3.push_back(Flap_lb);
    pnts3.push_back(Flap_lt);
    pnts3.push_back(Flap_rt);
    pnts3.push_back(Flap_rb);
    for (int i = 1; i <= 10; i++)
    {
        Real angle = Real(i) * PI / (11.0);
        Real x = Flap_rb[0] - 0.5 * Flap_width * (1.0 - cos(angle));
        Real y = Flap_rb[1] - 0.5 * Flap_width * sin(angle) - 0.5 * particle_spacing_ref;
        pnts3.push_back(Vecd(x, y));
    }
    pnts3.push_back(Flap_lb);

    return pnts3;
}

MultiPolygon createDampingBufferShape()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(DL - 5.0, 0.356 - BW));
    pnts.push_back(Vecd(DL - 5.0, DH));
    pnts.push_back(Vecd(DL + BW, DH));
    pnts.push_back(Vecd(DL + BW, 0.356 - BW));
    pnts.push_back(Vecd(DL - 5.0, 0.356 - BW));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}

std::vector<Vecd> createOuterWallShape()
{
    std::vector<Vecd> pnts1;
    pnts1.push_back(Vecd(-DL_Extra - BW, -BW));
    pnts1.push_back(Vecd(-DL_Extra - BW, DH + BW));
    pnts1.push_back(Vecd(DL + BW, DH + BW));
    pnts1.push_back(Vecd(DL + BW, 0.35 - BW));
    pnts1.push_back(Water_slope_1 + Vec2d(0.0, -BW));
    pnts1.push_back(Water_slope_2 + Vec2d(0.0, -BW));
    pnts1.push_back(Water_slope_3 + Vec2d(0.0, -BW));
    pnts1.push_back(Water_slope_4 + Vec2d(0.0, -BW));
    pnts1.push_back(Vecd(-DL_Extra - BW, -BW));

    return pnts1;
}

std::vector<Vecd> createInnerWallShape01()
{
    std::vector<Vecd> pnts2;
    pnts2.push_back(Water_lb);
    pnts2.push_back(Vecd(0.0, DH + BW));
    pnts2.push_back(Vecd(DL, DH + BW));
    pnts2.push_back(Water_rb);
    pnts2.push_back(Water_slope_1);
    pnts2.push_back(Water_slope_2);
    pnts2.push_back(Base_rb);
    pnts2.push_back(Base_rt);
    pnts2.push_back(Base_lt);
    pnts2.push_back(Base_lb);
    pnts2.push_back(Water_slope_3);
    pnts2.push_back(Water_slope_4);
    pnts2.push_back(Water_lb);

    return pnts2;
}

std::vector<Vecd> createInnerWallShape02()
{
    std::vector<Vecd> pnts3;
    pnts3.push_back(Vecd(-DL_Extra, 0.0));
    pnts3.push_back(Vecd(-DL_Extra, DH + BW));
    pnts3.push_back(Vecd(-BW, DH + BW));
    pnts3.push_back(Vecd(-BW, 0.0));
    pnts3.push_back(Vecd(-DL_Extra, 0.0));

    return pnts3;
}

MultiPolygon createWaveMakerShape()
{
    std::vector<Vecd> wave_make_shape;
    wave_make_shape.push_back(Vecd(-BW, 0.0));
    wave_make_shape.push_back(Vecd(-BW, DH + BW));
    wave_make_shape.push_back(Vecd(0.0, DH + BW));
    wave_make_shape.push_back(Vecd(0.0, 0.0));
    wave_make_shape.push_back(Vecd(-BW, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(wave_make_shape, ShapeBooleanOps::add);
    return multi_polygon;
}
//------------------------------------------------------------------------------
// geometric shapes for the bodies used in the case
//------------------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createFlapShape(), ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(createFlapConstrainShape(), ShapeBooleanOps::sub);
    }
};

class WallBoundary : public MultiPolygonShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createFlapConstrainShape(), ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(createInnerWallShape01(), ShapeBooleanOps::sub);
        multi_polygon_.addAPolygon(createInnerWallShape02(), ShapeBooleanOps::sub);
    }
};

class Flap : public MultiPolygonShape
{
  public:
    explicit Flap(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createFlapShape(), ShapeBooleanOps::add);
    }
};
//------------------------------------------------------------------------------
// Body parts used in the case
//------------------------------------------------------------------------------
MultiPolygon createFlapSimbodyConstrainShape()
{
    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(createFlapShape(), ShapeBooleanOps::add);
    return multi_polygon;
};

class FlapSystemForSimbody : public SolidBodyPartForSimbody
{
  public:
    FlapSystemForSimbody(SPHBody &sph_body, SharedPtr<Shape> shape_ptr)
        : SolidBodyPartForSimbody(sph_body, shape_ptr)
    {
        // Vecd mass_center = Vecd(7.92, 0.355); // 0.3355
        // initial_mass_center_ = SimTKVec3(mass_center[0], mass_center[1], 0.0);
        /** UnitInertia_ (const RealP &xx, const RealP &yy, const RealP &zz)
         * 	Create a principal unit inertia matrix (only non-zero on diagonal).
         */
        Real Iz = 1.84 / 33.04;
        body_part_mass_properties_ =
            mass_properties_ptr_keeper_
                .createPtr<SimTK::MassProperties>(33.04, SimTKVec3(0.0), SimTK::UnitInertia(0.0, 0.0, Iz));
    }
};

class WaveMaking : public solid_dynamics::BaseMotionConstraint<BodyPartByParticle>
{
    Real model_scale_;
    Real gravity_;
    Real water_depth_;
    Real wave_height_;
    Real wave_period_;
    Real wave_freq_;
    Real wave_stroke_;

    Vecd getDisplacement(const Real &time)
    {
        Vecd displacement{Vecd::Zero()};
        displacement[0] = 0.5 * wave_stroke_ * sin(wave_freq_ * time);
        return displacement;
    }

    Vec2d getVelocity(const Real &time)
    {
        Vec2d velocity{Vecd::Zero()};
        velocity[0] = 0.5 * wave_stroke_ * wave_freq_ * cos(wave_freq_ * time);
        return velocity;
    }

    Vec2d getAcceleration(const Real &time)
    {
        Vec2d acceleration{Vecd::Zero()};
        acceleration[0] = -0.5 * wave_stroke_ * wave_freq_ * wave_freq_ * sin(wave_freq_ * time);
        return acceleration;
    }

    void computeWaveStrokeAndFrequency()
    {
        Real scaled_wave_height = wave_height_ / model_scale_;
        Real scaled_wave_period = wave_period_ / sqrt(model_scale_);
        Real scaled_wave_freq = 2.0 * PI / scaled_wave_period;
        Real scaled_wave_amp = 0.5 * scaled_wave_height;

        int iterator = 20;
        Real Tol = 1.0e-6;
        Real scaled_wave_number = 1.0;
        for (int i = 1; i < iterator; i++)
        {
            Real term1 = tanh(scaled_wave_number * water_depth_);
            Real term2 = scaled_wave_freq * scaled_wave_freq / gravity_;
            Real term3 = scaled_wave_number * term1 - term2;
            Real term4 = term1 + scaled_wave_number * water_depth_ * (1.0 - term1 * term1);
            Real wave_number_old = scaled_wave_number;
            scaled_wave_number = wave_number_old - term3 / term4;
            Real error = abs(scaled_wave_number - wave_number_old) / abs(scaled_wave_number);
            if (error <= Tol)
                break;
        }

        Real term_1 = gravity_ / scaled_wave_freq / scaled_wave_freq;
        Real term_2 = 2.0 * scaled_wave_number * water_depth_;
        Real term_3 = scaled_wave_number * water_depth_;
        Real scaled_wave_stroke = 0.5 * scaled_wave_amp * scaled_wave_number * term_1 *
                                  (term_2 + sinh(term_2)) / (cosh(term_3) * sinh(term_3));

        wave_stroke_ = scaled_wave_stroke;
        wave_freq_ = scaled_wave_freq;
        std::cout << "Wave stroke: " << wave_stroke_ << " Wave frequency: " << wave_freq_ << std::endl;
    }

  public:
    WaveMaking(BodyPartByParticle &body_part)
        : solid_dynamics::BaseMotionConstraint<BodyPartByParticle>(body_part),
          model_scale_(25.0), gravity_(gravity_g), water_depth_(Water_H), wave_height_(5.0),
          wave_period_(10.0)
    {
        computeWaveStrokeAndFrequency();
    }

    void update(size_t index_i, Real dt = 0.0)
    {
        Real time = GlobalStaticVariables::physical_time_;
        pos_[index_i] = pos0_[index_i] + getDisplacement(time);
        vel_[index_i] = getVelocity(time);
        acc_[index_i] = getAcceleration(time);
    };
};

Real h = 1.3 * particle_spacing_ref;
MultiPolygon createWaveProbeShape4()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(3.99 - h, 0.0));
    pnts.push_back(Vecd(3.99 - h, 1.0));
    pnts.push_back(Vecd(3.99 + h, 1.0));
    pnts.push_back(Vecd(3.99 + h, 0.0));
    pnts.push_back(Vecd(3.99 - h, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}
MultiPolygon createWaveProbeShape5()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(7.02 - h, 0.155));
    pnts.push_back(Vecd(7.02 - h, 1.0));
    pnts.push_back(Vecd(7.02 + h, 1.0));
    pnts.push_back(Vecd(7.02 + h, 0.155));
    pnts.push_back(Vecd(7.02 - h, 0.155));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}

MultiPolygon createWaveProbeShape12()
{
    std::vector<Vecd> pnts;
    pnts.push_back(Vecd(8.82 - h, 0.155));
    pnts.push_back(Vecd(8.82 - h, 1.0));
    pnts.push_back(Vecd(8.82 + h, 1.0));
    pnts.push_back(Vecd(8.82 + h, 0.155));
    pnts.push_back(Vecd(8.82 - h, 0.155));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(pnts, ShapeBooleanOps::add);
    return multi_polygon;
}
//------------------------------------------------------------------------------
// Case-dependent observer particle generator
//------------------------------------------------------------------------------
class FlapObserverParticleGenerator : public ObserverParticleGenerator
{
  public:
    explicit FlapObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
    {
        /** the measuring particle with zero volume */
        positions_.push_back(Vecd(7.862, 0.645));
        positions_.push_back(Vecd(7.862, 0.741));
        positions_.push_back(Vecd(7.862, 0.391));
        positions_.push_back(Vecd(7.862, 0.574));
        positions_.push_back(Vecd(7.862, 0.716));
        positions_.push_back(Vecd(7.862, 0.452));
    }
};
#endif // TEST_2D_OWSC_CASE_H
