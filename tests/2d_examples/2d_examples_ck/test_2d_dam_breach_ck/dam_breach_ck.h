#ifndef DAM_BREACH_H
#define DAM_BREACH_H

#include "sphinxsys_ck.h"
using namespace SPH; // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.0;					/**< Water tank length. */
Real DH = 1.0;					/**< Water tank height. */
Real WL = 1.4;					/**< Water column length. */
Real WH = 0.2;					/**< Water column height. */
Real SH = 0.2;					/**< Sediment column height. */
Real particle_spacing_ref = 0.01;	/**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;  /**< Extending width for boundary conditions. */
Real DL_sponge = particle_spacing_ref * 20.0; /**< Sponge region to impose inflow condition. */
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
// observer location
StdVec<Vecd> observation_location = {Vecd(DL, 0.2)};
//----------------------------------------------------------------------
//	Material properties of the soil.
//----------------------------------------------------------------------
/*Soil parameters*/
Real rho0_s = 2040;                                                       // reference density of soil
Real gravity_g = 9.8;                                                     // gravity force of soil
Real Youngs_modulus = 5.84e6;                                             // reference Youngs modulus
Real poisson = 0.3;                                                       // Poisson ratio
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3.0 * (1.0 - 2.0 * poisson))); // sound speed
Real friction_angle = 30.0 * Pi / 180;
Real cohesion = 400.0;
Real U_max_s = 2.0 * sqrt(gravity_g * SH);
/*Fluid parameters*/
Real rho0_f = 1000;
Real U_max_f = 2.0 * sqrt(gravity_g * WH);
Real c_f = 10.0 * U_max_f;
Real mu_f = 1e-6;

Real U_ref = SMAX(U_max_s, U_max_f);

//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
//dike shape
Vec2d offset(0.0,0.0);
Vec2d dike_lb = Vec2d(1.66, 0.5) + offset;	  // left bottom
Vec2d dike_lt = Vec2d(2.06, 0.7) + offset; // left top
Vec2d dike_rt = Vec2d(2.16, 0.7) + offset;  // right top
Vec2d dike_rb = Vec2d(2.56, 0.5) + offset;	  // right bottom

//inflow buffer
Vec2d emitter_halfsize = Vec2d(0.5 * 0.66, 1.0 * BW);
Vec2d emitter_translation = Vec2d(0.5 * 0.66, 1.0 * BW);
Vec2d emitter_buffer_halfsize = Vec2d(0.5 * 0.66, 2.0 * DL_sponge);
Vec2d emitter_buffer_translation = Vec2d(0.5 * 0.66, -0.5 * DL_sponge);

//outflow buffer
Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(4.0, 0.5 + DH + 0.25 * DH) - disposer_halfsize;
//Entire buffer
Vec2d entire_disposer_halfsize = Vec2d(DL, DH);
Vec2d entire_disposer_translation = disposer_halfsize;
//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
/*Wall Boundary*/
std::vector<Vecd> createOuterWallShape()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(-BW, -BW));
	pnts.push_back(Vecd(-BW, 1.0 + BW));
	pnts.push_back(Vecd(4.2 + BW, 1.0 + BW));
	pnts.push_back(Vecd(4.2 + BW, 0.5 - BW));
	pnts.push_back(Vecd(0.66 + BW, 0.5 - BW));
	pnts.push_back(Vecd(0.66 + BW, -BW));
	pnts.push_back(Vecd(-BW, -BW));

	return pnts;
}
std::vector<Vecd> createInnerWallShape()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(0.0, -BW));
	pnts.push_back(Vecd(0.0, 1.0));
	pnts.push_back(Vecd(4.2, 1.0));
	pnts.push_back(Vecd(4.2, 0.5));
	pnts.push_back(Vecd(0.66, 0.5));
	pnts.push_back(Vecd(0.66, -BW));
	pnts.push_back(Vecd(0.0, -BW));

	return pnts;
}

/*Dam shape*/
std::vector<Vecd> createDikeBlockShape()
{
	std::vector<Vecd> pnts;
	pnts.push_back(dike_lb);
	pnts.push_back(dike_lt);
	pnts.push_back(dike_rt);
	pnts.push_back(dike_rb);
	pnts.push_back(dike_lb);

	return pnts;
}
/*Water shape*/
std::vector<Vecd> createWaterBlockShape()
{
	std::vector<Vecd> pnts;
	pnts.push_back(Vecd(0.0, 0.0));
	pnts.push_back(Vecd(0.0, 0.6));
	pnts.push_back(Vecd(2.06, 0.6));
	pnts.push_back(Vecd(2.06, 0.5));
	pnts.push_back(Vecd(0.66, 0.5));
	pnts.push_back(Vecd(0.66, 0.0));
	pnts.push_back(Vecd(0.0, 0.0));

	return pnts;
}

//----------------------------------------------------------------------
//	Complex for different sections
//----------------------------------------------------------------------
class WallBoundary : public MultiPolygonShape
{
public:
	explicit WallBoundary(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createOuterWallShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createInnerWallShape(), ShapeBooleanOps::sub);
	}
};

class WaterBlock : public MultiPolygonShape
{
public:
	explicit WaterBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(createDikeBlockShape(), ShapeBooleanOps::sub);
	}
};

class SoilBlock : public MultiPolygonShape
{
public:
	explicit SoilBlock(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createDikeBlockShape(), ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	Inlet inflow condition
//----------------------------------------------------------------------
class InletInflowCondition : public BaseStateCondition
{
  public:
    InletInflowCondition(BaseParticles *particles)
        : BaseStateCondition(particles) {};

    class ComputingKernel : public BaseStateCondition::ComputingKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        ComputingKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : BaseStateCondition::ComputingKernel(ex_policy, encloser){};

        void operator()(AlignedBox *aligned_box, UnsignedInt index_i, Real /*time*/)
        {
            vel_[index_i] = Vec2d(0.0, 0.048);
        };
    };
};
//----------------------------------------------------------------------
//	Outflow  condition
//----------------------------------------------------------------------

#endif //DAM_BREACH_H