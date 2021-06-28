/**
 * @file 	1d_shock_tube.h
 * @brief 	Numerical parameters and body definition for 1D shock_tube.
 * @author 	Zhentong Wang and Xiangyu Hu
 */
 /**
  * @brief 	SPHinXsys Library.
  */
#include "sphinxsys.h"
  /**
 * @brief Namespace cite here.
 */
using namespace SPH;
/**
 * @brief Basic geometry parameters and numerical setup.
 */
Real DL = 5.0; 								/**< Tube length. */
Real particle_spacing_ref = 1.0/200.0; 		/**< Initial reference particle spacing. */
Real DH = particle_spacing_ref * 4; 		/**< Tube height. */
/** Domain bounds of the system. */
BoundingBox system_domain_bounds(Vec2d(-2.0 / 5.0*DL, 0.0), Vec2d(3.0 / 5.0*DL, DH));
/**
 * @brief Material properties of the fluid.
 */
Real rho0_l = 1.0;						    /**< initial density of left state. */
Real rho0_r = 0.125;						/**< initial density of right state. */
Vecd velocity_l(0.0);						/**< initial velocity of left state. */
Vecd velocity_r(0.0);						/**< initial velocity of right state. */
Real p_l = 1.0;						        /**< initial pressure of left state. */
Real p_r = 0.1;						        /**< initial pressure of right state. */
Real heat_capacity_ratio = 1.4;             /**< heat capacity ratio. */
/** create the left wave block shape */
std::vector<Vecd> CreatWavesBlockShape()
{
	//geometry
	std::vector<Vecd> waves_block_shape;
	waves_block_shape.push_back(Vecd(-2.0 / 5.0*DL, 0.0));
	waves_block_shape.push_back(Vecd(-2.0 / 5.0*DL, DH));
	waves_block_shape.push_back(Vecd(3.0 / 5.0*DL, DH));
	waves_block_shape.push_back(Vecd(3.0 / 5.0*DL, 0.0));
	waves_block_shape.push_back(Vecd(-2.0 / 5.0*DL, 0.0));
	return waves_block_shape;
}
/**
*@brief 	Left Wave body definition.
*/
class WaveBlock : public EulerianFluidBody
{
public:
	WaveBlock(SPHSystem& sph_system, std::string body_name)
		: EulerianFluidBody(sph_system, body_name)
	{
		/** Geomtry definition. */
		std::vector<Vecd> waves_block_shape = CreatWavesBlockShape();
		body_shape_ = new ComplexShape(body_name);
		body_shape_->addAPolygon(waves_block_shape, ShapeBooleanOps::add);
	}
};
/**
 * @brief 	Case dependent material properties definition.
 */
class WaveMaterial : public CompressibleFluid
{
public:
	WaveMaterial() : CompressibleFluid()
	{
		/** Basic material parameter*/
		gamma_ = heat_capacity_ratio;

		assignDerivedMaterialParameters();
	}
};
/**
 * @brief 	Case dependent material properties definition.
 */
/**
 * application dependent initial condition
 */
class WavesInitialCondition
	: public eulerian_fluid_dynamics::CompressibleFluidInitialCondition
{
public:
	WavesInitialCondition(EulerianFluidBody* water)
		: eulerian_fluid_dynamics::CompressibleFluidInitialCondition(water) {};
protected:
	void Update(size_t index_i, Real dt) override
	{
		if (pos_n_[index_i][0] < DL / 10.0)
		{
			/** initial left wave pressure,momentum and energy profile */
			rho_n_[index_i] = rho0_l;
			p_[index_i] = p_l;
			Real rho_e = p_[index_i] / (gamma_ - 1.0);
			vel_n_[index_i] = velocity_l;
			mom_[index_i] = rho0_l * velocity_l;
			E_[index_i] = rho_e + 0.5 * rho_n_[index_i] * vel_n_[index_i].normSqr();
		}
		if (pos_n_[index_i][0] > DL / 10.0)
		{
			/** initial right wave pressure,momentum and energy profile */
			rho_n_[index_i] = rho0_r;
			p_[index_i] = p_r;
			Real rho_e = p_[index_i] / (gamma_ - 1.0);
			vel_n_[index_i] = velocity_r;
			mom_[index_i] = rho0_r * velocity_r;
			E_[index_i] = rho_e + 0.5 * rho_n_[index_i] * vel_n_[index_i].normSqr();
		}
	}
};
