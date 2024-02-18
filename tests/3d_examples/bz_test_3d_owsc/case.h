/* -------------------------------------------------------------------------*
*								SPHinXsys									*
* --------------------------------------------------------------------------*
* SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle	*
* Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
* physical accurate simulation and aims to model coupled industrial dynamic *
* systems including fluid, solid, multi-body dynamics and beyond with SPH	*
* (smoothed particle hydrodynamics), a meshless computational method using	*
* particle discretization.													*
*																			*
* SPHinXsys is partially funded by German Research Foundation				*
* (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1				*
* and HU1527/12-1.															*
*                                                                           *
* Portions copyright (c) 2017-2020 Technical University of Munich and		*
* the authors' affiliations.												*
*                                                                           *
* Licensed under the Apache License, Version 2.0 (the "License"); you may   *
* not use this file except in compliance with the License. You may obtain a *
* copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
*                                                                           *
* --------------------------------------------------------------------------*/
/**
* @file 	case.h
* @brief 	This is the case file for the test of Oscillating Wave Surge Converter (OWSC) in 3D wave tank.
* @author   Chi Zhang and Xiangyu Hu
* @version  0.3.0
* @note  	Observer, moving with mobile solid body, can find template in this case.
*			-- Chi ZHANG
*/
#pragma once

#include "sphinxsys.h"
using namespace SPH;

#define PI 3.1415926
Real total_physical_time = 12.0;
/** Set the file path to the stl file. */
Real length_scale = 0.001;
Real Water_H = 0.691; /**< Water height. */
Real Flap_width = 0.12;
Real dp_0 		= Flap_width / 4.0;
Real BW = 4.0 * dp_0;
BoundingBox system_domain_bounds(Vec3d(-1.0 - BW, 0.0 - BW, -2.29 - BW), Vec3d(18.42 + BW, 1.0 + BW, 2.29 + BW));

std::string path_to_wall_stl = "./input/tank.stl";
std::string path_to_wave_maker_stl = "./input/wave_maker.stl";
std::string path_to_water_stl = "./input/water.stl";
std::string path_to_flap_stl = "./input/flap.stl";

Vecd halfsize_damping(2.5, 0.175, 2.29);
Vecd translation_damping(15.92, 0.525, 0.0);

//WaveGauges 
Real h = 1.3 * dp_0;
Vecd halfsize_No4(h, 0.5, h);
Vecd translation_No4(3.99, 0.5, 0.0);
Vecd halfsize_No5(h, 0.5, h);
Vecd translation_No5(7.02, 0.5, 0.0);
Vecd halfsize_No12(h, 0.5, h);
Vecd translation_No12(8.82, 0.5, 0.0);

Vecd translation_wave_maker(0, 0, 0);
Vecd translation_flap(0, 0.5 * dp_0, 0);

/** Material parameters. */
Real gravity_g = 9.81;
/** material properties of the fluid.*/
Real rho0_f = 1000.0;
Real U_f    = 2.0 * sqrt(Water_H * gravity_g);
Real c_f    = 10.0 * U_f;
Real mu_f   = 1.0e-6;
/** properties of the OWSC*/
Real flap_mass = 33.04;
Real flap_vol  = 0.0599;
Real rho0_s = flap_mass / flap_vol;
Real poisson = 0.33;
Real Youngs_modulus = 7.8e6;

// the offset that the rubber flap shifted above the tank
// Real flap_off = Flap_x - 0.5 * Flap_width + DL_Extra + BW;
// Real off_set = particle_spacing_ref + floor(flap_off / particle_spacing_ref) * particle_spacing_ref - flap_off;
Vec3d offset = Vec3d::Zero();;

/** Define the wave tank geometry. */
TriangleMeshShape* CreateCADGeometryForTank()
{
	Vecd translation(0.0, 0.0, 0.0);
	TriangleMeshShape *tank_geometry = new TriangleMeshShapeSTL(path_to_wall_stl, translation, length_scale);

	return tank_geometry;
}
/** Define the wave maker geometry. */
TriangleMeshShape* CreateCADGeometryForWaveMaker()
{
	Vecd translation(0.0, 0.0, 0.0);
	TriangleMeshShape *tank_geometry = new TriangleMeshShapeSTL(path_to_wave_maker_stl, translation, length_scale);

	return tank_geometry;
}
/** Define the water geometry. */
TriangleMeshShape* CreateCADGeometryForWater()
{
	Vecd translation(0.0, 0.0, 0.0);
	TriangleMeshShape *tank_geometry = new TriangleMeshShapeSTL(path_to_water_stl, translation, length_scale);

	return tank_geometry;
}
/** Define the OWSC geometry. */
TriangleMeshShape* CreateCADGeometryForOWSC()
{
	Vecd translation(0.0, 0.5 * dp_0, 0.0);
	TriangleMeshShape *tank_geometry = new TriangleMeshShapeSTL(path_to_flap_stl, translation, length_scale);

	return tank_geometry;
}
/** Define Tank Body. */
class MyTank : public ComplexShape
{
public:
	explicit MyTank(const std::string &shape_name)
		: ComplexShape(shape_name)
	{
		Vecd translation(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(path_to_wall_stl, translation, length_scale);

		//ComplexShape original_body_shape;
		//original_body_shape.TriangleMeshShape(CreateCADGeometryForTank(), ShapeBooleanOps::add);
		//body_shape_ = new LevelSetComplexShape(this, original_body_shape);		
	}
};
/** Define the water body. */
class Water : public ComplexShape
{
public:
	explicit Water(const std::string &shape_name)
		: ComplexShape(shape_name)
	{
		Vecd translation(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(path_to_water_stl, translation, length_scale);
		//body_shape_ = new ComplexShape(body_name);
		//body_shape_->addTriangleMeshShape(CreateCADGeometryForWater(), ShapeBooleanOps::add);
	}
};
/** Define Tank Body. */
class MyWaveMaker : public ComplexShape
{
public:
	explicit MyWaveMaker(const std::string &shape_name)
		: ComplexShape(shape_name)
	{
		Vecd translation(0.0, 0.0, 0.0);
		add<TriangleMeshShapeSTL>(path_to_wave_maker_stl, translation, length_scale);
		//body_shape_ = new ComplexShape(body_name);
		//body_shape_->addTriangleMeshShape(CreateCADGeometryForWaveMaker(), ShapeBooleanOps::add);
	}
};

/** Define Tank Body. */
class MyFlap : public ComplexShape
{
public:
	explicit MyFlap(const std::string &shape_name)
		: ComplexShape(shape_name)
	{
		Vecd translation(0.0, 0.5 * dp_0, 0.0);
		add<TriangleMeshShapeSTL>(path_to_flap_stl, translation, length_scale);
		//ComplexShape original_body_shape;
		//original_body_shape.addTriangleMeshShape(CreateCADGeometryForOWSC(), ShapeBooleanOps::add);
		//body_shape_ = new LevelSetComplexShape(this, original_body_shape, true);
	}
};
/**
*/

class FlapSystemForSimbody : public SolidBodyPartForSimbody
{
public:
	FlapSystemForSimbody(SPHBody &sph_body, SharedPtr<Shape> shape_ptr)
		: SolidBodyPartForSimbody(sph_body, shape_ptr)
	{
		// Vecd mass_center = Vecd(7.92, 0.355); // 0.3355
		// initial_mass_center_ = SimTK::Vec3(mass_center[0], mass_center[1], 0.0);
		/** UnitInertia_ (const RealP &xx, const RealP &yy, const RealP &zz)
		 * 	Create a principal unit inertia matrix (only non-zero on diagonal).
		 */
		Real Iz = 1.84 / 33.04;
		body_part_mass_properties_ =
			mass_properties_ptr_keeper_
			.createPtr<SimTK::MassProperties>(33.04, SimTK::Vec3(0.0), SimTK::UnitInertia(0.0, 0.0, Iz));
	}
};
/**
* @brief making the wave
*/
class WaveMaking : public solid_dynamics::BaseMotionConstraint<SolidBody>
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
		Vecd displacement{ Vecd::Zero() };
		displacement[0] = 0.5 * wave_stroke_ * sin(wave_freq_ * time);
		return displacement;
	}

	Vec3d getVelocity(const Real &time)
	{
		Vec3d velocity{ Vec3d::Zero() };
		velocity[0] = 0.5 * wave_stroke_ * wave_freq_ * cos(wave_freq_ * time);
		return velocity;
	}

	Vec3d getAcceleration(const Real &time)
	{
		Vec3d acceleration{ Vec3d::Zero() };
		acceleration[0] = - 0.5 * wave_stroke_ * wave_freq_ * wave_freq_ * sin(wave_freq_ * time);
		return acceleration;
	}
	void computeWaveStrokeAndFrequency( )
	{
		Real scaled_wave_height = wave_height_ / model_scale_;
		Real scaled_wave_period = wave_period_ / sqrt(model_scale_);
		Real scaled_wave_freq = 2.0 * PI / scaled_wave_period;
		Real scaled_wave_amp  = 0.5 * scaled_wave_height;

		int iterator = 20; 
		Real Tol = 1.0e-6;
		Real scaled_wave_number = 1.0;
		for(int i = 1; i < iterator; i++)
		{
			Real term1 = tanh(scaled_wave_number * water_depth_);
			Real term2 = scaled_wave_freq * scaled_wave_freq / gravity_;
			Real term3 = scaled_wave_number * term1 - term2;
			Real term4 = term1 + scaled_wave_number * water_depth_ * (1.0 - term1 * term1);
			Real wave_number_old = scaled_wave_number;
			scaled_wave_number = wave_number_old - term3 / term4;
			Real error = abs(scaled_wave_number - wave_number_old) / abs(scaled_wave_number);
			if(error <= Tol) 
				break;
		}

		Real term_1 = gravity_ / scaled_wave_freq / scaled_wave_freq;
		Real term_2 = 2.0 * scaled_wave_number * water_depth_;
		Real term_3 = scaled_wave_number * water_depth_;
		Real scaled_wave_stroke = 0.5 * scaled_wave_amp * scaled_wave_number * term_1*                  
			(term_2 + sinh(term_2))/ (cosh(term_3) * sinh(term_3));

		wave_stroke_ = scaled_wave_stroke;
		wave_freq_ = scaled_wave_freq;
		std::cout<< "Wave number: " << scaled_wave_number << " Wave stroke: " << wave_stroke_ << " Wave frequency: " << wave_freq_ << std::endl;
	}
public:
	WaveMaking(SolidBody &solid_body)
		: solid_dynamics::BaseMotionConstraint<SolidBody>(solid_body)
	{
		model_scale_ = 25.0;
		wave_height_ = 5.0;
		wave_period_ = 10.0;
		gravity_ 	 = gravity_g;
		water_depth_ = Water_H;
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

/** Flap observer body */
class FlapObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	explicit FlapObserverParticleGenerator(SPHBody &sph_body) : ObserverParticleGenerator(sph_body)
	{
		/** the measuring particle with zero volume */
		positions_.push_back(Vecd(7.872, 0.645, 0.468));
		positions_.push_back(Vecd(7.872, 0.741, 0.364));
		positions_.push_back(Vecd(7.872, 0.391, 0.364));
		positions_.push_back(Vecd(7.872, 0.574, 0.156));
		positions_.push_back(Vecd(7.872, 0.716, 0.052));
		positions_.push_back(Vecd(7.872, 0.452, 0.052));
	}
};