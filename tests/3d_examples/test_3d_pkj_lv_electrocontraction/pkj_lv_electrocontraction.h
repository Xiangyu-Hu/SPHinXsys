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
 * @file 	pkj_lv_electrocontraction.h
 * @brief 	Electro-contraction of left ventricle heart model.
 * @author 	Chi ZHANG and Xiangyu HU
 * @version  0.3
 * @version 0.2.1
 * 			Chi Zhang
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 *@version 0.3
 *			Here, the coupling with Purkinje network will be conducted.
 */
#pragma once
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
/** Set the file path to the stl file. */
std::string full_path_to_lv = "./input/leftventricle.stl";
Real length_scale = 1.0;
Real time_scale = 1.0 / 12.9;
Real stress_scale = 1.0e-6;
/** Parameters and physical properties. */
Vec3d domain_lower_bound(-90.0 * length_scale, -40.0 * length_scale, -80.0 * length_scale);
Vec3d domain_upper_bound(40.0 * length_scale, 30.0 * length_scale, 50.0 * length_scale);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 100.0;
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);

/** Material properties. */
Real rho0_s = 1.06e-3;
/** Active stress factor */
Real k_a = 150 * stress_scale;
std::array<Real, 4> a0 = {Real(496.0 * stress_scale), Real(15196.0 * stress_scale), Real(3283.0 * stress_scale), Real(662.0 * stress_scale)};
std::array<Real, 4> b0 = {Real(7.209), Real(20.417), Real(11.176), Real(9.466)};
/** reference stress to achieve weakly compressible condition */
Real poisson = 0.4995;
Real bulk_modulus = 2.0 * a0[0] * (1.0 + poisson) / (3.0 * (1.0 - 2.0 * poisson));
/** Electrophysiology parameters. */
Real acceleration_factor = 27.5; /** Acceleration factor for fast diffusion on purkinje network. */
Real diffusion_coeff = 0.8;
Real bias_coeff = 0.0;
/** Electrophysiology parameters. */
std::string diffusion_species_name = "Phi";
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.01;
Real b = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.002;
/** Fibers and sheet. */
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);
/** Purkinje Network. */
Vec3d starting_point(-21.9347 * length_scale, 4.0284 * length_scale, 0.0 * length_scale);
Vec3d second_point(-21.9347 * length_scale, 4.0284 * length_scale, -1.1089 * length_scale);

namespace SPH
{
//----------------------------------------------------------------------
//	Define heart shape
//----------------------------------------------------------------------
class Heart : public ComplexShape
{
  public:
    explicit Heart(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd translation(0.0, 0.0, 0.0);
        add<TriangleMeshShapeSTL>(full_path_to_lv, translation, length_scale);
    }
};
/** Set diffusion relaxation. */
using FiberDirectionDiffusionRelaxation =
    DiffusionRelaxationRK2<DiffusionRelaxation<Inner<KernelGradientInner>, IsotropicDiffusion>>;
/** Imposing diffusion boundary condition */
class DiffusionBCs : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit DiffusionBCs(BodyPartByParticle &body_part, const std::string &species_name)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(species_name)) {};
    virtual ~DiffusionBCs() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd displacement = sph_body_->getInitialShape().findNormalDirection(pos_[index_i]);
        Vecd face_norm = displacement / (displacement.norm() + 1.0e-15);

        Vecd center_norm = pos_[index_i] / (pos_[index_i].norm() + 1.0e-15);

        Real angle = face_norm.dot(center_norm);
        if (angle >= 0.0)
        {
            phi_[index_i] = 1.0;
        }
        else
        {
            if (pos_[index_i][1] < -getSPHAdaptation().ReferenceSpacing())
                phi_[index_i] = 0.0;
        }
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};

/** Compute Fiber and Sheet direction after diffusion */
class ComputeFiberAndSheetDirections : public LocalDynamics
{
  protected:
    LocallyOrthotropicMuscle &muscle_material_;
    Vecd *pos_;
    Real *phi_;
    Real beta_epi_, beta_endo_;
    Vecd center_line_vector_; // parallel to the ventricular centerline and pointing  apex-to-base

  public:
    explicit ComputeFiberAndSheetDirections(SPHBody &sph_body, const std::string &species_name)
        : LocalDynamics(sph_body),
          muscle_material_(DynamicCast<LocallyOrthotropicMuscle>(this, sph_body_->getBaseMaterial())),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(species_name))
    {
        center_line_vector_ = Vecd(0.0, 1.0, 0.0);
        beta_epi_ = -(70.0 / 180.0) * M_PI;
        beta_endo_ = (80.0 / 180.0) * M_PI;
    };
    virtual ~ComputeFiberAndSheetDirections() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        /**
         * Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
         * 		Present  doi.org/10.1016/j.cma.2016.05.031
         */
        /** Probe the face norm from Level set field. */
        Vecd displacement = sph_body_->getInitialShape().findNormalDirection(pos_[index_i]);
        Vecd face_norm = displacement / (displacement.norm() + 1.0e-15);
        Vecd center_norm = pos_[index_i] / (pos_[index_i].norm() + 1.0e-15);
        if (face_norm.dot(center_norm) <= 0.0)
        {
            face_norm = -face_norm;
        }
        /** Compute the centerline's projection on the plane orthogonal to face norm. */
        Vecd circumferential_direction = getCrossProduct(center_line_vector_, face_norm);
        Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
        /** The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo */
        Real beta = (beta_epi_ - beta_endo_) * phi_[index_i] + beta_endo_;
        /** Compute the rotation matrix through Rodrigues rotation formulation. */
        Vecd f_0 = cos(beta) * cd_norm + sin(beta) * getCrossProduct(face_norm, cd_norm) +
                   face_norm.dot(cd_norm) * (1.0 - cos(beta)) * face_norm;

        if (pos_[index_i][2] < 2.0 * getSPHAdaptation().ReferenceSpacing())
        {
            muscle_material_.local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
            muscle_material_.local_s0_[index_i] = face_norm;
        }
        else
        {
            muscle_material_.local_f0_[index_i] = Vecd::Zero();
            muscle_material_.local_s0_[index_i] = Vecd::Zero();
        }
    };
};
//	define shape parameters which will be used for the constrained body part.
class MuscleBaseShapeParameters : public TriangleMeshShapeBrick::ShapeParameters
{
  public:
    MuscleBaseShapeParameters() : TriangleMeshShapeBrick::ShapeParameters()
    {
        Real l = domain_upper_bound[0] - domain_lower_bound[0];
        Real w = domain_upper_bound[1] - domain_lower_bound[1];
        Real h = domain_upper_bound[2];
        halfsize_ = Vec3d(0.5 * l, 0.5 * w, 0.5 * h);
        resolution_ = 20;
        translation_ = Vec3d(-25.0 * length_scale, -5.0 * length_scale, 0.5 * h * length_scale);
    }
};
/**
 * application dependent initial condition
 */
class ApplyStimulusCurrentToMyocardium : public LocalDynamics
{
  public:
    explicit ApplyStimulusCurrentToMyocardium(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          voltage_(particles_->registerStateVariableData<Real>("Voltage")) {};

    void update(size_t index_i, Real dt)
    {
        if (-32.0 * length_scale <= pos_[index_i][0] && pos_[index_i][0] <= -20.0 * length_scale)
        {
            if (-5.0 * length_scale <= pos_[index_i][1] && pos_[index_i][1] <= 5.0)
            {
                if (-10.0 * length_scale <= pos_[index_i][2] && pos_[index_i][2] <= 0.0 * length_scale)
                {
                    voltage_[index_i] = 0.92;
                }
            }
        }
    };

  protected:
    Vecd *pos_;
    Real *voltage_;
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    observation_points.push_back(Vecd(-45.0 * length_scale, -30.0 * length_scale, 0.0));
    observation_points.push_back(Vecd(0.0, -30.0 * length_scale, 26.0 * length_scale));
    observation_points.push_back(Vecd(-30.0 * length_scale, -50.0 * length_scale, 0.0));
    observation_points.push_back(Vecd(0.0, -50.0 * length_scale, 20.0 * length_scale));
    observation_points.push_back(Vecd(0.0, -70.0 * length_scale, 0.0));
    return observation_points;
};

/**
 * application dependent initial condition
 */
class ApplyStimulusCurrentToPKJ : public LocalDynamics
{
  public:
    explicit ApplyStimulusCurrentToPKJ(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          voltage_(particles_->registerStateVariableData<Real>("Voltage")) {};

    void update(size_t index_i, Real dt)
    {
        if (index_i <= 10)
        {
            voltage_[index_i] = 1.0;
        }
    };

  protected:
    Vecd *pos_;
    Real *voltage_;
};

/**
 * Derived network particle generator.
 */
class NetworkWithExtraCheck;
template <>
class ParticleGenerator<BaseParticles, NetworkWithExtraCheck>
    : public ParticleGenerator<BaseParticles, Network>
{
  protected:
    bool extraCheck(const Vecd &new_position) override
    {
        bool no_generation = false;
        if (new_position[2] > 0)
            no_generation = true;
        return no_generation;
    };

  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles,
                      Vecd starting_pnt, Vecd second_pnt, int iterator, Real grad_factor)
        : ParticleGenerator<BaseParticles, Network>(
              sph_body, base_particles, starting_pnt, second_pnt, iterator, grad_factor) {};
};
} // namespace SPH
