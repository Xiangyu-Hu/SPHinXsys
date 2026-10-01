/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4                                              *
 *                                                                           *
 * Portions copyright (c) 2017-2022 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	simbody_system.h
 * @brief 	tbd.
 * @author	Xiangyu Hu
 */
#ifndef SIMBODY_SYSTEM_H
#define SIMBODY_SYSTEM_H

#include "base_body_part.h"
#include "base_data_type.h"

namespace SimTK
{
class MobilizedBody;
class State;
} // namespace SimTK
namespace SPH
{
class RealBody;
class Shape;
class SimbodySystem;
class SolidBodyPartForSimbody;

class SolidBodyPartForSimbody : public BodyRegionByParticle
{
  public:
    SolidBodyPartForSimbody(SPHBody &body, Shape &body_part_shape);
    virtual ~SolidBodyPartForSimbody() {};

    Vec3d getInitialMassCenter() const { return initial_mass_center_; };
    Vec3d getInertiaMoments() const { return inertia_moments_; };
    Vec3d getInertiaProducts() const { return inertia_products_; };
    Real getTotalMass() const { return total_mass_; };

  protected:
    Real min_spacing_sqr_;
    Real total_mass_;
    Vec3d initial_mass_center_;
    Vec3d inertia_moments_;
    Vec3d inertia_products_;
    Real rho0_;
    Real *Vol_;
    Vecd *pos_;
};

struct SimbodyState
{
    Vec3d initial_origin_location_;
    Vec3d origin_location_, origin_velocity_, origin_acceleration_;
    Vec3d angular_velocity_, angular_acceleration_;
    Mat3d rotation_;

    SimbodyState();
    SimbodyState(SimTK::MobilizedBody &mobod, const SimTK::State &state);

    // implemented according to the Simbody API function with the same name
    void findStationLocationVelocityAndAccelerationInGround(
        const Vec3d &initial_location, const Vec3d &initial_normal,
        Vec3d &locationOnGround, Vec3d &velocityInGround,
        Vec3d &accelerationInGround, Vec3d &normalInGround)
    {
        Vec3d temp_location = rotation_ * (initial_location - initial_origin_location_);
        locationOnGround = origin_location_ + temp_location;

        Vec3d temp_velocity = angular_velocity_.cross(temp_location);
        velocityInGround = origin_velocity_ + temp_velocity;
        accelerationInGround = origin_acceleration_ +
                               angular_acceleration_.cross(temp_location) +
                               angular_velocity_.cross(temp_velocity);
        normalInGround = rotation_ * initial_normal;
    }

    void printSimbodyState();
};

template <>
struct ZeroData<SimbodyState>
{
    static inline const SimbodyState value = SimbodyState();
};

class SimbodySystem
{
  public:
    SimbodySystem();
    virtual ~SimbodySystem();
    void writeStateToXml(UnsignedInt iteration_step);
    void readStateFromXml(UnsignedInt iteration_step);
    std::string createRigidBody(SolidBodyPartForSimbody &simbody_part);
    std::string createFirstMobilizedPlanar(const std::string &body_name);
    std::string createFirstMobilizedPin(const std::string &body_name);
    void setUForMobilizedPlanar(const std::string &body_name, const Vec2d &velocity, Real angular_velocity);
    void setUForMobilizedPin(const std::string &body_name, Real angular_velocity);
    UnsignedInt getBodyIndexByName(const std::string &body_name);
    void realizeState();
    void initializeStateForIntegrator();
    void checkSimbodyState(const std::string &body_name);
    SimbodyState getSimbodyState(UnsignedInt body_index);
    Vec3d getSimbodyOriginLocation(UnsignedInt body_index);
    Real getSimbodySystemTime();
    void stepSimbodySystemTo(Real time);

  protected:
    class Impl;
    std::unique_ptr<Impl> impl_;
};
} // namespace SPH
#endif // SIMBODY_SYSTEM_H
