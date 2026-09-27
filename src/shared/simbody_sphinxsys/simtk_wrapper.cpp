#include "simtk_wrapper.h"

#include "io_environment.h"
#include "sph_system.h"
#include "state_engine.h"
namespace SPH
{
//=================================================================================================//
SimbodyState::SimbodyState()
    : initial_origin_location_(Vec3d::Zero()), origin_location_(Vec3d::Zero()),
      origin_velocity_(Vec3d::Zero()), origin_acceleration_(Vec3d::Zero()),
      angular_velocity_(Vec3d::Zero()), angular_acceleration_(Vec3d::Zero()),
      rotation_(Mat3d::Identity()) {}
//=================================================================================================//
SimbodyState::SimbodyState(
    const SimTKVec3 &sim_tk_initial_origin_location,
    SimTK::MobilizedBody &mobod, const SimTK::State &state)
    : initial_origin_location_(SimTKToEigen(sim_tk_initial_origin_location)),
      origin_location_(SimTKToEigen(mobod.getBodyOriginLocation(state))),
      origin_velocity_(SimTKToEigen(mobod.getBodyOriginVelocity(state))),
      origin_acceleration_(SimTKToEigen(mobod.getBodyOriginAcceleration(state))),
      angular_velocity_(SimTKToEigen(mobod.getBodyAngularVelocity(state))),
      angular_acceleration_(SimTKToEigen(mobod.getBodyAngularAcceleration(state))),
      rotation_(SimTKToEigen(mobod.getBodyRotation(state))) {}
//=================================================================================================//
void SimbodyState::printSimbodyState()
{
    std::cout << "initial_origin_location_: " << initial_origin_location_.transpose() << std::endl;
    std::cout << "origin_location_: " << origin_location_.transpose() << std::endl;
    std::cout << "origin_velocity_: " << origin_velocity_.transpose() << std::endl;
    std::cout << "origin_acceleration_: " << origin_acceleration_.transpose() << std::endl;
    std::cout << "angular_velocity_: " << angular_velocity_.transpose() << std::endl;
    std::cout << "angular_acceleration_: " << angular_acceleration_.transpose() << std::endl;
    std::cout << "rotation_: \n"
              << rotation_ << std::endl;
}
//=================================================================================================//
SimbodySystem::SimbodySystem()
    : simbody_matter_(MBsystem_), integ_(MBsystem_)
{
    state_engine_keeper_.createPtr<SimbodyStateEngine>(MBsystem_);
}
//=================================================================================================//
SimbodyStateEngine &SimbodySystem::getSimbodyStateEngine()
{
    return *state_engine_keeper_.getPtr();
}
//=================================================================================================//
SimTK::Body::Rigid &SimbodySystem::createRigidBody(
    const std::string &name, const SimTK::MassProperties &mass_properties)
{
    SimTK::Body::Rigid *rigid_body =
        rigid_bodies_keeper_.createPtr<SimTK::Body::Rigid>(mass_properties);
    rigid_bodies_.push_back(std::make_pair(name, rigid_body));
    return *rigid_body;
}
//=================================================================================================//
SimTK::Body::Rigid &SimbodySystem::getRigidBody(const std::string &name)
{
    for (size_t i = 0; i < rigid_bodies_.size(); ++i)
    {
        if (rigid_bodies_[i].first == name)
        {
            return *rigid_bodies_[i].second;
        }
    }
    throw std::runtime_error("SimbodySystem::getRigidBody: Rigid body with name " + name + " not found.");
}
//=============================================================================================//
} // namespace SPH
