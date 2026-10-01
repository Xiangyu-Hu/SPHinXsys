#include "simbody_system.h"

#include "adaptation.h"
#include "base_body.h"
#include "base_body_part.h"
#include "base_geometry.h"
#include "base_material.h"
#include "base_particles.hpp"
#include "ownership.h"
#include "simtk_wrapper.h"
#include "sphinxsys_entity.h"
#include "state_engine.h"
#include "vector_functions.h"

namespace SPH
{
//=================================================================================================//
SolidBodyPartForSimbody::
    SolidBodyPartForSimbody(SPHBody &body, Shape &body_part_shape)
    : BodyRegionByParticle(body, body_part_shape),
      min_spacing_sqr_(pow(body.getSPHAdaptation().MinimumSpacing(), 2)),
      total_mass_(0), initial_mass_center_(ZeroData<Vec3d>::value),
      inertia_moments_(ZeroData<Vec3d>::value), inertia_products_(ZeroData<Vec3d>::value),
      rho0_(DynamicCast<Solid>(this, body.getMatterMaterial()).ReferenceDensity()),
      Vol_(base_particles_.getVariableDataByName<Real>("VolumetricMeasure")),
      pos_(base_particles_.getVariableDataByName<Vecd>("Position"))
{
    Real body_part_volume(0);
    Vec3d mass_moment = Vec3d::Zero();
    for (size_t i = 0; i < body_part_particles_.size(); ++i)
    {
        size_t index_i = body_part_particles_[i];
        Vec3d particle_position = upgradeToVec3d(pos_[index_i]);
        Real particle_volume = Vol_[index_i];

        mass_moment += particle_volume * particle_position;
        body_part_volume += particle_volume;
    }
    initial_mass_center_ = mass_moment / body_part_volume;
    total_mass_ = rho0_ * body_part_volume;

    // computing unit inertia
    Vec3d total_inertia_moments = Vec3d::Zero();
    Vec3d total_inertia_products = Vec3d::Zero();
    for (size_t i = 0; i < body_part_particles_.size(); ++i)
    {
        size_t index_i = body_part_particles_[i];
        Vec3d particle_position = upgradeToVec3d(pos_[index_i]);
        Real particle_volume = Vol_[index_i];

        Vec3d displacement = particle_position - mass_moment;

        total_inertia_moments[0] +=
            particle_volume * (displacement[1] * displacement[1] + displacement[2] * displacement[2]);
        total_inertia_moments[1] +=
            particle_volume * (displacement[0] * displacement[0] + displacement[2] * displacement[2]);
        total_inertia_moments[2] +=
            particle_volume * (displacement[0] * displacement[0] + displacement[1] * displacement[1]);
        total_inertia_products[0] -= particle_volume * displacement[0] * displacement[1];
        total_inertia_products[1] -= particle_volume * displacement[0] * displacement[2];
        total_inertia_products[2] -= particle_volume * displacement[1] * displacement[2];
    }
    inertia_moments_ = total_inertia_moments / body_part_volume;
    inertia_products_ = total_inertia_products / body_part_volume;
    if constexpr (Dimensions == 2)
    {
        inertia_moments_[2] -= min_spacing_sqr_ / 6.0; // assume a small thickness
    }
}
//=================================================================================================//
class SimbodySystem::Impl
{
    SimTK::MultibodySystem MBsystem_;
    SimTK::SimbodyMatterSubsystem matter_{MBsystem_};
    SimTK::RungeKuttaMersonIntegrator integ_{MBsystem_};
    SimbodyStateEngine state_engine_{MBsystem_};
    SimTK::State state_;
    UniquePtrsKeeper<SimTK::Body::Rigid> rigid_bodies_keeper_;
    UniquePtrsKeeper<SimTK::MobilizedBody> mobilized_bodies_keeper_;
    StdVec<std::pair<std::string, SolidBodyPartForSimbody *>> body_parts_;
    StdVec<std::pair<std::string, SimTK::Body::Rigid *>> rigid_bodies_;
    StdVec<std::pair<std::string, SimTK::MobilizedBody *>> mobilized_bodies_;
    EntityManager config_manager_;

  public:
    SimTK::MultibodySystem &getMultibodySystem() { return MBsystem_; }
    SimTK::SimbodyMatterSubsystem &getSimbodyMatterSubsystem() { return matter_; }
    SimTK::RungeKuttaMersonIntegrator &getSimbodyIntegrator() { return integ_; }
    SimbodyStateEngine &getSimbodyStateEngine() { return state_engine_; }
    SimTK::State &getState() { return state_; }
    void setState(const SimTK::State &state) { state_ = state; }
    std::string createRigidBody(SolidBodyPartForSimbody &simbody_part);
    SimTK::Body::Rigid &getRigidBody(const std::string &name);
    SolidBodyPartForSimbody &getSimbodyPart(const std::string &name);
    UnsignedInt getBodyIndexByName(const std::string &body_name);

    template <class MobilizedBodyType, class ParentBodyType>
    MobilizedBodyType &createMobilizedBody(const std::string &name, ParentBodyType &parent_mobod);

    template <class MobilizedBodyType>
    MobilizedBodyType &getMobilizedBody(const std::string &name);

    SimTK::MobilizedBody &getMobilizedBody(const std::string &name);
    SimTK::MobilizedBody &getMobilizedBody(UnsignedInt body_index);
};
//=================================================================================================//
std::string SimbodySystem::Impl::createRigidBody(SolidBodyPartForSimbody &simbody_part)
{
    std::string part_name = simbody_part.Name();
    body_parts_.push_back(std::make_pair(part_name, &simbody_part));
    SimTK::UnitInertia inertia(
        EigenToSimTK(simbody_part.getInertiaMoments()),
        EigenToSimTK(simbody_part.getInertiaProducts()));
    SimTK::Body::Rigid *rigid_body = rigid_bodies_keeper_.createPtr<
        SimTK::Body::Rigid>(SimTK::MassProperties(
        simbody_part.getTotalMass(),
        EigenToSimTK(simbody_part.getInitialMassCenter()),
        inertia));
    rigid_bodies_.push_back(std::make_pair(part_name, rigid_body));
    return part_name;
}
//=================================================================================================//
SolidBodyPartForSimbody &SimbodySystem::Impl::getSimbodyPart(const std::string &name)
{
    for (size_t i = 0; i < body_parts_.size(); ++i)
    {
        if (body_parts_[i].first == name)
        {
            return *body_parts_[i].second;
        }
    }
    throw std::runtime_error("SimbodySystem::getSimbodyPart: Body part with name " + name + " not found.");
}
//=================================================================================================//
SimTK::Body::Rigid &SimbodySystem::Impl::getRigidBody(const std::string &name)
{
    for (size_t i = 0; i < rigid_bodies_.size(); ++i)
    {
        if (rigid_bodies_[i].first == name)
        {
            return *rigid_bodies_[i].second;
        }
    }
    throw std::runtime_error(
        "SimbodySystem::getRigidBody: Rigid body with name " + name + " not found.");
}
//=================================================================================================//
UnsignedInt SimbodySystem::Impl::getBodyIndexByName(const std::string &body_name)
{
    for (size_t i = 0; i < mobilized_bodies_.size(); ++i)
    {
        if (mobilized_bodies_[i].first == body_name)
        {
            return i;
        }
    }
    throw std::runtime_error(
        "SimbodySystem::getBodyIndexByName: Mobilized body with name " + body_name + " not found.");
}
//=================================================================================================//
template <class MobilizedBodyType, class ParentBodyType>
MobilizedBodyType &SimbodySystem::Impl::createMobilizedBody(
    const std::string &name, ParentBodyType &parent_mobod)
{
    SolidBodyPartForSimbody &simbody_part = getSimbodyPart(name);
    SimTK::Body::Rigid &rigid_body = getRigidBody(name);
    MobilizedBodyType *mobilized_body =
        mobilized_bodies_keeper_.createPtr<MobilizedBodyType>(
            getSimbodyMatterSubsystem().Ground(),
            EigenToSimTK(simbody_part.getInitialMassCenter()),
            rigid_body, SimTK::Transform());
    config_manager_.addEntity<MobilizedBodyType>(name, mobilized_body);
    mobilized_bodies_.push_back(std::make_pair(name, mobilized_body));
    setState(MBsystem_.realizeTopology());
    return *mobilized_body;
}
//=================================================================================================//
template <class MobilizedBodyType>
MobilizedBodyType &SimbodySystem::Impl::getMobilizedBody(const std::string &name)
{
    if (config_manager_.hasEntity<MobilizedBodyType>(name))
    {
        return config_manager_.getEntity<MobilizedBodyType>(name);
    }
    throw std::runtime_error(
        "SimbodySystem::getMobilizedBody: Mobilized body with name " + name + " not found.");
}
//=================================================================================================//
SimTK::MobilizedBody &SimbodySystem::Impl::getMobilizedBody(const std::string &name)
{
    for (size_t i = 0; i < mobilized_bodies_.size(); ++i)
    {
        if (mobilized_bodies_[i].first == name)
        {
            return *mobilized_bodies_[i].second;
        }
    }
    std::cerr << "\n Error: Mobilized body " << name << " not found!" << std::endl;
    throw std::runtime_error("Error: Mobilized body not found!");
}
//=================================================================================================//
SimTK::MobilizedBody &SimbodySystem::Impl::getMobilizedBody(UnsignedInt body_index)
{
    if (body_index < mobilized_bodies_.size())
    {
        return *mobilized_bodies_[body_index].second;
    }
    std::cerr << "\n Error: Mobilized body index " << body_index << " out of range!" << std::endl;
    throw std::runtime_error("Error: Mobilized body index out of range!");
}
//=================================================================================================//
SimbodyState::SimbodyState()
    : initial_origin_location_(Vec3d::Zero()), origin_location_(Vec3d::Zero()),
      origin_velocity_(Vec3d::Zero()), origin_acceleration_(Vec3d::Zero()),
      angular_velocity_(Vec3d::Zero()), angular_acceleration_(Vec3d::Zero()),
      rotation_(Mat3d::Identity()) {}
//=================================================================================================//
SimbodyState::SimbodyState(SimTK::MobilizedBody &mobod, const SimTK::State &state)
    : initial_origin_location_(SimTKToEigen(mobod.getBodyOriginLocation(state))),
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
SimbodySystem::SimbodySystem() : impl_(std::make_unique<Impl>()) {}
//=================================================================================================//
SimbodySystem::~SimbodySystem() = default;
//=================================================================================================//
void SimbodySystem::writeStateToXml(UnsignedInt iteration_step)
{
    impl_->getSimbodyStateEngine().writeStateToXml(iteration_step, impl_->getSimbodyIntegrator());
}
//=================================================================================================//
void SimbodySystem::readStateFromXml(UnsignedInt iteration_step)
{
    impl_->getSimbodyStateEngine().readStateFromXml(iteration_step, impl_->getState());
}
//=================================================================================================//
std::string SimbodySystem::createRigidBody(SolidBodyPartForSimbody &simbody_part)
{
    return impl_->createRigidBody(simbody_part);
}
//=================================================================================================//
std::string SimbodySystem::createFirstMobilizedPlanar(const std::string &body_name)
{
    auto &parent_body = impl_->getSimbodyMatterSubsystem().Ground();
    impl_->createMobilizedBody<SimTK::MobilizedBody::Planar>(body_name, parent_body);
    return body_name;
}
//=================================================================================================//
std::string SimbodySystem::createFirstMobilizedPin(const std::string &body_name)
{
    auto &parent_body = impl_->getSimbodyMatterSubsystem().Ground();
    impl_->createMobilizedBody<SimTK::MobilizedBody::Pin>(body_name, parent_body);
    return body_name;
}
//=================================================================================================//
void SimbodySystem::setUForMobilizedPlanar(
    const std::string &body_name, const Vec2d &velocity, Real angular_velocity)
{
    auto &mobilized_body = impl_->getMobilizedBody<SimTK::MobilizedBody::Planar>(body_name);
    SimTK::State &state = impl_->getState();
    mobilized_body.setU(state, SimTKVec3(angular_velocity, velocity[0], velocity[1]));
}
//=================================================================================================//
void SimbodySystem::setUForMobilizedPin(const std::string &body_name, Real angular_velocity)
{
    auto &mobilized_body = impl_->getMobilizedBody<SimTK::MobilizedBody::Pin>(body_name);
    SimTK::State &state = impl_->getState();
    mobilized_body.setU(state, angular_velocity);
}
//=================================================================================================//
void SimbodySystem::realizeState()
{
    SimTK::State &state = impl_->getState();
    impl_->getMultibodySystem().realize(state);
}
//=================================================================================================//
void SimbodySystem::initializeStateForIntegrator()
{
    SimTK::State &state = impl_->getState();
    impl_->getSimbodyIntegrator().initialize(state);
}
//=================================================================================================//
SimbodyState SimbodySystem::getSimbodyState(UnsignedInt body_index)
{
    auto &MBsystem = impl_->getMultibodySystem();
    const SimTK::State &state = impl_->getSimbodyIntegrator().getState();
    MBsystem.realize(state);
    auto &mobilized_body = impl_->getMobilizedBody(body_index);
    return SimbodyState(mobilized_body, state);
}
//=================================================================================================//
UnsignedInt SimbodySystem::getBodyIndexByName(const std::string &body_name)
{
    return impl_->getBodyIndexByName(body_name);
}
//=================================================================================================//
Vec3d SimbodySystem::getSimbodyOriginLocation(UnsignedInt body_index)
{
    const SimTK::State &state = impl_->getSimbodyIntegrator().getState();
    auto &mobilized_body = impl_->getMobilizedBody(body_index);
    return SimTKToEigen(mobilized_body.getBodyOriginLocation(state));
}
//=================================================================================================//
Real SimbodySystem::getSimbodySystemTime()
{
    return impl_->getSimbodyIntegrator().getState().getTime();
}
//=================================================================================================//
void SimbodySystem::stepSimbodySystemTo(Real time)
{
    auto &integ = impl_->getSimbodyIntegrator();
    integ.stepTo(time);
}
//=================================================================================================//
void SimbodySystem::checkSimbodyState(const std::string &body_name)
{
    auto &MBsystem = impl_->getMultibodySystem();
    auto &mobilized_body = impl_->getMobilizedBody(body_name);
    auto &simbody_body = impl_->getRigidBody(body_name);
    auto &mass_properties = simbody_body.getDefaultRigidBodyMassProperties();
    std::cout << "\n------------------------------------------------------------" << std::endl;
    std::cout << "Simbody constraint information: " << std::endl;
    std::cout << "Mass: " << mass_properties.getMass() << std::endl;
    std::cout << "UnitInertia Moments: " << mass_properties.getUnitInertia().getMoments() << std::endl;
    std::cout << "UnitInertia Products: " << mass_properties.getUnitInertia().getProducts() << std::endl;
    std::cout << "------------------------------------------------------------" << std::endl;

    auto &integ = impl_->getSimbodyIntegrator();
    SimTK::State state = integ.getState(); // copy to allow cache invalidation
    MBsystem.realize(state);
    SimbodyState test_simbody_state(mobilized_body, state);
    test_simbody_state.printSimbodyState();
}
//=================================================================================================//
} // namespace SPH
