#include "simbody_system.h"

#include "adaptation.h"
#include "base_body.h"
#include "base_body_part.h"
#include "base_geometry.h"
#include "base_material.h"
#include "base_particles.hpp"
#include "ownership.h"
#include "simtk_wrapper.h"
#include "state_engine.h"
#include "vector_functions.h"

#include "body_part_for_simbody.h"

namespace SPH
{
//=================================================================================================//
SolidBodyPartForSimbodyCK::
    SolidBodyPartForSimbodyCK(SPHBody &body, Shape &body_part_shape)
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
    SimTK::GeneralForceSubsystem force_system_{MBsystem_};
    SimTK::Force::DiscreteForces force_on_bodies_{force_system_, matter_};
    SimTK::RungeKuttaMersonIntegrator integ_{MBsystem_};
    SimbodyStateEngine state_engine_{MBsystem_};
    SimTK::State initial_state_for_integrator_;
    UniquePtrsKeeper<SimTK::Body::Rigid> rigid_bodies_keeper_;
    UniquePtrsKeeper<SimTK::Force> forces_keeper_;
    UniquePtrsKeeper<SimTK::MobilizedBody> mobilized_bodies_keeper_;
    StdVec<std::pair<std::string, SolidBodyPartForSimbodyCK *>> body_parts_;
    StdVec<std::pair<std::string, SimTK::Force *>> forces_;
    StdVec<std::pair<std::string, SimTK::Body::Rigid *>> rigid_bodies_;
    StdVec<std::pair<std::string, SimTK::MobilizedBody *>> mobilized_bodies_;

  public:
    Impl();
    SimTK::MultibodySystem &getMultibodySystem() { return MBsystem_; }
    SimTK::SimbodyMatterSubsystem &getSimbodyMatterSubsystem() { return matter_; }
    SimTK::GeneralForceSubsystem &getSimbodyForceSubsystem() { return force_system_; }
    SimTK::Force::DiscreteForces &getSimbodyForceOnBodies() { return force_on_bodies_; }
    SimTK::RungeKuttaMersonIntegrator &getSimbodyIntegrator() { return integ_; }
    SimbodyStateEngine &getSimbodyStateEngine() { return state_engine_; }
    SimTK::State &getInitialStateForIntegrator() { return initial_state_for_integrator_; }
    void setInitialStateForIntegrator(const SimTK::State &state) { initial_state_for_integrator_ = state; }

    SimTK::Body &createRigidBody(const std::string &name, SimTK::MassProperties mass_properties);
    SimTK::Body::Rigid &getRigidBody(const std::string &name);
    SolidBodyPartForSimbodyCK &getSimbodyPart(const std::string &name);
    UnsignedInt getBodyIndexByName(const std::string &body_name);

    template <class ForceType, typename... Args>
    ForceType &addForce(const std::string &name, Args &&...args);
    void updateForceOnBody(UnsignedInt body_index, const SimTK::SpatialVec &spatial_force);

    template <class MobilizedBodyType>
    MobilizedBodyType &createMobilizedBody(
        const std::string &name, SimTK::MobilizedBody &parent_mobod,
        const SimTK::Transform &X_PF, const SimTK::Transform &X_BM);

    SimTK::MobilizedBody &getMobilizedBody(const std::string &name);
    SimTK::MobilizedBody &getMobilizedBody(UnsignedInt body_index);
};
//=================================================================================================//
SimbodySystem::Impl::Impl() : initial_state_for_integrator_(MBsystem_.realizeTopology())
{
    mobilized_bodies_.push_back(std::make_pair("Ground", &matter_.Ground()));
}
//=================================================================================================//
SimTK::Body &SimbodySystem::Impl::createRigidBody(
    const std::string &name, SimTK::MassProperties mass_properties)
{
    SimTK::Body::Rigid *rigid_body = rigid_bodies_keeper_.createPtr<
        SimTK::Body::Rigid>(mass_properties);
    rigid_bodies_.push_back(std::make_pair(name, rigid_body));
    return *rigid_body;
}
//=================================================================================================//
template <class ForceType, typename... Args>
ForceType &SimbodySystem::Impl::addForce(const std::string &name, Args &&...args)
{
    ForceType *force = forces_keeper_.createPtr<ForceType>(
        force_system_, matter_, std::forward<Args>(args)...);
    forces_.push_back(std::make_pair(name, force));
    initial_state_for_integrator_ = MBsystem_.realizeTopology();
    return *force;
}
//=================================================================================================//
void SimbodySystem::Impl::updateForceOnBody(
    UnsignedInt body_index, const SimTK::SpatialVec &spatial_force)
{
    SimTK::State &state_for_update = integ_.updAdvancedState();
    auto &mobilized_body = getMobilizedBody(body_index);
    force_on_bodies_.setOneBodyForce(state_for_update, mobilized_body, spatial_force);
}
//=================================================================================================//
SolidBodyPartForSimbodyCK &SimbodySystem::Impl::getSimbodyPart(const std::string &name)
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
template <class MobilizedBodyType>
MobilizedBodyType &SimbodySystem::Impl::createMobilizedBody(
    const std::string &name, SimTK::MobilizedBody &parent_mobod,
    const SimTK::Transform &X_PF, const SimTK::Transform &X_BM)
{
    SimTK::Body::Rigid &rigid_body = getRigidBody(name);
    MobilizedBodyType *mobilized_body = mobilized_bodies_keeper_.createPtr<
        MobilizedBodyType>(parent_mobod, X_PF, rigid_body, X_BM);
    mobilized_bodies_.push_back(std::make_pair(name, mobilized_body));
    initial_state_for_integrator_ = MBsystem_.realizeTopology();
    return *mobilized_body;
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
SimbodyState::SimbodyState(
    const Vec3d &initial_origin_location, SimTK::MobilizedBody &mobod, const SimTK::State &state)
    : initial_origin_location_(initial_origin_location),
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
    SimTK::State state = impl_->getMultibodySystem().realizeTopology();
    impl_->getSimbodyStateEngine().readStateFromXml(iteration_step, state);
    impl_->getMultibodySystem().realize(state);
    impl_->getSimbodyIntegrator().initialize(state);
}
//=================================================================================================//
std::string SimbodySystem::createRigidBody(SolidBodyPartForSimbodyCK &simbody_part)
{
    std::string part_name = simbody_part.Name();
    config_manager_.addEntity<SolidBodyPartForSimbodyCK>(part_name, &simbody_part);
    SimTK::MassProperties mass_properties(
        simbody_part.getTotalMass(), EigenToSimTK(simbody_part.getInitialMassCenter()),
        SimTK::UnitInertia(EigenToSimTK(simbody_part.getInertiaMoments()),
                           EigenToSimTK(simbody_part.getInertiaProducts())));
    impl_->createRigidBody(part_name, mass_properties);
    return part_name;
}
//=================================================================================================//
void SimbodySystem::addUniformGravity(const Vec3d &gravity_vector)
{
    std::string force_name = "UniformGravity";
    auto &gravity = impl_->addForce<SimTK::Force::UniformGravity>(force_name, EigenToSimTK(gravity_vector));
    config_manager_.addEntity<SimTK::Force::UniformGravity>(force_name, &gravity);
}
//=================================================================================================//
SimTK::MobilizedBody &SimbodySystem::createMobilizedBody(
    const std::string &body_name, const std::string &parent_name, const std::string &mobilizer_type)
{
    auto &parent_mobod = impl_->getMobilizedBody(parent_name);
    auto &body_part = config_manager_.getEntity<SolidBodyPartForSimbodyCK>(body_name);
    SimTK::Vec3 initial_mass_center = EigenToSimTK(body_part.getInitialMassCenter());
    
    if (mobilizer_type == "Planar")
    {
        auto &mobilized_body = impl_->createMobilizedBody<SimTK::MobilizedBody::Planar>(
            body_name, parent_mobod, SimTK::Transform(initial_mass_center), SimTK::Transform());
        config_manager_.addEntity<SimTK::MobilizedBody::Planar>(body_name, &mobilized_body);
        return mobilized_body;
    }
    else if (mobilizer_type == "Pin")
    {
        auto &mobilized_body = impl_->createMobilizedBody<SimTK::MobilizedBody::Pin>(
            body_name, parent_mobod, SimTK::Transform(initial_mass_center), SimTK::Transform());
        config_manager_.addEntity<SimTK::MobilizedBody::Pin>(body_name, &mobilized_body);
        return mobilized_body;
    }
    else
    {
        throw std::runtime_error(
            "SimbodySystem::createMobilizedBody: Unsupported mobilizer type " + mobilizer_type);
    }
}
//=================================================================================================//
void SimbodySystem::updateForceOnBody(
    SimTK::MobilizedBody &mobilized_body, const TorqueAndForce &torque_and_force)
{
    SimTK::State &state_for_update = impl_->getSimbodyIntegrator().updAdvancedState();
    auto &force_on_bodies = impl_->getSimbodyForceOnBodies();
    SimTK::SpatialVec spatial_force(
        EigenToSimTK(torque_and_force.first), EigenToSimTK(torque_and_force.second));
    force_on_bodies.setOneBodyForce(state_for_update, mobilized_body, spatial_force);
}
//=================================================================================================//
void SimbodySystem::updateForceOnBody(UnsignedInt body_index, const TorqueAndForce &torque_and_force)
{
    SimTK::SpatialVec spatial_force(
        EigenToSimTK(torque_and_force.first), EigenToSimTK(torque_and_force.second));
    impl_->updateForceOnBody(body_index, spatial_force);
}
//=================================================================================================//
void SimbodySystem::setUForMobilizedPlanar(
    const std::string &body_name, const Vec2d &velocity, Real angular_velocity)
{
    auto &mobilized_body = config_manager_.getEntity<SimTK::MobilizedBody::Planar>(body_name);
    SimTK::State &state = impl_->getInitialStateForIntegrator();
    mobilized_body.setU(state, SimTKVec3(angular_velocity, velocity[0], velocity[1]));
}
//=================================================================================================//
void SimbodySystem::setUForMobilizedPin(const std::string &body_name, Real angular_velocity)
{
    auto &mobilized_body = config_manager_.getEntity<SimTK::MobilizedBody::Pin>(body_name);
    SimTK::State &state = impl_->getInitialStateForIntegrator();
    mobilized_body.setU(state, angular_velocity);
}
//=================================================================================================//
void SimbodySystem::initializeStateForIntegrator(Real accuracy, bool allow_interpolation)
{
    SimTK::State state = impl_->getInitialStateForIntegrator();
    impl_->getMultibodySystem().realize(state);
    auto &integ = impl_->getSimbodyIntegrator();
    integ.initialize(state);
    integ.setAccuracy(accuracy);
    integ.setAllowInterpolation(allow_interpolation);
}
//=================================================================================================//
SimbodyState SimbodySystem::getSimbodyState(
    const Vec3d &initial_origin_location, UnsignedInt body_index)
{
    auto &MBsystem = impl_->getMultibodySystem();
    const SimTK::State &state = impl_->getSimbodyIntegrator().getState();
    MBsystem.realize(state);
    auto &mobilized_body = impl_->getMobilizedBody(body_index);
    return SimbodyState(initial_origin_location, mobilized_body, state);
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
Vec3d SimbodySystem::getInitialSimbodyOriginLocation(UnsignedInt body_index)
{
    auto &MBsystem = impl_->getMultibodySystem();
    const SimTK::State &state = MBsystem.getDefaultState();
    MBsystem.realize(state);
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
void SimbodySystem::stepSimbodySystemBy(Real dt)
{
    auto &integ = impl_->getSimbodyIntegrator();
    integ.stepBy(dt);
}
//=================================================================================================//
void SimbodySystem::checkInitialSimbodyState(const std::string &body_name)
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
    Vec3d initial_origin_location = SimTKToEigen(mobilized_body.getBodyOriginLocation(state));
    SimbodyState test_simbody_state(initial_origin_location, mobilized_body, state);
    test_simbody_state.printSimbodyState();
}
//=================================================================================================//
} // namespace SPH