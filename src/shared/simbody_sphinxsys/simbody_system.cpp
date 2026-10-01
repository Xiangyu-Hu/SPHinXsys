#include "simbody_system.h"

#include "base_body.h"
#include "base_geometry.h"
#include "body_part_for_simbody.h"
#include "io_environment.h"
#include "simtk_wrapper.h"
#include "sph_system.h"
#include "sphinxsys_entity.h"
#include "state_engine.h"

namespace SPH
{
//=================================================================================================//
class SimbodySystem::Impl
{
    SimTK::MultibodySystem MBsystem_;
    SimTK::SimbodyMatterSubsystem matter_{MBsystem_};
    SimTK::RungeKuttaMersonIntegrator integ_{MBsystem_};
    SimbodyStateEngine state_engine_{MBsystem_};
    SimTK::State state_;
    UniquePtrsKeeper<SolidBodyPartForSimbody> body_parts_keeper_;
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

    std::string createRigidBody(RealBody &real_body, Shape &shape);

    SimTK::Body::Rigid &getRigidBody(const std::string &name);
    SolidBodyPartForSimbody &getSimbodyPart(const std::string &name);

    template <class MobilizedBodyType, class ParentBodyType>
    MobilizedBodyType &createMobilizedBody(
        const std::string &name, ParentBodyType &parent_mobod,
        const SimTK::Transform &X_PF, const SimTK::Body::Rigid &body,
        const SimTK::Transform &X_BM);

    template <class MobilizedBodyType>
    MobilizedBodyType &getMobilizedBody(const std::string &name);
    SimTK::MobilizedBody &getBaseMobilizedBody(const std::string &name);
};
//=================================================================================================//
std::string SimbodySystem::Impl::createRigidBody(RealBody &real_body, Shape &shape)
{
    SolidBodyPartForSimbody &body_part = real_body.addBodyPart<SolidBodyPartForSimbody>(shape);
    std::string part_name = body_part.Name();
    body_parts_.push_back(std::make_pair(part_name, &body_part));

    SimTK::Body::Rigid *rigid_body = rigid_bodies_keeper_.createPtr<
        SimTK::Body::Rigid>(body_part.getSimTKMassProperties());
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
    throw std::runtime_error("SimbodySystem::getRigidBody: Rigid body with name " + name + " not found.");
}
//=================================================================================================//
template <class MobilizedBodyType, class ParentBodyType>
MobilizedBodyType &SimbodySystem::Impl::createMobilizedBody(
    const std::string &name, ParentBodyType &parent_mobod,
    const SimTK::Transform &X_PF, const SimTK::Body::Rigid &body,
    const SimTK::Transform &X_BM)
{
    MobilizedBodyType *mobilized_body =
        mobilized_bodies_keeper_.createPtr<MobilizedBodyType>(
            parent_mobod, X_PF, body, X_BM);
    config_manager_.addEntity<MobilizedBodyType>(name, mobilized_body);
    mobilized_bodies_.push_back(std::make_pair(name, mobilized_body));
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
SimTK::MobilizedBody &SimbodySystem::Impl::getBaseMobilizedBody(const std::string &name)
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
std::string SimbodySystem::createRigidBody(RealBody &real_body, Shape &shape)
{
    return impl_->createRigidBody(real_body, shape);
}
//=================================================================================================//
std::string SimbodySystem::createFirstMobilizedPlanar(const std::string &body_name)
{
    SolidBodyPartForSimbody &simbody_part = impl_->getSimbodyPart(body_name);
    SimTK::Body::Rigid &rigid_body = impl_->getRigidBody(body_name);
    impl_->createMobilizedBody<SimTK::MobilizedBody::Planar>(
        body_name, impl_->getSimbodyMatterSubsystem().Ground(),
        simbody_part.getSimTKMassCenter(), rigid_body, simbody_part.getSimTKTransform());
    SimTK::State state = impl_->getMultibodySystem().realizeTopology();
    impl_->setState(state);
    return body_name;
}
//=================================================================================================//
void SimbodySystem::setUForMobilizedPlanar(
    const std::string &body_name, const Vec2d &velocity, Real angular_velocity)
{
    auto &mobilized_body =
        impl_->getMobilizedBody<SimTK::MobilizedBody::Planar>(body_name);
    SimTK::State &state = impl_->getState();
    mobilized_body.setU(state, SimTKVec3(angular_velocity, velocity[0], velocity[1]));
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
void SimbodySystem::checkSimbodyState(const std::string &body_name)
{
    auto &MBsystem = impl_->getMultibodySystem();
    auto &mobilized_body = impl_->getBaseMobilizedBody(body_name);
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
    SimbodyState test_simbody_state(mobilized_body.getBodyOriginLocation(state), mobilized_body, state);
    test_simbody_state.printSimbodyState();
}
//=================================================================================================//
} // namespace SPH
