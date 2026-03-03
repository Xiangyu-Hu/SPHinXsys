#include "state_engine.h"

#include "sph_system.h"
namespace SPH
{
//=================================================================================================//
SimbodyStateEngine::
    SimbodyStateEngine(SPHSystem &sph_system, SimTK::MultibodySystem &system)
    : simbody_xml_engine_("state_xml", "mbsystem")
{
    mbsystem_ = system;
    restart_folder_ = sph_system.getIOEnvironment().RestartFolder();
    if (!fs::exists(restart_folder_))
    {
        fs::create_directory(restart_folder_);
    }
}
//=============================================================================================//
void SimbodyStateEngine::InitializeState()
{
    /** Clear cached list of all related
        StateVariables if any from a previousSystem.
     */
    allstatevariables_.clear();
    getMultibodySystem().invalidateSystemTopologyCache();
    getMultibodySystem().realizeTopology();
    /** Set the model's operating state (internal member variable) to the
        default state that is stored inside the System.
     */
    working_state_ = getMultibodySystem().getDefaultState();
    /** Process the modified modeling option. */
    getMultibodySystem().realizeModel(working_state_);
    /** Realize instance variables that may have been set above. This
     *  means floating point parameters such as mass properties and
     * geometry placements are frozen.
     */
    getMultibodySystem().realize(working_state_, SimTK::Stage::Instance);
    /** Realize the initial configuration in preparation. This
     * initial configuration does not necessarily satisfy constraints.
     */
    getMultibodySystem().realize(working_state_, SimTK::Stage::Position);
}
//=============================================================================================//
SimTK::MultibodySystem &SimbodyStateEngine::getMultibodySystem()
{
    return mbsystem_.getRef();
}
//=============================================================================================//
void SimbodyStateEngine::addStateVariable(std::string statevariablename,
                                          SimTK::Stage invalidatestage)
{
    if ((invalidatestage < SimTK::Stage::Position) ||
        (invalidatestage > SimTK::Stage::Dynamics))
    {
        std::stringstream msg;
        msg << "SimbodyStateEngine::addStateVariable: invalidatestage "
               "must be Position, Velocity or Dynamics. "
            << __FILE__ << __LINE__;
        throw(msg.str());
    }
    /** Allocate space for a new state variable. */
    AddedStateVariable *asv =
        added_state_variable_keeper_
            .createPtr<AddedStateVariable>(statevariablename, *this, invalidatestage);
    // Add it to the Component and let it take ownership
    addStateVariable(asv);
}
//=============================================================================================//
void SimbodyStateEngine::addStateVariable(SimbodyStateEngine::StateVariable *statevariable)
{
    std::string &statevariablename = statevariable->getName();
    /** don't add state if there is another state variable with the same name. */
    std::map<std::string, StateVariableInfo>::const_iterator it;
    it = namedstatevariableinfo_.find(statevariablename);
    if (it != namedstatevariableinfo_.end())
    {
        std::stringstream msg;
        msg << "SimbodyStateEngine::addStateVariable: State variable " << statevariablename << " already exists." << __FILE__ << __LINE__;
        throw(msg.str());
    }
    int order = (int)namedstatevariableinfo_.size();
    /** assign a "slot" for a state variable by name
        state variable index will be invalid by default
        upon allocation during realizeTopology the index will be set
     */
    namedstatevariableinfo_[statevariablename] = StateVariableInfo(statevariable, order);
}
//=============================================================================================//
SimbodyStateEngine::StateVariable *SimbodyStateEngine::
    traverseToStateVariable(std::string &pathname)
{
    auto it = namedstatevariableinfo_.find(pathname);
    if (it != namedstatevariableinfo_.end())
    {
        return it->second.statevariable_.get();
    }
    else
    {
        return nullptr;
    }
}
//=============================================================================================//
StdVec<std::string> SimbodyStateEngine::getStateVariableNames()
{
    std::map<std::string, StateVariableInfo>::const_iterator it;
    it = namedstatevariableinfo_.begin();

    StdVec<std::string> names; //("",(int)namedstatevariableinfo_.size());

    while (it != namedstatevariableinfo_.end())
    {
        names[it->second.order] = it->first;
        it++;
    }
    return names;
}
//=============================================================================================//
int SimbodyStateEngine::getNumOfStateVariables()
{
    return getNumStateVariablesAddedByEngine();
}
//=============================================================================================//
bool SimbodyStateEngine::isAllStatesVariablesListValid()
{
    int nsv = getNumOfStateVariables();
    /** Consider the list of all StateVariables to be valid if all of
        the following conditions are true:
        1. a System has been associated with the list of StateVariables
        2. The list of all StateVariables is correctly sized (initialized)
        3. The System associated with the StateVariables is the current System */
    bool valid =
        !statesassociatedsystem_.empty() &&
        (int)allstatevariables_.size() == nsv &&
        getMultibodySystem().isSameSystem(statesassociatedsystem_.getRef());

    return valid;
}
//=============================================================================================//
SimTK::Vector SimbodyStateEngine::getStateVariableValues()
{
    int nsv = getNumOfStateVariables();
    /** if the StateVariables are invalid, rebuild the list. */
    if (!isAllStatesVariablesListValid())
    {
        statesassociatedsystem_.reset(&getMultibodySystem());
        allstatevariables_.clear();
        allstatevariables_.resize(nsv);
        StdVec<std::string> names = getStateVariableNames();
        for (int i = 0; i < nsv; ++i)
            allstatevariables_[i].reset(traverseToStateVariable(names[i]));
    }

    SimTK::Vector statevariablevalues(nsv, SimTK::NaN);
    for (int i = 0; i < nsv; ++i)
    {
        statevariablevalues[i] = allstatevariables_[i]->getValue();
        std::cout << statevariablevalues[i] << std::endl;
    }
    return statevariablevalues;
}
//=============================================================================================//
Real SimbodyStateEngine::AddedStateVariable::getValue()
{
    SimTK::ZIndex zix(getVarIndex());
    if (getSubsysIndex().isValid() && zix.isValid())
    {
        const SimTK::Vector &z = getOwner().getDefaultSubsystem().getZ(getOwner().working_state_);
        return z[SimTK::ZIndex(zix)];
    }

    std::stringstream msg;
    msg << "SimbodyStateEngine::AddedStateVariable::getValue: ERR- variable '"
        << getName() << "' is invalid! " << __FILE__ << __LINE__;
    throw(msg.str());
    return SimTK::NaN;
}
//=============================================================================================//
void SimbodyStateEngine::AddedStateVariable::setValue(Real value)
{
    SimTK::ZIndex zix(getVarIndex());
    if (getSubsysIndex().isValid() && zix.isValid())
    {
        SimTK::Vector &z = getOwner().getDefaultSubsystem().updZ(getOwner().working_state_);
        z[SimTK::ZIndex(zix)] = value;
        return;
    }

    std::stringstream msg;
    msg << "SimbodyStateEngine::AddedStateVariable::setValue: ERR- variable '"
        << getName() << "' is invalid! " << __FILE__ << __LINE__;
    ;
    throw(msg.str());
}
//=============================================================================================//
Real SimbodyStateEngine::AddedStateVariable::
    getDerivative()
{
    // return getCacheVariableValue<Real>(state, getName()+"_deriv");
    return 0.0;
}
//=============================================================================================//
void SimbodyStateEngine::AddedStateVariable::
    setDerivative(Real deriv)
{
    // return setCacheVariableValue<Real>(state, getName()+"_deriv", deriv);
}
//=============================================================================================//
void SimbodyStateEngine::reporter(SimTK::State &state_)
{
    const SimTK::SimbodyMatterSubsystem &matter_ = getMultibodySystem().getMatterSubsystem();
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter_.getNumBodies(); ++mbx)
    {

        const SimTK::MobilizedBody &mobod = matter_.getMobilizedBody(mbx);

        int num_q_ = mobod.getNumQ(state_);
        for (int i = 0; i < num_q_; i++)
        {
            std::cout << num_q_ << " " << mobod.getOneQ(state_, SimTK::QIndex(i)) << std::endl;
        }
        int num_u_ = mobod.getNumU(state_);
        for (int i = 0; i < num_u_; i++)
        {
            std::cout << num_u_ << " " << mobod.getOneU(state_, SimTK::UIndex(i)) << std::endl;
        }
        std::cout << " Body Info : " << std::endl;
        std::cout << " Transform : " << mobod.getBodyTransform(state_) << std::endl;
        std::cout << " Rotation : " << mobod.getBodyRotation(state_) << std::endl;
        std::cout << " Origin : " << mobod.getBodyOriginLocation(state_) << std::endl;
    }
}
//=============================================================================================//
void SimbodyStateEngine::resizeXmlDocForSimbody(size_t input_size)
{
    size_t total_elements = simbody_xml_engine_.SizeOfXmlDoc();

    if (total_elements <= input_size)
    {
        for (size_t i = total_elements; i != input_size; ++i)
            simbody_xml_engine_.addElementToXmlDoc("mbbody");
    }
}
//=============================================================================================//
void SimbodyStateEngine::writeStateToXml(int ite_rst, SimTK::RungeKuttaMersonIntegrator &integ)
{
    const SimTK::State &state = integ.getState();
    const SimTK::SimbodyMatterSubsystem &matter = getMultibodySystem().getMatterSubsystem();
    resizeXmlDocForSimbody(matter.getNumBodies());
    SimTK::Xml::element_iterator ele_ite = simbody_xml_engine_.root_element_.element_begin();
    for (SimTK::MobilizedBodyIndex mbx(0); mbx != matter.getNumBodies(); ++mbx)
    {
        const SimTK::MobilizedBody &mobod = matter.getMobilizedBody(mbx);

        int num_q_ = mobod.getNumQ(state);
        for (int i = 0; i < num_q_; i++)
        {
            Real mobod_q = mobod.getOneQ(state, SimTK::QIndex(i));
            std::string ele_name = "QIndx_" + std::to_string(i);
            simbody_xml_engine_.setAttributeToElement(ele_ite, ele_name, mobod_q);
        }

        int num_u_ = mobod.getNumU(state);
        for (int i = 0; i < num_u_; i++)
        {
            Real mobod_u = mobod.getOneU(state, SimTK::UIndex(i));
            std::string ele_name = "UIndx_" + std::to_string(i);
            simbody_xml_engine_.setAttributeToElement(ele_ite, ele_name, mobod_u);
        }
        ++ele_ite;
    }
    std::string filefullpath = restart_folder_ + "/simbody_rst_" + std::to_string(ite_rst) + ".xml";
    simbody_xml_engine_.writeToXmlFile(filefullpath);
}
//=============================================================================================//
void SimbodyStateEngine::readStateFromXml(int ite_rst, SimTK::State &state)
{
    std::string filefullpath = restart_folder_ + "/simbody_rst_" + std::to_string(ite_rst) + ".xml";
    const SimTK::SimbodyMatterSubsystem &matter = getMultibodySystem().getMatterSubsystem();
    if (!fs::exists(filefullpath))
    {
        std::cout << "\n Error: the input file:" << filefullpath << " is not valid" << std::endl;
        std::cout << __FILE__ << ':' << __LINE__ << std::endl;
        exit(1);
    }
    else
    {
        int num_mobod = 0;
        simbody_xml_engine_.loadXmlFile(filefullpath);
        SimTK::Xml::element_iterator ele_ite_ = simbody_xml_engine_.root_element_.element_begin();
        for (; ele_ite_ != simbody_xml_engine_.root_element_.element_end(); ++ele_ite_)
        {
            const SimTK::MobilizedBody &mobod = matter.getMobilizedBody(SimTK::MobilizedBodyIndex(num_mobod));
            int num_q_ = mobod.getNumQ(state);
            Real q_tmp_ = 0.0;
            if (num_q_ != 0)
            {
                for (int i = 0; i < num_q_; i++)
                {
                    std::string attr_name = "QIndx_" + std::to_string(i);
                    simbody_xml_engine_.getRequiredAttributeValue(ele_ite_, attr_name, q_tmp_);
                    mobod.setOneQ(state, SimTK::QIndex(i), q_tmp_);
                }
            }
            int num_u_ = mobod.getNumU(state);
            Real u_tmp_ = 0.0;
            if (num_u_ != 0)
            {
                for (int i = 0; i < num_u_; i++)
                {
                    std::string attr_name = "UIndx_" + std::to_string(i);
                    simbody_xml_engine_.getRequiredAttributeValue(ele_ite_, attr_name, u_tmp_);
                    mobod.setOneU(state, SimTK::UIndex(i), u_tmp_);
                }
            }
            num_mobod++;
        }
    }
}
//------------------------------------------------------------------------------
//          REALIZE THE SYSTEM TO THE REQUIRED COMPUTATIONAL STAGE
//------------------------------------------------------------------------------
//=============================================================================================//
void SimbodyStateEngine::realizeTime()
{
    getMultibodySystem().realize(working_state_, SimTK::Stage::Time);
}
//=============================================================================================//
void SimbodyStateEngine::realizePosition()
{
    getMultibodySystem().realize(working_state_, SimTK::Stage::Position);
}
//=============================================================================================//
void SimbodyStateEngine::realizeVelocity()
{
    getMultibodySystem().realize(working_state_, SimTK::Stage::Velocity);
}
//=============================================================================================//
void SimbodyStateEngine::realizeDynamics()
{
    getMultibodySystem().realize(working_state_, SimTK::Stage::Dynamics);
}
//=============================================================================================//
void SimbodyStateEngine::realizeAcceleration()
{
    getMultibodySystem().realize(working_state_, SimTK::Stage::Acceleration);
}
//=============================================================================================//
void SimbodyStateEngine::realizeReport()
{
    getMultibodySystem().realize(working_state_, SimTK::Stage::Report);
}
//=============================================================================================//
} // namespace SPH
