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
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	state_engine.h
 * @details The SimbodyStateEngine class defines the interface used to add computational
 *          elements to the underlying SimTK::System (MultibodySystem). It specifies
 *          the interface that simbody states must satisfy in order to be part of the system
 *          and provides a series of helper methods for adding variables
 *          (state, discrete, cache, ...) to the underlying system. As such, SimbodyState
 *          handles all of the bookkeeping of system indices and provides convenience
 *          access to variable values (incl. derivatives) via their names as strings.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef STATE_ENGINE_SIMBODY_H
#define STATE_ENGINE_SIMBODY_H

#define _SILENCE_EXPERIMENTAL_FILESYSTEM_DEPRECATION_WARNING

#include "base_data_type_package.h"
#include "exception.h"

#include <cstdio>
#include <iostream>
#include <stdexcept>
#include <string>

#include "simbody_middle.h"
#include "xml_engine.h"
namespace SPH
{

class SPHSystem; // forward declaration

class SimbodyStateEngine
{
  protected:
    XmlEngine simbody_xml_engine_;
    void resizeXmlDocForSimbody(size_t input_size);

    /**
     * @class StateVariable
     * @details Derived simbody sate must create concrete StateVariables to expose their state
     *      variables. When exposing state variables allocated by the underlying Simbody
     *      (MobilizedBody, Constraint, Force, etc...) use its interface to
     *      implement the virtual methods below.
     */
    class StateVariable
    {
      public:
        /** Constructor and destructor. */
        StateVariable()
            : name_(""), owner_(nullptr), subsysindex_(SimTK::InvalidIndex), varindex_(SimTK::InvalidIndex), sysyindex_(SimTK::InvalidIndex)
        {
        }

        explicit StateVariable(std::string &name, SimbodyStateEngine &owner, SimTK::SubsystemIndex subsys, int varindex)
            : name_(name), owner_(&owner), subsysindex_(subsys), varindex_(varindex), sysyindex_(SimTK::InvalidIndex)
        {
        }

        virtual ~StateVariable() {}

        std::string &getName() { return name_; }
        SimbodyStateEngine &getOwner() { return *owner_; }
        /** Get the index of simbody state variable. */
        int &getVarIndex() { return varindex_; }
        /** Return the index of the subsystem used to make resource allocations. */
        SimTK::SubsystemIndex &getSubsysIndex() { return subsysindex_; }
        /** Return the index in the global list of continuous state variables, Y. */
        SimTK::SystemYIndex &getSystemYIndex() { return sysyindex_; }
        /** Set the index of simbody variable. */
        void setVarIndex(int index) { varindex_ = index; }
        /** Set the index of the subsystem used to resource allocations. */
        void setSubsystemIndex(SimTK::SubsystemIndex &subsysindx)
        {
            subsysindex_ = subsysindx;
        }
        /** Concrete StateEngines implement how the state variable value is evaluated. */
        virtual Real getValue() = 0;
        virtual void setValue(Real value) = 0;
        virtual Real getDerivative() = 0;
        /** The derivative a state should be a cache entry and thus does not change the state. */
        virtual void setDerivative(Real deriv) = 0;

      private:
        std::string name_;
        SimTK::ReferencePtr<SimbodyStateEngine> owner_;
        /**
         *  Identify which subsystem this state variable belongs to, which should
         *  be determined and set at creation time
         */
        SimTK::SubsystemIndex subsysindex_;
        /**
         *  The local variable index in the subsystem also provided at creation
         * (e.g. can be QIndex, UIndex, or Zindex type)
         */
        int varindex_;
        /**
         *  Once allocated a state in the system will have a global index
         *  and that can be stored here as well
         */
        SimTK::SystemYIndex sysyindex_;
    };
    /**
     * @class AddedStateVariable
     * @brief Class for handling state variable added (allocated) by this SimbodyStateEngine.
     */
    class AddedStateVariable : public StateVariable
    {
      public:
        /** Constructors adn destructors. */
        AddedStateVariable() : StateVariable(), invalidatestage_(SimTK::Stage::Empty) {}

        /** Convenience constructor for defining a SimbodyStateEngine added state variable */
        explicit AddedStateVariable(std::string &name, /**< state var name. */
                                    SimbodyStateEngine &owner,
                                    SimTK::Stage invalidatestage) : /**< stage this variable invalidates. */
                                                                    StateVariable(name, owner,
                                                                                  SimTK::SubsystemIndex(SimTK::InvalidIndex),
                                                                                  SimTK::InvalidIndex),
                                                                    invalidatestage_(SimTK::Stage::Empty)
        {
        }

        /** override virtual methods. */
        Real getValue() override;
        void setValue(Real value) override;
        Real getDerivative() override;
        void setDerivative(Real deriv) override;

      private:
        /** Changes in state variables trigger recalculation of appropriate cache
         * variables by automatically invalidating the realization stage specified
         * upon allocation of the state variable.
         */
        SimTK::Stage invalidatestage_;
    };

    UniquePtrsKeeper<AddedStateVariable> added_state_variable_keeper_;
    /**
     * @struct StateVariableInfo
     * @brief   To hold related info about discrete variables.
     */
    struct StateVariableInfo
    {
        StateVariableInfo() {}
        explicit StateVariableInfo(SimbodyStateEngine::StateVariable *sv, int order) : statevariable_(sv), order(order) {}

        /** Need empty copy constructor because default compiler generated
            will fail since it cannot copy a unique_ptr. */
        StateVariableInfo(const StateVariableInfo &) {}
        /** Now handle assignment by moving ownership of the unique pointer. */
        StateVariableInfo &operator=(const StateVariableInfo &svi)
        {
            if (this != &svi)
            {
                /** assignment has to be const but cannot swap const
                    want to keep unique pointer to guarantee no multiple reference
                    so use const_cast to swap under the covers. */
                StateVariableInfo *mutableSvi = const_cast<StateVariableInfo *>(&svi);
                statevariable_.swap(mutableSvi->statevariable_);
            }
            order = svi.order;
            return *this;
        }

        // State variable
        std::unique_ptr<SimbodyStateEngine::StateVariable> statevariable_;
        // order of allocation
        int order;
    };

  public:
    /** Add a continuous system state variable belonging to this Engine,
        and assign a name by which to refer to it. Changing the value of this state
        variable will automatically invalidate everything at and above its
        \a invalidatesStage, which is normally Stage::Dynamics meaning that there
        are forces that depend on this variable. If you define one or more
        of these variables you must also override computeStateVariableDerivatives()
        to provide time derivatives for them. Note, all corresponding system
        indices are automatically determined using this interface. As an advanced
        option you may choose to hide the state variable from being accessed outside
        of this component, in which case it is considered to be "hidden".
        You may also want to create an Output for this state variable; see
        #OpenSim_DECLARE_OUTPUT_FOR_STATE_VARIABLE for more information. Reporters
        should use such an Output to get the StateVariable's value (instead of using
        getStateVariableValue()).

        @param[in] stateVariableName     string value to access variable by name
        @param[in] invalidatesStage      the system realization stage that is
                                 invalidated when variable value is changed
        @param[in] isHidden              flag (bool) to optionally hide this state
                                 variable from being accessed outside this
                                 component as an Output
    */
    void addStateVariable(std::string statevariablename,
                          SimTK::Stage invalidatestage);

    /** The above method provides a convenient interface to this method, which
        automatically creates an 'AddedStateVariable' and allocates resources in the
        SimTK::State for this variable.  This interface allows the creator to
        add/expose state variables that are allocated by underlying Simbody
        components and specify how the state variable value is accessed by
        implementing a concrete StateVariable and adding it to the SimbodyStateEngine using
        this method.
        You may also want to create an Output for this state variable; see
        #OpenSim_DECLARE_OUTPUT_FOR_STATE_VARIABLE for more information. Reporters
        should use such an Output to get the StateVariable's value (instead of
        using getStateVariableValue()).
    */
    void addStateVariable(SimbodyStateEngine::StateVariable *statevariable);

    SimTK::DefaultSystemSubsystem &getDefaultSubsystem()
    {
        return const_cast<SimTK::DefaultSystemSubsystem &>(getMultibodySystem().getDefaultSubsystem());
    }
    SimTK::DefaultSystemSubsystem &updDefaultSubsystem()
    {
        return getMultibodySystem().updDefaultSubsystem();
    }

    /**
     * Get a StateVariable anywhere in the state engine, given a
     * StateVariable path.
     * This returns nullptr if a StateVariable does not exist at the specified
     * path or if the path is invalid.
     */
    StateVariable *traverseToStateVariable(std::string &pathname);
    /** Map names of continuous state variables of the Engine to their
        underlying SimTK indices. */
    mutable std::map<std::string, StateVariableInfo> namedstatevariableinfo_;
    /** Check that the list of _allStateVariables is valid. */
    bool isAllStatesVariablesListValid();

    /** Array of all state variables for fast access during simulation. */
    mutable SimTK::Array_<SimTK::ReferencePtr<StateVariable>>
        allstatevariables_;
    /** A handle the System associated with the above state variables. */
    mutable SimTK::ReferencePtr<SimTK::System> statesassociatedsystem_;

    /** Default constructor **/
    SimbodyStateEngine(SPHSystem &sph_system, SimTK::MultibodySystem &system);

    /** Reference pointer to the system that this engine manage. */
    SimTK::ReferencePtr<SimTK::MultibodySystem> mbsystem_;
    /** This is the internal 'writable' state of the engine.
     * working_state_ will be set to the system default state when
     * initializeState() is called.
     */
    SimTK::State working_state_;

    /** Destructor is virtual to allow concrete SimbodyStateEngine to cleanup. **/
    virtual ~SimbodyStateEngine() {};
    /** Set up the working state in present engine */
    void InitializeState();
    /**
     * Get the underlying MultibodySystem that this SimbodyStateEngine is connected to.
     * Make sure you have called Model::initSystem() prior to accessing the System.
     * Throws an Exception if the System has not been created or the SimbodyStateEngine
     * has not added itself to the System.
     * @see hasSystem().  */
    SimTK::MultibodySystem &getMultibodySystem();
    /**
     * Get writable reference to the MultibodySystem that this component is
     * connected to.
     */
    SimTK::MultibodySystem &updMultibodySystem();
    /**
     * Get the number of "continuous" state variables maintained by the
     * State Engine.
     */
    int getNumOfStateVariables();
    /** Get the number of continuous states that the State Engine added to the
        underlying computational system.*/
    int getNumStateVariablesAddedByEngine()
    {
        return (int)namedstatevariableinfo_.size();
    }
    /**
     * Get the names of "continuous" state variables maintained by the Engine
     */
    StdVec<std::string> getStateVariableNames();
    /**
     * Get all values of the state variables allocated by this SimbodyStateEngine.
     *
     * @param state   the State for which to get the value
     * @return Vector of state variable values of length getNumStateVariables()
     *                in the order returned by getStateVariableNames()
     * @throws StateEngineHasNoSystem if this object has not been added to a
     *         System (i.e., if initSystem has not been called)
     */
    SimTK::Vector getStateVariableValues();
    /**
     * report the state info by requested.
     */
    void reporter(SimTK::State &state_);
    /**
     * Write the state data to xml file.
     * For all bodies in the matter system, their generalized coordinates,
     * generalized velocities and transformations of the origin points are written in
     * the output file
     */
    std::string restart_folder_;
    void writeStateToXml(int ite_rst, SimTK::RungeKuttaMersonIntegrator &integ);
    /**
     * read state info from xml and set it to sate.
     * For all bodies in the matter system, their generalized coordinates,
     * generalized velocities and transformations of the origin points are read from
     * the restart file
     */
    void readStateFromXml(int ite_rst, SimTK::State &state);
    /**@name  Realize the Simbody System and State to Computational Stage.
            Methods in this section enable advanced and scripting users access to
            realize the Simbody MultibodySystem and the provided state to a desired
            computational (realization) Stage.
    */

    /**
     * Perform computations that depend only on time and earlier stages.
     */
    void realizeTime();
    /**
     * Perform computations that depend only on position-level state
     * variables and computations performed in earlier stages (including time).
     */
    void realizePosition();
    /**
     * Perform computations that depend only on velocity-level state
     * variables and computations performed in earlier stages (including position,
     * and time).
     */
    void realizeVelocity();
    /**
     * Perform computations (typically forces) that may depend on
     * dynamics-stage state variables, and on computations performed in earlier
     * stages (including velocity, position, and time), but not on other forces,
     * accelerations, constraint multipliers, or reaction forces.
     */
    void realizeDynamics();
    /**
     * Perform computations that may depend on applied forces.
     */
    void realizeAcceleration();
    /**
     * Perform computations that may depend on anything but are only used
     * for reporting and cannot affect subsequent simulation behavior.
     */
    void realizeReport();
};
} // namespace SPH
#endif // STATE_ENGINE_SIMBODY_H