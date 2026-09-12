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
 * @file base_particles.h
 * @brief This is the base class of SPH particles containing data
 * and operation for all types of particles.
 * Note that there is no class of single particle.
 * TODO: It seems that I need to transfer the IO related functions to the IO classes.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef BASE_PARTICLES_H
#define BASE_PARTICLES_H

#include "base_data_type_package.h"
#include "sphinxsys_variable.h"

namespace tinyxml2
{
class XMLElement; // forward declaration
}

namespace SPH
{
class SPHBody;
class SPHAdaptation;
class BodyPartByParticle;
class XmlParser;
using TinyXMLElement = tinyxml2::XMLElement;
class GroupManager;

/** Generalized particle data type */
typedef DataContainerAssemble<AllocatedData> ParticleData;
/**
 * @class SubdomainExchangeInterface
 * @brief What the shared dynamics and the I/O need from the decomposition of a body:
 *        the halo refresh of a set of variables, the assembly of the global particle
 *        set in the host arrays for output, and the owner of a position. Implemented
 *        by BodyDecomposition; declared here so that the shared code can use it
 *        without depending on the decomposition headers.
 */
class SubdomainExchangeInterface
{
  public:
    virtual ~SubdomainExchangeInterface() {};
    virtual void refreshHalo(DiscreteVariables &variables) = 0;
    /** Assemble the owned particles of all subdomains in the host arrays, in subdomain
     *  order, and publish the global particle count to the host thread. */
    virtual void gatherToHost(DiscreteVariables &variables) = 0;
    /** Undo the count publication of gatherToHost(); required before the next fan-out. */
    virtual void finishHostAccess() = 0;
    /** The subdomain owning a position. */
    virtual int subdomainOf(const Vecd &position) const = 0;
};

/**
 * @class BaseParticles
 * @brief Particles with essential (geometric and matter) data.
 * There are three groups of particles，all particles of a same type are saved with continuous memory segments.
 * The first is for real particles whose states are updated by particle dynamics.
 * One is reserved particles allocated after the real particles.
 * The global value of total_real_particles_ separate the real and reserved particles.
 * Reserved particles may be switched from real particles or switch to real particles.
 * As the memory for both particles are continuous, such switch is achieved at the memory boundary applying atomic operations.
 * The basic idea is swap the data of the last real particle with the one will be switched particle,
 * and then switch this swapped last particle as reserved particle by decrease the total_real_particles_ by one.
 * Switch from reserved particle to real particle is easy. One just need to assign expect state to
 * the first reserved particle and increase total_real_particles_ by one.
 * The third group is for ghost particles whose states are updated according to
 * boundary condition if their indices are included in the neighbor particle list.
 * The ghost particles are temporarily defined behind the real particles after the latter is determined.
 * All particles whose states are updated either by particle dynamics or boundary conditions
 * are separated from the reserved particles by total_particles_.
 * All real and reserved particles are bounded by particle_bound_, which is the limit of particle in all types.
 * In SPHinXsys, the discrete variables (state of each particle) registered belong to a hierarchy of two layers.
 * The first is for the basic geometric and matter properties.
 * These variables are created after particles are generated.
 * The second is for the local, dynamics-method-related variables,
 * which are defined by its data type and a unique name in specific methods,
 * and are only used by the relevant methods.
 */
class BaseParticles
{
  private:
    DataContainerUniquePtrAssemble<DiscreteVariable> all_discrete_variable_ptrs_;
    DataContainerUniquePtrAssemble<SingleVariable> all_singular_variable_ptrs_;
    UniquePtrsKeeper<Quantity> unique_variable_ptrs_;
    UniquePtrKeeper<XmlParser> xml_parser_ptr_;
    UniquePtrKeeper<GroupManager> group_manager_ptr_;

  public:
    explicit BaseParticles(SPHBody &sph_body);
    virtual ~BaseParticles();
    SPHBody &getSPHBody() { return sph_body_; };
    SPHAdaptation &getSPHAdaptation();
    std::string getBodyName();
    //----------------------------------------------------------------------
    // Global information for defining particle groups
    // total_real_particles_ gives the run-time total number of real particles.
    // particles_bound_ gives the total number of particles in all groups.
    //----------------------------------------------------------------------
  protected:
    SingleVariable<UnsignedInt> *sv_total_real_particles_;
    /** Real particles plus the halo copies received from neighboring subdomains.
     *  Equal to sv_total_real_particles_ outside a domain decomposed run. The two
     *  counts differ in what they are used for: physics is integrated over the real
     *  (owned) particles only, whereas the cell linked list and hence the neighbor
     *  search must also see the halo, or particles near a cut plane would lose part
     *  of their support. */
    SingleVariable<UnsignedInt> *sv_total_local_particles_;
    UnsignedInt particles_bound_;

  public:
    /** initialize basic variables after the particles generated by particle generator */
    virtual void initializeBasicDiscreteVariables();
    //----------------------------------------------------------------------
    // Generalized particle manipulation
    //----------------------------------------------------------------------
    SingleVariable<UnsignedInt> *svTotalRealParticles() { return sv_total_real_particles_; };
    UnsignedInt TotalRealParticles() { return sv_total_real_particles_->getValue(); };
    SingleVariable<UnsignedInt> *svTotalLocalParticles() { return sv_total_local_particles_; };
    UnsignedInt TotalLocalParticles() { return sv_total_local_particles_->getValue(); };
    UnsignedInt ParticlesBound() { return particles_bound_; };
    GroupManager &getParticleGroupManager();
    void initializeAllParticlesBounds(UnsignedInt total_real_particles);
    void initializeAllParticlesBoundsFromReload();
    void increaseParticlesBounds(UnsignedInt extra_size);
    void checkEnoughReserve();
    //----------------------------------------------------------------------
    // Parameterized management on particle variables and data
    //----------------------------------------------------------------------
  public:
    template <typename DataType>
    DiscreteVariable<DataType> *getVariableByName(const std::string &name);
    template <class DataType, typename... Args>
    DiscreteVariable<DataType> *addUniqueDiscreteVariable(const std::string &name, UnsignedInt size, Args &&...args);
    template <typename DataType, typename... Args>
    DiscreteVariable<DataType> *registerDiscreteVariable(const std::string &name, UnsignedInt size, Args &&...args);
    template <typename DataType, typename... Args>
    DiscreteVariable<DataType> *registerStateVariable(const std::string &name, Args &&...args);
    template <typename DataType>
    DiscreteVariable<DataType> *registerStateVariableFrom(const std::string &new_name, const std::string &old_name);
    template <typename DataType>
    DiscreteVariable<DataType> *registerStateVariableFrom(const std::string &name, const StdVec<DataType> &geometric_data);
    template <typename DataType>
    DiscreteVariable<DataType> *registerStateVariableFromReload(const std::string &name);
    template <typename DataType>
    SingleVariable<DataType> *addUniqueSingleVariable(const std::string &name, DataType initial_value = ZeroData<DataType>::value);
    template <typename DataType>
    SingleVariable<DataType> *registerSingleVariable(const std::string &name, DataType initial_value = ZeroData<DataType>::value);
    template <typename DataType>
    SingleVariable<DataType> *getSingleVariableByName(const std::string &name);
    //----------------------------------------------------------------------
    // Manage subsets of particle variables
    //----------------------------------------------------------------------
    template <typename DataType>
    DiscreteVariable<DataType> *addDiscreteVariableToList(DiscreteVariables &variable_set, const std::string &name);
    template <typename DataType>
    DiscreteVariable<DataType> *addDiscreteVariableToList(DiscreteVariables &variable_set, DiscreteVariable<DataType> *variable);
    //----------------------------------------------------------------------
    // Domain decomposed run
    //----------------------------------------------------------------------
    /** Installed by the BodyDecomposition of this body, if any. Interaction algorithms
     *  call refreshHalo() with the variables they read at the neighbors right before
     *  their interaction step, and the recorders bracket their host side access with
     *  gatherToHost() and finishHostAccess(), so that neither needs placing by hand in
     *  the case file. All of them are no-ops when the body is not decomposed. */
    void setSubdomainExchange(SubdomainExchangeInterface *subdomain_exchange) { subdomain_exchange_ = subdomain_exchange; };
    SubdomainExchangeInterface *getSubdomainExchange() { return subdomain_exchange_; };
    void refreshHalo(DiscreteVariables &variables)
    {
        if (subdomain_exchange_ != nullptr)
            subdomain_exchange_->refreshHalo(variables);
    };
    void gatherToHost(DiscreteVariables &variables)
    {
        if (subdomain_exchange_ != nullptr)
            subdomain_exchange_->gatherToHost(variables);
    };
    void finishHostAccess()
    {
        if (subdomain_exchange_ != nullptr)
            subdomain_exchange_->finishHostAccess();
    };
    //----------------------------------------------------------------------
    // Particle data for sorting
    //----------------------------------------------------------------------
  protected:
    SubdomainExchangeInterface *subdomain_exchange_ = nullptr; /**< see setSubdomainExchange() */
    UnsignedInt *original_id_;             /**< the original ids assigned just after particle is generated. */
    UnsignedInt *sorted_id_;               /**< the current sorted particle ids of particles from original ids. */
    DiscreteVariables evolving_variables_; // particle variables which evolving during simulation

  public:
    DiscreteVariables &VariablesToWrite() { return variables_to_write_; };
    DiscreteVariables &EvolvingVariables() { return evolving_variables_; };
    DiscreteVariables &ParticleAttributesToWrite() { return particle_attributes_to_write_; };
    StdVec<std::string> &ParticleGroupsToWrite() { return particle_groups_to_write_; };
    void addParticleGroupToWrite(const std::string &name);
    template <typename DataType, typename... Args>
    void addEvolvingVariable(Args &&...args);
    template <typename DataType, typename... Args>
    void addVariableToWrite(Args &&...args);
    //----------------------------------------------------------------------
    // Particle data ouput functions
    //----------------------------------------------------------------------
    void resizeXmlDocForParticles(XmlParser &xml_parser);
    void resetTotalRealParticlesFromXmlDoc(XmlParser &xml_parser);
    void writeParticlesToXml(XmlParser &xml_parser, TinyXMLElement *body_element);
    void readParticlesFromXml(XmlParser &xml_parser, TinyXMLElement *body_element);
    void readReloadXmlFile(const std::string &filefullpath, const std::string &body_name);
    template <typename DataType>
    BaseParticles &reloadExtraVariable(const std::string &name);
    //----------------------------------------------------------------------
    // Function related to geometric variables and their relations
    //----------------------------------------------------------------------
    void registerPositionAndVolumetricMeasure(StdVec<Vecd> &pos, StdVec<Real> &Vol);
    void registerPositionAndVolumetricMeasureFromReload();
    DiscreteVariable<Vecd> *dvParticlePosition() { return dv_pos_; }

  protected:
    DiscreteVariable<Vecd> *dv_pos_; /**< Discrete variable position */
    Real *Vol_;                      /**< Volumetric measure, also area and length of surface and linear particle */
    SPHBody &sph_body_;
    XmlParser &reload_xml_parser_;
    DiscreteVariables all_discrete_variables_;
    SingleVariables all_singular_variables_;
    DiscreteVariables variables_to_write_; // position is included
    DiscreteVariables particle_attributes_to_write_; // position is excluded
    StdVec<std::string> particle_groups_to_write_;

  protected:
    int total_body_parts_;                                /**< total number of body parts indicated particle groups*/
    StdVec<BodyPartByParticle *> body_parts_by_particle_; /**< all body parts by particle */

  public:
    int getNewBodyPartID();
    void addBodyPartByParticle(BodyPartByParticle *body_part) { body_parts_by_particle_.push_back(body_part); };
    StdVec<BodyPartByParticle *> getBodyPartsByParticle() { return body_parts_by_particle_; };

  protected:
    //----------------------------------------------------------------------
    // Small structs for generalize particle operations on
    // assembled variables and data sets
    //----------------------------------------------------------------------
    struct CopyParticleState
    {
        template <typename DataType>
        void operator()(DataContainerKeeper<AllocatedData<DataType>> &data_keeper, UnsignedInt index, UnsignedInt another_index);
    };

    struct WriteParticleVariableToXmlElement
    {
        TinyXMLElement *element_;
        WriteParticleVariableToXmlElement(TinyXMLElement *element) : element_(element) {}
        template <typename DataType>
        void operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, XmlParser &xml_parser);
    };

    struct ReadParticleVariableFromXmlElement
    {
        TinyXMLElement *element_;
        ReadParticleVariableFromXmlElement(TinyXMLElement *element) : element_(element) {}
        template <typename DataType>
        void operator()(DataContainerAddressKeeper<DiscreteVariable<DataType>> &variables, BaseParticles *base_particles, XmlParser &xml_parser);
    };

    OperationOnDataAssemble<ParticleData, CopyParticleState> copy_particle_state_;
    //----------------------------------------------------------------------
    // Functions for old CPU code compatibility
    //----------------------------------------------------------------------
  public:
    void copyFromAnotherParticle(UnsignedInt index, UnsignedInt another_index);
    UnsignedInt allocateGhostParticles(UnsignedInt ghost_size);
    void updateGhostParticle(UnsignedInt ghost_index, UnsignedInt index);
    void switchToBufferParticle(UnsignedInt index);
    UnsignedInt createRealParticleFrom(UnsignedInt index);

    template <typename DataType>
    DataType *getVariableDataByName(const std::string &name);
    template <class DataType, typename... Args>
    DataType *addUniqueDiscreteVariableData(const std::string &name, UnsignedInt size, Args &&...args);
    template <typename DataType, typename... Args>
    DataType *registerDiscreteVariableData(const std::string &name, UnsignedInt size, Args &&...args);
    template <typename DataType, typename... Args>
    DataType *registerStateVariableData(const std::string &name, Args &&...args);
    template <typename DataType, typename... Args>
    DataType *registerStateVariableDataFrom(const std::string &new_name, Args &&...args);
    template <typename DataType>
    DataType *registerStateVariableDataFromReload(const std::string &name);

    ParticleData evolving_variables_data_;
    ParticleData all_state_data_; /**< all discrete variable data except those on particle IDs  */
    ParticleData &EvolvingVariablesData() { return evolving_variables_data_; };

    Vecd *ParticlePositions() { return dv_pos_->Data(); }
    Real *VolumetricMeasures() { return Vol_; }
    virtual Real ParticleVolume(UnsignedInt index) { return Vol_[index]; }
    virtual Real ParticleSpacing(UnsignedInt index) { return std::pow(Vol_[index], 1.0 / Real(Dimensions)); }
    UnsignedInt *ParticleOriginalIds() { return original_id_; };
    UnsignedInt *ParticleSortedIds() { return sorted_id_; };

  protected:
    void incrementTotalRealParticles(UnsignedInt increment = 1) { sv_total_real_particles_->incrementValue(increment); };
    void decrementTotalRealParticles(UnsignedInt decrement = 1) { sv_total_real_particles_->incrementValue(-decrement); };
};
} // namespace SPH
#endif // BASE_PARTICLES_H
