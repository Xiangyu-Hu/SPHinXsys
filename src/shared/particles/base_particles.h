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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	base_particles.h
 * @brief 	This is the base class of SPH particles. The basic data of the particles
 *			is saved in separated large vectors. Each derived class will introduce several extra
 * 			vectors for the new data. Note that there is no class of single particle.
 * @author	Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef BASE_PARTICLES_H
#define BASE_PARTICLES_H

#include "base_data_package.h"
#include "base_material.h"
#include "base_variable.h"
#include "particle_sorting.h"
#include "sph_data_containers.h"
#include "xml_engine.h"

#include <fstream>

namespace SPH
{

class SPHBody;
class BaseMaterial;
class ParticleGenerator;
class BodySurface;
template <class ReturnType>
class BaseDynamics;

/**
 * @class BaseParticles
 * @brief Particles with essential (geometric and kinematic) data.
 * 		  There are three types of particles， all particles of a same type are saved with continuous memory segments.
 * 		  The first type is real particles whose states are updated by particle dynamics.
 * 		  One is buffer particles whose state are not updated by particle dynamics.
 * 		  Buffer particles are saved behind real particles.
 * 		  The global value of total_real_particles_ separate the real and buffer particles.
 * 		  They may be switched from real particles or switch to real particles.
 * 		  As the memory for both particles are continuous, such switch is achieved at the memory boundary sequentially.
 * 		  The basic idea is swap the data of the last real particle with the one will be switched particle,
 * 		  and then switch this swapped last particle as buffer particle by decrease the total_real_particles_ by one.
 * 		  Switch from buffer particle to real particle is easy. One just need to assign expect state to
 * 		  the first buffer particle and increase total_real_particles_ by one.
 * 		  The other is ghost particles whose states are updated according to
 * 		  boundary condition if their indices are included in the neighbor particle list.
 * 		  The ghost particles are saved behind the buffer particles.
 * 		  The global value of real_particles_bound_ separate the sum of real and buffer particles with ghost particles.
 * 		  The global value of total_ghost_particles_ indicates the total number of ghost particles in use.
 * 		  It will be initialized to zero before a time step.
 * 		  In SPHinXsys, the discrete variables (state of each particle) registered in general particle data
 *      (ParticleData) belong to a hierarchy of two layers.
 * 		  The first is for the basic physical states to describe the physical process.
 * 		  These variables are defined within the classes of particles.
 * 		  The second is for the local, dynamics-method-related variables, which are defined in specific methods,
 * 		  and are only used by the relevant methods. Generally, a discrete variable is defined
 *      and the corresponding data owned by one object so that other objects can use it by the function
 *      getVariableByName. A shared discrete variable can also be defined by several objects.
 *      In this case, the data is owned by BaseParticles within all_shared_data_ptrs_.
 */
class BaseParticles
{
  private:
    DataContainerUniquePtrAssemble<DiscreteVariable> all_discrete_variable_ptrs_;
    DataContainerUniquePtrAssemble<StdLargeVec> shared_particle_data_ptrs_;
    DataContainerUniquePtrAssemble<GlobalVariable> all_global_variable_ptrs_;
    UniquePtrsKeeper<BaseDynamics<void>> derived_particle_data_;

  public:
    explicit BaseParticles(SPHBody &sph_body, BaseMaterial *base_material);
    virtual ~BaseParticles(){};

    StdLargeVec<Vecd> pos_;       /**< Position */
    StdLargeVec<Vecd> vel_;       /**< Velocity */
    StdLargeVec<Vecd> acc_;       /**< Acceleration induced by pressure- or stress */
    StdLargeVec<Vecd> acc_prior_; /**< Other, such as gravity and viscous, accelerations computed before acc_ */

    StdLargeVec<Real> Vol_;              /**< Volumetric measure, also area and length of surface and linear particle */
    StdLargeVec<Real> rho_;              /**< Density */
    StdLargeVec<Real> mass_;             /**< Massive measure, also mass per-unit thickness and per-unit cross-section of surface and linear particle */
    StdLargeVec<int> surface_indicator_; /**< free surface indicator */
    //----------------------------------------------------------------------
    // Global information for defining particle groups
    //----------------------------------------------------------------------
    size_t total_real_particles_;
    size_t real_particles_bound_; /**< Maximum possible number of real particles. Also start index of ghost particles. */
    size_t total_ghost_particles_;

    SPHBody &getSPHBody() { return sph_body_; };
    BaseMaterial &getBaseMaterial() { return base_material_; };
    ParticleData &getAllParticleData() { return all_particle_data_; };
    /** initialize other variables after the particles are generated */
    virtual void initializeOtherVariables();
    //----------------------------------------------------------------------
    //		Generalized particle manipulation
    //----------------------------------------------------------------------
    void addBufferParticles(size_t buffer_size);
    void copyFromAnotherParticle(size_t index, size_t another_index);
    void updateFromAnotherParticle(size_t index, size_t another_index);
    size_t insertAGhostParticle(size_t index);
    void switchToBufferParticle(size_t index);
    //----------------------------------------------------------------------
    //		Parameterized management on generalized particle data
    //----------------------------------------------------------------------
    template <typename DataType>
    void registerVariable(StdLargeVec<DataType> &variable_addrs, const std::string &variable_name,
                          DataType initial_value = ZeroData<DataType>::value);
    template <typename DataType, class InitializationFunction>
    void registerVariable(StdLargeVec<DataType> &variable_addrs, const std::string &variable_name,
                          const InitializationFunction &initialization);
    template <typename DataType>
    StdLargeVec<DataType> *registerSharedVariable(
        const std::string &variable_name, const DataType &default_value = ZeroData<DataType>::value);
    template <typename DataType>
    StdLargeVec<DataType> *getVariableByName(const std::string &variable_name);
    ParticleVariables &AllDiscreteVariables() { return all_discrete_variables_; };

    template <typename DataType>
    DataType *registerGlobalVariable(const std::string &variable_name,
                                     DataType initial_value = ZeroData<DataType>::value);
    template <typename DataType>
    DataType *getGlobalVariableByName(const std::string &variable_name);
    //----------------------------------------------------------------------
    //		Manage subsets of particle variables
    //----------------------------------------------------------------------
    template <typename DataType>
    void addVariableToList(ParticleVariables &variable_set, const std::string &variable_name);
    template <typename DataType>
    void addVariableToWrite(const std::string &variable_name);
    template <typename DataType>
    void addVariableToRestart(const std::string &variable_name);
    inline const ParticleVariables &getVariablesToRestart() const { return variables_to_restart_; }
    template <typename DataType>
    void addVariableToReload(const std::string &variable_name);
    inline const ParticleVariables &getVariablesToReload() const { return variables_to_reload_; }

    template <class DerivedVariableMethod, class... Ts>
    void addDerivedVariableToWrite(Ts &&...);
    //----------------------------------------------------------------------
    //		Particle data for sorting
    //----------------------------------------------------------------------
    StdLargeVec<size_t> unsorted_id_; /**< the ids assigned just after particle generated. */
    StdLargeVec<size_t> sorted_id_;   /**< the sorted particle ids of particles from unsorted ids. */
    StdLargeVec<size_t> sequence_;    /**< the sequence referred for sorting. */
    ParticleData sortable_data_;
    ParticleVariables sortable_variables_;
    ParticleSorting particle_sorting_;

    template <typename DataType>
    void registerSortableVariable(const std::string &variable_name);
    template <typename SequenceMethod>
    void sortParticles(SequenceMethod &sequence_method);
    //----------------------------------------------------------------------
    //		Particle data ouput functions
    //----------------------------------------------------------------------
    template <typename OutStreamType>
    void writeParticlesToVtk(OutStreamType &output_stream);
    void writeParticlesToPltFile(std::ofstream &output_file);
    virtual void writeSurfaceParticlesToVtuFile(std::ostream &output_file, BodySurface &surface_particles);
    void resizeXmlDocForParticles(XmlEngine &xml_engine);
    void writeParticlesToXmlForRestart(std::string &filefullpath);
    void readParticleFromXmlForRestart(std::string &filefullpath);
    void writeToXmlForReloadParticle(std::string &filefullpath);
    void readFromXmlForReloadParticle(std::string &filefullpath);
    XmlEngine *getReloadXmlEngine() { return &reload_xml_engine_; };
    virtual BaseParticles *ThisObjectPtr() { return this; };
    //----------------------------------------------------------------------
    //		Relation relate volume, surface and linear particles
    //----------------------------------------------------------------------
    virtual Real ParticleVolume(size_t index) { return Vol_[index]; }
    virtual Real ParticleMass(size_t index) { return mass_[index]; }

  protected:
    SPHBody &sph_body_;
    std::string body_name_;
    BaseMaterial &base_material_;
    XmlEngine restart_xml_engine_;
    XmlEngine reload_xml_engine_;
    ParticleData all_particle_data_;
    ParticleVariables all_discrete_variables_;
    GlobalVariables all_global_variables_;
    ParticleVariables variables_to_write_;
    ParticleVariables variables_to_restart_;
    ParticleVariables variables_to_reload_;
    StdVec<BaseDynamics<void> *> derived_variables_;

    void addAParticleEntry(); /**< Add a particle entry to the particle array. */
    virtual void writePltFileHeader(std::ofstream &output_file);
    virtual void writePltFileParticleData(std::ofstream &output_file, size_t index);
    //----------------------------------------------------------------------
    //		Small structs for generalize particle operations
    //----------------------------------------------------------------------
    template <typename DataType>
    struct resizeParticleData
    {
        void operator()(ParticleData &particle_data, size_t new_size) const;
    };

    /** Add a particle data with default value. */
    template <typename DataType>
    struct addParticleDataWithDefaultValue
    {
        void operator()(ParticleData &particle_data) const;
    };

    template <typename DataType>
    struct copyParticleData
    {
        void operator()(ParticleData &particle_data, size_t index, size_t another_index) const;
    };
    //----------------------------------------------------------------------
    //		Assemble based generalize particle operations
    //----------------------------------------------------------------------
    DataAssembleOperation<resizeParticleData> resize_particle_data_;
    DataAssembleOperation<addParticleDataWithDefaultValue> add_particle_data_with_default_value_;
    DataAssembleOperation<copyParticleData> copy_particle_data_;
};

/**
 * @struct WriteAParticleVariableToXml
 * @brief Define a operator for writing particle variable to XML format.
 */
struct WriteAParticleVariableToXml
{
    XmlEngine &xml_engine_;
    size_t &total_real_particles_;
    WriteAParticleVariableToXml(XmlEngine &xml_engine, size_t &total_real_particles)
        : xml_engine_(xml_engine), total_real_particles_(total_real_particles){};

    template <typename DataType>
    void operator()(const std::string &variable_name, StdLargeVec<DataType> &variable) const;
};

/**
 * @struct ReadAParticleVariableFromXml
 * @brief Define a operator for reading particle variable to XML format.
 */
struct ReadAParticleVariableFromXml
{
    XmlEngine &xml_engine_;
    size_t &total_real_particles_;
    ReadAParticleVariableFromXml(XmlEngine &xml_engine, size_t &total_real_particles)
        : xml_engine_(xml_engine), total_real_particles_(total_real_particles){};

    template <typename DataType>
    void operator()(const std::string &variable_name, StdLargeVec<DataType> &variable) const;
};

/**
 * @class BaseDerivedVariable
 * @brief computing displacement from current and initial particle position
 */
template <typename DataType>
class BaseDerivedVariable
{
  public:
    using DerivedDataType = DataType;
    std::string variable_name_;

    BaseDerivedVariable(SPHBody &sph_body, const std::string &variable_name);
    virtual ~BaseDerivedVariable(){};

  protected:
    StdLargeVec<DataType> derived_variable_;
};
} // namespace SPH
#endif // BASE_PARTICLES_H
