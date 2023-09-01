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
 * @file 	base_body.h
 * @brief 	These are the base classes of SPH bodies. The RealBody is for
 *			that with cell linked list and the SPHBody does not.
 * 			Before the definition of the SPH bodies, the shapes with complex
 *			geometries, i.e. those are produced by advanced binary operations,
 * 			such as intersection, should be produced first.
 * 			Then, all shapes used in body definition should be either contain
 * 			or not contain each other. Partial overlap between them are not permitted.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef BASE_BODY_H
#define BASE_BODY_H

#include "adaptation.h"
#include "all_geometries.h"
#include "base_data_package.h"
#include "base_material.h"
#include "base_particles.h"
#include "cell_linked_list.h"
#include "particle_sorting.h"
#include "sph_data_containers.h"
#include "sph_system.h"

#include <string>

namespace SPH
{
class SPHRelation;
class BodySurface;

/**
 * @class SPHBody
 * @brief SPHBody is a base body with basic data and functions.
 *		  Its derived class can be a real fluid body, a real deformable solid body,
 *        a static or moving solid body or an observer body.
 * 		  Note that only real bodies have cell linked list.
 */
class SPHBody
{
  private:
    SharedPtrKeeper<Shape> shape_ptr_keeper_;
    UniquePtrKeeper<SPHAdaptation> sph_adaptation_ptr_keeper_;
    UniquePtrKeeper<BaseParticles> base_particles_ptr_keeper_;
    UniquePtrKeeper<BaseMaterial> base_material_ptr_keeper_;

  protected:
    SPHSystem &sph_system_;
    std::string body_name_;
    bool newly_updated_;            /**< whether this body is in a newly updated state */
    BaseParticles *base_particles_; /**< Base particles for dynamic cast DataDelegate  */

  public:
    Shape *body_shape_;                    /**< volumetric geometry enclosing the body */
    SPHAdaptation *sph_adaptation_;        /**< numerical adaptation policy */
    BaseMaterial *base_material_;          /**< base material for dynamic cast in DataDelegate */
    StdVec<SPHRelation *> body_relations_; /**< all contact relations centered from this body **/

    SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr, const std::string &body_name);
    SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr);
    virtual ~SPHBody(){};

    std::string getName() { return body_name_; };
    SPHSystem &getSPHSystem();
    SPHBody &getSPHBody() { return *this; };
    BaseParticles &getBaseParticles() { return *base_particles_; };
    size_t &LoopRange() { return base_particles_->total_real_particles_; };
    size_t SizeOfLoopRange() { return base_particles_->total_real_particles_; };
    Real getSPHBodyResolutionRef() { return sph_adaptation_->ReferenceSpacing(); };
    void setNewlyUpdated() { newly_updated_ = true; };
    void setNotNewlyUpdated() { newly_updated_ = false; };
    bool checkNewlyUpdated() { return newly_updated_; };
    BoundingBox getBodyShapeBounds();
    BoundingBox getSPHSystemBounds();
    void allocateConfigurationMemoriesForBufferParticles();
    //----------------------------------------------------------------------
    //		Object factory template functions
    //----------------------------------------------------------------------
    virtual void defineAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio = 1.0);

    template <class AdaptationType, typename... ConstructorArgs>
    void defineAdaptation(ConstructorArgs &&...args)
    {
        sph_adaptation_ = sph_adaptation_ptr_keeper_
                              .createPtr<AdaptationType>(*this, std::forward<ConstructorArgs>(args)...);
    };

    template <typename... ConstructorArgs>
    LevelSetShape *defineComponentLevelSetShape(const std::string &shape_name, ConstructorArgs &&...args)
    {
        ComplexShape *complex_shape = DynamicCast<ComplexShape>(this, body_shape_);
        return complex_shape->defineLevelSetShape(*this, shape_name, std::forward<ConstructorArgs>(args)...);
    };

    template <typename... ConstructorArgs>
    LevelSetShape *defineBodyLevelSetShape(ConstructorArgs &&...args)
    {
        LevelSetShape *levelset_shape =
            shape_ptr_keeper_.resetPtr<LevelSetShape>(*this, *body_shape_, std::forward<ConstructorArgs>(args)...);

        body_shape_ = levelset_shape;
        return levelset_shape;
    };

    /** partial construct particles with an already constructed material */
    template <class ParticleType = BaseParticles, class MaterialType = BaseMaterial>
    void defineParticlesWithMaterial(MaterialType *material)
    {
        base_material_ = material;
        base_particles_ = base_particles_ptr_keeper_.createPtr<ParticleType>(*this, material);
    };

    /** partial construct particles with material informaiton. note that particle data not initialized yet */
    template <class ParticleType = BaseParticles, class MaterialType = BaseMaterial, typename... ConstructorArgs>
    void defineParticlesAndMaterial(ConstructorArgs &&...args)
    {
        MaterialType *material = base_material_ptr_keeper_.createPtr<MaterialType>(std::forward<ConstructorArgs>(args)...);
        defineParticlesWithMaterial<ParticleType>(material);
    };

    /** initialize particle data using a particle generator for geometric data.
     * the local material parameters are also initialized. */
    template <class ParticleGeneratorType, typename... ConstructorArgs>
    void generateParticles(ConstructorArgs &&...args)
    {
        sph_adaptation_->registerAdaptationVariables(*base_particles_);
        ParticleGeneratorType particle_generator(*this, std::forward<ConstructorArgs>(args)...);
        particle_generator.generateParticlesWithBasicVariables();
        base_particles_->initializeOtherVariables();
        base_material_->initializeLocalParameters(base_particles_);
    };

    template <typename DataType>
    void addBodyState(StdLargeVec<DataType> &variable_addrs, const std::string &variable_name)
    {
        base_particles_->template registerVariable<DataType>(variable_addrs, variable_name);
    };

    template <typename VariableType>
    void addBodyStateForRecording(const std::string &variable_name)
    {
        base_particles_->template addVariableToWrite<VariableType>(variable_name);
    };

    template <class DerivedVariableMethod, typename... Args>
    void addDerivedBodyStateForRecording(Args &&...args)
    {
        base_particles_->template addDerivedVariableToWrite<DerivedVariableMethod>(std::forward<Args>(args)...);
    };

    template <typename VariableType>
    void addBodyStateToRestart(const std::string &variable_name)
    {
        base_particles_->template addVariableToRestart<VariableType>(variable_name);
    };

    virtual void writeParticlesToVtuFile(std::ostream &output_file);
    virtual void writeParticlesToVtpFile(std::ofstream &output_file);
    virtual void writeParticlesToPltFile(std::ofstream &output_file);
    virtual void writeSurfaceParticlesToVtuFile(std::ofstream &output_file, BodySurface &surface_particles);
    virtual void writeParticlesToXmlForRestart(std::string &filefullpath);
    virtual void readParticlesFromXmlForRestart(std::string &filefullpath);
    virtual void writeToXmlForReloadParticle(std::string &filefullpath);
    virtual void readFromXmlForReloadParticle(std::string &filefullpath);
    virtual SPHBody *ThisObjectPtr() { return this; };
};

/**
 * @class RealBody
 * @brief Derived body with inner particle configuration or inner interactions.
 * After construction, the particle and material must be specified.
 */
class RealBody : public SPHBody
{
  private:
    UniquePtr<BaseCellLinkedList> cell_linked_list_ptr_;
    /**
     * @brief particle by cells lists is for parallel splitting algorithm.
     * All particles in each cell are collected together.
     * If two particles each belongs two different cell entries,
     * they have no interaction because they are too far.
     */
    SplitCellLists split_cell_lists_;
    bool use_split_cell_lists_;
    size_t iteration_count_;
    bool cell_linked_list_created_;

  public:
    template <typename... ConstructorArgs>
    RealBody(ConstructorArgs &&...args)
        : SPHBody(std::forward<ConstructorArgs>(args)...),
          use_split_cell_lists_(false), iteration_count_(1),
          cell_linked_list_created_(false)
    {
        this->getSPHSystem().real_bodies_.push_back(this);
        size_t number_of_split_cell_lists = pow(3, Dimensions);
        split_cell_lists_.resize(number_of_split_cell_lists);
    };
    virtual ~RealBody(){};
    BaseCellLinkedList &getCellLinkedList();
    void setUseSplitCellLists() { use_split_cell_lists_ = true; };
    bool getUseSplitCellLists() { return use_split_cell_lists_; };
    SplitCellLists &getSplitCellLists() { return split_cell_lists_; };
    void updateCellLinkedList();
    void updateCellLinkedListWithParticleSort(size_t particle_sort_period);
};
} // namespace SPH
#endif // BASE_BODY_H
