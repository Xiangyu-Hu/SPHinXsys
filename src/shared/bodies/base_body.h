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
#include "base_particle_generator.h"
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
    bool is_bound_set_;             /**< whether the bounding box is set */
    BoundingBox bound_;             /**< bounding box of the body */
    Shape *initial_shape_;          /**< initial volumetric geometry enclosing the body */

  public:
    SPHAdaptation *sph_adaptation_;        /**< numerical adaptation policy */
    BaseMaterial *base_material_;          /**< base material for dynamic cast in DataDelegate */
    StdVec<SPHRelation *> body_relations_; /**< all contact relations centered from this body **/

    SPHBody(SPHSystem &sph_system, Shape &shape, const std::string &name);
    SPHBody(SPHSystem &sph_system, Shape &shape);
    SPHBody(SPHSystem &sph_system, const std::string &name);
    SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr, const std::string &name);
    SPHBody(SPHSystem &sph_system, SharedPtr<Shape> shape_ptr);
    virtual ~SPHBody(){};

    std::string getName() { return body_name_; };
    SPHSystem &getSPHSystem();
    SPHBody &getSPHBody() { return *this; };
    Shape &getInitialShape() { return *initial_shape_; };
    void assignBaseParticles(BaseParticles *base_particles) { base_particles_ = base_particles; };
    BaseParticles &getBaseParticles();
    BaseMaterial &getBaseMaterial();
    StdVec<SPHRelation *> &getBodyRelations() { return body_relations_; };
    IndexRange LoopRange() { return IndexRange(0, base_particles_->TotalRealParticles()); };
    size_t SizeOfLoopRange() { return base_particles_->TotalRealParticles(); };
    Real getSPHBodyResolutionRef() { return sph_adaptation_->ReferenceSpacing(); };
    void setNewlyUpdated() { newly_updated_ = true; };
    void setNotNewlyUpdated() { newly_updated_ = false; };
    bool checkNewlyUpdated() { return newly_updated_; };
    void setSPHBodyBounds(const BoundingBox &bound);
    BoundingBox getSPHBodyBounds();
    BoundingBox getSPHSystemBounds();
    //----------------------------------------------------------------------
    //		Object factory template functions
    //----------------------------------------------------------------------
    virtual void defineAdaptationRatios(Real h_spacing_ratio, Real new_system_refinement_ratio = 1.0);

    template <class AdaptationType, typename... Args>
    void defineAdaptation(Args &&...args)
    {
        sph_adaptation_ = sph_adaptation_ptr_keeper_
                              .createPtr<AdaptationType>(sph_system_.ReferenceResolution(), std::forward<Args>(args)...);
    };

    template <typename... Args>
    LevelSetShape *defineComponentLevelSetShape(const std::string &shape_name, Args &&...args)
    {
        ComplexShape *complex_shape = DynamicCast<ComplexShape>(this, initial_shape_);
        return complex_shape->defineLevelSetShape(*this, shape_name, std::forward<Args>(args)...);
    };

    template <typename... Args>
    LevelSetShape *defineBodyLevelSetShape(Args &&...args)
    {
        LevelSetShape *level_set_shape =
            shape_ptr_keeper_.resetPtr<LevelSetShape>(*this, *initial_shape_, std::forward<Args>(args)...);

        initial_shape_ = level_set_shape;
        return level_set_shape;
    };

    template <class MaterialType>
    void assignMaterial(MaterialType *material)
    {
        base_material_ = material;
    };

    template <class MaterialType = BaseMaterial, typename... Args>
    MaterialType *defineMaterial(Args &&...args)
    {
        MaterialType *material = base_material_ptr_keeper_.createPtr<MaterialType>(std::forward<Args>(args)...);
        assignMaterial(material);
        return material;
    };
    //----------------------------------------------------------------------
    // Particle generating methods
    // Initialize particle data using a particle generator for geometric data.
    // The local material parameters are also initialized.
    //----------------------------------------------------------------------
    template <class ParticleType, class... Parameters, typename... Args>
    void generateParticles(Args &&...args)
    {
        ParticleType *particles = base_particles_ptr_keeper_.createPtr<ParticleType>(*this, base_material_);
        ParticleGenerator<ParticleType, Parameters...> particle_generator(*this, *particles, std::forward<Args>(args)...);
        particle_generator.generateParticlesWithGeometricVariables();
        particles->initializeBasicParticleVariables();
        sph_adaptation_->initializeAdaptationVariables(*particles);
        base_material_->setLocalParameters(sph_system_.ReloadParticles(), particles);
    };

    // Buffer or ghost particles can be generated together with real particles
    template <class ParticleType, typename... Parameters, class ReserveType, typename... Args>
    void generateParticlesWithReserve(ReserveType &particle_reserve, Args &&...args)
    {
        generateParticles<ParticleType, ReserveType, Parameters...>(particle_reserve, std::forward<Args>(args)...);
    };

    virtual void writeParticlesToVtuFile(std::ostream &output_file);
    virtual void writeParticlesToVtpFile(std::ofstream &output_file);
    virtual void writeParticlesToPltFile(std::ofstream &output_file);
    virtual void writeParticlesToXmlForRestart(std::string &filefullpath);
    virtual void readParticlesFromXmlForRestart(std::string &filefullpath);
    virtual void writeToXmlForReloadParticle(std::string &filefullpath);
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
    size_t iteration_count_;
    bool cell_linked_list_created_;

  public:
    template <typename... Args>
    RealBody(Args &&...args)
        : SPHBody(std::forward<Args>(args)...),
          iteration_count_(1), cell_linked_list_created_(false)
    {
        this->getSPHSystem().addRealBody(this);
    };
    virtual ~RealBody(){};
    BaseCellLinkedList &getCellLinkedList();
    void updateCellLinkedList();
    void updateCellLinkedListWithParticleSort(size_t particle_sort_period);
};
} // namespace SPH
#endif // BASE_BODY_H
