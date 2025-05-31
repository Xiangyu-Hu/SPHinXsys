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
 * @file 	body_partition.h
 * @brief 	A part of the body for adaptation.
 * @details	TBD.
 * @author	Xiangyu Hu
 */

#ifndef BODY_PARTITION_H
#define BODY_PARTITION_H

#include "base_body.h"

namespace SPH
{
class BodyPartition
{
  public:
    typedef BodyPartition BaseIdentifier;
    BodyPartition(SPHBody &sph_body, UnsignedInt present_adapt_level_);
    virtual ~BodyPartition() {};
    SPHBody &getSPHBody() { return sph_body_; };
    SPHSystem &getSPHSystem() { return sph_body_.getSPHSystem(); };
    std::string getName();
    SPHAdaptation &getSPHAdaptation() { return sph_adaptation_; };
    BaseParticles &getBaseParticles() { return base_particles_; };
    UnsignedInt PresentAdaptationLevel() { return present_adapt_level_; };
    DiscreteVariable<int> *dvAdaptationLevel() { return dv_adapt_level_; };
    virtual BaseCellLinkedList &getCellLinkedList() = 0;

    class SourceParticleMask
    {
      public:
        template <class ExecutionPolicy, typename EnclosureType>
        SourceParticleMask(ExecutionPolicy &ex_policy, EnclosureType &encloser)
            : present_adapt_level_(encloser.present_adapt_level_),
              adapt_level_(encloser.dv_adapt_level_->DelegatedData(ex_policy)) {}
        ~SourceParticleMask() {}

        bool operator()(UnsignedInt source_index)
        {
            return adapt_level_[source_index] == present_adapt_level_;
        }

      protected:
        int present_adapt_level_;
        int *adapt_level_;
    };
    using ListedParticleMask = SourceParticleMask;

    template <typename TargetCriterion>
    class TargetParticleMask : public TargetCriterion
    {
      public:
        template <class ExecutionPolicy, typename EnclosureType, typename... Args>
        TargetParticleMask(ExecutionPolicy &ex_policy, EnclosureType &encloser, Args &&...args)
            : TargetCriterion(std::forward<Args>(args)...) {}
        virtual ~TargetParticleMask() {}
    };

  protected:
    SPHBody &sph_body_;
    SPHAdaptation &sph_adaptation_;
    BaseParticles &base_particles_;
    UniquePtr<BaseCellLinkedList> cell_linked_list_ptr_;
    bool cell_linked_list_created_;
    UnsignedInt present_adapt_level_;
    DiscreteVariable<int> *dv_adapt_level_;
};

class BodyPartitionTemporal : public BodyPartition
{
  public:
    BodyPartitionTemporal(SPHBody &sph_body, UnsignedInt adapt_level);
    virtual ~BodyPartitionTemporal() {};
    virtual BaseCellLinkedList &getCellLinkedList() override;
};

class BodyPartitionSpatial : public BodyPartition
{
  public:
    BodyPartitionSpatial(SPHBody &sph_body, UnsignedInt adapt_level);
    virtual ~BodyPartitionSpatial() {};
    Real getReferenceSmoothingLength() { return sph_adaptation_.SmoothingLengthByLevel(present_adapt_level_); };
    virtual BaseCellLinkedList &getCellLinkedList() override;
};
} // namespace SPH
#endif // BODY_PARTITION_H
