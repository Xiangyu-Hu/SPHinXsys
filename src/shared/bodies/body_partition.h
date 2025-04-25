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

#include "base_body_part.h"

namespace SPH
{
class BodyPartition : public BodyPartByID
{
  public:
    BodyPartition(SPHBody &sph_body, UnsignedInt adaptation_level);
    virtual ~BodyPartition() {};
    UnsignedInt AdaptationLevel() { return adaptation_level_; };
    virtual BaseCellLinkedList &getCellLinkedList() = 0;

  protected:
    UniquePtr<BaseCellLinkedList> cell_linked_list_ptr_;
    bool cell_linked_list_created_;
    UnsignedInt adaptation_level_;
    SPHAdaptation &sph_adaptation_;
};

class BodyPartitionTemporal : public BodyPartition
{
  public:
    BodyPartitionTemporal(SPHBody &sph_body, UnsignedInt adaptation_level);
    virtual ~BodyPartitionTemporal() {};
    virtual BaseCellLinkedList &getCellLinkedList() override;
};

class BodyPartitionSpatial : public BodyPartition
{
  public:
    BodyPartitionSpatial(SPHBody &sph_body, UnsignedInt adaptation_level);
    virtual ~BodyPartitionSpatial() {};
    virtual BaseCellLinkedList &getCellLinkedList() override;
};
} // namespace SPH
#endif // BODY_PARTITION_H
