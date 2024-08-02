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
 * @file    update_cell_linked_list_sycl.h
 * @brief   TBD
 * @author	Xiangyu Hu
 */

#ifndef UPDATE_CELL_LINKED_LIST_SYCL_H
#define UPDATE_CELL_LINKED_LIST_SYCL_H

#include "update_cell_linked_list.h"

namespace SPH
{
template <typename MeshType>
class UpdateCellLinkedList<MeshType, ParallelDevicePolicy>
    : public LocalDynamics, public BaseDynamics<void>
{
  public:
    UpdateCellLinkedList(SPHBody &sph_body);
    virtual ~UpdateCellLinkedList(){};
    virtual void exec(Real dt = 0.0) override;

    class ComputingKernel
    {
      public:
      protected:
        friend class UpdateCellLinkedList<MeshType, ParallelDevicePolicy>;
        MeshType mesh_;
    };

  protected:
    MeshType mesh_;
    Vecd *pos_;
    UnsignedInt *particle_id_list_;
    UnsignedInt *particle_offset_list_;
};
} // namespace SPH
#endif // UPDATE_CELL_LINKED_LIST_SYCL_H