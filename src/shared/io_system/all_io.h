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
 * @file all_io.h
 * @brief This is the header file that user code should include to pick up all
 * IO class used in SPHinXsys.
 * @author Xiangyu Hu
 */

#ifndef ALL_IO_H
#define ALL_IO_H

#include "io_base.h"
#include "io_base_ck.h"
#include "io_log.h"
#include "io_observation.h"
#include "io_observation_ck.h"
#include "io_plt.hpp"
#include "io_simbody.h"
#include "io_vtk.h"
#include "io_vtk_mesh.h"
#include "io_vtk_mesh_ck.h"

#endif // ALL_IO_H
