/* -------------------------------------------------------------------------*
 *								SPHinXsys									*
 * -------------------------------------------------------------------------*
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle*
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for	*
 * physical accurate simulation and aims to model coupled industrial dynamic*
 * systems including fluid, solid, multi-body dynamics and beyond with SPH	*
 * (smoothed particle hydrodynamics), a meshless computational method using	*
 * particle discretization.													*
 *																			*
 * SPHinXsys is partially funded by German Research Foundation				*
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,			*
 *  HU1527/12-1 and HU1527/12-4												*
 *                                                                          *
 * Portions copyright (c) 2017-2022 Technical University of Munich and		*
 * the authors' affiliations.												*
 *                                                                          *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may  *
 * not use this file except in compliance with the License. You may obtain a*
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.       *
 *                                                                          *
 * ------------------------------------------------------------------------*/
/**
 * @file 	exceution_policy.h
 * @brief 	This is for the base functions for particle iterator.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef EXECUTION_POLICY_H
#define EXECUTION_POLICY_H

namespace SPH
{
    namespace execution
    {
        class SequencedPolicy
        {
        };

        class ParallelPolicy
        {
        };

        class ParallelUnsequencedPolicy
        {
        };

        inline constexpr SequencedPolicy seq{};
        inline constexpr ParallelPolicy par{};
        inline constexpr ParallelUnsequencedPolicy par_unseq{};
    }
}
#endif // EXECUTION_POLICY_H
