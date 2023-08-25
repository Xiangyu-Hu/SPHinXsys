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
 * @file 	my_memory_pool.h
 * @brief 	A class template for scalable memory allocation from memory blocks provided by an underlying allocator.
 * @details A memory_pool allocates and frees memory in a way that scales with the number of processors.
 *			The memory is obtained as big chunks from an underlying allocator specified by the template argument.
 *			The latter must satisfy the subset of the allocator requirements from the [allocator.requirements] ISO C++ Standard section.
 *			A memory_pool meet the Memory Pool named requirement.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef MY_MEMORY_POOL_H
#define MY_MEMORY_POOL_H

#define TBB_PREVIEW_MEMORY_POOL 1

#include "tbb/enumerable_thread_specific.h"
#include "tbb/memory_pool.h"

#include <list>

/**
 * @class MyMemoryPool
 * @brief Note that the data package T should has a default constructor.
 */
template <class T>
class MyMemoryPool
{

#ifdef __EMSCRIPTEN__
    std::allocator<T> my_pool; /**< memory pool. */
    std::list<T> data_list;    /**< list of all nodes allocated. */
#else
    tbb::memory_pool<std::allocator<T>> my_pool;            /**< memory pool. */
    typedef tbb::memory_pool_allocator<T> pool_allocator_t; /**< memory allocator. */
    std::list<T, pool_allocator_t> data_list;               /**< list of all nodes allocated. */
#endif
    std::list<T *> free_list; /**< list of all free nodes. */

  public:
#ifdef __EMSCRIPTEN__
    MyMemoryPool() : data_list(my_pool){};
#else
    MyMemoryPool() : data_list(pool_allocator_t(my_pool)){};
#endif

    ~MyMemoryPool(){
        // my_pool.recycle();
    };
    /**  Prepare an available node. */
    template <typename... Args>
    T *malloc(Args &&...args)
    {
        if (free_list.empty())
        {
            data_list.emplace_back(std::forward<Args>(args)...);
            return (&data_list.back());
        }
        else
        {
            T *result = free_list.front();
            free_list.pop_front();
            return result;
        }
    };
    /** Relinquish an unused node. */
    void free(T *ptr)
    {
        free_list.push_back(ptr);
    };
    /** Return the total number of nodes allocated. */
    int capacity()
    {
        return data_list.size();
    };
    /** Return the number of current available nodes. */
    int available_node()
    {
        return free_list.size();
    };
};

#endif // MY_MEMORY_POOL_H
