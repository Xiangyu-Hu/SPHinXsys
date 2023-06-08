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
 *  HU1527/12-1 and HU1527/12-4													*
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
 * @file 	particle_iterators.h
 * @brief 	This is for the base functions for particle iterator.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef PARTICLE_ITERATORS_H
#define PARTICLE_ITERATORS_H

#include "base_data_package.h"
#include "sph_data_containers.h"
#include "execution_policy.h"
#include "execution_unit/execution_proxy.hpp"
#include "execution_unit/execution_queue.hpp"

namespace SPH
{
	using namespace execution;

	template <class ExecutionPolicy, typename DynamicsRange, class LocalDynamicsFunction>
	void particle_for(const ExecutionPolicy &execution_policy, const DynamicsRange &dynamics_range,
					  const LocalDynamicsFunction &local_dynamics_function)
	{
		std::cout << "\n Error: ExecutionPolicy, DynamicsRange or LocalDynamicsFunction not defined for particle dynamics !" << std::endl;
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		exit(1);
	};

	/**
	 * Body-wise iterators (for sequential and parallel computing).
	 */

	template <class LocalDynamicsFunction>
	inline void particle_for(const SequencedPolicy &seq, const size_t &all_real_particles,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i < all_real_particles; ++i)
			local_dynamics_function(i);
	};

	template <class LocalDynamicsFunction>
	inline void particle_for(const ParallelPolicy &par, const size_t &all_real_particles,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		parallel_for(
			IndexRange(0, all_real_particles),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					local_dynamics_function(i);
				}
			},
			ap);
	};

    template <class LocalDynamicsFunction, class Proxy>
    inline void particle_for(const ParallelPolicy &par_policy, const size_t &all_real_particles,
                             const LocalDynamicsFunction &local_dynamics_function, Proxy& proxy)
    {
        auto& kernel = *proxy.get(par_policy);
        parallel_for(
                IndexRange(0, all_real_particles),
                [&](const IndexRange &r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        local_dynamics_function(i, kernel);
                    }
                },
                ap);
    }

	/**
	 * Bodypart By Particle-wise iterators (for sequential and parallel computing).
	 */
	template <class LocalDynamicsFunction>
	inline void particle_for(const SequencedPolicy &seq, const IndexVector &body_part_particles,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i < body_part_particles.size(); ++i)
			local_dynamics_function(body_part_particles[i]);
	};

	template <class LocalDynamicsFunction>
	inline void particle_for(const ParallelPolicy &par, const IndexVector &body_part_particles,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		parallel_for(
			IndexRange(0, body_part_particles.size()),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					local_dynamics_function(body_part_particles[i]);
				}
			},
			ap);
	};
	/**
	 * Bodypart By Cell-wise iterators (for sequential and parallel computing).
	 */
	template <class LocalDynamicsFunction>
	inline void particle_for(const SequencedPolicy &seq, const ConcurrentCellLists &body_part_cells,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i != body_part_cells.size(); ++i)
		{
			ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
			for (size_t num = 0; num < particle_indexes.size(); ++num)
			{
				local_dynamics_function(particle_indexes[num]);
			}
		}
	}

	template <class LocalDynamicsFunction>
	inline void particle_for(const ParallelPolicy &par, const ConcurrentCellLists &body_part_cells,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		parallel_for(
			IndexRange(0, body_part_cells.size()),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
					for (size_t num = 0; num < particle_indexes.size(); ++num)
					{
						local_dynamics_function(particle_indexes[num]);
					}
				}
			},
			ap);
	};

    template <class LocalDynamicsFunction, class Proxy>
    inline void particle_for(const ParallelPolicy &par_policy, const ConcurrentCellLists &body_part_cells,
                             const LocalDynamicsFunction &local_dynamics_function, Proxy& proxy)
    {
        auto& kernel = *proxy.get(par_policy);
        parallel_for(
                IndexRange(0, body_part_cells.size()),
                [&](const IndexRange &r)
                {
                    for (size_t i = r.begin(); i < r.end(); ++i)
                    {
                        ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
                        for (size_t num = 0; num < particle_indexes.size(); ++num)
                        {
                            local_dynamics_function(particle_indexes[num], kernel);
                        }
                    }
                },
                ap);
    }

	/**
	 * BodypartByCell-wise iterators on cells (for sequential and parallel computing).
	 */
	template <class LocalDynamicsFunction>
	inline void particle_for(const SequencedPolicy &seq, const DataListsInCells &body_part_cells,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i != body_part_cells.size(); ++i)
			local_dynamics_function(body_part_cells[i]);
	};

	template <class LocalDynamicsFunction>
	inline void particle_for(const ParallelPolicy &par, const DataListsInCells &body_part_cells,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		parallel_for(
			IndexRange(0, body_part_cells.size()),
			[&](const IndexRange &r)
			{
				for (size_t i = r.begin(); i < r.end(); ++i)
				{
					local_dynamics_function(body_part_cells[i]);
				}
			},
			ap);
	};
	/**
	 * Splitting algorithm (for sequential and parallel computing).
	 */
	template <class LocalDynamicsFunction>
	inline void particle_for(const SequencedPolicy &seq, const SplitCellLists &split_cell_lists,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		// forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k)
		{
			const ConcurrentCellLists &cell_lists = split_cell_lists[k];
			for (size_t l = 0; l != cell_lists.size(); ++l)
			{
				const ConcurrentIndexVector &particle_indexes = *cell_lists[l];
				for (size_t i = 0; i != particle_indexes.size(); ++i)
				{
					local_dynamics_function(particle_indexes[i]);
				}
			}
		}

		// backward sweeping
		for (size_t k = split_cell_lists.size(); k != 0; --k)
		{
			const ConcurrentCellLists &cell_lists = split_cell_lists[k - 1];
			for (size_t l = 0; l != cell_lists.size(); ++l)
			{
				const ConcurrentIndexVector &particle_indexes = *cell_lists[l];
				for (size_t i = particle_indexes.size(); i != 0; --i)
				{
					local_dynamics_function(particle_indexes[i - 1]);
				}
			}
		}
	}

	template <class LocalDynamicsFunction>
	inline void particle_for(const ParallelPolicy &par, const SplitCellLists &split_cell_lists,
							 const LocalDynamicsFunction &local_dynamics_function)
	{
		// forward sweeping
		for (size_t k = 0; k != split_cell_lists.size(); ++k)
		{
			const ConcurrentCellLists &cell_lists = split_cell_lists[k];
			parallel_for(
				IndexRange(0, cell_lists.size()),
				[&](const IndexRange &r)
				{
					for (size_t l = r.begin(); l < r.end(); ++l)
					{
						const ConcurrentIndexVector &particle_indexes = *cell_lists[l];
						for (size_t i = 0; i < particle_indexes.size(); ++i)
						{
							local_dynamics_function(particle_indexes[i]);
						}
					}
				},
				ap);
		}

		// backward sweeping
		for (size_t k = split_cell_lists.size(); k != 0; --k)
		{
			const ConcurrentCellLists &cell_lists = split_cell_lists[k - 1];
			parallel_for(
				IndexRange(0, cell_lists.size()),
				[&](const IndexRange &r)
				{
					for (size_t l = r.begin(); l < r.end(); ++l)
					{
						const ConcurrentIndexVector &particle_indexes = *cell_lists[l];
						for (size_t i = particle_indexes.size(); i != 0; --i)
						{
							local_dynamics_function(particle_indexes[i - 1]);
						}
					}
				},
				ap);
		}
	}

    template <class LocalDynamicsFunction, class Proxy>
    inline void particle_for(const ParallelSYCLDevicePolicy & sycl_policy, const size_t &all_real_particles,
                             const LocalDynamicsFunction &local_dynamics_function, Proxy& proxy)
    {
        try {
            auto &sycl_queue = ExecutionQueue::getQueue();
            auto* kernel = proxy.get(sycl_policy);
			auto kernel_buffer = sycl::buffer<typename Proxy::Kernel>(kernel, 1);
            sycl_queue.submit([&](sycl::handler &cgh) {
                auto kernel_accessor = kernel_buffer.get_access(cgh, sycl::read_write);
                cgh.parallel_for(all_real_particles, [=](sycl::item<1> index) {
                    auto i = index.get_id();
                    local_dynamics_function(i, kernel_accessor[0]);
                });
            });
            sycl_queue.wait_and_throw();
        } catch (const sycl::exception &error) {
            std::cerr << error.what() << std::endl;
        }
    }

	template <class ExecutionPolicy, typename DynamicsRange, class ReturnType,
			  typename Operation, class LocalDynamicsFunction>
	void particle_reduce(const ExecutionPolicy &execution_policy, const DynamicsRange &dynamics_range,
						 const LocalDynamicsFunction &local_dynamics_function)
	{
		std::cout << "\n Error: ExecutionPolicy, DynamicsRange or LocalDynamicsFunction not defined for particle dynamics !" << std::endl;
		std::cout << __FILE__ << ':' << __LINE__ << std::endl;
		exit(1);
	};

	/**
	 * Body-wise reduce iterators (for sequential and parallel computing).
	 */
	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	inline ReturnType particle_reduce(const SequencedPolicy &seq, const size_t &all_real_particles,
									  ReturnType temp, Operation &&operation,
									  const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i < all_real_particles; ++i)
		{
			temp = operation(temp, local_dynamics_function(i));
		}
		return temp;
	}

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	inline ReturnType particle_reduce(const ParallelPolicy &par, const size_t &all_real_particles,
									  ReturnType temp, Operation &&operation,
									  const LocalDynamicsFunction &local_dynamics_function)
	{
		return parallel_reduce(
			IndexRange(0, all_real_particles),
			temp, [&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					temp0 = operation(temp0, local_dynamics_function(i));
				}
				return temp0; },
			[&](const ReturnType &x, const ReturnType &y) -> ReturnType
			{
				return operation(x, y);
			});
	};

    template <class ReturnType, typename Operation, class LocalDynamicsFunction, class Proxy>
    inline ReturnType particle_reduce(const ParallelPolicy &par_policy, const size_t &all_real_particles,
                                      ReturnType identity, Operation &&operation,
                                      const LocalDynamicsFunction &local_dynamics_function, Proxy& proxy)
    {
        auto& kernel = *proxy.get(par_policy);
        return parallel_reduce(
                IndexRange(0, all_real_particles),
                identity, [&](const IndexRange &r, ReturnType temp0) -> ReturnType
                {
                    for (size_t i = r.begin(); i != r.end(); ++i)
                    {
                        temp0 = operation(temp0, local_dynamics_function(i, kernel));
                    }
                    return temp0; },
                [&](const ReturnType &x, const ReturnType &y) -> ReturnType
                {
                    return operation(x, y);
                });
    };

    template <class ReturnType, typename Operation, class LocalDynamicsFunction, class Proxy>
    inline ReturnType particle_reduce(const ParallelSYCLDevicePolicy &sycl_policy, const size_t &all_real_particles,
                                      ReturnType identity, Operation&&,
                                      const LocalDynamicsFunction &local_dynamics_function, Proxy& proxy)
    {
        ReturnType result = identity;
        auto* kernel = proxy.get(sycl_policy);
        auto kernel_buffer = sycl::buffer<typename Proxy::Kernel>(kernel, 1);
        try {
            auto &sycl_queue = ExecutionQueue::getQueue();
            sycl::buffer<ReturnType> buffer_result(&result, 1);
            sycl_queue.submit([&](sycl::handler &cgh) {
                auto kernel_accessor = kernel_buffer.get_access(cgh, sycl::read_write);
                auto reduction_operator = sycl::reduction(buffer_result, cgh, typename std::remove_reference_t<Operation>::SYCLOp());
                cgh.parallel_for(sycl::range(all_real_particles), reduction_operator, [=](sycl::id<1> idx, auto& reduction) {
                    reduction.combine(local_dynamics_function(idx, kernel_accessor[0]));
                });
            });
            sycl_queue.wait_and_throw();
        } catch (const sycl::exception &error) {
            std::cerr << error.what() << std::endl;
        }
        return result;
    }

	/**
	 * BodypartByParticle-wise reduce iterators (for sequential and parallel computing).
	 */
	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	inline ReturnType particle_reduce(const SequencedPolicy &seq, const IndexVector &body_part_particles,
									  ReturnType temp, Operation &&operation,
									  const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i < body_part_particles.size(); ++i)
		{
			temp = operation(temp, local_dynamics_function(body_part_particles[i]));
		}
		return temp;
	}

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	inline ReturnType particle_reduce(const ParallelPolicy &par, const IndexVector &body_part_particles,
									  ReturnType temp, Operation &&operation,
									  const LocalDynamicsFunction &local_dynamics_function)
	{
		return parallel_reduce(
			IndexRange(0, body_part_particles.size()),
			temp,
			[&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t n = r.begin(); n != r.end(); ++n)
				{
					temp0 = operation(temp0, local_dynamics_function(body_part_particles[n]));
				}
				return temp0;
			},
			[&](const ReturnType &x, const ReturnType &y) -> ReturnType
			{
				return operation(x, y);
			});
	};
	/**
	 * BodypartByCell-wise reduce iterators (for sequential and parallel computing).
	 */
	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	inline ReturnType particle_reduce(const SequencedPolicy &seq, const ConcurrentCellLists &body_part_cells,
									  ReturnType temp, Operation &&operation,
									  const LocalDynamicsFunction &local_dynamics_function)
	{
		for (size_t i = 0; i != body_part_cells.size(); ++i)
		{
			ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
			for (size_t num = 0; num < particle_indexes.size(); ++num)
			{
				temp = operation(temp, local_dynamics_function(particle_indexes[num]));
			}
		}

		return temp;
	}

	template <class ReturnType, typename Operation, class LocalDynamicsFunction>
	inline ReturnType particle_reduce(const ParallelPolicy &par, const ConcurrentCellLists &body_part_cells,
									  ReturnType temp, Operation &&operation,
									  const LocalDynamicsFunction &local_dynamics_function)
	{
		return parallel_reduce(
			IndexRange(0, body_part_cells.size()),
			temp,
			[&](const IndexRange &r, ReturnType temp0) -> ReturnType
			{
				for (size_t i = r.begin(); i != r.end(); ++i)
				{
					ConcurrentIndexVector &particle_indexes = *body_part_cells[i];
					for (size_t num = 0; num < particle_indexes.size(); ++num)
					{
						temp0 = operation(temp0, local_dynamics_function(particle_indexes[num]));
					}
				}
				return temp0;
			},
			[&](const ReturnType &x, const ReturnType &y) -> ReturnType
			{ return operation(x, y); });
	}
}
#endif // PARTICLE_ITERATORS_H
