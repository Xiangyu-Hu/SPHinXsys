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
 * @file 	sph_data_containers.h
 * @brief 	Set up of basic data structure.
 * @author	Chi ZHang and Xiangyu Hu
 */

#ifndef SPH_DATA_CONTAINERS_H
#define SPH_DATA_CONTAINERS_H

#include "base_data_package.h"
#include "base_data_type.h"

namespace SPH
{
	/**
	 * @brief Pre-claimed classes.
	 */
	class BaseMaterial;
	class SPHBody;
	class RealBody;
	class SolidBody;
	class BodyPart;
	class BaseParticles;

	using MaterialVector = StdVec<BaseMaterial *>;
	using SPHBodyVector = StdVec<SPHBody *>;
	using SolidBodyVector = StdVec<SolidBody *>;
	using RealBodyVector = StdVec<RealBody *>;
	using BodyPartVector = StdVec<BodyPart *>;

	using IndexVector = StdVec<size_t>;
	using ConcurrentIndexVector = ConcurrentVec<size_t>;

	/** List data pair: first for indexes, second for particle position. */
	using ListData = std::tuple<size_t, Vecd, Real>;
	using ListDataVector = StdLargeVec<ListData>;
	using ConcurrentIndexesInCells = StdLargeVec<ConcurrentIndexVector *>;
	using DataListsInCells = StdLargeVec<ListDataVector *>;
	using CellLists = std::pair<ConcurrentIndexesInCells, DataListsInCells>;

	using ConcurrentCellLists = ConcurrentVec<ConcurrentIndexVector *>;
	/** Cell list for splitting algorithms. */
	using SplitCellLists = StdVec<ConcurrentCellLists>;

	/** Generalized particle data type */
	typedef DataContainerAddressAssemble<StdLargeVec> ParticleData;
	constexpr size_t ParticleDataSize = std::tuple_size<ParticleData>::value;
	/** Generalized particle variable to index map */
	typedef std::array<std::map<std::string, size_t>, ParticleDataSize> ParticleDataMap;
	/** Generalized particle variable list */
	typedef std::array<StdVec<std::pair<std::string, size_t>>, ParticleDataSize> ParticleVariableList;

	/** operation by looping or going through a particle data map */
	template <typename VariableType>
	struct loopParticleDataMap
	{
		template <typename VariableOperation>
		void operator()(ParticleData &particle_data,
						ParticleDataMap &particle_data_map, VariableOperation &variable_operation) const
		{
			constexpr int type_index = DataTypeIndex<VariableType>::value;
			for (auto const &name_index : particle_data_map[type_index])
			{
				std::string variable_name = name_index.first;
				StdLargeVec<VariableType> &variable = *(std::get<type_index>(particle_data)[name_index.second]);
				variable_operation(variable_name, variable);
			}
		};
	};

	/** operation by looping or going through a variable name list */
	template <typename VariableType>
	struct loopVariableNameList
	{
		template <typename VariableOperation>
		void operator()(ParticleData &particle_data,
						ParticleVariableList &variable_name_list, VariableOperation &variable_operation) const
		{
			constexpr int type_index = DataTypeIndex<VariableType>::value;
			for (std::pair<std::string, size_t> &name_index : variable_name_list[type_index])
			{
				std::string variable_name = name_index.first;
				StdLargeVec<VariableType> &variable = *(std::get<type_index>(particle_data)[name_index.second]);
				variable_operation(variable_name, variable);
			}
		};
	};
}
#endif // SPH_DATA_CONTAINERS_H
