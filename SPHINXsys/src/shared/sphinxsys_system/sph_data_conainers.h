/**
 * @file 	sph_data_conainers.h
 * @brief 	Set up of basic data structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */
#pragma once

#include "base_data_package.h"

namespace SPH {
	/**
	 * @brief Preclaimed classes.
	 */
	class BaseMaterial;
	class SPHBody;
	class RealBody;
	class BodyPart;
	class FictitiousBody;
	class CellList;
	class BaseParticles;

	/** Bounding box for system, body, body part and shape, first: lower bound, second: upper bound. */
	typedef std::pair<Vecd, Vecd> BoundingBox;
	/** Generalized particle data type */
	typedef std::tuple<StdVec<StdLargeVec<Real>*>, StdVec<StdLargeVec<Vecd>*>, StdVec<StdLargeVec<Matd>*>,
		StdVec<StdLargeVec<int>*>> ParticleData;
	/** Generalized particle variable to index map */
	typedef std::array<std::map<std::string, size_t>, 4> ParticleDataMap;
	/** Generalized particle variable list */
	typedef std::array<StdVec<std::pair<std::string, size_t>>, 4> ParticleVariableList;
	/** Vector of Material. Note that vector of references are not allowed in c++.*/
	using MaterialVector = StdVec<BaseMaterial*>;
	/** Vector of bodies */
	using SPHBodyVector = StdVec<SPHBody*>;
	using RealBodyVector = StdVec<RealBody*>;
	using BodyPartVector = StdVec<BodyPart*>;
	using FictitiousBodyVector = StdVec<FictitiousBody*>;

	/** Index container with elements of size_t. */
	using IndexVector = StdVec<size_t>;
	/** Concurrent particle indexes .*/
	using ConcurrentIndexVector = LargeVec<size_t>;

	/** List data pair*/
	using ListData = std::pair<size_t, Vecd>;
	/** Cell list vector data. */
	using CellListDataVector = StdLargeVec<ListData>;
	/** Cell lists*/
	using CellLists = StdLargeVec<CellList*>;

	/** Concurrent vector .*/
	template<class DataType>
	using ConcurrentVector = LargeVec<DataType>;
	/** concurrent cell lists*/
	using ConcurrentCellLists = LargeVec<CellList*>;
	/** Split cell list for split algorithms. */
	using SplitCellLists = StdVec<ConcurrentCellLists>;
	/** Pair of point and volume. */
	using PositionsAndVolumes = StdVec<std::pair<Vecd, Real>>;

	/** register a variable to generalized particle data and particle map */
	template<int DataTypeIndex, typename VariableType>
	void registerAVariableToParticleData(ParticleData& particle_data, 
		ParticleDataMap& particle_map, StdLargeVec<VariableType>& variable_addrs, std::string variable_name)
	{
		if (particle_map[DataTypeIndex].find(variable_name) == particle_map[DataTypeIndex].end())
		{
			std::get<DataTypeIndex>(particle_data).push_back(&variable_addrs);
			particle_map[DataTypeIndex].insert(make_pair(variable_name, std::get<DataTypeIndex>(particle_data).size() - 1));
		}
		else
		{
			std::cout << "\n Error: the variable '" << variable_name << "' has already been registered!" << std::endl;
			exit(1);
		}
	};

	/** loop particle data with operations */
	template<template<int DataTypeIndex, typename VariableType> typename OperationType, 
		typename... ParticleArgs>
	void loopParticleData(ParticleData& particle_data, ParticleArgs... particle_args)
	{
		OperationType<indexScalar, Real>	scalar_operation;
		OperationType<indexVector, Vecd>	vector_operation;
		OperationType<indexMatrix, Matd>	matrix_operation;
		OperationType<indexInteger, int>	integer_operation;

		scalar_operation(particle_data, particle_args...);
		vector_operation(particle_data, particle_args...);
		matrix_operation(particle_data, particle_args...);
		integer_operation(particle_data, particle_args...);
	};

		/** operation by looping or going through a particle data map */
		template<int DataTypeIndex, typename VariableType>
		struct loopParticleDataMap
		{
			template<typename VariableOperation>
			void operator () (ParticleData& particle_data,
				ParticleDataMap& particle_data_map, VariableOperation& variable_operation) const
			{
				for (auto const& name_index : particle_data_map[DataTypeIndex])
				{
					std::string variable_name = name_index.first;
					StdLargeVec<VariableType>& variable = *(std::get<DataTypeIndex>(particle_data)[name_index.second]);
					variable_operation(variable_name, variable);
				}
			};
		};
}
