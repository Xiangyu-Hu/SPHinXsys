/**
 * @file 	sph_data_containers.h
 * @brief 	Set up of basic data structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#ifndef SPH_DATA_CONTAINERS_H
#define SPH_DATA_CONTAINERS_H

#include "base_data_package.h"
#include "base_data_type.h"

namespace SPH
{
	/**
	 * @brief Preclaimed classes.
	 */
	class BaseMaterial;
	class SPHBody;
	class RealBody;
	class SolidBody;
	class BodyPart;
	class FictitiousBody;
	class CellList;
	class BaseParticles;

	/** Bounding box for system, body, body part and shape, first: lower bound, second: upper bound. */
	typedef std::pair<Vecd, Vecd> BoundingBox;
	/** Check if a point is inside the bounding box */
	bool checkIfPointInBoundingBox(Vec3d point, BoundingBox &bbox);
	bool checkIfPointInBoundingBox(Vec2d point, BoundingBox &bbox);

	/** Generalized particle data type */
	typedef std::tuple<StdVec<StdLargeVec<Real> *>, StdVec<StdLargeVec<Vecd> *>, StdVec<StdLargeVec<Matd> *>,
					   StdVec<StdLargeVec<int> *>>
		ParticleData;
	/** Generalized particle variable to index map */
	typedef std::array<std::map<std::string, size_t>, 4> ParticleDataMap;
	/** Generalized particle variable list */
	typedef std::array<StdVec<std::pair<std::string, size_t>>, 4> ParticleVariableList;
	/** Vector of Material. Note that vector of references are not allowed in c++.*/
	using MaterialVector = StdVec<BaseMaterial *>;
	/** Vector of bodies */
	using SPHBodyVector = StdVec<SPHBody *>;
	using SolidBodyVector = StdVec<SolidBody *>;
	using RealBodyVector = StdVec<RealBody *>;
	using BodyPartVector = StdVec<BodyPart *>;
	using FictitiousBodyVector = StdVec<FictitiousBody *>;

	/** Index container with elements of size_t. */
	using IndexVector = StdVec<size_t>;
	/** Concurrent particle indexes .*/
	using ConcurrentIndexVector = LargeVec<size_t>;

	/** List data pair */
	using ListData = std::pair<size_t, Vecd>;
	/** Vector of list data pair */
	using ListDataVector = StdLargeVec<ListData>;
	/** Cell lists*/
	using CellLists = StdLargeVec<CellList *>;

	/** Concurrent vector .*/
	template <class DataType>
	using ConcurrentVector = LargeVec<DataType>;
	/** concurrent cell lists*/
	using ConcurrentCellLists = LargeVec<CellList *>;
	/** Split cell list for split algorithms. */
	using SplitCellLists = StdVec<ConcurrentCellLists>;
	/** Pair of point and volume. */
	using PositionsVolumes = StdVec<std::pair<Vecd, Real>>;

	/** an operation on particle data of all types */
	template <template <int DataTypeIndex, typename VariableType> typename OperationType>
	struct ParticleDataOperation
	{
		OperationType<indexScalar, Real> scalar_operation;
		OperationType<indexVector, Vecd> vector_operation;
		OperationType<indexMatrix, Matd> matrix_operation;
		OperationType<indexInteger, int> integer_operation;

		template <typename... ParticleArgs>
		void operator()(ParticleData &particle_data, ParticleArgs... particle_args)
		{
			scalar_operation(particle_data, particle_args...);
			vector_operation(particle_data, particle_args...);
			matrix_operation(particle_data, particle_args...);
			integer_operation(particle_data, particle_args...);
		}
	};

	template <int DataTypeIndex, typename VariableType>
	struct swapParticleDataValue
	{
		void operator()(ParticleData &particle_data, size_t index_a, size_t index_b) const
		{
			StdVec<StdLargeVec<VariableType> *> variables = std::get<DataTypeIndex>(particle_data);
			for (size_t i = 0; i != variables.size(); ++i)
			{
				StdLargeVec<VariableType> &variable = *variables[i];
				std::swap(variable[index_a], variable[index_b]);
			}
		};
	};
	
	/** operation by looping or going through a particle data map */
	template <int DataTypeIndex, typename VariableType>
	struct loopParticleDataMap
	{
		template <typename VariableOperation>
		void operator()(ParticleData &particle_data,
						ParticleDataMap &particle_data_map, VariableOperation &variable_operation) const
		{
			for (auto const &name_index : particle_data_map[DataTypeIndex])
			{
				std::string variable_name = name_index.first;
				StdLargeVec<VariableType> &variable = *(std::get<DataTypeIndex>(particle_data)[name_index.second]);
				variable_operation(variable_name, variable);
			}
		};
	};

	/** operation by looping or going through a variable name list */
	template <int DataTypeIndex, typename VariableType>
	struct loopVariabaleNameList
	{
		template <typename VariableOperation>
		void operator()(ParticleData &particle_data,
						ParticleVariableList &variable_name_list, VariableOperation &variable_operation) const
		{
			for (std::pair<std::string, size_t> &name_index : variable_name_list[DataTypeIndex])
			{
				std::string variable_name = name_index.first;
				StdLargeVec<VariableType> &variable = *(std::get<DataTypeIndex>(particle_data)[name_index.second]);
				variable_operation(variable_name, variable);
			}
		};
	};
}
#endif //SPH_DATA_CONTAINERS_H