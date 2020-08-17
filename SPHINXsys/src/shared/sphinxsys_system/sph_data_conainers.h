/**
 * @file 	sph_data_conainers.h
 * @brief 	Set up of basic data structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_data_package.h"

using namespace std;

namespace SPH {
	/**
	 * @brief Preclaimed classes.
	 */
	class BaseMaterial;
	class SPHBody;
	class CellList;
	class BaseParticleData;
	class NeighborRelation;

	/**< Vector of Material. Note that vector of references are not allowed in c++.*/
	using MaterialVector = StdVec<BaseMaterial*>;
	/** Vector of SPH body. Note that vector of references are not allowed in c++.*/
	using SPHBodyVector = StdVec<SPHBody*>;
	typedef pair<SPHBody*, SPHBodyVector> SPHBodyContactMap;

	/** Index container with elements of size_t. */
	using IndexVector = StdVec<size_t>;
	/** Cell container with elements of Vecu. */
	using CellVector = StdVec<Vecu>;

	/** Concurrent particle indexes .*/
	using ConcurrentIndexVector = LargeVec<size_t>;
	/** Concurrent cell indexes.*/
	using ConcurrentCellVector = LargeVec<Vecu>;
	/** List data pair*/
	using ListData = pair<size_t, Vecd>;
	/** Cell list concurrent vector data. */
	using CellListDataVector = StdLargeVec<ListData>;
	/** Concurrent vector .*/
	template<class DataType>
	using ConcurrentVector = LargeVec<DataType>;

	/** Cell lists*/
	using CellLists = StdLargeVec<CellList*>;
	/** concurrent cell lists*/
	using ConcurrentCellLists = LargeVec<CellList*>;
	/** Split cell list for split algorithms. */
	using SplitCellLists = StdVec<ConcurrentCellLists>;
	/** Pair of point and volume. */
	using PositionsAndVolumes = vector<pair<Point, Real>>;
}
