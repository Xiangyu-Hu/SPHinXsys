/**
 * @file 	sph_data_conainers.h
 * @brief 	Set up of basic data structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 */

#ifndef SPH_DATA_CONTAINERS_H
#define SPH_DATA_CONTAINERS_H



#include "base_data_package.h"

using namespace std;

namespace SPH {
	/**
	 * @brief Preclaimed classes.
	 */
	class BaseMaterial;
	class SPHBody;
	class RealBody;
	class FictitiousBody;
	class CellList;
	class BaseParticles;

	/** Bounding box for system, body, body part and shape, first: lower bound, second: upper bound. */
	typedef pair<Vecd, Vecd> BoundingBox;
	/** Vector of Material. Note that vector of references are not allowed in c++.*/
	using MaterialVector = StdVec<BaseMaterial*>;
	/** Vector of bodys */
	using SPHBodyVector = StdVec<SPHBody*>;
	using RealBodyVector = StdVec<RealBody*>;
	using FictitiousBodyVector = StdVec<FictitiousBody*>;

	/** Index container with elements of size_t. */
	using IndexVector = StdVec<size_t>;
	/** Concurrent particle indexes .*/
	using ConcurrentIndexVector = LargeVec<size_t>;

	/** List data pair*/
	using ListData = pair<size_t, Vecd>;
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
	using PositionsAndVolumes = vector<pair<Vecd, Real>>;
}
#endif //SPH_DATA_CONTAINERS_H