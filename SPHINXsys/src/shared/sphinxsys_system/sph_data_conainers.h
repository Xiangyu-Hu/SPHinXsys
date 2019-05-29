/**
 * @file 	sph_data_conainers.h
 * @brief 	Set up of basic data structure.
 * @author	Luhui Han, Chi ZHang and Xiangyu Hu
 * @version	0.1
 */
#pragma once

#include "base_data_package.h"

namespace SPH {
	/**
	 * @brief Friend classes.
	 */
	class ListData;
	class SPHBody;
	
	using IndexVector = StdVec<size_t>;		/**< Index containner with elements of size_t. */
	using CellVector = StdVec<Vecu>;		/**< Cell containner with elements of Vecu. */
	using SPHBodyVector = StdVec<SPHBody*>;	/**< Vector of SPH body. */

	using ListIndexVector = LargeVec<size_t>;	/**< Concurrent body index vector for cell linked lists. */

	using ConcurrentIndexVector = LargeVec<size_t>;	/**< Concurrent particle index which can be used to store particle position of a container .*/
	using ListDataVector = LargeVec<ListData>;		/**< List of vector data. */

	// index set
	using IndexSet = LargeSet<size_t>;				/**< Index set. */
	using ContactParticleList = ListIndexVector *;	/**< List of contact interacing partilces. */
	using ByCellLists = StdVec<IndexVector> *;		/**< Index vector for by cell list. */
	using PositionsAndVolumes = std::vector<std::pair<Point, Real>> ; /**< Pair of point and volume. */

}