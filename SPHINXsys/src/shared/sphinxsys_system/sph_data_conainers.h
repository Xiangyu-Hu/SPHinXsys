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
	class NeighboringParticle;
	class ReferenceNeighboringParticle;

	/**< Vector of Material. Note that vector of references are not allowed in c++.*/
	using MaterialVector = StdVec<BaseMaterial*>;
	/** Vector of SPH body. Note that vector of references are not allowed in c++.*/
	using SPHBodyVector = StdVec<SPHBody*>;	
	typedef pair<SPHBody*, SPHBodyVector> SPHBodyContactMap;
	typedef vector<SPHBodyContactMap> SPHBodyTopology;
	
	/** Index containner with elements of size_t. */
	using IndexVector = StdVec<size_t>;
	/** Cell containner with elements of Vecu. */
	using CellVector = StdVec<Vecu>;		

	/** Concurrent particle indexes .*/
	using ConcurrentIndexVector = LargeVec<size_t>;
	/** Concurrent cell indexes.*/
	using ConcurrentCellVector = LargeVec<Vecu>;
	/** List data pair*/
	using ListData = pair<size_t, Vecd>;
	/** Cell list concurrent vector data. */
	using ConcurrentListDataVector = LargeVec<ListData>;

	/** Split cell list for split algorithms. */
	using SplitCellLists = StdVec<StdLargeVec<CellList*>>;
	/** Pair of point and volume. */
	using PositionsAndVolumes =vector<pair<Point, Real>> ; 

	/** Neighboring particles list. 
	  * Using pointer for overloading derived neighbor relations.*/
	using NeighborList = StdVec<NeighboringParticle*>; 
	/** Neighboring particle list, the current and the previous number of neighbors. */	
	using Neighborhood = tuple<NeighborList, size_t, size_t>;
	/** A neighborhoods for all particles in a body. */
	using ParticleConfiguration = StdLargeVec<Neighborhood>;
	/** Inner neighborhoods for all particles in a body. */
	using InnerParticleConfiguration = ParticleConfiguration;
	/** All contact neighborhoods for all particles in a body. */
	using ContatcParticleConfiguration = StdVec<ParticleConfiguration>;
	/** Interacting neighborhoods for all particles in a body. */
	using InteractingParticleConfiguration = StdVec<ParticleConfiguration*>;

	/** List of partilces contact to another body. */
	using ContactParticleList = ConcurrentIndexVector;
	/** All contact particles lists. **/
	using ContactParticles = StdVec<ContactParticleList>;
	/** Contact particles to interacting bodies. **/
	using InteractingParticles = StdVec<ContactParticleList*>;

}