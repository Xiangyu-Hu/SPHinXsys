#pragma once

#include "small_vectors.h"
#include "large_data_containers.h"

namespace SPH {

	//for 3d build
	using Veci = Vec3i;
	using Vecu = Vec3u;
	using Vecd = Vec3d;
	using Point = Vec3d;
	using Index = Vec3i;
	using Matd = Mat3d;
	using SymMatd = SymMat3d;

	using Transformd = SimTK::Transform;

	class CellList;
	using matrix_cell = CellList ***;

	template<class DataType, int ARRAY_SIZE>
	using PackageDataMatrix = std::array<std::array<std::array<DataType, ARRAY_SIZE>, ARRAY_SIZE>, ARRAY_SIZE>;

	template<class DataType>
	using MeshDataMatrix = DataType***;

}