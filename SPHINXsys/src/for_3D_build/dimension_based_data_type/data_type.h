#pragma once
#include "small_vectors.h"
#include "large_data_containers.h"

namespace SPH {

	class CellList;
	class BackgroundData;

	//for 3d build
	using Veci = Vec3i;
	using Vecu = Vec3u;
	using Vecd = Vec3d;
	using Point = Vec3d;
	using Index = Vec3i;
	using Matd = Mat3d;
	using SymMatd = SymMat3d;

	using matrix_i = int ***;
	using Rotationd = SimTK::Real;
	using Transformd = SimTK::Transform;

	using matrix_cell = CellList ***;
	using matrix_grid = BackgroundData ***;

	template<class DataType>
	using MeshDataMatrix = DataType***;
	template<size_t PACK_SIZE, class DataType>
	using MeshData = MeshData3<PACK_SIZE, DataType>;
}