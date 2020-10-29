#pragma once
#include "small_vectors.h"
#include "large_data_containers.h"

namespace SPH {

	//for 2d build
	using Veci = Vec2i;
	using Vecu = Vec2u;
	using Vecd = Vec2d;
	using Point = Vec2d;
	using Index = Vec2i;
	using Matd = Mat2d;
	using SymMatd = SymMat2d;

	using Transformd = Transform2d;

	class CellList;
	using matrix_cell = CellList **;

	template<class DataType, int ARRAY_SIZE>
	using PackageDataMatrix = std::array<std::array<DataType, ARRAY_SIZE>, ARRAY_SIZE>;

	template<class DataType>
	using MeshDataMatrix = DataType**;

	/** only works for smoothing length ratio less or equal than 1.3*/
	const int MaximumNeighborhoodSize = int(Pi * 9);
}