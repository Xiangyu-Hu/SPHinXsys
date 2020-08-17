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

	template<class DataType>
	using MeshDataMatrix = DataType**;
}