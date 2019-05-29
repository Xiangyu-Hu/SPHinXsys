#pragma once
#include "small_vectors.h"
#include "large_data_containers.h"

namespace SPH {

	class CellList;
	class BackgroundData;

	//for 2d build
	using Veci = Vec2i;
	using Vecu = Vec2u;
	using Vecd = Vec2d;
	using Point = Vec2d;
	using Index = Vec2i;
	using Matd = Mat2d;
	using SymMatd = SymMat2d;

	using matrix_i = int **;
	using Angular = SimTK::Real;

	using matrix_cell = CellList **;
	using matrix_grid = BackgroundData **;

}