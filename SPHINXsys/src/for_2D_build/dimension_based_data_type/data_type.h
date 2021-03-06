#pragma once
#include "small_vectors.h"
#include "large_data_containers.h"

namespace SPH {

	//for 2d build
	using Veci = Vec2i;
	using Vecu = Vec2u;
	using Vecd = Vec2d;
	using Matd = Mat2d;
	using SymMatd = SymMat2d;
	using AngularVecd = Real;
	const int indexAngularVector = 0;


	using Transformd = Transform2d;

	class CellList;
	using matrix_cell = CellList **;

	template<class DataType, int ARRAY_SIZE>
	using PackageDataMatrix = std::array<std::array<DataType, ARRAY_SIZE>, ARRAY_SIZE>;

	template<class DataType>
	using MeshDataMatrix = DataType**;

	/** only works for smoothing length ratio less or equal than 1.3*/
	constexpr int MaximumNeighborhoodSize = int(M_PI * 9);
	const int Dimensions = 2;

	/** correction matrix, only works for thin structure dynamics. */
	const Matd correction_matrix = { 1, 0, 0, 0 };

	/** initial locald normal, only works for thin structure dynamics. */
	const Vecd n_local_0 = Vecd(0.0, 1.0);
}