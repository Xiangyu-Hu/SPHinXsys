
#ifndef DATA_TYPE_2D_H
#define DATA_TYPE_2D_H


#include "base_data_type.h"

namespace SPH {

	//for 2d build
	using Veci = Vec2i;
	using Vecu = Vec2u;
	using Vecd = Vec2d;
	using Matd = Mat2d;
	using SymMatd = SymMat2d;
	using AngularVecd = Real;

	using Transformd = Transform2d;

	template<class DataType, int ARRAY_SIZE>
	using PackageDataMatrix = std::array<std::array<DataType, ARRAY_SIZE>, ARRAY_SIZE>;

	template<class DataType>
	using MeshDataMatrix = DataType**;

	/** only works for smoothing length ratio less or equal than 1.3*/
	constexpr int MaximumNeighborhoodSize = int(M_PI * 9);
	const int Dimensions = 2;

	/** correction matrix, only works for thin structure dynamics. */
	const Matd reduced_unit_matrix = { 1, 0, 0, 0 };

	/** initial local normal, only works for thin structure dynamics. */
	const Vecd local_pseudo_n_0 = Vecd(0.0, 1.0);
}

#endif //DATA_TYPE_2D_H