
#ifndef DATA_TYPE_3D_H
#define DATA_TYPE_3D_H



#include "base_data_type.h"

namespace SPH {

	//for 3d build
	using Veci = Vec3i;
	using Vecu = Vec3u;
	using Vecd = Vec3d;
	using Matd = Mat3d;
	using SymMatd = SymMat3d;
	using AngularVecd = Vec3d;

	using Transformd = Transform3d;

	template<class DataType, int ARRAY_SIZE>
	using PackageDataMatrix = std::array<std::array<std::array<DataType, ARRAY_SIZE>, ARRAY_SIZE>, ARRAY_SIZE>;

	template<class DataType>
	using MeshDataMatrix = DataType***;

	/** only works for smoothing length ratio less or equal than 1.3*/
	constexpr int MaximumNeighborhoodSize = int(1.33 * M_PI * 27);

	const int Dimensions = 3;

	/** correction matrix, only works for thin structure dynamics. */
	const Matd reduced_unit_matrix = { 1, 0, 0, 0, 1, 0, 0, 0, 0 };

	/** initial local normal, only works for thin structure dynamics. */
	const Vecd local_pseudo_n_0 = Vecd(0.0, 0.0, 1.0);
}
#endif //DATA_TYPE_3D_H