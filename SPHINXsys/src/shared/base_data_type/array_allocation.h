#ifndef ARRAYALLOCATE_H
#define ARRAYALLOCATE_H

#include "small_vectors.h"

namespace SPH {
	//-------------------------------------------------------------------------------------------------
	//Allocate and deallocate 3d array
	//-------------------------------------------------------------------------------------------------
	template<class T>
	void Allocate3dArray(T*** &matrix, Vec3u res)
	{
		matrix = new T**[res[0]];
		for (size_t i = 0; i < res[0]; i++) {
			matrix[i] = new T*[res[1]];
			for (size_t j = 0; j < res[1]; j++) {
				matrix[i][j] = new T[res[2]];
			}
		}
	}

	template<class T>
	void Delete3dArray(T*** matrix, Vec3u res)
	{
		for (size_t i = 0; i < res[0]; i++) {
			for (size_t j = 0; j < res[1]; j++) {
				delete[] matrix[i][j];
			}
			delete[] matrix[i];
		}
		delete[] matrix;
	}
	//-------------------------------------------------------------------------------------------------
	//							Allocate 2d array
	//-------------------------------------------------------------------------------------------------
	template<class T>
	void Allocate2dArray(T** &matrix, Vec2u res)
	{
		matrix = new T*[res[0]];
		for (size_t i = 0; i < res[0]; i++) {
			matrix[i] = new T[res[1]];
		}
	}
	template<class T>
	void Delete2dArray(T** matrix, Vec2u res)
	{
		for (size_t i = 0; i < res[0]; i++) {
			delete[] matrix[i];
		}
		delete[] matrix;
	}

	//-------------------------------------------------------------------------------------------------
	//							Allocate 1d array
	//-------------------------------------------------------------------------------------------------
	template<class T>
	void Allocate1dArray(T* &matrix, size_t res)
	{
		matrix = new T[res];
	}
	template<class T>
	void Delete1dArray(T* matrix, size_t res)
	{
		delete[] matrix;
	}

}


#endif //ARRAYALLOCATE_H