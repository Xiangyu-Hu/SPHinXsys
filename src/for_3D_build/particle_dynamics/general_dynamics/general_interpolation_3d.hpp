/**
 * @file 	general_interpolation_3d.cpp
 * @brief
 * @author	Bo Zhang and Xiangyu Hu
 */

#ifndef GENERAL_INTERPOLATION_3D_HPP
#define GENERAL_INTERPOLATION_3D_HPP

#include "base_general_dynamics.h"
#include "general_interpolation.h"

namespace SPH
{
    template <typename DataType>
    void BaseInterpolation<DataType>::interaction(size_t index_i, Real dt)
    {
        DataType observed_quantity = ZeroData<DataType>::value;
        DataType gradientXprediction = ZeroData<DataType>::value;
        DataType gradientYprediction = ZeroData<DataType>::value;
        DataType gradientZprediction = ZeroData<DataType>::value;
        Mat4d restoring_matrix = Eps * Mat4d::Identity();
        Mat4d restoring_matrix_inverse = Mat4d::Identity();

        Real matrix_element_1 = 0.0;
        Vec3d matrix_element_2 = Vec3d::Zero();
        Vec3d matrix_element_3 = Vec3d::Zero();
        Mat3d matrix_element_4 = Mat3d::Zero();

        for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
        {
            Real* Vol_k = contact_Vol_[k];
            DataType* data_k = contact_data_[k];
            Neighborhood& contact_neighborhood = (*this->contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                size_t index_j = contact_neighborhood.j_[n];
                Vecd e_ij = contact_neighborhood.e_ij_[n];
                Vecd r_ji = contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
                Real W_ij = contact_neighborhood.W_ij_[n];
                Real dW_ij = contact_neighborhood.dW_ij_[n];

                Real element1 = W_ij * Vol_k[index_j];
                Vec3d element2 = W_ij * Vol_k[index_j] * r_ji;
                Vec3d element3 = dW_ij * Vol_k[index_j] * e_ij;
                Mat3d element4 = dW_ij * Vol_k[index_j] * r_ji * e_ij.transpose();

                observed_quantity += element1 * data_k[index_j];
                gradientXprediction += element3[0] * data_k[index_j];
                gradientYprediction += element3[1] * data_k[index_j];
                gradientZprediction += element3[2] * data_k[index_j];

                matrix_element_1 += element1;
                matrix_element_2 -= element2;
                matrix_element_3 += element3;
                matrix_element_4 -= element4;
            }
        }
        restoring_matrix(0, 0) = matrix_element_1;
        restoring_matrix(0, 1) = matrix_element_2(0);
        restoring_matrix(0, 2) = matrix_element_2(1);
        restoring_matrix(0, 3) = matrix_element_2(2);
        restoring_matrix(1, 0) = matrix_element_3(0);
        restoring_matrix(1, 1) = matrix_element_4(0, 0);
        restoring_matrix(1, 2) = matrix_element_4(0, 1);
        restoring_matrix(1, 3) = matrix_element_4(0, 2);
        restoring_matrix(2, 0) = matrix_element_3(1);
        restoring_matrix(2, 1) = matrix_element_4(1, 0);
        restoring_matrix(2, 2) = matrix_element_4(1, 1);
        restoring_matrix(2, 3) = matrix_element_4(1, 2);
        restoring_matrix(3, 0) = matrix_elemant_3(2);
        restoring_matrix(3, 1) = matrix_element_4(2, 0);
        restoring_matrix(3, 2) = matrix_element_4(2, 1);
        restoring_matrix(3, 3) = matrix_element_4(2, 2);
        restoring_matrix_inverse = restoring_matrix.inverse();

        interpolated_quantities_[index_i] = restoring_matrix_inverse(0, 0) * observed_quantity +
            restoring_matrix_inverse(0, 1) * gradientXprediction +
            restoring_matrix_inverse(0, 2) * gradientYprediction +
            restoring_matrix_inverse(0, 3) * gradientZprediction;
//=================================================================================================//
    } // namespace SPH
}//=================================================================================================//
#endif // GENERAL_INTERPOLATION_3D_HPP
