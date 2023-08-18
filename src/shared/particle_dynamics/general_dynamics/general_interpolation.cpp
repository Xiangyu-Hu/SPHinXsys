/**
 * @file 	general_interpolation.cpp
 * @brief
 * @author	Chi Zhang and Xiangyu Hu
 */

#include "general_interpolation.h"

namespace SPH
{
//=================================================================================================//
CorrectInterpolationKernelWeights::
    CorrectInterpolationKernelWeights(BaseContactRelation &contact_relation)
    : LocalDynamics(contact_relation.getSPHBody()),
      InterpolationContactData(contact_relation)
{
    for (size_t k = 0; k != contact_particles_.size(); ++k)
    {
        contact_Vol_.push_back(&(contact_particles_[k]->Vol_));
    }
}
//=================================================================================================//
} // namespace SPH
