/**
 * @file 	common_compressible_FVM_classes.h
 * @brief 	Here, we define the common compressible classes for fluid dynamics in FVM.
 * @author	Zhentong Wang and Xiangyu Hu
 */

#ifndef COMMON_COMPRESSIBLE_FVM_CLASSES_H
#define COMMON_COMPRESSIBLE_FVM_CLASSES_H

#include "sphinxsys.h"

namespace SPH
{
/**
 * @class CompressibleAcousticTimeStepSizeInFVM
 * @brief Computing the acoustic time step size
 */
class CompressibleAcousticTimeStepSizeInFVM : public fluid_dynamics::AcousticTimeStep
{
  protected:
    Real *rho_, *p_;
    Vecd *vel_;
    Real min_distance_between_nodes_;

  public:
    explicit CompressibleAcousticTimeStepSizeInFVM(SPHBody &sph_body, Real min_distance_between_nodes, Real acousticCFL = 0.6);
    virtual ~CompressibleAcousticTimeStepSizeInFVM(){};
    Real reduce(size_t index_i, Real dt = 0.0);
    virtual Real outputResult(Real reduced_value) override;
    Fluid &compressible_fluid_;
    Real acousticCFL_;
};
} // namespace SPH
#endif // COMMON_COMPRESSIBLE_FVM_CLASSES_H
