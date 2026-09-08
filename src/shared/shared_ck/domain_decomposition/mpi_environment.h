/**
 * @file    mpi_environment.h
 * @brief   Thin wrapper over MPI initialisation and the collectives SPHinXsys needs.
 * @details When SPHINXSYS_USE_MPI is off, every member degrades to the single-rank answer
 *          (rank 0 of 1, reductions return their input). A case written against this API
 *          therefore still compiles and runs serially without #ifdefs at the call site.
 * @author  SPHinXsys multi-GPU draft
 */

#ifndef MPI_ENVIRONMENT_H
#define MPI_ENVIRONMENT_H

#include "data_type.h"

#if SPHINXSYS_USE_MPI
#include <mpi.h>
#endif // SPHINXSYS_USE_MPI

namespace SPH
{
/**
 * @class MPIEnvironment
 * @brief Owns MPI_Init/MPI_Finalize for the lifetime of the object.
 *
 * Construct once, near the top of main(), before any SPHSystem is built.
 */
class MPIEnvironment
{
  public:
    MPIEnvironment(int &argc, char **&argv);
    ~MPIEnvironment();

    MPIEnvironment(const MPIEnvironment &) = delete;
    MPIEnvironment &operator=(const MPIEnvironment &) = delete;

    int Rank() const { return rank_; };
    int WorldSize() const { return world_size_; };
    bool isMainRank() const { return rank_ == 0; };

    /** Whether this build actually has MPI compiled in. */
    static bool isDistributedBuild();

    void barrier() const;

    /**
     * Collective reductions. The CFL time step is the motivating case: every rank computes
     * a local dt from its own particles, and without an all-reduce the ranks would advance
     * by different amounts and desynchronise.
     */
    Real allReduceMin(Real local_value) const;
    Real allReduceMax(Real local_value) const;
    Real allReduceSum(Real local_value) const;
    UnsignedInt allReduceSum(UnsignedInt local_value) const;

  protected:
    int rank_ = 0;
    int world_size_ = 1;
    bool owns_mpi_ = false; /**< false if MPI was already initialised by the host application */
};
} // namespace SPH
#endif // MPI_ENVIRONMENT_H
