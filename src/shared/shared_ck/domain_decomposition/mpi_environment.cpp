#include "mpi_environment.h"

namespace SPH
{
#if SPHINXSYS_USE_MPI
namespace
{
/** Real is float or double depending on SPHINXSYS_USE_FLOAT; pick the matching MPI type. */
MPI_Datatype mpiRealType()
{
    return sizeof(Real) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT;
}
/** UnsignedInt is uint32_t or size_t depending on the build; pick the matching MPI type. */
MPI_Datatype mpiUnsignedIntType()
{
    return sizeof(UnsignedInt) == sizeof(uint64_t) ? MPI_UINT64_T : MPI_UINT32_T;
}
} // namespace
#endif // SPHINXSYS_USE_MPI
//=================================================================================================//
MPIEnvironment::MPIEnvironment(int &argc, char **&argv)
{
#if SPHINXSYS_USE_MPI
    int already_initialized = 0;
    MPI_Initialized(&already_initialized);
    if (!already_initialized)
    {
        MPI_Init(&argc, &argv);
        owns_mpi_ = true;
    }
    MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size_);
#else
    (void)argc;
    (void)argv;
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
MPIEnvironment::~MPIEnvironment()
{
#if SPHINXSYS_USE_MPI
    if (owns_mpi_)
    {
        MPI_Finalize();
    }
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
bool MPIEnvironment::isDistributedBuild()
{
#if SPHINXSYS_USE_MPI
    return true;
#else
    return false;
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
void MPIEnvironment::barrier() const
{
#if SPHINXSYS_USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
Real MPIEnvironment::allReduceMin(Real local_value) const
{
#if SPHINXSYS_USE_MPI
    Real global_value = local_value;
    MPI_Allreduce(&local_value, &global_value, 1, mpiRealType(), MPI_MIN, MPI_COMM_WORLD);
    return global_value;
#else
    return local_value;
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
Real MPIEnvironment::allReduceMax(Real local_value) const
{
#if SPHINXSYS_USE_MPI
    Real global_value = local_value;
    MPI_Allreduce(&local_value, &global_value, 1, mpiRealType(), MPI_MAX, MPI_COMM_WORLD);
    return global_value;
#else
    return local_value;
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
Real MPIEnvironment::allReduceSum(Real local_value) const
{
#if SPHINXSYS_USE_MPI
    Real global_value = local_value;
    MPI_Allreduce(&local_value, &global_value, 1, mpiRealType(), MPI_SUM, MPI_COMM_WORLD);
    return global_value;
#else
    return local_value;
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
UnsignedInt MPIEnvironment::allReduceSum(UnsignedInt local_value) const
{
#if SPHINXSYS_USE_MPI
    UnsignedInt global_value = local_value;
    MPI_Allreduce(&local_value, &global_value, 1, mpiUnsignedIntType(), MPI_SUM, MPI_COMM_WORLD);
    return global_value;
#else
    return local_value;
#endif // SPHINXSYS_USE_MPI
}
//=================================================================================================//
} // namespace SPH
