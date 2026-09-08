/**
 * @file    halo_exchange_ck.h
 * @brief   Exchange of boundary ("halo") particles between neighbouring MPI ranks.
 * @details Structurally this is the distributed-memory analogue of
 *          PeriodicConditionUsingGhostParticles: the same split into an occasional
 *          "establish the pattern" phase and a per-step "refresh the payload" phase.
 *
 *              periodic (single rank)         |  halo (this class)
 *              -------------------------------|---------------------------------
 *              CreatPeriodicGhostParticles    |  buildSendLists() + negotiateCounts()
 *              UpdatePeriodicGhostParticles   |  exchangeState()
 *
 *          Halo particles are appended above the owned particles and counted as real
 *          particles, so UpdateCellLinkedList and the relation update pick them up with
 *          no modification (see MULTI_GPU_DESIGN.md section 6, option A). Physics
 *          dynamics must therefore loop over svTotalOwnedParticles(), not
 *          svTotalRealParticles(), or they will integrate halo slots.
 * @author  SPHinXsys multi-GPU draft
 */

#ifndef HALO_EXCHANGE_CK_H
#define HALO_EXCHANGE_CK_H

#include "base_body.h"
#include "base_particles.h"
#include "domain_partition.h"
#include "mpi_environment.h"
#include "sphinxsys_variable.h"

#include <memory>
#include <vector>

namespace SPH
{
/** MPI message tags used by the halo exchange. */
constexpr int tag_halo_count = 7101;
constexpr int tag_halo_state = 7102;

/**
 * @class HaloVariablePacker
 * @brief Type-erased view of one particle variable that takes part in the exchange.
 *
 * The set of variables to exchange is case-dependent (a fluid needs Position, Velocity
 * and Density; a solid may need more), so it is registered at setup time by name rather
 * than hardcoded, mirroring BodyStatesRecordingToVtpCK::addToWrite.
 */
class HaloVariablePacker
{
  public:
    virtual ~HaloVariablePacker() {};
    /** Bytes occupied by one particle's worth of this variable. */
    virtual size_t bytesPerParticle() const = 0;
    /** Gather @p count particles listed in @p indices into @p buffer. */
    virtual void pack(const UnsignedInt *indices, UnsignedInt count, void *buffer) = 0;
    /** Scatter @p count particles from @p buffer into slots starting at @p first_slot. */
    virtual void unpack(UnsignedInt first_slot, UnsignedInt count, const void *buffer) = 0;
    virtual std::string Name() const = 0;
};

template <typename DataType>
class TypedHaloVariablePacker : public HaloVariablePacker
{
  public:
    explicit TypedHaloVariablePacker(DiscreteVariable<DataType> *dv_variable)
        : dv_variable_(dv_variable) {};

    size_t bytesPerParticle() const override { return sizeof(DataType); };
    void pack(const UnsignedInt *indices, UnsignedInt count, void *buffer) override;
    void unpack(UnsignedInt first_slot, UnsignedInt count, const void *buffer) override;
    std::string Name() const override { return dv_variable_->Name(); };

  protected:
    DiscreteVariable<DataType> *dv_variable_;
};

/**
 * @class HaloExchangeCK
 * @brief Per-body halo exchange across the neighbour ranks of a DomainPartition.
 */
template <class ExecutionPolicy>
class HaloExchangeCK
{
  public:
    HaloExchangeCK(RealBody &real_body, const DomainPartition &partition,
                   const MPIEnvironment &mpi_env, Real halo_reserve_factor = 0.25);
    virtual ~HaloExchangeCK() {};

    /** Register a variable to be shipped every step. Call at setup, before exec(). */
    template <typename DataType>
    void addToExchange(const std::string &variable_name);

    /**
     * Re-derive which owned particles each neighbour rank needs, and agree on counts.
     * Must be re-run whenever particles may have moved across the halo shell, i.e. at
     * the same cadence as the cell linked list rebuild.
     */
    void buildSendLists();

    /** Ship the registered variables for the particles in the send lists. */
    void exchangeState();

    /** buildSendLists() followed by exchangeState(). */
    void exec(Real dt = 0.0);

    UnsignedInt TotalHaloParticles() const { return total_halo_particles_; };
    UnsignedInt TotalOwnedParticles() const { return total_owned_particles_; };
    SingleVariable<UnsignedInt> *svTotalOwnedParticles() { return sv_total_owned_particles_; };

  protected:
    RealBody &real_body_;
    BaseParticles &particles_;
    const DomainPartition &partition_;
    const MPIEnvironment &mpi_env_;

    DiscreteVariable<Vecd> *dv_pos_;
    /** Loop bound for physics: the owned particles only, excluding the halo block. */
    SingleVariable<UnsignedInt> *sv_total_owned_particles_;

    UnsignedInt halo_capacity_ = 0;   /**< reserved halo slots */
    UnsignedInt halo_first_slot_ = 0; /**< index of the first halo slot */
    UnsignedInt total_halo_particles_ = 0;
    UnsignedInt total_owned_particles_ = 0;

    StdVec<int> neighbour_ranks_;
    /** Per neighbour: owned particle indices that neighbour needs. */
    StdVec<StdVec<UnsignedInt>> send_indices_;
    StdVec<UnsignedInt> send_counts_;
    StdVec<UnsignedInt> recv_counts_;
    /** Per neighbour: first halo slot holding that neighbour's particles. */
    StdVec<UnsignedInt> recv_first_slot_;

    StdVec<std::unique_ptr<HaloVariablePacker>> packers_;
    size_t bytes_per_particle_ = 0;
    StdVec<std::vector<char>> send_buffers_;
    StdVec<std::vector<char>> recv_buffers_;

    /** Exchange send_counts_ <-> recv_counts_ and lay out the halo block. */
    void negotiateCounts();
    void checkHaloCapacity() const;
};
} // namespace SPH
#endif // HALO_EXCHANGE_CK_H
