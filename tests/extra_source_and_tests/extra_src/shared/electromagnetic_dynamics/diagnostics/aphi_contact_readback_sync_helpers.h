#ifndef APHI_CONTACT_READBACK_SYNC_HELPERS_H
#define APHI_CONTACT_READBACK_SYNC_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_contact_gmres_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

namespace SPH
{
namespace electromagnetics
{
namespace test
{

static constexpr const char *AphiContactProbeScalarName = "ContactProbeScalar";
static constexpr const char *AphiContactReadNeighborProbeName = "ContactReadNeighborProbe";

class RegisterContactProbeVariablesCK : public LocalDynamics
{
  public:
    RegisterContactProbeVariablesCK(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          dv_probe_(particles_->registerStateVariable<Real>(AphiContactProbeScalarName)),
          dv_read_(particles_->registerStateVariable<Real>(AphiContactReadNeighborProbeName))
    {
    }

  protected:
    DiscreteVariable<Real> *dv_probe_;
    DiscreteVariable<Real> *dv_read_;
};

class WriteConstantContactProbeCK : public LocalDynamics
{
  public:
    WriteConstantContactProbeCK(SPHBody &sph_body, Real value)
        : LocalDynamics(sph_body),
          dv_probe_(particles_->template getVariableByName<Real>(AphiContactProbeScalarName)), value_(value)
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : probe_(encloser.dv_probe_->DelegatedData(ex_policy)), value_(encloser.value_)
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            probe_[index_i] = value_;
        }

      protected:
        Real *probe_;
        Real value_;
    };

  protected:
    DiscreteVariable<Real> *dv_probe_;
    Real value_;
};

class ZeroContactProbeCK : public LocalDynamics
{
  public:
    explicit ZeroContactProbeCK(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          dv_probe_(particles_->template getVariableByName<Real>(AphiContactProbeScalarName))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : probe_(encloser.dv_probe_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            probe_[index_i] = 0.0;
        }

      protected:
        Real *probe_;
    };

  protected:
    DiscreteVariable<Real> *dv_probe_;
};

template <typename... RelationTypes>
class ContactMaxNeighborProbeCK;

template <typename... Parameters>
class ContactMaxNeighborProbeCK<Contact<Parameters...>> : public Interaction<Contact<Parameters...>>
{
    using BaseInteraction = Interaction<Contact<Parameters...>>;

  public:
    explicit ContactMaxNeighborProbeCK(Contact<Parameters...> &contact_relation)
        : BaseInteraction(contact_relation),
          dv_read_(this->particles_->template getVariableByName<Real>(AphiContactReadNeighborProbeName))
    {
        for (auto *contact_particles : this->contact_particles_)
        {
            dv_contact_probe_.push_back(
                contact_particles->template getVariableByName<Real>(AphiContactProbeScalarName));
        }
    }

    explicit ContactMaxNeighborProbeCK(DynamicsArgs<Contact<Parameters...>> parameters)
        : ContactMaxNeighborProbeCK(parameters.identifier_){};

    virtual ~ContactMaxNeighborProbeCK() = default;

    class InteractKernel : public BaseInteraction::InteractKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser, size_t contact_index);

        void interact(size_t index_i, Real dt = 0.0);

      protected:
        Real *read_;
        Real *contact_probe_;
    };

  protected:
    DiscreteVariable<Real> *dv_read_;
    StdVec<DiscreteVariable<Real> *> dv_contact_probe_;
};

template <typename... Parameters>
template <class ExecutionPolicy, class EncloserType>
inline ContactMaxNeighborProbeCK<Contact<Parameters...>>::InteractKernel::InteractKernel(const ExecutionPolicy &ex_policy,
                                                                                EncloserType &encloser,
                                                                                size_t contact_index)
    : BaseInteraction::InteractKernel(ex_policy, encloser, contact_index),
      read_(encloser.dv_read_->DelegatedData(ex_policy)),
      contact_probe_(encloser.dv_contact_probe_[contact_index]->DelegatedData(ex_policy))
{
}

template <typename... Parameters>
inline void ContactMaxNeighborProbeCK<Contact<Parameters...>>::InteractKernel::interact(size_t index_i, Real dt)
{
    (void)dt;
    Real max_probe = read_[index_i];
    for (UnsignedInt n = this->FirstNeighbor(index_i); n != this->LastNeighbor(index_i); ++n)
    {
        const UnsignedInt index_j = this->neighbor_index_[n];
        max_probe = std::max(max_probe, contact_probe_[index_j]);
    }
    read_[index_i] = max_probe;
}

class ZeroContactReadNeighborProbeCK : public LocalDynamics
{
  public:
    explicit ZeroContactReadNeighborProbeCK(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          dv_read_(particles_->template getVariableByName<Real>(AphiContactReadNeighborProbeName))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : read_(encloser.dv_read_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            read_[index_i] = 0.0;
        }

      protected:
        Real *read_;
    };

  protected:
    DiscreteVariable<Real> *dv_read_;
};

inline Real maxContactReadValueOnLeftBody(AphiTwoBodyInterfaceCase &case_setup)
{
    BaseParticles &particles = case_setup.left_body.getBaseParticles();
    syncVariableToHost<Real>(particles, AphiContactReadNeighborProbeName);
    const Real *read_values = particles.getVariableDataByName<Real>(AphiContactReadNeighborProbeName);
    const size_t total_real_particles = particles.TotalRealParticles();
    Real max_value = 0.0;
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        max_value = std::max(max_value, read_values[i]);
    }
    return max_value;
}

inline Real readNeighborProbeAfterDeviceWrite(AphiTwoBodyInterfaceCase &case_setup, Real probe_value)
{
    StateDynamics<MainExecutionPolicy, ZeroContactReadNeighborProbeCK> zero_read(case_setup.left_body);
    StateDynamics<MainExecutionPolicy, WriteConstantContactProbeCK> write_probe(case_setup.right_body, probe_value);
    InteractionDynamicsCK<MainExecutionPolicy, ContactMaxNeighborProbeCK<Contact<>>> read_neighbor(
        DynamicsArgs(case_setup.left_contact()));
    zero_read.exec();
    write_probe.exec();
    case_setup.updateRelations();
    read_neighbor.exec();
    return maxContactReadValueOnLeftBody(case_setup);
}

inline Real readNeighborProbeAfterHostWrite(AphiTwoBodyInterfaceCase &case_setup, Real probe_value, bool sync_to_device)
{
    BaseParticles &right_particles = case_setup.right_body.getBaseParticles();
    const size_t total_real_particles = right_particles.TotalRealParticles();
    Real *probe = right_particles.getVariableDataByName<Real>(AphiContactProbeScalarName);
    for (size_t i = 0; i != total_real_particles; ++i)
    {
        probe[i] = probe_value;
    }
    if (sync_to_device)
    {
        syncVariableToDevice<Real>(right_particles, AphiContactProbeScalarName);
    }
    StateDynamics<MainExecutionPolicy, ZeroContactReadNeighborProbeCK> zero_read(case_setup.left_body);
    InteractionDynamicsCK<MainExecutionPolicy, ContactMaxNeighborProbeCK<Contact<>>> read_neighbor(
        DynamicsArgs(case_setup.left_contact()));
    zero_read.exec();
    case_setup.updateRelations();
    read_neighbor.exec();
    return maxContactReadValueOnLeftBody(case_setup);
}

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_READBACK_SYNC_HELPERS_H
