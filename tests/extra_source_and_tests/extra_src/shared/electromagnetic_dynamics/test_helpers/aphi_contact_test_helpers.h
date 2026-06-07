#ifndef APHI_CONTACT_TEST_HELPERS_H
#define APHI_CONTACT_TEST_HELPERS_H

#include "electromagnetic_dynamics/test_helpers/aphi_lhs_test_helpers.h"
#include "electromagnetic_dynamics/test_helpers/aphi_test_device_sync.h"

#include <map>

namespace SPH
{
namespace electromagnetics
{
namespace test
{

struct AphiContactPositionKey
{
    Vecd position;
    bool operator<(const AphiContactPositionKey &other) const
    {
        for (int d = 0; d < 3; ++d)
        {
            if (position[d] + TinyReal < other.position[d])
                return true;
            if (position[d] > other.position[d] + TinyReal)
                return false;
        }
        return false;
    }
};

inline AphiContactPositionKey makeContactPositionKey(const Vecd &position)
{
    AphiContactPositionKey key;
    key.position = position;
    return key;
}

inline bool isAwayFromInterface(const Vecd &position, Real x_interface, Real dp_0)
{
    return std::abs(position[0] - x_interface) > dp_0 + TinyReal;
}

/** Particles whose scalar-grad / graddiv PC stencil should not reach the duplicated split-interface layer. */
inline bool isGradStencilSafeFromInterface(const Vecd &position, Real x_interface, Real dp_0)
{
    const Real support_margin = 3.5 * dp_0;
    return std::abs(position[0] - x_interface) > support_margin + TinyReal;
}

struct AphiContactBlockByPosition
{
    Vecd a_real = Vecd::Zero();
    Vecd a_imag = Vecd::Zero();
    Real phi_real = 0.0;
    Real phi_imag = 0.0;
};

using AphiBlockMapByPosition = std::map<AphiContactPositionKey, AphiContactBlockByPosition>;

inline void collectCoreBlockByPosition(AphiBlockMapByPosition &block_map, BaseParticles &particles,
                                       const AphiBlockNames &block, Real body_length, Real body_height, Real body_width,
                                       Real core_shell, Real x_interface, Real dp_0)
{
    syncAphiBlockToHost(particles, block);
    syncVariableToHost<Vecd>(particles, "Position");

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *a_real = particles.getVariableDataByName<Vecd>(block.a_real);
    const Vecd *a_imag = particles.getVariableDataByName<Vecd>(block.a_imag);
    const Real *phi_real = particles.getVariableDataByName<Real>(block.phi_real);
    const Real *phi_imag = particles.getVariableDataByName<Real>(block.phi_imag);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!isAwayFromInterface(positions[i], x_interface, dp_0))
        {
            continue;
        }

        AphiContactBlockByPosition value;
        value.a_real = a_real[i];
        value.a_imag = a_imag[i];
        value.phi_real = phi_real[i];
        value.phi_imag = phi_imag[i];
        block_map[makeContactPositionKey(positions[i])] = value;
    }
}

inline Real maxAbsBlockDifference(const AphiBlockMapByPosition &reference, const AphiBlockMapByPosition &candidate,
                                  size_t &matched_particles, size_t &missing_particles)
{
    matched_particles = 0;
    missing_particles = 0;
    Real max_diff = 0.0;

    for (const auto &entry : reference)
    {
        const auto it = candidate.find(entry.first);
        if (it == candidate.end())
        {
            missing_particles += 1;
            continue;
        }
        matched_particles += 1;
        const AphiContactBlockByPosition &lhs = entry.second;
        const AphiContactBlockByPosition &rhs = it->second;
        max_diff = std::max(max_diff, (lhs.a_real - rhs.a_real).norm());
        max_diff = std::max(max_diff, (lhs.a_imag - rhs.a_imag).norm());
        max_diff = std::max(max_diff, std::abs(lhs.phi_real - rhs.phi_real));
        max_diff = std::max(max_diff, std::abs(lhs.phi_imag - rhs.phi_imag));
    }

    return max_diff;
}

inline void collectCoreScalarByPosition(std::map<AphiContactPositionKey, Real> &value_map, BaseParticles &particles,
                                        const std::string &variable_name, Real body_length, Real body_height,
                                        Real body_width, Real core_shell, Real x_interface, Real dp_0)
{
    syncVariableToHost<Real>(particles, variable_name);
    syncVariableToHost<Vecd>(particles, "Position");

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *values = particles.getVariableDataByName<Real>(variable_name);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!isAwayFromInterface(positions[i], x_interface, dp_0))
        {
            continue;
        }
        value_map[makeContactPositionKey(positions[i])] = values[i];
    }
}

inline Real maxAbsScalarDifference(const std::map<AphiContactPositionKey, Real> &reference,
                                   const std::map<AphiContactPositionKey, Real> &candidate, size_t &matched_particles,
                                   size_t &missing_particles)
{
    matched_particles = 0;
    missing_particles = 0;
    Real max_diff = 0.0;

    for (const auto &entry : reference)
    {
        const auto it = candidate.find(entry.first);
        if (it == candidate.end())
        {
            missing_particles += 1;
            continue;
        }
        matched_particles += 1;
        max_diff = std::max(max_diff, std::abs(entry.second - it->second));
    }

    return max_diff;
}

inline void collectMatchedRealByPosition(std::map<AphiContactPositionKey, Real> &value_map, BaseParticles &particles,
                                           const std::string &variable_name, Real body_length, Real body_height,
                                           Real body_width, Real core_shell, Real x_interface, Real dp_0,
                                           bool require_core, bool require_away_from_interface)
{
    syncVariableToHost<Real>(particles, variable_name);
    syncVariableToHost<Vecd>(particles, "Position");

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Real *values = particles.getVariableDataByName<Real>(variable_name);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (require_core && !isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (require_away_from_interface && !isAwayFromInterface(positions[i], x_interface, dp_0))
        {
            continue;
        }
        value_map[makeContactPositionKey(positions[i])] = values[i];
    }
}

inline void collectCoreRealByPosition(std::map<AphiContactPositionKey, Real> &value_map, BaseParticles &particles,
                                      const std::string &variable_name, Real body_length, Real body_height,
                                      Real body_width, Real core_shell, Real x_interface, Real dp_0)
{
    collectMatchedRealByPosition(value_map, particles, variable_name, body_length, body_height, body_width, core_shell,
                                 x_interface, dp_0, true, true);
}

inline Real maxAbsVecdMapDifference(const std::map<AphiContactPositionKey, Vecd> &reference,
                                    const std::map<AphiContactPositionKey, Vecd> &candidate, size_t &matched_particles,
                                    size_t &missing_particles)
{
    matched_particles = 0;
    missing_particles = 0;
    Real max_diff = 0.0;

    for (const auto &entry : reference)
    {
        const auto it = candidate.find(entry.first);
        if (it == candidate.end())
        {
            missing_particles += 1;
            continue;
        }
        matched_particles += 1;
        max_diff = std::max(max_diff, (entry.second - it->second).norm());
    }

    return max_diff;
}

inline void collectCoreVecdByPosition(std::map<AphiContactPositionKey, Vecd> &value_map, BaseParticles &particles,
                                      const std::string &variable_name, Real body_length, Real body_height,
                                      Real body_width, Real core_shell, Real x_interface, Real dp_0)
{
    syncVariableToHost<Vecd>(particles, variable_name);
    syncVariableToHost<Vecd>(particles, "Position");

    const size_t total_real_particles = particles.TotalRealParticles();
    const Vecd *positions = particles.getVariableDataByName<Vecd>("Position");
    const Vecd *values = particles.getVariableDataByName<Vecd>(variable_name);

    for (size_t i = 0; i != total_real_particles; ++i)
    {
        if (!isCoreParticle(positions[i], body_length, body_height, body_width, core_shell))
        {
            continue;
        }
        if (!isAwayFromInterface(positions[i], x_interface, dp_0))
        {
            continue;
        }
        value_map[makeContactPositionKey(positions[i])] = values[i];
    }
}

class AphiHalfSpaceBoxShape : public ComplexShape
{
  public:
    AphiHalfSpaceBoxShape(const std::string &shape_name, const Vecd &center, const Vecd &halfsize) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(center), halfsize);
    }
};

class AssignSeparableAphiFieldsCK : public LocalDynamics
{
  public:
    explicit AssignSeparableAphiFieldsCK(SPHBody &sph_body, const AphiVariableNames &names)
        : LocalDynamics(sph_body),
          dv_position_(particles_->template getVariableByName<Vecd>("Position")),
          dv_a_real_(particles_->template getVariableByName<Vecd>(names.solution.a_real)),
          dv_a_imag_(particles_->template getVariableByName<Vecd>(names.solution.a_imag)),
          dv_phi_real_(particles_->template getVariableByName<Real>(names.solution.phi_real)),
          dv_phi_imag_(particles_->template getVariableByName<Real>(names.solution.phi_imag))
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : position_(encloser.dv_position_->DelegatedData(ex_policy)),
              a_real_(encloser.dv_a_real_->DelegatedData(ex_policy)),
              a_imag_(encloser.dv_a_imag_->DelegatedData(ex_policy)),
              phi_real_(encloser.dv_phi_real_->DelegatedData(ex_policy)),
              phi_imag_(encloser.dv_phi_imag_->DelegatedData(ex_policy))
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            const Vecd &position = position_[index_i];
            const Real x = position[0];
            const Real y = position[1];
            const Real z = position[2];
            const Real pi = Pi;

            a_real_[index_i] = Vecd(std::sin(pi * x), std::sin(pi * y), std::sin(pi * z));
            a_imag_[index_i] = Vecd(std::cos(pi * x), std::cos(pi * y), std::cos(pi * z));
            phi_real_[index_i] = std::sin(pi * x) * std::cos(pi * y);
            phi_imag_[index_i] = std::cos(pi * x) * std::sin(pi * z);
        }

      protected:
        Vecd *position_;
        Vecd *a_real_;
        Vecd *a_imag_;
        Real *phi_real_;
        Real *phi_imag_;
    };

  protected:
    DiscreteVariable<Vecd> *dv_position_;
    DiscreteVariable<Vecd> *dv_a_real_;
    DiscreteVariable<Vecd> *dv_a_imag_;
    DiscreteVariable<Real> *dv_phi_real_;
    DiscreteVariable<Real> *dv_phi_imag_;
};

class AssignConstantMaterialSigmaCK : public LocalDynamics
{
  public:
    AssignConstantMaterialSigmaCK(SPHBody &sph_body, Real sigma, const AphiMaterialNames &material_names)
        : LocalDynamics(sph_body),
          dv_sigma_(particles_->template getVariableByName<Real>(material_names.sigma)),
          sigma_(sigma)
    {
    }

    class UpdateKernel
    {
      public:
        template <class ExecutionPolicy, class EncloserType>
        UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser)
            : sigma_(encloser.dv_sigma_->DelegatedData(ex_policy)), sigma_value_(encloser.sigma_)
        {
        }

        void update(size_t index_i, Real dt = 0.0)
        {
            (void)dt;
            sigma_[index_i] = sigma_value_;
        }

      protected:
        Real *sigma_;
        Real sigma_value_;
    };

  protected:
    DiscreteVariable<Real> *dv_sigma_;
    Real sigma_;
};

} // namespace test
} // namespace electromagnetics
} // namespace SPH

#endif // APHI_CONTACT_TEST_HELPERS_H
