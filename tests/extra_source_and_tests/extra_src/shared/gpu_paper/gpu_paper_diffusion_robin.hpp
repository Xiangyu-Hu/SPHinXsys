/**
 * @file gpu_paper_diffusion_robin.hpp
 * @brief Implementations for GPU-paper Robin BC (A3).
 */
#ifndef SPHINXSYS_GPU_PAPER_DIFFUSION_ROBIN_HPP
#define SPHINXSYS_GPU_PAPER_DIFFUSION_ROBIN_HPP

#include "gpu_paper_diffusion_robin.h"

namespace SPH
{
//=================================================================================================//
template <class DiffusionType>
template <class DiffusionDynamics>
Robin<DiffusionType>::Robin(DiffusionDynamics &diffusion_dynamics, BaseParticles *contact_particles)
    : dv_n_(diffusion_dynamics.getBaseParticles().template getVariableByName<Vecd>("NormalDirection")),
      dv_contact_n_(contact_particles->template getVariableByName<Vecd>("NormalDirection")),
      dv_species_(diffusion_dynamics.getSpeciesByName(
          &diffusion_dynamics.getBaseParticles(), diffusion_dynamics.getSpeciesNames())),
      dv_contact_convection_(diffusion_dynamics.registerSpecies(
          contact_particles, diffusion_dynamics.getSpeciesNames(), "Convection")),
      sv_contact_species_infinity_(contact_particles->template registerSingleVariable<Real>(
          diffusion_dynamics.getSpeciesNames()[0] + "Infinity"))
{
    contact_particles->template addVariableToWrite<Real>(dv_contact_convection_);
}
//=================================================================================================//
template <class DiffusionType>
template <class ExecutionPolicy, class EncloserType>
Robin<DiffusionType>::ComputingKernel::ComputingKernel(const ExecutionPolicy &ex_policy,
                                                       EncloserType &encloser)
    : n_(encloser.dv_n_->DelegatedData(ex_policy)),
      contact_n_(encloser.dv_contact_n_->DelegatedData(ex_policy)),
      species_(encloser.dv_species_->DelegatedMultiEntryView(ex_policy)),
      contact_convection_(encloser.dv_contact_convection_->DelegatedMultiEntryView(ex_policy)),
      contact_species_infinity_(encloser.sv_contact_species_infinity_->DelegatedData(ex_policy)) {}
//=================================================================================================//
template <class DiffusionType>
Vecd Robin<DiffusionType>::ComputingKernel::operator()(UnsignedInt m, UnsignedInt index_i,
                                                       UnsignedInt index_j, const Vecd &e_ij,
                                                       const Vecd &vec_r_ij)
{
    Real phi_ij = *contact_species_infinity_ - species_[index_i][m];
    return contact_convection_[index_j][m] * phi_ij * (n_[index_i] - contact_n_[index_j]);
}
//=================================================================================================//
} // namespace SPH

#endif
