/* ------------------------------------------------------------------------- *
 *                                SPHinXsys                                  *
 * ------------------------------------------------------------------------- *
 * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
 * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
 * physical accurate simulation and aims to model coupled industrial dynamic *
 * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
 * (smoothed particle hydrodynamics), a meshless computational method using  *
 * particle discretization.                                                  *
 *                                                                           *
 * SPHinXsys is partially funded by German Research Foundation               *
 * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
 *  HU1527/12-1 and HU1527/12-4.                                             *
 *                                                                           *
 * Portions copyright (c) 2017-2025 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	weakly_compressible_fluid.h
 * @brief 	Describe the weakly compressible fluid which is used
 * 			model incompressible fluids. Here, we have included several equation of states.
 * 			Furthermore, A typical non-newtonian fluid model is included.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef WEAKLY_COMPRESSIBLE_FLUID_H
#define WEAKLY_COMPRESSIBLE_FLUID_H

#include "base_material.h"
#include "common_functors.h"
#include "sphinxsys_constant.h"

namespace SPH
{
class WeaklyCompressibleFluid : public Fluid
{
  protected:
    Real rho0_; /**< reference density. */
    Real c0_;   /**< reference sound speed. */
    Real p0_;   /**< reference pressure */

  public:
    explicit WeaklyCompressibleFluid(Real rho0, Real c0);
    explicit WeaklyCompressibleFluid(ConstructArgs<Real, Real> args);
    virtual ~WeaklyCompressibleFluid() {};
    virtual Real ReferenceDensity() const override { return rho0_; };
    virtual Real ReferenceSoundSpeed() const override { return c0_; };
    virtual Real getPressure(Real rho) override;
    virtual Real DensityFromPressure(Real p) override;
    virtual Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;

    class EosKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        EosKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser);
        Real PressureFromDensity(UnsignedInt, Real rho);
        Real DensityFromPressure(UnsignedInt, Real p);
        Real getSoundSpeed(UnsignedInt, Real, Real);
        Real getReferenceDensity(UnsignedInt);

      protected:
        Real rho0_, c0_, p0_;
    };
};

class WeaklyCompressibleMixture : public Fluid
{
  public:
    explicit WeaklyCompressibleMixture(Real c0);
    virtual ~WeaklyCompressibleMixture();
    virtual Real ReferenceSoundSpeed() const override { return c0_; };
    DiscreteVariable<Real> *dvReferenceDensity() const { return dv_rho0_; };
    // the following virtual functions are as they are deprecated for computing kernel
    // based implementations, but we keep them for backward compatibility
    [[deprecated("Use WeaklyCompressibleMultiSpecies::EosKernel::PressureFromDensity() instead.")]]
    Real getPressure(Real rho) override;
    [[deprecated("Use WeaklyCompressibleMultiSpecies::EosKernel::DensityFromPressure() instead.")]]
    Real DensityFromPressure(Real p) override;
    [[deprecated("Use WeaklyCompressibleMultiSpecies::EosKernel::getSoundSpeed() instead.")]]
    Real getSoundSpeed(Real p = 0.0, Real rho = 1.0) override;

    class EosKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        EosKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser);
        Real PressureFromDensity(UnsignedInt index_i, Real rho);
        Real DensityFromPressure(UnsignedInt index_i, Real p);
        Real getSoundSpeed(UnsignedInt, Real, Real);
        Real getReferenceDensity(UnsignedInt index_i);

      protected:
        Real c0_, c0_sq_;
        DataView<Real> rho0_;
    };

  protected:
    Real c0_;                         /**< reference sound speed. */
    DiscreteVariable<Real> *dv_rho0_; /**< local reference density. */
};

using NamesAndDensities = StdVec<std::pair<std::string, Real>>;
class WeaklyCompressibleMultiSpecies : public WeaklyCompressibleMixture
{
  protected:
    StdVec<std::string> species_name_list_; /**< species name list. */
    StdVec<Real> rho0_list_;                /**< reference density list. */
    ConstantArray<Real> *ca_inv_rho0_list_; /**< inverse reference density list. */
    DiscreteVariable<Real> *dv_Y_list_;     /**< species mass fraction list. */

  public:
    WeaklyCompressibleMultiSpecies(const NamesAndDensities &species_data, Real c0);
    WeaklyCompressibleMultiSpecies(const std::string &name, const NamesAndDensities &species_data, Real c0);
    virtual ~WeaklyCompressibleMultiSpecies();
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    virtual Real ReferenceDensity() const override { return rho0_list_[0]; };
    StdVec<std::string> getSpeciesNameList() const { return species_name_list_; };
    StdVec<Real> getReferenceDensityList() const { return rho0_list_; };
    DiscreteVariable<Real> *dvMassFraction() const { return dv_Y_list_; };
    ConstantArray<Real> *caInvReferenceDensity() const { return ca_inv_rho0_list_; };
    UnsignedInt NumberOfMixtures() const { return species_name_list_.size(); };

    class EosKernel : public WeaklyCompressibleMixture::EosKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        EosKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser);
        template <typename FractionType>
        void setMixtureFractions(UnsignedInt index_i, const FractionType &mass_fractions);
        Real computeReferenceDensity(UnsignedInt index_i);

      protected:
        ArrayView<Real> inv_rho0_list_;
        MultiEntryView<Real> Y_list_;
        DataView<Real> rho0_;
    };
};

class WeaklyCompressibleMultiPhase : public WeaklyCompressibleMixture
{
    UniquePtrsKeeper<Fluid> fluid_ptrs_;
    StdVec<WeaklyCompressibleFluid *> pure_phase_list_;                 /**< pure phase list. */
    StdVec<WeaklyCompressibleMultiSpecies *> multi_species_phase_list_; /**< multi-species phase list. */
    using PureEosKernel = typename WeaklyCompressibleFluid::EosKernel;
    using MultiSpeciesEosKernel = typename WeaklyCompressibleMultiSpecies::EosKernel;
    bool is_phases_set_ = false;

  protected:
    StdVec<std::string> phase_name_list_;      /**< phase name list. */
    DiscreteVariable<Real> *dv_phi_list_;      /**< phase volume fraction list. */
    DiscreteVariable<Vecd> *dv_velocity_list_; /**< phase velocity list. */
    ComputingKernelArray<WeaklyCompressibleFluid, PureEosKernel> *pure_eos_kernels_;
    ComputingKernelArray<WeaklyCompressibleMultiSpecies, MultiSpeciesEosKernel> *multi_species_eos_kernels_;

  public:
    WeaklyCompressibleMultiPhase(Real c0);
    virtual ~WeaklyCompressibleMultiPhase();
    void addPurePhases(const NamesAndDensities &pure_phases);
    void addMultiSpeciesPhases(const StdVec<std::pair<std::string, NamesAndDensities>> &multi_species_phases);
    void setPhases() { is_phases_set_ = true; };
    virtual void initializeLocalParameters(BaseParticles *base_particles) override;
    virtual Real ReferenceDensity() const override;
    StdVec<std::string> getPhaseNameList() const { return phase_name_list_; };
    DiscreteVariable<Real> *dvVolumeFraction() const { return dv_phi_list_; };
    UnsignedInt NumberOfMixtures() const;

    class EosKernel : public WeaklyCompressibleMixture::EosKernel
    {
      public:
        template <class ExecutionPolicy, class EnclosureType>
        EosKernel(const ExecutionPolicy &ex_policy, EnclosureType &encloser);
        template <typename FractionType>
        void setMixtureFractions(UnsignedInt index_i, const FractionType &volume_fractions);
        Real computeReferenceDensity(UnsignedInt index_i);
        template <typename VelocityType>
        Vecd computeMixtureVelocity(UnsignedInt index_i, const VelocityType &velocities);

      protected:
        ArrayView<WeaklyCompressibleFluid::EosKernel> pure_eos_;
        ArrayView<WeaklyCompressibleMultiSpecies::EosKernel> multi_species_eos_;
        MultiEntryView<Real> phi_list_;
        MultiEntryView<Vecd> velocity_list_;
    };
};
} // namespace SPH
#endif // WEAKLY_COMPRESSIBLE_FLUID_H