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
 * Portions copyright (c) 2017-2023 Technical University of Munich and       *
 * the authors' affiliations.                                                *
 *                                                                           *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
 * not use this file except in compliance with the License. You may obtain a *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
 *                                                                           *
 * ------------------------------------------------------------------------- */
/**
 * @file 	electro_physiology.h
 * @brief 	In is file, we declaim the dynamics relevant to electrophysiology,
 * 			including diffusion, reaction and muscle activation.
 * @author	Chi Zhang and Xiangyu Hu
 */

#ifndef ELECTRO_PHYSIOLOGY_H
#define ELECTRO_PHYSIOLOGY_H

#include "all_diffusion_reaction_dynamics.h"
#include "solid_particles.h"

namespace SPH
{
class ElectroPhysiologyReaction : public BaseReactionModel<3>
{
  protected:
    Real k_a_;
    size_t voltage_;
    size_t gate_variable_;
    size_t active_contraction_stress_;

    virtual Real getProductionRateIonicCurrent(LocalSpecies &species) = 0;
    virtual Real getLossRateIonicCurrent(LocalSpecies &species) = 0;
    virtual Real getProductionRateGateVariable(LocalSpecies &species) = 0;
    virtual Real getLossRateGateVariable(LocalSpecies &species) = 0;
    virtual Real getProductionActiveContractionStress(LocalSpecies &species);
    virtual Real getLossRateActiveContractionStress(LocalSpecies &species);

  public:
    explicit ElectroPhysiologyReaction(Real k_a)
        : BaseReactionModel<3>({"Voltage", "GateVariable", "ActiveContractionStress"}),
          k_a_(k_a), voltage_(species_indexes_map_["Voltage"]),
          gate_variable_(species_indexes_map_["GateVariable"]),
          active_contraction_stress_(species_indexes_map_["ActiveContractionStress"])
    {
        reaction_model_ = "ElectroPhysiologyReaction";
        initializeElectroPhysiologyReaction();
    };
    virtual ~ElectroPhysiologyReaction(){};

    void initializeElectroPhysiologyReaction();
};

/**
 * @class AlievPanfilowModel
 * @brief The simplest Electrophysiology Reaction model,
 * which reduces the complex of array of ion currents to two variables that
 * describe excitation and recovery.
 */
class AlievPanfilowModel : public ElectroPhysiologyReaction
{
  protected:
    /** Parameters for two variable cell model. */
    Real k_, a_, b_, mu_1_, mu_2_, epsilon_, c_m_;

    virtual Real getProductionRateIonicCurrent(LocalSpecies &species) override;
    virtual Real getLossRateIonicCurrent(LocalSpecies &species) override;
    virtual Real getProductionRateGateVariable(LocalSpecies &species) override;
    virtual Real getLossRateGateVariable(LocalSpecies &species) override;

  public:
    explicit AlievPanfilowModel(Real k_a, Real c_m, Real k, Real a, Real b, Real mu_1, Real mu_2, Real epsilon)
        : ElectroPhysiologyReaction(k_a), k_(k), a_(a), b_(b), mu_1_(mu_1), mu_2_(mu_2),
          epsilon_(epsilon), c_m_(c_m)
    {
        reaction_model_ = "AlievPanfilowModel";
    };
    virtual ~AlievPanfilowModel(){};
};

// type trait for pass type template constructor
// This is a C++17 replacement to the C++20 https://en.cppreference.com/w/cpp/types/type_identity.
template <typename T>
struct TypeIdentity
{
};
/**
 * @class MonoFieldElectroPhysiology
 * @brief material class for electro_physiology.
 */
class MonoFieldElectroPhysiology : public DiffusionReaction<Solid, 3>
{
  public:
    template <class DiffusionType>
    MonoFieldElectroPhysiology(SharedPtr<ElectroPhysiologyReaction> electro_physiology_reaction_ptr,
                               TypeIdentity<DiffusionType> empty_object,
                               Real diff_cf, Real bias_diff_cf, Vecd bias_direction)
        : DiffusionReaction<Solid, 3>({"Voltage", "GateVariable", "ActiveContractionStress"},
                                      electro_physiology_reaction_ptr)
    {
        material_type_name_ = "MonoFieldElectroPhysiology";
        initializeAnDiffusion<DiffusionType>("Voltage", "Voltage", diff_cf, bias_diff_cf, bias_direction);
    };
    virtual ~MonoFieldElectroPhysiology(){};
};

/**
 * @class ElectroPhysiologyParticles
 * @brief A group of particles with electrophysiology particle data.
 */
class ElectroPhysiologyParticles
    : public DiffusionReactionParticles<SolidParticles, MonoFieldElectroPhysiology>
{
  public:
    ElectroPhysiologyParticles(SPHBody &sph_body, MonoFieldElectroPhysiology *mono_field_electro_physiology)
        : DiffusionReactionParticles<SolidParticles, MonoFieldElectroPhysiology>(sph_body, mono_field_electro_physiology){};
    virtual ~ElectroPhysiologyParticles(){};
    virtual ElectroPhysiologyParticles *ThisObjectPtr() override { return this; };
};

/**
 * @class ElectroPhysiologyReducedParticles
 * @brief A group of reduced particles with electrophysiology particle data.
 */
class ElectroPhysiologyReducedParticles : public ElectroPhysiologyParticles
{
  public:
    /** Constructor. */
    ElectroPhysiologyReducedParticles(SPHBody &sph_body, MonoFieldElectroPhysiology *mono_field_electro_physiology)
        : ElectroPhysiologyParticles(sph_body, mono_field_electro_physiology){};
    /** Destructor. */
    virtual ~ElectroPhysiologyReducedParticles(){};
    virtual ElectroPhysiologyReducedParticles *ThisObjectPtr() override { return this; };
};

namespace electro_physiology
{
typedef DiffusionReactionSimpleData<ElectroPhysiologyParticles> ElectroPhysiologyDataDelegateSimple;
typedef DiffusionReactionInnerData<ElectroPhysiologyParticles> ElectroPhysiologyDataDelegateInner;
/**
 * @class ElectroPhysiologyInitialCondition
 * @brief  set initial condition for a muscle body
 * This is a abstract class to be override for case specific initial conditions.
 */
class ElectroPhysiologyInitialCondition : public LocalDynamics,
                                          public ElectroPhysiologyDataDelegateSimple
{
  public:
    explicit ElectroPhysiologyInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          ElectroPhysiologyDataDelegateSimple(sph_body),
          pos_(particles_->pos_), all_species_(particles_->all_species_){};
    virtual ~ElectroPhysiologyInitialCondition(){};

  protected:
    StdLargeVec<Vecd> &pos_;
    StdVec<StdLargeVec<Real>> &all_species_;
};
/**
 * @class GetElectroPhysiologyTimeStepSize
 * @brief Computing the time step size from diffusion criteria
 */
class GetElectroPhysiologyTimeStepSize : public GetDiffusionTimeStepSize<ElectroPhysiologyParticles>
{
  public:
    explicit GetElectroPhysiologyTimeStepSize(RealBody &real_body)
        : GetDiffusionTimeStepSize<ElectroPhysiologyParticles>(real_body){};
    virtual ~GetElectroPhysiologyTimeStepSize(){};
};

using ElectroPhysiologyDiffusionRelaxationInner = DiffusionRelaxationInner<ElectroPhysiologyParticles, CorrectedKernelGradientInner>;
/**
 * @class ElectroPhysiologyDiffusionInnerRK2
 * @brief Compute the diffusion relaxation process
 */
class ElectroPhysiologyDiffusionInnerRK2
    : public DiffusionRelaxationRK2<ElectroPhysiologyDiffusionRelaxationInner>
{
  public:
    explicit ElectroPhysiologyDiffusionInnerRK2(BaseInnerRelation &inner_relation)
        : DiffusionRelaxationRK2(inner_relation){};
    virtual ~ElectroPhysiologyDiffusionInnerRK2(){};
};

using DiffusionRelaxationWithDirichletContact = DiffusionRelaxationDirichlet<ElectroPhysiologyParticles, ElectroPhysiologyParticles>;
/**
 * @class ElectroPhysiologyDiffusionRelaxationComplex
 * @brief Compute the diffusion relaxation process
 */
class ElectroPhysiologyDiffusionRelaxationComplex
    : public DiffusionRelaxationRK2<
          ComplexInteraction<ElectroPhysiologyDiffusionRelaxationInner, DiffusionRelaxationWithDirichletContact>>
{
  public:
    explicit ElectroPhysiologyDiffusionRelaxationComplex(BaseInnerRelation &inner_relation, BaseContactRelation &contact_relation)
        : DiffusionRelaxationRK2<ComplexInteraction<ElectroPhysiologyDiffusionRelaxationInner, DiffusionRelaxationWithDirichletContact>>(inner_relation, contact_relation){};
    virtual ~ElectroPhysiologyDiffusionRelaxationComplex(){};
};

/** Solve the reaction ODE equation of trans-membrane potential	using forward sweeping */
using ElectroPhysiologyReactionRelaxationForward =
    SimpleDynamics<ReactionRelaxationForward<ElectroPhysiologyParticles>>;
/** Solve the reaction ODE equation of trans-membrane potential	using backward sweeping */
using ElectroPhysiologyReactionRelaxationBackward =
    SimpleDynamics<ReactionRelaxationBackward<ElectroPhysiologyParticles>>;
} // namespace electro_physiology
} // namespace SPH
#endif // ELECTRO_PHYSIOLOGY_H