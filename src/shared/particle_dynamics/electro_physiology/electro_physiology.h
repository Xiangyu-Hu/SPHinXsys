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
    virtual ~ElectroPhysiologyReaction() {};

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
    virtual ~AlievPanfilowModel() {};
};

/**
 * @class MonoFieldElectroPhysiology
 * @brief material class for electro_physiology.
 */
template <class DirectionalDiffusionType>
class MonoFieldElectroPhysiology
    : public ReactionDiffusion<ElectroPhysiologyReaction, DirectionalDiffusionType>
{
  public:
    template <typename... Args, size_t... Is>
    MonoFieldElectroPhysiology(ElectroPhysiologyReaction *electro_physiology_reaction,
                               ConstructArgs<Args...> args, std::index_sequence<Is...>)
        : ReactionDiffusion<ElectroPhysiologyReaction, DirectionalDiffusionType>(electro_physiology_reaction)
    {
        this->addDiffusion("Voltage", "Voltage", std::get<Is>(args)...);
    };
    template <class ElectroPhysiologyReactionType, typename... OtherArgs>
    explicit MonoFieldElectroPhysiology(ConstructArgs<ElectroPhysiologyReactionType *, ConstructArgs<OtherArgs...>> args)
        : MonoFieldElectroPhysiology(std::get<0>(args), std::get<1>(args), std::index_sequence_for<OtherArgs...>{}){};
    virtual ~MonoFieldElectroPhysiology() {};
};

namespace electro_physiology
{
template <class DirectionalDiffusionType>
using ElectroPhysiologyDiffusionRelaxationInner =
    DiffusionRelaxation<Inner<CorrectedKernelGradientInner>, DirectionalDiffusionType>;
/**
 * @class ElectroPhysiologyDiffusionInnerRK2
 * @brief Compute the diffusion relaxation process
 */
template <class DirectionalDiffusionType>
using ElectroPhysiologyDiffusionInnerRK2 =
    DiffusionRelaxationRK2<ElectroPhysiologyDiffusionRelaxationInner<DirectionalDiffusionType>>;

/**
 * @class ElectroPhysiologyDiffusionNetworkRK2
 * @brief Compute the diffusion relaxation process on network
 */
using ElectroPhysiologyDiffusionNetworkRK2 =
    DiffusionRelaxationRK2<DiffusionRelaxation<Inner<KernelGradientInner>, IsotropicDiffusion>>;

using DiffusionRelaxationWithDirichletContact =
    DiffusionRelaxation<Dirichlet<KernelGradientContact>, IsotropicDiffusion>;

template <class DirectionalDiffusionType, template <typename...> typename... ContactInteractionTypes>
using ElectroPhysiologyDiffusionRelaxationComplex =
    DiffusionBodyRelaxationComplex<DirectionalDiffusionType,
                                   KernelGradientInner, KernelGradientContact,
                                   ContactInteractionTypes...>;

/** Solve the reaction ODE equation of trans-membrane potential	using forward sweeping */
using ElectroPhysiologyReactionRelaxationForward =
    SimpleDynamics<ReactionRelaxationForward<ElectroPhysiologyReaction>>;
/** Solve the reaction ODE equation of trans-membrane potential	using backward sweeping */
using ElectroPhysiologyReactionRelaxationBackward =
    SimpleDynamics<ReactionRelaxationBackward<ElectroPhysiologyReaction>>;
} // namespace electro_physiology
} // namespace SPH
#endif // ELECTRO_PHYSIOLOGY_H