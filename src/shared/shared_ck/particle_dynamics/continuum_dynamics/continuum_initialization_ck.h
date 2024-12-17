// /* ------------------------------------------------------------------------- *
//  *                                SPHinXsys                                  *
//  * ------------------------------------------------------------------------- *
//  * SPHinXsys (pronunciation: s'finksis) is an acronym from Smoothed Particle *
//  * Hydrodynamics for industrial compleX systems. It provides C++ APIs for    *
//  * physical accurate simulation and aims to model coupled industrial dynamic *
//  * systems including fluid, solid, multi-body dynamics and beyond with SPH   *
//  * (smoothed particle hydrodynamics), a meshless computational method using  *
//  * particle discretization.                                                  *
//  *                                                                           *
//  * SPHinXsys is partially funded by German Research Foundation               *
//  * (Deutsche Forschungsgemeinschaft) DFG HU1527/6-1, HU1527/10-1,            *
//  *  HU1527/12-1 and HU1527/12-4.                                             *
//  *                                                                           *
//  * Portions copyright (c) 2017-2023 Technical University of Munich and       *
//  * the authors' affiliations.                                                *
//  *                                                                           *
//  * Licensed under the Apache License, Version 2.0 (the "License"); you may   *
//  * not use this file except in compliance with the License. You may obtain a *
//  * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.        *
//  *                                                                           *
//  * ------------------------------------------------------------------------- */
// /**
//  * @file 	continuum_integration.h
//  * @brief 	Here, we define the algorithm classes for continuum dynamics within the body.
//  * @details We consider here weakly compressible assumption to model elastic and
//  * 			plastic materials with the updated Lagrangian framework.
//  * @author	Shuaihao Zhang and Xiangyu Hu
//  */
// #ifndef CONTINUUM_INITIALIZATION_CK_H
// #define CONTINUUM_INITIALIZATION_CK_H

// #include "continuum_integration_1st_ck.h"

// namespace SPH
// {
// namespace continuum_dynamics
// {

// class StressDiffusionCK : public BasePlasticIntegration<DataDelegateInner>
// {
//   public:
//     explicit StressDiffusionCK(BaseInnerRelation &inner_relation);
//     virtual ~StressDiffusionCK(){};
//     void interaction(size_t index_i, Real dt = 0.0);

//   protected:
//     Real zeta_ = 0.1, phi_; /*diffusion coefficient*/
//     Real smoothing_length_, sound_speed_;
// };


// template <class RiemannSolverType, class KernelCorrectionType, typename... Parameters>
// class PlasticAcousticStep1stHalf<Inner<OneLevel, RiemannSolverType, KernelCorrectionType, Parameters...>>
//     : public PlasticAcousticStep<Interaction<Inner<Parameters...>>>
// {
//     using EosKernel = typename WeaklyCompressibleFluid::EosKernel;
//     using BaseInteraction = PlasticAcousticStep<Interaction<Inner<Parameters...>>>;

//   public:
//     explicit PlasticAcousticStep1stHalf(Relation<Inner<Parameters...>> &inner_relation);
//     virtual ~PlasticAcousticStep1stHalf(){};

//     class InitializeKernel
//     {
//       public:
//         template <class ExecutionPolicy, class EncloserType>
//         InitializeKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
//         void initialize(size_t index_i, Real dt = 0.0);

//       protected:
//         EosKernel eos_;
//         Real *rho_, *p_, *drho_dt_;
//         Vecd *vel_, *dpos_;
//         Mat3d *stress_tensor_3D_;
//     };

//     class InteractKernel : public BaseInteraction::InteractKernel
//     {
//       public:
//         template <class ExecutionPolicy, class EncloserType>
//         InteractKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
//         void interact(size_t index_i, Real dt = 0.0);

//       protected:
//         KernelCorrectionType correction_;
//         RiemannSolverType riemann_solver_;
//         Real *Vol_, *rho_, *p_, *drho_dt_, *mass_;
//         Vecd *force_;

//         //add
//         Mat3d *stress_tensor_3D_;
//     };

//     class UpdateKernel
//     {
//       public:
//         template <class ExecutionPolicy, class EncloserType>
//         UpdateKernel(const ExecutionPolicy &ex_policy, EncloserType &encloser);
//         void update(size_t index_i, Real dt = 0.0);

//       protected:
//         Real *mass_;
//         Vecd *vel_, *force_, *force_prior_;
//     };

//   protected:
//     KernelCorrectionType correction_;
//     RiemannSolverType riemann_solver_;
// };

// } // namespace continuum_dynamics
// } // namespace SPH
// #endif // CONTINUUM_INITIALIZATION_CK_H
