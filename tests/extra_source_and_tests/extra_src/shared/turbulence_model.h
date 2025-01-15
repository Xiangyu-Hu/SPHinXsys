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
 * @file 	turbulence_model.h
 * @brief 	Here, we implement the std k-epsilon model with the std wall function .
 * @author	Yash Mandaokar, Feng Wang and Xiangyu Hu
 */

#ifndef TURBULENCEMODEL_H
#define TURBULENCEMODEL_H
#include "fluid_integration.hpp"
#include "unstructured_mesh.h"
#include "all_particle_dynamics.h"
#include "fvm_ghost_boundary.h"
#include "rans_turbulence_dynamics.h"
#include "extended_eulerian_riemann_solver.h"

namespace SPH
{
    namespace fluid_dynamics
    {
    //=================================================================================================//
        class BaseTurbulence : public BaseIntegration<DataDelegateInner>
        {
            public:
            explicit BaseTurbulence(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~BaseTurbulence() {};
            

            protected:
            Vecd *mom_, *dmom_dt_;
            Real *dmass_dt_;
            Real *K_prod_p_, *K_prod_, *Eps_p_, *K_adv_, *K_lap_;
            Real *Eps_adv_, *Eps_lap_, *Eps_prod_, *Eps_destruction_, *Tau_wall_;
            Real C_mu_, sigma_k_;
            Real sigma_eps_, C1_eps_, C2_eps_, *K_, *Eps_, *mu_t_;
            GhostCreationFromMesh& ghost_creator_;
            Viscosity &viscosity_;
        };
        //=================================================================================================//
        class WallAdjacentCells : public BaseTurbulence
        {
        public:
            explicit WallAdjacentCells(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~WallAdjacentCells() {};
            

        protected:
          Vecd *wall_normal_;
          Real *wall_adjacent_cell_flag_, *yp_;
          StdLargeVec<Real> wall_adjacent_index_,  wall_ghost_index_;
          StdLargeVec<Vecd> wall_eij_;
          Real ymax_;

           void walladjacentcellyp();
           void update(size_t index_i, Real dt);
        };  
        //=================================================================================================//
        class TurbuleceVariablesGradient : public BaseTurbulence
        {
         public:
           explicit TurbuleceVariablesGradient(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator);
           virtual ~TurbuleceVariablesGradient(){};

         protected:
           Vecd *K_grad_, *Eps_grad_;

           void update(size_t index_i, Real dt);
        };
        //=================================================================================================//
        class StdWallFunctionFVM : public BaseTurbulence
        {
         public:
           explicit StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator);
           virtual ~StdWallFunctionFVM(){};
           void nearwallquantities(size_t index_i);

         protected:
           Real *y_star_, *yp_;
           Vecd *wall_normal_;
           Matd *vel_gradient_mat_;
           Real von_kar_, E_;
        };
        //=================================================================================================//
        template <class RiemannSolverType>
        class KEpsilonStd1stHalf : public StdWallFunctionFVM
        {
            public:
           explicit KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator, Real limiter_parameter = 0.0);
            virtual ~KEpsilonStd1stHalf() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            Real *dK_dt_, *wall_adjacent_cell_flag_, *strain_rate_;
            Real *dudx_, *dudy_, *dvdx_, *dvdy_;
            RiemannSolverType riemann_solver_;
            Vecd *K_grad_, *Eps_grad_;
        };
        using KEpsilonStd1stHalfExtendedHLLCRiemannSolver = KEpsilonStd1stHalf<ExtendedHLLCRiemannSolver>;
        using KEpsilonStd1stHalfSecondOrderUpwind = KEpsilonStd1stHalf<SecondOrderUpwind>;
        //=================================================================================================//
        template <class RiemannSolverType>
        class KEpsilonStd2ndHalf : public BaseTurbulence
        {
            public:
            explicit KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator, Real limiter_parameter = 0.0);
            virtual ~KEpsilonStd2ndHalf() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            Real *dEps_dt_, *wall_adjacent_cell_flag_;
            RiemannSolverType riemann_solver_;
            Vecd *K_grad_, *Eps_grad_;
            Matd *vel_gradient_mat_;
        };
        using KEpsilonStd2ndHalfExtendedHLLCRiemannSolver = KEpsilonStd2ndHalf<ExtendedHLLCRiemannSolver>;
        using KEpsilonStd2ndHalfSecondOrderUpwind = KEpsilonStd2ndHalf<SecondOrderUpwind>;
        //=================================================================================================//
        
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // TURBULENCEMODEL_H
