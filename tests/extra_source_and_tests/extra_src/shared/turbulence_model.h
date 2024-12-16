#ifndef TURBULENCEMODEL_H
#define TURBULENCEMODEL_H
#include "fluid_integration.hpp"
#include "unstructured_mesh.h"
#include "all_particle_dynamics.h"
#include "fvm_ghost_boundary.h"
#include "rans_dynamics.h"
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
          SPHBody &bounds_;

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
        using KEpsilonStd1stHalfAcousticRiemannSolver = KEpsilonStd1stHalf<AcousticRiemannSolver>;
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
        };
        using KEpsilonStd2ndHalfExtendedHLLCRiemannSolver = KEpsilonStd2ndHalf<ExtendedHLLCRiemannSolver>;
        using KEpsilonStd2dnHalfAcousticRiemannSolver = KEpsilonStd2ndHalf<AcousticRiemannSolver>;
        //=================================================================================================//
        
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // TURBULENCEMODEL_H
