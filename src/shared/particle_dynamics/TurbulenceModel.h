#ifndef TURBULENCEMODEL_H
#define TURBULENCEMODEL_H
#include "fluid_integration.hpp"
#include "unstructured_mesh.h"
#include "all_particle_dynamics.h"
#include "riemann_solver.h"
#include "rans_dynamics.h"
#include "fvm_ghost_boundary.h"
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
            Real *K_prod_p_, *K_prod_, *Eps_p_, *Eps_sum_, *K_adv_, *K_lap_;
            Real *Eps_adv_, *Eps_lap_, *Eps_prodscalar_, *Eps_scalar_;
            Real Cmu_, sigmak_;
            Real sigmaeps_, C1eps_, C2eps_;
            Real *Tau_wall_;
            Real *K_, *Eps_, *mu_t_;
            GhostCreationFromMesh& ghost_creator_;
        };
        //=================================================================================================//
        class WallAdjacentCells : public BaseTurbulence
        {
        public:
            explicit WallAdjacentCells(BaseInnerRelation& inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~WallAdjacentCells() {};
            

        protected:
          Vecd *wallnormal_;
          Real *walladjacentcellflag_, *yp_;
          StdLargeVec<Real> walladjacentindex_,  wallghostindex_;
          StdLargeVec<Vecd> walleij_;
          Real ymax_;
          SPHBody &bounds_;

           void walladjacentcellyp();
           void update(size_t index_i, Real dt);
        };
        //=================================================================================================//
        /*
        class WallAdjacentCells : public BaseTurbulence
        {
         public:
           explicit WallAdjacentCells(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator);
           virtual ~WallAdjacentCells(){};

         protected:
           StdLargeVec<Real> walladjacentindex_, &walladjacentcellflag_, wallghostindex_, &wallfacearea_;
           StdLargeVec<Vecd> walleij_;

           void walladjacentcellyp();
        };
        */ 
        //=================================================================================================//
        /*
        class StdWallFunctionFVM : public BaseTurbulence
        {
         public:
           explicit StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator);
           virtual ~StdWallFunctionFVM(){};
           void nearwallquantities(size_t index_i);

         protected:
           StdLargeVec<Real> &walladjacentcellflag_, &yp_;
           StdLargeVec<Vecd> &wallnormal_;
           Real vonkar_, E_;
        };*/ 
        //=================================================================================================//
        class StdWallFunctionFVM : public WallAdjacentCells
        {
         public:
           explicit StdWallFunctionFVM(BaseInnerRelation &inner_relation, GhostCreationFromMesh &ghost_creator);
           virtual ~StdWallFunctionFVM(){};
           void nearwallquantities(size_t index_i);

         protected:
           //StdLargeVec<Real> &walladjacentcellflag_, &yp_;
           //StdLargeVec<Vecd> &wallnormal_;
           Real vonkar_, E_;
           Real *ystar_;
           Matd *vel_gradient_mat_;
        };
        //=================================================================================================//

        class KEpsilonStd1stHalf : public StdWallFunctionFVM
        {
            public:
            explicit KEpsilonStd1stHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~KEpsilonStd1stHalf() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            Real *dK_dt_, *walladjacentcellflag_, *strain_rate_;
            Matd *vel_gradient_mat_;
            Real *dudx_, *dudy_, *dvdx_, *dvdy_;
        };
        //=================================================================================================//
        class KEpsilonStd2ndHalf : public BaseTurbulence
        {
            public:
            explicit KEpsilonStd2ndHalf(BaseInnerRelation &inner_relation, GhostCreationFromMesh& ghost_creator);
            virtual ~KEpsilonStd2ndHalf() {};
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

            protected:
            Real *dEps_dt_, *walladjacentcellflag_;
        };
        //=================================================================================================//
        
    }// namespace fluid_dynamics
    

}// namespace SPH
#endif // TURBULENCEMODEL_H
