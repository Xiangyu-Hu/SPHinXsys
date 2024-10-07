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
            Real *K_prod_p_, *K_prod_, *Eps_p_, *K_adv_, *K_lap_;
            Real *Eps_adv_, *Eps_lap_, *Eps_prodscalar_, *Eps_scalar_, *Tau_wall_;
            Real Cmu_, sigmak_;
            Real sigmaeps_, C1eps_, C2eps_, *K_, *Eps_, *mu_t_;
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
          Real *walladjacentcellflag_, *yp_, *cornercellflag_, *boundary_type_;
          StdLargeVec<Real> walladjacentindex_,  wallghostindex_;
          StdLargeVec<Vecd> walleij_;
          Real ymax_;
          SPHBody &bounds_;

           void walladjacentcellyp();
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
           Real *ystar_, *yp_;
           Vecd *wallnormal_;
           Matd *vel_gradient_mat_;
           Real vonkar_, E_;
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
            Real *dK_dt_, *walladjacentcellflag_, *strain_rate_, *cornercellflag_, *boundary_type_;
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
