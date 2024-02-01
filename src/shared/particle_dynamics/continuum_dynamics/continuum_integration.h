#ifndef CONTINUUM_INTEGRATION_H
#define CONTINUUM_INTEGRATION_H
#include "constraint_dynamics.h"
#include "continuum_particles.h"
#include "fluid_integration.hpp"
#include "general_continuum.h"
#include "base_continuum_dynamics.h"
namespace SPH
{
    namespace continuum_dynamics
    {
        typedef DataDelegateSimple<ContinuumParticles> ContinuumDataSimple;
        typedef DataDelegateInner<ContinuumParticles> ContinuumDataInner;
        typedef DataDelegateSimple<PlasticContinuumParticles> PlasticContinuumDataSimple;
        typedef DataDelegateInner<PlasticContinuumParticles> PlasticContinuumDataInner;

        class ContinuumInitialCondition : public LocalDynamics, public PlasticContinuumDataSimple
        {
        public:
            explicit ContinuumInitialCondition(SPHBody& sph_body);
            virtual ~ContinuumInitialCondition() {};

        protected:
            StdLargeVec<Vecd>& pos_, & vel_;
            StdLargeVec<Mat3d>& stress_tensor_3D_;
        };

        template <class FluidDynamicsType>
        class BaseIntegration1stHalf : public FluidDynamicsType
        {
        public:
            explicit BaseIntegration1stHalf(BaseInnerRelation& inner_relation);
            virtual ~BaseIntegration1stHalf() {};
            void update(size_t index_i, Real dt = 0.0);

        protected:
            StdLargeVec<Vecd>& acc_shear_;
        };
        using Integration1stHalf = BaseIntegration1stHalf<fluid_dynamics::Integration1stHalfInnerNoRiemann>;
        using Integration1stHalfRiemann = BaseIntegration1stHalf<fluid_dynamics::Integration1stHalfInnerRiemann>;

        class ShearAccelerationRelaxation : public fluid_dynamics::BaseIntegration<ContinuumDataInner>
        {
        public:
            explicit ShearAccelerationRelaxation(BaseInnerRelation& inner_relation);
            virtual ~ShearAccelerationRelaxation() {};
            void interaction(size_t index_i, Real dt = 0.0);

        protected:
            GeneralContinuum& continuum_;
            Real G_, smoothing_length_;
            StdLargeVec<Matd>& shear_stress_;
            StdLargeVec<Vecd>& acc_shear_;
        };

        class ShearStressRelaxation : public fluid_dynamics::BaseIntegration<ContinuumDataInner>
        {
        public:
            explicit ShearStressRelaxation(BaseInnerRelation& inner_relation);
            virtual ~ShearStressRelaxation() {};
            void initialization(size_t index_i, Real dt = 0.0);
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

        protected:
            GeneralContinuum& continuum_;
            StdLargeVec<Matd>& shear_stress_, & shear_stress_rate_, & velocity_gradient_, & strain_tensor_, & strain_tensor_rate_;
            StdLargeVec<Real>& von_mises_stress_, & von_mises_strain_, & Vol_;
            StdLargeVec<Matd>& B_;
        };

        using FixBodyPartConstraint = solid_dynamics::FixConstraint<BodyPartByParticle, ContinuumDataSimple>;
        using FixedInAxisDirection = solid_dynamics::BaseFixedInAxisDirection<ContinuumDataSimple>;
        using ConstrainSolidBodyMassCenter = solid_dynamics::BaseConstrainSolidBodyMassCenter<ContinuumDataSimple>;

        template <class DataDelegationType>
        class BasePlasticIntegration : public fluid_dynamics::BaseIntegration<DataDelegationType>
        {
        public:
            template <class BaseRelationType>
            explicit BasePlasticIntegration(BaseRelationType& base_relation);
            virtual ~BasePlasticIntegration() {};
            Matd reduceTensor(Mat3d tensor_3d);
            Mat3d increaseTensor(Matd tensor_2d);
        protected:
            PlasticContinuum& plastic_continuum_;
            StdLargeVec<Mat3d>& stress_tensor_3D_, & strain_tensor_3D_, & stress_rate_3D_, & strain_rate_3D_;
            StdLargeVec<Mat3d>& elastic_strain_tensor_3D_, & elastic_strain_rate_3D_;
            StdLargeVec<Matd>& velocity_gradient_;
        };

        template <typename... InteractionTypes>
        class PlasticIntegration1stHalf;

        template <class RiemannSolverType>
        class PlasticIntegration1stHalf<Inner<>, RiemannSolverType>
            : public BasePlasticIntegration<PlasticContinuumDataInner>
        {
        public:
            explicit PlasticIntegration1stHalf(BaseInnerRelation& inner_relation);
            virtual ~PlasticIntegration1stHalf() {};
            void initialization(size_t index_i, Real dt = 0.0);
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);
            virtual Vecd computeNonConservativeForce(size_t index_i);

        protected:
            RiemannSolverType riemann_solver_;
        };
        using PlasticIntegration1stHalfInnerNoRiemann = PlasticIntegration1stHalf<Inner<>, NoRiemannSolver>;
        using PlasticIntegration1stHalfInnerRiemann = PlasticIntegration1stHalf<Inner<>, AcousticRiemannSolver>;

        using BaseIntegrationWithWall = InteractionWithWall<BasePlasticIntegration>;

        template <class RiemannSolverType>
        class PlasticIntegration1stHalf<Contact<Wall>, RiemannSolverType>
            : public BaseIntegrationWithWall
        {
        public:
            explicit PlasticIntegration1stHalf(BaseContactRelation& wall_contact_relation);
            virtual ~PlasticIntegration1stHalf() {};
            inline void interaction(size_t index_i, Real dt = 0.0);
            virtual Vecd computeNonConservativeForce(size_t index_i);

        protected:
            RiemannSolverType riemann_solver_;
        };
        
        template <class RiemannSolverType>
        using PlasticIntegration1stHalfWithWall = ComplexInteraction<PlasticIntegration1stHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;
        using PlasticIntegration1stHalfWithWallNoRiemann = PlasticIntegration1stHalfWithWall<NoRiemannSolver>;
        using PlasticIntegration1stHalfWithWallRiemann = PlasticIntegration1stHalfWithWall<AcousticRiemannSolver>;

        template <typename... InteractionTypes>
        class PlasticIntegration2ndHalf;

        template <class RiemannSolverType>
        class PlasticIntegration2ndHalf<Inner<>, RiemannSolverType>
            : public BasePlasticIntegration<PlasticContinuumDataInner>
        {
        public:
            explicit PlasticIntegration2ndHalf(BaseInnerRelation& inner_relation);
            virtual ~PlasticIntegration2ndHalf() {};
            void initialization(size_t index_i, Real dt = 0.0);
            void interaction(size_t index_i, Real dt = 0.0);
            void update(size_t index_i, Real dt = 0.0);

        protected:
            RiemannSolverType riemann_solver_;
            StdLargeVec<Real>& acc_deviatoric_plastic_strain_, & vertical_stress_;
            StdLargeVec<Real>& Vol_, & mass_;
            Real E_, nu_;
        };
        using PlasticIntegration2ndHalfInnerNoRiemann = PlasticIntegration2ndHalf<Inner<>, NoRiemannSolver>;
        using PlasticIntegration2ndHalfInnerRiemann = PlasticIntegration2ndHalf<Inner<>, AcousticRiemannSolver>;

        template <class RiemannSolverType>
        class PlasticIntegration2ndHalf<Contact<Wall>, RiemannSolverType>
            : public BaseIntegrationWithWall
        {
        public:
            explicit PlasticIntegration2ndHalf(BaseContactRelation& wall_contact_relation);
            virtual ~PlasticIntegration2ndHalf() {};
            inline void interaction(size_t index_i, Real dt = 0.0);

        protected:
            RiemannSolverType riemann_solver_;
        };

        template <class RiemannSolverType>
        using PlasticIntegration2ndHalfWithWall = ComplexInteraction<PlasticIntegration2ndHalf<Inner<>, Contact<Wall>>, RiemannSolverType>;
        using PlasticIntegration2ndHalfWithWallNoRiemann = PlasticIntegration2ndHalfWithWall<NoRiemannSolver>;
        using PlasticIntegration2ndHalfWithWallRiemann = PlasticIntegration2ndHalfWithWall<AcousticRiemannSolver>;

        class StressDiffusion : public BasePlasticIntegration<PlasticContinuumDataInner>
        {
        public:
            explicit StressDiffusion(BaseInnerRelation& inner_relation);
            virtual ~StressDiffusion() {};
            void interaction(size_t index_i, Real dt = 0.0);

        protected:
            Real zeta_ = 0.1, fai_; /*diffusion coefficient*/
            Real smoothing_length_, sound_speed_;
        };
    } // namespace continuum_dynamics
} // namespace SPH
#endif // CONTINUUM_INTEGRATION_H