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
 * @file    diffusion_dynamics.h
 * @brief   These are particle dynamics applicable for all type of particles.
 * @author  Chi Zhang, Chenxi Zhao and Xiangyu Hu
 */

#ifndef DIFFUSION_DYNAMICS_H
#define DIFFUSION_DYNAMICS_H

#include "general_diffusion_reaction_dynamics.h"
#include "contact_diffusivity.h"

namespace SPH
{
    /**
     * Kernel gradient operator, viz. first order or renormalized with correct matrix. 
    */
    class KernelGradientInner
    {
    public:
        explicit KernelGradientInner(BaseParticles *inner_particles){};

        Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return dW_ijV_j * e_ij;
        };
    };

    class CorrectedKernelGradientInner
    {
        StdLargeVec<Matd> &B_;

    public:
        explicit CorrectedKernelGradientInner(BaseParticles *inner_particles)
            : B_(*inner_particles->getVariableByName<Matd>("KernelCorrectionMatrix"))
        {};

        Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return 0.5 * dW_ijV_j * (B_[index_i] + B_[index_j]) * e_ij;
        };
    };

    class KernelGradientContact
    {
    public:
        KernelGradientContact(BaseParticles *inner_particles, BaseParticles *contact_particles){};

        Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return dW_ijV_j * e_ij;
        };
    };

    class CorrectedKernelGradientContact
    {
        StdLargeVec<Matd> &B_;
        StdLargeVec<Matd> &contact_B_;

    public:
        CorrectedKernelGradientContact(BaseParticles *inner_particles, BaseParticles *contact_particles)
            : B_(*inner_particles->getVariableByName<Matd>("KernelCorrectionMatrix"))
            , contact_B_(*contact_particles->getVariableByName<Matd>("KernelCorrectionMatrix"))
        {};

        Vecd operator()(size_t index_i, size_t index_j, Real dW_ijV_j, const Vecd &e_ij)
        {
            return 0.5 * dW_ijV_j * (B_[index_i] + contact_B_[index_j]) * e_ij;
        };
    };

    /**
     * @class GetDiffusionTimeStepSize
     * @brief Computing the time step size based on diffusion coefficient and particle smoothing length
     */
    template <class ParticlesType>
    class GetDiffusionTimeStepSize : public BaseDynamics<Real>
                                   , public DiffusionReactionSimpleData<ParticlesType>
    {
    public:
        explicit GetDiffusionTimeStepSize(SPHBody &sph_body);
        virtual ~GetDiffusionTimeStepSize(){};

        virtual Real exec(Real dt = 0.0) override { return diff_time_step_; };

    protected:
        Real diff_time_step_;
    };

    /**
     * @class BaseDiffusionRelaxation
     * @brief Base for compute the diffusion of all species
     */
    template <class ParticlesType>
    class BaseDiffusionRelaxation : public LocalDynamics
                                  , public DiffusionReactionSimpleData<ParticlesType>
    {
    protected:
        typedef typename ParticlesType::DiffusionReactionMaterial Material;
        Material &material_;
        StdVec<BaseDiffusion *> &all_diffusions_;
        StdVec<StdLargeVec<Real> *> &diffusion_species_;
        StdVec<StdLargeVec<Real> *> &gradient_species_;
        StdVec<StdLargeVec<Real> *> diffusion_dt_;
        StdVec<Real> all_thermal_capcity_reciprocal_;

    public:
        typedef ParticlesType InnerParticlesType;
        explicit BaseDiffusionRelaxation(SPHBody &sph_body);
        virtual ~BaseDiffusionRelaxation(){};
    
        StdVec<BaseDiffusion *> &AllDiffusions() { return material_.AllDiffusions(); };

        void initialization(size_t index_i, Real dt = 0.0);
        void update(size_t index_i, Real dt = 0.0);
    };

    /**
     * @class BaseThermalDiffusionRelaxation
     * @brief Base for compute the diffusion of all species in thermal dynamics.
     */
    template <class ParticlesType>
    class BaseThermalDiffusionRelaxation : public BaseDiffusionRelaxation<ParticlesType>
    {
    protected:
        StdLargeVec<Real> &rho_;

    public:
        explicit BaseThermalDiffusionRelaxation(SPHBody &sph_body);
        virtual ~BaseThermalDiffusionRelaxation(){};

        void update(size_t index_i, Real dt = 0.0);
    };

    /**
     * @class BaseDiffusionRelaxationInner
     * @brief Compute the diffusion relaxation process of all species
     */
    template <class ParticlesType, class KernelGradientType = KernelGradientInner, template<class> class BaseDiffusionRelaxationType = BaseDiffusionRelaxation>
    class BaseDiffusionRelaxationInner : public BaseDiffusionRelaxationType<ParticlesType>
                                   , public DataDelegateInner<ParticlesType, DataDelegateEmptyBase>
    {
    protected:
        KernelGradientType kernel_gradient_;
        void getDiffusionChangeRate(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij);

    public:
        typedef BaseInnerRelation BodyRelationType;
        explicit BaseDiffusionRelaxationInner(BaseInnerRelation &inner_relation);
        virtual ~BaseDiffusionRelaxationInner(){};

        inline void interaction(size_t index_i, Real dt = 0.0);
    };

    template <class ParticlesType, class KernelGradientType = KernelGradientInner>
    using DiffusionRelaxationInner = BaseDiffusionRelaxationInner<ParticlesType, KernelGradientType, BaseDiffusionRelaxation>;

    template <class ParticlesType, class KernelGradientType = KernelGradientInner>
    using ThermalDiffusionRelaxationInner = BaseDiffusionRelaxationInner<ParticlesType, KernelGradientType, BaseThermalDiffusionRelaxation>;
    /**
     * @class DiffusionRelaxationContact
     * @brief Base class for diffusion relaxation process between two contact bodies.
     */
    template <class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact, template<class> class BaseDiffusionRelaxationType = BaseDiffusionRelaxation>
    class BaseDiffusionRelaxationContact : public BaseDiffusionRelaxationType<ParticlesType>
                                         , public DataDelegateContact<ParticlesType, ContactParticlesType, DataDelegateEmptyBase>
    {
    protected:
        StdVec<StdVec<std::string>> contact_gradient_species_names_;
        StdVec<KernelGradientType> contact_kernel_gradients_;

    public:
        typedef BaseContactRelation BodyRelationType;

        explicit BaseDiffusionRelaxationContact(BaseContactRelation &contact_relation);
        virtual ~BaseDiffusionRelaxationContact(){};
    };

    /**
     * @class DiffusionRelaxationContact
     * @brief Contact diffusion relaxation with Dirichlet boundary condition.
     */
    template < class DiffusionType, class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact, template<class> class BaseDiffusionRelaxationType = BaseDiffusionRelaxation>
    class DiffusionRelaxationContact
        : public BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType, BaseDiffusionRelaxationType>
    {
        StdVec<StdVec<BaseDiffusion *>> all_contact_diffusions_;
        StdVec<StdVec< ContactDiffusivity<DiffusionType> >> contact_thermal_diffusivities_;
        StdVec<StdVec<StdLargeVec<Real> *>> contact_gradient_species_;

     protected: 
        void getDiffusionChangeRateMultimaterialContact(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij, size_t contact_k);

    public:
        explicit DiffusionRelaxationContact(BaseContactRelation &contact_relation);
        virtual ~DiffusionRelaxationContact(){};

        inline void interaction(size_t index_i, Real dt = 0.0);
    };

    template <class DiffusionType, class ParticlesType, class KernelGradientType = KernelGradientContact>
    using ThermalDiffusionRelaxationContact = DiffusionRelaxationContact<DiffusionType, ParticlesType, KernelGradientType, BaseThermalDiffusionRelaxation<ParticlesType>>;

    /**
     * @class DiffusionRelaxationDirichlet
     * @brief Contact diffusion relaxation with Dirichlet boundary condition.
     */
    template <class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact, template<class> class BaseDiffusionRelaxationType = BaseDiffusionRelaxation>
    class DiffusionRelaxationDirichlet
        : public BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType, BaseDiffusionRelaxationType>
    {
    protected:
        StdVec<StdVec<StdLargeVec<Real> *>> contact_gradient_species_;
        void getDiffusionChangeRateDirichlet(size_t particle_i, size_t particle_j, Vecd &e_ij, Real surface_area_ij, size_t k);

    public:
        explicit DiffusionRelaxationDirichlet(BaseContactRelation &contact_relation);
        virtual ~DiffusionRelaxationDirichlet(){};
        inline void interaction(size_t index_i, Real dt = 0.0);
    };

    template <class ParticlesType, class KernelGradientType = KernelGradientContact>
    using ThermalDiffusionRelaxationDirichlet = DiffusionRelaxationDirichlet<ParticlesType, KernelGradientType, BaseThermalDiffusionRelaxation<ParticlesType>>;

    /**
     * @class DiffusionRelaxationNeumann
     * @brief Contact diffusion relaxation with Neumann boundary condition.
     */
    template <class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact, template<class> class BaseDiffusionRelaxationType = BaseDiffusionRelaxation>
    class DiffusionRelaxationNeumann
        : public BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType, BaseDiffusionRelaxationType>
    {
        StdLargeVec<Vecd> &n_;
        StdVec<StdLargeVec<Real> *> contact_heat_flux_;
        StdVec<StdLargeVec<Vecd> *> contact_n_;

    protected:
        void getDiffusionChangeRateNeumann(size_t particle_i, size_t particle_j, Real surface_area_ij_Neumann, size_t k);

    public:
        explicit DiffusionRelaxationNeumann(BaseContactRelation &contact_relation);
        virtual ~DiffusionRelaxationNeumann(){};

        inline void interaction(size_t index_i, Real dt = 0.0);
    };

    template <class ParticlesType, class KernelGradientType = KernelGradientContact>
    using ThermalDiffusionRelaxationNeumann = DiffusionRelaxationNeumann<ParticlesType, KernelGradientType, BaseThermalDiffusionRelaxation<ParticlesType>>;

    /**
     * @class RelaxationOfAllDiffusionSpeciesRobinContact
     * @brief Contact diffusion relaxation with Robin boundary condition.
     */
    template <class ParticlesType, class ContactParticlesType, class KernelGradientType = KernelGradientContact, template<class> class BaseDiffusionRelaxationType = BaseDiffusionRelaxation>
    class DiffusionRelaxationRobin
        : public BaseDiffusionRelaxationContact<ParticlesType, ContactParticlesType, KernelGradientType, BaseDiffusionRelaxationType>
    {
        StdLargeVec<Vecd> &n_;
        StdVec<StdLargeVec<Real> *> contact_convection_;
        StdVec<Real *> contact_T_infinity_;
        StdVec<StdLargeVec<Vecd> *> contact_n_;

    protected:
        void getDiffusionChangeRateRobin(size_t particle_i, size_t particle_j, Real surface_area_ij_Robin, size_t k);

    public:
        explicit DiffusionRelaxationRobin(BaseContactRelation &contact_relation);
        virtual ~DiffusionRelaxationRobin(){};

        inline void interaction(size_t index_i, Real dt = 0.0);
    };

    template <class ParticlesType, class KernelGradientType = KernelGradientContact>
    using ThermalDiffusionRelaxationRobin = DiffusionRelaxationRobin<ParticlesType, KernelGradientType, BaseThermalDiffusionRelaxation<ParticlesType>>;

    /**
     * @class InitializationRK
     * @brief Initialization of a runge-kutta integration scheme.
     */
    template <class ParticlesType>
    class InitializationRK : public BaseDiffusionRelaxation<ParticlesType>
    {
    protected:
        StdVec<StdLargeVec<Real>> &diffusion_species_s_;

    public:
        InitializationRK(SPHBody &sph_body, StdVec<StdLargeVec<Real>> &diffusion_species_s);
        virtual ~InitializationRK(){};

        void update(size_t index_i, Real dt = 0.0);
    };

    /**
     * @class SecondStageRK2
     * @brief The second stage of the 2nd-order Runge-Kutta scheme.
     */
    template <class FirstStageType>
    class SecondStageRK2 : public FirstStageType
    {
    protected:
        StdVec<StdLargeVec<Real>> &diffusion_species_s_;

    public:
        template <typename... ContactArgsType>
        SecondStageRK2(typename FirstStageType::BodyRelationType &body_relation
                     , StdVec<StdLargeVec<Real>> &diffusion_species_s
                     , ContactArgsType &&...contact_args)
            : FirstStageType(body_relation, std::forward<ContactArgsType>(contact_args)...)
            , diffusion_species_s_(diffusion_species_s)
        {};
        virtual ~SecondStageRK2(){};

        void update(size_t index_i, Real dt = 0.0);
    };

    /**
     * @class DiffusionRelaxationRK2
     * @brief The 2nd-order runge-kutta integration scheme.
     *        A intermediate state for species is introduced here to achieve multi-step integration.
     */
    template <class FirstStageType>
    class DiffusionRelaxationRK2 : public BaseDynamics<void>
    {
    protected:
        StdVec<StdLargeVec<Real>> diffusion_species_s_;                         /**< Intermediate state */

        SimpleDynamics<InitializationRK<typename FirstStageType::InnerParticlesType>> rk2_initialization_;
        Dynamics1Level<FirstStageType> rk2_1st_stage_;
        Dynamics1Level<SecondStageRK2<FirstStageType>> rk2_2nd_stage_;
        StdVec<BaseDiffusion *> all_diffusions_;

    public:
        template <typename... ContactArgsType>
        explicit DiffusionRelaxationRK2(typename FirstStageType::BodyRelationType &body_relation, ContactArgsType &&...contact_args);
        virtual ~DiffusionRelaxationRK2(){};

        virtual void exec(Real dt = 0.0) override;
    };
}         // namespace SPH
#endif    // DIFFUSION_DYNAMICS_H
