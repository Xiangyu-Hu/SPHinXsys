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
 * @file 	  contact_diffusivity.h
 * @brief 	This is the collection of handling the thermal diffusivity, 
 *          viz. a = k / rho/ c (https://en.wikipedia.org/wiki/Thermal_diffusivity), 
 *          at the particle interface, where different 
 *          materials are involved. 
 * @author	Chi Zhang
 */

#ifndef CONTACT_DIFFUSIVITY_H
#define CONTACT_DIFFUSIVITY_H

#include "base_data_package.h"
#include "diffusion_reaction.h"

namespace SPH
{
    /**
     * @struct IsotropicContactDiffusivity
     * @brief  Isotropic thermal diffusivity between different materials. 
     */
    template<class DiffusionType>
    class ContactDiffusivity
    {
    public:
        ContactDiffusivity(DiffusionType &diffusion_i, DiffusionType &diffusion_j){};
        virtual ~ContactDiffusivity(){};
    };

    template<>
    class ContactDiffusivity<BaseDiffusion>
    {
    public:
        ContactDiffusivity(BaseDiffusion &diffusion_i, BaseDiffusion &diffusion_j){};
        virtual ~ContactDiffusivity(){};

        Real getContactDiffusivity(size_t particle_i, size_t particle_j, const Vecd &r_ij){return 0;};
    };

    template<>
    class ContactDiffusivity<IsotropicDiffusion>
    {
    public:
        ContactDiffusivity(IsotropicDiffusion &diffusion_i, IsotropicDiffusion &diffusion_j)
        {
            contact_diff_cf_ = 2.0 * diffusion_i.getReferenceDiffusivity()  * diffusion_j.getReferenceDiffusivity() 
                                    / (diffusion_i.getReferenceDiffusivity()  + diffusion_j.getReferenceDiffusivity()) ;
        };
        virtual ~ContactDiffusivity(){};

        Real getContactDiffusivity(size_t particle_i, size_t particle_j, const Vecd &r_ij)
        {
            return contact_diff_cf_;
        };

    private:
        Real contact_diff_cf_;   /**< contact diffusion coefficient 2 * diff_cf__i * diff_cf__j / (diff_cf__i + diff_cf__j). */
    };
    
    template<>
    class ContactDiffusivity<IsotropicThermalDiffusion>
    {
    public:
        ContactDiffusivity(IsotropicThermalDiffusion &diffusion_i, IsotropicThermalDiffusion &diffusion_j)
        {
            contact_diffusivity_ = 2.0 * diffusion_i.getThermalConductivity()  * diffusion_j.getThermalConductivity() 
                                    / ( diffusion_i.getThermalConductivity()  + diffusion_j.getThermalConductivity() ) ;
        };
        virtual ~ContactDiffusivity(){};

        Real getContactDiffusivity(size_t particle_i, size_t particle_j, const Vecd &r_ij)
        {
            return contact_diffusivity_;
        };
        
    private:
        Real contact_diffusivity_;   /**< contact thermal diffusivity 2 * k_i * k_j / (k_i + k_j) . */
    };
}       // namespace SPH
#endif  //CONTACT_DIFFUSIVITY
