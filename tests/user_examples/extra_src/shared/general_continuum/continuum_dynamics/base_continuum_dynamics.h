#ifndef BASE_CONTINUUM_DYNAMICS_H
#define BASE_CONTINUUM_DYNAMICS_H

#include "base_fluid_dynamics.h"
#include "continuum_particles.h"

namespace SPH
{
    namespace continuum_dynamics
    {
        typedef DataDelegateContact<PlasticContinuumParticles, SolidParticles> FSIContactData;
        /**
         * @class InteractionWithWall
         * @brief Base class adding interaction with wall to general relaxation process
         */
        template <template <typename...> class BaseInteractionType>
        class InteractionWithWall : public BaseInteractionType<FSIContactData>
        {
        public:
            explicit InteractionWithWall(BaseContactRelation& wall_contact_relation)
                : BaseInteractionType<FSIContactData>(wall_contact_relation)
            {
                for (size_t k = 0; k != this->contact_particles_.size(); ++k)
                {
                    wall_vel_ave_.push_back(this->contact_particles_[k]->AverageVelocity());
                    wall_force_ave_.push_back(this->contact_particles_[k]->AverageForce());
                    wall_n_.push_back(&(this->contact_particles_[k]->n_));
                    wall_mass_.push_back(&(this->contact_particles_[k]->mass_));
                }
            };
            virtual ~InteractionWithWall() {};

        protected:
            StdVec<StdLargeVec<Vecd>*> wall_vel_ave_, wall_force_ave_, wall_n_;
            StdVec<StdLargeVec<Real>*> wall_mass_;
        };
    } // namespace continuum_dynamics
} // namespace SPH
#endif // BASE_CONTINUUM_DYNAMICS_H
