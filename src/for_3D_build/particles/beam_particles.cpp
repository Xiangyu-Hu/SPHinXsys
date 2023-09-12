#include "solid_particles.h"
#include "solid_particles_variable.h"

#include "base_body.h"
#include "elastic_solid.h"
#include "inelastic_solid.h"
#include "xml_engine.h"
#include "beam_particles.h"

namespace SPH
{
//=============================================================================================//



        BarParticles::BarParticles(SPHBody &sph_body, ElasticSolid *elastic_solid)
            : ShellParticles(sph_body, elastic_solid), width_ref_(1.0)
        {
                //----------------------------------------------------------------------
                //		modify kernel function for surface particles
                //----------------------------------------------------------------------
                // sph_body.sph_adaptation_->getKernel()->reduceTwice();
                //----------------------------------------------------------------------
                //		register geometric data only
                //----------------------------------------------------------------------
                registerVariable(b_n_, "BinormalDirection");
                registerVariable(width_, "Width");
                /**
                 * add particle reload data
                 */
                //addVariableNameToList<Vecd>(variables_to_reload_, "BinormalDirection");
                //addVariableNameToList<Real>(variables_to_reload_, "Width");
        }
        //=================================================================================================//
        void BarParticles::initializeOtherVariables()
        {
                ShellParticles::initializeOtherVariables();
                /**
                 * register particle data
                 */
                registerVariable(b_n0_, "InitialBinormalDirection",
                                 [&](size_t i) -> Vecd
                                 { return b_n_[i]; });
                registerVariable(pseudo_b_n_, "PseudoBinormal",
                                 [&](size_t i) -> Vecd
                                 { return pseudo_b_n_[i]; });
                registerVariable(dpseudo_b_n_dt_, "PseudoBinormalChangeRate");
                registerVariable(dpseudo_b_n_d2t_, "PseudoBinormal2ndOrderTimeDerivative");
                registerVariable(rotation_b_, "Rotation_b");
                registerVariable(angular_b_vel_, "AngularVelocity_b");
                registerVariable(dangular_b_vel_dt_, "AngularAccelerationofBinormal");
                registerVariable(F_b_bending_, "b_BendingDeformationGradient");
                registerVariable(dF_b_bending_dt_, "b_BendingDeformationGradientChangeRate");

                registerVariable(global_b_shear_stress_, "Global_b_ShearStress");
                registerVariable(global_b_stress_, "Global_b_Stress");
                registerVariable(global_b_moment_, "Global_b_Moment");

                /**
                 * for rotation.
                 */

                addVariableToRestart<Vecd>("PseudoBinormal");
                addVariableToRestart<Vecd>("Rotation_b");
                addVariableToRestart<Vecd>("AngularVelocity_b");
                /**
                 * add basic output particle data
                 */
                addVariableToWrite<Vecd>("BinormalDirection");
                addVariableToWrite<Vecd>("Rotation_b");
                /**
                 * initialize transformation matrix
                 */
                for (size_t i = 0; i != real_particles_bound_; ++i)
                {
        
                        transformation_matrix_[i] = getTransformationMatrix(n_[i], b_n_[i]);
                   
                }
        }
        }
