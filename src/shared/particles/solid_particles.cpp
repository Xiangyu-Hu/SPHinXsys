#include "solid_particles.h"
#include "solid_particles_variable.h"

#include "base_body.h"
#include "elastic_solid.h"
#include "inelastic_solid.h"
#include "xml_engine.h"

namespace SPH
{
//=============================================================================================//
SolidParticles::SolidParticles(SPHBody &sph_body, Solid *solid)
		: BaseParticles(sph_body, solid), solid_(*solid) {}
//=================================================================================================//
void SolidParticles::initializeOtherVariables()
{
		BaseParticles::initializeOtherVariables();
    registerVariable(pos0_, "InitialPosition", [&](size_t i) -> Vecd
                     { return pos_[i]; });
		registerVariable(n_, "NormalDirection");
    registerVariable(n0_, "InitialNormalDirection", [&](size_t i) -> Vecd
                     { return n_[i]; });
    registerVariable(B_, "CorrectionMatrix", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
}
//=============================================================================================//
ElasticSolidParticles::
		ElasticSolidParticles(SPHBody &sph_body, ElasticSolid *elastic_solid)
		: SolidParticles(sph_body, elastic_solid),
		  elastic_solid_(*elastic_solid) {}
//=================================================================================================//
void ElasticSolidParticles::initializeOtherVariables()
{
		SolidParticles::initializeOtherVariables();
		/**
		 *	register particle data
		 */
    registerVariable(F_, "DeformationGradient", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
		registerVariable(dF_dt_, "DeformationRate");
		/**
		 *	register FSI data
		 */
		registerVariable(vel_ave_, "AverageVelocity");
		registerVariable(acc_ave_, "AverageAcceleration");
		/**
		 *	add restart output particle data
		 */
		addVariableToRestart<Matd>("DeformationGradient");
		/**
		 *	add restart output particle data
     */
		addVariableToWrite<Vecd>("NormalDirection");
		addDerivedVariableToWrite<Displacement>();
		addDerivedVariableToWrite<VonMisesStress>();
		addDerivedVariableToWrite<VonMisesStrain>();
		addVariableToRestart<Matd>("DeformationGradient");
		// get which stress measure is relevant for the material
		stress_measure_ = elastic_solid_.getRelevantStressMeasureName();
}
//=================================================================================================//
Matd ElasticSolidParticles::getGreenLagrangeStrain(size_t particle_i)
{
		Matd F = F_[particle_i];
		return 0.5 * (F.transpose() * F - Matd::Identity());
}
//=================================================================================================//
Vecd ElasticSolidParticles::getPrincipalStrains(size_t particle_i)
{
    Matd epsilon = getGreenLagrangeStrain(particle_i);
		return getPrincipalValuesFromMatrix(epsilon);
}
//=================================================================================================//
Matd ElasticSolidParticles::getStressCauchy(size_t particle_i)
{
		Matd F = F_[particle_i];
		Matd stress_PK2 = elastic_solid_.StressPK2(F, particle_i);
		return (1.0 / F.determinant()) * F * stress_PK2 * F.transpose();
}
//=================================================================================================//
Matd ElasticSolidParticles::getStressPK2(size_t particle_i)
{
		return elastic_solid_.StressPK2(F_[particle_i], particle_i);
}
//=================================================================================================//
Vecd ElasticSolidParticles::getPrincipalStresses(size_t particle_i)
{
		Matd sigma;
		if (stress_measure_ == "Cauchy")
		{
			sigma = getStressCauchy(particle_i); // Cauchy stress
		}
		else if (stress_measure_ == "PK2")
		{
			sigma = getStressPK2(particle_i); // Second Piola-Kirchhoff stress
		}
		else
		{
			throw std::runtime_error("get_Principal_stresses: wrong input");
		}

		return getPrincipalValuesFromMatrix(sigma);
}
//=================================================================================================//
Real ElasticSolidParticles::getVonMisesStress(size_t particle_i)
{
		Matd sigma;
		if (stress_measure_ == "Cauchy")
		{
			sigma = getStressCauchy(particle_i); // Cauchy stress
		}
		else if (stress_measure_ == "PK2")
		{
			sigma = getStressPK2(particle_i); // Second Piola-Kirchhoff stress
		}
		else
		{
			throw std::runtime_error("get_von_Mises_stress: wrong input");
		}

		return getVonMisesStressFromMatrix(sigma);
}
//=================================================================================================//
StdLargeVec<Real> ElasticSolidParticles::getVonMisesStrainVector(std::string strain_measure)
{
		StdLargeVec<Real> strain_vector = {};
		for (size_t index_i = 0; index_i < pos0_.size(); index_i++)
		{
			Real strain = 0.0;
			if (strain_measure == "static")
			{
				strain = getVonMisesStrain(index_i);
			}
			else if (strain_measure == "dynamic")
			{
				strain = getVonMisesStrainDynamic(index_i, elastic_solid_.PoissonRatio());
			}
			else
			{
				throw std::runtime_error("getVonMisesStrainVector: wrong input");
			}

			strain_vector.push_back(strain);
		}

		return strain_vector;
}
//=================================================================================================//
Real ElasticSolidParticles::getVonMisesStrainMax(std::string strain_measure)
{
		Real strain_max = 0;
		for (size_t index_i = 0; index_i < pos0_.size(); index_i++)
		{
			Real strain = 0.0;
			if (strain_measure == "static")
			{
				strain = getVonMisesStrain(index_i);
			}
			else if (strain_measure == "dynamic")
			{
				strain = getVonMisesStrainDynamic(index_i, elastic_solid_.PoissonRatio());
			}
			else
			{
				throw std::runtime_error("getVonMisesStrainMax: wrong input");
			}
			if (strain_max < strain)
				strain_max = strain;
		}
		return strain_max;
}
//=================================================================================================//
Real ElasticSolidParticles::getPrincipalStressMax()
{
		Real stress_max = 0.0;
		for (size_t index_i = 0; index_i < pos0_.size(); index_i++)
		{
			/** take the max. component, which is the first one, this represents the max. tension. */
        Real stress = getPrincipalStresses(index_i)[0];
			if (stress_max < stress)
				stress_max = stress;
		}
		return stress_max;
};
//=================================================================================================//
StdLargeVec<Real> ElasticSolidParticles::getVonMisesStressVector()
{
		StdLargeVec<Real> stress_vector = {};
		for (size_t index_i = 0; index_i < pos0_.size(); index_i++)
		{
			Real stress = getVonMisesStress(index_i);
			stress_vector.push_back(stress);
		}
		return stress_vector;
}
//=================================================================================================//
Real ElasticSolidParticles::getVonMisesStressMax()
{
		Real stress_max = 0.0;
		for (size_t index_i = 0; index_i < pos0_.size(); index_i++)
		{
			Real stress = getVonMisesStress(index_i);
			if (stress_max < stress)
				stress_max = stress;
		}
		return stress_max;
}
//=================================================================================================//
Vecd ElasticSolidParticles::displacement(size_t particle_i)
{
		return pos_[particle_i] - pos0_[particle_i];
}
//=================================================================================================//
Vecd ElasticSolidParticles::normal(size_t particle_i)
{
		return n_[particle_i];
}
//=================================================================================================//
StdLargeVec<Vecd> ElasticSolidParticles::getDisplacement()
{
		StdLargeVec<Vecd> displacement_vector = {};
		for (size_t index_i = 0; index_i < pos0_.size(); index_i++)
		{
			displacement_vector.push_back(displacement(index_i));
		}
		return displacement_vector;
}
//=================================================================================================//
Real ElasticSolidParticles::getMaxDisplacement()
{
		Real displ_max = 0.0;
		for (size_t index_i = 0; index_i < pos0_.size(); index_i++)
		{
			Real displ = displacement(index_i).norm();
			if (displ_max < displ)
				displ_max = displ;
		}
		return displ_max;
};
//=================================================================================================//
StdLargeVec<Vecd> ElasticSolidParticles::getNormal()
{
		StdLargeVec<Vecd> normal_vector = {};
		for (size_t index_i = 0; index_i < pos0_.size(); index_i++)
		{
			normal_vector.push_back(normal(index_i));
		}
		return normal_vector;
}
//=============================================================================================//
ShellParticles::ShellParticles(SPHBody &sph_body, ElasticSolid *elastic_solid)
		: ElasticSolidParticles(sph_body, elastic_solid), thickness_ref_(1.0)
{
		//----------------------------------------------------------------------
		//		modify kernel function for surface particles
		//----------------------------------------------------------------------
		sph_body.sph_adaptation_->getKernel()->reduceOnce();
		//----------------------------------------------------------------------
		//		register geometric data only
		//----------------------------------------------------------------------
		registerVariable(n_, "NormalDirection");
		registerVariable(thickness_, "Thickness");
		/**
		 * add particle reload data
		 */
    addVariableToList<Vecd>(variables_to_reload_, "NormalDirection");
    addVariableToList<Real>(variables_to_reload_, "Thickness");
}
//=================================================================================================//
void ShellParticles::initializeOtherVariables()
{
    BaseParticles::initializeOtherVariables();
    /**
     * register particle data
     */
    registerVariable(pos0_, "InitialPosition",
                     [&](size_t i) -> Vecd
                     { return pos_[i]; });
    registerVariable(n0_, "InitialNormalDirection",
                     [&](size_t i) -> Vecd
                     { return n_[i]; });
    registerVariable(transformation_matrix_, "TransformationMatrix");
    registerVariable(B_, "CorrectionMatrix", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
    registerVariable(F_, "DeformationGradient", [&](size_t i) -> Matd
                     { return Matd::Identity(); });
    registerVariable(dF_dt_, "DeformationRate");
    registerVariable(pseudo_n_, "PseudoNormal",
                     [&](size_t i) -> Vecd
                     { return n_[i]; });
    registerVariable(dpseudo_n_dt_, "PseudoNormalChangeRate");
    registerVariable(dpseudo_n_d2t_, "PseudoNormal2ndOrderTimeDerivative");
    registerVariable(rotation_, "Rotation");
    registerVariable(angular_vel_, "AngularVelocity");
    registerVariable(dangular_vel_dt_, "AngularAcceleration");
    registerVariable(F_bending_, "BendingDeformationGradient");
    registerVariable(dF_bending_dt_, "BendingDeformationGradientChangeRate");
    registerVariable(global_shear_stress_, "GlobalShearStress");
    registerVariable(global_stress_, "GlobalStress");
    registerVariable(global_moment_, "GlobalMoment");
    registerVariable(mid_surface_cauchy_stress_, "MidSurfaceCauchyStress");
    registerVariable(numerical_damping_scaling_, "NemrticalDampingScaling_",
                     [&](size_t i) -> Matd
                     { return Matd::Identity() * sph_body_.sph_adaptation_->ReferenceSmoothingLength(); });
    /**
     * for FSI
     */
    registerVariable(vel_ave_, "AverageVelocity");
    registerVariable(acc_ave_, "AverageAcceleration");
    /**
     * for rotation.
     */
    addVariableToRestart<Matd>("DeformationGradient");
    addVariableToRestart<Vecd>("PseudoNormal");
    addVariableToRestart<Vecd>("Rotation");
    addVariableToRestart<Vecd>("AngularVelocity");
    /**
     * add basic output particle data
     */
    addVariableToWrite<Vecd>("NormalDirection");
    addDerivedVariableToWrite<Displacement>();
    addDerivedVariableToWrite<VonMisesStress>();
    addDerivedVariableToWrite<VonMisesStrain>();
    addVariableToRestart<Matd>("DeformationGradient");
    addVariableToWrite<Vecd>("Rotation");
    addDerivedVariableToWrite<MidSurfaceVonMisesStressofShells>();
    /**
     * initialize transformation matrix
     */
    for (size_t i = 0; i != real_particles_bound_; ++i)
    {
                        transformation_matrix_[i] = getTransformationMatrix(n_[i]);
                        numerical_damping_scaling_[i](Dimensions - 1, Dimensions - 1) =
                            thickness_[i] < sph_body_.sph_adaptation_->ReferenceSmoothingLength() ? thickness_[i] : sph_body_.sph_adaptation_->ReferenceSmoothingLength();
    }
}

        }
