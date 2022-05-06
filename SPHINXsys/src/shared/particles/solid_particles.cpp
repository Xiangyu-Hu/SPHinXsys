/**
 * @file solid_particles.cpp
 * @brief Definition of functions declared in solid_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 */
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
		: BaseParticles(sph_body, solid) {}
	//=================================================================================================//
	void SolidParticles::initializeOtherVariables()
	{
		BaseParticles::initializeOtherVariables();

		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(pos_0_, "InitialPosition", "Position");
		registerAVariable(n_, "NormalDirection");
		registerAVariable(n_0_, "InitialNormalDirection");
		registerAVariable(B_, "CorrectionMatrix", Matd(1.0));
		//----------------------------------------------------------------------
		//		for FSI
		//----------------------------------------------------------------------
		registerAVariable(vel_ave_, "AverageVelocity");
		registerAVariable(dvel_dt_ave_, "AverageAcceleration");
		registerAVariable(force_from_fluid_, "ForceFromFluid");
		//----------------------------------------------------------------------
		//		For solid-solid contact
		//----------------------------------------------------------------------
		registerAVariable(contact_density_, "ContactDensity");
		registerAVariable(contact_force_, "ContactForce");
	}
	//=================================================================================================//
	Vecd SolidParticles::normalizeKernelGradient(size_t particle_index_i, Vecd &kernel_gradient)
	{
		return B_[particle_index_i] * kernel_gradient;
	}
	//=================================================================================================//
	Vecd SolidParticles::getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd &e_ij)
	{
		return 0.5 * dW_ij * (B_[particle_index_i] + B_[particle_index_j]) * e_ij;
	}
	//=============================================================================================//
	ElasticSolidParticles::
		ElasticSolidParticles(SPHBody &sph_body, ElasticSolid *elastic_solid)
		: SolidParticles(sph_body, elastic_solid) {}
	//=================================================================================================//
	void ElasticSolidParticles::initializeOtherVariables()
	{
		SolidParticles::initializeOtherVariables();

		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(F_, "DeformationGradient", Matd(1.0));
		registerAVariable(dF_dt_, "DeformationRate");
		registerAVariable(stress_PK1_, "FirstPiolaKirchhoffStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableToRestart<Matd>("DeformationGradient");
		//----------------------------------------------------------------------
		//		add basic output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<Vecd>("NormalDirection");
		addDerivedVariableToWrite<Displacement>();
		addDerivedVariableToWrite<VonMisesStress>();
		addDerivedVariableToWrite<VonMisesStrain>();
		addAVariableToRestart<Matd>("DeformationGradient");
	}
	//=================================================================================================//
	StdLargeVec<Real> ElasticSolidParticles::getVonMisesStress()
	{
		StdLargeVec<Real> von_Mises_stress_vector = {};
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			von_Mises_stress_vector.push_back(von_Mises_stress(index_i));
		}
		return von_Mises_stress_vector;
	}
	//=================================================================================================//
	Real ElasticSolidParticles::getMaxVonMisesStress()
	{
		Real von_Mises_stress_max = 0;
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			Real von_Mises_stress_i = von_Mises_stress(index_i);
			if (von_Mises_stress_max < von_Mises_stress_i)
			{
				von_Mises_stress_max = von_Mises_stress_i;
			}
		}
		return von_Mises_stress_max;
	}
	//=================================================================================================//
	StdLargeVec<Real> ElasticSolidParticles::getVonMisesStrain()
	{
		StdLargeVec<Real> von_Mises_strain_vector = {};
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			von_Mises_strain_vector.push_back(von_Mises_strain(index_i));
		}
		return von_Mises_strain_vector;
	}
	//=================================================================================================//
	Real ElasticSolidParticles::getMaxVonMisesStrain()
	{
		Real von_Mises_strain_max = 0;
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			Real von_Mises_strain_i = von_Mises_strain(index_i);
			if (von_Mises_strain_max < von_Mises_strain_i)
			{
				von_Mises_strain_max = von_Mises_strain_i;
			}
		}
		return von_Mises_strain_max;
	}
	//=============================================================================================//
	ShellParticles::ShellParticles(SPHBody &sph_body, ElasticSolid *elastic_solid)
		: ElasticSolidParticles(sph_body, elastic_solid), thickness_ref_(1.0)
	{
		registerAVariable(n_, "NormalDirection");
		registerAVariable(thickness_, "Thickness");
		registerAVariable(transformation_matrix_, "TransformationMatrix");
	}
	//=================================================================================================//
	void ShellParticles::initializeOtherVariables()
	{
		BaseParticles::initializeOtherVariables();
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable(pos_0_, "InitialPosition", "Position");
		registerAVariable(n_0_, "InitialNormalDirection", "NormalDirection");
		registerAVariable(B_, "CorrectionMatrix", Matd(1.0));
		registerAVariable(F_, "DeformationGradient", Matd(1.0));
		registerAVariable(dF_dt_, "DeformationRate");
		registerAVariable(stress_PK1_, "FirstPiolaKirchhoffStress");
		registerAVariable(pseudo_n_, "PseudoNormal", "NormalDirection");
		registerAVariable(dpseudo_n_dt_, "PseudoNormalChangeRate");
		registerAVariable(dpseudo_n_d2t_, "PseudoNormal2ndOrderTimeDerivative");
		registerAVariable(rotation_, "Rotation");
		registerAVariable(angular_vel_, "AngularVelocity");
		registerAVariable(dangular_vel_dt_, "AngularAcceleration");
		registerAVariable(F_bending_, "BendingDeformationGradient");
		registerAVariable(dF_bending_dt_, "BendingDeformationGradientChangeRate");
		registerAVariable(global_shear_stress_, "GlobalShearStress");
		registerAVariable(global_stress_, "GlobalStress");
		registerAVariable(global_moment_, "GlobalMoment");
		//----------------------------------------------------------------------
		//		for FSI
		//----------------------------------------------------------------------
		registerAVariable(vel_ave_, "AverageVelocity");
		registerAVariable(dvel_dt_ave_, "AverageAcceleration");
		registerAVariable(force_from_fluid_, "ForceFromFluid");
		//----------------------------------------------------------------------
		//		For solid-solid contact
		//----------------------------------------------------------------------
		registerAVariable(contact_density_, "ContactDensity");
		registerAVariable(contact_force_, "ContactForce");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableToRestart<Matd>("DeformationGradient");
		addAVariableToRestart<Vecd>("PseudoNormal");
		addAVariableToRestart<Vecd>("Rotation");
		addAVariableToRestart<Vecd>("AngularVelocity");
		//----------------------------------------------------------------------
		//		add basic output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<Vecd>("NormalDirection");
		addDerivedVariableToWrite<Displacement>();
		addDerivedVariableToWrite<VonMisesStress>();
		addDerivedVariableToWrite<VonMisesStrain>();
		addAVariableToRestart<Matd>("DeformationGradient");
		addAVariableToWrite<Vecd>("Rotation");
	}
	//=================================================================================================//
}
