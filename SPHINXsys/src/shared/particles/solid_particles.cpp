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
		addDerivedVariableToWrite<Displacement>();
		addDerivedVariableToWrite<VonMisesStress>();
		addDerivedVariableToWrite<VonMisesStrain>();
		addAVariableToRestart<Matd>("DeformationGradient");
	}
	//=================================================================================================//
	Real ElasticSolidParticles::getPrincipalStressMax()
	{
		Real stress_max = 0.0;
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			Real stress = get_Principal_stresses(index_i)[0]; // take the max. component, which is the first one, this represents the max. tension
			if (stress_max < stress) stress_max = stress;
		}
		return stress_max;
	};
	//=================================================================================================//
	Vecd ElasticSolidParticles::displacement(size_t particle_i)
	{
		return pos_n_[particle_i]-pos_0_[particle_i];
	}
	//=================================================================================================//
	StdLargeVec<Vecd> ElasticSolidParticles::getDisplacement()
	{
		StdLargeVec<Vecd> displacement_vector = {};
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			displacement_vector.push_back(displacement(index_i));
		}
		return displacement_vector;
	}
	//=================================================================================================//
	Real ElasticSolidParticles::getMaxDisplacement()
	{
		Real displ_max = 0.0;
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			Real displ = displacement(index_i).norm();
			if (displ_max < displ) displ_max = displ;
		}
		return displ_max;
	};
	//=============================================================================================//
	ShellParticles::ShellParticles(SPHBody &sph_body, ElasticSolid *elastic_solid)
		: ElasticSolidParticles(sph_body, elastic_solid), thickness_ref_(1.0)
	{
		//----------------------------------------------------------------------
		//		register geometric data only
		//----------------------------------------------------------------------
		registerAVariable(n_, "NormalDirection");
		registerAVariable(thickness_, "Thickness");
		//----------------------------------------------------------------------
		//		add particle reload data
		//----------------------------------------------------------------------
		addAVariableNameToList<Vecd>(variables_to_reload_, "NormalDirection");
		addAVariableNameToList<Real>(variables_to_reload_, "Thickness");
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
		registerAVariable(transformation_matrix_, "TransformationMatrix");
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
		//----------------------------------------------------------------------
		//		initialize transformation matrix
		//----------------------------------------------------------------------
		for (size_t i = 0; i != real_particles_bound_; ++i)
		{
			transformation_matrix_[i] = getTransformationMatrix(n_[i]);
		}
	}
	//=================================================================================================//
}
