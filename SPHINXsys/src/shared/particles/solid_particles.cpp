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
	SolidParticles::SolidParticles(SPHBody &sph_body,
								   SharedPtr<Solid> shared_solid_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr)
		: BaseParticles(sph_body, shared_solid_ptr, particle_generator_ptr)
	{
		shared_solid_ptr->assignSolidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<Vecd>(n_, "NormalDirection");
		registerAVariable<Vecd>(n_0_, "InitialNormalDirection");
		registerAVariable<Matd>(B_, "CorrectionMatrix", Matd(1.0));
		//----------------------------------------------------------------------
		//		for FSI
		//----------------------------------------------------------------------
		registerAVariable<Vecd>(vel_ave_, "AverageVelocity");
		registerAVariable<Vecd>(dvel_dt_ave_, "AverageAcceleration");
		registerAVariable<Vecd>(force_from_fluid_, "ForceFromFluid");
		//----------------------------------------------------------------------
		//		For solid-solid contact
		//----------------------------------------------------------------------
		registerAVariable<Real>(contact_density_, "ContactDensity");
		registerAVariable<Vecd>(contact_force_, "ContactForce");
		//-----------------------------------------------------------------------------------------
		//		register sortable particle data before building up particle configuration
		//-----------------------------------------------------------------------------------------
		registerASortableVariable<Vecd>("Position");
		registerASortableVariable<Vecd>("InitialPosition");
		registerASortableVariable<Real>("Volume");
		//sorting particle once
		//DynamicCast<RealBody>(this, body)->sortParticleWithCellLinkedList();
	}
	//=============================================================================================//
	void SolidParticles::offsetInitialParticlePosition(Vecd offset)
	{
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			pos_n_[i] += offset;
			pos_0_[i] += offset;
		}
	}
	//=================================================================================================//
	void SolidParticles::initializeNormalDirectionFromBodyShape()
	{
		ComplexShape &body_shape = sph_body_->body_shape_;
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			Vecd normal_direction = body_shape.findNormalDirection(pos_n_[i]);
			n_[i] = normal_direction;
			n_0_[i] = normal_direction;
		}
	}
	//=================================================================================================//
	void SolidParticles::initializeNormalDirectionFromShapeAndOp(const std::string &shape_name)
	{
		ComplexShape &body_shape = sph_body_->body_shape_;
		ShapeAndOp *shape_and_op = body_shape.getShapeAndOpByName(shape_name);
		Real switch_sign = shape_and_op->second == ShapeBooleanOps::add ? 1.0 : -1.0;
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			Vecd normal_direction = switch_sign * shape_and_op->first->findNormalDirection(pos_n_[i]);
			n_[i] = normal_direction;
			n_0_[i] = normal_direction;
		}
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
	ElasticSolidParticles::ElasticSolidParticles(SPHBody &sph_body,
												 SharedPtr<ElasticSolid> shared_elastic_solid_ptr,
												 SharedPtr<ParticleGenerator> particle_generator_ptr)
		: SolidParticles(sph_body, shared_elastic_solid_ptr, particle_generator_ptr)
	{
		shared_elastic_solid_ptr->assignElasticSolidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<Matd>(F_, "DeformationGradient", Matd(1.0));
		registerAVariable<Matd>(dF_dt_, "DeformationRate");
		registerAVariable<Matd>(stress_PK1_, "FirstPiolaKirchhoffStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<Matd>(variables_to_restart_, "DeformationGradient");
		//----------------------------------------------------------------------
		//		add basic output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<Vecd>("NormalDirection");
		addDerivedVariableToWrite<Displacement>();
		addDerivedVariableToWrite<VonMisesStress>();
		addDerivedVariableToWrite<VonMisesStrain>();
		addAVariableNameToList<Matd>(variables_to_restart_, "DeformationGradient");
		// get which stress measure is relevant for the material
		stress_measure_ = shared_elastic_solid_ptr->getRelevantStressMeasureName();
	}
	//=================================================================================================//
	StdLargeVec<Real> ElasticSolidParticles::getVonMisesStrainVector(std::string strain_measure, Real poisson)
	{
		StdLargeVec<Real> strain_vector = {};
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			Real strain = 0.0;
			if (strain_measure == "static") {
				strain = von_Mises_strain_static(index_i);
			} else if (strain_measure == "dynamic") {
				strain = von_Mises_strain_dynamic(index_i, poisson);
			} else {
				throw std::runtime_error("getVonMisesStrainVector: wrong input");
			}
			strain_vector.push_back(strain);
		}
		return strain_vector;
	}
	//=================================================================================================//
	Real ElasticSolidParticles::getVonMisesStrainMax(std::string strain_measure, Real poisson)
	{
		Real strain_max = 0;
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			Real strain = 0.0;
			if (strain_measure == "static") {
				strain = von_Mises_strain_static(index_i);
			} else if (strain_measure == "dynamic") {
				strain = von_Mises_strain_dynamic(index_i, poisson);
			} else {
				throw std::runtime_error("getVonMisesStrainMax: wrong input");
			}
			if (strain_max < strain) strain_max = strain;
		}
		return strain_max;
	}
	//=================================================================================================//
	StdLargeVec<Real> ElasticSolidParticles::getVonMisesStressVector()
	{	
		StdLargeVec<Real> stress_vector = {};
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			Real stress = get_von_Mises_stress(index_i);
			stress_vector.push_back(stress);
		}
		return stress_vector;
	}
	//=================================================================================================//
	Real ElasticSolidParticles::getVonMisesStressMax()
	{
		Real stress_max = 0.0;
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			Real stress = get_von_Mises_stress(index_i);
			if (stress_max < stress) stress_max = stress;
		}
		return stress_max;
	}
	//=================================================================================================//
	Vecd ElasticSolidParticles::displacement(size_t particle_i)
	{
		return pos_n_[particle_i]-pos_0_[particle_i];
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
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			displacement_vector.push_back(displacement(index_i));
		}
		return displacement_vector;
	}
	//=================================================================================================//
	StdLargeVec<Vecd> ElasticSolidParticles::getNormal()
	{
		StdLargeVec<Vecd> normal_vector = {};
		for (size_t index_i = 0; index_i < pos_0_.size(); index_i++)
		{
			normal_vector.push_back(normal(index_i));
		}
		return normal_vector;
	}
	//=============================================================================================//
	void ActiveMuscleParticles::initializeActiveMuscleParticleData()
	{
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<Matd>(active_stress_, "ActiveStress");
		registerAVariable<Real>(active_contraction_stress_, "ActiveContractionStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<Real>(variables_to_restart_, "ActiveContractionStress");
	}
	//=============================================================================================//
	ShellParticles::ShellParticles(SPHBody &sph_body,
								   SharedPtr<ElasticSolid> shared_elastic_solid_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr, Real thickness)
		: ElasticSolidParticles(sph_body, shared_elastic_solid_ptr, particle_generator_ptr)
	{
		shared_elastic_solid_ptr->assignElasticSolidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<Matd>(transformation_matrix_, "TransformationMatrix", Matd(1.0));
		registerAVariable<Real>(shell_thickness_, "Thickness", thickness);
		registerAVariable<Vecd>(pseudo_n_, "PseudoNormal");
		registerAVariable<Vecd>(dpseudo_n_dt_, "PseudoNormalChangeRate");
		registerAVariable<Vecd>(dpseudo_n_d2t_, "PseudoNormal2ndOrderTimeDerivative");
		registerAVariable<Vecd>(rotation_, "Rotation");
		registerAVariable<Vecd>(angular_vel_, "AngularVelocity");
		registerAVariable<Vecd>(dangular_vel_dt_, "AngularAcceleration");
		registerAVariable<Matd>(F_bending_, "BendingDeformationGradient");
		registerAVariable<Matd>(dF_bending_dt_, "BendingDeformationGradientChangeRate");
		registerAVariable<Vecd>(global_shear_stress_, "GlobalShearStress");
		registerAVariable<Matd>(global_stress_, "GlobalStress");
		registerAVariable<Matd>(global_moment_, "GlobalMoment");
		//----------------------------------------------------------------------
		//		add basic output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<Vecd>("Rotation");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<Vecd>(variables_to_restart_, "PseudoNormal");
		addAVariableNameToList<Vecd>(variables_to_restart_, "Rotation");
		addAVariableNameToList<Vecd>(variables_to_restart_, "AngularVelocity");
	}
	//=================================================================================================//
}
