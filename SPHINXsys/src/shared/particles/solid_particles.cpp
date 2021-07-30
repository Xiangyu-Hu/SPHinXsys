/**
 * @file solid_particles.cpp
 * @brief Definition of functions declared in solid_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 */
#include "solid_particles.h"

#include "geometry.h"
#include "base_body.h"
#include "elastic_solid.h"
#include "inelastic_solid.h"
#include "xml_engine.h"

namespace SPH {
	//=============================================================================================//
	SolidParticles::SolidParticles(SPHBody* body)
		: SolidParticles(body, new Solid()){}
	//=============================================================================================//
	SolidParticles::SolidParticles(SPHBody* body, Solid* solid)
		: BaseParticles(body, solid)
	{
		solid->assignSolidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexVector, Vecd>(pos_0_, "InitialPosition");
		registerAVariable<indexVector, Vecd>(n_, "NormalDirection");
		registerAVariable<indexVector, Vecd>(n_0_, "InitialNormalDirection");
		registerAVariable<indexMatrix, Matd>(B_, "CorrectionMatrix", Matd(1.0));
		//----------------------------------------------------------------------
		//		for FSI
		//----------------------------------------------------------------------
		registerAVariable<indexVector, Vecd>(vel_ave_, "AverageVelocity");
		registerAVariable<indexVector, Vecd>(dvel_dt_ave_, "AverageAcceleration");
		registerAVariable<indexVector, Vecd>(force_from_fluid_, "ForceFromFluid");
		//----------------------------------------------------------------------
		//		For solid-solid contact
		//----------------------------------------------------------------------
		registerAVariable<indexScalar, Real>(contact_density_, "ContactDensity");
		registerAVariable<indexVector, Vecd>(contact_force_, "ContactForce");
		//-----------------------------------------------------------------------------------------
		//		register sortable particle data before building up particle configuration
		//-----------------------------------------------------------------------------------------
		registerASortableVariable<indexVector, Vecd>("Position");
		registerASortableVariable<indexVector, Vecd>("InitialPosition");
		registerASortableVariable<indexScalar, Real>("Volume");
		//set the initial value for initial particle position
		for (size_t i = 0; i != pos_n_.size(); ++i) pos_0_[i] =  pos_n_[i];
		//sorting particle once
		//dynamic_cast<RealBody*>(body)->sortParticleWithCellLinkedList();
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
	void SolidParticles::initializeNormalDirectionFromGeometry()
	{
		ComplexShape* body_shape = body_->body_shape_;
		for (size_t i = 0; i != total_real_particles_; ++i)
		{
			Vecd normal_direction = body_shape->findNormalDirection(pos_n_[i]);
			n_[i] = normal_direction;
			n_0_[i] = normal_direction;
		}
	}
	//=================================================================================================//	
	Vecd SolidParticles::normalizeKernelGradient(size_t particle_index_i, Vecd& kernel_gradient) 
	{
		return  B_[particle_index_i] * kernel_gradient;
	}
	//=================================================================================================//
	Vecd SolidParticles::getKernelGradient(size_t particle_index_i, size_t particle_index_j, Real dW_ij, Vecd& e_ij) 
	{
		return 0.5 * dW_ij * (B_[particle_index_i] + B_[particle_index_j]) * e_ij;
	}
	//=============================================================================================//
	ElasticSolidParticles::ElasticSolidParticles(SPHBody* body, ElasticSolid* elastic_solid)
		: SolidParticles(body, elastic_solid)
	{
		elastic_solid->assignElasticSolidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexMatrix, Matd>(F_, "DeformationGradient", Matd(1.0));
		registerAVariable<indexMatrix, Matd>(dF_dt_, "DeformationRate");
		registerAVariable<indexMatrix, Matd>(stress_PK1_, "FirstPiolaKirchhoffStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<indexMatrix, Matd>(variables_to_restart_, "DeformationGradient");
	}
	//=================================================================================================//
	void ElasticSolidParticles::writeParticlesToVtuFile(std::ofstream& output_file)
	{
		SolidParticles::writeParticlesToVtuFile(output_file);

		size_t total_real_particles = total_real_particles_;

		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			output_file << std::fixed << std::setprecision(9) << von_Mises_stress(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
	}
	//=================================================================================================//
	void ElasticSolidParticles::writePltFileHeader(std::ofstream& output_file)
	{
		SolidParticles::writePltFileHeader(output_file);

		output_file << ",\" von Mises stress \"";
	}
	//=================================================================================================//
	void ElasticSolidParticles::writePltFileParticleData(std::ofstream& output_file, size_t index_i)
	{
		SolidParticles::writePltFileParticleData(output_file, index_i);
		
		output_file << von_Mises_stress(index_i) << " ";
	}
	//=============================================================================================//
	void ActiveMuscleParticles::initializeActiveMuscleParticleData()
	{
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexMatrix, Matd>(active_stress_, "ActiveStress");
		registerAVariable<indexScalar, Real>(active_contraction_stress_, "ActiveContractionStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<indexScalar, Real>(variables_to_restart_, "ActiveContractionStress");
	}
	//=============================================================================================//
	ShellParticles::ShellParticles(SPHBody* body, ElasticSolid* elastic_solid, Real thickness)
		: ElasticSolidParticles(body, elastic_solid)
	{
		elastic_solid->assignElasticSolidParticles(this);
		//----------------------------------------------------------------------
		//		register particle data
		//----------------------------------------------------------------------
		registerAVariable<indexMatrix, Matd>(transformation_matrix_, "TransformationMatrix", Matd(1.0));
		registerAVariable<indexScalar, Real>(shell_thickness_, "Thickness", thickness);
		registerAVariable<indexVector, Vecd>(pseudo_n_, "PseudoNormal");
		registerAVariable<indexVector, Vecd>(dpseudo_n_dt_, "PseudoNormalChangeRate");
		registerAVariable<indexVector, Vecd>(dpseudo_n_d2t_, "PseudoNormal2ndOrderTimeDerivative");
		registerAVariable<indexVector, Vecd>(rotation_, "Rotation");
		registerAVariable<indexVector, Vecd>(angular_vel_, "AngularVelocity");
		registerAVariable<indexVector, Vecd>(dangular_vel_dt_, "AngularAcceleration");
		registerAVariable<indexMatrix, Matd>(F_bending_, "BendingDeformationGradient");
		registerAVariable<indexMatrix, Matd>(dF_bending_dt_, "BendingDeformationGradientChangeRate");
		registerAVariable<indexVector, Vecd>(global_shear_stress_, "GlobalShearStress");
		registerAVariable<indexMatrix, Matd>(global_stress_, "GlobalStress");
		registerAVariable<indexMatrix, Matd>(global_moment_, "GlobalMoment");
		//----------------------------------------------------------------------
		//		add basic output particle data
		//----------------------------------------------------------------------
		addAVariableToWrite<indexVector, Vecd>("Rotation");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "PseudoNormal");
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "Rotation");
		addAVariableNameToList<indexVector, Vecd>(variables_to_restart_, "AngularVelocity");
	}
	//=================================================================================================//
}
