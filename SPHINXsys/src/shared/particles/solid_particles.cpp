/**
 * @file solid_particles.cpp
 * @brief Definition of functions declared in solid_particles.h
 * @author	Xiangyu Hu and Chi Zhang
 */
#include "solid_particles.h"

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
		for (size_t i = 0; i != pos_n_.size(); ++i)
			pos_0_[i] = pos_n_[i];
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
		registerAVariable<indexMatrix, Matd>(F_, "DeformationGradient", Matd(1.0));
		registerAVariable<indexMatrix, Matd>(dF_dt_, "DeformationRate");
		registerAVariable<indexMatrix, Matd>(stress_PK1_, "FirstPiolaKirchhoffStress");
		//----------------------------------------------------------------------
		//		add restart output particle data
		//----------------------------------------------------------------------
		addAVariableNameToList<indexMatrix, Matd>(variables_to_restart_, "DeformationGradient");
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
	void ElasticSolidParticles::writeParticlesToVtpFile(std::ostream &output_file)
	{
		SolidParticles::writeParticlesToVtpFile(output_file);

		size_t total_real_particles = total_real_particles_;

		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i)
		{
			output_file << std::fixed << std::setprecision(9) << get_von_Mises_stress(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
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
	//=================================================================================================//
	void ElasticSolidParticles::writeParticlesToVtuFile(std::ostream &output_file)
	{
		SolidParticles::writeParticlesToVtuFile(output_file);

		size_t total_real_particles = total_real_particles_;

		//write von Mises stress
		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			output_file << std::fixed << std::setprecision(9) << get_von_Mises_stress(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write Displacement
		output_file << "    <DataArray Name=\"Displacement\" type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			Vecd displacement_vector = displacement(i);
			output_file << displacement_vector[0] << " " << displacement_vector[1] << " " << displacement_vector[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write Normal Vectors
		output_file << "    <DataArray Name=\"Normal Vector\" type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			Vecd normal_vector = normal(i);
			output_file << normal_vector[0] << " " << normal_vector[1] << " " << normal_vector[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		//write von Mises strain
		output_file << "    <DataArray Name=\"von Mises strain\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_real_particles; ++i) {
			output_file << std::fixed << std::setprecision(9) << von_Mises_strain_static(i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
	}
	//=================================================================================================//
	void ElasticSolidParticles::writeSurfaceParticlesToVtuFile(std::ostream& output_file, BodySurface& surface_particles)
	{
		SolidParticles::writeSurfaceParticlesToVtuFile(output_file, surface_particles);

		size_t total_surface_particles = surface_particles.body_part_particles_.size();

		/** write Min Principal stress */
		/** precision: 6 - higher precision because it depends on the E modulus */
		output_file << "    <DataArray Name=\"Principal stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			Vecd stress = get_Principal_stresses(particle_i);
			output_file << std::fixed << std::setprecision(6) << stress[0] << " "; // take the max. component, which is the first one, this represents the max. tension
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		/** write von Mises stress */
		/** precision: 6 - higher precision because it depends on the E modulus */
		output_file << "    <DataArray Name=\"von Mises stress\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			output_file << std::fixed << std::setprecision(6) << get_von_Mises_stress(particle_i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		/** write Min Principal strain */
		/** precision: 3 - 0.1% accuracy */
		output_file << "    <DataArray Name=\"Principal strain\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			Vecd strain = get_Principal_strains(particle_i);
			output_file << std::fixed << std::setprecision(3) << strain[0] << " "; // take the max. component, which is the first one, this represents the max. tension
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		/** write von Mises strain */
		/** precision: 3 - 0.1% accuracy */
		output_file << "    <DataArray Name=\"von Mises strain\" type=\"Float32\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			output_file << std::fixed << std::setprecision(3) << von_Mises_strain_static(particle_i) << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		/** write Displacement */
		/** precision: 3 - 0.1 mm accuracy */
		output_file << "    <DataArray Name=\"Displacement\" type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			Vecd displacement_vector = displacement(particle_i);
			output_file << std::fixed << std::setprecision(4) << displacement_vector[0] << " "
						<< std::fixed << std::setprecision(4) << displacement_vector[1] << " "
						<< std::fixed << std::setprecision(4) << displacement_vector[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";

		/** write Normal Vectors  */
		/*
		// removed for production
		output_file << "    <DataArray Name=\"Normal Vector\" type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
		output_file << "    ";
		for (size_t i = 0; i != total_surface_particles; ++i) {
			size_t particle_i = surface_particles.body_part_particles_[i];
			Vecd normal_vector = normal(particle_i);
			output_file << normal_vector[0] << " " << normal_vector[1] << " " << normal_vector[2] << " ";
		}
		output_file << std::endl;
		output_file << "    </DataArray>\n";
		*/
	}
	//=================================================================================================//
	void ElasticSolidParticles::writePltFileHeader(std::ofstream &output_file)
	{
		SolidParticles::writePltFileHeader(output_file);

		output_file << ",\" von Mises stress \"";
		output_file << ",\" Displacement \"";
		output_file << ",\" Normal Vectors \"";
		output_file << ",\" von Mises strain \"";
	}
	//=================================================================================================//
	void ElasticSolidParticles::writePltFileParticleData(std::ofstream &output_file, size_t index_i)
	{
		SolidParticles::writePltFileParticleData(output_file, index_i);
		
		output_file << get_von_Mises_stress(index_i) << " ";
		Vecd displacement_vector = displacement(index_i);
		output_file << displacement_vector[0] << " " << displacement_vector[1] << " " << displacement_vector[2] << " "
			<< index_i << " ";
		output_file << von_Mises_strain_static(index_i) << " ";
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
	ShellParticles::ShellParticles(SPHBody &sph_body,
								   SharedPtr<ElasticSolid> shared_elastic_solid_ptr,
								   SharedPtr<ParticleGenerator> particle_generator_ptr, Real thickness)
		: ElasticSolidParticles(sph_body, shared_elastic_solid_ptr, particle_generator_ptr)
	{
		shared_elastic_solid_ptr->assignElasticSolidParticles(this);
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
