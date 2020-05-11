/**
 * @file 	heart_reader.cpp
 * @brief 	This is the first example of reading a simple heart geometry form stl file
 * @author 	Chi Zhang and Xiangyu Hu
 * @version 0.2.0
 */
/** header file and namespace. */
#include "sphinxsys.h"
using namespace SPH;
#define PI 3.1415926
/** Set the file path to the stl file. */
std::string full_path_to_stl_file = "./input/heart-new.stl";
/** Paremeters and physical properties. */
Vec3d domain_lower_bound(-55.0, -75.0, -35.0);
Vec3d domain_upper_bound(35.0, 5.0, 35.0);			
/** reference particle spacing. */
Real dp_0 	= (domain_upper_bound[0] - domain_lower_bound[0]) / 45.0;	
/** Material properties. */
Real diffusion_coff = 1.0;
/** Define the geometry. */
Geometry *CreateHeart()
{
	Vecd translation(-53.5, -70.0, -32.5);
	Geometry *geometry_myocardium = new Geometry(full_path_to_stl_file, translation, 1.0);

	return geometry_myocardium;
}
/** Define the myheart body. */
class MyHeart : public SolidBody
{
public:
	MyHeart(SPHSystem &system, string body_name, int refinement_level, ParticlesGeneratorOps op)
		: SolidBody(system, body_name, refinement_level, op)
	{
		body_region_.add_geometry(CreateHeart(), RegionBooleanOps::add);
		body_region_.done_modeling();
		/** Initilize the background mesh, e.g., Level set data. */
		addBackgroundMesh();
	}
};
/**
 * Setup diffusion material properties for mapping the fiber direction
 */
class DiffusionMaterial
	: public DiffusionReactionMaterial<ElasticSolidParticles, LocallyOrthotropicMuscle>
{
public:
	DiffusionMaterial()
		: DiffusionReactionMaterial<ElasticSolidParticles, LocallyOrthotropicMuscle>()
	{
		insertASpecies("Phi");
		assignDerivedMaterialParameters();
		initializeDiffusion();
	}

	/** Initialize diffusion reaction material. */
	virtual void initializeDiffusion() override 
	{
		IsotropicDiffusion* phi_diffusion
			= new IsotropicDiffusion(species_indexes_map_["Phi"], species_indexes_map_["Phi"], diffusion_coff);
		species_diffusion_.push_back(phi_diffusion);
	};
};
/** Set diffusion relaxation. */
class DiffusionRelaxation
	: public RelaxationOfAllDifussionSpeciesRK2<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>
{
public:
	DiffusionRelaxation(SolidBody* body)
		: RelaxationOfAllDifussionSpeciesRK2<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>(body){};
	virtual ~DiffusionRelaxation() {};
};
/** Imposing diffusion boundary condition */
class DiffusionBCs
	: public DiffusionReactionConstraint<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>
{
protected:
	size_t phi_;
	virtual void ConstraintAParticle(size_t index_particle_i, Real dt = 0.0) override
	{
		BaseParticleData &base_particle_data_i 			= particles_->base_particle_data_[index_particle_i];
		DiffusionReactionData& diffusion_data_i 		= particles_->diffusion_reaction_data_[index_particle_i];

		Vecd dist_2_face = body_->mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);	
		Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);

		Vecd center_norm = base_particle_data_i.pos_n_ / (base_particle_data_i.pos_n_.norm() + 1.0e-15);

		Real angle = dot(face_norm, center_norm);
		if (angle >= 0.0) 
		{
				diffusion_data_i.species_n_[phi_] = 1.0;
		}
		else 
		{
				if(base_particle_data_i.pos_n_[1] < - body_->particle_spacing_)
					diffusion_data_i.species_n_[phi_] = 0.0;
		}
	};
public:
	DiffusionBCs(SolidBody* body, BodySurface* body_part)
		: DiffusionReactionConstraint<SolidBody, ElasticSolidParticles, BodySurface, LocallyOrthotropicMuscle>(body, body_part)
		{
			phi_ = material_->getSpeciesIndexMap()["Phi"];
		};
	virtual ~DiffusionBCs() {};
};
/** Compute FiberandSheet direction after diffuision */
class ComputeFiberandSheetDirections
	: public DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>
{
protected:
	size_t phi_;
	Real beta_epi_, beta_endo_;
	/** We define the centerline vector, which is parallel to the ventricular centerline and pointing  apex-to-base.*/
	Vecd center_line_; 
	virtual void Update(size_t index_particle_i, Real dt = 0.0) override
	{
			BaseParticleData &base_particle_data_i 			= particles_->base_particle_data_[index_particle_i];
			DiffusionReactionData& diffusion_data_i 		= particles_->diffusion_reaction_data_[index_particle_i];
			/**
			 * Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
			 * 		Present  doi.org/10.1016/j.cma.2016.05.031
			 */
			/** Probe the face norm from Levelset field. */
			Vecd dist_2_face = body_->mesh_background_->ProbeNormalDirection(base_particle_data_i.pos_n_);	
			Vecd face_norm = dist_2_face / (dist_2_face.norm() + 1.0e-15);
			Vecd center_norm = base_particle_data_i.pos_n_ / (base_particle_data_i.pos_n_.norm() + 1.0e-15);
			if (dot(face_norm, center_norm) <= 0.0) 
			{
				face_norm = -face_norm;
			}
			/** Compute the centerline's projection on the plane orthogonal to face norm. */
			Vecd circumferential_direction = getCrossProduct(center_line_, face_norm); 
			Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
			/** The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo */
			Real beta = (beta_epi_ - beta_endo_) * diffusion_data_i.species_n_[phi_] + beta_endo_;
			/** Compute the rotation matrix through Rodrigues rotation formulation. */
			Vecd f_0 = cos(beta) * cd_norm + sin(beta) * getCrossProduct(face_norm, cd_norm) + 
				       dot(face_norm, cd_norm) * (1.0 - cos(beta)) * face_norm;

			if (base_particle_data_i.pos_n_[1] < -body_->particle_spacing_) {
				material_->local_f0_[index_particle_i] = f_0 / (f_0.norm() + 1.0e-15);
				material_->local_s0_[index_particle_i] = face_norm;
			} 
			else {
				material_->local_f0_[index_particle_i] = Vecd(0);
				material_->local_s0_[index_particle_i] = Vecd(0);
			}
	};
public:
	ComputeFiberandSheetDirections(SolidBody* body)
		: DiffusionBasedMapping<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle>(body)
		{
			phi_ = material_->getSpeciesIndexMap()["Phi"];
			center_line_ = Vecd(0.0, 1.0, 0.0);
			beta_epi_  = -(70.0 / 180.0) * PI;
			beta_endo_ = (80.0 / 180.0) * PI;
		};
	virtual ~ComputeFiberandSheetDirections() {};
};
/**
 *  The main program
 */
int main()
{
	/** Setup the system. */
	SPHSystem system(domain_lower_bound, domain_upper_bound, dp_0);
	/** Creat a Heart body, corresponding material, particles and reaction model. */
	MyHeart *heart_body = new MyHeart(system, "SimpleHeart", 0, ParticlesGeneratorOps::lattice);
	DiffusionMaterial *diffusion_material = new DiffusionMaterial();
	DiffusionReactionParticles<ElasticSolidParticles, LocallyOrthotropicMuscle>	diffusion_particles(heart_body, diffusion_material);
	/** 
	 * Set body contact map
	 * The contact map gives the data conntections between the bodies
	 * basically the the range of bidies to build neighbor particle lists
	 */
	SPHBodyTopology body_topology = { { heart_body,{} } };
	system.SetBodyTopology(&body_topology);
	/** Setting up the simulation. */
	system.SetupSPHSimulation();
	/**
	 * @brief 	Methods used for updating data structure.
	 */
	 /** Update the cell linked list system. */
	ParticleDynamicsCellLinkedList		update_cell_list(heart_body);
	/** Update the configuration of bodies when neccessary. */
	ParticleDynamicsInnerConfiguration 		update_inner_configuration(heart_body);
	/** Random reset the relax solid particle position. */
	RandomizePartilePosition  			random_particles(heart_body);
	/**
	 * @brief 	Algorithms for particle relaxation.
	 */
	relax_dynamics::BodySurfaceBounding
		body_surface_bounding(heart_body, new NearBodySurface(heart_body));
	/** Compute the time step for physics relaxation. */
	relax_dynamics::GetTimeStepSize get_relax_timestep(heart_body);
	/** Physics relax algorith without contact interactions. */
	relax_dynamics::PhysicsRelaxationInner	relax_process(heart_body);
	/**
	 * Diffusion process.
	 */
	/** Time step for diffusion. */
	GetDiffusionTimeStepSize<SolidBody, ElasticSolidParticles, LocallyOrthotropicMuscle> get_time_step_size(heart_body);
	/** Diffusion process for diffusion body. */
	DiffusionRelaxation 			diffusion_relaxation(heart_body);
	/** Compute the fiber and sheet after diffusion. */
	ComputeFiberandSheetDirections compute_fiber_sheet(heart_body);
	/** 
	 * Output 
	 */
	In_Output in_output(system);
	WriteBodyStatesToVtu 		write_states(in_output, system.real_bodies_);
	/** Write the particle reload files. */
	WriteReloadParticle 		write_particle_reload_files(in_output, system.real_bodies_);
	/** Write mesh data. */
	WriteBodyMeshToPlt 	write_background_mesh(in_output, heart_body);
	/** Write material property to xml file. */
	WriteReloadMaterialProperty write_material_property(in_output, diffusion_material);
	/** Mesh output */
	write_background_mesh.WriteToFile(0.0);
	/**
	 * @brief 	Physics relaxation starts here.
	 */
	 /** Relax the elastic structure. */
	random_particles.parallel_exec(0.25);
	body_surface_bounding.parallel_exec();
	write_states.WriteToFile(0.0);
	/**
	 * From here the time stepping begines.
	 * Set the starting time.
	 */
	int ite = 0;
	int relax_step = 1000;
	int diffusion_step = 100;
	Real dt = 0.0;
	while (ite < relax_step)
	{

		relax_process.parallel_exec(dt);
		body_surface_bounding.parallel_exec();
		dt = get_relax_timestep.parallel_exec();
		ite++;

		update_cell_list.parallel_exec();
		update_inner_configuration.parallel_exec();
		if(ite % 100 ==0)
		{
			cout << fixed << setprecision(9) << "Relaxation steps N = " << ite << "\n";
			write_states.WriteToFile(Real(ite)*1.0e-4);
		}
	}

	BodySurface* surface_part = new BodySurface(heart_body);
	/** constraint boundary contidtion for diffusion. */
	DiffusionBCs impose_diffusion_bc(heart_body, surface_part);
	impose_diffusion_bc.parallel_exec();

	write_states.WriteToFile(Real(ite)*1.0e-4);

	dt = 0.0;
	while (ite <= diffusion_step + relax_step)
	{
		dt = get_time_step_size.parallel_exec();
		diffusion_relaxation.parallel_exec(dt);
		impose_diffusion_bc.parallel_exec();
		if (ite % 10 == 0)
		{
			cout << "Diffusion steps N=" << ite - relax_step << "	dt: " << dt << "\n";
			write_states.WriteToFile(Real(ite) * 1.0e-4);
		}
		ite++;
	}
	compute_fiber_sheet.exec();
	ite++;
	write_states.WriteToFile(Real(ite)*1.0e-4);
	compute_fiber_sheet.parallel_exec();
	write_material_property.WriteToFile(0);
	write_particle_reload_files.WriteToFile(0);

	return 0;
}
