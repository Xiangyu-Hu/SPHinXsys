/**
 * @file 	3d_diffusion_with_anisotropic_resolution.cpp
 * @brief 	This is a test to validate our anisotropic diffusion algorithm with anisotropic resolution.
 * @author Xiaojing Tang and Xiangyu Hu 
 */
 
 
#include "anisotropic_diffusion_relaxation.h"
#include "base_particles.h"
#include "contact_body_relation.h"
#include "inner_body_relation.h"
#include "kernel_correction.h"
#include "observer_particles.h"
#include "anisotropic_diffusion_relaxation_3d.hpp"


#include "sphinxsys.h"  
using namespace SPH;   // Namespace cite here

//------------------------------------------------------------------------------
// global parameters for the case
//------------------------------------------------------------------------------
Real PL = 1.0;	/**< X-direction domain. */    
Real PH = 0.1; /**< Z-direction domain. */ 
Real BC = 0.3 * PL; 
/** Domain bounds of the system. */
Vec3d domain_lower_bound(-PL, -PL, -PH);
Vec3d domain_upper_bound(3.0 * PL, 3.0 * PL, 3.0 * PH);
BoundingBox system_domain_bounds(domain_lower_bound, domain_upper_bound);
int resolution(20);

int z_num = 6;
Real ratio_ = 4.0;

Vec3d scaling_vector = Vec3d(1.0, 1.0, 1.0 / ratio_);
Real scaling_factor = 1.0 / ratio_;

//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1.0;
Real rho0 = 1.0;


// reference particle spacing
Real resolution_ref = PH / z_num;
Real resolution_ref_large = ratio_ * resolution_ref;
int BP = 3; 
 
Real boundary_x = BP * resolution_ref_large;
Real boundary_z = BP * resolution_ref;

 
int x_num = PL / resolution_ref_large ;
int y_num = x_num;

Real V_j = resolution_ref_large * resolution_ref_large * resolution_ref;

//----------------------------------------------------------------------
//	Geometric shapes used in the case.
//----------------------------------------------------------------------
Vec3d halfsize_membrane(0.5 * PL , 0.5 * PL, 0.5 * PH);
Vec3d translation_membrane(0.5 * PL + boundary_x , 0.5 * PL + boundary_x, 0.5 * PH + boundary_z );


Vec3d halfsize_boundary(0.5 * PL + boundary_x, 0.5 * PL + boundary_x, 0.5 * PH + boundary_z );
Vec3d translation_boundary(0.5 * PL + boundary_x , 0.5 * PL + boundary_x, 0.5 * PH + boundary_z );

namespace SPH
{

class DiffusionBlock : public ComplexShape
{
public:
	explicit DiffusionBlock(const std::string &shape_name) : ComplexShape(shape_name)
	{
		add<TriangleMeshShapeBrick>(halfsize_membrane, resolution, translation_membrane);
	}
};

class DiffusionBoundary : public ComplexShape
{
public:
	explicit DiffusionBoundary(const std::string &shape_name) : ComplexShape(shape_name)
	{ 
		add<TriangleMeshShapeBrick>(halfsize_boundary, resolution, translation_boundary);
	}
};

 
//particle generator for domain particles
template <>
class ParticleGenerator<BaseParticles, DiffusionBlock>: public ParticleGenerator<BaseParticles>
{
  public:
   ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles){};

    virtual void prepareGeometricData() override
    {
        // set particles directly
		for (int i = 0; i < x_num; i++)
		{
			for (int j = 0; j < y_num; j++)
			{
				for (int k = 0; k < z_num; k++)
				{
					Real x = (i + 0.5 + BP) * resolution_ref_large;
					Real y = (j + 0.5 + BP) * resolution_ref_large;
					Real z = (k + 0.5 + BP) * resolution_ref;
					addPositionAndVolumetricMeasure(Vec3d(x, y, z),V_j);
				}
			}
		}
    }
};


//particle generator for boundary particles
template <>
class ParticleGenerator<BaseParticles, DiffusionBoundary>: public ParticleGenerator<BaseParticles>
{
  public:
    ParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles)
        : ParticleGenerator<BaseParticles>(sph_body, base_particles){};

    virtual void prepareGeometricData() override
    {
         // set front and back particles directly
		 for (int i = 0; i < BP; i++)
		 {
			 for (int j = 0; j < (y_num + 2*BP); j++)
			 {
				 for (int k = 0; k < (z_num + 2*BP ); k++)
				 {
					 Real x = (i + 0.5 ) * resolution_ref_large;
					 Real y = (j + 0.5 ) * resolution_ref_large;
					 Real z = (k + 0.5 ) * resolution_ref;
					 addPositionAndVolumetricMeasure(Vec3d(x, y, z), V_j);

					 Real x_top = (x_num + i + 0.5 + BP) * resolution_ref_large;
					addPositionAndVolumetricMeasure(Vec3d(x_top, y, z), V_j);
				 }
				 
			 }
		 }


   // set  left and right particles directly
	   for (int i = 0; i < x_num; i++)
	   {
		   for (int j = 0; j < BP; j++)
		   {
			   for (int k = 0; k < (z_num + 2*BP ); k++)
			   {
				   Real x = (i + 0.5 + 3.0) * resolution_ref_large;
				   Real y = (j + 0.5 ) * resolution_ref_large;
				   Real z = (k + 0.5) * resolution_ref;
				   addPositionAndVolumetricMeasure(Vec3d(x, y, z),V_j);

				   Real y_top= (y_num + j + 0.5 + BP) * resolution_ref_large;
				 
				   addPositionAndVolumetricMeasure(Vec3d(x, y_top, z),V_j);
			   }
		   }
	   }
 
   // set up and down particles directly
	   for (int i = 0; i < (x_num + 2* BP ); i++)
	   {
		   for (int j = 0; j < (y_num + 2* BP ); j++)
		   {
			   for (int k = 0; k < BP; k++)
			   {
				   Real x = (i + 0.5 ) * resolution_ref_large;
				   Real y = (j + 0.5 ) * resolution_ref_large;
				   Real z = (k + 0.5) * resolution_ref;
				   addPositionAndVolumetricMeasure(Vec3d(x, y, z),V_j);

				   Real z_top = (z_num + k + 0.5+BP) * resolution_ref;
				   addPositionAndVolumetricMeasure(Vec3d(x, y, z_top), V_j);
			   }
		   }
	   }
 
    }
};
 
 

 
class DiffusionInitialCondition : public LocalDynamics
{
  public:
    DiffusionInitialCondition(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()),
          pos_(particles_->getVariableDataByName<Vec3d>("Position")),
          species_(particles_->getVariableDataByName<Real>("Species")){};
    virtual ~DiffusionInitialCondition(){};

    Vec3d *pos_;
    Real  *species_;

  protected:
    void update(size_t index_i, Real dt = 0.0)
    { 
  species_[index_i] = pos_[index_i][0] *pos_[index_i][0] + pos_[index_i][1] * pos_[index_i][1] + pos_[index_i][2] * pos_[index_i][2];  
      
    };
};
 
 
class GetLaplacianTimeStepSize : public LocalDynamicsReduce<ReduceMin>                    
{
  protected:
    Real smoothing_length;
   
    AnisotropicDiffusionSolid &anisotropic_diffusion_solid_;

  public:
    GetLaplacianTimeStepSize(SPHBody &sph_body)
        : LocalDynamicsReduce<ReduceMin>(sph_body),
       anisotropic_diffusion_solid_(DynamicCast<AnisotropicDiffusionSolid>(this, sph_body.getBaseMaterial())) 
    {
        smoothing_length = sph_body.sph_adaptation_->ReferenceSmoothingLength();
    };

    Real reduce(size_t index_i, Real dt)
    {
        return 0.5 * smoothing_length * smoothing_length  * smoothing_length  
             / anisotropic_diffusion_solid_.DiffusivityCoefficient() / Dimensions;
    }

    virtual ~GetLaplacianTimeStepSize(){};


};
 } // namespace SPH
   
//------------------------------------------------------------------------------
// the main program
//------------------------------------------------------------------------------

int main(int ac, char *av[])
{
   
    //----------------------------------------------------------------------
    //	Build up the environment of a SPHSystem.
    //----------------------------------------------------------------------
    SPHSystem sph_system(system_domain_bounds, resolution_ref_large);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody boundary_body(sph_system, makeShared<DiffusionBoundary>("DiffusionBoundary"));
    boundary_body.sph_adaptation_->resetKernel<AnisotropicKernel<KernelWendlandC2>>(scaling_vector);
    boundary_body.defineMaterial<Solid>();
    boundary_body.generateParticles<BaseParticles, DiffusionBoundary>();

    SolidBody diffusion_body(sph_system, makeShared<DiffusionBlock>("DiffusionBlock"));
    diffusion_body.sph_adaptation_->resetKernel<AnisotropicKernel<KernelWendlandC2>>(scaling_vector);
    diffusion_body.defineMaterial<AnisotropicDiffusionSolid>(rho0, diffusion_coeff);
    diffusion_body.generateParticles<BaseParticles, DiffusionBlock>();

    //----------------------------------------------------------------------
    //	Define body relation map.
    //	The contact map gives the topological connections between the bodies.
    //	Basically the the range of bodies to build neighbor particle lists.
    //  Generally, we first define all the inner relations, then the contact relations.
    //  At last, we define the complex relaxations by combining previous defined
    //  inner and contact relations.
    //----------------------------------------------------------------------
    InnerRelation diffusion_body_inner_relation(diffusion_body); 
    ContactRelation diffusion_block_contact(diffusion_body, {&boundary_body});

   
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------
    InteractionWithUpdate<AnisotropicLinearGradientCorrectionMatrixComplex> correct_configuration(diffusion_body_inner_relation, diffusion_block_contact);
   	InteractionWithUpdate<AnisotropicKernelCorrectionMatrixComplex> correct_second_configuration(diffusion_body_inner_relation, diffusion_block_contact);
    ReduceDynamics<GetLaplacianTimeStepSize> get_time_step_size(diffusion_body);
    InteractionWithUpdate<AnisotropicDiffusionRelaxationComplex> diffusion_relaxation(diffusion_body_inner_relation, diffusion_block_contact);
    SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body_inner_relation);

    
    BodyStatesRecordingToVtp write_states_vtp(sph_system);

    write_states_vtp.addToWrite<Real>(diffusion_body, "Species");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_x");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_y");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_z");
	  write_states_vtp.addToWrite<Mat3d>(diffusion_body,"LinearGradientCorrectionMatrix");
 
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
    sph_system.initializeSystemConfigurations();
    correct_configuration.exec();
    correct_second_configuration.exec();
    setup_diffusion_initial_condition.exec();
   //----------------------------------------------------------------------
    //	Setup for time-stepping control
    //----------------------------------------------------------------------
    Real &physical_time = *sph_system.getSystemVariableDataByName<Real>("PhysicalTime");
    size_t number_of_iterations = 0;
    
    Real end_time = 1.0;
    Real output_interval = 0.01 * end_time;
    Real Observe_time = 0.01 * output_interval;
    Real dt = 0.0;

     //----------------------------------------------------------------------
    //	Statistics for CPU time
    //----------------------------------------------------------------------
    TickCount t1 = TickCount::now();
    TimeInterval interval;
    //----------------------------------------------------------------------
    //	First output before the main loop.
    //----------------------------------------------------------------------
    write_states_vtp.writeToFile();

 
	//----------------------------------------------------------------------
    //	Main loop starts here.
    //----------------------------------------------------------------------
    while (physical_time < end_time)
    {
        Real integration_time = 0.0;
        while (integration_time < output_interval)
        {
            Real relaxation_time = 0.0;
            while (relaxation_time < Observe_time)
            {
              dt = 0.1 * scaling_factor * get_time_step_size.exec();
               diffusion_relaxation.exec(dt);
               if (number_of_iterations <3) //for comparing the results with the theoretical results
                {
                  write_states_vtp.writeToFile(number_of_iterations);
                }  
           
                if (number_of_iterations % 1000 == 0)
                {
                    std::cout << "N=" << number_of_iterations << " Time: "
                              << physical_time << "	dt: "
                              << dt << "\n";
                }

                number_of_iterations++;
                relaxation_time += dt;
                integration_time += dt;
                physical_time += dt;
            } 
			
        }

        TickCount t2 = TickCount::now();
        write_states_vtp.writeToFile(number_of_iterations);
        TickCount t3 = TickCount::now();
        interval += t3 - t2;
    }
    TickCount t4 = TickCount::now();

    TimeInterval tt;
    tt = t4 - t1 - interval;
    std::cout << "Total wall time for computation: " << tt.seconds() << " seconds." << std::endl;

    return 0;
}
