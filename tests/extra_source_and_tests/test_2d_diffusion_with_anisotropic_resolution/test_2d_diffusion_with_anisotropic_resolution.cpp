/**
 * @file 	diffusion_with_anisotropic_resolution.cpp
 * @brief 	This is a test to validate our anisotropic diffusion algorithm with anisotropic resolution.
 * @author Xiaojing Tang and Xiangyu Hu 
 */
#include "base_particles.h"
#include "contact_body_relation.h"
#include "inner_body_relation.h"
#include "kernel_correction.h"
#include "observer_particles.h"
#include "anisotropic_diffusion_relaxation_2d.hpp"


#include "sphinxsys.h"  
using namespace SPH;   // Namespace cite here
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 0.1;
 
int y_num = 10;  // Particle number in y direction
Real ratio_ = 4.0;   //The anisotropic resolution ratio 
Real resolution_ref = H / y_num;
Real resolution_ref_large = ratio_ * resolution_ref;
int x_num = L / resolution_ref_large;
Vec2d scaling_vector = Vec2d(1.0, 1.0 / ratio_); // scaling_vector for defining the anisotropic kernel
Real scaling_factor = 1.0 / ratio_;              // scaling factor to calculate the time step


Real V_j = resolution_ref_large * resolution_ref;
BoundingBox system_domain_bounds(Vec2d(-L, -H), Vec2d(2.0 * L, 2.0 * H));
Real BL = 3.0 * resolution_ref_large;
Real BH = 3.0 * resolution_ref;
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coeff = 1.0;
Real rho0 = 1.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the case.
//----------------------------------------------------------------------
std::vector<Vec2d> diffusion_shape{Vec2d(0.0, 0.0), Vec2d(0.0, H), Vec2d(L, H), Vec2d(L, 0.0), Vec2d(0.0, 0.0)};
 std::vector<Vec2d> boundary{Vec2d(-BL, -BH),Vec2d(-BL, H + BH), Vec2d(L + BL, H + BH),
  Vec2d(L + BL, -BH), Vec2d(-BL, -BH)};

namespace SPH
{
class DiffusionBlock : public MultiPolygonShape
{
public:
    explicit DiffusionBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(diffusion_shape, ShapeBooleanOps::add);
    }
};

class DiffusionBoundary: public MultiPolygonShape
{
public:
    explicit DiffusionBoundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
    
        multi_polygon_.addAPolygon(boundary, ShapeBooleanOps::add);
        multi_polygon_.addAPolygon(diffusion_shape, ShapeBooleanOps::sub);
    }
};

//----------------------------------------------------------------------
//	An observer particle generator.
//----------------------------------------------------------------------

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the entrance of the channel. */
    size_t number_of_observation_points = 11;
    Real range_of_measure = 0.9 * L;
    Real start_of_measure = 0.05 * L;

    /** the measuring locations */
     for (size_t i = 0; i < number_of_observation_points; ++i)
        {
            Vec2d point_coordinate(range_of_measure * (Real)i / (Real)(number_of_observation_points - 1) + start_of_measure, 0.5 * H);
            observation_points.push_back(point_coordinate);
        }

    return observation_points;
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
                Real x = (i + 0.5) * resolution_ref_large;
                Real y = (j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vecd(x, y), V_j);
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
        // set particles directly
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < y_num + 6; j++)
            {
                Real x = -BL + (i + 0.5) * resolution_ref_large;
                Real y = -BH + (j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vec2d(x, y), V_j);
            }
        }

        // set particles directly
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < y_num + 6; j++)
            {
                Real x = (x_num + i + 0.5) * resolution_ref_large;
                Real y = -BH + (j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vec2d(x, y), V_j);
            }
        }

        // set particles directly
        for (int i = 0; i < x_num; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Real x = (i + 0.5) * resolution_ref_large;
                Real y = -BH + (j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vec2d(x, y), V_j);
            }
        }

        // set particles directly
        for (int i = 0; i < x_num; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Real x = (i + 0.5) * resolution_ref_large;
                Real y = (y_num + j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vec2d(x, y), V_j);
            }
        }
    }
};

class DiffusionInitialCondition : public LocalDynamics
{
  public:
    DiffusionInitialCondition(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          species_(particles_->getVariableDataByName<Real>("Species")){};
    virtual ~DiffusionInitialCondition(){};

    Vec2d *pos_;
    Real  *species_;

  protected:
    void update(size_t index_i, Real dt = 0.0)
    {
        /*if (pos_[index_i][0] >= 0.4* L && pos_[index_i][0] <= 0.6 * L)
        {
            species_[index_i] = 1.0;     
        }*/
       species_[index_i] = pos_[index_i][0] *pos_[index_i][0] + pos_[index_i][1] * pos_[index_i][1];  
      
    };
};

/*
This class is used to check if the gradient is correct or not with specific function given.*/
class GradientCheck : public LocalDynamics, public DataDelegateInner
{
  public:
    GradientCheck(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()), DataDelegateInner(inner_relation),
        B_(this->particles_->template getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
        species_(particles_->getVariableDataByName<Real>("Species")),
        Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure"))
    ,show_neighbor_(particles_->registerStateVariable<Real>("ShowingNeighbor", Real(0.0)))
 {
       Gradient_x = particles_->registerStateVariable<Real>("Gradient_x", [&](size_t i) -> Real { return Real(0.0); });
       Gradient_y = particles_->registerStateVariable<Real>("Gradient_y", [&](size_t i) -> Real { return Real(0.0); });
    };

    virtual ~GradientCheck(){};
  
        
  protected:
    Mat2d *B_;
    Real  *species_;
    Real  *Vol_;
    Real *Gradient_x;
    Real *Gradient_y;
   Real *show_neighbor_;
    void interaction(size_t index_i, Real dt = 0.0)
    {
        Vec2d rate_ = Vec2d::Zero();
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ij
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vec2d gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            

             if (index_i == 59)
            {
                show_neighbor_[index_j] = 1.0;
            };

            rate_ += (species_[index_j] - species_[index_i]) * (B_[index_i].transpose() * gradW_ijV_j);
        }

        Gradient_x[index_i] = rate_[0];
        Gradient_y[index_i] = rate_[1];
    };

    void update(size_t index_i, Real dt = 0.0){};
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
        return 0.5 * smoothing_length * smoothing_length  
             / anisotropic_diffusion_solid_.DiffusivityCoefficient() / Dimensions;
    }

    virtual ~GetLaplacianTimeStepSize(){};


};
} // namespace SPH

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
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
    //	Particle and body creation of fluid observers.
    //----------------------------------------------------------------------
    ObserverBody temperature_observer(sph_system, "TemperatureObserver");
    temperature_observer.generateParticles<ObserverParticles>(createObservationPoints());

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
    InteractionWithUpdate<GradientCheck> gradient_check(diffusion_body_inner_relation);
    
    BodyStatesRecordingToVtp write_states_vtp(sph_system);

    write_states_vtp.addToWrite<Real>(diffusion_body, "Species");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_x");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_y");
    //write_states_vtp.addToWrite<Real>(diffusion_body,"Gradient_x");
   // write_states_vtp.addToWrite<Real>(diffusion_body,"Gradient_y");
    write_states_vtp.addToWrite<Matd>(diffusion_body,"LinearGradientCorrectionMatrix");
	//write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_xy");
	//write_states_vtp.addToWrite<Real>(diffusion_body,"diffusion_dt");

    write_states_vtp.addToWrite<Real>(diffusion_body, "ShowingNeighbor");
    write_states_vtp.addToWrite<Vecd>(diffusion_body, "FirstOrderCorrectionVector1");
    write_states_vtp.addToWrite<Vecd>(diffusion_body, "FirstOrderCorrectionVector2");
    write_states_vtp.addToWrite<Vecd>(diffusion_body, "FirstOrderCorrectionVector3");


    
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
    
    Real end_time = 0.4;
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
               gradient_check.exec(dt);  // for checking the first order accuracy
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
