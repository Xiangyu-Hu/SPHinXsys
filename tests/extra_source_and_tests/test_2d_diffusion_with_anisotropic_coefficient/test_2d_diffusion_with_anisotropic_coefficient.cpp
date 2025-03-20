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
Real L = 20.0;
Real H = 20.0;
  
int y_num = 100;
Real resolution_ref = H / y_num;
 
BoundingBox system_domain_bounds(Vec2d(-L, -H), Vec2d(2.0 * L, 2.0 * H));
Real BL = 3.0 * resolution_ref;
Real BH = 3.0 * resolution_ref;
//If this is anisotropic diffusion tensor
Mat2d transform_tensor{ 
    {0.1, 0.03},  
    {0.03, 0.03},
}; 
Mat2d inverse_decomposed_transform_tensor = inverseCholeskyDecomposition(transform_tensor);
Mat2d  decomposed_transform_tensor = inverse_decomposed_transform_tensor.inverse();

 
 //----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real rho0 = 1.0;


//----------------------------------------------------------------------
//	Geometric shapes used in the case.
//----------------------------------------------------------------------
std::vector<Vec2d> diffusion_shape{Vec2d(0.0, 0.0), Vec2d(0.0, H), Vec2d(L, H), Vec2d(L, 0.0), Vec2d(0.0, 0.0)};



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
        std::vector<Vec2d> boundary;
        boundary.push_back(Vec2d(-BL, -BH));
        boundary.push_back(Vec2d(-BL, H + BH));
        boundary.push_back(Vec2d(L + BL, H + BH));
        boundary.push_back(Vec2d(L + BL, -BH));
        boundary.push_back(Vec2d(-BL, -BH));
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
    size_t number_of_observation_points = 41;
    Real range_of_measure = 1.0 * H;
    Real start_of_measure = 0.0 * H;

    /** the measuring locations */
    for (size_t i = 0; i < number_of_observation_points; ++i)
    {
        Vec2d point_coordinate(0.5*L, range_of_measure * (Real)i / (Real)(number_of_observation_points - 1) + start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
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
        if (pos_[index_i][0] >= (0.5 * L- 0.5 ) && pos_[index_i][0] <=  (0.5 * L + 0.5 ))
        {
            if (pos_[index_i][1] >=  (0.5 * H - 0.5) && pos_[index_i][1] <= (0.5* H + 0.5))
            {
                species_[index_i] = 1.0;
            }
          
        }  
   
    };
};
 

class GetLaplacianTimeStepSize : public LocalDynamicsReduce<ReduceMin>                    
{
  protected:
    Real smoothing_length;
    Real diff_coeff; //should be the minimum vone in this tensor


  public:
    GetLaplacianTimeStepSize(SPHBody &sph_body)
        : LocalDynamicsReduce<ReduceMin>(sph_body)
    {
        smoothing_length =this->sph_body_.getSPHAdaptation().ReferenceSmoothingLength();
        diff_coeff = decomposed_transform_tensor.trace() /Dimensions ;

    };

    Real reduce(size_t index_i, Real dt)
    {
        return 0.5 * smoothing_length * smoothing_length / diff_coeff / Dimensions;
        
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
    SPHSystem sph_system(system_domain_bounds, resolution_ref);
    sph_system.handleCommandlineOptions(ac, av);
    IOEnvironment io_environment(sph_system);
    //----------------------------------------------------------------------
    //	Creating body, materials and particles.
    //----------------------------------------------------------------------
    SolidBody boundary_body(sph_system, makeShared<DiffusionBoundary>("DiffusionBoundary"));
    boundary_body.defineMaterial<Solid>();
    boundary_body.generateParticles<BaseParticles, Lattice>();

    SolidBody diffusion_body(sph_system, makeShared<DiffusionBlock>("DiffusionBlock"));
    diffusion_body.defineMaterial<AnisotropicDiffusionSolid>(rho0, decomposed_transform_tensor);
    diffusion_body.generateParticles<BaseParticles, Lattice>();


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
    
    BodyStatesRecordingToVtp write_states_vtp(sph_system);

    write_states_vtp.addToWrite<Real>(diffusion_body, "Species");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_x");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_y");
    write_states_vtp.addToWrite<Matd>(diffusion_body,"LinearGradientCorrectionMatrix");
    
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
    
    Real end_time = 100;
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
              dt = 0.1  * get_time_step_size.exec();
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
