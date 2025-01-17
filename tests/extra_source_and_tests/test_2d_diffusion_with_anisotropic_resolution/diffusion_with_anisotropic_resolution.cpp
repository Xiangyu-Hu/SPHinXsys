/**
 * @file 	diffusion_with_anisotropic_resolution.cpp
 * @brief 	This is a test to validate our anisotropic diffusion algorithm with anisotropic resolution.
 * @author Xiaojing Tang and Xiangyu Hu 
 */
#include "base_particles.h"
#include "observer_particles.h"
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
Real youngs_modulus = 1.0;
Real poisson_ratio = 1.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the case.
//----------------------------------------------------------------------
std::vector<Vec2d> diffusion_shape{Vec2d(0.0, 0.0), Vec2d(0.0, H), Vec2d(L, H), Vec2d(L, 0.0), Vec2d(0.0, 0.0)};

AnisotropicKernel<KernelWendlandC2>
    wendland(1.15 * resolution_ref_large, scaling_vector, Vecd(0.0, 0.0)); // no rotation introduced

Mat2d transform_tensor = wendland.getCoordinateTransformationTensorG(scaling_vector, Vecd(0.0, 0.0)); // tensor

class DiffusionBlock : public MultiPolygonShape
{
  public:
    explicit DiffusionBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(diffusion_shape, ShapeBooleanOps::add);
    }
};

class Boundary: public MultiPolygonShape
{
  public:
    explicit Boundary(const std::string &shape_name) : MultiPolygonShape(shape_name)
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
class AnisotropicParticleGenerator : public ParticleGenerator<BaseParticles>
{
  public:
  explicit  AnisotropicParticleGenerator(SPHBody &sph_body, BaseParticles &base_particles) 
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
                addPositionAndVolumetricMeasure(Vecd(x, y), (resolution_ref * resolution_ref_large));
            }
        }
    }
};

//particle generator for boundary particles
class AnisotropicParticleGeneratorBoundary : public  ParticleGenerator<BaseParticles>
{
  public:
    AnisotropicParticleGeneratorBoundary(SPHBody &sph_body,BaseParticles &base_particles) 
     :  ParticleGenerator<BaseParticles>(sph_body, base_particles){};

    virtual void prepareGeometricData() override
    {
        // set particles directly
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < y_num + 6; j++)
            {
                Real x = -BL + (i + 0.5) * resolution_ref_large;
                Real y = -BH + (j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vec2d(x, y), (resolution_ref * resolution_ref_large));
            }
        }

        // set particles directly
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < y_num + 6; j++)
            {
                Real x = (x_num + i + 0.5) * resolution_ref_large;
                Real y = -BH + (j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vec2d(x, y), (resolution_ref * resolution_ref_large));
            }
        }

        // set particles directly
        for (int i = 0; i < x_num; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Real x = (i + 0.5) * resolution_ref_large;
                Real y = -BH + (j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vec2d(x, y), (resolution_ref * resolution_ref_large));
            }
        }

        // set particles directly
        for (int i = 0; i < x_num; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                Real x = (i + 0.5) * resolution_ref_large;
                Real y = (y_num + j + 0.5) * resolution_ref;
                addPositionAndVolumetricMeasure(Vec2d(x, y), (resolution_ref * resolution_ref_large));
            }
        }
    }
};

//material definition
class LaplacianDiffusionSolid : public LinearElasticSolid
{
  public:
    LaplacianDiffusionSolid(Real rho0, Real coeff, Real youngs_modulus, Real poisson_ratio)
        : LinearElasticSolid(rho0, youngs_modulus, poisson_ratio), diffusion_coeff(coeff)
    {
        material_type_name_ = "LaplacianDiffusionSolid";
    };
    virtual ~LaplacianDiffusionSolid(){};

    Real diffusion_coeff;
    Real DiffusivityCoefficient() { return diffusion_coeff; };
};

class LaplacianDiffusionSolid;

class LaplacianDiffusionParticles : public BaseParticles
{
  public:
    LaplacianDiffusionSolid &laplacian_solid_;

    LaplacianDiffusionParticles(SPHBody &sph_body, LaplacianDiffusionSolid *laplacian_solid)
        : BaseParticles(sph_body, laplacian_solid), laplacian_solid_(*laplacian_solid),
           phi_(nullptr),  A1_(nullptr),  A2_(nullptr),  A3_(nullptr){};
        
    virtual ~LaplacianDiffusionParticles(){};


  void registerAnisotropicDiffusionProperties
        (StdLargeVec<Real> &phi, StdLargeVec<Vecd> &A1,  
              StdLargeVec<Vecd> &A2, StdLargeVec<Vecd> &A3)
{
    phi_ = registerStateVariableFrom<Real >("Phi", phi);
    A1_ = registerStateVariableFrom<Vecd>("FirstOrderCorrectionVectorA1", A1);
    A2_ = registerStateVariableFrom<Vecd>("FirstOrderCorrectionVectorA2", A2);
    A3_ = registerStateVariableFrom<Vecd>("FirstOrderCorrectionVectorA3", A3);
  
    addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA1");
    addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA2");
    addVariableToWrite<Vec2d>("FirstOrderCorrectionVectorA3");
}

        Real *phi_;
        Vec2d *A1_;
        Vec2d *A2_;
        Vec2d *A3_;

/*

    StdLargeVec<Real> *phi_;
    StdLargeVec<Vec2d> *A1_;
    StdLargeVec<Vec2d> *A2_;
    StdLargeVec<Vec2d> *A3_;*/
};

typedef DataDelegateSimple<LaplacianDiffusionParticles> LaplacianSolidDataSimple;
typedef DataDelegateInner<LaplacianDiffusionParticles> LaplacianSolidDataInner;
typedef DataDelegateContact<LaplacianDiffusionParticles, BaseParticles, DataDelegateEmptyBase> LaplacianSolidDataContactOnly;
typedef DataDelegateComplex<LaplacianDiffusionParticles, BaseParticles> LaplacianSolidDataComplex;
typedef DataDelegateComplex<BaseParticles, BaseParticles>GeneralDataDelegateComplex;
 

class NonisotropicKernelCorrectionMatrixComplex : 
    public LinearGradientCorrectionMatrix<DataDelegateContact>
{
  public:
    NonisotropicKernelCorrectionMatrixComplex(ComplexRelation &complex_relation, Real alpha = Real(0))
     : LinearGradientCorrectionMatrix<DataDelegateContact>(complex_relation), 
		B_(*particles_->registerSharedVariable<Mat2d>("KernelCorrectionMatrix")) {};

    virtual ~NonisotropicKernelCorrectionMatrixComplex(){};

  protected:
	  StdLargeVec<Mat2d> &B_;
 
	  void initialization(size_t index_i, Real dt = 0.0)
	  {
		  Mat2d local_configuration = Eps * Matd::Identity();
		  const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
         {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd gradW_ij = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];
            local_configuration -= r_ji * gradW_ij.transpose();
    
		  }
		  B_[index_i] = local_configuration;

	  };

    void interaction(size_t index_i, Real dt = 0.0)
    { 
		 Mat2d local_configuration = Eps * Mat2d::Identity();
		for (size_t k = 0; k < contact_configuration_.size(); ++k)
		{
			Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
			for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
			{
                 size_t index_j = contact_neighborhood.j_[n];
				Vec2d r_ji = contact_neighborhood.r_ij_[n] *contact_neighborhood.e_ij_[n] ;
				Vec2d gradW_ij = contact_neighborhood.dW_ij_[n] * Vol_[index_j]  * contact_neighborhood.e_ij_[n];
				local_configuration -= r_ji * gradW_ij.transpose();
			}
		}
		B_[index_i] += local_configuration; 
     }; 

	void update(size_t index_i, Real dt)
	{
		Mat2d inverse = B_[index_i].inverse();
		B_[index_i] = inverse;	
	}
};



class NonisotropicKernelCorrectionMatrixComplexAC : public LocalDynamics, public LaplacianSolidDataComplex
{
  public:
    NonisotropicKernelCorrectionMatrixComplexAC(ComplexRelation &complex_relation): 
                LocalDynamics(complex_relation.getInnerRelation().getSPHBody()), LaplacianSolidDataComplex(complex_relation),
                 B_(particles_->B_), A1_(particles_->A1_), A2_(particles_->A2_), A3_(particles_->A3_) {};

    virtual ~NonisotropicKernelCorrectionMatrixComplexAC(){};

  protected:
    StdLargeVec<Mat2d> &B_;
    StdLargeVec<Vec2d> &A1_,&A2_,&A3_;

    void initialization(size_t index_i, Real dt = 0.0)
    {
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ik
        {
            Vec2d gradW_ikV_k = inner_neighborhood.dW_ijV_j_[n] * inner_neighborhood.e_ij_[n];
            Vec2d r_ik = -inner_neighborhood.r_ij_vector_[n];

            A1_[index_i] += r_ik[0] * r_ik[0] * (B_[index_i].transpose() * gradW_ikV_k);
            A2_[index_i] += r_ik[1] * r_ik[1] * (B_[index_i].transpose() * gradW_ikV_k);
            A3_[index_i] += r_ik[0] * r_ik[1] * (B_[index_i].transpose() * gradW_ikV_k);
        }
    };

    void interaction(size_t index_i, Real dt = 0.0)
    {
         for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {
                Vec2d r_ik = -contact_neighborhood.r_ij_vector_[n];
                Vec2d gradW_ikV_k = contact_neighborhood.dW_ijV_j_[n] * contact_neighborhood.e_ij_[n];

                A1_[index_i] += r_ik[0] * r_ik[0] * (B_[index_i].transpose() * gradW_ikV_k);
                A2_[index_i] += r_ik[1] * r_ik[1] * (B_[index_i].transpose() * gradW_ikV_k);
                A3_[index_i] += r_ik[0] * r_ik[1] * (B_[index_i].transpose() * gradW_ikV_k);
            }
        }

    };

	void update(size_t index_i, Real dt = 0.0) {};
};

class LaplacianBodyRelaxation : public LocalDynamics, public LaplacianSolidDataComplex
{
  public:
    LaplacianBodyRelaxation(ComplexRelation &complex_relation): 
                LocalDynamics(complex_relation.getInnerRelation().getSPHBody()), LaplacianSolidDataComplex(complex_relation),
          pos_(particles_->pos_), B_(particles_->B_), phi_(particles_->phi_), A1_(particles_->A1_), A2_(particles_->A2_), A3_(particles_->A3_)
    {
        particles_->registerVariable(SC_, "FirstOrderCorrectionMatrixSC", [&](size_t i) -> Mat3d { return Eps * Mat3d::Identity(); });
        particles_->registerVariable(E_, "FirstOrderCorrectionVectorE", [&](size_t i) -> Vec2d { return Eps * Vec2d::Identity(); });
        particles_->registerVariable(G_, "FirstOrderCorrectionVectorG", [&](size_t i) -> Vec3d { return Eps * Vec3d::Identity(); });
        particles_->registerVariable(Laplacian_, "Laplacian", [&](size_t i) -> Vec3d { return Vec3d::Zero(); });

        particles_->registerVariable(Laplacian_x, "Laplacian_x", [&](size_t i) -> Real { return Real(0.0); });
        particles_->registerVariable(Laplacian_y, "Laplacian_y", [&](size_t i) -> Real { return Real(0.0); });
        particles_->registerVariable(Laplacian_xy, "Laplacian_xy", [&](size_t i) -> Real { return Real(0.0); });
		particles_->registerVariable(diffusion_dt_, "diffusion_dt", [&](size_t i) -> Real { return Real(0.0); });
		
        diffusion_coeff_ = particles_->laplacian_solid_.DiffusivityCoefficient();

	     /*StdLargeVec<Real>  contact_phi_;  
		contact_phi_.resize(this->contact_particles_.size());

		 
		for (size_t k = 0; k != this->contact_particles_.size(); ++k)
		{
			 
			contact_phi_[k].push_back(&all_contact_species_k[contact_species_index_k_m]);
		}
		 */ 
 
    };
    virtual ~LaplacianBodyRelaxation(){};

    StdLargeVec<Vec2d> &pos_;
    StdLargeVec<Mat2d> &B_;
    StdLargeVec<Real> &phi_;

    StdLargeVec<Vec2d> &A1_, &A2_, &A3_;

    StdLargeVec<Mat3d> SC_;
    StdLargeVec<Vec2d> E_;
    StdLargeVec<Vec3d> G_;
    StdLargeVec<Vec3d> Laplacian_;

    StdLargeVec<Real> Laplacian_x, Laplacian_y, Laplacian_xy, diffusion_dt_;

    Real diffusion_coeff_;

  protected:
    void initialization(size_t index_i, Real dt = 0.0)
    {   
        Neighborhood &inner_neighborhood = inner_configuration_[index_i];
        Vec2d E_rate = Vec2d::Zero();
        for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ik
        {
              
            size_t index_k = inner_neighborhood.j_[n];
            Vec2d gradW_ikV_k = inner_neighborhood.dW_ij_[n] * Vol_[index_k] * inner_neighborhood.e_ij_[n];

            E_rate += (phi_[index_k] - phi_[index_i]) * (B_[index_i].transpose() * gradW_ikV_k); // HOW TO DEFINE IT???
        }
        E_[index_i] = E_rate;

         
        Vec3d G_rate = Vec3d::Zero();
        Mat3d SC_rate = Mat3d::Zero();
        Real H_rate = 1.0;
         for (size_t n = 0; n != inner_neighborhood.current_size_; ++n) // this is ik
        {
            size_t index_j = inner_neighborhood.j_[n];
            Vecd r_ji = inner_neighborhood.r_ij_[n] * inner_neighborhood.e_ij_[n];

            Vec2d gradW_ijV_j = inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            Vec3d S_ = Vec3d(r_ij[0] * r_ij[0], r_ij[1] * r_ij[1], r_ij[0] * r_ij[1]);
            H_rate = r_ij.dot(B_[index_i].transpose() * gradW_ijV_j) / pow(r_ij.norm(), 4.0);
			 
            Real FF_ = 2.0 * (phi_[index_j] - phi_[index_i] - r_ij.dot(E_[index_i]));
            G_rate += S_ *H_rate * FF_;
            
             //TO DO
            Vec3d C_ = Vec3d::Zero();
            C_[0] = (r_ij[0] * r_ij[0]- r_ij.dot(A1_[index_i]));
            C_[1] = (r_ij[1] * r_ij[1]- r_ij.dot(A2_[index_i]));
            C_[2] = (r_ij[0] * r_ij[1]- r_ij.dot(A3_[index_i]));
            SC_rate += S_ *H_rate* C_.transpose();   

        }
        
        G_[index_i] = G_rate;
        SC_[index_i] = SC_rate;
  
    };

    void interaction(size_t index_i, Real dt = 0.0)
    {
        Mat3d SC_rate_contact = Mat3d::Zero();
        Vec3d G_rate_contact = Vec3d::Zero();
        Real H_rate_contact = 1.0;
         for (size_t k = 0; k < contact_configuration_.size(); ++k)
        {
            Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
            for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
            {   
                 size_t index_j = contact_neighborhood.j_[n];
                Vec2d r_ij = -contact_neighborhood.r_ij_[n] * contact_neighborhood.e_ij_[n];
                Vec2d gradW_ijV_j = contact_neighborhood.dW_ij_[n] * Vol_[index_j]* contact_neighborhood.e_ij_[n];

            
                Vec3d S_ = Vec3d(r_ij[0] * r_ij[0], r_ij[1] * r_ij[1], r_ij[0] * r_ij[1]);
                Real FF_ = 2.0 * (0.0 - r_ij.dot(E_[index_i])); ///here when it is periodic boundary condition, should notice the 0.0
                H_rate_contact = r_ij.dot(B_[index_i].transpose() * gradW_ijV_j) / pow(r_ij.norm(), 4.0);
		
                //TO DO
                Vec3d C_ = Vec3d::Zero();
                  C_[0] = (r_ij[0] * r_ij[0]- r_ij.dot(A1_[index_i]));
                C_[1] = (r_ij[1] * r_ij[1]- r_ij.dot(A2_[index_i]));
                C_[2] = (r_ij[0] * r_ij[1]- r_ij.dot(A3_[index_i]));
          

                SC_rate_contact += S_ * H_rate_contact * C_.transpose();
                G_rate_contact += S_ * H_rate_contact * FF_;

            }
			SC_[index_i] += SC_rate_contact;
            G_[index_i] += G_rate_contact;
        }

        Laplacian_[index_i] = diffusion_coeff_ * SC_[index_i].inverse() * G_[index_i];

        Laplacian_x[index_i] = Laplacian_[index_i][0];
        Laplacian_y[index_i] = Laplacian_[index_i][1];
        Laplacian_xy[index_i] = Laplacian_[index_i][2];
		diffusion_dt_[index_i] = Laplacian_[index_i][0] + Laplacian_[index_i][1];
    };

    void update(size_t index_i, Real dt = 0.0)
    {
        phi_[index_i] += dt * (Laplacian_[index_i][0] + Laplacian_[index_i][1]);
    };
};

class DiffusionInitialCondition : public LocalDynamics, public LaplacianSolidDataInner
{
  public:
    DiffusionInitialCondition(BaseInnerRelation &inner_relation)
        : LocalDynamics(inner_relation.getSPHBody()), LaplacianSolidDataInner(inner_relation),
          pos_(particles_->pos_), phi_(particles_->phi_){};
    virtual ~DiffusionInitialCondition(){};

    StdLargeVec<Vec2d> &pos_;
    StdLargeVec<Real> &phi_;

  protected:
    void update(size_t index_i, Real dt = 0.0)
    {
        if (pos_[index_i][0] >= 0.4* L && pos_[index_i][0] <= 0.6 * L)
        {
            phi_[index_i] = 1.0;     
        }
      
    };
};



class BoundaryCondition : public LocalDynamics, public LaplacianSolidDataInner
{
public:
	BoundaryCondition(BaseInnerRelation &inner_relation)
		: LocalDynamics(inner_relation.getSPHBody()), LaplacianSolidDataInner(inner_relation),
		pos_(particles_->pos_), phi_(particles_->phi_) {};
	virtual ~BoundaryCondition() {};

	StdLargeVec<Vec2d> &pos_;
	StdLargeVec<Real> &phi_;

protected:
	void update(size_t index_i, Real dt = 0.0)
	{
		phi_[index_i] = 1.0;
	};
};

class GetLaplacianTimeStepSize : public LocalDynamicsReduce<Real, ReduceMin>,
                                 public LaplacianSolidDataSimple
{
  protected:
    Real smoothing_length;

  public:
    GetLaplacianTimeStepSize(SPHBody &sph_body)
        : LocalDynamicsReduce<Real, ReduceMin>(sph_body, Real(MaxRealNumber)),
          LaplacianSolidDataSimple(sph_body)
    {
        smoothing_length = sph_body.sph_adaptation_->ReferenceSmoothingLength();
    };

    Real reduce(size_t index_i, Real dt)
    {
        return 0.5 * smoothing_length * smoothing_length / particles_->laplacian_solid_.DiffusivityCoefficient() / Dimensions;
    }

    virtual ~GetLaplacianTimeStepSize(){};
};

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
    SolidBody diffusion_body(sph_system, makeShared<DiffusionBlock>("DiffusionBlock"));
    diffusion_body.sph_adaptation_->resetKernel<AnisotropicKernel<KernelWendlandC2>>(scaling_vector);
    diffusion_body.defineMaterial<LaplacianDiffusionSolid>(rho0, diffusion_coeff, youngs_modulus, poisson_ratio);
    diffusion_body.generateParticles<AnisotropicParticleGenerator>();

    SolidBody boundary_body(sph_system, makeShared<Boundary>("Boundary"));
    boundary_body.sph_adaptation_->resetKernel<AnisotropicKernel<KernelWendlandC2>>(scaling_vector);
    boundary_body.defineMaterial<Solid>();
    boundary_body.generateParticles<AnisotropicParticleGeneratorBoundary>();

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
	InnerRelation boundary_body_inner_relation(boundary_body);

    ContactRelation temperature_observer_contact(temperature_observer, {&diffusion_body});
    //contact realtion???? todo or ComplexRelation
    ContactRelation diffusion_block_complex(diffusion_body, {&boundary_body});
    //----------------------------------------------------------------------
    //	Define the main numerical methods used in the simulation.
    //	Note that there may be data dependence on the constructors of these methods.
    //----------------------------------------------------------------------

	Dynamics1Level<NonisotropicKernelCorrectionMatrixComplex> correct_configuration(diffusion_block_complex);
	Dynamics1Level<NonisotropicKernelCorrectionMatrixComplexAC> correct_second_configuration(diffusion_block_complex);
    ReduceDynamics<GetLaplacianTimeStepSize> get_time_step_size(diffusion_body);
    Dynamics1Level<LaplacianBodyRelaxation> diffusion_relaxation(diffusion_block_complex);

    SimpleDynamics<DiffusionInitialCondition> setup_diffusion_initial_condition(diffusion_body_inner_relation);

/*  
    write_states_vtp.addToWrite<Real>(diffusion_body, "Phi");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_x");
    write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_y");
    write_states_vtp.addToWrite<Mat2d>(diffusion_body,"KernelCorrectionMatrix");
	write_states_vtp.addToWrite<Real>(diffusion_body,"Laplacian_xy");
	write_states_vtp.addToWrite<Real>(diffusion_body,"diffusion_dt");
    write_states_vtp.addToWrite<Vec2d>(diffusion_body,"FirstOrderCorrectionVectorE");
*/ 	

    PeriodicAlongAxis periodic_along_x(diffusion_body.getSPHBodyBounds(), xAxis);
    PeriodicAlongAxis periodic_along_y(diffusion_body.getSPHBodyBounds(), yAxis);

	PeriodicConditionUsingCellLinkedList periodic_condition_y(diffusion_body, periodic_along_y);
	PeriodicConditionUsingCellLinkedList periodic_condition_x(diffusion_body,periodic_along_x);

    //----------------------------------------------------------------------
    //	Define the methods for I/O operations and observations of the simulation.
    //----------------------------------------------------------------------
    BodyStatesRecordingToVtp write_states_vtp(sph_system);
   /* RegressionTestEnsembleAverage<ObservedQuantityRecording<Real>>
        write_solid_temperature("Phi", io_environment, temperature_observer_contact);*/
    //----------------------------------------------------------------------
    //	Prepare the simulation with cell linked list, configuration
    //	and case specified initial condition if necessary.
    //----------------------------------------------------------------------
    sph_system.initializeSystemCellLinkedLists();
	periodic_condition_y.update_cell_linked_list_.exec();
	periodic_condition_x.update_cell_linked_list_.exec();
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

    std::cout << "transform_tensor: "<< transform_tensor << std::endl;
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
