/**
 * @file 	heart_volume_change.cpp
 * @brief 	This is the case studying the electromechanics on a biventricular heart model in 3D, including the volume change of the ventricles.
 * @author 	John Benjamin, Chi Zhang and Xiangyu Hu
 * 			Unit :
 *			time t = ms = 12.9 [-]
 * 			length l = mm
 * 			mass m = g
 *			density rho = g * (mm)^(-3)
 *			Pressure pa = g * (mm)^(-1) * (ms)^(-2)
 *			diffusion d = (mm)^(2) * (ms)^(-2)
 */
#include "heart_volume_change.h"
using namespace SPH; // Namespace cite here.
/** Geometry parameter. */
/** Set the file path to the stl file. */
std::string full_path_to_myocardium = "./input/myocardium_simple.stl";
std::string full_path_to_lv = "./input/myocardium_simple_lv.stl";
std::string full_path_to_rv = "./input/myocardium_simple_rv.stl";
Real length_scale = 1.0;
Real time_scale = 1.0 / 12.9;
Real stress_scale = 1.0e-6;
Vecd translation(-53.5 * length_scale, -70.0 * length_scale, -32.5 * length_scale);
/** Parameters and physical properties. */
Vec3d domain_lower_bound(-55.0 * length_scale, -75.0 * length_scale, -35.0 * length_scale);
Vec3d domain_upper_bound(35.0 * length_scale, 5.0 * length_scale, 35.0 * length_scale);
Real dp_0 = (domain_upper_bound[0] - domain_lower_bound[0]) / 45.0; /**< Initial particle spacing. */
/** Domain bounds of the system. */
BoundingBoxd system_domain_bounds(domain_lower_bound, domain_upper_bound);

/** Material properties. */
Real rho0_s = 1.06e-3;
/** Active stress factor */
Real k_a = 100 * stress_scale;
std::array<Real, 4> a0 = {Real(496.0) * stress_scale, Real(15196.0) * stress_scale, Real(3283.0) * stress_scale, Real(662.0) * stress_scale};
std::array<Real, 4> b0 = {Real(7.209), Real(20.417), Real(11.176), Real(9.466)};
/** reference stress to achieve weakly compressible condition */
Real poisson = 0.4995;
Real bulk_modulus = 2.0 * a0[0] * (1.0 + poisson) / (3.0 * (1.0 - 2.0 * poisson));
/** Electrophysiology parameters. */
std::string diffusion_species_name = "Phi";
Real diffusion_coeff = 0.8;
Real bias_coeff = 0.0;
/** Electrophysiology parameters. */
Real c_m = 1.0;
Real k = 8.0;
Real a = 0.01;
Real b = 0.15;
Real mu_1 = 0.2;
Real mu_2 = 0.3;
Real epsilon = 0.002;
/** Fibers and sheet. */
Vec3d fiber_direction(1.0, 0.0, 0.0);
Vec3d sheet_direction(0.0, 1.0, 0.0);

namespace SPH
{
/**
 * Define heart geometry
 */
class Heart : public ComplexShape
{
  public:
    explicit Heart(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<TriangleMeshShapeSTL>(full_path_to_myocardium, translation, length_scale);
    }
};
/** Set diffusion relaxation method. */
using FiberDirectionDiffusionRelaxation =
    DiffusionRelaxationRK2<DiffusionRelaxation<Inner<KernelGradientInner>, IsotropicDiffusion>>;
/** Imposing diffusion boundary condition */
class DiffusionBCs : public BaseLocalDynamics<BodyPartByParticle>
{
  public:
    explicit DiffusionBCs(BodyPartByParticle &body_part, const std::string &species_name)
        : BaseLocalDynamics<BodyPartByParticle>(body_part),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(species_name)) {};
    virtual ~DiffusionBCs() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        Vecd displacement = sph_body_->getInitialShape().findNormalDirection(pos_[index_i]);
        Vecd face_norm = displacement / (displacement.norm() + 1.0e-15);

        Vecd center_norm = pos_[index_i] / (pos_[index_i].norm() + 1.0e-15);

        Real angle = face_norm.dot(center_norm);
        if (angle >= 0.0)
        {
            phi_[index_i] = 1.0;
        }
        else
        {
            if (pos_[index_i][1] < -getSPHAdaptation().ReferenceSpacing())
                phi_[index_i] = 0.0;
        }
    };

  protected:
    Vecd *pos_;
    Real *phi_;
};
/** Compute Fiber and Sheet direction after diffusion */
class ComputeFiberAndSheetDirections : public LocalDynamics
{
  protected:
    LocallyOrthotropicMuscle &muscle_material_;
    Vecd *pos_;
    Real *phi_;
    Real beta_epi_, beta_endo_;
    Vecd center_line_vector_; // parallel to the ventricular centerline and pointing  apex-to-base

  public:
    explicit ComputeFiberAndSheetDirections(SPHBody &sph_body, const std::string &species_name)
        : LocalDynamics(sph_body),
          muscle_material_(DynamicCast<LocallyOrthotropicMuscle>(this, sph_body_->getBaseMaterial())),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(species_name))
    {
        center_line_vector_ = Vecd(0.0, 1.0, 0.0);
        beta_epi_ = -(70.0 / 180.0) * M_PI;
        beta_endo_ = (80.0 / 180.0) * M_PI;
    };
    virtual ~ComputeFiberAndSheetDirections() {};

    void update(size_t index_i, Real dt = 0.0)
    {
        /**
         * Ref: original doi.org/10.1016/j.euromechsol.2013.10.009
         * 		Present  doi.org/10.1016/j.cma.2016.05.031
         */
        /** Probe the face norm from Level set field. */
        Vecd displacement = sph_body_->getInitialShape().findNormalDirection(pos_[index_i]);
        Vecd face_norm = displacement / (displacement.norm() + 1.0e-15);
        Vecd center_norm = pos_[index_i] / (pos_[index_i].norm() + 1.0e-15);
        if (face_norm.dot(center_norm) <= 0.0)
        {
            face_norm = -face_norm;
        }
        /** Compute the centerline's projection on the plane orthogonal to face norm. */
        Vecd circumferential_direction = getCrossProduct(center_line_vector_, face_norm);
        Vecd cd_norm = circumferential_direction / (circumferential_direction.norm() + 1.0e-15);
        /** The rotation angle is given by beta = (beta_epi - beta_endo) phi + beta_endo */
        Real beta = (beta_epi_ - beta_endo_) * phi_[index_i] + beta_endo_;
        /** Compute the rotation matrix through Rodrigues rotation formulation. */
        Vecd f_0 = cos(beta) * cd_norm + sin(beta) * getCrossProduct(face_norm, cd_norm) +
                   face_norm.dot(cd_norm) * (1.0 - cos(beta)) * face_norm;

        if (pos_[index_i][1] < -getSPHAdaptation().ReferenceSpacing())
        {
            muscle_material_.local_f0_[index_i] = f_0 / (f_0.norm() + 1.0e-15);
            muscle_material_.local_s0_[index_i] = face_norm;
        }
        else
        {
            muscle_material_.local_f0_[index_i] = Vecd::Zero();
            muscle_material_.local_s0_[index_i] = Vecd::Zero();
        }
    };
};

//	define shape parameters which will be used for the constrained body part.
class MuscleBaseShapeParameters : public TriangleMeshShapeBrick::ShapeParameters
{
  public:
    MuscleBaseShapeParameters() : TriangleMeshShapeBrick::ShapeParameters()
    {
        Real l = domain_upper_bound[0] - domain_lower_bound[0];
        Real w = domain_upper_bound[2] - domain_lower_bound[2];
        halfsize_ = Vec3d(0.5 * l, 1.0 * dp_0, 0.5 * w);
        resolution_ = 20;
        translation_ = Vec3d(-10.0 * length_scale, -1.0 * dp_0, 0.0);
    }
};

//	application dependent initial condition
class ApplyStimulusCurrentSI : public LocalDynamics
{
  public:
    explicit ApplyStimulusCurrentSI(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          voltage_(particles_->registerStateVariableData<Real>("Voltage")) {};

    void update(size_t index_i, Real dt)
    {
        if (-30.0 * length_scale <= pos_[index_i][0] && pos_[index_i][0] <= -15.0 * length_scale)
        {
            if (-2.0 * length_scale <= pos_[index_i][1] && pos_[index_i][1] <= 0.0)
            {
                if (-3.0 * length_scale <= pos_[index_i][2] && pos_[index_i][2] <= 3.0 * length_scale)
                {
                    voltage_[index_i] = 0.92;
                }
            }
        }
    };

  protected:
    Vecd *pos_;
    Real *voltage_;
};
/**
 * application dependent initial condition
 */
class ApplyStimulusCurrentSII : public LocalDynamics
{
  public:
    explicit ApplyStimulusCurrentSII(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          voltage_(particles_->registerStateVariableData<Real>("Voltage")) {};

    void update(size_t index_i, Real dt)
    {
        if (0.0 <= pos_[index_i][0] && pos_[index_i][0] <= 6.0 * length_scale)
        {
            if (-6.0 * length_scale <= pos_[index_i][1])
            {
                if (12.0 * length_scale <= pos_[index_i][2])
                {
                    voltage_[index_i] = 0.95;
                }
            }
        }
    };

  protected:
    Vecd *pos_;
    Real *voltage_;
};

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    observation_points.push_back(Vecd(-45.0 * length_scale, -30.0 * length_scale, 0.0));
    observation_points.push_back(Vecd(0.0, -30.0 * length_scale, 26.0 * length_scale));
    observation_points.push_back(Vecd(-30.0 * length_scale, -50.0 * length_scale, 0.0));
    observation_points.push_back(Vecd(0.0, -50.0 * length_scale, 20.0 * length_scale));
    observation_points.push_back(Vecd(0.0, -70.0 * length_scale, 0.0));
    return observation_points;
};
} // namespace SPH
