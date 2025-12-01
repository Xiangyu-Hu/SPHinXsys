/**
 * @file 	cohesive_soil_failure.h
 * @brief 	2D cohesive soil failure.
 * @details This case employs the SPH to simulate cohesive granular materials.
 *          The transport velocity formulation is adopted to address the tensile instability.
 * @author Shuaihao Zhang
 */
#include "sphinxsys.h" //SPHinXsys Library.
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real LL = 4.0;                        /**< Soil column length. */
Real LH = 2.0;                        /**< Soil column height. */
Real DL = 5.0;                        /**< Tank length. */
Real DH = 2.5;                        /**< Tank height. */
Real particle_spacing_ref = LL / 100; /**< Initial reference particle spacing. */
Real BW = particle_spacing_ref * 4;   /**< Extending width for boundary conditions. */
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(DL + BW, DH + BW));
//----------------------------------------------------------------------
//	Material properties of the soil.
//----------------------------------------------------------------------
Real rho0_s = 1850;                                                       // reference density of soil
Real gravity_g = 9.8;                                                     // gravity force of soil
Real Youngs_modulus = 1.8e6;                                              // reference Youngs modulus
Real poisson = 0.3;                                                       // Poisson ratio
Real c_s = sqrt(Youngs_modulus / (rho0_s * 3.0 * (1.0 - 2.0 * poisson))); // sound speed
Real cohesion = 5.0e3;
Real friction_angle = 25.0 * Pi / 180;
//----------------------------------------------------------------------
//	Geometric shapes used in this case.
//----------------------------------------------------------------------
Vec2d soil_block_halfsize = Vec2d(0.5 * LL, 0.5 * LH); // local center at origin:
Vec2d soil_block_translation = soil_block_halfsize;
Vec2d outer_wall_halfsize = Vec2d(0.5 * DL + BW, 0.5 * DH + BW);
Vec2d outer_wall_translation = Vec2d(-BW, -BW) + outer_wall_halfsize;
Vec2d inner_wall_halfsize = Vec2d(0.5 * DL, 0.5 * DH);
Vec2d inner_wall_translation = inner_wall_halfsize;
//----------------------------------------------------------------------
//	Complex for wall boundary
//----------------------------------------------------------------------
class WallBoundary : public ComplexShape
{
  public:
    explicit WallBoundary(const std::string &shape_name) : ComplexShape(shape_name)
    {
        add<GeometricShapeBox>(Transform(outer_wall_translation), outer_wall_halfsize);
        subtract<GeometricShapeBox>(Transform(inner_wall_translation), inner_wall_halfsize);
    }
};
std::vector<Vecd> soil_shape{
    Vecd(0, 0), Vecd(0, LH), Vecd(LL, LH), Vecd(LL, 0), Vecd(0, 0)};

class Soil : public MultiPolygonShape
{
  public:
    explicit Soil(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(soil_shape, ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	application dependent initial condition
//----------------------------------------------------------------------
class SoilInitialCondition : public continuum_dynamics::ContinuumInitialCondition
{
  public:
    explicit SoilInitialCondition(RealBody &granular_column)
        : continuum_dynamics::ContinuumInitialCondition(granular_column) {};

  protected:
    void update(size_t index_i, Real dt)
    {
        /** initial stress */
        Real y = pos_[index_i][1];
        Real gama = 1 - sin(friction_angle);
        Real stress_yy = -rho0_s * gravity_g * y;
        stress_tensor_3D_[index_i](1, 1) = stress_yy;
        stress_tensor_3D_[index_i](0, 0) = stress_yy * gama;
        stress_tensor_3D_[index_i](2, 2) = stress_yy * gama;
    };
};
//----------------------------------------------------------------------
//	Unified transport velocity correction
//----------------------------------------------------------------------
template <typename... T>
class TransportVelocityCorrection;

template <class AdaptationType, class LimiterType, typename... CommonControlTypes>
class TransportVelocityCorrection<Inner<AdaptationType, LimiterType>, CommonControlTypes...>
    : public fluid_dynamics::TransportVelocityCorrection<Base, DataDelegateInner, CommonControlTypes...>
{
    using SmoothingRatioType = typename AdaptationType::SmoothingRatioType;

  public:
    explicit TransportVelocityCorrection(BaseInnerRelation &inner_relation, Real coefficient = 0.2)
        : fluid_dynamics::TransportVelocityCorrection<Base, DataDelegateInner, CommonControlTypes...>(inner_relation),
          h_ref_(this->getSPHAdaptation().ReferenceSmoothingLength()),
          correction_scaling_(coefficient * h_ref_ * h_ref_),
          Vol_(this->particles_->template getVariableDataByName<Real>("VolumetricMeasure")),
          pos_div_(this->particles_->template registerStateVariableData<Real>("PositionDivergence")),
          pos_(this->particles_->template getVariableDataByName<Vecd>("Position")),
          h_ratio_(DynamicCast<AdaptationType>(this, this->getSPHAdaptation())), limiter_(h_ref_ * h_ref_),
          indicator_(this->particles_->template registerStateVariableData<int>("Indicator")),
          corner_indicator_(this->particles_->template registerStateVariableData<int>("CornerIndicator")),
          surface_normal_(this->particles_->template registerStateVariableData<Vecd>("SurfaceNormal"))
    {
        static_assert(std::is_base_of<Limiter, LimiterType>::value,
                      "Limiter is not the base of LimiterType!");
    }
    virtual ~TransportVelocityCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (this->within_scope_(index_i))
        {
            Vecd inconsistency = Vecd::Zero();
            const Neighborhood &inner_neighborhood = this->inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                // acceleration for transport velocity
                inconsistency -= (this->kernel_correction_(index_i) + this->kernel_correction_(index_j)) *
                                 inner_neighborhood.dW_ij_[n] * this->Vol_[index_j] * inner_neighborhood.e_ij_[n];
            }
            this->kernel_gradient_integral_[index_i] = inconsistency;
        }
    };
    void update(size_t index_i, Real dt = 0.0)
    {
        if (this->within_scope_(index_i))
        {
            Real inv_h_ratio = 1.0 / h_ratio_(index_i);
            Real squared_norm = this->kernel_gradient_integral_[index_i].squaredNorm();
            Vecd pos_transport = correction_scaling_ * limiter_(squared_norm) *
                                 this->kernel_gradient_integral_[index_i] * inv_h_ratio * inv_h_ratio;
            if (this->indicator_[index_i])
            {
                pos_transport = pos_transport - pos_transport.dot(this->surface_normal_[index_i]) * this->surface_normal_[index_i];
                if (this->pos_div_[index_i] < 0.6 * Dimensions)
                    pos_transport = Vecd::Zero();
            }
            pos_[index_i] += pos_transport;
        }
    };

  protected:
    const Real h_ref_, correction_scaling_;
    Real *Vol_, *pos_div_;
    Vecd *pos_;
    SmoothingRatioType h_ratio_;
    LimiterType limiter_;
    int *indicator_, *corner_indicator_;
    Vecd *surface_normal_;
};
template <class LimiterType, class ParticleScope>
using TransportVelocityCorrectionInner =
    TransportVelocityCorrection<Inner<SPHAdaptation, LimiterType>, NoKernelCorrection, ParticleScope>;

template <typename... CommonControlTypes>
class TransportVelocityCorrection<Contact<Boundary>, CommonControlTypes...>
    : public fluid_dynamics::TransportVelocityCorrection<Base, DataDelegateContact, CommonControlTypes...>
{
  public:
    explicit TransportVelocityCorrection(BaseContactRelation &contact_relation)
        : fluid_dynamics::TransportVelocityCorrection<Base, DataDelegateContact, CommonControlTypes...>(contact_relation)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            wall_Vol_.push_back(this->contact_particles_[k]->template getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };
    virtual ~TransportVelocityCorrection() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (this->within_scope_(index_i))
        {
            Vecd inconsistency = Vecd::Zero();
            for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
            {
                Real *wall_Vol_k = wall_Vol_[k];
                Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    // acceleration for transport velocity
                    inconsistency -= 2.0 * this->kernel_correction_(index_i) * contact_neighborhood.dW_ij_[n] *
                                     wall_Vol_k[index_j] * contact_neighborhood.e_ij_[n];
                }
            }
            this->kernel_gradient_integral_[index_i] += inconsistency;
        }
    };

  protected:
    StdVec<Real *> wall_Vol_;
};
template <class AdaptationType, class LimiterType, typename... CommonControlTypes>
using BaseTransportVelocityCorrectionComplex =
    ComplexInteraction<TransportVelocityCorrection<Inner<AdaptationType, LimiterType>, Contact<Boundary>>, CommonControlTypes...>;

template <class ParticleScope>
using TransportVelocityCorrectionComplex =
    BaseTransportVelocityCorrectionComplex<SPHAdaptation, NoLimiter, NoKernelCorrection, ParticleScope>;

//----------------------------------------------------------------------
//	Free surface normal direction
//----------------------------------------------------------------------
template <typename... InteractionTypes>
class FreeSurfaceNormal;

template <class DataDelegationType>
class FreeSurfaceNormal<DataDelegationType>
    : public LocalDynamics, public DataDelegationType
{
  public:
    template <class BaseRelationType>
    explicit FreeSurfaceNormal(BaseRelationType &base_relation)
        : LocalDynamics(base_relation.getSPHBody()), DataDelegationType(base_relation),
          surface_normal_(particles_->registerStateVariableData<Vecd>("SurfaceNormal")),
          color_gradient_(particles_->registerStateVariableData<Vecd>("ColorGradient")),
          B_(particles_->getVariableDataByName<Matd>("LinearGradientCorrectionMatrix")),
          indicator_(particles_->registerStateVariableData<int>("Indicator")),
          Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure"))
    {
        particles_->addEvolvingVariable<Vecd>("SurfaceNormal");
        particles_->addEvolvingVariable<Vecd>("ColorGradient");
    };
    virtual ~FreeSurfaceNormal() {};

  protected:
    Vecd *surface_normal_, *color_gradient_;
    Matd *B_;
    int *indicator_;
    Real *Vol_;
};

template <>
class FreeSurfaceNormal<Inner<>>
    : public FreeSurfaceNormal<DataDelegateInner>
{
  public:
    explicit FreeSurfaceNormal(BaseInnerRelation &inner_relation) : FreeSurfaceNormal<DataDelegateInner>(inner_relation) {};
    virtual ~FreeSurfaceNormal() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i])
        {
            Vecd color_gradient = ZeroData<Vecd>::value;
            const Neighborhood &inner_neighborhood = inner_configuration_[index_i];
            for (size_t n = 0; n != inner_neighborhood.current_size_; ++n)
            {
                size_t index_j = inner_neighborhood.j_[n];
                color_gradient -= inner_neighborhood.dW_ij_[n] * Vol_[index_j] * inner_neighborhood.e_ij_[n];
            }
            color_gradient_[index_i] = color_gradient;
        }
    };
    void update(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i])
        {
            surface_normal_[index_i] = B_[index_i] * color_gradient_[index_i] / (B_[index_i] * color_gradient_[index_i]).norm();
        }
        else
        {
            surface_normal_[index_i] = ZeroData<Vecd>::value;
        }
    };
};
using FreeSurfaceNormalInner = FreeSurfaceNormal<Inner<>>;

template <>
class FreeSurfaceNormal<Contact<>>
    : public FreeSurfaceNormal<DataDelegateContact>
{
  public:
    explicit FreeSurfaceNormal(BaseContactRelation &contact_relation)
        : FreeSurfaceNormal<DataDelegateContact>(contact_relation)
    {
        for (size_t k = 0; k != this->contact_particles_.size(); ++k)
        {
            contact_Vol_.push_back(contact_particles_[k]->getVariableDataByName<Real>("VolumetricMeasure"));
        }
    };
    virtual ~FreeSurfaceNormal() {};
    void interaction(size_t index_i, Real dt = 0.0)
    {
        if (indicator_[index_i])
        {
            for (size_t k = 0; k < contact_configuration_.size(); ++k)
            {
                Vecd color_gradient = ZeroData<Vecd>::value;
                Real *Vol_k = contact_Vol_[k];
                Neighborhood &contact_neighborhood = (*contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    color_gradient -= contact_neighborhood.dW_ij_[n] * Vol_k[index_j] * contact_neighborhood.e_ij_[n];
                }
                color_gradient_[index_i] += color_gradient;
            }
        }
    };

  protected:
    StdVec<Real *> contact_Vol_;
};
using FreeSurfaceNormalComplex =
    ComplexInteraction<FreeSurfaceNormal<Inner<>, Contact<>>>;