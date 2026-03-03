/**
 * @file 	windows_frame_diffusion_D4.h
 * @brief 	Head files used by windows_frame_diffusion_D4.cpp.
 * @author	Haotian Ji, Dong Wu, Chi Zhang and Xiangyu Hu
 */
#ifndef WINDOWS_FRAME_DIFFUSION_D4_H
#define WINDOWS_FRAME_DIFFUSION_D4_H

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 0.3; // Unit in m.
Real H = 0.093;
Real global_resolution = 0.001;
Real BW = global_resolution * 2.0;
BoundingBoxd system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
std::string diffusion_species_name = "Phi";
// Const for calculating air cavities conductivity
Real C1 = 0.025;        // Unit W/(m*K)
Real C3 = 1.57;         // Unit W/(m2*K)
Real C4 = 2.11;         // Unit W/(m*K)
Real wood_cond = 0.13;  // Unit W/(m*k), soft wood conductivity
Real epdm_cond = 0.25;  // epdm conductivity
Real pane_cond = 0.035; // insulation panel conductivity

Real getACConductivity(Real b, Real d, Real A) // calculate air cavities conductivity
{
    Real b_equal = sqrt(A * b / d);
    Real d_equal = sqrt(A * d / b);

    Real ha = 0.0;
    if (b_equal < 0.005)
    {
        ha = C1 / d_equal;
    }
    else
    {
        ha = SMAX(C1 / d_equal, C3);
    }

    Real hr = C4 * (1 - d_equal / b_equal + sqrt(1 + pow(d_equal / b_equal, 2)));
    Real Rs = 1 / (ha + hr);
    Real cond = d_equal / Rs;

    return cond;
}

/*---Geometric parameter of air cavities, unit m and m2, d is the cavity dimension
in the heat flow rate direction, b is the cavity dimension perpendicular to the heat
flow rate direction. The Area parameters for non-rectangular are given, for rectangular
formulas are given.---*/
Real d1 = 0.054;
Real b1 = 0.006;
Real A1 = d1 * b1;
Real d2 = 0.034;
Real b2 = 0.005;
Real A2 = d2 * b2;
Real od1 = 0.018;
Real ob1 = 0.005;
Real oA1 = od1 * ob1;

// unventilated air cavities conductivity
Real ac1_cond = getACConductivity(b1, d1, A1);
Real ac2_cond = getACConductivity(b2, d2, A2);

// slightly ventilated air conductivity
Real ac1_open_cond = 2 * getACConductivity(ob1, od1, oA1);
//----------------------------------------------------------------------
//	Initial and boundary condition parameters.
//----------------------------------------------------------------------
// Temperature initialization,unit Celsius
Real initial_temperature = 10.0;
Real T_infinity_e = 0.0;
Real T_infinity_i = 20.0;

// Surface resistance, unit (K*m2)/W
Real rs_i = 0.13;
Real rs_i_increased = 0.20;
Real rs_e = 0.04;

// Convection coefficient, unit W/(K*m2)
Real convection_e = 1 / rs_e;
Real convection_i = 1 / rs_i;
Real convection_i_decreased = 1 / rs_i_increased;

namespace SPH
{
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
MultiPolygon createOverallStructureBody()
{
    std::vector<Vecd> overallStructureDomainShape;
    overallStructureDomainShape.push_back(Vecd(0.0, 0.005));
    overallStructureDomainShape.push_back(Vecd(0.0, 0.071));
    overallStructureDomainShape.push_back(Vecd(0.026, 0.071));
    overallStructureDomainShape.push_back(Vecd(0.026, 0.088));
    overallStructureDomainShape.push_back(Vecd(0.11, 0.088));
    overallStructureDomainShape.push_back(Vecd(0.11, 0.051));
    overallStructureDomainShape.push_back(Vecd(0.3, 0.051));
    overallStructureDomainShape.push_back(Vecd(0.3, 0.023));
    overallStructureDomainShape.push_back(Vecd(0.11, 0.023));
    overallStructureDomainShape.push_back(Vecd(0.11, 0.005));
    overallStructureDomainShape.push_back(Vecd(0.0, 0.005));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(overallStructureDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createInternalAirBody()
{
    std::vector<Vecd> internalAirDomainShape;
    internalAirDomainShape.push_back(Vecd(0.0, 0.071));
    internalAirDomainShape.push_back(Vecd(0.0, 0.076));
    internalAirDomainShape.push_back(Vecd(0.009, 0.076));
    internalAirDomainShape.push_back(Vecd(0.009, 0.093));
    internalAirDomainShape.push_back(Vecd(0.140, 0.093));
    internalAirDomainShape.push_back(Vecd(0.140, 0.056));
    internalAirDomainShape.push_back(Vecd(0.3, 0.056));
    internalAirDomainShape.push_back(Vecd(0.3, 0.051));
    internalAirDomainShape.push_back(Vecd(0.11, 0.051));
    internalAirDomainShape.push_back(Vecd(0.11, 0.088));
    internalAirDomainShape.push_back(Vecd(0.026, 0.088));
    internalAirDomainShape.push_back(Vecd(0.026, 0.071));
    internalAirDomainShape.push_back(Vecd(0.0, 0.071));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(internalAirDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createDecreasedInternalConvectionBody()
{
    std::vector<Vecd> decreasedInternalConvectionDomainShape1;
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.009, 0.071));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.026, 0.088));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.026, 0.071));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.009, 0.071));

    std::vector<Vecd> decreasedInternalConvectionDomainShape2;
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.110, 0.051));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.110, 0.088));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.140, 0.051));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.110, 0.051));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(decreasedInternalConvectionDomainShape1, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(decreasedInternalConvectionDomainShape2, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createExternalAirBody()
{
    std::vector<Vecd> externalAirDomainShape;
    externalAirDomainShape.push_back(Vecd(0.0, 0.0));
    externalAirDomainShape.push_back(Vecd(0.0, 0.005));
    externalAirDomainShape.push_back(Vecd(0.110, 0.005));
    externalAirDomainShape.push_back(Vecd(0.110, 0.023));
    externalAirDomainShape.push_back(Vecd(0.300, 0.023));
    externalAirDomainShape.push_back(Vecd(0.300, 0.018));
    externalAirDomainShape.push_back(Vecd(0.115, 0.018));
    externalAirDomainShape.push_back(Vecd(0.115, 0.0));
    externalAirDomainShape.push_back(Vecd(0.0, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(externalAirDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createWoodBody()
{
    std::vector<Vecd> woodDomainShape1;
    woodDomainShape1.push_back(Vecd(0.000, 0.005));
    woodDomainShape1.push_back(Vecd(0.000, 0.071));
    woodDomainShape1.push_back(Vecd(0.042, 0.071));
    woodDomainShape1.push_back(Vecd(0.042, 0.020));
    woodDomainShape1.push_back(Vecd(0.063, 0.020));
    woodDomainShape1.push_back(Vecd(0.063, 0.005));
    woodDomainShape1.push_back(Vecd(0.000, 0.005));

    std::vector<Vecd> woodDomainShape2;
    woodDomainShape2.push_back(Vecd(0.068, 0.005));
    woodDomainShape2.push_back(Vecd(0.068, 0.023));
    woodDomainShape2.push_back(Vecd(0.048, 0.023));
    woodDomainShape2.push_back(Vecd(0.048, 0.074));
    woodDomainShape2.push_back(Vecd(0.026, 0.074));
    woodDomainShape2.push_back(Vecd(0.026, 0.088));
    woodDomainShape2.push_back(Vecd(0.110, 0.088));
    woodDomainShape2.push_back(Vecd(0.110, 0.054));
    woodDomainShape2.push_back(Vecd(0.090, 0.054));
    woodDomainShape2.push_back(Vecd(0.090, 0.020));
    woodDomainShape2.push_back(Vecd(0.110, 0.020));
    woodDomainShape2.push_back(Vecd(0.110, 0.005));
    woodDomainShape2.push_back(Vecd(0.068, 0.005));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(woodDomainShape1, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(woodDomainShape2, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createEPDMBody()
{
    std::vector<Vecd> epdmDomainShape1;
    epdmDomainShape1.push_back(Vecd(0.026, 0.071));
    epdmDomainShape1.push_back(Vecd(0.026, 0.074));
    epdmDomainShape1.push_back(Vecd(0.042, 0.074));
    epdmDomainShape1.push_back(Vecd(0.042, 0.071));
    epdmDomainShape1.push_back(Vecd(0.026, 0.071));

    std::vector<Vecd> epdmDomainShape2;
    epdmDomainShape2.push_back(Vecd(0.048, 0.020));
    epdmDomainShape2.push_back(Vecd(0.048, 0.023));
    epdmDomainShape2.push_back(Vecd(0.063, 0.023));
    epdmDomainShape2.push_back(Vecd(0.063, 0.020));
    epdmDomainShape2.push_back(Vecd(0.048, 0.020));

    std::vector<Vecd> epdmDomainShape3;
    epdmDomainShape3.push_back(Vecd(0.095, 0.051));
    epdmDomainShape3.push_back(Vecd(0.095, 0.054));
    epdmDomainShape3.push_back(Vecd(0.11, 0.054));
    epdmDomainShape3.push_back(Vecd(0.11, 0.051));
    epdmDomainShape3.push_back(Vecd(0.095, 0.051));

    std::vector<Vecd> epdmDomainShape4;
    epdmDomainShape4.push_back(Vecd(0.095, 0.020));
    epdmDomainShape4.push_back(Vecd(0.095, 0.023));
    epdmDomainShape4.push_back(Vecd(0.11, 0.023));
    epdmDomainShape4.push_back(Vecd(0.11, 0.020));
    epdmDomainShape4.push_back(Vecd(0.095, 0.020));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(epdmDomainShape1, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(epdmDomainShape2, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(epdmDomainShape3, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(epdmDomainShape4, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createPanelBody()
{
    std::vector<Vecd> panelDomainShape;
    panelDomainShape.push_back(Vecd(0.095, 0.023));
    panelDomainShape.push_back(Vecd(0.095, 0.051));
    panelDomainShape.push_back(Vecd(0.3, 0.051));
    panelDomainShape.push_back(Vecd(0.3, 0.023));
    panelDomainShape.push_back(Vecd(0.095, 0.023));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(panelDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody1()
{
    std::vector<Vecd> acDomainShape1;
    acDomainShape1.push_back(Vecd(0.042, 0.020));
    acDomainShape1.push_back(Vecd(0.042, 0.074));
    acDomainShape1.push_back(Vecd(0.048, 0.074));
    acDomainShape1.push_back(Vecd(0.048, 0.020));
    acDomainShape1.push_back(Vecd(0.042, 0.020));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape1, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody2()
{
    std::vector<Vecd> acDomainShape2;
    acDomainShape2.push_back(Vecd(0.090, 0.020));
    acDomainShape2.push_back(Vecd(0.090, 0.054));
    acDomainShape2.push_back(Vecd(0.095, 0.054));
    acDomainShape2.push_back(Vecd(0.095, 0.020));
    acDomainShape2.push_back(Vecd(0.090, 0.020));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape2, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACOpenBody1()
{
    std::vector<Vecd> acOpenDomainShape1;
    acOpenDomainShape1.push_back(Vecd(0.063, 0.005));
    acOpenDomainShape1.push_back(Vecd(0.063, 0.023));
    acOpenDomainShape1.push_back(Vecd(0.068, 0.023));
    acOpenDomainShape1.push_back(Vecd(0.068, 0.005));
    acOpenDomainShape1.push_back(Vecd(0.063, 0.005));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acOpenDomainShape1, ShapeBooleanOps::add);

    return multi_polygon;
}

//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
template <class DynamicsIdentifier>
class LocalQuantityDefinition
    : public BaseLocalDynamics<DynamicsIdentifier>
{
  public:
    LocalQuantityDefinition(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier) {};
    virtual ~LocalQuantityDefinition() {};
};

class LocalDiffusivityDefinition : public LocalQuantityDefinition<BodyPartByParticle>
{
  public:
    explicit LocalDiffusivityDefinition(BodyPartByParticle &body_part, Real local_diff)
        : LocalQuantityDefinition<BodyPartByParticle>(body_part),
          thermal_conductivity(particles_->getVariableDataByName<Real>("ThermalConductivity")),
          local_diff(local_diff) {};

    void update(size_t index_i, Real dt)
    {
        thermal_conductivity[index_i] = local_diff;
    };

  protected:
    Real *thermal_conductivity;
    Real local_diff;
};

class DiffusionInitialCondition : public LocalDynamics
{
  public:
    explicit DiffusionInitialCondition(SPHBody &sph_body)
        : LocalDynamics(sph_body),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = initial_temperature;
    };

  protected:
    Real *phi_;
};

class RobinBoundaryDefinition : public LocalDynamics
{
  public:
    explicit RobinBoundaryDefinition(SolidBody &diffusion_body)
        : LocalDynamics(diffusion_body),
          pos_(particles_->getVariableDataByName<Vecd>("Position")),
          phi_(particles_->registerStateVariableData<Real>(diffusion_species_name)),
          phi_convection_(particles_->template getVariableDataByName<Real>(diffusion_species_name + "Convection")),
          phi_infinity_(*(this->particles_->template getSingularVariableByName<Real>(diffusion_species_name + "Infinity")->Data())) {};

    void update(size_t index_i, Real dt)
    {
        phi_[index_i] = -0.0;

        if (pos_[index_i][1] >= 0.051)
        {
            phi_convection_[index_i] = convection_i;
            phi_infinity_ = T_infinity_i;
        }
        if (pos_[index_i][1] <= 0.023)
        {
            phi_convection_[index_i] = convection_e;
            phi_infinity_ = T_infinity_e;
        }
    };

  protected:
    Vecd *pos_;
    Real *phi_, *phi_convection_;
    Real &phi_infinity_;
};

class LocalConvectionDefinition : public LocalQuantityDefinition<BodyPartByParticle>
{
  public:
    explicit LocalConvectionDefinition(BodyPartByParticle &body_part, Real local_convection)
        : LocalQuantityDefinition<BodyPartByParticle>(body_part),
          phi_convection_(particles_->template getVariableDataByName<Real>(diffusion_species_name + "Convection")),
          local_convection_(local_convection) {};

    void update(size_t index_i, Real dt)
    {
        phi_convection_[index_i] = local_convection_;
    };

  protected:
    Real *phi_convection_;
    Real local_convection_;
};

using DiffusionBodyRelaxation =
    DiffusionBodyRelaxationComplex<BaseDiffusion, KernelGradientInner, KernelGradientContact, Robin>;

StdVec<Vecd> createObservationPoints()
{
    StdVec<Vecd> observation_points;
    /** A line of measuring points at the given position. */
    size_t number_of_observation_points = 5;
    Real range_of_measure = H - 0.02;
    Real start_of_measure = 0.01;

    for (size_t i = 0; i < number_of_observation_points; ++i)
    {
        Vec2d point_coordinate(
            0.028, range_of_measure * Real(i) / Real(number_of_observation_points - 1) + start_of_measure);
        observation_points.push_back(point_coordinate);
    }
    return observation_points;
};
} // namespace SPH
#endif // WINDOWS_FRAME_DIFFUSION_D4_H
