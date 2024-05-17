/**
 * @file 	windows_frame_diffusion_D7.h
 * @brief 	Head files used by windows_frame_diffusion_D7.cpp.
 * @author	Haotian Ji, Dong Wu, Chi Zhang and Xiangyu Hu
 */
#ifndef WINDOWS_FRAME_DIFFUSION_D7_H
#define WINDOWS_FRAME_DIFFUSION_D7_H

#include "sphinxsys.h"
using namespace SPH;
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 0.238; // Unit in m.
Real H = 0.109;
Real resolution_ref = 0.001;
Real BW = resolution_ref * 2.0;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
// Const for calculating air cavities conductivity
Real C1 = 0.025;        // Unit W/(m*K)
Real C3 = 1.57;         // Unit W/(m2*K)
Real C4 = 2.11;         // Unit W/(m*K)
Real poly_cond = 0.25;  // Unit W/(m*k), polyamide conductivity
Real epdm_cond = 0.25;  // epdm conductivity
Real pvc_cond = 0.17;  // pvc conductivity
Real pane_cond = 0.035; // insulation panel conductivity

Real getACConductivity(Real b, Real d, Real A)//calculate air cavities conductivity
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
Real d1 = 0.031;Real b1 = 0.025;Real A1 = 0.00058;
Real d2 = 0.009;Real b2 = 0.010;Real A2 = 0.000048;
Real d3 = 0.019;Real b3 = 0.012;Real A3 = d3 * b3;
Real d4 = 0.019;Real b4 = 0.025;Real A4 = 0.000367;
Real d5 = 0.03;Real b5 = 0.005;Real A5 = d5 * b5;
Real d6 = 0.035;Real b6 = 0.015;Real A6 = 0.000417;
Real d7 = 0.037;Real b7 = 0.036;Real A7 = 0.0006735;
Real od1 = 0.018;Real ob1 = 0.005;Real oA1 = od1 * ob1;

// unventilated air cavities conductivity
Real ac1_cond = getACConductivity(b1, d1, A1); 
Real ac2_cond = getACConductivity(b2, d2, A2);
Real ac3_cond = getACConductivity(b3, d3, A3);
Real ac4_cond = getACConductivity(b4, d4, A4);
Real ac5_cond = getACConductivity(b5, d5, A5);
Real ac6_cond = getACConductivity(b6, d6, A6);
Real ac7_cond = getACConductivity(b7, d7, A7);

// slightly ventilated air conductivity
Real acopen1_cond = 2 * getACConductivity(ob1, od1, oA1);
//----------------------------------------------------------------------
//	Initial and boundary condition parameters.
//----------------------------------------------------------------------
// Temperature initialization,unit Celsius
Real initial_temperature = 10.0;
Real T_infinity_e = 0.0;
Real T_infinity_i = 20.0;

//Surface resistance, unit (K*m2)/W
Real rs_i = 0.13;
Real rs_i_increased = 0.20;
Real rs_e = 0.04;

// Convection coefficient, unit W/(K*m2)
Real convection_e = 1 / rs_e;
Real convection_i = 1 / rs_i;
Real convection_i_decreased = 1 / rs_i_increased;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
MultiPolygon createOverallStructureBody()
{
    std::vector<Vecd> overallStructureDomainShape;
    overallStructureDomainShape.push_back(Vecd(0.0, 0.005));
    overallStructureDomainShape.push_back(Vecd(0.0, 0.104));
    overallStructureDomainShape.push_back(Vecd(0.031, 0.104));
    overallStructureDomainShape.push_back(Vecd(0.031, 0.092));
    overallStructureDomainShape.push_back(Vecd(0.048, 0.092));
    overallStructureDomainShape.push_back(Vecd(0.048, 0.064));
    overallStructureDomainShape.push_back(Vecd(0.238, 0.064));
    overallStructureDomainShape.push_back(Vecd(0.238, 0.04));
    overallStructureDomainShape.push_back(Vecd(0.048, 0.04));
    overallStructureDomainShape.push_back(Vecd(0.048, 0.034));
    overallStructureDomainShape.push_back(Vecd(0.031, 0.005));
    overallStructureDomainShape.push_back(Vecd(0.0, 0.005));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(overallStructureDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createInternalAirBody()
{
    std::vector<Vecd> internalAirDomainShape;
    internalAirDomainShape.push_back(Vecd(0.0, 0.104));
    internalAirDomainShape.push_back(Vecd(0.0, 0.109));
    internalAirDomainShape.push_back(Vecd(0.043, 0.109));
    internalAirDomainShape.push_back(Vecd(0.043, 0.097));
    internalAirDomainShape.push_back(Vecd(0.076, 0.097));
    internalAirDomainShape.push_back(Vecd(0.076, 0.069));
    internalAirDomainShape.push_back(Vecd(0.238, 0.069));
    internalAirDomainShape.push_back(Vecd(0.238, 0.064));
    internalAirDomainShape.push_back(Vecd(0.048, 0.064));
    internalAirDomainShape.push_back(Vecd(0.048, 0.092));
    internalAirDomainShape.push_back(Vecd(0.031, 0.092));
    internalAirDomainShape.push_back(Vecd(0.031, 0.104));
    internalAirDomainShape.push_back(Vecd(0.0, 0.104));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(internalAirDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createDecreasedInternalConvectionBody()
{
    std::vector<Vecd> decreasedInternalConvectionDomainShape1;
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.031, 0.092));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.031, 0.104));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.043, 0.092));
    decreasedInternalConvectionDomainShape1.push_back(Vecd(0.031, 0.092));

    std::vector<Vecd> decreasedInternalConvectionDomainShape2;
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.048, 0.064));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.048, 0.092));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.076, 0.064));
    decreasedInternalConvectionDomainShape2.push_back(Vecd(0.048, 0.064));

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
    externalAirDomainShape.push_back(Vecd(0.031, 0.005));
    externalAirDomainShape.push_back(Vecd(0.048, 0.034));
    externalAirDomainShape.push_back(Vecd(0.048, 0.04));
    externalAirDomainShape.push_back(Vecd(0.238, 0.04));
    externalAirDomainShape.push_back(Vecd(0.238, 0.035));
    externalAirDomainShape.push_back(Vecd(0.053, 0.035));
    externalAirDomainShape.push_back(Vecd(0.031, 0.0));
    externalAirDomainShape.push_back(Vecd(0.0, 0.0));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(externalAirDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createPolyamideBody()
{
    std::vector<Vecd> polyamideDomainShape1;
    polyamideDomainShape1.push_back(Vecd(0.031, 0.067));
    polyamideDomainShape1.push_back(Vecd(0.031, 0.073));
    polyamideDomainShape1.push_back(Vecd(0.021, 0.073));
    polyamideDomainShape1.push_back(Vecd(0.021, 0.079));
    polyamideDomainShape1.push_back(Vecd(0.031, 0.079));
    polyamideDomainShape1.push_back(Vecd(0.031, 0.092));
    polyamideDomainShape1.push_back(Vecd(0.048, 0.092));
    polyamideDomainShape1.push_back(Vecd(0.048, 0.067));
    polyamideDomainShape1.push_back(Vecd(0.031, 0.067));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(polyamideDomainShape1, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createEPDMBody()
{
    std::vector<Vecd> epdmDomainShape1;
    epdmDomainShape1.push_back(Vecd(0.036, 0.037));
    epdmDomainShape1.push_back(Vecd(0.036, 0.04));
    epdmDomainShape1.push_back(Vecd(0.048, 0.04));
    epdmDomainShape1.push_back(Vecd(0.048, 0.037));
    epdmDomainShape1.push_back(Vecd(0.036, 0.037));

    std::vector<Vecd> epdmDomainShape2;
    epdmDomainShape2.push_back(Vecd(0.036, 0.064));
    epdmDomainShape2.push_back(Vecd(0.036, 0.067));
    epdmDomainShape2.push_back(Vecd(0.048, 0.067));
    epdmDomainShape2.push_back(Vecd(0.048, 0.064));
    epdmDomainShape2.push_back(Vecd(0.036, 0.064));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(epdmDomainShape1, ShapeBooleanOps::add);
    multi_polygon.addAPolygon(epdmDomainShape2, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createPanelBody()
{
    std::vector<Vecd> panelDomainShape;
    panelDomainShape.push_back(Vecd(0.036, 0.04));
    panelDomainShape.push_back(Vecd(0.036, 0.064));
    panelDomainShape.push_back(Vecd(0.238, 0.064));
    panelDomainShape.push_back(Vecd(0.238, 0.04));
    panelDomainShape.push_back(Vecd(0.036, 0.04));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(panelDomainShape, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody1()
{
    std::vector<Vecd> acDomainShape1;
    acDomainShape1.push_back(Vecd(0.003, 0.070));
    acDomainShape1.push_back(Vecd(0.003, 0.101));
    acDomainShape1.push_back(Vecd(0.028, 0.101));
    acDomainShape1.push_back(Vecd(0.028, 0.085));
    acDomainShape1.push_back(Vecd(0.015, 0.085));
    acDomainShape1.push_back(Vecd(0.015, 0.070));
    acDomainShape1.push_back(Vecd(0.003, 0.070));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape1, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody2()
{
    std::vector<Vecd> acDomainShape2;
    acDomainShape2.push_back(Vecd(0.018, 0.073));
    acDomainShape2.push_back(Vecd(0.018, 0.082));
    acDomainShape2.push_back(Vecd(0.028, 0.082));
    acDomainShape2.push_back(Vecd(0.028, 0.079));
    acDomainShape2.push_back(Vecd(0.021, 0.079));
    acDomainShape2.push_back(Vecd(0.021, 0.073));
    acDomainShape2.push_back(Vecd(0.018, 0.073));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape2, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody3()
{
    std::vector<Vecd> acDomainShape3;
    acDomainShape3.push_back(Vecd(0.033, 0.070));
    acDomainShape3.push_back(Vecd(0.033, 0.089));
    acDomainShape3.push_back(Vecd(0.045, 0.089));
    acDomainShape3.push_back(Vecd(0.045, 0.070));
    acDomainShape3.push_back(Vecd(0.033, 0.070));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape3, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody4()
{
    std::vector<Vecd> acDomainShape4;
    acDomainShape4.push_back(Vecd(0.003, 0.054));
    acDomainShape4.push_back(Vecd(0.003, 0.067));
    acDomainShape4.push_back(Vecd(0.028, 0.067));
    acDomainShape4.push_back(Vecd(0.028, 0.048));
    acDomainShape4.push_back(Vecd(0.021, 0.048));
    acDomainShape4.push_back(Vecd(0.021, 0.054));
    acDomainShape4.push_back(Vecd(0.003, 0.054));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape4, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody5()
{
    std::vector<Vecd> acDomainShape5;
    acDomainShape5.push_back(Vecd(0.031, 0.037));
    acDomainShape5.push_back(Vecd(0.031, 0.067));
    acDomainShape5.push_back(Vecd(0.036, 0.067));
    acDomainShape5.push_back(Vecd(0.036, 0.037));
    acDomainShape5.push_back(Vecd(0.031, 0.037));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape5, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody6()
{
    std::vector<Vecd> acDomainShape6;
    acDomainShape6.push_back(Vecd(0.003, 0.016));
    acDomainShape6.push_back(Vecd(0.003, 0.051));
    acDomainShape6.push_back(Vecd(0.018, 0.051));
    acDomainShape6.push_back(Vecd(0.018, 0.04));
    acDomainShape6.push_back(Vecd(0.009, 0.016));
    acDomainShape6.push_back(Vecd(0.003, 0.016));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape6, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACBody7()
{
    std::vector<Vecd> acDomainShape7;
    acDomainShape7.push_back(Vecd(0.009, 0.008));
    acDomainShape7.push_back(Vecd(0.009, 0.013));
    acDomainShape7.push_back(Vecd(0.012, 0.013));
    acDomainShape7.push_back(Vecd(0.021, 0.04));
    acDomainShape7.push_back(Vecd(0.021, 0.045));
    acDomainShape7.push_back(Vecd(0.028, 0.045));
    acDomainShape7.push_back(Vecd(0.028, 0.034));
    acDomainShape7.push_back(Vecd(0.045, 0.034));
    acDomainShape7.push_back(Vecd(0.028, 0.008));
    acDomainShape7.push_back(Vecd(0.009, 0.008));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acDomainShape7, ShapeBooleanOps::add);

    return multi_polygon;
}
MultiPolygon createACOpenBody1()
{
    std::vector<Vecd> acOpenDomainShape1;
    acOpenDomainShape1.push_back(Vecd(0.003, 0.005));
    acOpenDomainShape1.push_back(Vecd(0.003, 0.013));
    acOpenDomainShape1.push_back(Vecd(0.006, 0.013));
    acOpenDomainShape1.push_back(Vecd(0.006, 0.005));
    acOpenDomainShape1.push_back(Vecd(0.003, 0.005));

    MultiPolygon multi_polygon;
    multi_polygon.addAPolygon(acOpenDomainShape1, ShapeBooleanOps::add);

    return multi_polygon;
}
//----------------------------------------------------------------------
//	Setup diffusion material properties.
//----------------------------------------------------------------------
class DiffusionMaterial : public DiffusionReaction<Solid>
{
  public:
    DiffusionMaterial() : DiffusionReaction<Solid>({"Phi"}, SharedPtr<NoReaction>())
    {
        initializeAnDiffusion<LocalIsotropicDiffusion>("Phi", "Phi", pvc_cond);
    }
};

class WallMaterial : public DiffusionReaction<Solid>
{
  public:
    WallMaterial() : DiffusionReaction<Solid>({"Phi"}, SharedPtr<NoReaction>())
    {
        initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", poly_cond);
        //polyamide conductivity is the largest, we use it to define the time step.
    }
};
using DiffusionParticles = DiffusionReactionParticles<SolidParticles, DiffusionMaterial>;
using WallParticles = DiffusionReactionParticles<SolidParticles, WallMaterial>;
//----------------------------------------------------------------------
//	Application dependent initial condition.
//----------------------------------------------------------------------
template <class DynamicsIdentifier, class ParticlesType>
class LocalQuantityDefinition
    : public BaseLocalDynamics<DynamicsIdentifier>,
      public DiffusionReactionSimpleData<ParticlesType>
{
  public:
    LocalQuantityDefinition(DynamicsIdentifier &identifier)
        : BaseLocalDynamics<DynamicsIdentifier>(identifier),
          DiffusionReactionSimpleData<ParticlesType>(identifier.getSPHBody()){};
    virtual ~LocalQuantityDefinition(){};
};
class ThermalConductivityInitialization
    : public LocalQuantityDefinition<BodyPartByParticle, DiffusionParticles>
{
  protected:
    StdLargeVec<Real> &thermal_conductivity;
    Real local_diff;

  public:
    explicit ThermalConductivityInitialization(BodyPartByParticle &body_part, Real local_diff)
        : LocalQuantityDefinition<BodyPartByParticle, DiffusionParticles>(body_part),
          thermal_conductivity(*(particles_->getVariableByName<Real>("ThermalConductivity"))),
          local_diff(local_diff){};

    void update(size_t index_i, Real dt)
    {
        thermal_conductivity[index_i] = local_diff;
    }
};
class LocalConvectionInitialization
    : public LocalQuantityDefinition<BodyPartByParticle, WallParticles>
{
  protected:
    StdLargeVec<Real> &convection_;
    Real local_convection;

  public:
    explicit LocalConvectionInitialization(BodyPartByParticle &body_part, Real local_convection)
        : LocalQuantityDefinition<BodyPartByParticle, WallParticles>(body_part),
          convection_(*(this->particles_->template getVariableByName<Real>("Convection"))),
          local_convection(local_convection){};

    void update(size_t index_i, Real dt)
    {
        convection_[index_i] = local_convection;
    }
};
class LocalHeatTransferConvection
    : public LocalQuantityDefinition<BodyPartByParticle, WallParticles>
{
  protected:
    StdLargeVec<Real> &ht_convection_;
    Real local_ht_convection;

  public:
    explicit LocalHeatTransferConvection(BodyPartByParticle &body_part, Real local_ht_convection)
        : LocalQuantityDefinition<BodyPartByParticle, WallParticles>(body_part),
          ht_convection_(*(this->particles_->template getVariableByName<Real>("HT_Convection"))),
          local_ht_convection(local_ht_convection){};

    void update(size_t index_i, Real dt)
    {
        ht_convection_[index_i] = local_ht_convection;
    }
};
class DiffusionInitialCondition
    : public DiffusionReactionInitialCondition<DiffusionParticles>
{
  protected:
    size_t phi_;

  public:
    explicit DiffusionInitialCondition(SPHBody &sph_body)
        : DiffusionReactionInitialCondition<DiffusionParticles>(sph_body)
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    };

    void update(size_t index_i, Real dt)
    {
        all_species_[phi_][index_i] = initial_temperature;
    };
};
class RobinWallBoundaryInitialCondition
    : public DiffusionReactionInitialCondition<WallParticles>
{
  protected:
    size_t phi_;
    StdLargeVec<Real> &convection_;
    Real &T_infinity_;

  public:
    explicit RobinWallBoundaryInitialCondition(SolidBody &diffusion_body)
        : DiffusionReactionInitialCondition<WallParticles>(diffusion_body),
          convection_(*(this->particles_->template getVariableByName<Real>("Convection"))),
          T_infinity_(*(this->particles_->template getSingleVariableByName<Real>("T_infinity")))
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    }

    void update(size_t index_i, Real dt)
    {
        all_species_[phi_][index_i] = -0.0;

        if (pos_[index_i][1] >= 0.064)
        {
            convection_[index_i] = convection_i;
            T_infinity_ = T_infinity_i;
        }
        if (pos_[index_i][1] <= 0.04)
        {
            convection_[index_i] = convection_e;
            T_infinity_ = T_infinity_e;
        }
    }
};
class ExternalHeatTransferInitialCondition
    : public DiffusionReactionInitialCondition<WallParticles>
{
  protected:
    size_t phi_;
    StdLargeVec<Real> &ht_convection_;
    Real &ht_T_infinity_;

  public:
    explicit ExternalHeatTransferInitialCondition(SolidBody &diffusion_body)
        : DiffusionReactionInitialCondition<WallParticles>(diffusion_body),
          ht_convection_(*(this->particles_->template getVariableByName<Real>("HT_Convection"))),
          ht_T_infinity_(*(this->particles_->template getSingleVariableByName<Real>("HT_T_infinity")))
    {
        phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    }

    void update(size_t index_i, Real dt)
    {
       ht_convection_[index_i] = convection_e;
       ht_T_infinity_ = T_infinity_e;
    }
};
class InternalHeatTransferInitialCondition
    : public DiffusionReactionInitialCondition<WallParticles>
{
  protected:
    size_t phi_;
    StdLargeVec<Real> &ht_convection_;
    Real &ht_T_infinity_;

  public:
    explicit InternalHeatTransferInitialCondition(SolidBody &diffusion_body)
        : DiffusionReactionInitialCondition<WallParticles>(diffusion_body),
          ht_convection_(*(this->particles_->template getVariableByName<Real>("HT_Convection"))),
          ht_T_infinity_(*(this->particles_->template getSingleVariableByName<Real>("HT_T_infinity")))
    {
       phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
    }

    void update(size_t index_i, Real dt)
    {
       ht_convection_[index_i] = convection_i;
       ht_T_infinity_ = T_infinity_i;
    }
};

template <typename... ControlTypes>
class RobinFlux; /*Calculate heat transfer flux of Robin boundary condition*/
template <class ParticlesType, class ContactParticlesType, class ContactKernelGradientType>
class DiffusionRelaxation<RobinFlux<ParticlesType, ContactParticlesType, ContactKernelGradientType>>
    : public DiffusionRelaxation<Contact<Base>, ParticlesType, ContactParticlesType, ContactKernelGradientType>
{
    StdLargeVec<Vecd> &n_;
    StdVec<StdLargeVec<Vecd> *> ht_n_;
    StdVec<Real *> ht_T_infinity_;
    StdVec<StdLargeVec<Real> *> contact_Vol_;
    StdVec<StdLargeVec<Real> *> ht_flux_;
    StdVec<StdLargeVec<Real> *> ht_convection_; // ht for heat transfer

  public:
    explicit DiffusionRelaxation(BaseContactRelation &contact_relation)
        : DiffusionRelaxation<Contact<Base>, ParticlesType, ContactParticlesType, ContactKernelGradientType>(contact_relation),
          n_(this->particles_->n_)
    {
       ht_flux_.resize(this->all_diffusions_.size());
       ht_T_infinity_.resize(this->all_diffusions_.size());
       ht_convection_.resize(this->contact_particles_.size());

       for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
       {
            ht_flux_[m] = this->particles_->template registerSharedVariable<Real>("HT_Flux");
            for (size_t k = 0; k != this->contact_particles_.size(); ++k)
            {
                ht_n_.push_back(&(this->contact_particles_[k]->n_));
                contact_Vol_.push_back(&(this->contact_particles_[k]->Vol_));
                ht_convection_[k] = this->contact_particles_[k]->template registerSharedVariable<Real>("HT_Convection");
                ht_T_infinity_[m] = this->contact_particles_[k]->template registerSingleVariable<Real>("HT_T_infinity");
            }
       }
    };
    virtual ~DiffusionRelaxation(){};

    void update(size_t index_i, Real dt)
    {
       for (size_t m = 0; m < this->all_diffusions_.size(); ++m)
       {
            for (size_t k = 0; k < this->contact_configuration_.size(); ++k)
            {
                Real ht_flux_ = 0.0;

                StdLargeVec<Vecd> &n_k = *(ht_n_[k]);
                StdLargeVec<Real> &Vol_k = *(contact_Vol_[k]);
                StdLargeVec<Real> &ht_convection_k = *(ht_convection_[k]);
                Real &ht_T_infinity_k = *(ht_T_infinity_[k]);

                Neighborhood &contact_neighborhood = (*this->contact_configuration_[k])[index_i];
                for (size_t n = 0; n != contact_neighborhood.current_size_; ++n)
                {
                    size_t index_j = contact_neighborhood.j_[n];
                    Real dW_ijV_j_ = contact_neighborhood.dW_ij_[n] * Vol_k[index_j];
                    Vecd &e_ij = contact_neighborhood.e_ij_[n];

                    const Vecd &grad_ijV_j = this->contact_kernel_gradients_[k](index_i, index_j, dW_ijV_j_, e_ij);
                    Vecd n_ij = n_[index_i] - n_k[index_j];
                    Real area_ij_Robin = grad_ijV_j.dot(n_ij);

                    Real phi_ij = ht_T_infinity_k - (*this->diffusion_species_[m])[index_i];
                    ht_flux_ += ht_convection_k[index_j] * phi_ij * area_ij_Robin * Vol_k[index_j];
                }
                (*this->ht_flux_[m])[index_i] = ht_flux_;
            }
       }
    };
};

using DiffusionBodyRelaxation = DiffusionBodyRelaxationComplex<
    DiffusionParticles, WallParticles, KernelGradientInner, KernelGradientContact, Robin, Robin>;
using RobinFluxCalculation = SimpleDynamics<DiffusionRelaxation<
    RobinFlux<DiffusionParticles, WallParticles, KernelGradientContact>>>;
//----------------------------------------------------------------------
//	An observer body to measure temperature at given positions.
//----------------------------------------------------------------------
class ParticleGeneratorTemperatureObserver : public ParticleGenerator<Observer>
{
  public:
    explicit ParticleGeneratorTemperatureObserver(SPHBody &sph_body) : ParticleGenerator<Observer>(sph_body)
    {
       /** A line of measuring points at the given position. */
       size_t number_of_observation_points = 5;
       Real range_of_measure = H - 0.02;
       Real start_of_measure = 0.01;

       for (size_t i = 0; i < number_of_observation_points; ++i)
       {
            Vec2d point_coordinate(
                0.028, range_of_measure * Real(i) / Real(number_of_observation_points - 1) + start_of_measure);
            positions_.push_back(point_coordinate);
       }
    }
};
#endif // WINDOWS_FRAME_DIFFUSION_D7_H