/**
 * @file 	2d_diffusion_test_with_NeumannBC.h
 * @brief 	This is the head files used by 2d_diffusion_test_with_NeumannBC.cpp.
 * @author	Chenxi Zhao, Bo Zhang, Chi Zhang and Xiangyu Hu
 */
#ifndef	DIFFUSION_TEST_WITH_NEUMANNBC_H
#define DIFFUSION_TEST_WITH_NEUMANNBC_H

#include "sphinxsys.h"
#include "particle_dynamics_diffusion_reaction_with_boundary.h"
#include "particle_dynamics_diffusion_reaction_with_boundary.hpp"

using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real L = 1.0;
Real H = 1.0;
Real resolution_ref = H / 100.0;
Real BW = resolution_ref * 2.0;
BoundingBox system_domain_bounds(Vec2d(-BW, -BW), Vec2d(L + BW, H + BW));
//----------------------------------------------------------------------
//	Basic parameters for material properties.
//----------------------------------------------------------------------
Real diffusion_coff = 1;
std::array<std::string, 1> species_name_list{ "Phi" };
//----------------------------------------------------------------------
//	Initial and boundary conditions.
//----------------------------------------------------------------------
Real initial_temperature = 0.0;
Real left_temperature = 300.0;
Real right_temperature = 350.0;
Real heat_flux = 900.0;
//----------------------------------------------------------------------
//	Geometric shapes used in the system.
//----------------------------------------------------------------------
std::vector<Vecd> createThermalDomain()
{
	std::vector<Vecd> thermalDomainShape;
	thermalDomainShape.push_back(Vecd(0.0, 0.0));
	thermalDomainShape.push_back(Vecd(0.0, H));
	thermalDomainShape.push_back(Vecd(L, H));
	thermalDomainShape.push_back(Vecd(L, 0.0));
	thermalDomainShape.push_back(Vecd(0.0, 0.0));

	return thermalDomainShape;
}

std::vector<Vecd> left_temperature_region
{
	Vecd(0.3 * L, H), Vecd(0.3 * L, H + BW), Vecd(0.4 * L, H + BW),
	Vecd(0.4 * L, H), Vecd(0.3 * L, H)
};

std::vector<Vecd> right_temperature_region
{
	Vecd(0.6 * L, H), Vecd(0.6 * L, H + BW), Vecd(0.7 * L, H + BW),
	Vecd(0.7 * L, H), Vecd(0.6 * L, H)
};

std::vector<Vecd> heat_flux_region
{
	Vecd(0.45 * L, -BW), Vecd(0.45 * L, 0), Vecd(0.55 * L, 0),
	Vecd(0.55 * L, -BW), Vecd(0.45 * L, -BW)
};

//----------------------------------------------------------------------
//	Define SPH bodies. 
//----------------------------------------------------------------------
class DiffusionBody : public MultiPolygonShape
{
public:
	explicit DiffusionBody(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(createThermalDomain(), ShapeBooleanOps::add);
	}
};

class WallBoundaryDirichlet : public MultiPolygonShape
{
public:
	explicit WallBoundaryDirichlet(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(left_temperature_region, ShapeBooleanOps::add);
		multi_polygon_.addAPolygon(right_temperature_region, ShapeBooleanOps::add);
	}
};

class WallBoundaryNeumann : public MultiPolygonShape
{
public:
	explicit WallBoundaryNeumann(const std::string& shape_name) : MultiPolygonShape(shape_name)
	{
		multi_polygon_.addAPolygon(heat_flux_region, ShapeBooleanOps::add);
	}
};

//----------------------------------------------------------------------
//	Setup diffusion material properties. 
//----------------------------------------------------------------------
class DiffusionMaterial : public DiffusionReaction<Solid>
{
public:
	DiffusionMaterial() : DiffusionReaction<Solid>({ "Phi" }, SharedPtr<NoReaction>())
	{
		initializeAnDiffusion<IsotropicDiffusion>("Phi", "Phi", diffusion_coff);
	}
};
using DiffusionParticles = DiffusionReactionParticles<SolidParticles, DiffusionMaterial>;
using WallParticles = DiffusionReactionParticles<SolidParticles, DiffusionMaterial>;
//----------------------------------------------------------------------
//	Application dependent initial condition. 
//----------------------------------------------------------------------
class DiffusionInitialCondition
	: public DiffusionReactionInitialConditionWithBoundary<DiffusionParticles>
{
protected:
	size_t phi_;

public:
	explicit DiffusionInitialCondition(SPHBody& sph_body)
		: DiffusionReactionInitialConditionWithBoundary<DiffusionParticles>(sph_body)
	{
		phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
	};

	void update(size_t index_i, Real dt)
	{
		all_species_[phi_][index_i] = 200 + 50 * (double)rand() / RAND_MAX;
	};
};

class WallBoundaryInitialCondition
	: public DiffusionReactionInitialConditionWithBoundary<WallParticles>
{
protected:
	size_t phi_;

public:
	WallBoundaryInitialCondition(SolidBody& diffusion_body) :
		DiffusionReactionInitialConditionWithBoundary<WallParticles>(diffusion_body)
	{
		phi_ = particles_->diffusion_reaction_material_.AllSpeciesIndexMap()["Phi"];
	}

	void update(size_t index_i, Real dt)
	{
		all_species_[phi_][index_i] = -0.0;
		if (pos_[index_i][1] > H && pos_[index_i][0] > 0.3 * L && pos_[index_i][0] < 0.4 * L)
		{
			all_species_[phi_][index_i] = left_temperature;
		}
		if (pos_[index_i][1] > H && pos_[index_i][0] > 0.6 * L && pos_[index_i][0] < 0.7 * L)
		{
			all_species_[phi_][index_i] = right_temperature;
		}
		if (pos_[index_i][1] < 0 && pos_[index_i][0] > 0.45 * L && pos_[index_i][0] < 0.55 * L)
		{
			heat_flux_[index_i] = heat_flux;
		}
	}
};

using DiffusionRelaxationInner = RelaxationOfAllDiffusionSpeciesInner<DiffusionParticles>;
using DiffusionRelaxationWithDirichletContact = RelaxationOfAllDiffusionSpeciesDirichletContact<DiffusionParticles, WallParticles>;
using DiffusionRelaxationWithNeumannContact = RelaxationOfAllDiffusionSpeciesNeumannContact<DiffusionParticles, WallParticles>;
//----------------------------------------------------------------------
//	Specify diffusion relaxation method. 
//----------------------------------------------------------------------
class DiffusionBodyRelaxation
	: public RelaxationOfAllDiffusionSpeciesRK2Complex<ComplexInteraction<DiffusionRelaxationInner, DiffusionRelaxationWithDirichletContact, DiffusionRelaxationWithNeumannContact>>
{
public:
	explicit DiffusionBodyRelaxation(InnerRelation& inner_relation, ContactRelation& body_contact_relation_Dirichlet, ContactRelation& body_contact_relation_Neumann)
		: RelaxationOfAllDiffusionSpeciesRK2Complex<ComplexInteraction<DiffusionRelaxationInner, DiffusionRelaxationWithDirichletContact, DiffusionRelaxationWithNeumannContact>>(inner_relation, body_contact_relation_Dirichlet, body_contact_relation_Neumann) {};
	virtual ~DiffusionBodyRelaxation() {};
};
//----------------------------------------------------------------------
//	An observer body to measure temperature at given positions. 
//----------------------------------------------------------------------
class TemperatureObserverParticleGenerator : public ObserverParticleGenerator
{
public:
	TemperatureObserverParticleGenerator(SPHBody& sph_body) : ObserverParticleGenerator(sph_body)
	{
		/** A line of measuring points at the middle line. */
		size_t number_of_observation_points = 5;
		Real range_of_measure = L;
		Real start_of_measure = 0;

		for (size_t i = 0; i < number_of_observation_points; ++i)
		{
			Vec2d point_coordinate(0.5 * L, range_of_measure * Real(i) /
				Real(number_of_observation_points - 1) + start_of_measure);
			positions_.push_back(point_coordinate);
		}
	}
};

#endif //DIFFUSION_TEST_WITH_NEUMANNBC_H