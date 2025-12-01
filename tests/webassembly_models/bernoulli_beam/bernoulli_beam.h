/**
 * @file total_artificial_heart.cpp
 * @brief This is the example of total artificial heart implantation path simulation
 * @author John Benjamin, Bence Rochlitz - Virtonomy GmbH
 */
#ifndef SIM_BROWSER_BEAM_H
#define SIM_BROWSER_BEAM_H

#include "structural_simulation_class.h"
using namespace SPH;

#ifdef __EMSCRIPTEN__
#include <emscripten/val.h>
#endif

struct BernoulliBeamInput
{
    Real scale_stl;
    std::vector<Real> resolution;
    Real rho_0;
    Real poisson;
    Real Youngs_modulus;
    Real physical_viscosity;
    std::array<Real, 3> translation;
    StlList stls;
    std::string relative_input_path;
};

StructuralSimulationInput createSimulationInput(const BernoulliBeamInput &input, std::SharedPtr<SaintVenantKirchhoffSolid> material)
{
    StlList imported_stl_list = input.stls;
    std::vector<Vec3d> translation_list = {
        Vec3d(input.translation[0], input.translation[1], input.translation[2])};
    std::vector<Real> resolution_list = input.resolution;

    std::vector<SharedPtr<SaintVenantKirchhoffSolid>> material_model_list = {
        material};

    /** INPUT DECLARATION */
    StructuralSimulationInput inputStructuralSim = {
        input.relative_input_path,
        imported_stl_list,
        input.scale_stl,
        translation_list,
        resolution_list,
        material_model_list,
        StdVec<Real>(1, input.physical_viscosity),
        {}};
    inputStructuralSim.non_zero_gravity_ = std::vector<GravityPair>{GravityPair(0, Vec3d(0.0, -100.0, 0.0))}; // gravity

    BoundingBoxd fixation(-0.1 * Vec3d::Ones(), Vec3d(0, 0.1, 0.1));
    inputStructuralSim.body_indices_fixed_constraint_region_ = StdVec<ConstrainedRegionPair>{ConstrainedRegionPair(0, fixation)};
    inputStructuralSim.particle_relaxation_list_ = {true};

    return inputStructuralSim;
}

class BernoulliBeam
{
  public:
    BernoulliBeam(const BernoulliBeamInput &input) : material_(make_shared<SaintVenantKirchhoffSolid>(input.rho_0, input.Youngs_modulus, input.poisson))
    {
        sim.reset(new StructuralSimulation(createSimulationInput(input, material_)));
    }

  public: // C++ Backend functions
    void runCompleteSimulation(Real endTime) { sim->runSimulation(SPH::Real(endTime)); };

  private:
    std::SharedPtr<SaintVenantKirchhoffSolid> material_;
    std::unique_ptr<StructuralSimulation> sim;
};

#ifdef __EMSCRIPTEN__

class BernoulliBeamJS
{
  public:
    BernoulliBeamJS(const BernoulliBeamInput &input) :

                                                       material_(make_shared<SaintVenantKirchhoffSolid>(input.rho_0, input.Youngs_modulus, input.poisson))
    {
        sim_js_.reset(new StructuralSimulationJS(createSimulationInput(input, material_)));
    }

    void runSimulation(int number_of_steps)
    {
        try
        {
            sim_js_->runSimulationFixedDuration(number_of_steps);
        }
        catch (const std::exception &e)
        {
            if (on_error_)
                on_error_(e.what());
        }
    }

    VtuStringData getVtuData() const { return sim_js_->getVtuData(); }

    void onError(emscripten::val on_error)
    {
        on_error_ = [on_error](const std::string &error_message)
        { on_error(error_message); };
    }

  private:
    std::SharedPtr<SaintVenantKirchhoffSolid> material_;
    std::unique_ptr<StructuralSimulationJS> sim_js_;
    std::function<void(const std::string &)> on_error_;
};

#endif //__EMSCRIPTEN__

#endif // SIM_BROWSER_BEAM_H