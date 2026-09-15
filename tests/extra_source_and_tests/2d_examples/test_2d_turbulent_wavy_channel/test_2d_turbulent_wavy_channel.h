#include "bidirectional_buffer.h"
#include "udf_common_turbulence_model.cpp"
#include "density_correction.h"
#include "density_correction.hpp"
#include "udf_k-omega_turbulent_model.cpp"
#include "kernel_summation.h"
#include "kernel_summation.hpp"
#include "pressure_boundary.h"
#include "udf_P_refinement.h"
#include "sphinxsys.h"
using namespace SPH;

Real DH = 1.0;
Real DL = 1.0;
Real num_fluid_cross_section = 40.0;

constexpr Real wave_amplitude = 0.1;
constexpr Real wave_length = 1.0;
constexpr Real pi = 3.14159265358979323846;
Real lowerWallHeight(Real x)
{
    return -wave_amplitude *
        std::sin(2.0 * pi * x / wave_length);
}

Vecd external_acc = Vecd(0.01687141, 0.0);
Real external_acc_gradually_impose_t = 2.0;
int is_blended = 0;
int is_AMRD = 0;
bool is_constrain_normal_velocity_in_P_region = false;
bool is_source_term_linearisation = false;
static constexpr int num_node_sublayer_model = 5;
static constexpr int type_tdma_sublayer_model = 5;
Real turbulent_module_activate_time = 2.0;
StdVec<Real> initial_turbu_values = {0.01, 2.056, 0.02};
Real y_p_constant = DH / 2.0 / num_fluid_cross_section;
Real resolution_ref = DH / num_fluid_cross_section;

Real BW = resolution_ref * 4;

BoundingBoxd system_domain_bounds(Vecd(-2.0 * BW, -wave_amplitude - BW), Vecd(DL + 2.0 * BW, DH + BW));

Real U_inlet = 0.816;
Real U_f = U_inlet;
Real U_max = 1.5 * U_inlet;
Real c_f = 10.0 * U_max;
Real rho0_f = 1.0;
Real mu_f = 1.0e-4;
Real Re_calculated = U_f * DH * rho0_f / mu_f;

Real x_observe_center = 0.5 * DL;
Real y_observe_center =
0.5 * (DH + lowerWallHeight(x_observe_center));

StdVec<Vecd> observer_location_center_point = {
    Vecd(x_observe_center, y_observe_center) };

std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    constexpr size_t number_of_segments = 200;
    water_block_shape.push_back(Vecd(0.0, lowerWallHeight(0.0)));
    water_block_shape.push_back(Vecd(0.0, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    for (size_t i = 0; i <= number_of_segments; ++i)
    {
        Real x = DL *
            static_cast<Real>(number_of_segments - i) /
            static_cast<Real>(number_of_segments);

        water_block_shape.push_back(Vecd(x, lowerWallHeight(x)));
    }
    return water_block_shape;
}
class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name)
        : ComplexShape(shape_name)
    {
        MultiPolygon computational_domain(createWaterBlockShape());
        add<MultiPolygonShape>(
            computational_domain, "ComputationalDomain");
    }
};
std::vector<Vecd> createUpperWallShape()
{
    Real wall_x_min = -2.0 * BW;
    Real wall_x_max = DL + 2.0 * BW;

    std::vector<Vecd> upper_wall_shape;

    upper_wall_shape.push_back(Vecd(wall_x_min, DH));
    upper_wall_shape.push_back(Vecd(wall_x_min, DH + BW));
    upper_wall_shape.push_back(Vecd(wall_x_max, DH + BW));
    upper_wall_shape.push_back(Vecd(wall_x_max, DH));
    upper_wall_shape.push_back(Vecd(wall_x_min, DH));

    return upper_wall_shape;
}
std::vector<Vecd> createLowerWallShape()
{
    std::vector<Vecd> lower_wall_shape;

    Real wall_x_min = -2.0 * BW;
    Real wall_x_max = DL + 2.0 * BW;

    constexpr size_t number_of_segments = 240;
    for (size_t i = 0; i <= number_of_segments; ++i)
    {
        Real x = wall_x_min +
            (wall_x_max - wall_x_min) *
            static_cast<Real>(i) /
            static_cast<Real>(number_of_segments);

        lower_wall_shape.push_back(
            Vecd(x, lowerWallHeight(x)));
    }
    for (size_t i = 0; i <= number_of_segments; ++i)
    {
        Real x = wall_x_max -
            (wall_x_max - wall_x_min) *
            static_cast<Real>(i) /
            static_cast<Real>(number_of_segments);

        lower_wall_shape.push_back(
            Vecd(x, lowerWallHeight(x) - BW));
    }
    lower_wall_shape.push_back(
        Vecd(wall_x_min, lowerWallHeight(wall_x_min)));

    return lower_wall_shape;
}
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string& shape_name)
        : ComplexShape(shape_name)
    {
        MultiPolygon upper_wall(createUpperWallShape());
        add<MultiPolygonShape>(upper_wall, "UpperWall");

        MultiPolygon lower_wall(createLowerWallShape());
        add<MultiPolygonShape>(lower_wall, "LowerWavyWall");
    }
};
namespace SPH
{

    class UpdateVolume : public LocalDynamics
    {
    public:
        explicit UpdateVolume(SPHBody& sph_body)
            : LocalDynamics(sph_body),
            rho_(particles_->getVariableDataByName<Real>("Density")),
            mass_(particles_->getVariableDataByName<Real>("Mass")),
            Vol_(particles_->getVariableDataByName<Real>("VolumetricMeasure")),
            indicator_(particles_->getVariableDataByName<int>("Indicator")) {}
        virtual ~UpdateVolume() {};

        void update(size_t index_i, Real dt = 0.0)
        {
            Vol_[index_i] = mass_[index_i] / rho_[index_i];
        };

    protected:
        Real* rho_, * mass_, * Vol_;
        int* indicator_;
    };

}
