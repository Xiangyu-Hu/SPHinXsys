/**
 * @file 	2d_turbulent_channel_PBC.h
 * @brief 	This is the case file for the test of flow passing by a cylinder.
 * @details  We consider a flow passing by a cylinder in 2D.
 * @author 	Xiangyu Hu, Chi Zhang and Luhui Han
 */

#include "sphinxsys.h"
//#include "k-epsilon_turbulent_model.cpp"
using namespace SPH;

//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DH = 2.0;                         /**< Channel height. */
Real DL = 15.0;                         /**< Channel length. */
Real num_fluid_cross_section = 40.0;
Real amplitude = 0.1;
Real wave_length = 1;


//----------------------------------------------------------------------
//	Unique parameters for turbulence. 
//----------------------------------------------------------------------
Real characteristic_length = DH; /**<It needs characteristic Length to calculate turbulent length and the inflow turbulent epsilon>*/
//** Intial values for K, Epsilon and Mu_t *
StdVec<Real> initial_turbu_values = { 0.000180001 ,3.326679e-5 ,1.0e-9 };  

Real y_p_constant = 0.05;
Real resolution_ref = (DH - 2.0 * y_p_constant) / (num_fluid_cross_section - 1.0); /**< Initial reference particle spacing. */
Real offset_distance = y_p_constant - resolution_ref / 2.0; //** Basically offset distance is large than or equal to 0 *

Real BW = resolution_ref * 4;         /**< Reference size of the emitter. */
Real DL_sponge = resolution_ref * 20;
Real half_channel_height = DH / 2.0;
//----------------------------------------------------------------------
//	Domain bounds of the system. 
//----------------------------------------------------------------------
BoundingBox system_domain_bounds(Vec2d(-DL_sponge - 2.0 * BW, -BW), Vec2d(DL + 2.0 * BW, DH + 2.0 * BW));

//----------------------------------------------------------------------
//	Material properties of the fluid.
//----------------------------------------------------------------------
//** Laminar *
Real U_inlet = 0.5;
Real U_max = 0.75;
Real U_f = U_inlet; //*Characteristic velo is regarded as average velo here
Real c_f = 10.0 * U_max;                                        /**< Speed of sound. */
Real rho0_f = 1.0;                                            /**< Density. */
Real mu_f = 0.01;

Real Re_calculated = U_f * DH * rho0_f / mu_f; 
//----------------------------------------------------------------------
//	The emitter block. 
//----------------------------------------------------------------------
Real DH_C = DH - 2.0 * offset_distance;
Vec2d emitter_halfsize = Vec2d(0.5 * BW, 0.5 * DH_C);
Vec2d emitter_translation = Vec2d(-DL_sponge, 0.0) + emitter_halfsize + Vecd(0.0, offset_distance);
Vec2d inlet_buffer_halfsize = Vec2d(0.5 * DL_sponge, 0.5 * DH_C);
Vec2d inlet_buffer_translation = Vec2d(-DL_sponge, 0.0) + inlet_buffer_halfsize + Vecd(0.0, offset_distance);

Vec2d disposer_halfsize = Vec2d(0.5 * BW, 0.75 * DH);
Vec2d disposer_translation = Vec2d(DL, DH + 0.25 * DH) - disposer_halfsize;
//----------------------------------------------------------------------
// Observation.
//----------------------------------------------------------------------
Real x_observe = 0.90 * DL;
Real x_observe_start = 0.90 * DL;
Real observe_spacing_x = 0.02 * DL;
int num_observer_points_x = 1;
int num_observer_points = std::round(DH / resolution_ref); //**Evrey particle is regarded as a cell monitor* 
StdVec<Real> monitoring_bound = { 109 ,111 };
Real observe_spacing = DH / num_observer_points;
StdVec<Vecd> observation_locations;
//----------------------------------------------------------------------
//	Cases-dependent geometries 
//----------------------------------------------------------------------
std::vector<Vecd> createWaterBlockShape()
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance, DH));
    water_block_shape.push_back(Vecd(DL + offset_distance, DH));
    water_block_shape.push_back(Vecd(DL + offset_distance, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - offset_distance,  0.0));
    return water_block_shape;
    for (int k = 0; k <=0; ++k)
    {
       
    }
}
class WaterBlock : public ComplexShape
{
public:
    explicit WaterBlock(const std::string& shape_name) : ComplexShape(shape_name)
    {
        //MultiPolygon outer_boundary(createWaterBlockShape());
        //add<MultiPolygonShape>(outer_boundary, "OuterBoundary");
        MultiPolygon fluid_domain(createWaterBlockShape());
        add<ExtrudeShape<MultiPolygonShape>>(-offset_distance, fluid_domain, "FluidDomain");

    }
};

std::vector<Vecd> createChannelWallShape() //outerwall
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, DH));
    water_block_shape.push_back(Vecd(DL, DH));
    water_block_shape.push_back(Vecd(DL, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge, 0.0));

    return water_block_shape;
}
std::vector<Vecd> createComputationalDomainShape()// inner wall
{
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, DH));
    water_block_shape.push_back(Vecd(DL + 2.0 * BW, DH));
    water_block_shape.push_back(Vecd(DL + 2.0 * BW, 0.0));
    water_block_shape.push_back(Vecd(-DL_sponge - 2.0 * BW, 0.0));

    return water_block_shape;
}

/**
 * @brief 	Wall boundary body definition.
 */
class WallBoundary : public ComplexShape
{
public:
    explicit WallBoundary(const std::string& shape_name) : ComplexShape(shape_name)
    {
        MultiPolygon outer_boundary(createChannelWallShape());
        add<ExtrudeShape<MultiPolygonShape>>(-offset_distance + BW, outer_boundary, "OuterBoundary");

        MultiPolygon computational_domain(createComputationalDomainShape());
        subtract<ExtrudeShape<MultiPolygonShape>>(-offset_distance, computational_domain, "ComputationalDomain");
    }
};

//----------------------------------------------------------------------
//	Inflow velocity
//----------------------------------------------------------------------
struct InflowVelocity
{
    Real u_ref_, t_ref_;
    AlignedBoxShape& aligned_box_;
    Vecd halfsize_;

    template <class BoundaryConditionType>
    InflowVelocity(BoundaryConditionType& boundary_condition)
        : u_ref_(U_inlet), t_ref_(2.0),
        aligned_box_(boundary_condition.getAlignedBox()),
        halfsize_(aligned_box_.HalfSize()) {}

    Vecd operator()(Vecd& position, Vecd& velocity)
    {
        Vecd target_velocity = velocity;
        Real run_time = GlobalStaticVariables::physical_time_;
        Real u_ave = run_time < t_ref_ ? 0.5 * u_ref_ * (1.0 - cos(Pi * run_time / t_ref_)) : u_ref_;
        //target_velocity[0] = 1.5 * u_ave * SMAX(0.0, 1.0 - position[1] * position[1] / halfsize_[1] / halfsize_[1]);
        target_velocity[0] = 1.5 * u_ave * (1.0 - position[1] * position[1] / half_channel_height / half_channel_height);
        //target_velocity[0] = u_ave;
        if (position[1] > half_channel_height)
        {
            std::cout << "Particles out of domain, wrong inlet velocity." << std::endl;
            system("pause");
        }
        target_velocity[1] = 0.0;
        return target_velocity;
    }
};
//----------------------------------------------------------------------
//	Define time dependent acceleration in x-direction
//----------------------------------------------------------------------
class TimeDependentAcceleration : public Gravity
{
    Real du_ave_dt_, u_ref_, t_ref_;

public:
    explicit TimeDependentAcceleration(Vecd gravity_vector)
        : Gravity(gravity_vector), t_ref_(2.0), u_ref_(U_inlet), du_ave_dt_(0) {}

    virtual Vecd InducedAcceleration(const Vecd& position) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        du_ave_dt_ = 0.5 * u_ref_ * (Pi / t_ref_) * sin(Pi * run_time_ / t_ref_);
        return run_time_ < t_ref_ ? Vecd(du_ave_dt_, 0.0) : global_acceleration_;
    }
};