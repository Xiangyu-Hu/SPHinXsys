/**
 * @file 	static_confinement.cpp
 * @brief 	2D dambreak example in which the solid wall boundary are static confinement.
 * @details This is the one of the basic test cases.
 * @author 	Xiangyu Hu
 */
#include "sphinxsys.h" //SPHinXsys Library.
#include "body_part_by_cell_tracing.h"
#include "level_set_confinement.h"
#include "math.h"
#include <iostream>
#include <fstream>
using namespace SPH;   // Namespace cite here.
//----------------------------------------------------------------------
//	Basic geometry parameters and numerical setup.
//----------------------------------------------------------------------
Real DL = 5.366;              /**< Tank length. */
Real DH = 5.366;              /**< Tank height. */
Real LL = 5.366;                /**< Liquid column length. */
Real LH = 2.0;                /**< Liquid column height. */
Real resolution_ref = 0.05;  /**< Global reference resolution. */


// circle parameters
Vecd insert_circle_center (2.0, 1.0);
Real insert_circle_radius = 0.25;


std::vector<Vecd> createWaterBlockShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.0, 0.0));
    water_block_shape.push_back(Vecd(0.0, LH));
    water_block_shape.push_back(Vecd(LL, LH));
    water_block_shape.push_back(Vecd(LL, 0.0));
    water_block_shape.push_back(Vecd(0.0, 0.0));
    return water_block_shape;
}
/** create wall shape */
std::vector<Vecd> createWallShape()
{
    std::vector<Vecd> inner_wall_shape;
    inner_wall_shape.push_back(Vecd(0.0, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, DH));
    inner_wall_shape.push_back(Vecd(DL, DH));
    inner_wall_shape.push_back(Vecd(DL, 0.0));
    inner_wall_shape.push_back(Vecd(0.0, 0.0));

    return inner_wall_shape;
}
/** create a structure shape */
std::vector<Vecd> createStructureShape()
{
    // geometry
    std::vector<Vecd> water_block_shape;
    water_block_shape.push_back(Vecd(0.5 * DL, 0.05 * DH));
    water_block_shape.push_back(Vecd(0.5 * DL + 0.5 * LL, 0.05 * DH + 0.5 * LH));
    water_block_shape.push_back(Vecd(0.5 * DL + 0.5 * LL, 0.05 * DH));
    water_block_shape.push_back(Vecd(0.5 * DL, 0.05 * DH));
    return water_block_shape;
}

std::vector<Vecd> creatSquare()
{
    //geometry
    std::vector<Vecd> square_shape;
    square_shape.push_back(Vecd(2.5, 0.75));
    square_shape.push_back(Vecd(2.5, 1.25));
    square_shape.push_back(Vecd(3.0, 1.25));
    square_shape.push_back(Vecd(3.0, 0.75));
    square_shape.push_back(Vecd(2.5, 0.75));
    return square_shape;
}

//----------------------------------------------------------------------
// Water body shape definition.
//----------------------------------------------------------------------
class WaterBlock : public MultiPolygonShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWaterBlockShape(), ShapeBooleanOps::add);
        //multi_polygon_.addAPolygon(createStructureShape(), ShapeBooleanOps::sub);
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::sub);
        //multi_polygon_.addAPolygon(creatSquare(), ShapeBooleanOps::sub);
    }
};
//----------------------------------------------------------------------
//	Shape for the wall.
//----------------------------------------------------------------------
class Wall : public MultiPolygonShape
{
  public:
    explicit Wall(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
        multi_polygon_.addAPolygon(createWallShape(), ShapeBooleanOps::add);
    }
};
//----------------------------------------------------------------------
//	Shape for a structure.
//----------------------------------------------------------------------
class Triangle : public MultiPolygonShape
{
  public:
    explicit Triangle(const std::string &shape_name) : MultiPolygonShape(shape_name)
    {
       // multi_polygon_.addAPolygon(createStructureShape(), ShapeBooleanOps::add);
        multi_polygon_.addACircle(insert_circle_center, insert_circle_radius, 100, ShapeBooleanOps::add);
        //multi_polygon_.addAPolygon(creatSquare(), ShapeBooleanOps::add);
    }
};

class HorizontalMovement: public BaseTracingMethod
{
 public:
     HorizontalMovement() {};
     virtual ~HorizontalMovement(){};

     virtual Vecd tracingPosition (Vecd previous_position, Real current_time = 0.0) override
     {
         Real run_time = GlobalStaticVariables::physical_time_;
         Vecd current_position (0.0, 0.0);
         current_position[0]= previous_position[0] - 0.2 * run_time;
         //current_position[0] = previous_position[0];
         if(run_time <= 10.0)
         {
             current_position[1] = previous_position[1] + 0.05 * run_time;
         }
         else{
             current_position[1] = previous_position[1] - 0.05 * (run_time-20.0);
         }
         
         

         return current_position;
     }
};

class CircleMovement : public BaseTracingMethod
{
public:
    CircleMovement() {};
    virtual ~CircleMovement() {};

   
    virtual Vecd tracingPosition(Vecd previous_position, Real current_time = 0.0) override
    {
        Real dt = 0.1;
        Vecd rotation_center(2.0, 1.0);
        Real rotation_v = 0.2*Pi;
        Real rho = (previous_position - rotation_center).norm();
        Real theta = atan2(previous_position[0] - rotation_center[0], previous_position[1] - rotation_center[1]);
        Real run_time = GlobalStaticVariables::physical_time_;
        Vecd current_position(0.0, 0.0);
        current_position[0] = rotation_center[0] + cos(theta + rotation_v * run_time) * rho;
        current_position[1] = rotation_center[1] + sin(theta + rotation_v * run_time) * rho;
       
        return current_position;
    }
};
//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    BoundingBox mesh_domain(Vecd(1.5, 0.5), Vecd(2.5, 1.5));
    SPHSystem sph_system(mesh_domain, resolution_ref);
    IOEnvironment io_environment(sph_system);
    Mesh test_region(mesh_domain, 0.2, 4);
    
    CircleMovement circle_movement;
    Real end_time = 10.0;
    
    while (GlobalStaticVariables::physical_time_ <= end_time)
    {
        std::string output_folder_ = "./output";
        std::string filePath = output_folder_ + "/" + "rotation_test_" + std::to_string(GlobalStaticVariables::physical_time_) + ".dat";
        std::ofstream outputFile(filePath.c_str(), std::ios::app);
        Real dt = 0.5;
        Array2i cells = test_region.AllGridPoints();
        size_t x = cells[0];
        size_t y = cells[1];
        size_t n = cells.size();
        std::vector<Vecd> cell_position;
        cell_position.clear();
        for(size_t i= 0; i != x; ++i)
            for(size_t j=0; j!= y; ++j)
        {
                //cell_position.push_back(test_region.CellPositionFromIndex(Array2i(i, j)));
                cell_position.push_back(circle_movement.tracingPosition(test_region.CellPositionFromIndex(Array2i(i, j))));
        }
        
        for(size_t n=0; n!= cell_position.size(); ++n)
        {
            outputFile << cell_position[n][0] << " " << cell_position[n][1] << std::endl;
        }
  
        GlobalStaticVariables::physical_time_ += dt;
    }
    return 0;
}
