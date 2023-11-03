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
// general parameters for geometry
Real resolution_ref = 0.5;   // particle spacing
Real BW = resolution_ref * 4; // boundary width
Real DL = 3.0;              // tank length
Real DH = 3.0;                // tank height
Real DW = 1.0;                // tank width
Real LL = 3.0;                // liquid length
Real LH = 3.0;                // liquid height
Real LW = 1.0;                // liquid width
Real CL = 1.0;                // liquid length
Real CH = 1.0;                // liquid height
Real CW = 1.0;                // liquid width

Vecd axis_point_1(0.0, 0.0, 0.0);
Vecd axis_point_2(0.0, 0.0, 1.0);

// for material properties of the fluid
Real rho0_f = 1.0;
Real gravity_g = 1.0;
Real U_f = 2.0 * sqrt(gravity_g * LH);
Real c_f = 10.0 * U_f;

class WaterBlock : public ComplexShape
{
  public:
    explicit WaterBlock(const std::string &shape_name) : ComplexShape(shape_name)
    {
        Vecd halfsize_water(0.5 * LL, 0.5 * LH, 0.5 * LW);
        Transform translation_water(halfsize_water);
        Transform translation_start_point(Vecd(0.0, 0.0, 0.0));
        add<GeometricShapeBox>(halfsize_water);
        //add<TransformShape<GeometricShapeBox>>(Transform(translation_water), halfsize_water);
        //subtract<TriangleMeshShapeSTL>(full_path_to_file, translation, scaling);
        //subtract<TransformShape<GeometricShapeBox>>(Transform(translation_start_point), halfsize_cubic);
    }
};

class RotationMovement : public BaseTracingMethod
{
public:
    RotationMovement(Vecd axis_point_1, Vecd axis_point_2, Real rotation_velocity) : axis_point_1_(axis_point_1), 
        axis_point_2_(axis_point_2), rotation_v_(rotation_velocity)
    {
        rotation_axis_[0] = axis_point_1_[0] - axis_point_2_[0];
        rotation_axis_[1] = axis_point_1_[1] - axis_point_2_[1];
        rotation_axis_[2] = axis_point_1_[2] - axis_point_2_[2];
    };
    virtual ~RotationMovement() {};
   
    virtual Vecd tracingPosition(Vecd previous_position, Real current_time = 0.0) override
    {
        Real rho = (previous_position - findProjectionPoint(previous_position)).norm();
        Real difference_0 = previous_position[0] - findProjectionPoint(previous_position)[0];
        Real difference_1 = previous_position[1] - findProjectionPoint(previous_position)[1];
        Real difference_2 = previous_position[2] - findProjectionPoint(previous_position)[2];
        Real theta = atan2(difference_0, difference_1);
        Real run_time = GlobalStaticVariables::physical_time_;
        Vecd current_position(0.0, 0.0, 0.0);
        current_position[0] = rotation_axis_[0] + cos(theta + rotation_v_ * run_time) * rho;
        current_position[2] = previous_position[2];
        current_position[1] = rotation_axis_[1] + sin(theta + rotation_v_ * run_time) * rho;
       
        return current_position;
    }

    virtual Vecd updateNormalForVector(Vecd previous_vector) override
    {
        Real run_time = GlobalStaticVariables::physical_time_;
        Real magnitude = (previous_vector - findProjectionPoint(previous_vector)).norm();
        Real theta = atan2(previous_vector[0], previous_vector[2]);
        Vecd current_vector(0.0, 0.0, 0.0);
        current_vector[0] = magnitude * cos(theta + run_time * rotation_v_);
        current_vector[1] = magnitude * sin(theta + run_time * rotation_v_);
        current_vector[2] = previous_vector[2];

        return current_vector;
    }

protected:
    Vecd findProjectionPoint(Vecd previous_position)
    {
        Vecd point_vector (0.0,0.0,0.0);
        Vecd projection_vector = Vecd::Zero();
        point_vector[0] = previous_position[0] - axis_point_1_[0];
        point_vector[1] = previous_position[1] - axis_point_1_[1];
        point_vector[2] = previous_position[2] - axis_point_1_[2];

        Real axisLengthSquared = projection_vector[0] * rotation_axis_[0] + rotation_axis_[1] * rotation_axis_[1] + rotation_axis_[2] * rotation_axis_[2];
        Real dotProduct = rotation_axis_[0] * point_vector[0] + rotation_axis_[1] * point_vector[1] + rotation_axis_[2] * point_vector[2];

        projection_vector[0] = rotation_axis_[0] * (dotProduct / axisLengthSquared);
        projection_vector[1] = rotation_axis_[1] * (dotProduct / axisLengthSquared);
        projection_vector[2] = rotation_axis_[2] * (dotProduct / axisLengthSquared);

        return projection_vector;
    }
    Vecd rotation_axis_; //two points axis
    Real rotation_v_;
    Vecd axis_point_1_;
    Vecd axis_point_2_;
};

class ThreeDRotation : public BaseTracingMethod
{
public:
    ThreeDRotation(Vecd axis_point_1, Vecd axis_point_2, Real rotation_velocity) :
        axis_point_A_(axis_point_1), angular_v_(rotation_velocity)
    {
        axis_ = (axis_point_1 - axis_point_2).normalized();
        
    }

    virtual Vecd tracingPosition(Vecd previous_position, Real current_time = 0.0) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        Eigen::Quaterniond rotation_position;
        rotation_position = Eigen::AngleAxisd(angular_v_ * run_time_, axis_);
        Eigen::Quaterniond point_quaternion_previos(0, previous_position.x() - axis_point_A_.x(), previous_position.y() - axis_point_A_.y(), previous_position.z() - axis_point_A_.z());
        Eigen::Quaterniond rotated_point_quaternion = rotation_position * point_quaternion_previos * rotation_position.conjugate();
        Vecd new_position(rotated_point_quaternion.x() + axis_point_A_.x(),
            rotated_point_quaternion.y() + axis_point_A_.y(),
            rotated_point_quaternion.z() + axis_point_A_.z());
        return new_position;
    }

    virtual Vecd updateNormalForVector(Vecd previous_vector) override
    {
        Real run_time_ = GlobalStaticVariables::physical_time_;
        Eigen::Quaterniond rotation_vector;
        rotation_vector = Eigen::AngleAxisd(angular_v_ * run_time_, axis_);
        Eigen::Quaterniond vector_quaternion_previous(0, previous_vector.x(), previous_vector.y(), previous_vector.z());
        Eigen::Quaterniond rotated_quaternion = rotation_vector * vector_quaternion_previous * rotation_vector.conjugate();
        Vecd new_vector(rotated_quaternion.x(), rotated_quaternion.y(), rotated_quaternion.z());
        return new_vector;
    }
protected:
    Vecd axis_;
    Vecd axis_point_A_;
    Real angular_v_;
    
    
};

//----------------------------------------------------------------------
//	Main program starts here.
//----------------------------------------------------------------------
int main(int ac, char *av[])
{
    BoundingBox system_domain(Vecd(-5.0, -5.0, -5.0), Vecd(5.0, 5.0, 5.0));
    SPHSystem sph_system(system_domain, resolution_ref);
    IOEnvironment io_environment(sph_system);
    /*BoundingBox mesh_domain(Vecd(-2.0, -2.0, -0.5), Vecd(2.0, 2.0, 0.5));
    Mesh test_region(mesh_domain, 0.5, 0);*/
    
    FluidBody water_block(sph_system, makeShared<WaterBlock>("WaterBody"));
    water_block.defineParticlesAndMaterial<BaseParticles, WeaklyCompressibleFluid>(rho0_f, c_f);
    water_block.generateParticles<ParticleGeneratorLattice>();

    RotationMovement circle_movement(axis_point_1, axis_point_2, 0.5* Pi);
    ThreeDRotation rotation(axis_point_1, axis_point_2, 0.5* Pi);
    Real end_time = 10.0;
    
    while (GlobalStaticVariables::physical_time_ <= end_time)
    {
        std::string output_folder_ = "./output";
        std::string filePath_1 = output_folder_ + "/" + "rotation_test_" + std::to_string(GlobalStaticVariables::physical_time_) + ".dat";
        std::ofstream outputFile_1(filePath_1.c_str(), std::ios::app);
        std::string filePath_2 = output_folder_ + "/" + "rotation_test_orginal_position_" + ".dat";
        std::ofstream outputFile_2(filePath_2.c_str(), std::ios::app);
        Real dt = 0.2;
        //Array3i cells = test_region.AllGridPoints();
        //size_t x = cells[0];
        //size_t y = cells[1];
        //size_t z = cells[2];
        //size_t n = cells.size();
        //std::vector<Vecd> cell_position;
        //cell_position.clear();
        //for(size_t i= 0; i != x; ++i)
        //    for(size_t j=0; j!= y; ++j)
        //        for(size_t k=0; k!= z; ++k)
        //{
        //        cell_position.push_back(test_region.CellPositionFromIndex(Array3i(i, j, k)));
        //        //cell_position.push_back(circle_movement.tracingPosition(test_region.CellPositionFromIndex(Array3i(i, j, k))));
        //}
        std::vector<Vecd> new_position;
        new_position.clear();
        /*for (size_t n = 0; n !=  water_block.getBaseParticles().pos_.size(); ++n)
        {
            new_position.push_back(circle_movement.tracingPosition(water_block.getBaseParticles().pos_[n], GlobalStaticVariables::physical_time_));
            outputFile_2 << water_block.getBaseParticles().pos_[n][0] << " " << water_block.getBaseParticles().pos_[n][1] << " " << water_block.getBaseParticles().pos_[n][2] << std::endl;
        };*/

        for (size_t n = 0; n !=  water_block.getBaseParticles().pos_.size(); ++n)
        {
            new_position.push_back(rotation.tracingPosition(water_block.getBaseParticles().pos_[n], GlobalStaticVariables::physical_time_));
            outputFile_2 << water_block.getBaseParticles().pos_[n][0] << " " << water_block.getBaseParticles().pos_[n][1] << " " << water_block.getBaseParticles().pos_[n][2] << std::endl;
        };


        for(size_t n=0; n!= water_block.getBaseParticles().pos_.size(); ++n)
        {
            outputFile_1 << new_position[n][0] << " " << new_position[n][1] << " " << new_position[n][2]<<std::endl;
            //outputFile << water_block.getBaseParticles().pos_[n][0] << " " << water_block.getBaseParticles().pos_[n][1] << " " << water_block.getBaseParticles().pos_[n][2]<<std::endl;
        }
  
        GlobalStaticVariables::physical_time_ += dt;
    }
    return 0;
}
